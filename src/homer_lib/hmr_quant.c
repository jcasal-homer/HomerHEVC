/*****************************************************************************
 * hmr_quant.c : homerHEVC encoding library
/*****************************************************************************
 * Copyright (C) 2014 homerHEVC project
 *
 * Juan Casal <jcasal@homerhevc.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *****************************************************************************/
/*
 * some of the work below is derived from HM HEVC reference code where 
 * the following license applies
/****************************************************************************
/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.  
 *
 * Copyright (c) 2010-2014, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *****************************************************************************/

#include <memory.h>
#include <math.h>
#include <limits.h>

#include "hmr_common.h"
#include "hmr_private.h"




#define SCAN_SET_SIZE		16
#define LOG2_SCAN_SET_SIZE	4
#define SBH_THRESHOLD		4

void sign_bit_hidding( short * dst, short * src, uint const *scan, short* deltaU, int width, int height )
{
	int lastCG = -1;
	int absSum = 0 ;
	int n ;
	int subSet;
	for( subSet = (width*height-1) >> LOG2_SCAN_SET_SIZE; subSet >= 0; subSet-- )
	{
		int  subPos     = subSet << LOG2_SCAN_SET_SIZE;
		int  firstNZPosInCG=SCAN_SET_SIZE , lastNZPosInCG=-1 ;
		absSum = 0 ;

		for(n = SCAN_SET_SIZE-1; n >= 0; --n )
		{
			if( dst[ scan[ n + subPos ]] )
			{
				lastNZPosInCG = n;
				break;
			}
		}

		for(n = 0; n <SCAN_SET_SIZE; n++ )
		{
			if( dst[ scan[ n + subPos ]] )
			{
				firstNZPosInCG = n;
				break;
			}
		}

		for(n = firstNZPosInCG; n <=lastNZPosInCG; n++ )
		{
			absSum += dst[ scan[ n + subPos ]];
		}

		if(lastNZPosInCG>=0 && lastCG==-1) 
			lastCG = 1 ; 


		if( lastNZPosInCG-firstNZPosInCG>=SBH_THRESHOLD )	
		{
			uint signbit = (dst[scan[subPos+firstNZPosInCG]]>0?0:1) ;
			if( signbit!=(absSum&0x1) )  
			{
				int minCostInc = INT_MAX,  minPos =-1, finalChange=0, curCost=INT_MAX, curChange=0;

				for( n = (lastCG==1?lastNZPosInCG:SCAN_SET_SIZE-1) ; n >= 0; --n )
				{
					uint blkPos   = scan[ n+subPos ];
					if(dst[ blkPos ] != 0 )
					{
						if(deltaU[blkPos]>0)
						{
							curCost = - deltaU[blkPos]; 
							curChange=1 ;
						}
						else 
						{
							//curChange =-1;
							if(n==firstNZPosInCG && abs(dst[blkPos])==1)
								curCost=INT_MAX ; 
							else
							{
								curCost = deltaU[blkPos]; 
								curChange =-1;
							}
						}
					}
					else //if(dst[ blkPos ] == 0 )
					{
						if(n<firstNZPosInCG)
						{
							uint thisSignBit = (src[blkPos]>=0?0:1);
							if(thisSignBit != signbit )
								curCost = INT_MAX;
							else
							{ 
								curCost = - (deltaU[blkPos])  ;
								curChange = 1 ;
							}
						}
						else
						{
							curCost = - (deltaU[blkPos])  ;
							curChange = 1 ;
						}
					}

					if( curCost<minCostInc)
					{
						minCostInc = curCost ;
						finalChange = curChange ;
						minPos = blkPos ;
					}
				} //CG loop

				if(dst[minPos] == 32767 || dst[minPos] == -32768)
					finalChange = -1;

				if(src[minPos]>=0)
					dst[minPos] += finalChange ; 
				else 
					dst[minPos] -= finalChange ;
			} // Hide
		}
		if(lastCG==1) 
			lastCG=0 ;
	} // TU loop
}


void quant(henc_thread_t* et, int16_t* src, int16_t* dst, int scan_mode, int depth, int comp, int cu_mode, int is_intra, int *ac_sum, int cu_size, int per, int rem)
{
	int iLevel, auxLevel;
	int  iSign;
	picture_t *currpict = &et->enc_engine->current_pict;	
	slice_t *currslice = &currpict->slice;		
	int inv_depth = (et->max_cu_size_shift - (depth + (comp!=Y_COMP)));//et->max_cu_size_shift
	uint *scan = et->enc_engine->hvenc->scan_pyramid[scan_mode][inv_depth-1];
	int scan_type = (is_intra?0:3) + comp;
	int *quant_coeff = et->enc_engine->hvenc->quant_pyramid[inv_depth-2][scan_type][rem];
	uint transform_shift = MAX_TR_DYNAMIC_RANGE - et->bit_depth - inv_depth;
	int qbits = QUANT_SHIFT + per + transform_shift;                
	int qbits8 = qbits-8;
	int add = 171<<(qbits-9);//(currslice->slice_type==I_SLICE ? 171 : 85) << (qbits-9);//
	short *deltaU = et->aux_buff;
	int n;

	*ac_sum = 0;
	for( n = 0; n < cu_size*cu_size; n++ )
	{	
		auxLevel  = abs(src[n]);
		iSign   = SIGN(src[n]);      
#if ADAPTIVE_QP_SELECTION
		Int64 tmpLevel = (Int64)abs(iLevel) * quant_coeff[n];
		if( m_bUseAdaptQpSelect )
		{
			piArlCCoef[n] = (int)((tmpLevel + iAddC ) >> iQBitsC);
		}
		iLevel = (int)((tmpLevel + add ) >> qbits);
		deltaU[n] = (int)((tmpLevel - (iLevel<<qbits) )>> qbits8);
#else
		iLevel = (auxLevel * quant_coeff[n] + add ) >> qbits;

		deltaU[n] = (int)( (auxLevel * quant_coeff[n] - (iLevel<<qbits) )>> qbits8 );
#endif
		*ac_sum += iLevel;
		iLevel *= iSign;        
		dst[n] = clip(iLevel ,-32768, 32767); 
	} // for n

	if( et->pps->sign_data_hiding_flag)
	{
		if(*ac_sum>=2)
		{
			sign_bit_hidding( dst, src, scan, deltaU, cu_size, cu_size ) ;
		}
	}

}



void iquant(henc_thread_t* et, int16_t* src, int16_t* dst, int depth, int comp, int is_intra, int cu_size, int per, int rem)
{
	int iLevel;
	int inv_depth = (et->max_cu_size_shift - (depth+(comp!=Y_COMP)));
	int scan_type = (is_intra?0:3) + comp;
	int *dequant_coeff = et->enc_engine->hvenc->dequant_pyramid[inv_depth-2][scan_type][rem];
	uint transform_shift = MAX_TR_DYNAMIC_RANGE - et->bit_depth - inv_depth;
	int iq_shift = QUANT_IQUANT_SHIFT - QUANT_SHIFT - transform_shift + 4;
	int iq_add = (iq_shift>per)? 1 << (iq_shift - per - 1): 0;

	if(iq_shift>per)
	{
		int n;
		iq_shift=iq_shift-per;
		for( n = 0; n < cu_size*cu_size; n++ )
		{	
			iLevel  = src[n];
			iLevel  = (iLevel*dequant_coeff[n]+iq_add)>>iq_shift;
			dst[n] = clip(iLevel ,-32768, 32767); 
		} // for n
	}
	else 
	{
		int n;
		iq_shift=(per-iq_shift);
		iq_add = 0;
		for( n = 0; n < cu_size*cu_size; n++ )
		{	
			iLevel  = src[n];
			iLevel  = (iLevel*dequant_coeff[n])<<iq_shift;
			dst[n] = clip(iLevel ,-32768, 32767); 
		} // for n
	}

}
