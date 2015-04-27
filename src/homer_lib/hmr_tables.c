/*****************************************************************************
 * hmr_init_tables.c : homerHEVC encoding library
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
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *****************************************************************************/
/*
 * some of the work below is derived from HM HEVC reference code where 
 * the following license apply
 */
/****************************************************************************
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

#include <stdio.h>
#include	<malloc.h>
#include	<memory.h>
#include	<math.h>
#include	<limits.h>

#include	"hmr_os_primitives.h"
#include 	"homer_hevc_enc_api.h"
#include 	"hmr_common.h"
#include	"hmr_tables.h"
#include	"hmr_profiler.h"
#include	"hmr_sse42_functions.h"


static int num_scaling_list[NUM_SCALING_MODES]={6,6,6,2};
extern uint g_sigLastScanCG32x32[];


//scanning buffers: zigzag, horizontal, vertical y diagonal - zigzag is not used anymore
void init_scan_pyramid(hvenc_enc_t* ed, uint* pBuffZ, uint* pBuffH, uint* pBuffV, uint* pBuffD, int iWidth, int iHeight, int iDepth)
{
	const uint  uiNumScanPos  = (uint32_t) iWidth * iWidth;
	uint        uiNextScanPos = 0;

	if( iWidth < 16 )
	{
		uint uiScanLine;
		uint* pBuffTemp = pBuffD;
		if( iWidth == 8 )
		{
			pBuffTemp = g_sigLastScanCG32x32;
		}
		for( uiScanLine = 0; uiNextScanPos < uiNumScanPos; uiScanLine++ )
		{
			int    iPrimDim  = (int) uiScanLine;
			int    iScndDim  = 0;
			while( iPrimDim >= iWidth )
			{
				iScndDim++;
				iPrimDim--;
			}
			while( iPrimDim >= 0 && iScndDim < iWidth )
			{
				pBuffTemp[ uiNextScanPos ] = iPrimDim * iWidth + iScndDim ;
				uiNextScanPos++;
				iScndDim++;
				iPrimDim--;
			}
		}
	}
	if( iWidth > 4 )//diagonal
	{
		uint uiNumBlkSide = iWidth >> 2;
		uint uiNumBlks    = uiNumBlkSide * uiNumBlkSide;
		uint log2Blk      = iDepth-2;//g_aucConvertToBit[ uiNumBlkSide ] + 1;
		uint uiBlk;

		for( uiBlk = 0; uiBlk < uiNumBlks; uiBlk++ )
		{
			uint initBlkPos = ed->scan_pyramid[ DIAG_SCAN ][ log2Blk ][ uiBlk ];
			uint offsetY, offsetX, offsetD, offsetScan;
			uint uiScanLine;
			uiNextScanPos   = 0;

			if( iWidth == 32 )
			{
				initBlkPos = g_sigLastScanCG32x32[ uiBlk ];
			}
			offsetY    = initBlkPos / uiNumBlkSide;
			offsetX    = initBlkPos - offsetY * uiNumBlkSide;
			offsetD    = 4 * ( offsetX + offsetY * iWidth );
			offsetScan = 16 * uiBlk;
			for( uiScanLine = 0; uiNextScanPos < 16; uiScanLine++ )
			{
				int    iPrimDim  = (int)uiScanLine;
				int    iScndDim  = 0;
				while( iPrimDim >= 4 )
				{
					iScndDim++;
					iPrimDim--;
				}
				while( iPrimDim >= 0 && iScndDim < 4 )
				{
					pBuffD[ uiNextScanPos + offsetScan ] = iPrimDim * iWidth + iScndDim + offsetD;
					uiNextScanPos++;
					iScndDim++;
					iPrimDim--;
				}
			}
		}
	}

	if( iWidth > 2 )
	{
		uint uiCnt = 0;
		int blkY, blkX, x, y;
		int numBlkSide = iWidth >> 2;
		for(blkY=0; blkY < numBlkSide; blkY++)//horizontal
		{
			for(blkX=0; blkX < numBlkSide; blkX++)
			{
				uint offset    = blkY * 4 * iWidth + blkX * 4;
				for(y=0; y < 4; y++)
				{
					for(x=0; x < 4; x++)
					{
						pBuffH[uiCnt] = y*iWidth + x + offset;
						uiCnt ++;
					}
				}
			}
		}

		uiCnt = 0;
		for(blkX=0; blkX < numBlkSide; blkX++)//vertical
		{
			for(blkY=0; blkY < numBlkSide; blkY++)
			{
				uint offset    = blkY * 4 * iWidth + blkX * 4;
				for(x=0; x < 4; x++)
				{
					for(y=0; y < 4; y++)
					{
						pBuffV[uiCnt] = y*iWidth + x + offset;
						uiCnt ++;
					}
				}
			}
		}
	}
	else
	{
		uint uiCnt = 0;
		int iY, iX;
		for(iY=0; iY < iHeight; iY++)
		{
			for(iX=0; iX < iWidth; iX++)
			{
				pBuffH[uiCnt] = iY*iWidth + iX;
				uiCnt ++;
			}
		}

		//en este cambia la iteracion (0,2,1,3) para 2x2
		uiCnt = 0;
		for(iX=0; iX < iWidth; iX++)
		{
			for(iY=0; iY < iHeight; iY++)
			{
				pBuffV[uiCnt] = iY*iWidth + iX;
				uiCnt ++;
			}
		}    
	}
}

short* get_default_qtable(int size_mode, int list_index)
{
  short *src = 0;
  switch(size_mode)
  {
    case SCALING_MODE_4x4:
      src = (short*)QUANT_DEFAULT_4x4;
      break;
    case SCALING_MODE_8x8:
    case SCALING_MODE_16x16:
      src = (list_index<3) ? (short*)QUANT_DEFAULT_INTRA_8x8 : (short*)QUANT_DEFAULT_INTER_8x8;
      break;
    case SCALING_MODE_32x32:
      src = (list_index<1) ? (short*)QUANT_DEFAULT_INTRA_8x8 : (short*)QUANT_DEFAULT_INTER_8x8;
      break;
    default:
      src = NULL;
      break;
  }
  return src;
}

void init_quant_pyramids(hvenc_enc_t* ed, int* quant_pyramid, int* dequant_pyramid, double* scaling_error_pyramid, 
						 short* quant_def_table, int width, int height, int ratio, uint sizuNum, uint dc, int inv_depth, int qp)
{
	uint quant_scale[6] =		{26214,23302,20560,18396,16384,14564};    
	uint inv_quant_scale[6] =	{40,45,51,57,64,72};
	int i, j;
	int nsqth = (height < width) ? 4: 1; //height ratio for NSQT(NSQT= non square quant table)
	int nsqtw = (width < height) ? 4: 1; //width ratio for NSQT
	int quant_scales = quant_scale[qp]<<4;
	int iquant_scales = inv_quant_scale[qp];
	int iTransformShift = MAX_TR_DYNAMIC_RANGE - ed->bit_depth - inv_depth;  // Represents scaling through forward transform
	double dErrScale = ((double)(1<<SCALE_BITS))*pow(2.0,-2.0*iTransformShift); 


	for(j=0;j<height;j++)
	{
		for(i=0;i<width;i++)
		{
			quant_pyramid[j*width + i] = quant_scales / quant_def_table[sizuNum * (j * nsqth / ratio) + i * nsqtw /ratio];
			dequant_pyramid[j*width + i] = iquant_scales * quant_def_table[sizuNum * (j / ratio) + i / ratio];

			scaling_error_pyramid[j*width + i] = dErrScale / quant_pyramid[j*width + i] / quant_pyramid[j*width + i] / (1<<(2*(ed->bit_depth-8)));
		}
	}
	if(ratio > 1)
	{
		quant_pyramid[0] = quant_scales / dc;
		dequant_pyramid[0] = iquant_scales * dc;
	}
}


//setFlatScalingListen HM
void init_flat_quant_pyramids( hvenc_engine_t* ed, uint* quant_pyramid, uint* dequant_pyramid, double* scaling_error_pyramid, uint size, int inv_depth, int qp)
{
	uint quant_scale[6] =		{26214,23302,20560,18396,16384,14564};    
	uint inv_quant_scale[6] =	{40,45,51,57,64,72};
	uint i;
	uint quant = quant_scale[qp];
	uint inv_quant = inv_quant_scale[qp]<<4;
	int iTransformShift = MAX_TR_DYNAMIC_RANGE - ed->bit_depth - inv_depth;  // Represents scaling through forward transform
	double dErrScale = ((double)(1<<SCALE_BITS))*pow(2.0,-2.0*iTransformShift); 
	dErrScale = dErrScale / quant / quant / (1<<(2*(ed->bit_depth-8))); //(1<<DISTORTION_PRECISION_ADJUSTMENT(2*(bitDepth-8)));

	for(i=0;i<size;i++)
	{ 
		quant_pyramid[i] = quant;
		dequant_pyramid[i] = inv_quant;
		scaling_error_pyramid[i] = dErrScale;
	}
}	


void create_abs2raster_tables( unsigned short **zigzag, int total_depth, int depth, int start_value)
{
	int stride = 1 << ( total_depth - 1 );


	if ( depth == total_depth )
	{
		*zigzag[0] = start_value;
		(*zigzag)++;
	}
	else
	{
		int step = stride >> depth;
		create_abs2raster_tables( zigzag, total_depth, depth+1, start_value);
		create_abs2raster_tables( zigzag, total_depth, depth+1, start_value+step);
		create_abs2raster_tables( zigzag, total_depth, depth+1, start_value+step*stride);
		create_abs2raster_tables( zigzag, total_depth, depth+1, start_value+step*stride+step);
	}
}


void create_raster2abs_tables( unsigned short *zigzag, unsigned short *inv_zigzag, int max_cu_width, int max_cu_height, int total_depth)
{
	int i;

	int min_cu_width  = max_cu_width  >> ( total_depth - 1 );
	int min_cu_height = max_cu_height >> ( total_depth - 1 );
  
	int num_part_in_width  = max_cu_width  / min_cu_width;
	int num_part_in_height = max_cu_height / min_cu_height;

	for ( i = 0; i < num_part_in_width*num_part_in_height; i++ )
	{
		inv_zigzag[zigzag[i]] = i;
	}
}


void hmr_rd_init(hvenc_engine_t* ed, slice_t *currslice)
{
#define SHIFT_QP	12
	int		bitdepth_luma_qp_scale = 0;
	double	qp_factor = 0.4624;
	double	qp_temp = (double) ed->current_pict.slice.qp /* pict_qp */+ bitdepth_luma_qp_scale - SHIFT_QP;//
	double	lambda_scale = 1.0 - clip(0.05*(double)(/*ed->mb_interlaced*/0 ? (ed->gop_size-1)/2 : (ed->gop_size-1)), 0.0, 0.5);
	double	lambda;
    int depth, poc = currslice->poc%ed->gop_size;

	if (poc == 0)
		depth = 0;
	else
	{
		int i, j;
		int step = ed->gop_size;
		depth = 0;
		for(i=step>>1; i>=1; i>>=1)
		{
			for (j=i; j<ed->gop_size; j+=step )
			{
				if (j == poc)
				{
					i=0;
					break;
				}
			}
			step >>= 1;
			depth++;
		}
	}

//	if(currslice->slice_type==I_SLICE)
	{
		qp_factor=0.57*lambda_scale;
	}

	lambda = qp_factor*pow( 2.0, qp_temp/3.0 );

    if ( depth>0 )
    {
        lambda *= clip((qp_temp / 6.0), 2.00, 4.00); // (j == B_SLICE && p_cur_frm->layer != 0 )
    }

	ed->rd.lambda = lambda;
	ed->rd.sqrt_lambda = sqrt(lambda);
	ed->rd.lambda_SAD = (uint32_t)floor(65536.0 * ed->rd.sqrt_lambda);
	ed->rd.lambda_SSE = (uint32_t)floor(65536.0 * ed->rd.sqrt_lambda);
}

int find_scan_mode(int is_intra, int is_luma, int width, int dir_mode, int up_left_luma_dir_mode)//up_left_luma_dir_mode solo vale para la chroma
{
	int scan_idx  = DIAG_SCAN;
	int ctx_idx;

	if (!is_intra) 
	{
		scan_idx = DIAG_SCAN;
		return scan_idx;
	}

	switch(width)
	{
		case  2: ctx_idx = 6; break;
		case  4: ctx_idx = 5; break;
		case  8: ctx_idx = 4; break;
		case 16: ctx_idx = 3; break;
		case 32: ctx_idx = 2; break;
		case 64: ctx_idx = 1; break;
		default: ctx_idx = 0; break;
	}

	if(is_luma)
	{
		if (ctx_idx >3 && ctx_idx < 6) 
			scan_idx = abs(dir_mode - VER_IDX) < 5 ? HOR_SCAN : (abs(dir_mode - HOR_IDX) < 5 ? VER_SCAN : DIAG_SCAN);
	}
	else
	{
		if( dir_mode == DM_CHROMA_IDX )	
		{
			dir_mode = up_left_luma_dir_mode;
		}
		if (ctx_idx >4 && ctx_idx < 7) //if multiple scans supported for transform size
			scan_idx = abs(dir_mode - VER_IDX) < 5 ? HOR_SCAN : (abs(dir_mode - HOR_IDX) < 5 ? VER_SCAN : DIAG_SCAN);
	}
	return scan_idx;
}