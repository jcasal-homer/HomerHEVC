/*****************************************************************************
 * hmr_deblocking_filter.c : homerHEVC encoding library
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *****************************************************************************/

#include <math.h>
#include <memory.h>
#include "hmr_common.h"
#include "hmr_private.h"

const uint8_t sm_tcTable[54] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,6,6,7,8,9,10,11,13,14,16,18,20,22,24
};

const uint8_t sm_betaTable[52] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64
};

/*uint xCalcBsIdx( TComDataCU* pcCU, UInt uiAbsZorderIdx, Int dir, Int iEdgeIdx, Int iBaseUnitIdx )
  {
    TComPic* const pcPic = pcCU->getPic();
    const UInt uiLCUWidthInBaseUnits = pcPic->getNumPartInWidth();
    if( dir == 0 )
    {
      return g_auiRasterToZscan[g_auiZscanToRaster[uiAbsZorderIdx] + iBaseUnitIdx * uiLCUWidthInBaseUnits + iEdgeIdx ];
    }
    else
    {
      return g_auiRasterToZscan[g_auiZscanToRaster[uiAbsZorderIdx] + iEdgeIdx * uiLCUWidthInBaseUnits + iBaseUnitIdx ];
    }
  } 
*/ 

__inline uint bx_index(henc_thread_t* et, cu_partition_info_t* curr_cu_info, int deblock_dir, int cu_width_units, int deblock_edge_idx, int index)
{
	if(deblock_dir==EDGE_VER)
		return et->enc_engine->raster2abs_table[et->enc_engine->abs2raster_table[curr_cu_info->abs_index]+index*cu_width_units+deblock_edge_idx];
	else
		return et->enc_engine->raster2abs_table[et->enc_engine->abs2raster_table[curr_cu_info->abs_index]+deblock_edge_idx*cu_width_units+index];

}


void set_edge_filter_0(henc_thread_t* et, cu_partition_info_t* curr_cu_info, uint8_t *edge_buff , uint8_t *bs_buff , int depht, int abs_index, int cu_width, int cu_height, int deblock_dir, int deblock_internal_edge, int num_elements, int max_cu_width_units)
{
	int i;
	int aux0;

	if(deblock_dir==EDGE_VER)
		aux0 = max_cu_width_units;
	else
		aux0 = 1;
	for( i = 0; i < num_elements; i++ )
	{
//			uint bs_idx = enc_engine->raster2abs_table[enc_engine->abs2raster_table[curr_cu_info->abs_index]+i*max_cu_width_units+deblock_edge_idx];
		uint bs_idx = et->enc_engine->raster2abs_table[et->enc_engine->abs2raster_table[curr_cu_info->abs_index]+i*aux0];//+aux1];
		edge_buff[bs_idx] = deblock_internal_edge;
		bs_buff[bs_idx] = deblock_internal_edge;
	}
}

void set_edge_filter_0_subdiv(henc_thread_t* et, cu_partition_info_t*	curr_cu_info, uint8_t *edge_buff , uint8_t *bs_buff , int depht, int abs_index, int cu_width, int cu_height, int deblock_dir, int deblock_internal_edge, int num_elements, int max_cu_width_units, int edge_idx)
{
	if(deblock_dir==EDGE_VER)
	{
		int i;
		for( i = 0; i < num_elements; i++ )
		{
			uint bs_idx = et->enc_engine->raster2abs_table[et->enc_engine->abs2raster_table[curr_cu_info->abs_index]+i*max_cu_width_units+edge_idx];
			edge_buff[bs_idx] = deblock_internal_edge;
			if(edge_idx==0)
				bs_buff[bs_idx] = deblock_internal_edge;
		}
	}
	else
	{
		int i;		
		for( i = 0; i < num_elements; i++ )
		{
			uint bs_idx = et->enc_engine->raster2abs_table[et->enc_engine->abs2raster_table[curr_cu_info->abs_index]+edge_idx*max_cu_width_units+i];
			edge_buff[bs_idx] = deblock_internal_edge;
			if(edge_idx==0)
				bs_buff[bs_idx] = deblock_internal_edge;
		}
	}

}



void set_edge_filter(hvenc_engine_t* enc_engine, cu_partition_info_t*	curr_cu_info, uint8_t *edge_buff , uint8_t *bs_buff , int depht, int abs_index, int cu_width, int cu_height, int deblock_dir, int deblock_internal_edge, int deblock_edge_idx, int num_elements, int max_cu_width_units)
{
	int i;
//	int num_elements = (max(cu_width, cu_height)>>2);
//	int max_cu_width_units = (enc_engine->max_cu_size>>2);
	int aux0, aux1;

	if(deblock_dir==EDGE_VER)
	{
		aux0 = max_cu_width_units;
		aux1 = deblock_edge_idx;
	}
	else
	{
		aux0 = 1;
		aux1 = deblock_edge_idx*max_cu_width_units;
	}

	for( i = 0; i < num_elements; i++ )
	{
//			uint bs_idx = enc_engine->raster2abs_table[enc_engine->abs2raster_table[curr_cu_info->abs_index]+i*max_cu_width_units+deblock_edge_idx];
		uint bs_idx = enc_engine->raster2abs_table[enc_engine->abs2raster_table[curr_cu_info->abs_index]+i*aux0+aux1];
		edge_buff[bs_idx] = deblock_internal_edge;
		if (deblock_edge_idx == 0)
		{
			bs_buff[bs_idx] = deblock_internal_edge;
		}
	}
}

//xGetBoundaryStrengthSingle
void get_boundary_strength_single(henc_thread_t* et, slice_t *currslice, ctu_info_t *ctu, cu_partition_info_t *curr_cu_info, cu_partition_info_t *sub_cu_info, int dir)
{
	int bs;
	ctu_info_t	*ctu_aux;
	uint		abs_idx, aux_abs_idx = 0;

	abs_idx = sub_cu_info->abs_index;
	//-- Calculate Block Index
	if (dir == EDGE_VER)
	{
		ctu_aux = get_pu_left(ctu, sub_cu_info, &aux_abs_idx);//ctu->ctu_left;
	}
	else  // (dir == EDGE_HOR)
	{
		ctu_aux = get_pu_top(ctu, sub_cu_info, &aux_abs_idx, FALSE);//ctu->ctu_top;//getPUAbove( aux_abs_idx, m_uiAbsIdxInLCU + abs_index, TRUE, TRUE );
	}

	//-- Set BS for Intra MB : BS = 4 or 3
	if(ctu_aux->pred_mode[aux_abs_idx] == INTRA_MODE || ctu->pred_mode[abs_idx] == INTRA_MODE)//if ( pcCUP->isIntra(uiPartP) || pcCUQ->isIntra(uiPartQ) )
	{
		bs = 2;
	}

	//-- Set BS for not Intra MB : BS = 2 or 1 or 0
	if(ctu_aux->pred_mode[aux_abs_idx] != INTRA_MODE && ctu->pred_mode[sub_cu_info->abs_index] != INTRA_MODE)//if ( !pcCUP->isIntra(uiPartP) && !pcCUQ->isIntra(uiPartQ) )
	{
		uint ns_part_curr = sub_cu_info->abs_index;//nsPartQ
		uint ns_part_aux = aux_abs_idx;//nsPartP

//		if ( m_aapucBS[dir][uiAbsPartIdx] && (pcCUQ->getCbf( nsPartQ, TEXT_LUMA, pcCUQ->getTransformIdx(nsPartQ)) != 0 || pcCUP->getCbf( nsPartP, TEXT_LUMA, pcCUP->getTransformIdx(nsPartP) ) != 0) )
		if(et->deblock_filter_strength_bs[dir][abs_idx] && ((CBF(ctu, ns_part_curr, Y_COMP, ctu->tr_idx[ns_part_curr]) != 0) || (CBF(ctu_aux, ns_part_aux, Y_COMP, ctu_aux->tr_idx[ns_part_aux]) != 0)))
		{
			bs = 1;
		}
		else
		{
			if (dir == EDGE_HOR)
			{
//				pcCUP = pcCUQ->getPUAbove(uiPartP, uiPartQ, !pcCU->getSlice()->getLFCrossSliceBoundaryFlag(), false, !m_bLFCrossTileBoundary);
				ctu_aux = get_pu_top(ctu, sub_cu_info, &aux_abs_idx, FALSE);
			}
			if (currslice->slice_type == B_SLICE)//(pcSlice->isInterB() || pcCUP->getSlice()->isInterB())
			{
/*				Int iRefIdx;
				TComPic *piRefP0, *piRefP1, *piRefQ0, *piRefQ1;
				iRefIdx = pcCUP->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartP);
				piRefP0 = (iRefIdx < 0) ? NULL : pcCUP->getSlice()->getRefPic(REF_PIC_LIST_0, iRefIdx);
				iRefIdx = pcCUP->getCUMvField(REF_PIC_LIST_1)->getRefIdx(uiPartP);
				piRefP1 = (iRefIdx < 0) ? NULL : pcCUP->getSlice()->getRefPic(REF_PIC_LIST_1, iRefIdx);
				iRefIdx = pcCUQ->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartQ);
				piRefQ0 = (iRefIdx < 0) ? NULL : pcSlice->getRefPic(REF_PIC_LIST_0, iRefIdx);
				iRefIdx = pcCUQ->getCUMvField(REF_PIC_LIST_1)->getRefIdx(uiPartQ);
				piRefQ1 = (iRefIdx < 0) ? NULL : pcSlice->getRefPic(REF_PIC_LIST_1, iRefIdx);

				TComMv pcMvP0 = pcCUP->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartP);
				TComMv pcMvP1 = pcCUP->getCUMvField(REF_PIC_LIST_1)->getMv(uiPartP);
				TComMv pcMvQ0 = pcCUQ->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartQ);
				TComMv pcMvQ1 = pcCUQ->getCUMvField(REF_PIC_LIST_1)->getMv(uiPartQ);

				if (piRefP0 == NULL) pcMvP0.setZero();
				if (piRefP1 == NULL) pcMvP1.setZero();
				if (piRefQ0 == NULL) pcMvQ0.setZero();
				if (piRefQ1 == NULL) pcMvQ1.setZero();

				if ( ((piRefP0==piRefQ0)&&(piRefP1==piRefQ1)) || ((piRefP0==piRefQ1)&&(piRefP1==piRefQ0)) )
				{
					if ( piRefP0 != piRefP1 )   // Different L0 & L1
					{
						if ( piRefP0 == piRefQ0 )
						{
							uiBs  = ((abs(pcMvQ0.getHor() - pcMvP0.getHor()) >= 4) ||
								(abs(pcMvQ0.getVer() - pcMvP0.getVer()) >= 4) ||
								(abs(pcMvQ1.getHor() - pcMvP1.getHor()) >= 4) ||
								(abs(pcMvQ1.getVer() - pcMvP1.getVer()) >= 4)) ? 1 : 0;
						}
						else
						{
							uiBs  = ((abs(pcMvQ1.getHor() - pcMvP0.getHor()) >= 4) ||
								(abs(pcMvQ1.getVer() - pcMvP0.getVer()) >= 4) ||
								(abs(pcMvQ0.getHor() - pcMvP1.getHor()) >= 4) ||
								(abs(pcMvQ0.getVer() - pcMvP1.getVer()) >= 4)) ? 1 : 0;
						}
					}
					else    // Same L0 & L1
					{
						uiBs  = ((abs(pcMvQ0.getHor() - pcMvP0.getHor()) >= 4) ||
							(abs(pcMvQ0.getVer() - pcMvP0.getVer()) >= 4) ||
							(abs(pcMvQ1.getHor() - pcMvP1.getHor()) >= 4) ||
							(abs(pcMvQ1.getVer() - pcMvP1.getVer()) >= 4)) &&
							((abs(pcMvQ1.getHor() - pcMvP0.getHor()) >= 4) ||
							(abs(pcMvQ1.getVer() - pcMvP0.getVer()) >= 4) ||
							(abs(pcMvQ0.getHor() - pcMvP1.getHor()) >= 4) ||
							(abs(pcMvQ0.getVer() - pcMvP1.getVer()) >= 4)) ? 1 : 0;
					}
				}
				else // for all different Ref_Idx
				{
					uiBs = 1;
				}
*/			}
			else  // pcSlice->isInterP()
			{
				int ref_idx = ctu_aux->mv_ref_idx[REF_PIC_LIST_0][aux_abs_idx];
				video_frame_t *ref_frame_aux, *ref_frame_curr;
				motion_vector_t mv_aux, mv_curr;
				ref_frame_aux  = (ref_idx<0)?NULL:(currslice->ref_pic_list[REF_PIC_LIST_0][ref_idx]);
				ref_idx = ctu->mv_ref_idx[REF_PIC_LIST_0][sub_cu_info->abs_index];
				ref_frame_curr  = (ref_idx<0)?NULL:(currslice->ref_pic_list[REF_PIC_LIST_0][ref_idx]);
				mv_aux = ctu_aux->mv_ref[REF_PIC_LIST_0][aux_abs_idx];
				mv_curr = ctu->mv_ref[REF_PIC_LIST_0][sub_cu_info->abs_index];

				if (ref_frame_aux == NULL) {mv_aux.hor_vector = mv_aux.ver_vector = 0;}
				if (ref_frame_curr == NULL) {mv_curr.hor_vector = mv_curr.ver_vector = 0;}

				bs  = ((ref_frame_aux != ref_frame_curr) ||
					(abs(mv_curr.hor_vector - mv_aux.hor_vector) >= 4) ||
					(abs(mv_curr.ver_vector - mv_aux.ver_vector) >= 4)) ? 1 : 0;
			}
		}   // enf of "if( one of BCBP == 0 )"
	}   // enf of "if( not Intra )"

	et->deblock_filter_strength_bs[dir][abs_idx] = bs;
//	m_aapucBS[dir][uiAbsPartIdx] = uiBs;
}

__inline int calc_dp( int16_t* src, int offset)
{
  return abs( src[-offset*3] - 2*src[-offset*2] + src[-offset] ) ;
}
  
__inline int calc_dq( int16_t* src, int offset)
{
  return abs( src[0] - 2*src[offset] + src[offset*2] );
}


__inline int use_strong_filter( int offset, int d, int beta, int tc, int16_t* src)
{
  int16_t m4  = src[0];
  int16_t m3  = src[-offset];
  int16_t m7  = src[ offset*3];
  int16_t m0  = src[-offset*4];

  int d_strong = abs(m0-m3) + abs(m7-m4);

  return ( (d_strong < (beta>>3)) && (d<(beta>>2)) && ( abs(m3-m4) < ((tc*5+1)>>1)) );
}

__inline void filter_luma( int16_t* src, int offset, int tc , int sw, int bPartPNoFilter, int bPartQNoFilter, int iThrCut, int bFilterSecondP, int bFilterSecondQ)
{
	int delta;

	int16_t m4  = src[0];
	int16_t m3  = src[-offset];
	int16_t m5  = src[ offset];
	int16_t m2  = src[-offset*2];
	int16_t m6  = src[ offset*2];
	int16_t m1  = src[-offset*3];
	int16_t m7  = src[ offset*3];
	int16_t m0  = src[-offset*4];

	if (sw)
	{
		src[-offset]   = clip(((m1 + 2*m2 + 2*m3 + 2*m4 + m5 + 4) >> 3), m3-2*tc, m3+2*tc);
		src[0]          = clip(((m2 + 2*m3 + 2*m4 + 2*m5 + m6 + 4) >> 3), m4-2*tc, m4+2*tc);
		src[-offset*2] = clip(((m1 + m2 + m3 + m4 + 2)>>2), m2-2*tc, m2+2*tc);
		src[ offset]   = clip(((m3 + m4 + m5 + m6 + 2)>>2), m5-2*tc, m5+2*tc);
		src[-offset*3] = clip(((2*m0 + 3*m1 + m2 + m3 + m4 + 4 )>>3), m1-2*tc, m1+2*tc);
		src[ offset*2] = clip(((m3 + m4 + m5 + 3*m6 + 2*m7 +4 )>>3), m6-2*tc, m6+2*tc);
	}
	else
	{
		/* Weak filter */
		delta = (9*(m4-m3) -3*(m5-m2) + 8)>>4 ;

		if ( abs(delta) < iThrCut )
		{
			int tc2;
			delta = clip(delta,-tc, tc);        
			src[-offset] = clip((m3+delta), 0, 255);
			src[0] = clip((m4-delta), 0, 255);

			tc2 = tc>>1;
			if(bFilterSecondP)
			{
				int delta1 = clip((( ((m1+m3+1)>>1)- m2+delta)>>1), -tc2, tc2);
				src[-offset*2] = clip((m2+delta1), 0, 255);
			}
			if(bFilterSecondQ)
			{
				int delta2 = clip((( ((m6+m4+1)>>1)- m5-delta)>>1), -tc2, tc2);
				src[ offset] = clip((m5+delta2), 0, 255);
			}
		}
	}

	if(bPartPNoFilter)
	{
		src[-offset] = m3;
		src[-offset*2] = m2;
		src[-offset*3] = m1;
	}
	if(bPartQNoFilter)
	{
		src[0] = m4;
		src[ offset] = m5;
		src[ offset*2] = m6;
	}
}


//Void xEdgeFilterLuma( TComDataCU* pcCU, UInt uiAbsZorderIdx, UInt uiDepth, Int iDir, Int iEdge  )
void deblock_filter_luma(henc_thread_t* et, ctu_info_t *ctu, cu_partition_info_t *curr_cu_info, wnd_t *decoded_wnd, slice_t *currslice, int dir, int edge)
{
	int offset, src_increment;
	int decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, Y_COMP);
	int curr_part_global_x = ctu->x[Y_COMP]+curr_cu_info->x_position;
	int curr_part_global_y = ctu->y[Y_COMP]+curr_cu_info->y_position;
	int partition_p_no_filter = 0;
	int partition_q_no_filter = 0; 
	int num_pels_in_partition = (et->max_cu_size>>et->max_cu_depth);
	int16_t *decoded_buff  = WND_POSITION_2D(int16_t *, *decoded_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, 0, et->ctu_width);
	int16_t *aux_decoded_buff = decoded_buff;
//	int	pixels_per_partition = enc_engine->max_cu_size >> enc_engine->max_cu_depth;
	int	idx;
	int num_partitions;
	int max_cu_width_units = (et->max_cu_size>>2);		

	if (dir == EDGE_VER)
	{
		offset = 1;
		src_increment = decoded_buff_stride;
		aux_decoded_buff += edge*num_pels_in_partition;
	}
	else  // (iDir == EDGE_HOR)
	{
		offset = decoded_buff_stride;
		src_increment = 1;
		aux_decoded_buff += edge*num_pels_in_partition*decoded_buff_stride;
	}

//	num_partitions = enc_engine->max_cu_size >> curr_cu_info->depth;
//	num_partitions /= num_pels_in_partition;

	num_partitions = curr_cu_info->size/num_pels_in_partition;

	for ( idx = 0; idx < num_partitions; idx++ )
	{
		int index_tc, index_b;
		int bs_abs_idx = bx_index(et, curr_cu_info, dir, max_cu_width_units, edge, idx);
		int bs = et->deblock_filter_strength_bs[dir][bs_abs_idx];

		if ( bs )
		{
			int block_idx;
			int qp;
			int partition_q_idx = bs_abs_idx;
			ctu_info_t	*ctu_aux;
			uint /*abs_idx, */aux_abs_idx = 0;
			int bit_depth_scale;
			int  beta_offset_div2 = currslice->slice_beta_offset_div2;
			int  tc_offset_div2 = currslice->slice_tc_offset_div2;
			int tc, beta, side_threshold, threshold_cut;
			int blocks_in_partition;// = max(num_pels_in_partition>>2, 1);
			int pcm_filter = (currslice->sps->pcm_enabled_flag && currslice->sps->pcm_loop_filter_disable_flag);//(pcCU->getSlice()->getSPS()->getUsePCM() && pcCU->getSlice()->getSPS()->getPCMFilterDisableFlag())? true : false;
			cu_partition_info_t *sub_cu_info = &ctu->partition_list[et->enc_engine->partition_depth_start[et->enc_engine->max_cu_depth]];//4x4 blocks

			sub_cu_info += bs_abs_idx;

			// Derive neighboring PU index
			if (dir == EDGE_VER)
			{
				ctu_aux = get_pu_left(ctu, sub_cu_info, &aux_abs_idx);
			}
			else  // (iDir == EDGE_HOR)
			{
				ctu_aux = get_pu_top(ctu, sub_cu_info, &aux_abs_idx, FALSE);
			}

			//qp = (ctu_aux->qp+ctu->qp+1)>>1;
			//if(ctu_aux->qp[aux_abs_idx]!=32)
			qp = (ctu_aux->qp[aux_abs_idx] + ctu->qp[curr_cu_info->abs_index]+1)>>1;
			//qp = (sub_cu_info->qp + sub_cu_info->qp)>>1;

			bit_depth_scale = 1<<(et->bit_depth-8);			

			index_tc =  clip((int)(qp + DEFAULT_INTRA_TC_OFFSET*(bs-1) + (tc_offset_div2 << 1)), 0, MAX_QP+DEFAULT_INTRA_TC_OFFSET);
			index_b = clip(qp + (beta_offset_div2 << 1), 0, MAX_QP);
			tc = sm_tcTable[index_tc];
			beta = sm_betaTable[index_b];
			side_threshold = (beta+(beta>>1))>>3;
			threshold_cut = tc*10;

			blocks_in_partition = max(num_pels_in_partition>>2, 1);
//			UInt  uiBlocksInPart = uiPelsInPart / 4 ? uiPelsInPart / 4 : 1;
			for (block_idx = 0; block_idx<blocks_in_partition; block_idx++)//This is always [0,1]
			{
				int dp0 = calc_dp( aux_decoded_buff+src_increment*(idx*num_pels_in_partition+block_idx*4+0), offset);
				int dq0 = calc_dq( aux_decoded_buff+src_increment*(idx*num_pels_in_partition+block_idx*4+0), offset);
				int dp3 = calc_dp( aux_decoded_buff+src_increment*(idx*num_pels_in_partition+block_idx*4+3), offset);
				int dq3 = calc_dq( aux_decoded_buff+src_increment*(idx*num_pels_in_partition+block_idx*4+3), offset);
				int d0 = dp0 + dq0;
				int d3 = dp3 + dq3;
				int d =  d0 + d3;
				int dp = dp0 + dp3;
				int dq = dq0 + dq3;

				if (pcm_filter || currslice->pps->transquant_bypass_enable_flag)
				{
					// Check if each of PUs is I_PCM with LF disabling
//					partition_p_no_filter = (pcm_filter && ctu_aux->getIPCMFlag(uiPartPIdx));
//					partition_q_no_filter = (pcm_filter && ctu->getIPCMFlag(uiPartQIdx));

					// check if each of PUs is lossless coded
//					partition_p_no_filter = partition_p_no_filter || (ctu_aux->isLosslessCoded(uiPartPIdx) );
//					partition_q_no_filter = partition_q_no_filter || (ctu->isLosslessCoded(uiPartQIdx) );
				}

				if (d < beta)
				{ 
					int i;
					int filter_p = (dp < side_threshold);
					int filter_q = (dq < side_threshold);

					int sw = use_strong_filter( offset, 2*d0, beta, tc, aux_decoded_buff+src_increment*(idx*num_pels_in_partition+block_idx*4+0)) && 
						use_strong_filter( offset, 2*d3, beta, tc, aux_decoded_buff+src_increment*(idx*num_pels_in_partition+block_idx*4+3));

					for ( i = 0; i < DEBLOCK_SMALLEST_BLOCK/2; i++)
					{
						filter_luma(aux_decoded_buff+src_increment*(idx*num_pels_in_partition+block_idx*4+i), offset, tc, sw, partition_p_no_filter, partition_q_no_filter, threshold_cut, filter_p, filter_q);
					}
				}
			}
		}
	}
}


//__inline Void TComLoopFilter::xPelFilterChroma( Pel* piSrc, Int iOffset, Int tc, Bool bPartPNoFilter, Bool bPartQNoFilter)
__inline void filter_chroma( int16_t* src, int offset, int tc , int bPartPNoFilter, int bPartQNoFilter)
{
  int delta;
  
  int16_t m4  = src[0];
  int16_t m3  = src[-offset];
  int16_t m5  = src[ offset];
  int16_t m2  = src[-offset*2];
  
  delta = clip((((( m4 - m3 ) << 2 ) + m2 - m5 + 4 ) >> 3), -tc, tc);
  src[-offset] = clip((m3+delta),0,255);
  src[0] = clip((m4-delta),0,255);

  if(bPartPNoFilter)
  {
    src[-offset] = m3;
  }
  if(bPartQNoFilter)
  {
    src[0] = m4;
  }
}


extern const uint8_t chroma_scale_conversion_table[];

void deblock_filter_chroma(henc_thread_t* et, ctu_info_t *ctu, cu_partition_info_t *curr_cu_info, wnd_t *decoded_wnd, slice_t *currslice, int dir, int edge)
{
	int offset, src_increment;
	int decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, CHR_COMP);
	int curr_part_global_x = ctu->x[CHR_COMP]+curr_cu_info->x_position_chroma;
	int curr_part_global_y = ctu->y[CHR_COMP]+curr_cu_info->y_position_chroma;
	int partition_p_no_filter = 0;
	int partition_q_no_filter = 0; 
	int num_pels_in_partition = (et->max_cu_size>>(et->max_cu_depth+1));
//	uint8_t *decoded_buff[3]  = {0,WND_POSITION_2D(uint8_t *, *decoded_wnd, U_COMP, curr_part_global_x, curr_part_global_y, 0, enc_engine->ctu_width), WND_POSITION_2D(uint8_t *, *decoded_wnd, V_COMP, curr_part_global_x, curr_part_global_y, 0, enc_engine->ctu_width)};
	int16_t *aux_decoded_buff;
//	int	pixels_per_partition = enc_engine->max_cu_size >> enc_engine->max_cu_depth;
	int	idx;
	int num_partitions;
	int max_cu_width_units = (et->max_cu_size>>2);		
	int ch_component;
	int raster_index = et->enc_engine->abs2raster_table[curr_cu_info->abs_index];
	int uiEdgeNumInLCUVert = raster_index % max_cu_width_units + edge;
	int uiEdgeNumInLCUHor = raster_index / max_cu_width_units + edge;


	num_partitions = curr_cu_info->size_chroma/num_pels_in_partition;//enc_engine->max_cu_size >> curr_cu_info->depth;

	if ( (num_pels_in_partition < DEBLOCK_SMALLEST_BLOCK) && (( (uiEdgeNumInLCUVert%(DEBLOCK_SMALLEST_BLOCK/num_pels_in_partition))&&(dir==EDGE_VER) ) || ( (uiEdgeNumInLCUHor%(DEBLOCK_SMALLEST_BLOCK/num_pels_in_partition))&& dir==EDGE_HOR ) ))
	{
		return;
	}


	if (dir == EDGE_VER)
	{
		offset = 1;
		src_increment = decoded_buff_stride;
//		aux_decoded_buff += edge*num_pels_in_partition;
	}
	else  // (iDir == EDGE_HOR)
	{
		offset = decoded_buff_stride;
		src_increment = 1;
//		aux_decoded_buff += edge*num_pels_in_partition*decoded_buff_stride;
	}

	for ( idx = 0; idx < num_partitions; idx++ )
	{
		int index_tc;
		int bs_abs_idx = bx_index(et, curr_cu_info, dir, max_cu_width_units, edge, idx);
		int bs = et->deblock_filter_strength_bs[dir][bs_abs_idx];

		if( bs > 1)
		{
			int block_idx;
			int qp;
			int partition_q_idx = bs_abs_idx;
			ctu_info_t	*ctu_aux;
			uint /*abs_idx, */aux_abs_idx = 0;
			int bit_depth_scale;
			int  beta_offset_div2 = currslice->slice_beta_offset_div2;
			int  tc_offset_div2 = currslice->slice_tc_offset_div2;
			int tc;//, beta, side_threshold, threshold_cut;
//			int blocks_in_partition;// = max(num_pels_in_partition>>2, 1);
			int pcm_filter = (currslice->sps->pcm_enabled_flag && currslice->sps->pcm_loop_filter_disable_flag);//(pcCU->getSlice()->getSPS()->getUsePCM() && pcCU->getSlice()->getSPS()->getPCMFilterDisableFlag())? true : false;
			cu_partition_info_t *sub_cu_info = &ctu->partition_list[et->enc_engine->partition_depth_start[et->max_cu_depth]];//4x4 blocks

			sub_cu_info += bs_abs_idx;

			// Derive neighboring PU index
			if (dir == EDGE_VER)
			{
				ctu_aux = get_pu_left(ctu, sub_cu_info, &aux_abs_idx);
			}
			else  // (iDir == EDGE_HOR)
			{
				ctu_aux = get_pu_top(ctu, sub_cu_info, &aux_abs_idx, FALSE);
			}

/*				if (bPCMFilter || pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
			{
				// Check if each of PUs is I_PCM with LF disabling
				bPartPNoFilter = (bPCMFilter && pcCUP->getIPCMFlag(uiPartPIdx));
				bPartQNoFilter = (bPCMFilter && pcCUQ->getIPCMFlag(uiPartQIdx));

				// check if each of PUs is lossless coded
				bPartPNoFilter = bPartPNoFilter || (pcCUP->isLosslessCoded(uiPartPIdx));
				bPartQNoFilter = bPartQNoFilter || (pcCUQ->isLosslessCoded(uiPartQIdx));
			}
*/

			for(ch_component = U_COMP;ch_component<=V_COMP;ch_component++)
			{
				int chr_qp_offset = (ch_component == U_COMP)?currslice->pps->cb_qp_offset:currslice->pps->cr_qp_offset;
				int chr_qp;
				aux_decoded_buff = WND_POSITION_2D(int16_t *, *decoded_wnd, ch_component, curr_part_global_x, curr_part_global_y, 0, et->ctu_width);
				aux_decoded_buff += edge*num_pels_in_partition*offset;

//				qp = ((ctu_aux->qp+ctu->qp+1)>>1)+chr_qp_offset;
//				qp = (ctu_aux->qp[aux_abs_idx] + sub_cu_info->qp)>>1;
//				qp = (sub_cu_info->qp + curr_cu_info->qp)>>1+chr_qp_offset;
				qp = ((ctu_aux->qp[aux_abs_idx] + ctu->qp[curr_cu_info->abs_index]+1)>>1)+chr_qp_offset;

				chr_qp = chroma_scale_conversion_table[clip(qp,0,57)];

				bit_depth_scale = 1<<(et->bit_depth-8);			

				index_tc =  clip((int)(chr_qp + DEFAULT_INTRA_TC_OFFSET*(bs-1) + (tc_offset_div2 << 1)), 0, MAX_QP+DEFAULT_INTRA_TC_OFFSET);
//				index_b = clip(qp + (beta_offset_div2 << 1), MIN_QP, MAX_QP);
				tc = sm_tcTable[index_tc];
//				beta = sm_betaTable[index_b];
//				side_threshold = (beta+(beta>>1))>>3;
//				threshold_cut = tc*10;

//				blocks_in_partition = max(num_pels_in_partition>>2, 1);
	//			UInt  uiBlocksInPart = uiPelsInPart / 4 ? uiPelsInPart / 4 : 1;
				for (block_idx = 0; block_idx<num_pels_in_partition; block_idx++)//This is always [0,1]
				{
					filter_chroma(aux_decoded_buff+src_increment*(block_idx+idx*num_pels_in_partition), offset, tc , partition_p_no_filter, partition_q_no_filter);
//					filter_luma(aux_decoded_buff+src_increment*(idx*num_pels_in_partition+block_idx*4+i), offset, tc, sw, partition_p_no_filter, partition_q_no_filter, threshold_cut, filter_p, filter_q);
				}
			}//for (chr_component..
		}
	}//while

}



void deblock_cu(henc_thread_t* et, slice_t *currslice, ctu_info_t* ctu, cu_partition_info_t *curr_cu_info, int dir)
{
	int part_idx;
	int abs_x = ctu->x[Y_COMP]+curr_cu_info->x_position;
	int abs_y = ctu->y[Y_COMP]+curr_cu_info->y_position;
//	int left_edge = (!(abs_x==0 || deblocking_filter_disable));//left_edge indica si hay que filtrar la columna izq del PU o no
//	int top_edge = (!(abs_y==0 || deblocking_filter_disable));//top_edge indica si hay que filtrar la fila top del PU o no
	int num_pels_in_partition = (et->max_cu_size>>et->max_cu_depth);
	int deblock_partition_idx_incr = max((DEBLOCK_SMALLEST_BLOCK / num_pels_in_partition), 1);
	int pu_size_in_deblock_units = curr_cu_info->size/num_pels_in_partition;
					
	cu_partition_info_t *sub_cu_info;

	int curr_cu_size = curr_cu_info->size;
	int abs_index = curr_cu_info->abs_index;

	sub_cu_info = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]];//4x4 blocks
	sub_cu_info += curr_cu_info->abs_index/sub_cu_info->num_part_in_cu;
	for( part_idx = abs_index; part_idx < (abs_index + curr_cu_info->num_part_in_cu); part_idx++,sub_cu_info++)
	{
		uint bs_check;
						
		if(num_pels_in_partition==4)//if( (g_uiMaxCUWidth >> g_uiMaxCUDepth) == 4 ) // I think this is always true. At least in homerHEVC
		{
			bs_check = ((dir==EDGE_VER && ((part_idx&0x1) == 0)) || (dir==EDGE_HOR && ((part_idx&0x2) == 0)));//bs_check signals whether the partition unit contains the right column or the bottom line
		}
		else
		{
			bs_check = 1;
		}

		if (et->deblock_edge_filter[dir][part_idx] && bs_check)
		{
			get_boundary_strength_single(et, currslice, ctu, curr_cu_info, sub_cu_info, dir);//xGetBoundaryStrengthSingle (enc_engine, ctu, curr_cu_info, sub_cu_info, currslice, dir);
//							xGetBoundaryStrengthSingle ( pcCU, dir, uiPartIdx );
		}
	}

	for(part_idx = 0; part_idx < pu_size_in_deblock_units; part_idx+=deblock_partition_idx_incr)
	{
		deblock_filter_luma(et, ctu, curr_cu_info, &et->enc_engine->curr_reference_frame->img, currslice, dir, part_idx);
		if ( (num_pels_in_partition>DEBLOCK_SMALLEST_BLOCK) || (part_idx % ( (DEBLOCK_SMALLEST_BLOCK<<1)/num_pels_in_partition ) ) == 0 )
		{
			deblock_filter_chroma(et, ctu, curr_cu_info, &et->enc_engine->curr_reference_frame->img, currslice, dir, part_idx);
//		  xEdgeFilterChroma   ( pcCU, uiAbsZorderIdx, uiDepth, iDir, iEdge );
		}
	}
}


//xSetEdgefilterPU
void set_edge_filter_pu(henc_thread_t* et, slice_t *currslice, ctu_info_t* ctu, cu_partition_info_t *curr_cu_info, int num_elements_width, int num_elements_height, int max_cu_width_units, int deblock_internal_edge)//, int dir)
{
//	int part_idx;
	int abs_x = ctu->x[Y_COMP]+curr_cu_info->x_position;
	int abs_y = ctu->y[Y_COMP]+curr_cu_info->y_position;
	int deblocking_filter_disable = currslice->deblocking_filter_disabled_flag;
	int left_edge = (!(abs_x==0 || deblocking_filter_disable));//left_edge indica si hay que filtrar la columna izq del PU o no
	int top_edge = (!(abs_y==0 || deblocking_filter_disable));//top_edge indica si hay que filtrar la fila top del PU o no
	int num_pels_in_partition = (et->max_cu_size>>et->max_cu_depth);
//	cu_partition_info_t *sub_cu_info;
	int curr_depth = curr_cu_info->size;
	int curr_cu_size = curr_cu_info->size;
	int abs_index = curr_cu_info->abs_index;

	//creo que esto solo habria que hacerlo si left_edge=!1
//	  xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, 0, m_stLFCUParam.bLeftEdge );
	set_edge_filter_0(et, curr_cu_info, et->deblock_edge_filter[EDGE_VER] , et->deblock_filter_strength_bs[EDGE_VER] , curr_depth, abs_index, curr_cu_size, curr_cu_size, EDGE_VER, left_edge, num_elements_width, max_cu_width_units);
	//creo que esto solo habria que hacerlo si top_edge=!1
	set_edge_filter_0(et, curr_cu_info, et->deblock_edge_filter[EDGE_HOR] , et->deblock_filter_strength_bs[EDGE_HOR] , curr_depth, abs_index, curr_cu_size, curr_cu_size, EDGE_HOR, top_edge, num_elements_height, max_cu_width_units);		

	switch (ctu->part_size_type[abs_index])
	{
		case SIZE_2Nx2N:
		{
			break;
		}
		case SIZE_2NxN:
		{
			break;							
		}
		case SIZE_NxN:
		{
			//creo que esto solo habria que hacerlo si left_edge=!1
			set_edge_filter_0_subdiv(et, curr_cu_info, et->deblock_edge_filter[EDGE_VER] , et->deblock_filter_strength_bs[EDGE_VER] , curr_depth, abs_index, curr_cu_size, curr_cu_size, EDGE_VER, deblock_internal_edge, num_elements_width, max_cu_width_units, num_elements_width>>1);
			//creo que esto solo habria que hacerlo si top_edge=!1
			set_edge_filter_0_subdiv(et, curr_cu_info, et->deblock_edge_filter[EDGE_HOR] , et->deblock_filter_strength_bs[EDGE_HOR] , curr_depth, abs_index, curr_cu_size, curr_cu_size, EDGE_HOR, deblock_internal_edge, num_elements_height, max_cu_width_units, num_elements_height>>1);		

//			xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, uiHWidthInBaseUnits, m_stLFCUParam.bInternalEdge );
//			xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, uiHHeightInBaseUnits, m_stLFCUParam.bInternalEdge );
			break;
		}
		//...
	}
}

void hmr_deblock_filter_cu(henc_thread_t* et, slice_t *currslice, ctu_info_t* ctu, int dir)
{
	int curr_depth = 0;
	cu_partition_info_t *curr_cu_info;//, *pred_cu_info;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int deblock_internal_edge = !currslice->deblocking_filter_disabled_flag;
	int curr_cu_size, abs_index;
	int num_elements_width, num_elements_height;// = (max(ctu->size, ctu->size)>>2);
	int max_cu_width_units = (et->max_cu_size>>2);
	int max_cu_height_units = (et->max_cu_size>>2);
//	int dir;

//	for(dir=EDGE_VER;dir<=EDGE_HOR;dir++)
	curr_cu_info = ctu->partition_list;//this should be mapped correctly

	//parent
	curr_depth = curr_cu_info->depth;
	memset(depth_state, 0, sizeof(depth_state));

	memset(et->deblock_filter_strength_bs[dir], 0, et->num_partitions_in_cu*sizeof(et->deblock_filter_strength_bs[dir][0]));
	memset(et->deblock_edge_filter[dir], 0, et->num_partitions_in_cu*sizeof(et->deblock_edge_filter[dir][0]));

	while(curr_depth!=0|| depth_state[curr_depth]!=1)
	{
		int pred_depth = ctu->pred_depth[curr_cu_info->abs_index];
		int tr_depth = ctu->tr_idx[curr_cu_info->abs_index];
		int total_depth = pred_depth + tr_depth;
		curr_cu_size = curr_cu_info->size;

		depth_state[curr_depth]++;
		if(curr_depth < total_depth && curr_cu_info->is_tl_inside_frame && curr_cu_size>DEBLOCK_SMALLEST_BLOCK)//go down one level
		{
			curr_depth++;
			curr_cu_info = curr_cu_info->children[depth_state[curr_depth]];
		}
		else//
		{
//			int l;
			if(curr_cu_info->is_tl_inside_frame)
			{
				abs_index = curr_cu_info->abs_index;
				num_elements_width = num_elements_height = (max(curr_cu_size, curr_cu_size)>>2);//(width, height)
				//num_elements_width = num_elements_height = (et->max_cu_size>>total_depth)>>2;//num_elements_height = (max(curr_cu_size, curr_cu_size)>>2);//(width, height)
				set_edge_filter_0(et, curr_cu_info, et->deblock_edge_filter[EDGE_VER] , et->deblock_filter_strength_bs[EDGE_VER] , curr_depth, abs_index, curr_cu_size, curr_cu_size, EDGE_VER, deblock_internal_edge, num_elements_width, max_cu_width_units);
				set_edge_filter_0(et, curr_cu_info, et->deblock_edge_filter[EDGE_HOR] , et->deblock_filter_strength_bs[EDGE_HOR] , curr_depth, abs_index, curr_cu_size, curr_cu_size, EDGE_HOR, deblock_internal_edge, num_elements_height, max_cu_height_units);				
	
				while(depth_state[curr_depth]==4 && curr_depth > pred_depth)
				{
					depth_state[curr_depth] = 0;
					curr_depth--;
					curr_cu_info = curr_cu_info->parent;
				}

				if(curr_depth == pred_depth)//xDeblockCU
				{
					//xSetEdgefilterPU
					abs_index = curr_cu_info->abs_index;
					num_elements_width = num_elements_height = (max(curr_cu_info->size, curr_cu_info->size)>>2);//(width, height)

					set_edge_filter_pu(et, currslice, ctu, curr_cu_info, num_elements_width, num_elements_height, max_cu_width_units, deblock_internal_edge);
					deblock_cu(et, currslice, ctu, curr_cu_info, dir);
				}
			}

			while(depth_state[curr_depth]==4)
			{
				depth_state[curr_depth] = 0;
				curr_depth--;
				curr_cu_info = curr_cu_info->parent;
			}
			
			if(curr_cu_info->parent != NULL)
				curr_cu_info = curr_cu_info->parent->children[depth_state[curr_depth]];
		}
	}
}


void hmr_deblock_filter(hvenc_engine_t* enc_engine, slice_t *currslice)
{
	int ctu_num;
	ctu_info_t* ctu;
//	int dir;

//	if(enc_engine->debug_file!=NULL)
//		wnd_write2file(&enc_engine->curr_reference_frame->img, enc_engine->debug_file);//debug

	//EDGE_VER = horizontal filter, EDGE_HOR = vertical filter

/*	for(dir=EDGE_VER;dir<=EDGE_HOR;dir++)
	{
		for(ctu_num = 0;ctu_num < enc_engine->pict_total_ctu;ctu_num++)
		{	
			ctu = &enc_engine->ctu_info[ctu_num];

			create_partition_ctu_neighbours(enc_engine->thread[0], ctu, ctu->partition_list);//this call should be removed

			hmr_deblock_filter_cu(enc_engine, currslice, ctu, dir);
		}
//		if(dir==EDGE_VER)
//			wnd_write2file(&enc_engine->curr_reference_frame->img, enc_engine->debug_file);//debug 		
	}
//	if(enc_engine->debug_file!=NULL)
//		wnd_write2file(&enc_engine->curr_reference_frame->img, enc_engine->debug_file);//debug 		

*/

//	int current_line=0;
	for(ctu_num = 0;ctu_num < enc_engine->pict_total_ctu;ctu_num++)
	{
		ctu = &enc_engine->ctu_info[ctu_num];
//		ctu->partition_list = enc_engine->thread[0]->deblock_partition_info;
//		create_partition_ctu_neighbours(enc_engine->thread[0], ctu, ctu->partition_list);//ctu->partition_list);//this call should be removed
		hmr_deblock_filter_cu(enc_engine->thread[0], currslice, ctu, EDGE_VER);
		if(ctu_num>=1*enc_engine->pict_width_in_ctu)
		{
			int ctu_num_horizontal = ctu_num-1*enc_engine->pict_width_in_ctu;
			ctu = &enc_engine->ctu_info[ctu_num_horizontal];
//			ctu->partition_list = enc_engine->thread[0]->deblock_partition_info;
//			create_partition_ctu_neighbours(enc_engine->thread[0], ctu, ctu->partition_list);//ctu->partition_list);//this call should be removed
			hmr_deblock_filter_cu(enc_engine->thread[0], currslice, ctu, EDGE_HOR);
		}
	}

	//finish the remaining ctus of the horizontal filter
	for(ctu_num = (enc_engine->pict_total_ctu-1*enc_engine->pict_width_in_ctu) ; ctu_num < enc_engine->pict_total_ctu ; ctu_num++)
	{
		ctu = &enc_engine->ctu_info[ctu_num];
//		ctu->partition_list = enc_engine->thread[0]->deblock_partition_info;
//		create_partition_ctu_neighbours(enc_engine->thread[0], ctu, ctu->partition_list);//ctu->partition_list);//this call should be removed
		hmr_deblock_filter_cu(enc_engine->thread[0], currslice, ctu, EDGE_HOR);		
	}
}




