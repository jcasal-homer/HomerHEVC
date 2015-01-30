/*****************************************************************************
 * hmr_motion_intra_chroma.c : homerHEVC encoding library
/*****************************************************************************
 * Copyright (C) 2014 homerHEVC project
 *
 * Juan Casal <jcasal.homer@gmail.com>
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
#include <math.h>
#include <limits.h>
#include <memory.h>
#include "hmr_common.h"
#include "hmr_private.h"


void synchronize_reference_buffs_chroma(henc_thread_t* et, cu_partition_info_t* curr_part, wnd_t *decoded_src, wnd_t * decoded_dst, int gcnt)
{
	int j, comp;
	int decoded_buff_stride;
	int16_t *  decoded_buff_src;
	int16_t *  decoded_buff_dst;

	for(comp=U_COMP;comp<=V_COMP;comp++)
	{
		decoded_buff_stride = WND_STRIDE_2D(*decoded_src, comp);//-curr_part->size;
		decoded_buff_src = WND_POSITION_2D(int16_t *, *decoded_src, comp, curr_part->x_position_chroma, curr_part->y_position_chroma, gcnt, et->ctu_width);
		decoded_buff_dst = WND_POSITION_2D(int16_t *, *decoded_dst, comp, curr_part->x_position_chroma, curr_part->y_position_chroma, gcnt, et->ctu_width);

		//bottom line
		memcpy(decoded_buff_dst+decoded_buff_stride*(curr_part->size_chroma-1), decoded_buff_src+decoded_buff_stride*(curr_part->size_chroma-1), curr_part->size_chroma*sizeof(decoded_buff_src[0]));

		//right column
		decoded_buff_src+=curr_part->size_chroma-1;
		decoded_buff_dst+=curr_part->size_chroma-1;
		for(j=0;j<curr_part->size_chroma-1;j++)
		{
			decoded_buff_dst[j*decoded_buff_stride] = decoded_buff_src[j*decoded_buff_stride];
		}
	}
}

//this function consolidate buffers from bottom to top
void synchronize_motion_buffers_chroma(henc_thread_t* et, cu_partition_info_t* curr_part, wnd_t * quant_src, wnd_t * quant_dst, wnd_t *decoded_src, wnd_t * decoded_dst, int gcnt)
{
	int j, comp;//, i;
	int decoded_buff_stride;
	int16_t *  decoded_buff_src;
	int16_t *  decoded_buff_dst;

	int quant_buff_stride = curr_part->size_chroma;//0;//es lineal
	int16_t *  quant_buff_src;
	int16_t *  quant_buff_dst;

	for(comp=U_COMP;comp<=V_COMP;comp++)
	{
		decoded_buff_stride = WND_STRIDE_2D(*decoded_src, comp);//-curr_part->size;
		decoded_buff_src = WND_POSITION_2D(int16_t *, *decoded_src, comp, curr_part->x_position_chroma, curr_part->y_position_chroma, gcnt, et->ctu_width);
		decoded_buff_dst = WND_POSITION_2D(int16_t *, *decoded_dst, comp, curr_part->x_position_chroma, curr_part->y_position_chroma, gcnt, et->ctu_width);

		quant_buff_src = WND_POSITION_1D(int16_t  *, *quant_src, comp, gcnt, et->ctu_width, (curr_part->abs_index<<et->num_partitions_in_cu_shift)>>2);
		quant_buff_dst = WND_POSITION_1D(int16_t  *, *quant_dst, comp, gcnt, et->ctu_width, (curr_part->abs_index<<et->num_partitions_in_cu_shift)>>2);

		for(j=0;j<curr_part->size_chroma;j++)
		{
			memcpy(quant_buff_dst, quant_buff_src, curr_part->size_chroma*sizeof(quant_buff_src[0]));
			memcpy(decoded_buff_dst, decoded_buff_src, curr_part->size_chroma*sizeof(decoded_buff_src[0]));
			quant_buff_src += quant_buff_stride;
			quant_buff_dst += quant_buff_stride;
			decoded_buff_src += decoded_buff_stride;
			decoded_buff_dst += decoded_buff_stride;
		}
	}
}

void create_chroma_dir_list(int* list, int luma_mode)
{
	int i;

	list[0] = PLANAR_IDX;
	list[1] = VER_IDX;
	list[2] = HOR_IDX;
	list[3] = DC_IDX;
	list[4] = DM_CHROMA_IDX;

	for( i = 0; i < NUM_CHROMA_MODE - 1; i++ )
	{
		if( luma_mode == list[i] )
		{
			list[i] = 34;
			break;
		}
	}
}	

extern const uint8_t chroma_scale_conversion_table[];

uint encode_intra_chroma(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position,  int part_size_type)
{
	int cu_mode, cu_mode_idx;
	uint distortion = 0, best_distortion=0, bit_cost, cost, best_cost = INT_MAX, best_mode, best_mode_idx;
	uint sum = 0, best_sum;
	picture_t *currpict = &et->ed->current_pict;
	slice_t *currslice = &currpict->slice;
	ctu_info_t* ctu_rd = et->ctu_rd;
	int pred_buff_stride, orig_buff_stride, residual_buff_stride, decoded_buff_stride;
	uint8_t * orig_buff;
	int16_t *pred_buff, *residual_buff, * quant_buff, * iquant_buff, *decoded_buff;
	uint8_t *cbf_buff[3] = {NULL,NULL,NULL};
	int num_candidates = 3;
	int best_pred_index = 0;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL;//, *resi_wnd = NULL;
	wnd_t *best_quant_wnd= NULL, *best_quant_wnd2= NULL, *best_decoded_wnd = NULL, *best_decoded_wnd2 = NULL;//, *best_resi_wnd = NULL, *best_resi_wnd2 = NULL;

	int best_pred_modes[3];
	double best_pred_cost[3]= {DOUBLE_MAX,DOUBLE_MAX,DOUBLE_MAX};
	uint best_cu_bitcost[3];

	int curr_part_size, curr_part_size_shift;
	int curr_adi_size;
	int curr_part_x, curr_part_y;
	int curr_depth;// = et->max_pred_partition_depth;
	cu_partition_info_t*	parent_part_info;
	cu_partition_info_t*	curr_partition_info;
	int curr_scan_mode;
	int curr_sum = 0;
	int partition_cost[4];
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int cbf_split[2][MAX_PARTITION_DEPTH] = {{0,0,0,0,0},{0,0,0,0,0}};
	int tr_depth_luma, original_depth;
	int  mode_list[NUM_CHROMA_MODE];
	int ch_component;
	double weight;
	int initial_state, end_state;
	int luma_mode;
	int qp_chroma, per, rem;// = chroma_scale_conversion_table[clip(curr_cu_info->qp,0,57)];
	int chr_qp_offset = et->ed->chroma_qp_offset;

	curr_partition_info = &ctu->partition_list[et->partition_depth_start[depth]]+part_position;

	weight = pow( 2.0, (currslice->qp-chroma_scale_conversion_table[clip(currslice->qp+chr_qp_offset,0,57)])/3.0 ); 
	qp_chroma = chroma_scale_conversion_table[clip(curr_partition_info->qp+chr_qp_offset,0,57)];
	per = qp_chroma/6;
	rem = qp_chroma%6;

//	ctu->qp_chroma = chroma_scale_conversion_table[clip(currslice->qp,0,57)];
//	ctu->per = ctu->qp_chroma/6;
//	ctu->rem = ctu->qp_chroma%6;

//	weight = pow( 2.0, (currslice->qp-ctu->qp_chroma)/3.0 ); 

	if(depth==0 && et->max_cu_size == MAX_CU_SIZE)
	{
		parent_part_info = &ctu->partition_list[et->partition_depth_start[depth]];
		curr_partition_info = parent_part_info->children[0];
	}
	else
	{
		curr_partition_info = &ctu->partition_list[et->partition_depth_start[depth]]+part_position;
		parent_part_info = curr_partition_info->parent;
	}


	luma_mode = et->intra_mode_buffs[Y_COMP][depth][curr_partition_info->abs_index] ;
	create_chroma_dir_list(mode_list, luma_mode);

#ifndef COMPUTE_AS_HM
	if(curr_partition_info->size_chroma == 2)//4 2x2 partitions are encoded as 1 4x4
	{
		curr_partition_info = parent_part_info;
		parent_part_info = curr_partition_info->parent;
	}

	decoded_wnd = &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1];//[curr_depth];
	for(cu_mode_idx=0;cu_mode_idx<NUM_CHROMA_MODE;cu_mode_idx++)
	{	
			distortion = cost = sum = 0;
			for(ch_component = U_COMP;ch_component<=V_COMP;ch_component++)
			{
				curr_part_x = curr_partition_info->x_position_chroma;
				curr_part_y = curr_partition_info->y_position_chroma;
				curr_depth = curr_partition_info->depth;
				curr_part_size = curr_partition_info->size_chroma;
				curr_part_size_shift = et->max_cu_size_shift-(curr_depth+1);
				curr_adi_size = 2*2*curr_part_size+1;

				pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd, ch_component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, ch_component);
				orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, ch_component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, ch_component);
				decoded_buff = WND_POSITION_2D(int16_t *, *decoded_wnd, ch_component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, ch_component);

				cu_mode = mode_list[cu_mode_idx];
				if(cu_mode == DM_CHROMA_IDX)//chroma direct mode infered from luma
				{
					cu_mode = luma_mode;
				}

				fill_reference_samples(et, ctu, curr_partition_info, curr_adi_size, decoded_buff-decoded_buff_stride-1, decoded_buff_stride, curr_part_size, CHR_COMP, FALSE);//don't create filtered adi samples

				//create prediction
				if(cu_mode== PLANAR_IDX)
					et->funcs->create_intra_planar_prediction(et, pred_buff, pred_buff_stride, et->adi_pred_buff, curr_adi_size, curr_part_size, curr_part_size_shift);//creamos el array de prediccion planar
				else
					et->funcs->create_intra_angular_prediction(et, ctu, pred_buff, pred_buff_stride, et->adi_pred_buff, curr_adi_size, curr_part_size, cu_mode, FALSE);//creamos el array de prediccion angular

				distortion += (double)et->funcs->sad(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride,curr_part_size);//R-D

				cost += distortion;
			}

			if(et->rd_mode==RD_FULL)			//info for rd
			{
//				int bit_cost;
				ctu_rd->intra_mode[CHR_COMP] = et->intra_mode_buffs[CHR_COMP][depth];

				memset(&et->intra_mode_buffs[CHR_COMP][depth][curr_partition_info->abs_index], mode_list[cu_mode_idx], curr_partition_info->num_part_in_cu*sizeof(et->intra_mode_buffs[CHR_COMP][depth][0]));

				bit_cost = (double)rd_estimate_bits_intra_mode(et, ctu_rd, curr_partition_info, depth-(part_size_type==SIZE_NxN), FALSE);
				cost += bit_cost*et->rd.sqrt_lambda;
			}
			else if(et->rd_mode != RD_FULL)
			{
//				int bit_cost = 0;

				if(mode_list[cu_mode_idx] == DM_CHROMA_IDX)
					bit_cost = 1;
				else
					bit_cost = 12;
				cost += bit_cost*et->rd.sqrt_lambda;
			}
			homer_update_cand_list( mode_list[cu_mode_idx], cost, bit_cost, best_pred_modes, best_pred_cost, best_cu_bitcost);
	}

	//depth = 0 keeps the best result
	cu_mode_idx=0;
//	for(cu_mode_idx=0;cu_mode_idx<1;cu_mode_idx++)
	{
		cu_mode = best_pred_modes[cu_mode_idx];
		bit_cost = best_cu_bitcost[cu_mode_idx];
#else	//#ifndef COMPUTE_AS_HM

	for(cu_mode_idx=0;cu_mode_idx<NUM_CHROMA_MODE;cu_mode_idx++)
	{
		cu_mode = mode_list[cu_mode_idx];

#endif //#ifndef COMPUTE_AS_HM
		memset(depth_state, 0, sizeof(depth_state));

		if(cu_mode == DM_CHROMA_IDX)//chroma direct mode infered from luma
		{
			cu_mode = luma_mode;
		}

		if(depth==0 && et->max_cu_size == MAX_CU_SIZE)
		{
			parent_part_info = &ctu->partition_list[et->partition_depth_start[depth]];
			curr_partition_info = parent_part_info->children[0];
		}
		else//(depth!=0 || et->max_cu_size != MAX_CU_SIZE)//if(part_size_type==SIZE_NxN)
		{
			curr_partition_info = &ctu->partition_list[et->partition_depth_start[depth]]+part_position;
			parent_part_info = curr_partition_info->parent;
		}

		curr_depth = curr_partition_info->depth;
		curr_part_size_shift = et->max_cu_size_shift-(curr_depth+1);
		distortion = 0;
		depth_state[curr_depth] = part_position&0x3;

		//get buffs
		quant_wnd = &et->transform_quant_wnd[NUM_QUANT_WNDS-1];//depth =4 buffs are not used to consolidate
		decoded_wnd = &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1];
		cbf_buff[U_COMP] = et->cbf_buffs_chroma[U_COMP];//temporal buffer
		cbf_buff[V_COMP] = et->cbf_buffs_chroma[V_COMP];//temporal buffer

		while(!(curr_depth==(depth-(part_size_type==SIZE_NxN)) && depth_state[curr_depth]==(part_position&0x3)+1))
		{
			curr_partition_info = (parent_part_info==NULL)?curr_partition_info:parent_part_info->children[depth_state[curr_depth]];

			tr_depth_luma = et->tr_idx_buffs[depth][curr_partition_info->abs_index]+depth-(part_size_type==SIZE_NxN);

			while(curr_depth<tr_depth_luma)
			{
				parent_part_info = curr_partition_info;
				curr_depth++;
				curr_partition_info = parent_part_info->children[depth_state[curr_depth]];
			}

			curr_scan_mode = find_scan_mode(TRUE, FALSE, curr_partition_info->size_chroma, cu_mode, 0);

			original_depth = curr_partition_info->depth;
			if(curr_partition_info->size_chroma == 2)//4 2x2 partitions are encoded as 1 4x4
			{
				curr_partition_info = parent_part_info;
				parent_part_info = curr_partition_info->parent;
			}

			curr_part_x = curr_partition_info->x_position_chroma;
			curr_part_y = curr_partition_info->y_position_chroma;
			curr_depth = curr_partition_info->depth;
			curr_part_size = curr_partition_info->size_chroma;
			curr_part_size_shift = et->max_cu_size_shift-(curr_depth+1);
			curr_adi_size = 2*2*curr_part_size+1;

			partition_cost[depth_state[curr_depth]] = 0;

			for(ch_component = U_COMP;ch_component<=V_COMP;ch_component++)
			{
				pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, ch_component);
				pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd, ch_component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, ch_component);
				orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, ch_component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, ch_component);
				residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, ch_component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, ch_component, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift)>>2);//420
				iquant_buff = WND_POSITION_1D(int16_t  *, et->itransform_iquant_wnd, ch_component, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift)>>2);//420
				decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, ch_component);
				decoded_buff = WND_POSITION_2D(int16_t *, *decoded_wnd, ch_component, curr_part_x, curr_part_y, gcnt, et->ctu_width);

				fill_reference_samples(et, ctu, curr_partition_info, curr_adi_size, decoded_buff-decoded_buff_stride-1, decoded_buff_stride, curr_part_size, CHR_COMP, FALSE);//don't create filtered adi samples

				if(cu_mode== PLANAR_IDX)
					et->funcs->create_intra_planar_prediction(et, pred_buff, pred_buff_stride, et->adi_pred_buff, curr_adi_size, curr_part_size, curr_part_size_shift);//creamos el array de prediccion planar
				else
					et->funcs->create_intra_angular_prediction(et, ctu, pred_buff, pred_buff_stride, et->adi_pred_buff, curr_adi_size, curr_part_size, cu_mode, FALSE);//creamos el array de prediccion angular

				//intra code
				et->funcs->predict(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride, residual_buff, residual_buff_stride, curr_part_size);
				et->funcs->transform(et->bit_depth, residual_buff, et->pred_aux_buff, residual_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, curr_part_size_shift, REG_DCT, quant_buff);//usamos quant buff como auxiliar
				et->funcs->quant(et, et->pred_aux_buff, quant_buff, curr_scan_mode, curr_depth, ch_component, cu_mode, 1, &curr_sum, curr_part_size, per, rem);//Si queremos quitar el bit de signo necesitamos hacerlo en dos arrays distintos

				sum+=curr_sum;
				//set cbf
				memset(&cbf_buff[ch_component][curr_partition_info->abs_index], ((curr_sum ? 1 : 0) << (original_depth-depth+(part_size_type==SIZE_NxN)))|((curr_sum ? 1 : 0) << (curr_depth-depth+(part_size_type==SIZE_NxN))), curr_partition_info->num_part_in_cu*sizeof(cbf_buff[ch_component][0]));//(width*width)>>4 num parts of 4x4 in partition
				cbf_split[ch_component-1][curr_depth] |= (curr_sum ? 1 : 0);

				if(curr_sum)
				{
					//intra decode 
					et->funcs->inv_quant(et, quant_buff, iquant_buff, curr_depth, ch_component, 1, curr_part_size, per, rem);
					et->funcs->itransform(et->bit_depth, residual_buff, iquant_buff, residual_buff_stride, curr_part_size, curr_part_size, REG_DCT, et->pred_aux_buff);
					et->funcs->reconst(pred_buff, pred_buff_stride, residual_buff, residual_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
				}
				else
				{
					et->funcs->reconst(pred_buff, pred_buff_stride, quant_buff, 0, decoded_buff, decoded_buff_stride, curr_part_size);
				}
				partition_cost[depth_state[curr_depth]] += (int)(weight*et->funcs->ssd(orig_buff, orig_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size));//R-D
			}//for(ch_component = U_COMP;ch_component<=V_COMP;ch_component++)

			curr_partition_info->sum += sum;
			distortion += partition_cost[depth_state[curr_depth]];

#ifndef COMPUTE_AS_HM
			if(distortion>best_cost)
			{
				distortion = best_cost+1;
				break;
			}

#endif
			depth_state[curr_depth]++;

			if(depth_state[curr_depth]==4)
			{	
				while(depth_state[curr_depth]==4 && curr_depth>(depth-(part_size_type==SIZE_NxN)))
				{
					int ll;
					for(ll=parent_part_info->abs_index;ll<parent_part_info->abs_index+parent_part_info->num_part_in_cu;ll++)
					{
						cbf_buff[U_COMP][ll] |= (cbf_split[U_COMP-1][curr_depth]<<(curr_depth-1-depth+(part_size_type==SIZE_NxN)));
						cbf_buff[V_COMP][ll] |= (cbf_split[V_COMP-1][curr_depth]<<(curr_depth-1-depth+(part_size_type==SIZE_NxN)));
					}
					cbf_split[U_COMP-1][curr_depth-1] |= cbf_split[U_COMP-1][curr_depth];
					cbf_split[V_COMP-1][curr_depth-1] |= cbf_split[V_COMP-1][curr_depth];

					cbf_split[U_COMP-1][curr_depth] = cbf_split[V_COMP-1][curr_depth] = 0;

					depth_state[curr_depth] = 0;
					curr_depth--;
					depth_state[curr_depth]++;
					if(curr_depth!=0)
						parent_part_info = parent_part_info->parent;
				}									
			}
		}//while(curr_depth!=1 || depth_state[curr_depth]!=4)

		//curr_partition_info top most element of the recursive tree. 
		if(depth==0)
			curr_partition_info = &ctu->partition_list[et->partition_depth_start[depth]];
		else
		{
			curr_partition_info = (&ctu->partition_list[et->partition_depth_start[depth]]+part_position);
			if(part_size_type==SIZE_NxN)
				curr_partition_info = curr_partition_info->parent;
		}

		cost = distortion;

		if(et->rd_mode == RD_FULL && cost<best_cost)
		{
			//rd
			ctu_rd->cbf[U_COMP] = cbf_buff[U_COMP];
			ctu_rd->cbf[V_COMP] = cbf_buff[V_COMP];
			ctu_rd->intra_mode[U_COMP] = et->intra_mode_buffs[Y_COMP][depth];
			ctu_rd->intra_mode[CHR_COMP] = et->intra_mode_buffs[CHR_COMP][depth];
			ctu_rd->tr_idx = et->tr_idx_buffs[depth];
			ctu_rd->coeff_wnd = quant_wnd;

#ifndef COMPUTE_AS_HM
			memset(&et->intra_mode_buffs[CHR_COMP][depth][curr_partition_info->abs_index], best_pred_modes[cu_mode_idx], curr_partition_info->num_part_in_cu*sizeof(et->intra_mode_buffs[CHR_COMP][depth][0]));
#else
			memset(&et->intra_mode_buffs[CHR_COMP][depth][curr_partition_info->abs_index], mode_list[cu_mode_idx], curr_partition_info->num_part_in_cu*sizeof(et->intra_mode_buffs[CHR_COMP][depth][0]));
#endif
			bit_cost = rd_get_intra_bits_qt(et, ctu_rd, curr_partition_info, depth, FALSE, gcnt);
			cost += bit_cost*et->rd.lambda+.5;
		}
		else if(et->rd_mode != RD_FULL && cost<best_cost)
		{
			double correction = calc_mv_correction(curr_partition_info->qp, et->ed->avg_dist);//.25+et->ed->avg_dist*et->ed->avg_dist/5000000.;
			//cost += bit_cost*curr_partition_info->qp/clip((3500000/(et->ed->avg_dist*et->ed->avg_dist)),.35,4.)+.5;
			cost += bit_cost*correction+.5;
		}

		if(cost < best_cost)
		{
			best_distortion = distortion;
			best_cost = cost;
			best_sum = sum;
#ifndef COMPUTE_AS_HM
			best_mode = best_pred_modes[cu_mode_idx];
#else
			best_mode = mode_list[cu_mode_idx];
#endif
			best_mode_idx = cu_mode_idx;

//			curr_partition_info->distortion_chroma = best_distortion;
//			curr_partition_info->cost_chroma = cost;
//			curr_partition_info->mode_chroma = cu_mode;

			//synchronize buffers for next iterations esto en verdad no seria necesario, quedaria todo en depth=1
			synchronize_motion_buffers_chroma(et, curr_partition_info, quant_wnd, &et->transform_quant_wnd[depth+1], decoded_wnd, &et->decoded_mbs_wnd[depth+1], gcnt);
			memcpy(&et->cbf_buffs[U_COMP][depth][curr_partition_info->abs_index], &cbf_buff[U_COMP][curr_partition_info->abs_index], curr_partition_info->num_part_in_cu*sizeof(et->cbf_buffs[U_COMP][depth][0]));
			memcpy(&et->cbf_buffs[V_COMP][depth][curr_partition_info->abs_index], &cbf_buff[V_COMP][curr_partition_info->abs_index], curr_partition_info->num_part_in_cu*sizeof(et->cbf_buffs[U_COMP][depth][0]));
		}
	}//for(cu_mode_idx=0;cu_mode_idx<NUM_CHROMA_MODE;cu_mode_idx++)

	memset(&et->intra_mode_buffs[CHR_COMP][depth][curr_partition_info->abs_index], best_mode, curr_partition_info->num_part_in_cu*sizeof(et->intra_mode_buffs[CHR_COMP][depth][0]));
	curr_partition_info->sum+=best_sum;
	return best_cost;//curr_partition_info->cost_chroma;
}


