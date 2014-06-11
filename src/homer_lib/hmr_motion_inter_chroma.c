/*****************************************************************************
 * hmr_motion_inter_chroma.c : homerHEVC encoding library
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


#include "hmr_private.h"
#include "hmr_common.h"
#include "hmr_profiler.h"

#include "hmr_sse42_functions.h"


#ifdef ENABLE___

int encode_inter_cu_chroma(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, int component, int depth, int cu_mode, PartSize part_size_type, int *curr_sum, motion_vector_t *mv, int gcnt)//depth = prediction depth
{		
	int ssd_;
	int pred_buff_stride, orig_buff_stride, residual_buff_stride, residual_dec_buff_stride, decoded_buff_stride;
	uint8_t *pred_buff, *orig_buff, *decoded_buff;
	int16_t*residual_buff, *residual_dec_buff, *quant_buff, *iquant_buff;
	uint8_t *cbf_buff = NULL;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL;
	int inv_depth, diff;//, is_filtered;
		
	int curr_depth = curr_partition_info->depth;
	int curr_part_x = curr_partition_info->x_position;
	int curr_part_y = curr_partition_info->y_position;
	int curr_part_size = curr_partition_info->size;
	int curr_part_size_shift = et->max_cu_size_shift-curr_depth;
	int curr_adi_size = 2*2*curr_part_size+1;

	int curr_scan_mode = find_scan_mode(TRUE, TRUE, curr_part_size, cu_mode, 0);

	quant_wnd = &et->transform_quant_wnd[curr_depth];
	decoded_wnd = &et->decoded_mbs_wnd[curr_depth];
	cbf_buff = et->cbf_buffs[component][curr_depth];

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, component);
	pred_buff = WND_POSITION_2D(uint8_t *, et->prediction_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, component);
	orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, component);
	residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_dec_buff_stride = WND_STRIDE_2D(et->residual_dec_wnd, Y_COMP);
	residual_dec_buff = WND_POSITION_2D(int16_t *, et->residual_dec_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, component, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));
	iquant_buff = WND_POSITION_1D(int16_t  *, et->itransform_iquant_wnd, component, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));
	decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, component);
	decoded_buff = WND_POSITION_2D(uint8_t *, *decoded_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, component, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));

	inv_depth = (et->max_cu_size_shift - curr_depth);
	diff = min(abs(cu_mode - HOR_IDX), abs(cu_mode - VER_IDX));

	//2d -> 1D buffer
	PROFILER_RESET(inter_luma_tr)
	et->funcs->transform(et->bit_depth, residual_buff, et->pred_aux_buff, residual_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, curr_part_size_shift, cu_mode, quant_buff);//usamos quant buff como auxiliar
	PROFILER_ACCUMULATE(inter_luma_tr)

	PROFILER_RESET(inter_luma_q)
	et->funcs->quant(et, ctu, et->pred_aux_buff, quant_buff, curr_scan_mode, curr_depth, component, cu_mode, 0, curr_sum, curr_part_size);//Si queremos quitar el bit de signo necesitamos hacerlo en dos arrays distintos
	curr_partition_info->sum = *curr_sum;
	PROFILER_ACCUMULATE(inter_luma_q)


	curr_partition_info->inter_cbf[component] = (( curr_partition_info->sum ? 1 : 0 ) << (curr_depth-depth+(part_size_type==SIZE_NxN)));
	curr_partition_info->inter_tr_idx[component]  = (curr_depth-depth+(part_size_type==SIZE_NxN));
//	memset(&et->cbf_buffs[component][curr_depth][curr_partition_info->abs_index], (( curr_partition_info->sum ? 1 : 0 ) << (curr_depth-depth+(part_size_type==SIZE_NxN))), curr_partition_info->num_part_in_cu*sizeof(cbf_buff[0]));//(width*width)>>4 num parts of 4x4 in partition
//	memset(&et->tr_idx_buffs[curr_depth][curr_partition_info->abs_index], (curr_depth-depth+(part_size_type==SIZE_NxN)), curr_partition_info->num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));//(width*width)>>4 num parts of 4x4 in partition
//	memset(&et->intra_mode_buffs[component][curr_depth][curr_partition_info->abs_index], cu_mode, curr_partition_info->num_part_in_cu*sizeof(et->intra_mode_buffs[component][0][0]));//(width*width)>>4 num parts of 4x4 in partition

	et->funcs->inv_quant(et, ctu, quant_buff, iquant_buff, curr_depth, component, 0, curr_part_size);

	//1D ->2D buffer
	et->funcs->itransform(et->bit_depth, residual_dec_buff, iquant_buff, residual_buff_stride, curr_part_size, curr_part_size, cu_mode, et->pred_aux_buff);
	ssd_ = ssd_16(residual_buff, residual_buff_stride, residual_dec_buff, residual_buff_stride, curr_part_size);

	et->funcs->reconst(pred_buff, pred_buff_stride, residual_dec_buff, residual_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
//	ssd_ = et->funcs->ssd(orig_buff, orig_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	return ssd_;
}


/*void hmr_motion_inter_uni(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint8_t *orig_buff, int orig_buff_stride, uint8_t *reference_buff, int reference_buff_stride, uint8_t *pred_buff, int pred_buff_stride,  int curr_part_global_x, int curr_part_global_y, int curr_part_size, int curr_part_size_shift, motion_vector_t *mv)
{
	int j;
	mv->hor_vector = 0;
	mv->ver_vector = 0;

	for(j=0;j<curr_part_size;j++)
	{
		memcpy(pred_buff, reference_buff, curr_part_size);
		pred_buff += pred_buff_stride;
		reference_buff += reference_buff_stride;
	}
}
*/
int encode_inter_chroma(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type)
{
	int k;
	int cu_mode;
	double distortion = 0.;
	uint cost, best_cost;
	slice_t *currslice = &et->ed->current_pict.slice;
	ctu_info_t *ctu_rd = et->ctu_rd;
	int pred_buff_stride, orig_buff_stride, decoded_buff_stride, reference_buff_stride, residual_buff_stride;
	uint8_t *pred_buff, *orig_buff, *decoded_buff, *reference_buff, *reference_buff_cu_position;
	uint8_t *cbf_buff = NULL;//, *best_cbf_buff = NULL, *best_cbf_buff2 = NULL;
	int16_t *residual_buff;
	wnd_t *decoded_wnd = NULL, *reference_wnd=NULL;//, *resi_wnd = NULL;

	int curr_part_size, curr_part_size_shift;
	int curr_adi_size;
	int curr_part_x, curr_part_y, curr_part_global_x, curr_part_global_y;
	int curr_depth = depth;
	cu_partition_info_t*	parent_part_info;
	cu_partition_info_t*	curr_partition_info;
	int curr_sum = 0, best_sum;
	int num_part_in_cu;
	int partition_cost;
	int cu_min_tu_size_shift;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int max_tr_depth, max_tr_processing_depth;
	int initial_state, end_state;
//	int cbf_split[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int acc_cost[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int bitcost_cu_mode;
	int log2cu_size;
	int cu_x_position;
	motion_vector_t mv;
	int component;

	//rate-control - en motion_intra_chroma se modifica
	ctu->qp = currslice->qp;
	ctu->per = ctu->qp/6;
	ctu->rem = ctu->qp%6;

	PROFILER_RESET(inter_chroma)
	for(component=U_COMP;component<=V_COMP;component++)
	{
		curr_partition_info = &ctu->partition_list[et->partition_depth_start[depth]]+part_position;

		parent_part_info = curr_partition_info->parent;

		curr_depth = curr_partition_info->depth;
		curr_part_x = curr_partition_info->x_position;
		curr_part_y = curr_partition_info->y_position;
		curr_part_global_x = ctu->x[CHR_COMP]+curr_part_x;
		curr_part_global_y = ctu->y[CHR_COMP]+curr_part_y;
		curr_part_size = curr_partition_info->size;
		curr_part_size_shift = et->max_cu_size_shift-curr_depth;
		curr_adi_size = 2*2*curr_part_size+1;

		pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, component);
		pred_buff = WND_POSITION_2D(uint8_t *, et->prediction_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
		orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, component);
		orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);

		decoded_wnd = &et->decoded_mbs_wnd[depth];	
		decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, component);
		decoded_buff = WND_POSITION_2D(uint8_t *, *decoded_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);

		residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, component);
		residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);

		reference_wnd = &currslice->ref_pic_list[REF_PIC_LIST_0][0]->img;//[0] up to now we only use one reference	
		reference_buff_stride = WND_STRIDE_2D(*reference_wnd, component);
		reference_buff = WND_POSITION_2D(uint8_t *, *reference_wnd, component, 0, 0, gcnt, et->ctu_width);
		reference_buff_cu_position = WND_POSITION_2D(uint8_t *, *reference_wnd, component, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);

		if(depth==0 && et->max_cu_size == MAX_CU_SIZE)
		{
			parent_part_info = &ctu->partition_list[et->partition_depth_start[depth]];
			curr_partition_info = parent_part_info->children[0];
			initial_state = part_position;
			end_state = initial_state+4;
		}
		else
		{
			curr_partition_info = &ctu->partition_list[et->partition_depth_start[curr_depth]]+part_position;
			parent_part_info = curr_partition_info->parent;
			initial_state = part_position & 0x3;
			end_state = initial_state+1;
		}

		best_cost = INT_MAX;
		curr_depth = curr_partition_info->depth;
		memset(depth_state, 0, sizeof(depth_state));

		PROFILER_RESET(inter_luma_bucle1)
		hmr_motion_inter_uni(et, ctu, curr_partition_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, pred_buff, pred_buff_stride,  curr_part_global_x, curr_part_global_y, curr_part_size, curr_part_size_shift, &mv);
		//hmr_motion_search_uni(et, ctu, curr_partition_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, pred_buff, pred_buff_stride, curr_part_global_x, curr_part_global_y, curr_part_size, curr_part_size_shift, &mv);
		PROFILER_ACCUMULATE(inter_luma_bucle1)

		PROFILER_RESET(inter_luma_predict)
		et->funcs->predict(orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, residual_buff, residual_buff_stride, curr_part_size);
		PROFILER_ACCUMULATE(inter_luma_predict)

		cu_mode = 0;//best_mode;

		curr_depth = curr_partition_info->depth;

		max_tr_depth = et->max_inter_tr_depth;

		log2cu_size = et->max_cu_size_shift-(depth-(part_size_type==SIZE_NxN));

		max_tr_processing_depth = 0;

		if(curr_part_size_shift > et->max_tu_size_shift)
		{
			while(curr_part_size_shift > (et->max_tu_size_shift+max_tr_processing_depth)) max_tr_processing_depth++;
		}
		max_tr_processing_depth++;

		while((curr_part_size_shift-max_tr_processing_depth) < (et->max_cu_size_shift-et->max_cu_depth)) max_tr_processing_depth--;

	//	max_tr_processing_depth = et->max_cu_size_shift-cu_min_tu_size_shift;

		memset(acc_cost, 0, sizeof(acc_cost));
		memset(depth_state, 0, sizeof(depth_state));

		depth_state[curr_depth] = initial_state;

		while(curr_depth!=depth || depth_state[curr_depth]!=end_state)
		{		
			uint bit_cost = 0;
			curr_partition_info = (parent_part_info==NULL)?curr_partition_info:parent_part_info->children[depth_state[curr_depth]];//if cu_size=64 we process 4 32x32 partitions, else just the curr_partition
			curr_depth = curr_partition_info->depth;

			curr_partition_info->distortion_chroma = encode_inter_cu_chroma(et, ctu, curr_partition_info, component, depth, cu_mode, part_size_type, &curr_sum, &mv, gcnt);//depth = prediction depth

			curr_partition_info->cost_chroma = curr_partition_info->distortion_chroma;

			acc_cost[curr_depth] += curr_partition_info->distortion_chroma;

			depth_state[curr_depth]++;

			//HM order
			if(curr_depth < max_tr_processing_depth)//go one level down
			{
				curr_depth++;
				parent_part_info = curr_partition_info;
			}
			else if(depth_state[curr_depth]==4)//consolidate 
			{	
				while(depth_state[curr_depth]==4 && (curr_depth > (depth)))
				{
					distortion = parent_part_info->children[0]->distortion+parent_part_info->children[1]->distortion+parent_part_info->children[2]->distortion+parent_part_info->children[3]->distortion;
					cost = distortion;

					depth_state[curr_depth] = 0;

					if(cost < parent_part_info->cost)
					{
						acc_cost[parent_part_info->depth]+= cost - parent_part_info->cost;
						parent_part_info->cost = cost;
						parent_part_info->distortion = distortion;

						//consolidate in parent
						memcpy(&et->cbf_buffs[component][parent_part_info->depth][parent_part_info->abs_index], &et->cbf_buffs[component][curr_depth][parent_part_info->abs_index], parent_part_info->num_part_in_cu*sizeof(cbf_buff[0]));
	//					memcpy(&et->intra_mode_buffs[component][parent_part_info->depth][parent_part_info->abs_index], &et->intra_mode_buffs[component][curr_depth][parent_part_info->abs_index], parent_part_info->num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));
						memcpy(&et->tr_idx_buffs[parent_part_info->depth][parent_part_info->abs_index], &et->tr_idx_buffs[curr_depth][parent_part_info->abs_index], parent_part_info->num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));

						//synchronize buffers for next iterations
	//					synchronize_motion_buffers_luma(et, parent_part_info, &et->transform_quant_wnd[curr_depth], &et->transform_quant_wnd[curr_depth-1], &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[curr_depth-1], gcnt);
					}
					else
					{
	//					synchronize_reference_buffs(et, parent_part_info, &et->decoded_mbs_wnd[curr_depth-1], &et->decoded_mbs_wnd[curr_depth], gcnt);	
					}

					acc_cost[curr_depth] = 0;
					curr_depth--;
					parent_part_info = parent_part_info->parent;
				}

				if(curr_depth+2 <= max_tr_processing_depth)//sync buffers of higher depth 
				{
					int aux_depth, aux_abs_index, aux_num_part_in_cu;
					cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
					aux_abs_index = aux_partition_info->abs_index;
					aux_num_part_in_cu = aux_partition_info->num_part_in_cu;
					for(aux_depth=curr_depth+2;aux_depth<=max_tr_processing_depth;aux_depth++)//+1 pq el nivel superior esta ya sincronizado y tenemos que sincronizar los siguientes
					{
	//					synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);	
	//					if(et->rd_mode==1)
						{
	//						memcpy(&et->intra_mode_buffs[component][aux_depth][aux_abs_index], &et->intra_mode_buffs[component][curr_depth][aux_abs_index], aux_num_part_in_cu*sizeof(et->intra_mode_buffs[component][0][0]));
							memcpy(&et->cbf_buffs[component][aux_depth][aux_abs_index], &et->cbf_buffs[component][curr_depth][aux_abs_index], aux_num_part_in_cu*sizeof(et->cbf_buffs[component][0][0]));
							memcpy(&et->tr_idx_buffs[aux_depth][aux_abs_index], &et->tr_idx_buffs[curr_depth][aux_abs_index], aux_num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));
						}
					}
				}
			}
		}
	}
	PROFILER_ACCUMULATE(inter_chroma)

	return curr_partition_info->cost;
}

#endif
