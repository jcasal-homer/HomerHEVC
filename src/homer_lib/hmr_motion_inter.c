/*****************************************************************************
 * hmr_motion_intra.c : homerHEVC encoding library
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
 * You should have received a copy of the GNU General Public License
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


int motion_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt)
{
	double cost, cost_luma, cost_chroma, best_cost;
	int position = 0;
	int curr_depth = 0;
	PartSize part_size_type = SIZE_2Nx2N;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL;
	ctu_info_t *ctu_rd = et->ctu_rd;
	cu_partition_info_t	*parent_part_info = NULL;
	cu_partition_info_t	*curr_partition_info = ctu->partition_list;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	uint cost_sum[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int abs_index;
	int num_part_in_cu;
	uint8_t cbf_split[NUM_PICT_COMPONENTS];
	int fast_skip = 0;
	int ll;

	//init rd auxiliar ctu
	if(et->rd_mode==1)
	{
		copy_ctu(ctu, ctu_rd);
		ctu_rd->pred_mode = INTRA_MODE;
	}

#ifndef COMPUTE_AS_HM
	if(et->performance_mode != 0)
		analyse_intra_recursive_info(et, ctu, gcnt);
#endif

	memset(cbf_split, 0, sizeof(cbf_split));
	while(curr_depth!=0 || depth_state[curr_depth]!=1)
	{
		fast_skip = 0;
		curr_depth = curr_partition_info->depth;
		num_part_in_cu = curr_partition_info->num_part_in_cu;
		abs_index = curr_partition_info->abs_index;
		part_size_type = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//
		position = curr_partition_info->list_index - et->partition_depth_start[curr_depth];
		quant_wnd = &et->transform_quant_wnd[curr_depth];
		decoded_wnd = &et->decoded_mbs_wnd[curr_depth];
		
		cost_luma = cost_chroma = 0;

/*		if(curr_partition_info->is_b_inside_frame && curr_partition_info->is_r_inside_frame)//if br (and tl) are inside the frame, process
		{
#ifndef COMPUTE_AS_HM

//			if(curr_partition_info->recursive_split)//con esto activado deberia computar solo las hojas del analisis de recursividad. Se ve bien pero el resultado es distinto con rd que sin el, no hace exactamente lo que quiero
//			if((et->performance_mode == 1 && curr_partition_info->size == 64 && curr_partition_info->recursive_split) || (et->performance_mode == 2 && curr_partition_info->recursive_split))// && curr_depth<ed->max_pred_partition_depth && curr_partition_info->children[0]->recursive_split))//we skip 64x64 as it is not used very often
			if((et->performance_mode == 1 && curr_partition_info->recursive_split && curr_partition_info->children[0] && curr_partition_info->children[0]->recursive_split && curr_partition_info->children[1]->recursive_split && curr_partition_info->children[2]->recursive_split && curr_partition_info->children[3]->recursive_split) ||
				(et->performance_mode == 2 && curr_partition_info->recursive_split))
			{
				cost_luma = UINT_MAX;
				cost_chroma = 0;//UINT_MAX;
			}
			else

#endif
			{
				//encode
				PROFILER_RESET(intra_luma)
				cost_luma = encode_intra_luma(et, ctu, gcnt, curr_depth, position, part_size_type);
				PROFILER_ACCUMULATE(intra_luma)

				cost_chroma = 0;
				if(part_size_type == SIZE_2Nx2N)
				{
					int iiiiii;
					PROFILER_RESET(intra_chroma)
					cost_chroma = encode_intra_chroma(et, ctu, gcnt, curr_depth, position, part_size_type);//encode_intra_chroma(et, gcnt, depth, position, part_size_type);	
					PROFILER_ACCUMULATE(intra_chroma)
					iiiiii=0;
				}
			}
		}
*/
		curr_partition_info->cost = cost_luma+cost_chroma;//+cost_bits*et->rd.sqrt_lambda;;

		depth_state[curr_depth]++;

		cost_sum[curr_depth]+=curr_partition_info->cost;

#ifndef COMPUTE_AS_HM
		//if this matches, it is useless to continue the recursion. the case where part_size_type != SIZE_NxN is checked at the end of the consolidation buffer)
		if(/*et->performance_mode == 0 && */part_size_type == SIZE_NxN && depth_state[curr_depth]!=4 && cost_sum[curr_depth] > parent_part_info->cost && ctu->partition_list[0].is_b_inside_frame && ctu->partition_list[0].is_r_inside_frame)//parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
		{
			depth_state[curr_depth]=4;
		}

		//stop recursion
		if(et->performance_mode>0 && (curr_partition_info->recursive_split==0 || curr_partition_info->cost == 0) && /*curr_depth && */part_size_type != SIZE_NxN && curr_partition_info->is_b_inside_frame && curr_partition_info->is_r_inside_frame)
		{
			int max_processing_depth;// = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

			//if we fill this in here we don't have to consolidate
			memset(&ctu->pred_depth[abs_index], curr_depth-(part_size_type==SIZE_NxN), num_part_in_cu*sizeof(ctu->pred_depth[0]));
			memset(&ctu->part_size_type[abs_index], part_size_type, num_part_in_cu*sizeof(ctu->part_size_type[0]));
			if(et->rd_mode==1)//rd
			{
				memset(&ctu_rd->pred_depth[abs_index], curr_depth-(part_size_type==SIZE_NxN), num_part_in_cu*sizeof(ctu_rd->pred_depth[0]));
				memset(&ctu_rd->part_size_type[abs_index], part_size_type, num_part_in_cu*sizeof(ctu_rd->part_size_type[0]));
			}
			cost_chroma = 0;//only!=0 if part_size_type == SIZE_NxN
			fast_skip = 1;

			max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);
			if(curr_depth+1 <= max_processing_depth)//el = es para cuando et->max_intra_tr_depth!=4
			{
				int aux_depth;
				cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
				abs_index = aux_partition_info->abs_index;
				num_part_in_cu  = aux_partition_info->num_part_in_cu;

				for(aux_depth=curr_depth+1;aux_depth<=max_processing_depth;aux_depth++)//+1 pq el nivel superior esta ya sincronizado y tenemos que sincronizar los siguientes
				{
					synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);	
					synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);
					//for rd
					consolidate_recursive_info_buffers(et, gcnt, parent_part_info, curr_depth, aux_depth, abs_index, num_part_in_cu);
				}
				synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
			}
		}
#endif 

		if(!fast_skip && curr_depth<et->max_pred_partition_depth && curr_partition_info->is_tl_inside_frame)//depth_state[curr_depth]!=4 is for fast skip//if tl is not inside the frame don't process the next depths
		{
			curr_depth++;
			parent_part_info = curr_partition_info;
		}
		else if(/*fast_skip || */depth_state[curr_depth]==4)//la depth =1 lo hemos consolidado antes del bucle
		{
			int max_processing_depth;
			cost_chroma = 0;

/*			if(cost_sum[curr_depth] < parent_part_info->cost &&  part_size_type == SIZE_NxN && (curr_partition_info->is_b_inside_frame && curr_partition_info->is_r_inside_frame))//esto se tiene que hacer aqui pq si sno se eligen cu_modes distintos para cada una de las 4 particiones
			{
				position = parent_part_info->children[0]->list_index - et->partition_depth_start[curr_depth];
				PROFILER_RESET(intra_chroma)
				cost_chroma = encode_intra_chroma(et, ctu, gcnt, curr_depth, position, part_size_type);//encode_intra_chroma(et, gcnt, depth, position, part_size_type);	
				PROFILER_ACCUMULATE(intra_chroma)
			}
*/			while(depth_state[curr_depth]==4 && curr_depth>0)//>0 pq consolidamos sobre el padre, 
			{
/*				cost_luma = parent_part_info->children[0]->cost + parent_part_info->children[1]->cost +parent_part_info->children[2]->cost+parent_part_info->children[3]->cost;

				depth_state[curr_depth] = 0;
				cost = cost_luma+cost_chroma;//cost_chroma=0 after the first iteration 
				cost_sum[curr_depth]+=cost_chroma;
				best_cost = parent_part_info->cost;

				abs_index = parent_part_info->abs_index;
				num_part_in_cu = parent_part_info->num_part_in_cu;

				//choose best
				if(cost<best_cost || !(curr_partition_info->is_b_inside_frame && curr_partition_info->is_r_inside_frame))//if we get here, tl should be inside the frame
				{
					//Aqui consolidamos los resultados bottom-up
					//here we consolidate the bottom-up results for being preferred to the top-down computation
					int part_size_type2 = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//

					cost_sum[parent_part_info->depth] -= parent_part_info->cost;
					cost_sum[parent_part_info->depth] += cost;
					parent_part_info->cost = cost;
					quant_wnd = &et->transform_quant_wnd[curr_depth];//&et->transform_quant_wnd[curr_depth+1];
					decoded_wnd = &et->decoded_mbs_wnd[curr_depth];//&et->decoded_mbs_wnd[curr_depth+1];
			
					//consolidate in parent
					synchronize_motion_buffers_luma(et, parent_part_info, quant_wnd, &et->transform_quant_wnd[curr_depth-1], decoded_wnd, &et->decoded_mbs_wnd[curr_depth-1], gcnt);
					synchronize_motion_buffers_chroma(et, parent_part_info, quant_wnd, &et->transform_quant_wnd[curr_depth-1], decoded_wnd, &et->decoded_mbs_wnd[curr_depth-1], gcnt);

					consolidate_recursive_info_buffers(et, gcnt, parent_part_info, curr_depth, curr_depth-1, abs_index, num_part_in_cu);

					if(part_size_type2==SIZE_NxN)
					{
						int num_part_in_sub_cu = parent_part_info->children[0]->num_part_in_cu;
						cbf_split[Y_COMP] = (et->cbf_buffs[Y_COMP][curr_depth][abs_index]&2) || (et->cbf_buffs[Y_COMP][curr_depth][abs_index+num_part_in_sub_cu]&2) || (et->cbf_buffs[Y_COMP][curr_depth][abs_index+2*num_part_in_sub_cu]&2) || (et->cbf_buffs[Y_COMP][curr_depth][abs_index+3*num_part_in_sub_cu]&2);
						//consolidate cbf flags
						for(ll=abs_index;ll<abs_index+num_part_in_cu;ll++)
						{
							et->cbf_buffs[Y_COMP][curr_depth-1][ll] |= cbf_split[Y_COMP];
						}
					}

					if(curr_depth==et->max_pred_partition_depth)
					{
						//if we fill this in here we don't have to consolidate
						memset(&ctu->pred_depth[abs_index], curr_depth-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu->pred_depth[0]));
						memset(&ctu->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu->part_size_type[0]));
						if(et->rd_mode==1)
						{
							memset(&ctu_rd->pred_depth[abs_index], curr_depth-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu_rd->pred_depth[0]));
							memset(&ctu_rd->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu_rd->part_size_type[0]));
						}
					}
				}
				else
				{
					//Aqui prevalecen los resultados de la computacion top-down
					//top-down computation results are prefered
					int part_size_type2 = (curr_depth-1<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//

					//if we fill this in here we don't have to consolidate
					memset(&ctu->pred_depth[abs_index], curr_depth-1-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu->pred_depth[0]));
					memset(&ctu->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu->part_size_type[0]));
					if(et->rd_mode==1)//rd
					{
						memset(&ctu_rd->pred_depth[abs_index], curr_depth-1-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu_rd->pred_depth[0]));
						memset(&ctu_rd->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu_rd->part_size_type[0]));
					}
				}
				cost_sum[curr_depth] = 0;
				curr_depth--;
				parent_part_info = parent_part_info->parent;

				if(et->performance_mode == 0 && curr_depth>0 && curr_depth<et->max_pred_partition_depth && depth_state[curr_depth]<4 && ctu->partition_list[0].is_b_inside_frame && ctu->partition_list[0].is_r_inside_frame)
				{
					double totalcost = 0;
					int h;
					for(h=0;h<depth_state[curr_depth];h++)
						totalcost += parent_part_info->children[h]->cost;

					//if(curr_depth>0 && curr_depth<et->max_pred_partition_depth && /*depth_state[curr_depth]!=4 && totalcost > parent_part_info->cost && et->partition_info[0].is_b_inside_frame && et->partition_info[0].is_r_inside_frame)//parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
					if(totalcost > parent_part_info->cost)//cost of childs is already higher than cost of parent - stop recursion
					{
						depth_state[curr_depth] = 4;
					}
				}
				cost_chroma = 0;//only!=0 if part_size_type == SIZE_NxN
				*/
			}

/*			max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

			if(curr_depth+1 <= max_processing_depth)//el = es para cuando et->max_intra_tr_depth!=4
			{
				int aux_depth;
				cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
				abs_index = aux_partition_info->abs_index;
				num_part_in_cu  = aux_partition_info->num_part_in_cu;

				for(aux_depth=curr_depth+1;aux_depth<=max_processing_depth;aux_depth++)//+1 pq el nivel superior esta ya sincronizado y tenemos que sincronizar los siguientes
				{
					synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);	
					synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);
					//for rd
					consolidate_recursive_info_buffers(et, gcnt, parent_part_info, curr_depth, aux_depth, abs_index, num_part_in_cu);
				}
				synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
			}
//			acc_cost = 0;
*/
		}

		if(parent_part_info!=NULL)
			curr_partition_info = parent_part_info->children[depth_state[curr_depth]];
	}
	
/*	curr_partition_info = &ctu->partition_list[0];
	abs_index = curr_partition_info->abs_index;
	curr_depth = curr_partition_info->depth;
	num_part_in_cu = curr_partition_info->num_part_in_cu;

	//esta informacion o parte de ella se necesita para la codificacion de otros ctus (por lo menos el intra_mode y el pred_mode. Creo que los cbfs y los tr_idx solo se necesitan los actuales, no los de ctus anteriores)
	memcpy(&ctu->cbf[Y_COMP][abs_index], &et->cbf_buffs[Y_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->cbf[Y_COMP][0]));
	memcpy(&ctu->cbf[U_COMP][abs_index], &et->cbf_buffs[U_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->cbf[U_COMP][0]));
	memcpy(&ctu->cbf[V_COMP][abs_index], &et->cbf_buffs[V_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->cbf[V_COMP][0]));
	memcpy(&ctu->intra_mode[Y_COMP][abs_index], &et->intra_mode_buffs[Y_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->intra_mode[Y_COMP][0]));
	memcpy(&ctu->intra_mode[CHR_COMP][abs_index], &et->intra_mode_buffs[CHR_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->intra_mode[CHR_COMP][0]));
	memcpy(&ctu->tr_idx[abs_index], &et->tr_idx_buffs[curr_depth][abs_index], num_part_in_cu*sizeof(ctu->tr_idx[0]));

//	memset(&ctu->pred_mode[abs_index], INTRA_MODE, num_part_in_cu*sizeof(ctu->pred_mode[0]));//indicamos que todas las codificaciones son intra
	ctu->pred_mode = INTRA_MODE;
	return best_cost;
*/
}

