/*****************************************************************************
* hmr_sao.c : homerHEVC encoding library
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





uint g_saoMaxOffsetQVal[NUM_PICT_COMPONENTS];
uint m_offsetStepLog2[NUM_PICT_COMPONENTS];


//para alojar
sao_stat_data_t stat_data[12][NUM_PICT_COMPONENTS][NUM_SAO_NEW_TYPES];//enc_engine or ctu
sao_blk_param_t recon_params[12];//enc_engine
sao_blk_param_t coded_params[12];//enc_engine


int8_t m_signLineBuf1[64+1]; //en thread
int8_t m_signLineBuf2[64+1]; //en thread

int calculate_preblock_stats = FALSE;//esto no se si deberia ser siempre false

static int skiped_lines_r[NUM_PICT_COMPONENTS] = {5,3,3};
static int skiped_lines_b[NUM_PICT_COMPONENTS] = {4,2,2};

static int num_lcu_sao_off[NUM_PICT_COMPONENTS];//not used ¿?
static int sao_disabled_rate[NUM_PICT_COMPONENTS];//[NUM_TEMP_LAYERS];//slice

int slice_enabled[NUM_PICT_COMPONENTS];//slice



//double m_lambda[3] = {57.908390375799925, 45.961919900164979, 45.961919900164979};



void sao_init(int bit_depth)
{
	int component;
	for(component =0; component < NUM_PICT_COMPONENTS; component++)
	{
		m_offsetStepLog2  [component] = max(bit_depth - MAX_SAO_TRUNCATED_BITDEPTH, 0);
		g_saoMaxOffsetQVal[component] = (1<<(min(bit_depth,MAX_SAO_TRUNCATED_BITDEPTH)-5))-1; //Table 9-32, inclusive
	}
}


void get_ctu_stats(hvenc_engine_t* enc_engine, slice_t *currslice, ctu_info_t* ctu, sao_stat_data_t stats[][NUM_PICT_COMPONENTS][NUM_SAO_NEW_TYPES])	
{
	int component;
	int l_available, r_available, t_available, b_available, tl_available, bl_available, tr_available, br_available;
	int ctu_idx = ctu->ctu_number;
	int height_luma = (ctu->y[Y_COMP] + ctu->size > enc_engine->pict_height[Y_COMP])?(enc_engine->pict_height[Y_COMP]-ctu->y[Y_COMP]):ctu->size;
	int width_luma = (ctu->x[Y_COMP] + ctu->size > enc_engine->pict_width[Y_COMP])?(enc_engine->pict_width[Y_COMP]-ctu->x[Y_COMP]):ctu->size;

	l_available = (ctu->x[Y_COMP]> 0);
	t_available = (ctu->y[Y_COMP]> 0);
	r_available = ((ctu->x[Y_COMP] + ctu->size) < enc_engine->pict_width[Y_COMP]);
	b_available = ((ctu->y[Y_COMP] + ctu->size) < enc_engine->pict_height[Y_COMP]);
	tl_available = t_available & l_available;
	tr_available = t_available & r_available;
	bl_available = b_available & l_available;
	br_available = b_available & r_available;

	memset(&stats[ctu_idx][0][0], 0, sizeof(stats[ctu_idx]));

	if(ctu->ctu_number==6)
	{
		int iiiiii=0;
	}

	for(component=Y_COMP; component < NUM_PICT_COMPONENTS; component++)
	{
		int skip_lines_r = skiped_lines_r[component];
		int skip_lines_b = skiped_lines_b[component];
		int decoded_buff_stride = WND_STRIDE_2D(enc_engine->curr_reference_frame->img, component);
		int16_t *decoded_buff  = WND_POSITION_2D(int16_t *, enc_engine->curr_reference_frame->img, component, ctu->x[component], ctu->y[component], 0, enc_engine->ctu_width);
		int orig_buff_stride = WND_STRIDE_2D(enc_engine->current_pict.img2encode->img, component);
		uint8_t *orig_buff = WND_POSITION_2D(uint8_t *, enc_engine->current_pict.img2encode->img, component, ctu->x[component], ctu->y[component], 0, enc_engine->ctu_width);
		int type_idx;
		int chroma_shift = (component==Y_COMP)?0:1;
		int height = height_luma>>chroma_shift;
		int width = width_luma>>chroma_shift;

		//getBlkStats
		for(type_idx=0; type_idx< NUM_SAO_NEW_TYPES; type_idx++)
		{
			sao_stat_data_t *curr_stat_data = &stats[ctu_idx][component][type_idx];
			int64_t *diff = curr_stat_data->diff;
			int64_t *count = curr_stat_data->count;
			int x,y, start_x, start_y, end_x, end_y, edge_type, first_line_start_x, first_line_end_x;
			int8_t sign_left, sign_right, sign_down;
			int16_t *src_line = decoded_buff;
			int src_stride = decoded_buff_stride;
			uint8_t *org_line = orig_buff;
			int org_stride = orig_buff_stride;

			switch(type_idx)
			{
			case SAO_TYPE_EO_0:
				{
					diff +=2;
					count+=2;
					end_y   = (b_available) ? (height - skip_lines_b) : height;
					start_x = (!calculate_preblock_stats) ? (l_available ? 0 : 1): (r_available ? (width - skip_lines_r) : (width - 1));
					end_x   = (!calculate_preblock_stats) ? (r_available ? (width - skip_lines_r) : (width - 1)): (r_available ? width : (width - 1));

					for (y=0; y<end_y; y++)
					{
						int aux = src_line[start_x] - src_line[start_x-1];
						sign_left = aux==0?0:(SIGN(aux));//(Char)m_sign[srcLine[startX] - srcLine[startX-1]];
						for (x=start_x; x<end_x; x++)
						{
							int aux = src_line[x] - src_line[x+1];
							sign_right =  aux==0?0:(SIGN(aux));//(Char)m_sign[srcLine[x] - srcLine[x+1]]; 
							edge_type  =  sign_right + sign_left;
							sign_left  = -sign_right;

							diff [edge_type] += (org_line[x] - src_line[x]);
							count[edge_type] ++;
						}
						src_line  += src_stride;
						org_line  += org_stride;
					}				
				}
				break;
			case SAO_TYPE_EO_90:
				{
					int8_t *signUpLine = m_signLineBuf1;
					int16_t* srcLineAbove, *srcLineBelow;
					diff +=2;
					count+=2;
					
					start_x = (!calculate_preblock_stats)?0:(r_available ? (width - skip_lines_r) : width);
					start_y = t_available ? 0 : 1;
					end_x   = (!calculate_preblock_stats) ? (r_available ? (width - skip_lines_r) : width): width;
					end_y   = b_available ? (height - skip_lines_b) : (height - 1);
					if (!t_available)
					{
						src_line  += src_stride;
						org_line  += org_stride;
					}

					srcLineAbove = src_line - src_stride;
					for (x=start_x; x<end_x; x++) 
					{
						int aux = src_line[x] - srcLineAbove[x];
						signUpLine[x] = aux==0?0:(SIGN(aux));//(Char)m_sign[srcLine[x] - srcLineAbove[x]];
					}

					for (y=start_y; y<end_y; y++)
					{
						srcLineBelow = src_line + src_stride;

						for (x=start_x; x<end_x; x++)
						{
							int aux = src_line[x] - srcLineBelow[x];
							sign_down  = aux==0?0:(SIGN(aux));//(Char)m_sign[srcLine[x] - srcLineBelow[x]]; 
							edge_type  = sign_down + signUpLine[x];
							signUpLine[x]= -sign_down;

							diff [edge_type] += (org_line[x] - src_line[x]);
							count[edge_type] ++;
						}
						src_line += src_stride;
						org_line += org_stride;
					}

				}
				break;
			case SAO_TYPE_EO_135:
				{
					int8_t *signUpLine, *signDownLine, *signTmpLine;
					int16_t* srcLineBelow, *srcLineAbove;
					diff +=2;
					count+=2;
					
					signUpLine = m_signLineBuf1;
					signDownLine = m_signLineBuf2;

					start_x = (!calculate_preblock_stats) ? (l_available ? 0 : 1) : (r_available ? (width - skip_lines_r) : (width - 1));

					end_x = (!calculate_preblock_stats) ? (r_available ? (width - skip_lines_r): (width - 1)) : (r_available ? width : (width - 1));
					end_y = b_available ? (height - skip_lines_b) : (height - 1);

					//prepare 2nd line's upper sign
					srcLineBelow = src_line + src_stride;
					for (x=start_x; x<end_x+1; x++)
					{
						int aux = srcLineBelow[x] - src_line[x-1];
						signUpLine[x] = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLineBelow[x] - srcLine[x-1]];
					}

					//1st line
					srcLineAbove = src_line - src_stride;
					first_line_start_x = (!calculate_preblock_stats) ? (tl_available ? 0 : 1) : start_x;
					first_line_end_x   = (!calculate_preblock_stats) ? (t_available ? end_x : 1) : end_x;
					for(x=first_line_start_x; x<first_line_end_x; x++)
					{
						int aux = src_line[x] - srcLineAbove[x-1];
						edge_type = ((aux==0)?0:(SIGN(aux))) - signUpLine[x+1];//m_sign[srcLine[x] - srcLineAbove[x-1]] - signUpLine[x+1];

						diff [edge_type] += (org_line[x] - src_line[x]);
						count[edge_type] ++;
					}
					src_line  += src_stride;
					org_line  += org_stride;

					//middle lines
					for (y=1; y<end_y; y++)
					{
						int aux;
						srcLineBelow = src_line + src_stride;

						for (x=start_x; x<end_x; x++)
						{
							int aux0 = src_line[x] - srcLineBelow[x+1];
							sign_down = (aux0==0)?0:(SIGN(aux0));//(Char)m_sign[srcLine[x] - srcLineBelow[x+1]] ;
							
							edge_type = sign_down + signUpLine[x];
							diff [edge_type] += (org_line[x] - src_line[x]);
							count[edge_type] ++;

							signDownLine[x+1] = -sign_down; 
						}
						aux = srcLineBelow[start_x] - src_line[start_x-1];
						signDownLine[start_x] = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLineBelow[startX] - srcLine[startX-1]];

						signTmpLine  = signUpLine;
						signUpLine   = signDownLine;
						signDownLine = signTmpLine;

						src_line += src_stride;
						org_line += org_stride;
					}

				}
				break;
			case SAO_TYPE_EO_45:
				{
					int8_t *signUpLine = m_signLineBuf1+1;
					int16_t *srcLineBelow, *srcLineAbove;
					diff +=2;
					count+=2;

					start_x = (!calculate_preblock_stats) ? (l_available ? 0 : 1) : (r_available ? (width - skip_lines_r) : (width - 1));
					end_x = (!calculate_preblock_stats) ? (r_available ? (width - skip_lines_r) : (width - 1)) : (r_available ? width : (width - 1));
					end_y = b_available ? (height - skip_lines_b) : (height - 1);

					//prepare 2nd line upper sign
					srcLineBelow = src_line + src_stride;
					for (x=start_x-1; x<end_x; x++)
					{
						int aux = srcLineBelow[x] - src_line[x+1];
						signUpLine[x] = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLineBelow[x] - srcLine[x+1]];
					}

					//first line
					srcLineAbove = src_line - src_stride;
					first_line_start_x = (!calculate_preblock_stats) ? (t_available ? start_x : end_x) : start_x;
					first_line_end_x   = (!calculate_preblock_stats) ? ((!r_available && tr_available) ? width : end_x) : end_x;
					for(x=first_line_start_x; x<first_line_end_x; x++)
					{
						int aux = src_line[x] - srcLineAbove[x+1];
						edge_type = ((aux==0)?0:(SIGN(aux))) - signUpLine[x-1];//m_sign[srcLine[x] - srcLineAbove[x+1]] - signUpLine[x-1];
						diff [edge_type] += (org_line[x] - src_line[x]);
						count[edge_type] ++;
					}

					src_line += src_stride;
					org_line += org_stride;

					//middle lines
					for (y=1; y<end_y; y++)
					{
						int aux;
						srcLineBelow = src_line + src_stride;

						for(x=start_x; x<end_x; x++)
						{
							int aux0 = src_line[x] - srcLineBelow[x-1];
							sign_down = (aux0==0)?0:(SIGN(aux0));//(Char)m_sign[srcLine[x] - srcLineBelow[x-1]] ;
							edge_type = sign_down + signUpLine[x];

							diff [edge_type] += (org_line[x] - src_line[x]);
							count[edge_type]++;

							signUpLine[x-1] = -sign_down; 
						}
						aux = srcLineBelow[end_x-1] - src_line[end_x];
						signUpLine[end_x-1] = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLineBelow[endX-1] - srcLine[endX]];
						src_line  += src_stride;
						org_line  += org_stride;
					}
				}
				break;
			case SAO_TYPE_BO:
				{
					int shiftBits = enc_engine->bit_depth - NUM_SAO_BO_CLASSES_LOG2;
					start_x = (!calculate_preblock_stats) ? 0:( r_available?(width- skip_lines_r):width);
					end_x = (!calculate_preblock_stats) ? (r_available ? (width - skip_lines_r) : (width)) : width;
					end_y = b_available ? (height - skip_lines_b) : height;
					
					for (y=0; y< end_y; y++)
					{
						for (x=start_x; x< end_x; x++)
						{
							int bandIdx= src_line[x] >> shiftBits; 
							diff [bandIdx] += (org_line[x] - src_line[x]);
							count[bandIdx] ++;
						}
						src_line += src_stride;
						org_line += org_stride;
					}
				}
				break;
			default:
				{
					printf("Not a supported SAO types\n");
//					assert(0);
					exit(-1);
				}
			}

		}

	}
}

void decide_pic_params(int *slice_enable)	
{
	int component;
	for(component=Y_COMP; component < NUM_PICT_COMPONENTS; component++)
	{		// reset flags & counters
		slice_enable[component] = TRUE;
	}
}


int get_merge_list(hvenc_engine_t* enc_engine, int ctu_idx, sao_blk_param_t *blk_params, sao_blk_param_t* merge_list[])
{
	int ctuX = ctu_idx % enc_engine->pict_width_in_ctu;
	int ctuY = ctu_idx / enc_engine->pict_width_in_ctu;
	int mergedCTUPos;
	int numValidMergeCandidates = 0;
	int merge_type;

	for(merge_type=0; merge_type< NUM_SAO_MERGE_TYPES; merge_type++)
	{
		sao_blk_param_t *mergeCandidate = NULL;

		switch(merge_type)
		{
		case SAO_MERGE_ABOVE:
			{
				if(ctuY > 0)
				{
					mergedCTUPos = ctu_idx - enc_engine->pict_width_in_ctu;
//					if( pic->getSAOMergeAvailability(ctu_idx, mergedCTUPos) )
					{
						mergeCandidate = &(blk_params[mergedCTUPos]);
					}
				}
			}
			break;
		case SAO_MERGE_LEFT:
			{
				if(ctuX > 0)
				{
					mergedCTUPos = ctu_idx - 1;
//					if( pic->getSAOMergeAvailability(ctu_idx, mergedCTUPos) )
					{
						mergeCandidate = &(blk_params[mergedCTUPos]);
					}
				}
			}
			break;
		default:
			{
				printf("not a supported merge type");
//				//assert(0);
				exit(-1);
			}
		}

//		mergeList.push_back(mergeCandidate);
		if (mergeCandidate != NULL)
		{
			merge_list[numValidMergeCandidates] = mergeCandidate;
			numValidMergeCandidates++;
		}
	}

	return numValidMergeCandidates;

}

__inline double xRoundIbdi2(int bit_depth, double x)
{
	return ((x)>0) ? (int)(((int)(x)+(1<<(bit_depth-8-1)))/(1<<(bit_depth-8))) : ((int)(((int)(x)-(1<<(bit_depth-8-1)))/(1<<(bit_depth-8))));
}

__inline double x_round_ibdi(int bit_depth, double x)
{
	return (bit_depth > 8 ? xRoundIbdi2(bit_depth, (x)) : ((x)>=0 ? ((int)((x)+0.5)) : ((int)((x)-0.5)))) ;
}


//inline Int64 TEncSampleAdaptiveOffset::estSaoDist(Int64 count, Int64 offset, Int64 diffSum, Int shift)
int64_t estSaoDist(int64_t count, int64_t offset, int64_t diffSum, int64_t shift)
{
	return (( count*offset*offset-diffSum*offset*2 ) >> shift);
}

//inline int TEncSampleAdaptiveOffset::estIterOffset(int type_idx, int class_idx, Double lambda, int offset_input, Int64 count, Int64 diff_sum, int shift, int bit_increase, Int64& best_dist, Double& bestCost, int offset_thrshl )
int est_iter_offset(int type_idx, int class_idx, double lambda, int offset_input, int64_t count, int64_t diff_sum, int shift, int bit_increase, int64_t *best_dist, double* best_cost, int offset_thrshl)
{
	int iterOffset, tempOffset;
	int64_t tempDist, tempRate;
	double tempCost, tempMinCost;
	int offsetOutput = 0;
	iterOffset = offset_input;
	// Assuming sending quantized value 0 results in zero offset and sending the value zero needs 1 bit. entropy coder can be used to measure the exact rate here. 
	tempMinCost = lambda; 
	while (iterOffset != 0)
	{
		// Calculate the bits required for signaling the offset
		tempRate = (type_idx == SAO_TYPE_BO) ? (abs((int)iterOffset)+2) : (abs((int)iterOffset)+1); 
		if (abs((int)iterOffset)==offset_thrshl) //inclusive 
		{  
			tempRate --;
		}
		// Do the dequantization before distortion calculation
		tempOffset  = iterOffset << bit_increase;
		tempDist    = estSaoDist( count, tempOffset, diff_sum, shift);
		tempCost    = ((double)tempDist + lambda * (double) tempRate);
		if(tempCost < tempMinCost)
		{
			tempMinCost = tempCost;
			offsetOutput = iterOffset;
			*best_dist = tempDist;
			*best_cost = tempCost;
		}
		iterOffset = (iterOffset > 0) ? (iterOffset-1):(iterOffset+1);
	}
	return offsetOutput;
}


//Void TEncSampleAdaptiveOffset::deriveOffsets(int ctu, int compIdx, int typeIdc, SAOStatData& statData, int* quant_offsets, int& typeAuxInfo)
void derive_offsets(hvenc_engine_t *enc_engine, int component, int type_idc, sao_stat_data_t *stats, int* quant_offsets, int* type_aux_info)
{
	int bit_depth = enc_engine->bit_depth;//(component== Y_COMP) ? g_bitDepthY : g_bitDepthC;
	int shift = 2 * DISTORTION_PRECISION_ADJUSTMENT(bit_depth-8);
	int offset_thrshld = g_saoMaxOffsetQVal[component];  //inclusive
	int num_classes, class_idx;
	double *lambdas = enc_engine->lambdas;

	memset(quant_offsets, 0, sizeof(int)*MAX_NUM_SAO_CLASSES);

	//derive initial offsets 
	num_classes = (type_idc == SAO_TYPE_BO)?((int)NUM_SAO_BO_CLASSES):((int)NUM_SAO_EO_CLASSES);
	for(class_idx=0; class_idx< num_classes; class_idx++)
	{
		if( (type_idc != SAO_TYPE_BO) && (class_idx==SAO_CLASS_EO_PLAIN)  ) 
		{
			continue; //offset will be zero
		}

		if(stats->count[class_idx] == 0)
		{
			continue; //offset will be zero
		}

		quant_offsets[class_idx] = (int) x_round_ibdi(bit_depth, (double)( stats->diff[class_idx]<<(bit_depth-8))/(double)(stats->count[class_idx]<< m_offsetStepLog2[component]));

		quant_offsets[class_idx] = clip(quant_offsets[class_idx], -offset_thrshld, offset_thrshld);
	}

	// adjust offsets
	switch(type_idc)
	{
	case SAO_TYPE_EO_0:
	case SAO_TYPE_EO_90:
	case SAO_TYPE_EO_135:
	case SAO_TYPE_EO_45:
		{
			int64_t class_dist;
			double class_cost;
			int class_idx;
			for(class_idx=0; class_idx<NUM_SAO_EO_CLASSES; class_idx++)  
			{         
				if(class_idx==SAO_CLASS_EO_FULL_VALLEY && quant_offsets[class_idx] < 0) quant_offsets[class_idx] =0;
				if(class_idx==SAO_CLASS_EO_HALF_VALLEY && quant_offsets[class_idx] < 0) quant_offsets[class_idx] =0;
				if(class_idx==SAO_CLASS_EO_HALF_PEAK   && quant_offsets[class_idx] > 0) quant_offsets[class_idx] =0;
				if(class_idx==SAO_CLASS_EO_FULL_PEAK   && quant_offsets[class_idx] > 0) quant_offsets[class_idx] =0;

				if( quant_offsets[class_idx] != 0 ) //iterative adjustment only when derived offset is not zero
				{
					quant_offsets[class_idx] = est_iter_offset( type_idc, class_idx, lambdas[component], quant_offsets[class_idx], stats->count[class_idx], stats->diff[class_idx], shift, m_offsetStepLog2[component], &class_dist , &class_cost , offset_thrshld);
				}
			}

			*type_aux_info = 0;
		}
		break;
	case SAO_TYPE_BO:
		{
			int64_t dist_bo_classes[NUM_SAO_BO_CLASSES];
			double cost_bo_classes[NUM_SAO_BO_CLASSES];
			int clear_quant_offset[NUM_SAO_BO_CLASSES];
			double min_cost, cost;
			int band, i;
			int class_idx;
			memset(dist_bo_classes, 0, sizeof(dist_bo_classes));
			for(class_idx=0; class_idx< NUM_SAO_BO_CLASSES; class_idx++)
			{         
				cost_bo_classes[class_idx]= lambdas[component];
				if( quant_offsets[class_idx] != 0 ) //iterative adjustment only when derived offset is not zero
				{
					quant_offsets[class_idx] = est_iter_offset( type_idc, class_idx, lambdas[component], quant_offsets[class_idx], stats->count[class_idx], stats->diff[class_idx], shift, m_offsetStepLog2[component], &dist_bo_classes[class_idx], &cost_bo_classes[class_idx], offset_thrshld);
				}
			}

			//decide the starting band index
			min_cost = MAX_COST;
			for(band=0; band< NUM_SAO_BO_CLASSES- 4+ 1; band++) 
			{
				cost  = cost_bo_classes[band  ];
				cost += cost_bo_classes[band+1];
				cost += cost_bo_classes[band+2];
				cost += cost_bo_classes[band+3];

				if(cost < min_cost)
				{
					min_cost = cost;
					*type_aux_info = band;
				}
			}
			//clear those unused classes

			memset(clear_quant_offset, 0, sizeof(int)*NUM_SAO_BO_CLASSES);
			for(i=0; i< 4; i++) 
			{
				int band = (*type_aux_info+i)%NUM_SAO_BO_CLASSES;
				clear_quant_offset[band] = quant_offsets[band];
			}
			memcpy(quant_offsets, clear_quant_offset, sizeof(clear_quant_offset));        
		}
		break;
	default:
		{
			printf("Not a supported type");
			//assert(0);
			exit(-1);
		}

	}
}


//Void TComSampleAdaptiveOffset::invertQuantOffsets(Int compIdx, Int typeIdc, Int typeAuxInfo, Int* dstOffsets, Int* srcOffsets)
void invert_quant_offsets(int component, int type_idc, int typeAuxInfo, int* dstOffsets, int* srcOffsets)
{
  int codedOffset[MAX_NUM_SAO_CLASSES];

  memcpy(codedOffset, srcOffsets, sizeof(int)*MAX_NUM_SAO_CLASSES);
  memset(dstOffsets, 0, sizeof(int)*MAX_NUM_SAO_CLASSES);

  if(type_idc == SAO_TYPE_START_BO)
  {
	int i;
    for(i=0; i< 4; i++)
    {
      dstOffsets[(typeAuxInfo+ i)%NUM_SAO_BO_CLASSES] = codedOffset[(typeAuxInfo+ i)%NUM_SAO_BO_CLASSES]*(1<<m_offsetStepLog2[component]);
    }
  }
  else //EO
  {
	int i;
    for(i=0; i< NUM_SAO_EO_CLASSES; i++)
    {
      dstOffsets[i] = codedOffset[i] *(1<<m_offsetStepLog2[component]);
    }
//    assert(dstOffsets[SAO_CLASS_EO_PLAIN] == 0); //keep EO plain offset as zero
  }

}

//inline Int64 TEncSampleAdaptiveOffset::estSaoDist(Int64 count, Int64 offset, Int64 diffSum, Int shift)
__inline int64_t est_sao_dist(int64_t count, int64_t offset, int64_t diffSum, int shift)
{
  return (( count*offset*offset-diffSum*offset*2 ) >> shift);
}


//Int64 TEncSampleAdaptiveOffset::getDistortion(Int ctu, Int compIdx, Int typeIdc, Int typeAuxInfo, Int* invQuantOffset, SAOStatData& statData)
int64_t get_distortion(int typeIdc, int typeAuxInfo, int* invQuantOffset, sao_stat_data_t *stats, int bit_depth)
{
  int64_t dist=0;
  int inputBitDepth = bit_depth;//(component == Y_COMP) ? g_bitDepthY : g_bitDepthC ;
  int shift = 2 * DISTORTION_PRECISION_ADJUSTMENT(inputBitDepth-8);

  switch(typeIdc)
  {
    case SAO_TYPE_EO_0:
    case SAO_TYPE_EO_90:
    case SAO_TYPE_EO_135:
    case SAO_TYPE_EO_45:
      {
		int offsetIdx;
        for (offsetIdx=0; offsetIdx<NUM_SAO_EO_CLASSES; offsetIdx++)
        {
          dist += estSaoDist( stats->count[offsetIdx], invQuantOffset[offsetIdx], stats->diff[offsetIdx], shift);
        }        
      }
      break;
    case SAO_TYPE_BO:
      {
		  int offsetIdx;
        for (offsetIdx=typeAuxInfo; offsetIdx<typeAuxInfo+4; offsetIdx++)
        {
          int bandIdx = offsetIdx % NUM_SAO_BO_CLASSES ; 
          dist += estSaoDist( stats->count[bandIdx], invQuantOffset[bandIdx], stats->diff[bandIdx], shift);
        }
      }
      break;
    default:
      {
        printf("Not a supported type");
//        assert(0);
        exit(-1);
      }
  }

  return dist;
}


//Void TEncSampleAdaptiveOffset::deriveModeNewRDO(int ctu, std::vector<SAOBlkParam*>& mergeList, Bool* sliceEnabled, SAOStatData*** blkStats, SAOBlkParam& mode_param, Double& modeNormCost, TEncSbac** cabacCoderRDO, int inCabacLabel)
void derive_mode_new_rdo(hvenc_engine_t* enc_engine, sao_stat_data_t stats[NUM_PICT_COMPONENTS][NUM_SAO_NEW_TYPES], sao_blk_param_t *mode_param, double *mode_cost, int slice_enabled[] )
{
	double minCost, cost;
	int rate;
	uint previousWrittenBits;
	int64_t dist[NUM_PICT_COMPONENTS], modeDist[NUM_PICT_COMPONENTS];
	sao_offset_t testOffset[NUM_PICT_COMPONENTS];
	int component;
	int invQuantOffset[MAX_NUM_SAO_CLASSES];
	int type_idc;

	modeDist[Y_COMP]= modeDist[U_COMP] = modeDist[V_COMP] = 0;

	//pre-encode merge flags
	mode_param->offsetParam[Y_COMP].modeIdc = SAO_MODE_OFF;
//	m_pcRDGoOnSbacCoder->load(cabacCoderRDO[inCabacLabel]);
//	m_pcRDGoOnSbacCoder->codeSAOBlkParam(mode_param, sliceEnabled, (mergeList[SAO_MERGE_LEFT]!= NULL), (mergeList[SAO_MERGE_ABOVE]!= NULL), true);
//	m_pcRDGoOnSbacCoder->store(cabacCoderRDO[SAO_CABACSTATE_BLK_MID]);

	//------ luma --------//
	component = Y_COMP;
	//"off" case as initial cost
	mode_param->offsetParam[component].modeIdc = SAO_MODE_OFF;
//	m_pcRDGoOnSbacCoder->resetBits();
//	m_pcRDGoOnSbacCoder->codeSAOOffsetParam(component, mode_param[component], sliceEnabled[component]);
	modeDist[component] = 0;
	minCost = MAX_COST;
//	minCost= m_lambda[component]*((Double)m_pcRDGoOnSbacCoder->getNumberOfWrittenBits());
//	m_pcRDGoOnSbacCoder->store(cabacCoderRDO[SAO_CABACSTATE_BLK_TEMP]);
	if(slice_enabled[component])
	{
		int type_idc;
		for(type_idc=0; type_idc< NUM_SAO_NEW_TYPES; type_idc++)
		{
			testOffset[component].modeIdc = SAO_MODE_NEW;
			testOffset[component].typeIdc = type_idc;

			//derive coded offset
//			deriveOffsets(ctu, component, type_idc, blkStats[ctu][component][type_idc], testOffset[component].offset, testOffset[component].typeAuxInfo);
			
			derive_offsets(enc_engine, component, type_idc, &stats[component][type_idc], testOffset[component].offset, &testOffset[component].typeAuxInfo);

			//inversed quantized offsets
			//invertQuantOffsets(component, type_idc, testOffset[component].typeAuxInfo, invQuantOffset, testOffset[component].offset);
			invert_quant_offsets(component, type_idc, testOffset[component].typeAuxInfo, invQuantOffset, testOffset[component].offset);

			//get distortion
			//dist[component] = getDistortion(ctu, component, testOffset[component].type_idc, testOffset[component].typeAuxInfo, invQuantOffset, blkStats[ctu][component][type_idc]);
			dist[component] = get_distortion(testOffset[component].typeIdc, testOffset[component].typeAuxInfo, invQuantOffset, &stats[component][type_idc], enc_engine->bit_depth);
//				int64_t get_distortion(int typeIdc, int typeAuxInfo, int* invQuantOffset, sao_stat_data_t *stats, int bit_depth)
			//get rate
//			m_pcRDGoOnSbacCoder->load(cabacCoderRDO[SAO_CABACSTATE_BLK_MID]);
//			m_pcRDGoOnSbacCoder->resetBits();
//			m_pcRDGoOnSbacCoder->codeSAOOffsetParam(component, testOffset[component], sliceEnabled[component]);
//			rate = m_pcRDGoOnSbacCoder->getNumberOfWrittenBits();
			cost = (double)dist[component];// + m_lambda[component]*((Double)rate);
			if(cost < minCost)
			{
				minCost = cost;
				modeDist[component] = dist[component];
				mode_param->offsetParam[component]= testOffset[component];
//				m_pcRDGoOnSbacCoder->store(cabacCoderRDO[SAO_CABACSTATE_BLK_TEMP]);
			}
		}
	}
//	m_pcRDGoOnSbacCoder->load(cabacCoderRDO[SAO_CABACSTATE_BLK_TEMP]);
//	m_pcRDGoOnSbacCoder->store(cabacCoderRDO[SAO_CABACSTATE_BLK_MID]);

	//------ chroma --------//
	//"off" case as initial cost
	cost = 0;
	previousWrittenBits = 0;
//	m_pcRDGoOnSbacCoder->resetBits();
	for (component = U_COMP; component < NUM_PICT_COMPONENTS; component++)
	{
		mode_param->offsetParam[component].modeIdc = SAO_MODE_OFF; 
		modeDist [component] = 0;

//		m_pcRDGoOnSbacCoder->codeSAOOffsetParam(component, mode_param[component], sliceEnabled[component]);

//		const UInt currentWrittenBits = m_pcRDGoOnSbacCoder->getNumberOfWrittenBits();
//		cost += m_lambda[component] * (currentWrittenBits - previousWrittenBits);
//		previousWrittenBits = currentWrittenBits;
	}

	minCost = cost;
	minCost = MAX_COST;
	//doesn't need to store cabac status here since the whole CTU parameters will be re-encoded at the end of this function

	for(type_idc=0; type_idc< NUM_SAO_NEW_TYPES; type_idc++)
	{
//		m_pcRDGoOnSbacCoder->load(cabacCoderRDO[SAO_CABACSTATE_BLK_MID]);
//		m_pcRDGoOnSbacCoder->resetBits();
		previousWrittenBits = 0;
		cost = 0;

		for(component= U_COMP; component< NUM_PICT_COMPONENTS; component++)
		{
			if(!slice_enabled[component])
			{
				testOffset[component].modeIdc = SAO_MODE_OFF;
				dist[component]= 0;
				continue;
			}
			testOffset[component].modeIdc = SAO_MODE_NEW;
			testOffset[component].typeIdc = type_idc;

			//derive offset & get distortion
//			deriveOffsets(ctu, component, typeIdc, blkStats[ctu][component][typeIdc], testOffset[component].offset, testOffset[component].typeAuxInfo);
			derive_offsets(enc_engine, component, type_idc, &stats[component][type_idc], testOffset[component].offset, &testOffset[component].typeAuxInfo);

//			invertQuantOffsets(component, typeIdc, testOffset[component].typeAuxInfo, invQuantOffset, testOffset[component].offset);
			invert_quant_offsets(component, type_idc, testOffset[component].typeAuxInfo, invQuantOffset, testOffset[component].offset);

//			dist[component]= getDistortion(ctu, component, typeIdc, testOffset[component].typeAuxInfo, invQuantOffset, blkStats[ctu][component][typeIdc]);
			dist[component] = get_distortion(testOffset[component].typeIdc, testOffset[component].typeAuxInfo, invQuantOffset, &stats[component][type_idc], enc_engine->bit_depth);

//			m_pcRDGoOnSbacCoder->codeSAOOffsetParam(component, testOffset[component], sliceEnabled[component]);

//			const UInt currentWrittenBits = m_pcRDGoOnSbacCoder->getNumberOfWrittenBits();
			cost += dist[component];// + (m_lambda[component] * (currentWrittenBits - previousWrittenBits));
//			previousWrittenBits = currentWrittenBits;
		}

		if(cost < minCost)
		{
			minCost = cost;
			for(component= U_COMP; component< NUM_PICT_COMPONENTS; component++)
			{
				modeDist [component] = dist      [component];
				mode_param->offsetParam[component] = testOffset[component];
			}
		}
	}

	//----- re-gen rate & normalized cost----//
	*mode_cost = (double)modeDist[Y_COMP] + (double)modeDist[U_COMP] + (double)modeDist[V_COMP];
//	for(component = Y_COMP; component < NUM_PICT_COMPONENTS; component++)
//	{
//		mode_cost += (double)modeDist[component];// / m_lambda[component];
//	}
//	m_pcRDGoOnSbacCoder->load(cabacCoderRDO[inCabacLabel]);
//	m_pcRDGoOnSbacCoder->resetBits();
//	m_pcRDGoOnSbacCoder->codeSAOBlkParam(mode_param, sliceEnabled, (mergeList[SAO_MERGE_LEFT]!= NULL), (mergeList[SAO_MERGE_ABOVE]!= NULL), false);
//	modeNormCost += (Double)m_pcRDGoOnSbacCoder->getNumberOfWrittenBits();

}

//Void TEncSampleAdaptiveOffset::deriveModeMergeRDO(Int ctu, std::vector<SAOBlkParam*>& mergeList, Bool* sliceEnabled, SAOStatData*** blkStats, SAOBlkParam& modeParam, Double& modeNormCost, TEncSbac** cabacCoderRDO, Int inCabacLabel)
void derive_mode_merge_rdo(hvenc_engine_t *enc_engine, sao_blk_param_t** merge_list, int merge_list_size, int* slice_enable, sao_stat_data_t stats[NUM_PICT_COMPONENTS][NUM_SAO_NEW_TYPES], sao_blk_param_t* mode_param, double *mode_cost)//, TEncSbac** cabacCoderRDO, Int inCabacLabel)
{
  double cost;
  sao_blk_param_t test_blk_param;
  int merge_type;
    double norm_dist=0;
	int component;
//  int mergeListSize = (Int)mergeList.size();
  *mode_cost = MAX_COST;

  for(merge_type=0; merge_type< merge_list_size; merge_type++)
  {
    if(merge_list[merge_type] == NULL)
    {
      continue;
    }

    test_blk_param = *(merge_list[merge_type]);
    //normalized distortion

    for(component=0; component< NUM_PICT_COMPONENTS; component++)
    {
		sao_offset_t *merged_offsetparam;
		test_blk_param.offsetParam[component].modeIdc = SAO_MODE_MERGE;
		test_blk_param.offsetParam[component].typeIdc = merge_type;

		merged_offsetparam = &(*(merge_list[merge_type])).offsetParam[component];

      if( merged_offsetparam->modeIdc != SAO_MODE_OFF)
      {
        //offsets have been reconstructed. Don't call inversed quantization function.
		  //get_distortion(testOffset[component].typeIdc, testOffset[component].typeAuxInfo, invQuantOffset, &stats[component][type_idc], enc_engine->bit_depth);        
		  norm_dist += (((double)get_distortion(merged_offsetparam->typeIdc, merged_offsetparam->typeAuxInfo, merged_offsetparam->offset, &stats[component][merged_offsetparam->typeIdc], enc_engine->bit_depth)))/enc_engine->lambdas[component];///m_lambda[component]);
      }
    }

    //rate
//    m_pcRDGoOnSbacCoder->load(cabacCoderRDO[inCabacLabel]);
//    m_pcRDGoOnSbacCoder->resetBits();
//    m_pcRDGoOnSbacCoder->codeSAOBlkParam(testBlkParam, sliceEnabled, (mergeList[SAO_MERGE_LEFT]!= NULL), (mergeList[SAO_MERGE_ABOVE]!= NULL), false);
//    Int rate = m_pcRDGoOnSbacCoder->getNumberOfWrittenBits();

    cost = norm_dist;//+(Double)rate;

    if(cost < *mode_cost)
    {
      *mode_cost = cost;
      *mode_param    = test_blk_param;
//      m_pcRDGoOnSbacCoder->store(cabacCoderRDO[SAO_CABACSTATE_BLK_TEMP]);
    }
  }

//  m_pcRDGoOnSbacCoder->load(cabacCoderRDO[SAO_CABACSTATE_BLK_TEMP]);


}



//Void TComSampleAdaptiveOffset::reconstructBlkSAOParam(SAOBlkParam& recParam, std::vector<SAOBlkParam*>& mergeList)
void reconstruct_blk_sao_param(sao_blk_param_t *rec_param, sao_blk_param_t* merge_list[], int merge_list_size)
{
	int component;
  for(component=0; component< NUM_PICT_COMPONENTS; component++)
  {
    sao_offset_t *offset_param = &rec_param->offsetParam[component];

    if(offset_param->modeIdc == SAO_MODE_OFF)
    {
      continue;
    }

    switch(offset_param->modeIdc)
    {
    case SAO_MODE_NEW:
      {
//        invertQuantOffsets(component, offsetParam.typeIdc, offsetParam.typeAuxInfo, offsetParam.offset, offsetParam.offset);
		  invert_quant_offsets(component, offset_param->typeIdc, offset_param->typeAuxInfo, offset_param->offset, offset_param->offset);
      }
      break;
    case SAO_MODE_MERGE:
      {
        sao_blk_param_t* merge_target = merge_list[offset_param->typeIdc];
//        assert(mergeTarget != NULL);

        offset_param = &(*merge_target).offsetParam[component];
      }
      break;
    default:
      {
        printf("Not a supported mode");
//        assert(0);
        exit(-1);
      }
    }
  }
}

//Void TComSampleAdaptiveOffset::offsetBlock(Int component, Int typeIdx, Int* offset  
//										   , Pel* srcBlk, Pel* resBlk, Int srcStride, Int resStride,  Int width, Int height
//										   , Bool isLeftAvail,  Bool isRightAvail, Bool isAboveAvail, Bool isBelowAvail, Bool isAboveLeftAvail, Bool isAboveRightAvail, Bool isBelowLeftAvail, Bool isBelowRightAvail)
void offset_block(int compIdx, int typeIdx, int* offset, int16_t* srcBlk, int16_t* resBlk, int srcStride, int resStride,  int width, int height,
				int isLeftAvail,  int isRightAvail, int isAboveAvail, int isBelowAvail, int isAboveLeftAvail, int isAboveRightAvail, int isBelowLeftAvail, int isBelowRightAvail, int bit_depth)
{

/*	if(m_lineBufWidth != m_maxCUWidth)
	{
		m_lineBufWidth = m_maxCUWidth;

		if (m_signLineBuf1) delete[] m_signLineBuf1; m_signLineBuf1 = NULL;
		m_signLineBuf1 = new Char[m_lineBufWidth+1];

		if (m_signLineBuf2) delete[] m_signLineBuf2; m_signLineBuf2 = NULL;
		m_signLineBuf2 = new Char[m_lineBufWidth+1];
	}
*/
	//int* offsetClip = m_offsetClip[compIdx];

	int x,y, startX, startY, endX, endY, edgeType;
	int firstLineStartX, firstLineEndX, lastLineStartX, lastLineEndX;
	char signLeft, signRight, signDown;

	int16_t* srcLine = srcBlk;
	int16_t* resLine = resBlk;

	switch(typeIdx)
	{
	case SAO_TYPE_EO_0:
		{
			offset += 2;
			startX = isLeftAvail ? 0 : 1;
			endX   = isRightAvail ? width : (width -1);
			for (y=0; y< height; y++)
			{
				int aux = srcLine[startX] - srcLine[startX-1];
				signLeft = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLine[startX] - srcLine[startX-1]];
				for (x=startX; x< endX; x++)
				{
					aux = srcLine[x] - srcLine[x+1];
					signRight = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLine[x] - srcLine[x+1]]; 
					edgeType =  signRight + signLeft;
					signLeft  = -signRight;

					resLine[x] = clip(srcLine[x] + offset[edgeType], 0, 255);//offsetClip[srcLine[x] + offset[edgeType]];
				}
				srcLine  += srcStride;
				resLine += resStride;
			}

		}
		break;
	case SAO_TYPE_EO_90:
		{
			int8_t* signUpLine = m_signLineBuf1;
			int16_t* srcLineAbove, *srcLineBelow;
			offset += 2;

			startY = isAboveAvail ? 0 : 1;
			endY   = isBelowAvail ? height : height-1;
			if (!isAboveAvail)
			{
				srcLine += srcStride;
				resLine += resStride;
			}

			srcLineAbove = srcLine- srcStride;
			for (x=0; x< width; x++)
			{
				int aux = srcLine[x] - srcLineAbove[x];
				signUpLine[x] = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLine[x] - srcLineAbove[x]];
			}

			
			for (y=startY; y<endY; y++)
			{
				srcLineBelow = srcLine+ srcStride;

				for (x=0; x< width; x++)
				{
					int aux = srcLine[x] - srcLineBelow[x];
					signDown  = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLine[x] - srcLineBelow[x]]; 
					edgeType = signDown + signUpLine[x];
					signUpLine[x]= -signDown;

					resLine[x] = clip(srcLine[x] + offset[edgeType], 0, 255);//offsetClip[srcLine[x] + offset[edgeType]];
				}
				srcLine += srcStride;
				resLine += resStride;
			}

		}
		break;
	case SAO_TYPE_EO_135:
		{
			int8_t *signUpLine, *signDownLine, *signTmpLine;
			int16_t *srcLineBelow, *srcLineAbove;

			offset += 2;
			signUpLine  = m_signLineBuf1;
			signDownLine= m_signLineBuf2;

			startX = isLeftAvail ? 0 : 1 ;
			endX   = isRightAvail ? width : (width-1);

			//prepare 2nd line's upper sign
			srcLineBelow= srcLine+ srcStride;
			for (x=startX; x< endX+1; x++)
			{
				int aux = srcLineBelow[x] - srcLine[x- 1];
				signUpLine[x] = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLineBelow[x] - srcLine[x- 1]];
			}

			//1st line
			srcLineAbove= srcLine- srcStride;
			firstLineStartX = isAboveLeftAvail ? 0 : 1;
			firstLineEndX   = isAboveAvail? endX: 1;
			for(x= firstLineStartX; x< firstLineEndX; x++)
			{
				int aux = srcLine[x] - srcLineAbove[x- 1];
				edgeType  =  ((aux==0)?0:(SIGN(aux))) - signUpLine[x+1];//m_sign[srcLine[x] - srcLineAbove[x- 1]] - signUpLine[x+1];
				resLine[x] = clip(srcLine[x] + offset[edgeType], 0, 255);//offsetClip[srcLine[x] + offset[edgeType]];
			}
			srcLine  += srcStride;
			resLine  += resStride;


			//middle lines
			for (y= 1; y< height-1; y++)
			{
				int aux;
				srcLineBelow= srcLine+ srcStride;

				for (x=startX; x<endX; x++)
				{
					aux = srcLine[x] - srcLineBelow[x+ 1];
					signDown = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLine[x] - srcLineBelow[x+ 1]] ;
					edgeType =  signDown + signUpLine[x];
					resLine[x] = clip(srcLine[x] + offset[edgeType], 0, 255);//offsetClip[srcLine[x] + offset[edgeType]];

					signDownLine[x+1] = -signDown; 
				}
				aux = srcLineBelow[startX] - srcLine[startX-1];
				signDownLine[startX] = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLineBelow[startX] - srcLine[startX-1]];

				signTmpLine  = signUpLine;
				signUpLine   = signDownLine;
				signDownLine = signTmpLine;

				srcLine += srcStride;
				resLine += resStride;
			}

			//last line
			srcLineBelow= srcLine+ srcStride;
			lastLineStartX = isBelowAvail ? startX : (width -1);
			lastLineEndX   = isBelowRightAvail ? width : (width -1);
			for(x= lastLineStartX; x< lastLineEndX; x++)
			{
				int aux = srcLine[x] - srcLineBelow[x+ 1];
				edgeType =  ((aux==0)?0:(SIGN(aux))) + signUpLine[x];//m_sign[srcLine[x] - srcLineBelow[x+ 1]] + signUpLine[x];
				resLine[x] = clip(srcLine[x] + offset[edgeType], 0, 255);//offsetClip[srcLine[x] + offset[edgeType]];
			}
		}
		break;
	case SAO_TYPE_EO_45:
		{
			int8_t *signUpLine = m_signLineBuf1+1;
			int16_t *srcLineBelow, *srcLineAbove;

			offset += 2;
			startX = isLeftAvail ? 0 : 1;
			endX   = isRightAvail ? width : (width -1);

			//prepare 2nd line upper sign
			srcLineBelow= srcLine+ srcStride;
			for (x=startX-1; x< endX; x++)
			{
				int aux = srcLineBelow[x] - srcLine[x+1]; 
				signUpLine[x] = ((aux==0)?0:(SIGN(aux)));//(Char)m_sign[srcLineBelow[x] - srcLine[x+1]];
			}


			//first line
			srcLineAbove= srcLine- srcStride;
			firstLineStartX = isAboveAvail ? startX : (width -1 );
			firstLineEndX   = isAboveRightAvail ? width : (width-1);
			for(x= firstLineStartX; x< firstLineEndX; x++)
			{
				int aux = srcLine[x] - srcLineAbove[x+1];
				edgeType = ((aux==0)?0:(SIGN(aux))) - signUpLine[x-1];//m_sign[srcLine[x] - srcLineAbove[x+1]] -signUpLine[x-1];
				resLine[x] = clip(srcLine[x] + offset[edgeType], 0, 255);//offsetClip[srcLine[x] + offset[edgeType]];
			}
			srcLine += srcStride;
			resLine += resStride;

			//middle lines
			for (y= 1; y< height-1; y++)
			{
				int aux;
				srcLineBelow= srcLine+ srcStride;

				for(x= startX; x< endX; x++)
				{
					aux = srcLine[x] - srcLineBelow[x-1];
					signDown =  (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLine[x] - srcLineBelow[x-1]] ;
					edgeType =  signDown + signUpLine[x];
					resLine[x] = clip(srcLine[x] + offset[edgeType], 0, 255);//offsetClip[srcLine[x] + offset[edgeType]];
					signUpLine[x-1] = -signDown; 
				}
				aux = srcLineBelow[endX-1] - srcLine[endX];
				signUpLine[endX-1] = (aux==0)?0:(SIGN(aux));//(Char)m_sign[srcLineBelow[endX-1] - srcLine[endX]];
				srcLine  += srcStride;
				resLine += resStride;
			}

			//last line
			srcLineBelow= srcLine+ srcStride;
			lastLineStartX = isBelowLeftAvail ? 0 : 1;
			lastLineEndX   = isBelowAvail ? endX : 1;
			for(x= lastLineStartX; x< lastLineEndX; x++)
			{
				int aux = srcLine[x] - srcLineBelow[x-1];
				edgeType = ((aux==0)?0:(SIGN(aux))) + signUpLine[x]; //m_sign[srcLine[x] - srcLineBelow[x-1]] + signUpLine[x];
				resLine[x] = clip(srcLine[x] + offset[edgeType], 0, 255);//offsetClip[srcLine[x] + offset[edgeType]];
			}
		}
		break;
	case SAO_TYPE_BO:
		{
			int shiftBits = bit_depth - NUM_SAO_BO_CLASSES_LOG2;//((compIdx == Y_COMP)?g_bitDepthY:g_bitDepthC)- NUM_SAO_BO_CLASSES_LOG2;
			for (y=0; y< height; y++)
			{
				for (x=0; x< width; x++)
				{
					resLine[x] = clip(srcLine[x] + offset[srcLine[x] >> shiftBits], 0, 255);//offsetClip[ srcLine[x] + offset[srcLine[x] >> shiftBits] ];
				}
				srcLine += srcStride;
				resLine += resStride;
			}
		}
		break;
	default:
		{
			printf("Not a supported SAO types\n");
//			assert(0);
			exit(-1);
		}
	}
}

//Void TComSampleAdaptiveOffset::offsetCTU(Int ctu, TComPicYuv* srcYuv, TComPicYuv* resYuv, SAOBlkParam& saoblkParam, TComPic* pPic)
void offset_ctu(hvenc_engine_t *enc_engine, ctu_info_t *ctu, sao_blk_param_t* sao_blk_param)
{
	int y_pos, x_pos, height, width;
	int pic_height = enc_engine->pict_height[Y_COMP];
	int pic_width = enc_engine->pict_width[Y_COMP];
	int max_cu_size = enc_engine->max_cu_size;
	//int isLeftAvail,isRightAvail,isAboveAvail,isBelowAvail,isAboveLeftAvail,isAboveRightAvail,isBelowLeftAvail,isBelowRightAvail;
	int l_available, r_available, t_available, b_available, tl_available, bl_available, tr_available, br_available;
	int component;


	l_available = (ctu->x[Y_COMP]> 0);
	t_available = (ctu->y[Y_COMP]> 0);
	r_available = ((ctu->x[Y_COMP] + ctu->size) < enc_engine->pict_width[Y_COMP]);
	b_available = ((ctu->y[Y_COMP] + ctu->size) < enc_engine->pict_height[Y_COMP]);
	tl_available = t_available & l_available;
	tr_available = t_available & r_available;
	bl_available = b_available & l_available;
	br_available = b_available & r_available;

	if((sao_blk_param->offsetParam[Y_COMP].modeIdc == SAO_MODE_OFF) && (sao_blk_param->offsetParam[U_COMP].modeIdc == SAO_MODE_OFF) && (sao_blk_param->offsetParam[V_COMP].modeIdc == SAO_MODE_OFF))
	{
		return;
	}

	//block boundary availability
//	pPic->getPicSym()->deriveLoopFilterBoundaryAvailibility(ctu, isLeftAvail,isRightAvail,isAboveAvail,isBelowAvail,isAboveLeftAvail,isAboveRightAvail,isBelowLeftAvail,isBelowRightAvail);

//	y_pos   = (ctu->ctu_number / m_numCTUInWidth)*m_maxCUHeight;
//	x_pos   = (ctu->ctu_number % m_numCTUInWidth)*m_maxCUWidth;
	y_pos   = ctu->y[Y_COMP];//(ctu->ctu_number / enc_engine->pict_width_in_ctu)*enc_engine->max_cu_size;
	x_pos   = ctu->x[Y_COMP];//(ctu->ctu_number % enc_engine->pict_width_in_ctu)*enc_engine->max_cu_size;

	height = (y_pos + max_cu_size > pic_height)?(pic_height- y_pos):max_cu_size;
	width  = (x_pos + max_cu_size > pic_width)?(pic_width - x_pos):max_cu_size;


	if(ctu->ctu_number == 3)
	{
		int iiiii=0;
	}

	for(component= Y_COMP; component < NUM_PICT_COMPONENTS; component++)
	{
		sao_offset_t* ctb_offset = &sao_blk_param->offsetParam[component];

		if(ctb_offset->modeIdc != SAO_MODE_OFF)
		{
			int isLuma     = (component == Y_COMP);
			int formatShift= isLuma?0:1;

			int  blkWidth   = (width  >> formatShift);
			int  blkHeight  = (height >> formatShift);
			//int  blkYPos    = (y_pos   >> formatShift);
			//int  blkXPos    = (x_pos   >> formatShift);

			int decoded_buff_stride = WND_STRIDE_2D(enc_engine->curr_reference_frame->img, component);
			int16_t *decoded_buff  = WND_POSITION_2D(int16_t *, enc_engine->curr_reference_frame->img, component, ctu->x[component], ctu->y[component], 0, enc_engine->ctu_width);
			int src_buff_stride = WND_STRIDE_2D(enc_engine->sao_aux_wnd, component);
			int16_t *src_buff = WND_POSITION_2D(int16_t *, enc_engine->sao_aux_wnd, component, ctu->x[component], ctu->y[component], 0, enc_engine->ctu_width);

//			int  srcStride = isLuma?srcYuv->getStride():srcYuv->getCStride();
//			uint8_t* srcBlk    = getPicBuf(srcYuv, compIdx)+ (yPos >> formatShift)*srcStride+ (xPos >> formatShift);

//			Int  resStride  = isLuma?resYuv->getStride():resYuv->getCStride();
//			Pel* resBlk     = getPicBuf(resYuv, compIdx)+ blkYPos*resStride+ blkXPos;



/*			offsetBlock( compIdx, ctbOffset.typeIdc, ctbOffset.offset
				, srcBlk, resBlk, srcStride, resStride, blkWidth, blkHeight
				, isLeftAvail, isRightAvail
				, isAboveAvail, isBelowAvail
				, isAboveLeftAvail, isAboveRightAvail
				, isBelowLeftAvail, isBelowRightAvail
				);
*/

			offset_block(component, ctb_offset->typeIdc, ctb_offset->offset, src_buff, decoded_buff, src_buff_stride,  decoded_buff_stride, blkWidth, blkHeight,
				l_available,  r_available, t_available, b_available, tl_available, tr_available, bl_available, br_available, enc_engine->bit_depth);
		}
	} //component

}


void decide_blk_params(hvenc_engine_t *enc_engine, slice_t *currslice, ctu_info_t *ctu, sao_stat_data_t stats[][NUM_PICT_COMPONENTS][NUM_SAO_NEW_TYPES], int *slice_enable)	
{
	sao_blk_param_t* merge_list[12]; 
	int is_all_blks_disabled = FALSE;
	sao_blk_param_t mode_param;
	double minCost, modeCost;
	int ctu_idx = ctu->ctu_number;
	int component;
	int num_lcus_for_sao_off[NUM_PICT_COMPONENTS];

	if(!slice_enable[Y_COMP] && !slice_enable[U_COMP] && !slice_enable[V_COMP])
	{
		is_all_blks_disabled = TRUE;
	}

	//  m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[ SAO_CABACSTATE_PIC_INIT ]);

//	for(ctu_idx=0; ctu_idx< enc_engine->pict_total_ctu; ctu_idx++)
	{
		int mode;
		int merge_list_size;
		if(is_all_blks_disabled)
		{
			sao_blk_param_t *cp = &coded_params[ctu_idx];
			for(component=0; component< 3; component++)
			{
				cp->offsetParam[component].modeIdc = SAO_MODE_OFF;
				cp->offsetParam[component].typeIdc = -1;
				cp->offsetParam[component].typeAuxInfo = -1;
				memset(cp->offsetParam[component].offset, 0, sizeof(cp->offsetParam[component].offset));
			}
			return;
		}

		merge_list_size = get_merge_list(enc_engine, ctu_idx, recon_params, merge_list);

		if(ctu_idx == 6)
		{
			int iiiii=0;
		}


		minCost = MAX_COST;
		for(mode=0; mode < NUM_SAO_MODES; mode++)
		{
			switch(mode)
			{
			case SAO_MODE_OFF:
				{
					continue; //not necessary, since all-off case will be tested in SAO_MODE_NEW case.
				}
				break;
			case SAO_MODE_NEW:
				{
					//deriveModeNewRDO(ctu, mergeList, sliceEnabled, blkStats, mode_param, modeCost, m_pppcRDSbacCoder, SAO_CABACSTATE_BLK_CUR);
//					derive_mode_new_rdo(hvenc_engine_t* enc_engine, sao_stat_data_t stats[][NUM_PICT_COMPONENTS][NUM_SAO_NEW_TYPES], sao_blk_param_t *mode_param, int slice_enabled[] )
					derive_mode_new_rdo(enc_engine, stats[ctu_idx], &mode_param, &modeCost, slice_enable);
				}
				break;
			case SAO_MODE_MERGE:
				{
//					deriveModeMergeRDO(ctu, mergeList, sliceEnabled, blkStats , mode_param, modeCost, m_pppcRDSbacCoder, SAO_CABACSTATE_BLK_CUR);
					derive_mode_merge_rdo(enc_engine, merge_list, merge_list_size, slice_enable, stats[ctu_idx], &mode_param, &modeCost);//, TEncSbac** cabacCoderRDO, Int inCabacLabel)
				}
				break;
			default:
				{
					printf("Not a supported SAO mode\n");
					//assert(0);
					exit(-1);
				}
			}

			if(modeCost < minCost)
			{
				minCost = modeCost;
				coded_params[ctu_idx] = mode_param;
//				m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[ SAO_CABACSTATE_BLK_NEXT ]);

			}
		} //mode
//		m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[ SAO_CABACSTATE_BLK_NEXT ]);

		if(ctu_idx==3)
		{
			int iiiii=0;
		}
		//apply reconstructed offsets
		recon_params[ctu_idx] = coded_params[ctu_idx];

		//reconstructBlkSAOParam(reconParams[ctu], mergeList);
		reconstruct_blk_sao_param(&recon_params[ctu_idx], merge_list, merge_list_size);

//		offsetCTU(ctu, srcYuv, resYuv, reconParams[ctu], pic);
		offset_ctu(enc_engine, ctu, &recon_params[ctu_idx]);
	} //ctu

//#if SAO_ENCODING_CHOICE 
/*	num_lcus_for_sao_off[Y_COMP ] = num_lcus_for_sao_off[U_COMP]= num_lcus_for_sao_off[V_COMP]= 0;

	for (component=Y_COMP; component<NUM_PICT_COMPONENTS; component++)
	{
		for(ctu_idx=0; ctu_idx< enc_engine->pict_total_ctu; ctu_idx++)
		{
			if( recon_params->offsetParam[component].modeIdc == SAO_MODE_OFF)
			{
				num_lcus_for_sao_off[component]++;
			}
		}
	}
//#if SAO_ENCODING_CHOICE_CHROMA
	for (component=Y_COMP; component<NUM_PICT_COMPONENTS; component++)
	{
		sao_disabled_rate[component] = (double)num_lcus_for_sao_off[component]/(double)enc_engine->pict_total_ctu;
	}
*/
//#endif

}


void hmr_sao_hm(hvenc_engine_t *enc_engine, slice_t *currslice)
{
	int ctu_num;
	ctu_info_t* ctu;	
	int num_lcus_for_sao_off[NUM_PICT_COMPONENTS];
	int component;

	sao_init(enc_engine->bit_depth);

	memset(num_lcu_sao_off, 0, sizeof(num_lcu_sao_off));
	memset(sao_disabled_rate, 0, sizeof(sao_disabled_rate));


	//slice on/off 
	decide_pic_params(slice_enabled);// decidePicParams(sliceEnabled, pPic->getSlice(0)->getDepth()); 

	memset(recon_params, 0, sizeof(recon_params));

	for(ctu_num = 0;ctu_num < enc_engine->pict_total_ctu;ctu_num++)
	{
		ctu = &enc_engine->ctu_info[ctu_num];
		//		ctu->partition_list = enc_engine->thread[0]->deblock_partition_info;
		//		create_partition_ctu_neighbours(enc_engine->thread[0], ctu, ctu->partition_list);//this call should be removed

		get_ctu_stats(enc_engine, currslice, ctu, stat_data);	
//	}

//	for(ctu_num = 0;ctu_num < enc_engine->pict_total_ctu;ctu_num++)
//	{
//		ctu = &enc_engine->ctu_info[ctu_num];
		decide_blk_params(enc_engine, currslice, ctu, stat_data, slice_enabled);// decidePicParams(sliceEnabled, pPic->getSlice(0)->getDepth()); 
		reference_picture_border_padding_ctu(&enc_engine->curr_reference_frame->img, ctu);
	}

	num_lcus_for_sao_off[Y_COMP ] = num_lcus_for_sao_off[U_COMP]= num_lcus_for_sao_off[V_COMP]= 0;

	for (component=Y_COMP; component<NUM_PICT_COMPONENTS; component++)
	{
		for(ctu_num=0; ctu_num< enc_engine->pict_total_ctu; ctu_num++)
		{
			if( recon_params[ctu_num].offsetParam[component].modeIdc == SAO_MODE_OFF)
			{
				num_lcus_for_sao_off[component]++;
			}
		}
	}
//#if SAO_ENCODING_CHOICE_CHROMA
	for (component=Y_COMP; component<NUM_PICT_COMPONENTS; component++)
	{
		sao_disabled_rate[component] = (double)num_lcus_for_sao_off[component]/(double)enc_engine->pict_total_ctu;
	}

}