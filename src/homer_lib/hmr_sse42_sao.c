/*****************************************************************************
* hmr_sse42_sao.c : homerHEVC encoding library
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

#include <math.h>

#include "hmr_private.h"
#include "hmr_common.h"
#include "hmr_sse42_primitives.h"
#include "hmr_sse42_macros.h"
#include "hmr_os_primitives.h"

extern int skiped_lines_r[NUM_PICT_COMPONENTS];
extern int skiped_lines_b[NUM_PICT_COMPONENTS];

//void sse_get_ctu_stats(int16_t *decoded, int decoded_stride, uint8_t *orig, int orig_stride, sao_stat_data_t stats[][NUM_SAO_NEW_TYPES], unsigned int avaliable_mask, int height_luma, int width_luma, int component, int type_idx, int bit_depth)	
void sse_sao_get_ctu_stats(henc_thread_t *wpp_thread, slice_t *currslice, ctu_info_t* ctu, sao_stat_data_t stats[][NUM_SAO_NEW_TYPES])	
{
	ALIGN(16) int16_t edge_type_buff[64+8+1];
	ALIGN(16) int16_t diff_coeff_buff[64+8+1];

	int component;
	int l_available, r_available, t_available, b_available, tl_available, bl_available, tr_available, br_available;
	int height_luma = (ctu->y[Y_COMP] + ctu->size > wpp_thread->pict_height[Y_COMP])?(wpp_thread->pict_height[Y_COMP]-ctu->y[Y_COMP]):ctu->size;
	int width_luma = (ctu->x[Y_COMP] + ctu->size > wpp_thread->pict_width[Y_COMP])?(wpp_thread->pict_width[Y_COMP]-ctu->x[Y_COMP]):ctu->size;
	int calculate_preblock_stats = wpp_thread->enc_engine->calculate_preblock_stats;

	l_available = (ctu->x[Y_COMP]> 0);
	t_available = (ctu->y[Y_COMP]> 0);
	r_available = ((ctu->x[Y_COMP] + ctu->size) < wpp_thread->pict_width[Y_COMP]);
	b_available = ((ctu->y[Y_COMP] + ctu->size) < wpp_thread->pict_height[Y_COMP]);
	tl_available = t_available & l_available;
	tr_available = t_available & r_available;
	bl_available = b_available & l_available;
	br_available = b_available & r_available;

	for(component=Y_COMP; component < NUM_PICT_COMPONENTS; component++)
	{
		int type_idx;
		int skip_lines_r = skiped_lines_r[component];
		int skip_lines_b = skiped_lines_b[component];
		int decoded_buff_stride = WND_STRIDE_2D(wpp_thread->enc_engine->curr_reference_frame->img, component);
		int16_t *decoded_buff  = WND_POSITION_2D(int16_t *, wpp_thread->enc_engine->curr_reference_frame->img, component, ctu->x[component], ctu->y[component], 0, wpp_thread->ctu_width);
		int orig_buff_stride = WND_STRIDE_2D(wpp_thread->enc_engine->current_pict.img2encode->img, component);
		int16_t *orig_buff = WND_POSITION_2D(int16_t *, wpp_thread->enc_engine->current_pict.img2encode->img, component, ctu->x[component], ctu->y[component], 0, wpp_thread->ctu_width);
		int chroma_shift = (component==Y_COMP)?0:1;
		int height = height_luma>>chroma_shift;
		int width = width_luma>>chroma_shift;

		//getBlkStats
		for(type_idx=0; type_idx< NUM_SAO_NEW_TYPES; type_idx++)
		{
			sao_stat_data_t *curr_stat_data = &stats[component][type_idx];
			int64_t *diff = curr_stat_data->diff;
			int64_t *count = curr_stat_data->count;
			int x,y, start_x, start_y, end_x, end_y, first_line_start_x, first_line_end_x;
			int16_t sign_left, sign_right, sign_down;
			int16_t *src_line = decoded_buff;
			int src_stride = decoded_buff_stride;
			int16_t *org_line = orig_buff;
			int org_stride = orig_buff_stride;
			__m128_i16 one = sse_128_vector_i16(1);
			__m128_i16 one_ = sse_128_vector_i16(-1);
			switch(type_idx)
			{
			case SAO_TYPE_EO_0:
				{
					diff +=2;
					count+=2;
					end_y   = (b_available) ? (height - skip_lines_b) : height;
					start_x = (!calculate_preblock_stats) ? (l_available ? 0 : 1): (r_available ? (width - skip_lines_r) : (width - 1));
					end_x   = (!calculate_preblock_stats) ? (r_available ? (width - skip_lines_r) : (width - 1)): (r_available ? width : (width - 1));

//					src_line+=start_x;
//					org_line+=start_x;
//					end_x-=start_x;
					for (y=0; y<end_y; y++)
					{
						for (x=start_x; (x)<end_x; x+=8)
						{
							__m128_i16 src = sse_128_load_vector_u(src_line+x);
							__m128_i16 src_1 = sse_128_load_vector_u(src_line+x-1);
							__m128_i16 src1 = sse_128_load_vector_u(src_line+x+1);
							__m128_i16 aux = sse_128_sub_i16(src, src_1);
							__m128_i16 sign_left = sse_128_sign_16(one, aux);
							__m128_i16 aux2 = sse_128_sub_i16(src, src1);
							__m128_i16 sign_right = sse_128_sign_16(one, aux2);
							__m128_i16 edge_type = sse_128_add_i16(sign_left, sign_right);
							__m128_i16 diff = sse_128_sub_i16(sse_128_load_vector_u(org_line+x), src);
							sse_128_store_vector_u(edge_type_buff+x, edge_type);
							sse_128_store_vector_u(diff_coeff_buff+x, diff);
						}
						for (x=start_x; (x)<end_x; x++)
						{
							int ll=edge_type_buff[x]; 
							diff [ll] += diff_coeff_buff[x];//(org_line[x] - src_line[x]);
							count[ll] ++;
						}


/*						for (x=0; (x+8)<end_x; x+=8)
						{
							int i, iter_end = (end_x-x)<8?(end_x-x):8;
							__m128_i16 src = sse_128_load_vector_u(src_line+x);
							__m128_i16 src_1 = sse_128_load_vector_u(src_line+x-1);
							__m128_i16 src1 = sse_128_load_vector_u(src_line+x+1);
							__m128_i16 aux_ = sse_128_sub_i16(src, src_1);
							__m128_i16 sign_left_ = sse_128_sign_16(one, aux_);
							__m128_i16 aux2_ = sse_128_sub_i16(src, src1);
							__m128_i16 sign_right_ = sse_128_sign_16(one, aux2_);
							__m128_i16 edge_type_ = sse_128_add_i16(sign_left_, sign_right_);
							__m128_i16 diff_ = sse_128_sub_i16(_mm_cvtepu8_epi16(sse_128_load_vector_u(org_line+x)), src);
							sse_128_store_vector_a(edge_type, edge_type_);
							sse_128_store_vector_a(diff_coeff, diff_);
							for(i=0;i<8;i++)
							{
								int ll=edge_type[i]; 
								diff [ll] += diff_coeff[i];//(org_line[x] - src_line[x]);
								count[ll] ++;
							}
						}
						if(x<end_x)
						{
							int i, iter_end = (end_x-x);
							__m128_i16 src = sse_128_load_vector_u(src_line+x);
							__m128_i16 src_1 = sse_128_load_vector_u(src_line+x-1);
							__m128_i16 src1 = sse_128_load_vector_u(src_line+x+1);
							__m128_i16 aux_ = sse_128_sub_i16(src, src_1);
							__m128_i16 sign_left_ = sse_128_sign_16(one, aux_);
							__m128_i16 aux2_ = sse_128_sub_i16(src, src1);
							__m128_i16 sign_right_ = sse_128_sign_16(one, aux2_);
							__m128_i16 edge_type_ = sse_128_add_i16(sign_left_, sign_right_);
							__m128_i16 diff_ = sse_128_sub_i16(_mm_cvtepu8_epi16(sse_128_load_vector_u(org_line+x)), src);
							sse_128_store_vector_a(edge_type, edge_type_);
							sse_128_store_vector_a(diff_coeff, diff_);
							for(i=0;i<iter_end;i++)
							{
								diff [edge_type[i]] += diff_coeff[i];//(org_line[x] - src_line[x]);
								count[edge_type[i]] ++;
							}
						}
*/
						src_line  += src_stride;
						org_line  += org_stride;
					}				
				}
				break;
			case SAO_TYPE_EO_90:
				{
					int16_t *signUpLine = wpp_thread->sao_sign_line_buff1;//m_signLineBuf1;
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
					for (x=start_x; x<end_x; x+=8) 
					{
						__m128_i16 src = sse_128_load_vector_u(src_line+x);
						__m128_i16 src_above = sse_128_load_vector_u(srcLineAbove+x);
						__m128_i16 aux = sse_128_sub_i16(src, src_above);
						__m128_i16 sign_above = sse_128_sign_16(one, aux);
						sse_128_store_vector_u(signUpLine+x, sign_above);
					}

					for (y=start_y; y<end_y; y++)
					{
						srcLineBelow = src_line + src_stride;

						for (x=start_x; x<end_x; x+=8)
						{
							__m128_i16 src = sse_128_load_vector_u(src_line+x);
							__m128_i16 src_below = sse_128_load_vector_u(srcLineBelow+x);
							__m128_i16 aux = sse_128_sub_i16(src, src_below);
							__m128_i16 sign_below = sse_128_sign_16(one, aux);
							__m128_i16 next_sign_above = sse_128_sign_16(one_, aux);
							__m128_i16 sign_above = sse_128_load_vector_u(signUpLine+x);
							__m128_i16 edge_type = sse_128_add_i16(sign_below, sign_above);
							__m128_i16 diff = sse_128_sub_i16(sse_128_load_vector_u(org_line+x), src);
							sse_128_store_vector_u(signUpLine+x, next_sign_above);
							sse_128_store_vector_u(edge_type_buff+x, edge_type);
							sse_128_store_vector_u(diff_coeff_buff+x, diff);
						}
						for (x=start_x; x<end_x; x++)
						{
							int ll=edge_type_buff[x]; 
							diff [ll] += diff_coeff_buff[x];//(org_line[x] - src_line[x]);
							count[ll] ++;
						}

						src_line += src_stride;
						org_line += org_stride;
					}

				}
				break;
			case SAO_TYPE_EO_135:
				{
					int16_t *signUpLine, *signDownLine, *signTmpLine;
					int16_t* srcLineBelow, *srcLineAbove;
					diff +=2;
					count+=2;
					
					signUpLine = wpp_thread->sao_sign_line_buff1;//m_signLineBuf1;
					signDownLine = wpp_thread->sao_sign_line_buff2;//m_signLineBuf2;

					start_x = (!calculate_preblock_stats) ? (l_available ? 0 : 1) : (r_available ? (width - skip_lines_r) : (width - 1));

					end_x = (!calculate_preblock_stats) ? (r_available ? (width - skip_lines_r): (width - 1)) : (r_available ? width : (width - 1));
					end_y = b_available ? (height - skip_lines_b) : (height - 1);

					//prepare 2nd line's upper sign
					srcLineBelow = src_line + src_stride;
					for (x=start_x; x<end_x+1; x+=8) 
					{
						__m128_i16 src = sse_128_load_vector_u(src_line+x-1);
						__m128_i16 src_below = sse_128_load_vector_u(srcLineBelow+x);
						__m128_i16 aux = sse_128_sub_i16(src_below, src);
						__m128_i16 sign_below = sse_128_sign_16(one, aux);
						sse_128_store_vector_u(signUpLine+x, sign_below);
					}

					//1st line
					srcLineAbove = src_line - src_stride;
					first_line_start_x = (!calculate_preblock_stats) ? (tl_available ? 0 : 1) : start_x;
					first_line_end_x   = (!calculate_preblock_stats) ? (t_available ? end_x : 1) : end_x;
					for(x=first_line_start_x; x<first_line_end_x; x+=8)
					{
						__m128_i16 src = sse_128_load_vector_u(src_line+x);
						__m128_i16 src_above_1 = sse_128_load_vector_u(srcLineAbove+x-1);
						__m128_i16 aux = sse_128_sub_i16(src, src_above_1);
						__m128_i16 sign_above = sse_128_sign_16(one, aux);
						__m128_i16 edge_type = sse_128_sub_i16(sign_above, sse_128_load_vector_u(signUpLine+x+1));
						__m128_i16 diff = sse_128_sub_i16(sse_128_load_vector_u(org_line+x), src);
						sse_128_store_vector_u(edge_type_buff+x, edge_type);
						sse_128_store_vector_u(diff_coeff_buff+x, diff);
					}
					for(x=first_line_start_x; x<first_line_end_x; x++)
					{
						int ll=edge_type_buff[x]; 
						diff [ll] += diff_coeff_buff[x];//(org_line[x] - src_line[x]);
						count[ll] ++;
					}
					src_line  += src_stride;
					org_line  += org_stride;

					//middle lines
					for (y=1; y<end_y; y++)
					{
						int aux;
						srcLineBelow = src_line + src_stride;

						for (x=start_x; x<end_x; x+=8)
						{
							__m128_i16 src = sse_128_load_vector_u(src_line+x);
							__m128_i16 src_below1 = sse_128_load_vector_u(srcLineBelow+x+1);
							__m128_i16 aux = sse_128_sub_i16(src, src_below1);
							__m128_i16 sign_down = sse_128_sign_16(one, aux);
							__m128_i16 next_sign_down = sse_128_sign_16(one_, aux);
							__m128_i16 edge_type = sse_128_add_i16(sign_down, sse_128_load_vector_u(signUpLine+x));
							__m128_i16 diff = sse_128_sub_i16(sse_128_load_vector_u(org_line+x), src);
							sse_128_store_vector_u(signDownLine+x+1, next_sign_down);
							sse_128_store_vector_u(edge_type_buff+x, edge_type);
							sse_128_store_vector_u(diff_coeff_buff+x, diff);
						}
						for (x=start_x; x<end_x; x++)
						{
							int ll=edge_type_buff[x]; 
							diff [ll] += diff_coeff_buff[x];//(org_line[x] - src_line[x]);
							count[ll] ++;
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
					int16_t *signUpLine = wpp_thread->sao_sign_line_buff1+1;//m_signLineBuf1+1;
					int16_t *srcLineBelow, *srcLineAbove;
					diff +=2;
					count+=2;

					start_x = (!calculate_preblock_stats) ? (l_available ? 0 : 1) : (r_available ? (width - skip_lines_r) : (width - 1));
					end_x = (!calculate_preblock_stats) ? (r_available ? (width - skip_lines_r) : (width - 1)) : (r_available ? width : (width - 1));
					end_y = b_available ? (height - skip_lines_b) : (height - 1);

					//prepare 2nd line upper sign
					srcLineBelow = src_line + src_stride;
					for (x=start_x-1; x<end_x; x+=8) 
					{
						__m128_i16 src1 = sse_128_load_vector_u(src_line+x+1);
						__m128_i16 src_below = sse_128_load_vector_u(srcLineBelow+x);
						__m128_i16 aux = sse_128_sub_i16(src_below, src1);
						__m128_i16 sign_below = sse_128_sign_16(one, aux);
						sse_128_store_vector_u(signUpLine+x, sign_below);
					}

					//first line
					srcLineAbove = src_line - src_stride;
					first_line_start_x = (!calculate_preblock_stats) ? (t_available ? start_x : end_x) : start_x;
					first_line_end_x   = (!calculate_preblock_stats) ? ((!r_available && tr_available) ? width : end_x) : end_x;
					for(x=first_line_start_x; x<first_line_end_x; x+=8)
					{
						__m128_i16 src = sse_128_load_vector_u(src_line+x);
						__m128_i16 src_above1 = sse_128_load_vector_u(srcLineAbove+x+1);
						__m128_i16 aux = sse_128_sub_i16(src, src_above1);
						__m128_i16 sign_above = sse_128_sign_16(one, aux);
						__m128_i16 edge_type = sse_128_sub_i16(sign_above, sse_128_load_vector_u(signUpLine+x-1));
						__m128_i16 diff = sse_128_sub_i16(sse_128_load_vector_u(org_line+x), src);
						sse_128_store_vector_u(edge_type_buff+x, edge_type);
						sse_128_store_vector_u(diff_coeff_buff+x, diff);
					}
					for(x=first_line_start_x; x<first_line_end_x; x++)
					{
						int ll=edge_type_buff[x]; 
						diff [ll] += diff_coeff_buff[x];//(org_line[x] - src_line[x]);
						count[ll] ++;
					}

					src_line += src_stride;
					org_line += org_stride;

					//middle lines
					for (y=1; y<end_y; y++)
					{
						int aux;
						srcLineBelow = src_line + src_stride;

						for (x=start_x; x<end_x; x+=8)
						{
							__m128_i16 src = sse_128_load_vector_u(src_line+x);
							__m128_i16 src_below1 = sse_128_load_vector_u(srcLineBelow+x-1);
							__m128_i16 aux = sse_128_sub_i16(src, src_below1);
							__m128_i16 sign_down = sse_128_sign_16(one, aux);
							__m128_i16 next_sign_down = sse_128_sign_16(one_, aux);
							__m128_i16 edge_type = sse_128_add_i16(sign_down, sse_128_load_vector_u(signUpLine+x));
							__m128_i16 diff = sse_128_sub_i16(sse_128_load_vector_u(org_line+x), src);
							sse_128_store_vector_u(signUpLine+x-1, next_sign_down);
							sse_128_store_vector_u(edge_type_buff+x, edge_type);
							sse_128_store_vector_u(diff_coeff_buff+x, diff);
						}
						for (x=start_x; x<end_x; x++)
						{
							int ll=edge_type_buff[x]; 
							diff [ll] += diff_coeff_buff[x];//(org_line[x] - src_line[x]);
							count[ll] ++;
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
					int shiftBits = wpp_thread->bit_depth - NUM_SAO_BO_CLASSES_LOG2;
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
			}
		}
	}
}

