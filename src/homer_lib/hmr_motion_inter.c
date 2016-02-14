/*****************************************************************************
 * hmr_motion_inter.c : homerHEVC encoding library
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
#include <limits.h>
#include <memory.h>


#include "hmr_private.h"
#include "hmr_common.h"
//#include "hmr_profiler.h"

#include "hmr_sse42_functions.h"

static const int16_t zero_buff[256] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
										0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
										0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
										0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


int encode_inter_cu(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_cu_info, int depth, PartSize part_size_type, int *curr_sum, int gcnt)//depth = prediction depth
{		
	int ssd_;
	int pred_buff_stride, orig_buff_stride, residual_buff_stride, residual_dec_buff_stride, decoded_buff_stride;
	int16_t *orig_buff;
	int16_t *pred_buff, *residual_buff, *residual_dec_buff, *quant_buff, *iquant_buff, *decoded_buff;
//	uint8_t *cbf_buff = NULL;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL;
	int inv_depth;//, diff;//, is_filtered;
	int cu_mode = REG_DCT;//if !IsIntra(abs_index)cu_mode = REG_DCT;
		
	int curr_depth = curr_cu_info->depth;
	int curr_part_x = curr_cu_info->x_position;
	int curr_part_y = curr_cu_info->y_position;
	int curr_part_size = curr_cu_info->size;
	int curr_part_size_shift = et->max_cu_size_shift-curr_depth;
	int curr_scan_mode = find_scan_mode(TRUE, TRUE, curr_part_size, cu_mode, 0);
	int	 per = curr_cu_info->qp/6;
	int	 rem = curr_cu_info->qp%6;
	double div = 2.5;//et->enc_engine->performance_mode == PERF_FULL_COMPUTATION?1.75:2.5;
	double offset = 5.;//et->enc_engine->performance_mode == PERF_FULL_COMPUTATION?20.:.5;
	quant_wnd = et->transform_quant_wnd[curr_depth + 1 + (part_size_type!=SIZE_2Nx2N)];
	decoded_wnd = et->decoded_mbs_wnd[curr_depth + 1  + (part_size_type!=SIZE_2Nx2N)];
//	cbf_buff = et->cbf_buffs[Y_COMP][curr_depth];

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd[0], Y_COMP);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
	orig_buff = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, Y_COMP);
	residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_dec_buff_stride = WND_STRIDE_2D(et->residual_dec_wnd, Y_COMP);
	residual_dec_buff = WND_POSITION_2D(int16_t *, et->residual_dec_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, Y_COMP, gcnt, et->ctu_width, (curr_cu_info->abs_index<<et->num_partitions_in_cu_shift));
	iquant_buff = WND_POSITION_1D(int16_t  *, et->itransform_iquant_wnd, Y_COMP, gcnt, et->ctu_width, (curr_cu_info->abs_index<<et->num_partitions_in_cu_shift));
	decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, Y_COMP);
	decoded_buff = WND_POSITION_2D(int16_t *, *decoded_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);

	inv_depth = (et->max_cu_size_shift - curr_depth);
//	diff = min(abs(cu_mode - HOR_IDX), abs(cu_mode - VER_IDX));

	//2d -> 1D buffer
	et->funcs->transform(et->bit_depth, residual_buff, et->pred_aux_buff, residual_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, curr_part_size_shift, cu_mode, quant_buff);

	et->funcs->quant(et, et->pred_aux_buff, quant_buff, curr_scan_mode, curr_depth, Y_COMP, cu_mode, 0, curr_sum, curr_part_size, per, rem);

	curr_cu_info->inter_cbf[Y_COMP] = (( *curr_sum ? 1 : 0 ) << (curr_depth - depth));// + (part_size_type == SIZE_NxN)));
	curr_cu_info->inter_tr_idx = (curr_depth - depth);// + (part_size_type == SIZE_NxN));

	if(*curr_sum>0)//curr_cu_info->size)
	{
		uint32_t ssd_zero;
//		int16_t zero_buff[256];
//		memset(zero_buff,0,sizeof(zero_buff));
		ssd_zero = et->funcs->ssd16b(residual_buff, residual_buff_stride, (int16_t*)zero_buff, 0, curr_part_size);

		et->funcs->inv_quant(et, quant_buff, iquant_buff, curr_depth, Y_COMP, 0, curr_part_size, per, rem);

		//1D ->2D buffer
		et->funcs->itransform(et->bit_depth, residual_dec_buff, iquant_buff, residual_buff_stride, curr_part_size, curr_part_size, cu_mode, et->pred_aux_buff);

		ssd_ = et->funcs->ssd16b(residual_buff, residual_buff_stride, residual_dec_buff, residual_buff_stride, curr_part_size);

#ifndef COMPUTE_AS_HM
//		if(ssd_zero < clip((200./et->enc_engine->avg_dist),1.01,1.25)*ssd_)
		//if(ssd_zero < ssd_+curr_part_size**curr_sum)// clip((1000./et->enc_engine->avg_dist),1.01,1.5)**curr_sum)
		if(ssd_zero <= ssd_+clip(et->enc_engine->avg_dist/div-offset,1.,20000.)*(*curr_sum))
//		if(ssd_zero < clip((200./((double)ssd_/curr_cu_info->num_part_in_cu)),1.01,1.25)*ssd_)
		{
			memset(quant_buff, 0, curr_part_size*curr_part_size*sizeof(quant_buff[0]));
			*curr_sum = 0;
			curr_cu_info->inter_cbf[Y_COMP] = 0;//
			et->funcs->reconst(pred_buff, pred_buff_stride, quant_buff, 0, decoded_buff, decoded_buff_stride, curr_part_size);//quant buff is full of zeros - a memcpy could do
		}
		else
#endif
		et->funcs->reconst(pred_buff, pred_buff_stride, residual_dec_buff, residual_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	}
	else
	{
		ssd_ = et->funcs->ssd16b(residual_buff, residual_buff_stride, quant_buff, 0, curr_part_size);

		et->funcs->reconst(pred_buff, pred_buff_stride, quant_buff, 0, decoded_buff, decoded_buff_stride, curr_part_size);//quant buff is full of zeros - a memcpy could do
	}
	curr_cu_info->sum = *curr_sum;
//	ssd_ = et->funcs->ssd(orig_buff, orig_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	return ssd_;
}


 
extern const uint8_t chroma_scale_conversion_table[];

int encode_inter_cu_chroma(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_cu_info, int component, int depth, PartSize part_size_type, int *curr_sum, int gcnt)//depth = prediction depth
{		
	uint32_t ssd_;
	picture_t *currpict = &et->enc_engine->current_pict;
	slice_t *currslice = &currpict->slice;
	int pred_buff_stride, orig_buff_stride, residual_buff_stride, residual_dec_buff_stride, decoded_buff_stride;
	int16_t *orig_buff;
	int16_t *pred_buff, *residual_buff, *residual_dec_buff, *quant_buff, *iquant_buff, *decoded_buff;
//	uint8_t *cbf_buff = NULL;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL;
	int diff;//, is_filtered;
	int cu_mode = REG_DCT;//if !IsIntra(abs_index)cu_mode = REG_DCT;
	int original_depth = curr_cu_info->depth;
	cu_partition_info_t* processing_partition_info = (curr_cu_info->size_chroma!=2)?curr_cu_info:curr_cu_info->parent;
	int curr_depth = processing_partition_info->depth;
	int curr_part_x = processing_partition_info->x_position_chroma;
	int curr_part_y = processing_partition_info->y_position_chroma;
	int curr_part_size = processing_partition_info->size_chroma;
	int curr_part_size_shift = et->max_cu_size_shift-curr_depth-1;//420
	int curr_scan_mode = find_scan_mode(TRUE, TRUE, curr_part_size, cu_mode, 0);
	int chr_qp_offset = et->enc_engine->chroma_qp_offset;
	int qp_chroma = chroma_scale_conversion_table[clip(curr_cu_info->qp+chr_qp_offset,0,57)];
	double weight = pow( 2.0, (currslice->qp-chroma_scale_conversion_table[clip(currslice->qp+chr_qp_offset,0,57)])/3.0 );
	double div = 2.5;//et->enc_engine->performance_mode == PERF_FULL_COMPUTATION?1.75:2.5;
	double offset = 5.;//et->enc_engine->performance_mode == PERF_FULL_COMPUTATION?20.:.5;

	int per = qp_chroma/6;
	int rem = qp_chroma%6;

	quant_wnd = et->transform_quant_wnd[original_depth+1+(part_size_type!=SIZE_2Nx2N)];
	decoded_wnd = et->decoded_mbs_wnd[original_depth+1+(part_size_type!=SIZE_2Nx2N)];
//	cbf_buff = et->cbf_buffs[component][original_depth];

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd[0], component);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, component);
	orig_buff = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, component);
	residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_dec_buff_stride = WND_STRIDE_2D(et->residual_dec_wnd, component);
	residual_dec_buff = WND_POSITION_2D(int16_t *, et->residual_dec_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, component, gcnt, et->ctu_width, (processing_partition_info->abs_index<<et->num_partitions_in_cu_shift)>>2);//420
	iquant_buff = WND_POSITION_1D(int16_t  *, et->itransform_iquant_wnd, component, gcnt, et->ctu_width, (processing_partition_info->abs_index<<et->num_partitions_in_cu_shift)>>2);//420
	decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, component);
	decoded_buff = WND_POSITION_2D(int16_t *, *decoded_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);

//	inv_depth = (et->max_cu_size_shift - curr_depth);
	diff = min(abs(cu_mode - HOR_IDX), abs(cu_mode - VER_IDX));

	//2d -> 1D buffer
	et->funcs->transform(et->bit_depth, residual_buff, et->pred_aux_buff, residual_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, curr_part_size_shift, cu_mode, quant_buff);

	et->funcs->quant(et, et->pred_aux_buff, quant_buff, curr_scan_mode, curr_depth, component, cu_mode, 0, curr_sum, curr_part_size, per, rem);


	curr_cu_info->inter_cbf[component] = (( *curr_sum ? 1 : 0 ) << (original_depth-depth));//+(part_size_type==SIZE_NxN)));

	if(*curr_sum>0)//curr_cu_info->size_chroma)
	{
		uint32_t ssd_zero;
//		int16_t zero_buff[256];
//		memset(zero_buff,0,sizeof(zero_buff));
		ssd_zero = (uint32_t)(weight*et->funcs->ssd16b(residual_buff, residual_buff_stride, (int16_t*)zero_buff, 0, curr_part_size));

		et->funcs->inv_quant(et, quant_buff, iquant_buff, curr_depth, component, 0, curr_part_size, per, rem);

		//1D ->2D buffer
		et->funcs->itransform(et->bit_depth, residual_dec_buff, iquant_buff, residual_buff_stride, curr_part_size, curr_part_size, cu_mode, et->pred_aux_buff);
		ssd_ = (uint32_t)(weight*et->funcs->ssd16b(residual_buff, residual_buff_stride, residual_dec_buff, residual_buff_stride, curr_part_size));

#ifndef COMPUTE_AS_HM
//		if(ssd_zero < clip((200./et->enc_engine->avg_dist),1.01,1.25)*ssd_)
//		if(ssd_zero < ssd_+clip((1000./et->enc_engine->avg_dist),1.01,1.5)**curr_sum)
//		if(ssd_zero < ssd_+curr_part_size**curr_sum)
//		if(ssd_zero <= ssd_+clip(et->enc_engine->avg_dist/3.-10,1.,20000.)*(*curr_sum))
//		if(ssd_zero < clip((200./((double)ssd_/curr_cu_info->num_part_in_cu)),1.01,1.25)*ssd_)
		if(ssd_zero <= ssd_+clip(et->enc_engine->avg_dist/div-offset,1.,20000.)*(*curr_sum))
		{
			memset(quant_buff, 0, curr_part_size*curr_part_size*sizeof(quant_buff[0]));
			*curr_sum = 0;
			curr_cu_info->inter_cbf[component] = 0;//
			et->funcs->reconst(pred_buff, pred_buff_stride, quant_buff, 0, decoded_buff, decoded_buff_stride, curr_part_size);//quant buff is full of zeros - a memcpy could do
		}
		else
#endif
			et->funcs->reconst(pred_buff, pred_buff_stride, residual_dec_buff, residual_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	}
	else
	{
		ssd_ = (uint32_t)(weight*et->funcs->ssd16b(residual_buff, residual_buff_stride, quant_buff, 0, curr_part_size));

		et->funcs->reconst(pred_buff, pred_buff_stride, quant_buff, 0, decoded_buff, decoded_buff_stride, curr_part_size);//quant buff is full of zeros - a memcpy could do
	}
	curr_cu_info->sum += *curr_sum;

//	ssd_ = et->funcs->ssd(orig_buff, orig_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	return ssd_;
}



#define NTAPS_LUMA        8 ///< Number of taps for luma
#define NTAPS_CHROMA      4 ///< Number of taps for chroma
#define IF_INTERNAL_PREC 14 ///< Number of bits for internal precision
#define IF_FILTER_PREC    6 ///< Log2 of sum of filter taps
#define IF_INTERNAL_OFFS (1<<(IF_INTERNAL_PREC-1)) ///< Offset used internally

const int16_t luma_filter_coeffs[4][NTAPS_LUMA] =
{
  {  0, 0,   0, 64,  0,   0, 0,  0 },
  { -1, 4, -10, 58, 17,  -5, 1,  0 },
  { -1, 4, -11, 40, 40, -11, 4, -1 },
  {  0, 1,  -5, 17, 58, -10, 4, -1 }
};

const int16_t chroma_filter_coeffs[8][NTAPS_CHROMA] =
{
  {  0, 64,  0,  0 },
  { -2, 58, 10, -2 },
  { -4, 54, 16, -2 },
  { -6, 46, 28, -4 },
  { -4, 36, 36, -4 },
  { -4, 28, 46, -6 },
  { -2, 16, 54, -4 },
  { -2, 10, 58, -2 }
};


//TComInterpolationFilter::filterCopy(int bitDepth, const Pel *src, int srcStride, Short *dst, int dstStride, int width, int height, Bool isFirst, Bool isLast)
void filter_copy(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int row, col;
  
	if ( is_first == is_last )
	{
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col++)
			{
				dst[col] = src[col];
			}
			src += src_stride;
			dst += dst_stride;
		}              
	}
	else if ( is_first )
	{
		int shift = IF_INTERNAL_PREC - bit_depth;
    
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col++)
			{
				int16_t val = src[col] << shift;
				dst[col] = val - (int16_t)IF_INTERNAL_OFFS;
			}
			src += src_stride;
			dst += dst_stride;
		}          
	}
	else
	{
		int shift = IF_INTERNAL_PREC - bit_depth;
		int16_t offset = IF_INTERNAL_OFFS + (shift?(1 << (shift - 1)):0);
		int16_t max_val = (1 << bit_depth) - 1;
		int16_t min_val = 0;
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col++)
			{
				int16_t val = src[col];
				val = ( val + offset ) >> shift;
				val = clip(val, min_val, max_val);
				dst[col] = val;
			}
			src += src_stride;
			dst += dst_stride;
		}              
	}
}


void hmr_interpolation_filter_luma(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int num_taps = NTAPS_LUMA;//argument
	int ref_stride = ( is_vertical ) ? src_stride : 1;

	int offset;
	short maxVal;
	int headRoom = IF_INTERNAL_PREC - bit_depth;
	int shift = IF_FILTER_PREC;
	int row, col;
	int16_t c[8];
	const int16_t	*coeffs = luma_filter_coeffs[fraction];

	src -= ( num_taps/2 - 1 ) * ref_stride;

	c[0] = coeffs[0];
	c[1] = coeffs[1];
	c[2] = coeffs[2];
	c[3] = coeffs[3];
	c[4] = coeffs[4];
	c[5] = coeffs[5];
	c[6] = coeffs[6];
	c[7] = coeffs[7];
	if ( is_last )
	{
		shift += (is_first) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (is_first) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
		maxVal = (1 << bit_depth) - 1;
	}
	else
	{
		shift -= (is_first) ? headRoom : 0;
		offset = (is_first) ? -IF_INTERNAL_OFFS << shift : 0;
		maxVal = 0;
	}

	for (row = 0; row < height; row++)
	{
		for (col = 0; col < width; col++)
		{
			int sum;
			short val;
			sum  = src[ col + 0 * ref_stride] * c[0];
			sum += src[ col + 1 * ref_stride] * c[1];
			sum += src[ col + 2 * ref_stride] * c[2];
			sum += src[ col + 3 * ref_stride] * c[3];
			sum += src[ col + 4 * ref_stride] * c[4];
			sum += src[ col + 5 * ref_stride] * c[5];
			sum += src[ col + 6 * ref_stride] * c[6];
			sum += src[ col + 7 * ref_stride] * c[7];
			val = ( sum + offset ) >> shift;
			if ( is_last)
			{
				val = clip(val,0,maxVal);
			}
			dst[col] = val;
		}

		src += src_stride;
		dst += dst_stride;
	}
}

void hmr_interpolate_luma(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	if(fraction==0)
	{
		filter_copy(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
	}
	else
	{
		hmr_interpolation_filter_luma(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
	}
}



void hmr_half_pixel_estimation_luma_hm(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, cu_partition_info_t* curr_cu_info, int width, int height, int curr_part_size_shift, motion_vector_t *mv)
{
	int is_bi_predict = 0;

//	int x_vect = mv->hor_vector>>3;
//	int y_vect = mv->ver_vector>>3;
	int filter_size = NTAPS_LUMA;
	int half_filter_size = (filter_size>>1);
	int src_stride = reference_buff_stride;
	int16_t *src = reference_buff - half_filter_size*src_stride - 1;
	int16_t *dst;
	int dst_stride;
	int curr_part_x = curr_cu_info->x_position, curr_part_y = curr_cu_info->y_position;

	//copy original samples
	dst_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width+1, height+filter_size, INTERPOLATE_HOR, TRUE, FALSE);

	//half pixel horizontal interpolation 
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width+1, height+filter_size, INTERPOLATE_HOR, TRUE, FALSE);

	//copy original samples
	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], Y_COMP);
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + half_filter_size*src_stride+1;
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[0][0], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);

	//half pixel vertical interpolation of original pixels
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride+1;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width, height+1, INTERPOLATE_VERT, FALSE, TRUE);
  
	//copy half pixel horizontal interpolation
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width+1, height, INTERPOLATE_VERT, FALSE, TRUE);
	
	//half pixel vertical interpolation of horizontaly interpolated pixels
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width+1, height+1, INTERPOLATE_VERT, FALSE, TRUE);
}


void hmr_quarter_pixel_estimation_luma_hm(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, cu_partition_info_t* curr_cu_info, int width, int height, int curr_part_size_shift, motion_vector_t *half_pel_mv)
{
	int filter_size = NTAPS_LUMA;
	int half_filter_size = (filter_size>>1);
	int src_stride = reference_buff_stride;
	int16_t *src;
	int16_t *dst;
	int dst_stride;
	int curr_part_x = curr_cu_info->x_position, curr_part_y = curr_cu_info->y_position;
	int ext_height = (half_pel_mv->ver_vector == 0) ? height + filter_size : height + filter_size-1;
	//horizontal filter 1,0
	src = reference_buff - half_filter_size*src_stride - 1;
	if(half_pel_mv->ver_vector>0)
		src += src_stride;
	if(half_pel_mv->hor_vector>=0)
		src += 1;
	dst_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);

	src = reference_buff - half_filter_size*src_stride - 1;
	//horizontal filter 3,0
	if(half_pel_mv->ver_vector>0)
		src += src_stride;
	if(half_pel_mv->hor_vector>0)
		src += 1;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);

	//vertical filter 1,1
	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[1][1], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	if(half_pel_mv->ver_vector==0)
		src += src_stride;
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

	//vertical filter 1,3
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);

	if(half_pel_mv->ver_vector != 0)
	{
		//vertical filter 1,2
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector==0)
			src += src_stride;
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//vertical filter 3,2
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector==0)
			src += src_stride;
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
	}
	else
	{
		//vertical filter 1,0
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//vertical filter 3,0
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);    		
	}


	if(half_pel_mv->hor_vector != 0)
	{
		//filter 2,1
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->hor_vector>0)
			src += 1;
		if(half_pel_mv->ver_vector>=0)
			src += src_stride;
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//filter 2,3
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->hor_vector>0)
			src += 1;
		if(half_pel_mv->ver_vector>0)
			src += src_stride;
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
	}
	else
	{
		//filter 0,1
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector>=0)
			src += src_stride;
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//filter 0,3
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector>0)
			src += src_stride;
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
	}

	//filter 3,1
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	if(half_pel_mv->ver_vector==0)
		src += src_stride;
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

	//filter 3,3
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
}


void hmr_quarter_pixel_estimation_luma_fast(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, cu_partition_info_t* curr_cu_info, int width, int height, int curr_part_size_shift, motion_vector_t *half_pel_mv, int dir_x, int dir_y, int zero_curr_best)
{
	int filter_size = NTAPS_LUMA;
	int half_filter_size = (filter_size>>1);
	int src_stride = reference_buff_stride;
	int16_t *src;
	int16_t *dst;
	int dst_stride;
	int processed = 0;
	int curr_part_x = curr_cu_info->x_position, curr_part_y = curr_cu_info->y_position;
	int ext_height = (half_pel_mv->ver_vector == 0) ? height + filter_size : height + filter_size-1;
	//horizontal filter 1,0
	src = reference_buff - half_filter_size*src_stride - 1;
	if(half_pel_mv->ver_vector>0)
		src += src_stride;
	if(half_pel_mv->hor_vector>=0)
		src += 1;
	dst_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);

	src = reference_buff - half_filter_size*src_stride - 1;
	//horizontal filter 3,0
	if(half_pel_mv->ver_vector>0)
		src += src_stride;
	if(half_pel_mv->hor_vector>0)
		src += 1;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);

	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[1][1], Y_COMP);

	if(zero_curr_best)
	{
		if(dir_x<0)
		{
			//vertical filter 3,0
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);    					
			processed++;
		}
		else if(dir_x>0)
		{
			//vertical filter 1,0
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);			
			processed++;		
		}

		if(dir_y>0)
		{
			//filter 0,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector>=0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}
		else if(dir_y<0)
		{
			//filter 0,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector>0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x<=0 && dir_y<=0)
		{
			//filter 3,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;		
		}

		if(dir_x<=0 && dir_y>=0)
		{
			//filter 3,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x>=0 && dir_y<=0)
		{
			//vertical filter 1,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x>=0 && dir_y>=0)
		{
			//vertical filter 1,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}
	}
	else
	{
		if(dir_y==0)
		{
			if(dir_x>0)
			{
				//vertical filter 1,0
				src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
				dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
				et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);			
				processed++;
			}
			else if(dir_x<0)
			{
				//vertical filter 3,0
				src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
				dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
				et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);    					
				processed++;
			}
		}
		
		if(dir_x==0)
		{
			if(dir_y>0)
			{
				//filter 0,1
				src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
				dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
				if(half_pel_mv->ver_vector>=0)
					src += src_stride;
				et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
				processed++;
			}
			else if(dir_y<0)
			{
				//filter 0,3
				src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
				dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
				if(half_pel_mv->ver_vector>0)
					src += src_stride;
				et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
				processed++;
			}
		}

		if(dir_x>0 && dir_y>0)
		{
			//vertical filter 1,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x>0 && dir_y<0)
		{
			//vertical filter 1,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x<0 && dir_y>0)
		{
			//filter 3,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x<0 && dir_y<0)
		{
			//filter 3,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}


		if(dir_x>=0 && dir_y!=0)
		{
			//vertical filter 1,2
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}
		
		if(dir_x<=0 && dir_y!=0)
		{
			//vertical filter 3,2
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
			processed++;
		}

		if(dir_y>=0 && dir_x!=0)
		{
			//filter 2,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->hor_vector>0)
				src += 1;
			if(half_pel_mv->ver_vector>=0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}
		
		if(dir_y<=0 && dir_x!=0)
		{
			//filter 2,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->hor_vector>0)
				src += 1;
			if(half_pel_mv->ver_vector>0)
				src += src_stride;
			et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
			processed++;
		}
	}
	if(processed!=3)
		printf("processed=%d", processed);

}


void hmr_half_pixel_estimation_luma_fast(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, cu_partition_info_t* curr_cu_info, int width, int height, int curr_part_size_shift, motion_vector_t *mv, int dir_x, int dir_y)
{
	int is_bi_predict = 0;

//	int x_vect = mv->hor_vector>>3;
//	int y_vect = mv->ver_vector>>3;
	int filter_size = NTAPS_LUMA;
	int half_filter_size = (filter_size>>1);
	int src_stride = reference_buff_stride;
	int16_t *src = reference_buff - half_filter_size*src_stride - 1;
	int16_t *dst;
	int dst_stride;
	int curr_part_x = curr_cu_info->x_position, curr_part_y = curr_cu_info->y_position;

	dst_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], Y_COMP);
	if(dir_y!=0)
	{
		//copy original samples
		src_stride = reference_buff_stride;
		src = reference_buff - half_filter_size*src_stride - 1;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width+1, height+filter_size, INTERPOLATE_HOR, TRUE, FALSE);

		//half pixel vertical interpolation of original pixels
		src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], Y_COMP);
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride+1;
		dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[2][0], Y_COMP);
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width, height+1, INTERPOLATE_VERT, FALSE, TRUE);
/*
		src_stride = reference_buff_stride;
		src = reference_buff - half_filter_size*src_stride - 1  +  (half_filter_size-1)*src_stride+1;
		dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[2][0], Y_COMP);
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width, height+1, INTERPOLATE_VERT, TRUE, TRUE);
*/	}

	//half pixel horizontal interpolation 
	src_stride = reference_buff_stride;
	src = reference_buff - half_filter_size*src_stride - 1;
	dst_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[2], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width+1, height+filter_size, INTERPOLATE_HOR, TRUE, FALSE);

	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[2], Y_COMP);
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[0][2], Y_COMP);

	//copy half pixel horizontal interpolation
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 0, width+1, height, INTERPOLATE_HOR, FALSE, TRUE);

	//half pixel vertical interpolation of horizontaly interpolated pixels
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma_m_compensation(src, src_stride, dst, dst_stride, 2, width+1, height+1, INTERPOLATE_VERT, FALSE, TRUE);
}

void hmr_interpolate_chroma(int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int num_taps = NTAPS_CHROMA;//argument
	int ref_stride = ( is_vertical ) ? reference_buff_stride : 1;

	int offset;
	short maxVal;
	int headRoom = IF_INTERNAL_PREC - bit_depth;
	int shift = IF_FILTER_PREC;
	int row, col;
	int16_t c[4];
	const int16_t	*coeffs = chroma_filter_coeffs[fraction];


	reference_buff -= ( num_taps/2 - 1 ) * ref_stride;

	c[0] = coeffs[0];
	c[1] = coeffs[1];
	c[2] = coeffs[2];
	c[3] = coeffs[3];

	if ( is_last )
	{
		shift += (is_first) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (is_first) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
		maxVal = (1 << bit_depth) - 1;
	}
	else
	{
		shift -= (is_first) ? headRoom : 0;
		offset = (is_first) ? -IF_INTERNAL_OFFS << shift : 0;
		maxVal = 0;
	}

	for (row = 0; row < height; row++)
	{
		for (col = 0; col < width; col++)
		{
			int sum;
			short val;
			sum  = reference_buff[ col + 0 * ref_stride] * c[0];
			sum += reference_buff[ col + 1 * ref_stride] * c[1];
			sum += reference_buff[ col + 2 * ref_stride] * c[2];
			sum += reference_buff[ col + 3 * ref_stride] * c[3];

			val = ( sum + offset ) >> shift;
			if ( is_last)
			{
				val = clip(val,0,maxVal);
			}
			pred_buff[col] = val;
		}

		reference_buff += reference_buff_stride;
		pred_buff += pred_buff_stride;
	}
}

float squareRoot(float x)
{
  unsigned int i = *(unsigned int*) &x;

  // adjust bias
  i  += 127 << 23;
  // approximation of square root
  i >>= 1;

  return *(float*) &i;
}


#define DISABLING_CLIP_FOR_BIPREDME	1

//TComYuv::removeHighFreq
void remove_high_freq(int16_t* src, int src_stride, int16_t* dst, int dst_stride, int height, int width)
{
	int x, y;

	for ( y = 0; y < height; y++ )
	{
		for ( x = 0; x < width; x ++)
		{
#if DISABLING_CLIP_FOR_BIPREDME
			dst[x ] = 2 * dst[x] - src[x];
#else
			dst[x ] = ClipY(2 * dst[x] - src[x]);
#endif
		}
		src += src_stride;
		dst += dst_stride;
	}
}


//xCheckBestMVP
uint32_t select_mv_candidate(henc_thread_t* et, cu_partition_info_t* curr_cu_info, mv_candiate_list_t	*search_candidate_list, motion_vector_t *mv, int *best_candidate)
{
//	mv_candiate_list_t	*search_candidate_list = &et->amvp_candidates[ref_pic_list];
	int idx;
	uint32_t best_cost = INT_MAX, best_idx = 0;

	for (idx = 0; idx < search_candidate_list->num_mv_candidates; idx++)
	{
#ifdef COMPUTE_AS_HM
		double cost_mvx = 30*sqrt((float)abs(search_candidate_list->mv_candidates[idx].mv.hor_vector - mv->hor_vector));
		double cost_mvy = 30*sqrt((float)abs(search_candidate_list->mv_candidates[idx].mv.ver_vector - mv->ver_vector));
#else
		double cost_mvx = curr_cu_info->qp*squareRoot((float)abs(search_candidate_list->mv_candidates[idx].mv.hor_vector - mv->hor_vector));
		double cost_mvy = curr_cu_info->qp*squareRoot((float)abs(search_candidate_list->mv_candidates[idx].mv.ver_vector - mv->ver_vector));
#endif
		uint32_t cost = (uint32_t) (3.+cost_mvx+cost_mvy+.5);

		if(best_cost>cost)
		{
			best_cost = cost;
			best_idx = idx;		
		}
	}

	*best_candidate = best_idx;
	return best_cost;
}


uint32_t select_mv_candidate_fast(henc_thread_t* et, cu_partition_info_t* curr_cu_info, mv_candiate_list_t	*search_candidate_list, motion_vector_t *mv, int *best_candidate)
{
#ifdef COMPUTE_AS_HM
	return 0;
#else
//	mv_candiate_list_t	*search_candidate_list = &et->amvp_candidates[ref_pic_list];
	int idx;
	uint32_t best_cost = INT_MAX, best_idx = 0;

	for (idx = 0; idx < search_candidate_list->num_mv_candidates; idx++)
	{
		double correction = calc_mv_correction(curr_cu_info->qp, et->enc_engine->avg_dist);//.25+et->enc_engine->avg_dist*et->enc_engine->avg_dist/5000000.;
		double cost_mvx = correction*((float)abs(search_candidate_list->mv_candidates[idx].mv.hor_vector - mv->hor_vector));
		double cost_mvy = correction*((float)abs(search_candidate_list->mv_candidates[idx].mv.ver_vector - mv->ver_vector));

		uint32_t cost = (uint32_t) (cost_mvx+cost_mvy+.5);

		if(best_cost>cost)
		{
			best_cost = cost;
			best_idx = idx;		
		}
	}

	*best_candidate = best_idx;
	return best_cost;
#endif
}



static const int s_acMvRefineH_HM[9][2] =
{
	{  0,  0 }, // 0
	{  0, -1 }, // 1
	{  0,  1 }, // 2
	{ -1,  0 }, // 3
	{  1,  0 }, // 4
	{ -1, -1 }, // 5
	{  1, -1 }, // 6
	{ -1,  1 }, // 7
	{  1,  1 }  // 8
};

static const int s_acMvRefineH[8][2] =
{
//	{  0,  0 }, // 0
	{ -1, -1 }, // 1
	{  0, -1 }, // 2
	{  1, -1 }, // 3
	{  1,  0 }, // 4
	{  1,  1 }, // 5
	{  0,  1 }, // 6
	{ -1,  1 }, // 7
	{ -1,  0 }  // 8
};


static const int s_acMvRefineQ[9][2] =
{
	{  0,  0 }, // 0
	{  0, -1 }, // 1
	{  0,  1 }, // 2
	{ -1, -1 }, // 5
	{  1, -1 }, // 6
	{ -1,  0 }, // 3
	{  1,  0 }, // 4
	{ -1,  1 }, // 7
	{  1,  1 }  // 8
};

#ifdef COMPUTE_AS_HM
static const int diamond_small[][2] = {{0,-1},{-1,0},{1,0},{0,1}};
static const int diamond_big[][2] = {{0,-2},{-1,-1},{1,-1},{-2,0},{2,0},{-1,1},{1,1},{0,2}};
#else
static const int diamond_small[][2] = {{-1,0},{0,-1},{1,0},{0,1}};
static const int diamond_big[][2] = {{-2,0},{-1,-1},{0,-2},{1,-1},{2,0},{1,1},{0,2},{-1,1}};//{{0,-2},{-1,-1},{1,-1},{-2,0},{2,0},{-1,1},{1,1},{0,2}};
#endif
static const int square[][2] = {{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0}};

uint32_t hmr_motion_estimation_HM(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_cu_info, int16_t *orig_buff, int orig_buff_stride, int16_t *reference_buff, int reference_buff_stride, int curr_part_global_x, 
						int curr_part_global_y, int init_x, int init_y, int curr_part_size, int curr_part_size_shift, int search_range_x, int search_range_y, int frame_size_x, int frame_size_y, motion_vector_t *mv)
{
	int i; 
	int xlow, xhigh, ylow, yhigh;
	uint32_t curr_best_sad, best_sad;
	int curr_best_x, curr_best_y, best_x, best_y;
	int dist = 1;

	xlow=((curr_part_global_x - search_range_x)<0)?-curr_part_global_x:-search_range_x;
	xhigh=((curr_part_global_x + search_range_x)>(frame_size_x-curr_part_size))?frame_size_x-curr_part_global_x-curr_part_size:search_range_x;
	ylow=((curr_part_global_y - search_range_y)<0)?-curr_part_global_y:-search_range_y;
	yhigh=((curr_part_global_y + search_range_y)>(frame_size_y-curr_part_size))?frame_size_y-curr_part_global_y-curr_part_size:search_range_y;

	curr_best_x = init_x;
	curr_best_y = init_y;

	if (curr_best_x<xlow)
		curr_best_x=xlow;
	if (curr_best_x>xhigh)
		curr_best_x=xhigh;
	if (curr_best_y<ylow)
		curr_best_y=ylow;
	if (curr_best_y>yhigh)
		curr_best_y=yhigh;

	curr_best_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_best_y*reference_buff_stride+curr_best_x, reference_buff_stride, curr_cu_info->size);

	best_sad = curr_best_sad;
	best_x = curr_best_x;
	best_y = curr_best_y;

	for(i=0;i<sizeof(diamond_small)/sizeof(diamond_small[0]);i++)
	{
		int curr_x = best_x+diamond_small[i][0];
		int curr_y = best_y+diamond_small[i][1];
		uint32_t curr_sad;

		if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
		{
			curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
			}
		}
	}

	best_x = curr_best_x-1;
	best_y = curr_best_y-1;
	for(dist = 1; 2*dist<search_range_x; )//dist*=2)
	{
		if(best_x != curr_best_x || best_y != curr_best_y)
		{
			best_sad = curr_best_sad;
			best_x = curr_best_x;
			best_y = curr_best_y;
		}
		else
			dist*=2;

		for(i=0;i<sizeof(diamond_big)/sizeof(diamond_big[0]);i++)
		{
			int curr_x = best_x+diamond_big[i][0]*dist;
			int curr_y = best_y+diamond_big[i][1]*dist;
			uint32_t curr_sad;

			if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
			{
				curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);
				if(curr_sad < curr_best_sad)
				{
					curr_best_sad = curr_sad;
					curr_best_x = curr_x;
					curr_best_y = curr_y;
				}
			}
		}
	}

	best_sad = curr_best_sad;
	best_x = curr_best_x;
	best_y = curr_best_y;

	for(i=0;i<sizeof(diamond_small)/sizeof(diamond_small[0]);i++)
	{
		int curr_x = best_x+diamond_small[i][0];
		int curr_y = best_y+diamond_small[i][1];
		uint32_t curr_sad;

		if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
		{
			curr_sad = sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
			}
		}
	}

	best_sad = curr_best_sad;
	best_x = curr_best_x;
	best_y = curr_best_y;

	//perform half-sample and quarter-sample motion estimation
	{
		int curr_best_idx;
		int curr_part_x = curr_cu_info->x_position;
		int curr_part_y = curr_cu_info->y_position;
		uint32_t curr_sad;
		motion_vector_t half_mv;
		hmr_half_pixel_estimation_luma_hm(et, reference_buff+best_y*reference_buff_stride+best_x, reference_buff_stride, curr_cu_info, curr_part_size, curr_part_size, curr_part_size_shift, mv);

		curr_best_x = 0;
		curr_best_y = 0;
		curr_best_idx = 0;
		for (i = 0; i < 9; i++)
		{
			int curr_x = s_acMvRefineH_HM[i][0]*2;
			int curr_y = s_acMvRefineH_HM[i][1]*2;
			int src_stride;
			int16_t *src;

			src_stride = WND_STRIDE_2D(et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP);
			src = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			
			if ( curr_x == 2 && ( curr_y & 1 ) == 0 )
			{
				src += 1;
			}
			if ( ( curr_x & 1 ) == 0 && curr_y == 2 )
			{
				src += src_stride;
			}

			curr_sad = sad(orig_buff, orig_buff_stride, src, src_stride, curr_cu_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
				curr_best_idx = i;
			}
		}
		half_mv.hor_vector = s_acMvRefineH_HM[curr_best_idx][0];
		half_mv.ver_vector = s_acMvRefineH_HM[curr_best_idx][1];

		hmr_quarter_pixel_estimation_luma_hm(et, reference_buff+best_y*reference_buff_stride+best_x, reference_buff_stride, curr_cu_info, curr_part_size, curr_part_size, curr_part_size_shift, &half_mv);

		curr_best_x = half_mv.hor_vector*2;
		curr_best_y = half_mv.ver_vector*2;
		for (i = 0; i < 9; i++)
		{
			int curr_x = half_mv.hor_vector*2+s_acMvRefineQ[i][0]*1;
			int curr_y = half_mv.ver_vector*2+s_acMvRefineQ[i][1]*1;
			int src_stride;
			int16_t *src;

			src_stride = WND_STRIDE_2D(et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP);
			src = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			
			if ( curr_x == 2 && ( curr_y & 1 ) == 0 )
			{
				src += 1;
			}
			if ( ( curr_x & 1 ) == 0 && curr_y == 2 )
			{
				src += src_stride;
			}

			curr_sad = sad(orig_buff, orig_buff_stride, src, src_stride, curr_cu_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
				curr_best_idx = i;
			}
		}
	}



	best_sad = curr_best_sad;
	mv->hor_vector = (best_x<<2) + curr_best_x;
	mv->ver_vector = (best_y<<2) + curr_best_y;

	return best_sad;
}

//hay que reducir el numero de parametros
uint32_t hmr_bi_motion_estimation_HM(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_cu_info, int16_t *orig_buff, int orig_buff_stride, int16_t *reference_buff, int reference_buff_stride, int curr_part_global_x, 
						int curr_part_global_y, int init_x, int init_y, int curr_part_size, int curr_part_size_shift, int search_range_x, int search_range_y, int frame_size_x, int frame_size_y, motion_vector_t *mv)
{
	int j, i; 
	int xlow, xhigh, ylow, yhigh;
	uint32_t curr_best_sad = MAX_COST, best_sad;
	int curr_best_x = 0, curr_best_y = 0, best_x, best_y;
	int16_t *ref_buff = reference_buff;

	xlow=((curr_part_global_x - search_range_x)<0)?-curr_part_global_x:-search_range_x;
	xhigh=((curr_part_global_x + search_range_x)>(frame_size_x-curr_part_size))?frame_size_x-curr_part_global_x-curr_part_size:search_range_x;
	ylow=((curr_part_global_y - search_range_y)<0)?-curr_part_global_y:-search_range_y;
	yhigh=((curr_part_global_y + search_range_y)>(frame_size_y-curr_part_size))?frame_size_y-curr_part_global_y-curr_part_size:search_range_y;

//	curr_best_x = init_x;
//	curr_best_y = init_y;
	ref_buff+=ylow*reference_buff_stride;
	for(j=ylow;j<yhigh;j++)
	{
		for(i=xlow;i<xhigh;i++)
		{
			uint32_t curr_sad;

			curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, ref_buff+i, reference_buff_stride, curr_cu_info->size);

			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = i;
				curr_best_y = j;
			}
		}
		ref_buff+=reference_buff_stride;
	}

	best_sad = curr_best_sad;
	best_x = curr_best_x;
	best_y = curr_best_y;

	//perform half-sample and quarter-sample motion estimation
	{
		int curr_best_idx;
		int curr_part_x = curr_cu_info->x_position;
		int curr_part_y = curr_cu_info->y_position;
		uint32_t curr_sad;
		motion_vector_t half_mv;
		hmr_half_pixel_estimation_luma_hm(et, reference_buff+best_y*reference_buff_stride+best_x, reference_buff_stride, curr_cu_info, curr_part_size, curr_part_size, curr_part_size_shift, mv);

		curr_best_x = 0;
		curr_best_y = 0;
		curr_best_idx = 0;
		for (i = 0; i < 9; i++)
		{
			int curr_x = s_acMvRefineH_HM[i][0]*2;
			int curr_y = s_acMvRefineH_HM[i][1]*2;
			int src_stride;
			int16_t *src;

			src_stride = WND_STRIDE_2D(et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP);
			src = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			
			if ( curr_x == 2 && ( curr_y & 1 ) == 0 )
			{
				src += 1;
			}
			if ( ( curr_x & 1 ) == 0 && curr_y == 2 )
			{
				src += src_stride;
			}

			curr_sad = sad(orig_buff, orig_buff_stride, src, src_stride, curr_cu_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
				curr_best_idx = i;
			}
		}
		half_mv.hor_vector = s_acMvRefineH_HM[curr_best_idx][0];
		half_mv.ver_vector = s_acMvRefineH_HM[curr_best_idx][1];

		hmr_quarter_pixel_estimation_luma_hm(et, reference_buff+best_y*reference_buff_stride+best_x, reference_buff_stride, curr_cu_info, curr_part_size, curr_part_size, curr_part_size_shift, &half_mv);

		curr_best_x = half_mv.hor_vector*2;
		curr_best_y = half_mv.ver_vector*2;
		for (i = 0; i < 9; i++)
		{
			int curr_x = half_mv.hor_vector*2+s_acMvRefineQ[i][0]*1;
			int curr_y = half_mv.ver_vector*2+s_acMvRefineQ[i][1]*1;
			int src_stride;
			int16_t *src;

			src_stride = WND_STRIDE_2D(et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP);
			src = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			
			if ( curr_x == 2 && ( curr_y & 1 ) == 0 )
			{
				src += 1;
			}
			if ( ( curr_x & 1 ) == 0 && curr_y == 2 )
			{
				src += src_stride;
			}

			curr_sad = sad(orig_buff, orig_buff_stride, src, src_stride, curr_cu_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
				curr_best_idx = i;
			}
		}
	}



	best_sad = curr_best_sad;
	mv->hor_vector = (best_x<<2) + curr_best_x;
	mv->ver_vector = (best_y<<2) + curr_best_y;

	return best_sad;
}

uint32_t hmr_motion_estimation(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_cu_info, int16_t *orig_buff, int orig_buff_stride, int16_t *reference_buff, int reference_buff_stride, int curr_part_global_x, 
						int curr_part_global_y, int init_x, int init_y, int curr_part_size, int curr_part_size_shift, int search_range_x, int search_range_y, int frame_size_x, int frame_size_y, motion_vector_t *mv, motion_vector_t *subpix_mv, mv_candiate_list_t *amvp_candidate_list, uint32_t threshold, unsigned int action)
{
	int i, l; 
	int xlow, xhigh, ylow, yhigh;
	uint32_t prev_best_sad = 0, prev_best_rd = 0, curr_best_sad = 0, curr_best_rd = 0, best_sad = MAX_COST, best_rd = 0;
	int prev_best_x, prev_best_y, curr_best_x, curr_best_y, best_x = 0, best_y = 0;
	int dist = 1;
	int end;
	mv_candiate_list_t	*search_candidate_list = &et->mv_search_candidates;//&et->mv_candidates[REF_PIC_LIST_0];
	int mv_cost = 0;
	motion_vector_t half_mv;
	int next_start;
	int search_size, diamond_size;
	int best_candidate;
	subpix_mv->hor_vector = 0;
	subpix_mv->ver_vector = 0;

	threshold = 0;

	xlow=((curr_part_global_x - search_range_x)<0)?-curr_part_global_x:-search_range_x;
	xhigh=((curr_part_global_x + search_range_x)>(frame_size_x-curr_part_size))?frame_size_x-curr_part_global_x-curr_part_size:search_range_x;
	ylow=((curr_part_global_y - search_range_y)<0)?-curr_part_global_y:-search_range_y;
	yhigh=((curr_part_global_y + search_range_y)>(frame_size_y-curr_part_size))?frame_size_y-curr_part_global_y-curr_part_size:search_range_y;

	if(action & MOTION_PEL_MASK)
	{
		curr_best_x = init_x;
		curr_best_y = init_y;

		if (curr_best_x<xlow)
			curr_best_x=xlow;
		if (curr_best_x>xhigh)
			curr_best_x=xhigh;
		if (curr_best_y<ylow)
			curr_best_y=ylow;
		if (curr_best_y>yhigh)
			curr_best_y=yhigh;

		curr_best_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_best_y*reference_buff_stride+curr_best_x, reference_buff_stride, curr_cu_info->size);

		mv->hor_vector = curr_best_x<<2;
		mv->ver_vector = curr_best_y<<2;

		subpix_mv->hor_vector = 0;
		subpix_mv->ver_vector = 0;

		mv_cost = select_mv_candidate_fast(et, curr_cu_info, amvp_candidate_list, mv, &best_candidate);
		curr_best_rd = curr_best_sad+mv_cost;

		best_sad = curr_best_sad;
		best_rd = curr_best_rd;
		best_x = curr_best_x;
		best_y = curr_best_y;

	//	if(best_sad<=2*curr_part_size*curr_part_size)
		if(best_sad<=threshold)
			goto last_search;
	//	curr_best_x = 0;
	//	curr_best_y = 0;

		for(i=0;i<search_candidate_list->num_mv_candidates;i++)
		{
			uint32_t curr_sad, curr_rd;
			int curr_x = search_candidate_list->mv_candidates[i].mv.hor_vector>>2;
			int curr_y = search_candidate_list->mv_candidates[i].mv.ver_vector>>2;

			if(curr_x == 0 && curr_y == 0)
			{
				continue;
			}
			if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
			{
				curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);
				mv->hor_vector = curr_x<<2;
				mv->ver_vector = curr_y<<2;
				mv_cost = select_mv_candidate_fast(et, curr_cu_info, amvp_candidate_list, mv, &best_candidate);
				curr_rd = curr_sad+mv_cost;
				if(curr_rd < curr_best_rd)
				{
					curr_best_sad = curr_sad;
					curr_best_rd = curr_rd;
					curr_best_x = curr_x;
					curr_best_y = curr_y;
				}
			}
		}

		best_sad = curr_best_sad;
		best_rd = curr_best_rd;
		best_x = curr_best_x;
		best_y = curr_best_y;

		//if(best_sad<=2*curr_part_size*curr_part_size)
		if(best_sad<=threshold)//curr_part_size*curr_part_size)
			goto last_search;

		for(i=0;i<sizeof(diamond_small)/sizeof(diamond_small[0]);i++)
		{
			int curr_x = best_x+diamond_small[i][0];
			int curr_y = best_y+diamond_small[i][1];
			uint32_t curr_sad, curr_rd;

			if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
			{
				curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);

				mv->hor_vector = curr_x<<2;
				mv->ver_vector = curr_y<<2;
				mv_cost = select_mv_candidate_fast(et, curr_cu_info, amvp_candidate_list, mv, &best_candidate);
				curr_rd = curr_sad+mv_cost;
				if(curr_rd < curr_best_rd)
				{
					curr_best_sad = curr_sad;
					curr_best_rd = curr_rd;
					curr_best_x = curr_x;
					curr_best_y = curr_y;
				}
			}
		}

		if(best_sad<=threshold)//curr_part_size*curr_part_size)
			goto last_search;

		dist = 2;
		if(best_x!=0 && best_y!=0)
			end = 4;
		else
			end = 8;

		diamond_size = sizeof(diamond_big)/sizeof(diamond_big[0]);
		for(l = 0; l < 1; l++)
		{
			next_start = 0;
			search_size = diamond_size;
 
			best_sad = curr_best_sad;
			best_rd = curr_best_rd;
			best_x = curr_best_x;
			best_y = curr_best_y;
		
			if(l==1)
			{
				dist=1;
				end=2;
			}

		//	for(dist = 1; dist<8; dist*=2)
			while(dist < end)
			{
//				for(i=0;i<diamond_size;i++)
				for(i=next_start;i<next_start+search_size;i++)
				{
					int index = i%diamond_size;
					int curr_x = best_x+diamond_big[index][0]*dist;
					int curr_y = best_y+diamond_big[index][1]*dist;
					uint32_t curr_sad, curr_rd;

					if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
					{
						curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);

						mv->hor_vector = curr_x<<2;
						mv->ver_vector = curr_y<<2;
						mv_cost = select_mv_candidate_fast(et, curr_cu_info, amvp_candidate_list, mv, &best_candidate);
						curr_rd = curr_sad+mv_cost;
						if(curr_rd < curr_best_rd)
						{
							curr_best_sad = curr_sad;
							curr_best_rd = curr_rd;
							curr_best_x = curr_x;
							curr_best_y = curr_y;
							next_start = (index-2+diamond_size)%diamond_size;
							search_size = diamond_size-3;
						}
					}
				}
				if(l==1 && (best_x != curr_best_x || best_y != curr_best_y))
				{
					best_sad = curr_best_sad;
					best_rd = curr_best_rd;
					best_x = curr_best_x;
					best_y = curr_best_y;
					dist = 2;
				}
				else
					dist*=2;
			}
		}

	last_search:

		best_sad = curr_best_sad;
		best_rd = curr_best_rd;
		best_x = curr_best_x;
		best_y = curr_best_y;
		prev_best_rd = MAX_COST;
		prev_best_sad = MAX_COST;
		prev_best_x = best_x;
		prev_best_y = best_y;

	//	if(best_sad<=threshold/2)
	//		goto end;

		next_start = 0;
		diamond_size = sizeof(diamond_small)/sizeof(diamond_small[0]);
		search_size = diamond_size;

		while(1)
		{
//			for(i=0;i<diamond_size;i++)
			for(i=next_start;i<next_start+search_size;i++)
			{
				int index = i%diamond_size;
				int curr_x = best_x+diamond_small[index][0];
				int curr_y = best_y+diamond_small[index][1];
				uint32_t curr_sad, curr_rd;

				if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
				{
					curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);
					mv->hor_vector = curr_x<<2;
					mv->ver_vector = curr_y<<2;
					mv_cost = select_mv_candidate_fast(et, curr_cu_info, amvp_candidate_list, mv, &best_candidate);
					curr_rd = curr_sad+mv_cost;
					if(curr_rd < curr_best_rd)
					{
						prev_best_sad = curr_best_sad;
						prev_best_rd = curr_best_rd;
						prev_best_x = curr_best_x;
						prev_best_y = curr_best_y;

						curr_best_sad = curr_sad;
						curr_best_rd = curr_rd;
						curr_best_x = curr_x;
						curr_best_y = curr_y;
						next_start = (index-1+diamond_size)%diamond_size;
						search_size = diamond_size-1;
					}
					else if(curr_rd < prev_best_rd)
					{
						prev_best_sad = curr_sad;
						prev_best_rd = curr_rd;
						prev_best_x = curr_x;
						prev_best_y = curr_y;			
					}
				}
			}
			if(best_x == curr_best_x && best_y == curr_best_y)
			{
				break;
			}
			else
			{
				best_sad = curr_best_sad;
				best_rd = curr_best_rd;
				best_x = curr_best_x;
				best_y = curr_best_y;
			}
		}
		best_sad = curr_best_sad;
		best_x = curr_best_x;
		best_y = curr_best_y;

		mv->hor_vector = (best_x<<2);// + curr_best_x;
		mv->ver_vector = (best_y<<2);// + curr_best_y;
		subpix_mv->hor_vector = 0;
		subpix_mv->ver_vector = 0;
	}

	if(action & MOTION_HALF_PEL_MASK)
	//perform half-sample and quarter-sample motion estimation
	{
		int curr_best_idx;
		int curr_part_x = curr_cu_info->x_position;
		int curr_part_y = curr_cu_info->y_position;
		uint32_t curr_sad;

		best_x = mv->hor_vector>>2;
		best_y = mv->ver_vector>>2;

		if(!(action & MOTION_PEL_MASK))
			curr_best_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+best_y*reference_buff_stride+best_x, reference_buff_stride, curr_cu_info->size);

		hmr_half_pixel_estimation_luma_hm(et, reference_buff+best_y*reference_buff_stride+best_x, reference_buff_stride, curr_cu_info, curr_part_size, curr_part_size, curr_part_size_shift, mv);

		curr_best_x = 0;
		curr_best_y = 0;
		curr_best_idx = 0;
		for (i = 0; i < 9; i++)
		{
			int curr_x = s_acMvRefineH_HM[i][0]*2;
			int curr_y = s_acMvRefineH_HM[i][1]*2;
			int src_stride;
			int16_t *src;

			src_stride = WND_STRIDE_2D(et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP);
			src = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			
			if ( curr_x == 2 && ( curr_y & 1 ) == 0 )
			{
				src += 1;
			}
			if ( ( curr_x & 1 ) == 0 && curr_y == 2 )
			{
				src += src_stride;
			}

			curr_sad = sad(orig_buff, orig_buff_stride, src, src_stride, curr_cu_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
				curr_best_idx = i;
			}
		}

		mv->hor_vector = (best_x<<2) + curr_best_x;
		mv->ver_vector = (best_y<<2) + curr_best_y;
		subpix_mv->hor_vector = curr_best_x;
		subpix_mv->ver_vector = curr_best_y;
		best_sad = curr_best_sad;

		if(action & MOTION_QUARTER_PEL_MASK)
		{

			half_mv.hor_vector = s_acMvRefineH_HM[curr_best_idx][0];
			half_mv.ver_vector = s_acMvRefineH_HM[curr_best_idx][1];

			hmr_quarter_pixel_estimation_luma_hm(et, reference_buff+best_y*reference_buff_stride+best_x, reference_buff_stride, curr_cu_info, curr_part_size, curr_part_size, curr_part_size_shift, &half_mv);

			curr_best_x = half_mv.hor_vector*2;
			curr_best_y = half_mv.ver_vector*2;
			for (i = 0; i < 9; i++)
			{
				int curr_x = half_mv.hor_vector*2+s_acMvRefineQ[i][0]*1;
				int curr_y = half_mv.ver_vector*2+s_acMvRefineQ[i][1]*1;
				int src_stride;
				int16_t *src;

				src_stride = WND_STRIDE_2D(et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP);
				src = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[curr_y&3][curr_x&3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			
				if ( curr_x == 2 && ( curr_y & 1 ) == 0 )
				{
					src += 1;
				}
				if ( ( curr_x & 1 ) == 0 && curr_y == 2 )
				{
					src += src_stride;
				}

				curr_sad = sad(orig_buff, orig_buff_stride, src, src_stride, curr_cu_info->size);
				if(curr_sad < curr_best_sad)
				{
					curr_best_sad = curr_sad;
					curr_best_x = curr_x;
					curr_best_y = curr_y;
					curr_best_idx = i;
				}
			}
			best_sad = curr_best_sad;
			mv->hor_vector = (best_x<<2) + curr_best_x;
			mv->ver_vector = (best_y<<2) + curr_best_y;
			subpix_mv->hor_vector = curr_best_x;
			subpix_mv->ver_vector = curr_best_y;
		}
	}

	return best_sad;
}



void hmr_motion_compensation_luma(henc_thread_t *et, cu_partition_info_t* curr_cu_info, int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int width, int height, int curr_part_size_shift, motion_vector_t *mv, int is_bi_predict)
{
//	int is_bi_predict = 0;
	int src_stride;
	int16_t *src;

//#ifdef COMPUTE_AS_HM
	int x_fraction = mv->hor_vector&0x3;
	int y_fraction = mv->ver_vector&0x3;

	int x_vect = mv->hor_vector>>2;
	int y_vect = mv->ver_vector>>2;
	src = reference_buff+y_vect*reference_buff_stride+x_vect;
	src_stride = reference_buff_stride;

	if (x_fraction == 0)
	{
		//vertical filter 
		et->funcs->interpolate_luma_m_compensation(src, src_stride, pred_buff, pred_buff_stride, y_fraction, width, height, INTERPOLATE_VERT, TRUE, !is_bi_predict);
	}
	else if (y_fraction  == 0)
	{
		et->funcs->interpolate_luma_m_compensation(src, src_stride, pred_buff, pred_buff_stride, x_fraction, width, height, INTERPOLATE_HOR, TRUE, !is_bi_predict);
	}
	else
	{
		int filter_size = NTAPS_LUMA;
		int half_filter_size = filter_size>>1;
		int16_t *temp_buff = WND_DATA_PTR(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP);
		int temp_buff_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], Y_COMP);

		//horizontal
		et->funcs->interpolate_luma_m_compensation(src - (half_filter_size-1)*src_stride, src_stride, temp_buff, temp_buff_stride, x_fraction, width, height+filter_size-1, INTERPOLATE_HOR, TRUE, FALSE);
		//vertical filter 
		et->funcs->interpolate_luma_m_compensation(temp_buff + (half_filter_size-1)*temp_buff_stride, temp_buff_stride, pred_buff, pred_buff_stride, y_fraction, width, height, INTERPOLATE_VERT, FALSE, !is_bi_predict);
	}
/*#else
	{
		int16_t *dst;
		int dst_stride;
		int x_fraction = curr_cu_info->subpix_mv->hor_vector;//    mv->hor_vector&0x3;
		int y_fraction = curr_cu_info->subpix_mv->ver_vector;	//mv->ver_vector&0x3;

		if (x_fraction == 0 && y_fraction==0)
		{
			int x_vect = mv->hor_vector>>2;
			int y_vect = mv->ver_vector>>2;
			reference_buff+= y_vect*reference_buff_stride+x_vect;

			src = reference_buff;
			src_stride = reference_buff_stride;
		}
		else
		{
			src_stride = WND_STRIDE_2D(et->filtered_block_wnd[y_fraction&3][x_fraction&3], Y_COMP);
			src = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[y_fraction&3][x_fraction&3], Y_COMP, curr_cu_info->x_position, curr_cu_info->y_position, 0, et->ctu_width);
			
			if ( x_fraction == 2 && ( y_fraction & 1 ) == 0 )
			{
				src += 1;
			}
			if ( ( x_fraction & 1 ) == 0 && y_fraction == 2 )
			{
				src += src_stride;
			}	
		}

		//use interpolate function to make a copy 
		dst_stride = pred_buff_stride;//WND_STRIDE_2D(et->filtered_block_wnd[0][0], Y_COMP);
		dst = pred_buff;//WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][0], Y_COMP, curr_cu_info->x_position, curr_cu_info->y_position, 0, et->ctu_width);
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, TRUE, TRUE);	
	}

#endif
*/
}





void hmr_motion_compensation_chroma(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int curr_part_size, int curr_part_size_shift, motion_vector_t *mv, int is_bi_predict)
{
//	int is_bi_predict = 0;
	int x_fraction = mv->hor_vector&0x7;
	int y_fraction = mv->ver_vector&0x7;

	int x_vect = mv->hor_vector>>3;
	int y_vect = mv->ver_vector>>3;
	reference_buff+= y_vect*reference_buff_stride+x_vect;


/*	if(x_fraction==0 && y_fraction==0)
	{
		int j, i;
		//copy samples
		for(j=0 ; j<curr_part_size ; j++)
		{

			for(i=0 ; i<curr_part_size ; i++)
			{
				pred_buff[i] = reference_buff[i];
			}
//			memcpy(pred_buff, reference_buff, curr_part_size);
			reference_buff+=reference_buff_stride;
			pred_buff+=pred_buff_stride;
		}
	}
	else 
*/	if(x_fraction == 0)
	{
		//vertical filter 
		et->funcs->interpolate_chroma_m_compensation(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, y_fraction, curr_part_size, curr_part_size, INTERPOLATE_VERT, TRUE, !is_bi_predict);
	}
	else if(y_fraction == 0)
	{
		//horizontal filter 
		et->funcs->interpolate_chroma_m_compensation(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, x_fraction, curr_part_size, curr_part_size, INTERPOLATE_HOR, TRUE, !is_bi_predict);	
	}
	else //if(x_fraction!=0 && y_fraction!=0)
	{
		int filter_size = NTAPS_CHROMA;
		int half_filter_size = filter_size>>1;
		int16_t *temp_buff = WND_DATA_PTR(int16_t *, et->filtered_block_temp_wnd[0], U_COMP);//this is only for temporal results so the component does not care
		int temp_buff_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], U_COMP);
		//horizontal
		et->funcs->interpolate_chroma_m_compensation(reference_buff - (half_filter_size-1)*reference_buff_stride, reference_buff_stride, temp_buff, temp_buff_stride, x_fraction, curr_part_size, curr_part_size+filter_size+1, INTERPOLATE_HOR, TRUE, FALSE);
		//vertical filter 
		et->funcs->interpolate_chroma_m_compensation(temp_buff + (half_filter_size-1)*temp_buff_stride, temp_buff_stride, pred_buff, pred_buff_stride, y_fraction, curr_part_size, curr_part_size, INTERPOLATE_VERT, FALSE, !is_bi_predict);
	}
}





int equal_motion(ctu_info_t* ctu_a, int abs_idx_a, ctu_info_t* ctu_b, int abs_idx_b)
{
	int ref_list_i;

	if(ctu_a->inter_mode[abs_idx_a] != ctu_b->inter_mode[abs_idx_b])
		return FALSE;

	for(ref_list_i=0;ref_list_i<2;ref_list_i++)
	{
		if((ctu_a->inter_mode[abs_idx_a]) & (1<<ref_list_i))
		{

			if(ctu_a->mv_ref[ref_list_i][abs_idx_a].hor_vector != ctu_b->mv_ref[ref_list_i][abs_idx_b].hor_vector || 
				ctu_a->mv_ref[ref_list_i][abs_idx_a].ver_vector != ctu_b->mv_ref[ref_list_i][abs_idx_b].ver_vector ||
				ctu_a->mv_ref_idx[ref_list_i][abs_idx_a] != ctu_b->mv_ref_idx[ref_list_i][abs_idx_b])
				return FALSE;	
		}	
	}
	return TRUE;
}

//getInterMergeCandidates
void get_merge_mvp_candidates(henc_thread_t* et, slice_t *currslice, ctu_info_t* ctu, cu_partition_info_t *curr_cu_info, PartSize part_size_type, uint8_t *inter_mode_neighbours)//get candidates for motion search from the neigbour CUs
{
	mv_candiate_list_t	*mv_candidate_list0 = &et->merge_mvp_candidates[REF_PIC_LIST_0];
	mv_candiate_list_t	*mv_candidate_list1 = &et->merge_mvp_candidates[REF_PIC_LIST_1];
	ctu_info_t	*ctu_left=NULL, *ctu_left_bottom=NULL, *ctu_top=NULL, *ctu_top_right=NULL, *ctu_top_left=NULL;
	uint	part_idx_lb, part_idx_l, part_idx_tr, part_idx_t, part_idx_tl;
//	int	added = FALSE;
	int	num_partitions_height = curr_cu_info->size>>2;
	int	num_partitions_width = curr_cu_info->size>>2;
	int	max_partitions_width = et->max_cu_size>>2;
	uint num_partitions_mask = curr_cu_info->num_part_in_cu-1;
	int abs_index_lb = et->enc_engine->raster2abs_table[curr_cu_info->raster_index + max_partitions_width*(num_partitions_height-1)];
	int abs_index_tl = et->enc_engine->raster2abs_table[curr_cu_info->raster_index];
	int abs_index_tr = et->enc_engine->raster2abs_table[curr_cu_info->raster_index + num_partitions_width-1];
	cu_partition_info_t *partition_info_lb = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_lb;
	cu_partition_info_t *partition_tl = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_tl;
	cu_partition_info_t *partition_tr = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_tr;
	int cand_is_inter[MERGE_MVP_MAX_NUM_CANDS];
	uint array_addr, cut_off;
	int num_ref_idx, r, ref_cnt;
	int l, cnt = 0;
	int avaliableA0 = FALSE, avaliableA1 = FALSE, avaliableB0 = FALSE, avaliableB1 = FALSE, avaliableB2 = FALSE;

	for(l=0;l<currslice->max_num_merge_candidates;l++)
	{
		cand_is_inter[l] = FALSE;
		mv_candidate_list0->mv_candidates[l].ref_idx = -1;
		mv_candidate_list1->mv_candidates[l].ref_idx = -1;
	}

	mv_candidate_list0->num_mv_candidates = 0;
	mv_candidate_list1->num_mv_candidates = 0;

	ctu_left = get_pu_left(ctu, partition_info_lb, &part_idx_l);//ctu->ctu_left;

	avaliableA1 = (ctu_left!=NULL && ctu_left->pred_mode[part_idx_l] != INTRA_MODE);
	if(avaliableA1)//ctu_left->mv_ref_idx[REF_PIC_LIST_0][part_idx_l]>=0)
	{
		cand_is_inter[cnt] = TRUE;
		inter_mode_neighbours[cnt] = ctu_left->inter_mode[part_idx_l];
		{
			mv_candidate_list0->mv_candidates[cnt].mv = ctu_left->mv_ref[REF_PIC_LIST_0][part_idx_l];
			mv_candidate_list0->mv_candidates[cnt].ref_idx = ctu_left->mv_ref_idx[REF_PIC_LIST_0][part_idx_l];
		}
		if(currslice->slice_type == B_SLICE)
		{
			mv_candidate_list1->mv_candidates[cnt].mv = ctu_left->mv_ref[REF_PIC_LIST_1][part_idx_l];
			mv_candidate_list1->mv_candidates[cnt].ref_idx = ctu_left->mv_ref_idx[REF_PIC_LIST_1][part_idx_l];
		}
		cnt++;
	}

	if(cnt >= currslice->max_num_merge_candidates)
	{
		currslice->max_num_merge_candidates = cnt;
		mv_candidate_list0->num_mv_candidates = mv_candidate_list1->num_mv_candidates = cnt;
		return;
	}

	//int equal_motion(ctu_info_t* ctu_a, int abs_idx_a, ctu_info_t* ctu_b, int abs_idx_b)
	ctu_top = get_pu_top(ctu, partition_tr, &part_idx_t, 0);

	avaliableB1 = (ctu_top!=NULL && ctu_top->pred_mode[part_idx_t] != INTRA_MODE);
	if(avaliableB1 && (!avaliableA1 || !equal_motion(ctu_left, part_idx_l, ctu_top, part_idx_t)))
	{
		cand_is_inter[cnt] = TRUE;
		inter_mode_neighbours[cnt] = ctu_top->inter_mode[part_idx_t];
//		if(ctu_top!=NULL && ctu_top->mv_ref_idx[REF_PIC_LIST_0][part_idx_t]>=0 && (ctu_left==NULL || ))
		{
			mv_candidate_list0->mv_candidates[cnt].mv = ctu_top->mv_ref[REF_PIC_LIST_0][part_idx_t];
			mv_candidate_list0->mv_candidates[cnt].ref_idx = ctu_top->mv_ref_idx[REF_PIC_LIST_0][part_idx_t];
		}
		if(currslice->slice_type == B_SLICE)// && ctu_top->mv_ref_idx[REF_PIC_LIST_1][part_idx_t]>=0 && (ctu_left==NULL || !equal_motion(ctu_left, part_idx_l, ctu_top, part_idx_t)))
		{
			mv_candidate_list1->mv_candidates[cnt].mv = ctu_top->mv_ref[REF_PIC_LIST_1][part_idx_t];
			mv_candidate_list1->mv_candidates[cnt].ref_idx = ctu_top->mv_ref_idx[REF_PIC_LIST_1][part_idx_t];
		}
		cnt++;
	}

	if(cnt >= currslice->max_num_merge_candidates)
	{
		currslice->max_num_merge_candidates = cnt;
		mv_candidate_list0->num_mv_candidates = mv_candidate_list1->num_mv_candidates = cnt;
		return;
	}

	partition_tr->top_right_neighbour = curr_cu_info->top_right_neighbour;
	ctu_top_right = get_pu_top_right(ctu, partition_tr, &part_idx_tr);

	avaliableB0 = (ctu_top_right!=NULL && ctu_top_right->pred_mode[part_idx_tr] != INTRA_MODE);
	if(avaliableB0 && (!avaliableB1 || !equal_motion(ctu_top, part_idx_t, ctu_top_right, part_idx_tr)))
	{
		cand_is_inter[cnt] = TRUE;
		inter_mode_neighbours[cnt] = ctu_top_right->inter_mode[part_idx_tr];
	//	ctu_top_right = get_pu_top_right(ctu, curr_cu_info, &aux_part_idx);
		//if(ctu_top_right->mv_ref_idx[REF_PIC_LIST_0][part_idx_tr]>=0 && (ctu_top==NULL || !equal_motion(ctu_top, part_idx_t, ctu_top_right, part_idx_tr)))
		{
			mv_candidate_list0->mv_candidates[cnt].mv = ctu_top_right->mv_ref[REF_PIC_LIST_0][part_idx_tr];
			mv_candidate_list0->mv_candidates[cnt].ref_idx = ctu_top_right->mv_ref_idx[REF_PIC_LIST_0][part_idx_tr];
		}
		if(currslice->slice_type == B_SLICE)// && ctu_top_right->mv_ref_idx[REF_PIC_LIST_1][part_idx_tr]>=0 && (ctu_top==NULL || !equal_motion(ctu_top, part_idx_t, ctu_top_right, part_idx_tr)))
		{
			mv_candidate_list1->mv_candidates[cnt].mv = ctu_top_right->mv_ref[REF_PIC_LIST_1][part_idx_tr];
			mv_candidate_list1->mv_candidates[cnt].ref_idx = ctu_top_right->mv_ref_idx[REF_PIC_LIST_1][part_idx_tr];
		}
		cnt++;
	}

	if(cnt >= currslice->max_num_merge_candidates)
	{
		currslice->max_num_merge_candidates = cnt;
		mv_candidate_list0->num_mv_candidates = mv_candidate_list1->num_mv_candidates = cnt;
		return;
	}

	//get spatial candidates
	partition_info_lb->left_bottom_neighbour = curr_cu_info->left_bottom_neighbour;
	ctu_left_bottom = get_pu_left_bottom(et, ctu, partition_info_lb, &part_idx_lb);
	avaliableA0 = (ctu_left_bottom!=NULL && ctu_left_bottom->pred_mode[part_idx_lb] != INTRA_MODE);
	if(avaliableA0 && (!avaliableA1 || !equal_motion(ctu_left, part_idx_l, ctu_left_bottom, part_idx_lb)))
	{
		cand_is_inter[cnt] = TRUE;
		inter_mode_neighbours[cnt] = ctu_left_bottom->inter_mode[part_idx_lb];
//		if(ctu_left_bottom->mv_ref_idx[REF_PIC_LIST_0][part_idx_lb]>=0 && (ctu_left==NULL || !equal_motion(ctu_left, part_idx_l, ctu_left_bottom, part_idx_lb)))
		{
			mv_candidate_list0->mv_candidates[cnt].mv = ctu_left_bottom->mv_ref[REF_PIC_LIST_0][part_idx_lb];
			mv_candidate_list0->mv_candidates[cnt].ref_idx = ctu_left_bottom->mv_ref_idx[REF_PIC_LIST_0][part_idx_lb];
		}
		if(currslice->slice_type == B_SLICE)// && ctu_left_bottom->mv_ref_idx[REF_PIC_LIST_1][part_idx_lb]>=0 && (ctu_left==NULL || !equal_motion(ctu_left, part_idx_l, ctu_left_bottom, part_idx_lb)))
		{
			mv_candidate_list1->mv_candidates[cnt].mv = ctu_left_bottom->mv_ref[REF_PIC_LIST_1][part_idx_lb];
			mv_candidate_list1->mv_candidates[cnt].ref_idx = ctu_left_bottom->mv_ref_idx[REF_PIC_LIST_1][part_idx_lb];
		}
		cnt++;
	}

	if(cnt >= currslice->max_num_merge_candidates)
	{
		currslice->max_num_merge_candidates = cnt;
		mv_candidate_list0->num_mv_candidates = mv_candidate_list1->num_mv_candidates = cnt;
		return;
	}

	if(cnt < 4)
	{
		ctu_top_left = get_pu_top_left(ctu, partition_tl, &part_idx_tl);
		avaliableB2 = (ctu_top_left!=NULL && ctu_top_left->pred_mode[part_idx_tl] != INTRA_MODE);
		if(avaliableB2 && (!avaliableA1 || !equal_motion(ctu_left, part_idx_l, ctu_top_left, part_idx_tl)) && (!avaliableB1 || !equal_motion(ctu_top, part_idx_t, ctu_top_left, part_idx_tl)))
		{
			cand_is_inter[cnt] = TRUE;
			inter_mode_neighbours[cnt] = ctu_top_left->inter_mode[part_idx_tl];
//			if(ctu_top_left->mv_ref_idx[REF_PIC_LIST_0][part_idx_tl]>=0 && (ctu_left==NULL || !equal_motion(ctu_left, part_idx_l, ctu_top_left, part_idx_tl))  && (ctu_top==NULL || !equal_motion(ctu_top, part_idx_t, ctu_top_left, part_idx_tl)))
			{
				mv_candidate_list0->mv_candidates[cnt].mv = ctu_top_left->mv_ref[REF_PIC_LIST_0][part_idx_tl];
				mv_candidate_list0->mv_candidates[cnt].ref_idx = ctu_top_left->mv_ref_idx[REF_PIC_LIST_0][part_idx_tl];
			}
			if(currslice->slice_type == B_SLICE)// if(ctu_top_left->mv_ref_idx[REF_PIC_LIST_0][part_idx_tl]>=0 && (ctu_left==NULL || !equal_motion(ctu_left, part_idx_l, ctu_top_left, part_idx_tl))  && (ctu_top==NULL || !equal_motion(ctu_top, part_idx_t, ctu_top_left, part_idx_tl)))
			{
				mv_candidate_list1->mv_candidates[cnt].mv = ctu_top_left->mv_ref[REF_PIC_LIST_1][part_idx_tl];
				mv_candidate_list1->mv_candidates[cnt].ref_idx = ctu_top_left->mv_ref_idx[REF_PIC_LIST_1][part_idx_tl];
			}
			cnt++;
		}
	}

	if(cnt >= currslice->max_num_merge_candidates)
	{
		currslice->max_num_merge_candidates = cnt;
		mv_candidate_list0->num_mv_candidates = mv_candidate_list1->num_mv_candidates = cnt;
		return;
	}

	array_addr = cnt;
	cut_off = array_addr;

	if ( currslice->slice_type == B_SLICE)
	{
		uint priority_list0[12] = {0 , 1, 0, 2, 1, 2, 0, 3, 1, 3, 2, 3};
		uint priority_list1[12] = {1 , 0, 2, 0, 2, 1, 3, 0, 3, 1, 3, 2};
		int idx;


		for (idx=0; idx<cut_off*(cut_off-1) && array_addr!= currslice->max_num_merge_candidates; idx++)
		{
			int i = priority_list0[idx]; 
			int j = priority_list1[idx];
			if (cand_is_inter[i] && cand_is_inter[j] && (inter_mode_neighbours[i] & 0x1) && (inter_mode_neighbours[j] & 0x2))
			{
				uint ref_poc0, ref_poc1;
				cand_is_inter[array_addr] = TRUE;
				inter_mode_neighbours[array_addr] = 3;

				// get Mv from cand[i] and cand[j]
				//pcMvFieldNeighbours[array_addr << 1].setMvField(pcMvFieldNeighbours[i<<1].getMv(), pcMvFieldNeighbours[i<<1].getRefIdx());
				//pcMvFieldNeighbours[( array_addr << 1 ) + 1].setMvField(pcMvFieldNeighbours[(j<<1)+1].getMv(), pcMvFieldNeighbours[(j<<1)+1].getRefIdx());
				mv_candidate_list0->mv_candidates[array_addr] = mv_candidate_list0->mv_candidates[i];
				mv_candidate_list1->mv_candidates[array_addr] = mv_candidate_list1->mv_candidates[j];

				//int iRefPOCL0 = m_pcSlice->getRefPOC( REF_PIC_LIST_0, pcMvFieldNeighbours[(array_addr<<1)].getRefIdx() );
				//int iRefPOCL1 = m_pcSlice->getRefPOC( REF_PIC_LIST_1, pcMvFieldNeighbours[(array_addr<<1)+1].getRefIdx() );
				ref_poc0 = currslice->ref_poc_list[REF_PIC_LIST_0][mv_candidate_list0->mv_candidates[array_addr].ref_idx];
				ref_poc1 = currslice->ref_poc_list[REF_PIC_LIST_1][mv_candidate_list1->mv_candidates[array_addr].ref_idx];

				if (ref_poc0 == ref_poc1 && mv_candidate_list0->mv_candidates[array_addr].mv.hor_vector == mv_candidate_list1->mv_candidates[array_addr].mv.hor_vector 
					&& mv_candidate_list0->mv_candidates[array_addr].mv.ver_vector == mv_candidate_list1->mv_candidates[array_addr].mv.ver_vector)
				{
					cand_is_inter[array_addr] = FALSE;
				}
				else
				{
					array_addr++;
				}
			}
		}
	}

	if(cnt >= currslice->max_num_merge_candidates)
	{
		currslice->max_num_merge_candidates = cnt;
		mv_candidate_list0->num_mv_candidates = mv_candidate_list1->num_mv_candidates = cnt;
		return;
	}

	num_ref_idx = currslice->num_ref_idx[REF_PIC_LIST_0];
		
	r = 0;
	ref_cnt = 0;
	while(array_addr < currslice->max_num_merge_candidates)
	{
		cand_is_inter[array_addr] = TRUE;
		inter_mode_neighbours[array_addr] = 1;

		mv_candidate_list0->mv_candidates[array_addr].mv.hor_vector = 0;
		mv_candidate_list0->mv_candidates[array_addr].mv.ver_vector = 0;
		mv_candidate_list0->mv_candidates[array_addr].ref_idx = r;

		if(currslice->slice_type == B_SLICE)// if(ctu_top_left->mv_ref_idx[REF_PIC_LIST_0][part_idx_tl]>=0 && (ctu_left==NULL || !equal_motion(ctu_left, part_idx_l, ctu_top_left, part_idx_tl))  && (ctu_top==NULL || !equal_motion(ctu_top, part_idx_t, ctu_top_left, part_idx_tl)))
		{
			inter_mode_neighbours[array_addr] = 3;
			mv_candidate_list1->mv_candidates[array_addr].mv.hor_vector = 0;
			mv_candidate_list1->mv_candidates[array_addr].mv.ver_vector = 0;
			mv_candidate_list1->mv_candidates[array_addr].ref_idx = r;
		}
		array_addr++;

		if (ref_cnt == num_ref_idx - 1)
		{
			r = 0;
		}
		else
		{
			r++;
			ref_cnt++;
		}
	}

	mv_candidate_list0->num_mv_candidates = mv_candidate_list1->num_mv_candidates = array_addr;
}


int add_amvp_cand(mv_candiate_list_t *search_candidate_list, slice_t *currslice, ctu_info_t* ctu, int ref_pic_list, int ref_idx, int part_idx)
{
//	int added = FALSE;
	int ref_pic_list_2nd = (ref_pic_list==REF_PIC_LIST_0)?REF_PIC_LIST_1:REF_PIC_LIST_0;

	if(ctu!=NULL && ctu->mv_ref_idx[ref_pic_list][part_idx]>=0 && currslice->ref_pic_list[ref_pic_list][ref_idx]->temp_info.poc == currslice->ref_poc_list[ref_pic_list][ctu->mv_ref_idx[ref_pic_list][part_idx]])
	{
		search_candidate_list->mv_candidates[search_candidate_list->num_mv_candidates++].mv = ctu->mv_ref[ref_pic_list][part_idx];
		return TRUE;
	}

	if(ctu!=NULL && ctu->mv_ref_idx[ref_pic_list_2nd][part_idx]>=0 && currslice->ref_pic_list[ref_pic_list][ref_idx]->temp_info.poc == currslice->ref_poc_list[ref_pic_list_2nd][ctu->mv_ref_idx[ref_pic_list_2nd][part_idx]])
	{
		search_candidate_list->mv_candidates[search_candidate_list->num_mv_candidates++].mv = ctu->mv_ref[ref_pic_list_2nd][part_idx];
		return TRUE;
	}
	return FALSE;
}


int xGetDistScaleFactor(int iCurrPOC, int iCurrRefPOC, int iColPOC, int iColRefPOC)
{
	int iDiffPocD = iColPOC - iColRefPOC;
	int iDiffPocB = iCurrPOC - iCurrRefPOC;

	if( iDiffPocD == iDiffPocB )
	{
		return 4096;
	}
	else
	{
		int iTDB      = clip(iDiffPocB, -128, 127);
		int iTDD      = clip(iDiffPocD, -128, 127);
		int iX        = (0x4000 + abs(iTDD/2)) / iTDD;
		int iScale    = clip((iTDB * iX + 32) >> 6, -4096, 4095);
		return iScale;
	}
}

motion_vector_t scale_mv( motion_vector_t mv, int iScale )
{
	motion_vector_t aux;
	int mvx = clip( (iScale * mv.hor_vector + 127 + (iScale * mv.hor_vector < 0)) >> 8, -32768, 32767);
	int mvy = clip( (iScale * mv.ver_vector + 127 + (iScale * mv.ver_vector < 0)) >> 8, -32768, 32767);
	aux.hor_vector = mvx;
	aux.ver_vector = mvy;
	return aux;
}

int add_amvp_cand_order(mv_candiate_list_t *search_candidate_list, slice_t *currslice, ctu_info_t* ctu, int ref_pic_list, int ref_idx, int part_idx)
{
	int ref_pic_list_2nd;
	int curr_poc = currslice->poc;
	int curr_ref_poc = currslice->ref_pic_list[ref_pic_list][ref_idx]->temp_info.poc;
	int neib_poc = curr_poc;
	int neib_ref_poc;

	int bIsCurrRefLongTerm = FALSE;//currslice->ref_pic_list[ref_pic_list][ref_idx]->;// = m_pcSlice->getRefPic( eRefPicList, iRefIdx)->getIsLongTerm();
	int bIsNeibRefLongTerm = FALSE;// = false;

	if(ctu==NULL)
		return FALSE;

	ref_pic_list_2nd = REF_PIC_LIST_0;
	if( ref_pic_list == REF_PIC_LIST_0 )
	{
		ref_pic_list_2nd = REF_PIC_LIST_1;
	}
	else if ( ref_pic_list == REF_PIC_LIST_1)
	{
		ref_pic_list_2nd = REF_PIC_LIST_0;
	}

	if(ctu->mv_ref_idx[ref_pic_list][part_idx]>=0)// pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) >= 0)
	{
		//currslice->ref_pic_list[ref_pic_list][ref_idx]->temp_info.poc == currslice->ref_poc_list[ref_pic_list][ctu->mv_ref_idx[ref_pic_list][part_idx]])
		motion_vector_t mv_pred = ctu->mv_ref[ref_pic_list][part_idx];//pcTmpCU->getCUMvField(eRefPicList)->getMv(uiIdx);
		motion_vector_t rc_mv;
		neib_ref_poc = currslice->ref_poc_list[ref_pic_list][ctu->mv_ref_idx[ref_pic_list][part_idx]];//pcTmpCU->getSlice()->getRefPOC( eRefPicList, pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) );

		bIsNeibRefLongTerm = FALSE;//pcTmpCU->getSlice()->getRefPic( eRefPicList, pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) )->getIsLongTerm();
		if ( bIsCurrRefLongTerm == bIsNeibRefLongTerm ) 
		{
			if ( bIsCurrRefLongTerm || bIsNeibRefLongTerm )
			{
				rc_mv = mv_pred;
			}
			else
			{
				int iScale = xGetDistScaleFactor( curr_poc, curr_ref_poc, neib_poc, neib_ref_poc);
				if ( iScale == 4096 )
				{
					rc_mv = mv_pred;
				}
				else
				{
					rc_mv = scale_mv(mv_pred, iScale); //cMvPred.scaleMv( iScale );
				}
			}
			search_candidate_list->mv_candidates[search_candidate_list->num_mv_candidates++].mv = rc_mv;
			return TRUE;
		}
	}

	if(ctu->mv_ref_idx[ref_pic_list_2nd][part_idx]>=0)//( pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) >= 0)
	{
		motion_vector_t mv_pred = ctu->mv_ref[ref_pic_list_2nd][part_idx];//pcTmpCU->getCUMvField(eRefPicList2nd)->getMv(uiIdx);
		motion_vector_t rc_mv;
		neib_ref_poc = currslice->ref_poc_list[ref_pic_list_2nd][ctu->mv_ref_idx[ref_pic_list_2nd][part_idx]];//pcTmpCU->getSlice()->getRefPOC( eRefPicList2nd, pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) );

		bIsNeibRefLongTerm = FALSE;//pcTmpCU->getSlice()->getRefPic( eRefPicList2nd, pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) )->getIsLongTerm();
		if ( bIsCurrRefLongTerm == bIsNeibRefLongTerm ) 
		{
			if ( bIsCurrRefLongTerm || bIsNeibRefLongTerm )
			{
				rc_mv = mv_pred;
			}
			else
			{
				int iScale = xGetDistScaleFactor( curr_poc, curr_ref_poc, neib_poc, neib_ref_poc);
				if ( iScale == 4096 )
				{
					rc_mv = mv_pred;
				}
				else
				{
					rc_mv = scale_mv(mv_pred, iScale); //cMvPred.scaleMv( iScale );
				}
			}
			search_candidate_list->mv_candidates[search_candidate_list->num_mv_candidates++].mv = rc_mv;
			return TRUE;
		}
	}

/*	if(ctu!=NULL && ctu->mv_ref_idx[ref_pic_list][part_idx]>=0 && currslice->ref_pic_list[ref_pic_list][ref_idx]->temp_info.poc == currslice->ref_poc_list[ref_pic_list][ctu->mv_ref_idx[ref_pic_list][part_idx]])
	{
		search_candidate_list->mv_candidates[search_candidate_list->num_mv_candidates++].mv = ctu->mv_ref[ref_pic_list][part_idx];
		added = TRUE;
	}
*/	return FALSE;
}


//fillMvpCand
void get_amvp_candidates(henc_thread_t* et, slice_t *currslice, ctu_info_t* ctu, cu_partition_info_t *curr_cu_info, mv_candiate_list_t	*search_candidate_list, int ref_pic_list, int ref_idx, PartSize part_size_type)//get candidates for motion search from the neigbour CUs
{
//	mv_candiate_list_t	*search_candidate_list = &et->amvp_candidates[ref_pic_list];
	ctu_info_t	*ctu_left=NULL, *ctu_left_bottom=NULL, *ctu_top=NULL, *ctu_top_right=NULL, *ctu_top_left=NULL;
	uint	part_idx_lb, part_idx_l, aux_part_idx;
	int	added = FALSE,  added_smvp = FALSE;
	int	num_partitions_height = curr_cu_info->size>>2;
	int	num_partitions_width = curr_cu_info->size>>2;
	int	max_partitions_width = et->max_cu_size>>2;
	uint num_partitions_mask = curr_cu_info->num_part_in_cu-1;
	int abs_index_lb = et->enc_engine->raster2abs_table[curr_cu_info->raster_index + max_partitions_width*(num_partitions_height-1)];
	int abs_index_tl = et->enc_engine->raster2abs_table[curr_cu_info->raster_index];
	int abs_index_tr = et->enc_engine->raster2abs_table[curr_cu_info->raster_index + num_partitions_width-1];
	cu_partition_info_t *partition_info_lb = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_lb;
	cu_partition_info_t *partition_tl = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_tl;
	cu_partition_info_t *partition_tr = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_tr;
	int left_pixel_x = ctu->x[Y_COMP]+curr_cu_info->x_position+curr_cu_info->size;

	search_candidate_list->num_mv_candidates = 0;

	partition_info_lb->left_bottom_neighbour = curr_cu_info->left_bottom_neighbour;

	ctu_left_bottom = get_pu_left_bottom(et, ctu, partition_info_lb, &part_idx_lb);

//	tmpCU = getPUBelowLeft(idx, uiPartIdxLB);
	added_smvp = (ctu_left_bottom != NULL) && (ctu_left_bottom->pred_mode[part_idx_lb] != INTRA_MODE);

	if (!added_smvp)
	{
		ctu_left = get_pu_left(ctu, partition_info_lb, &part_idx_l);//ctu->ctu_left;
		added_smvp = (ctu_left != NULL) && (ctu_left->pred_mode[part_idx_l] != INTRA_MODE);
	}

	//get spatial candidates
	added = add_amvp_cand(search_candidate_list, currslice, ctu_left_bottom, ref_pic_list, ref_idx, part_idx_lb);
	if(!added)
	{
		ctu_left = get_pu_left(ctu, partition_info_lb, &part_idx_l);//ctu->ctu_left;
		added = add_amvp_cand(search_candidate_list, currslice, ctu_left, ref_pic_list, ref_idx, part_idx_l);
	}

	if(!added)
	{
		added = add_amvp_cand_order(search_candidate_list, currslice, ctu_left_bottom, ref_pic_list, ref_idx, part_idx_lb);
		if(!added)
		{
			added = add_amvp_cand_order(search_candidate_list, currslice, ctu_left, ref_pic_list, ref_idx, part_idx_l);
		}
	}

	partition_tr->top_right_neighbour = curr_cu_info->top_right_neighbour;
	ctu_top_right = get_pu_top_right(ctu, partition_tr, &aux_part_idx);
	added = add_amvp_cand(search_candidate_list, currslice, ctu_top_right, ref_pic_list, ref_idx, aux_part_idx);

	if(!added)
	{
		ctu_top = get_pu_top(ctu, partition_tr, &aux_part_idx, 0);
		added = add_amvp_cand(search_candidate_list, currslice, ctu_top, ref_pic_list, ref_idx, aux_part_idx);
	}

	if(!added)
	{
		ctu_top_left = get_pu_top_left(ctu, partition_tl, &aux_part_idx);
		added = add_amvp_cand(search_candidate_list, currslice, ctu_top_left, ref_pic_list, ref_idx, aux_part_idx);
	}

	if (!added_smvp)
	{
		ctu_top_right = get_pu_top_right(ctu, partition_tr, &aux_part_idx);
		added = add_amvp_cand_order(search_candidate_list, currslice, ctu_top_right, ref_pic_list, ref_idx, aux_part_idx);

		if(!added)
		{
			ctu_top = get_pu_top(ctu, partition_tr, &aux_part_idx, 0);
			added = add_amvp_cand_order(search_candidate_list, currslice, ctu_top, ref_pic_list, ref_idx, aux_part_idx);
		}

		if(!added)
		{
			ctu_top_left = get_pu_top_left(ctu, partition_tl, &aux_part_idx);
			added = add_amvp_cand_order(search_candidate_list, currslice, ctu_top_left, ref_pic_list, ref_idx, aux_part_idx);
		}
	}


	if((ctu_left_bottom!=NULL && ctu_left_bottom->pred_mode[part_idx_lb] != INTRA_MODE)|| (ctu_left!=NULL && ctu_left->pred_mode[part_idx_l] != INTRA_MODE))
	{
		//reorder		�?
	}

	if(search_candidate_list->num_mv_candidates==2 && search_candidate_list->mv_candidates[0].mv.hor_vector == search_candidate_list->mv_candidates[1].mv.hor_vector && search_candidate_list->mv_candidates[0].mv.ver_vector == search_candidate_list->mv_candidates[1].mv.ver_vector)
	{
		search_candidate_list->num_mv_candidates = 1;
	}

	//get temporal candidates

	//.....
	if(search_candidate_list->num_mv_candidates>AMVP_MAX_NUM_CANDS)
		search_candidate_list->num_mv_candidates=AMVP_MAX_NUM_CANDS;

	while(search_candidate_list->num_mv_candidates<AMVP_MAX_NUM_CANDS)
	{
		search_candidate_list->mv_candidates[search_candidate_list->num_mv_candidates].mv.hor_vector = 0;
		search_candidate_list->mv_candidates[search_candidate_list->num_mv_candidates++].mv.ver_vector = 0;
	}
}


#define SET_ENC_INFO_BUFFS(et, cu_info, depth, abs_idx, num_partitions)																		\
{																																			\
	memset(&et->cbf_buffs[Y_COMP][depth][abs_idx], cu_info->inter_cbf[Y_COMP], num_partitions*sizeof(et->cbf_buffs[0][0][0]));				\
	memset(&et->cbf_buffs[U_COMP][depth][abs_idx], cu_info->inter_cbf[U_COMP], num_partitions*sizeof(et->cbf_buffs[0][0][0]));				\
	memset(&et->cbf_buffs[V_COMP][depth][abs_idx], cu_info->inter_cbf[V_COMP], num_partitions*sizeof(et->cbf_buffs[0][0][0]));				\
	memset(&et->tr_idx_buffs[depth][abs_idx], cu_info->inter_tr_idx, num_partitions*sizeof(et->tr_idx_buffs[0][0]));						\
}				


#define SET_INTER_MV_BUFFS(et, ctu, cu_info, abs_idx, num_partitions)																						\
{																																							\
	int i;																																					\
	for(i=0;i<num_partitions;i++)																															\
	{																																						\
		ctu->mv_ref[REF_PIC_LIST_0][abs_idx+i] = cu_info->inter_mv[REF_PIC_LIST_0];																			\
		ctu->mv_ref[REF_PIC_LIST_1][abs_idx+i] = cu_info->inter_mv[REF_PIC_LIST_1];																			\
	}																																						\
}


int hmr_cu_motion_estimation(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type, uint threshold, unsigned int action)
{
	double distortion = 0.;
	int i;
	uint32_t sad = 0, mv_total_cost = 0;
	picture_t *currpict = &et->enc_engine->current_pict;
	slice_t *currslice = &currpict->slice;
	ctu_info_t *ctu_rd = et->ctu_rd;
	int orig_buff_stride, reference_buff_stride;
	int orig_buff_stride_chroma, reference_buff_stride_chroma;
	int16_t  *orig_buff, *orig_buff_u, *orig_buff_v;
	int16_t  *reference_buff_cu_position, *reference_buff_cu_position_u, *reference_buff_cu_position_v;
	wnd_t *reference_wnd=NULL;//, *resi_wnd = NULL;
	int curr_part_size, curr_part_size_shift;
	int curr_part_size_chroma, curr_part_size_shift_chroma;
	int curr_part_x, curr_part_y, curr_part_global_x, curr_part_global_y;
	int curr_part_x_chroma, curr_part_y_chroma, curr_part_global_x_chroma, curr_part_global_y_chroma;
	int curr_depth = depth;
	cu_partition_info_t *parent_part_info;
	cu_partition_info_t *curr_cu_info;
//	int cbf_split[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
//	int acc_cost[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int num_partitions, npart, part_incr = 1;

	if(ctu->ctu_number == 2)
	{
		int iiiii=0;
	}

	curr_cu_info = &ctu->partition_list[et->partition_depth_start[curr_depth]]+part_position;

	if(part_size_type == SIZE_2Nx2N)
	{
		parent_part_info = curr_cu_info->parent;	
		num_partitions = 1;
	}
	else if(part_size_type == SIZE_NxN)
	{
		parent_part_info = curr_cu_info->parent;
		curr_cu_info = parent_part_info->children[0];
		num_partitions = 4;
		part_incr = 1;
	}

	sad = 0;
	for(npart=0;npart<num_partitions;npart+=part_incr)
	{
		int ref_list = 0;
		uint32_t best_cost[2] = {MAX_COST, MAX_COST};
		uint32_t mv_best_cost[2] = {0, 0};
		uint32_t best_cost_bi = MAX_COST;
		uint32_t mv_best_cost_bi;
		motion_vector_t mv_bi[2];
		motion_vector_t mv_best_cand_bi[2];
		motion_vector_t mvp_bi[2];
		int ref_index_bi[2];
		int mvp_idx[2];
		int best_cand_idx[2];

		curr_depth = curr_cu_info->depth;
		curr_part_x = curr_cu_info->x_position;
		curr_part_y = curr_cu_info->y_position;
		curr_part_global_x = ctu->x[Y_COMP]+curr_part_x;
		curr_part_global_y = ctu->y[Y_COMP]+curr_part_y;
		curr_part_size = curr_cu_info->size;
		curr_part_size_shift = et->max_cu_size_shift-curr_depth;
		curr_part_x_chroma = curr_cu_info->x_position_chroma;
		curr_part_y_chroma = curr_cu_info->y_position_chroma;
		curr_part_global_x_chroma = ctu->x[CHR_COMP]+curr_part_x_chroma;
		curr_part_global_y_chroma = ctu->y[CHR_COMP]+curr_part_y_chroma;
		curr_part_size_chroma = curr_cu_info->size_chroma;
		curr_part_size_shift_chroma = et->max_cu_size_shift-curr_depth-1;//420

		orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
		orig_buff = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
		orig_buff_stride_chroma = WND_STRIDE_2D(et->curr_mbs_wnd, CHR_COMP);
		orig_buff_u = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
		orig_buff_v = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

		if(et->enc_engine->num_encoded_frames==2 && ctu->ctu_number==5 && curr_cu_info->abs_index>=144 && curr_cu_info->depth==2)// && curr_cu_info->depth==3 && npart==2)// && curr_cu_info->abs_index>=32)// && curr_partition_info->abs_index>=192)// && curr_depth==2)
		{
			int iiiiii=0;
		}

		//unidirectional prediction
		for (ref_list = 0; ref_list <= (int)REF_PIC_LIST_1; ref_list++ )
		{
			int ref_idx;
			for ( ref_idx = 0; ref_idx < currslice->num_ref_idx[ref_list]; ref_idx++)
			{
				mv_candiate_list_t	amvp_candidate_list;
				int mvp_candidate_idx;
				uint32_t cost, mv_cost;
				uint32_t cost_temp_l0[32]; 
				motion_vector_t mv_temp_l0[32], mv, subpix_mv;
				reference_wnd = &currslice->ref_pic_list[ref_list][ref_idx]->img;//[0] up to now we only use one reference	
				reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
				reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
				reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, CHR_COMP);
				reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
				reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);

				get_amvp_candidates(et, currslice, ctu, curr_cu_info, &amvp_candidate_list, ref_list, ref_idx, part_size_type);//get candidates for advanced motion vector prediction from the neigbour CUs


				if(ref_list == REF_PIC_LIST_1 && currslice->list1_idx_to_list0_idx[ref_idx] >= 0)//equal reference frames in both lists
				{
					mv = mv_temp_l0[currslice->list1_idx_to_list0_idx[ref_idx]];
					cost = cost_temp_l0[currslice->list1_idx_to_list0_idx[ref_idx]];
					mv_cost = select_mv_candidate(et, curr_cu_info, &amvp_candidate_list, &mv, &mvp_candidate_idx);//&curr_cu_info->best_candidate_idx[ref_list]);

//					curr_cu_info->best_dif_mv[ref_list].hor_vector = mv.hor_vector-et->amvp_candidates[ref_list].mv_candidates[curr_cu_info->best_candidate_idx[ref_list]].mv.hor_vector;
//					curr_cu_info->best_dif_mv[ref_list].ver_vector = mv.ver_vector-et->amvp_candidates[ref_list].mv_candidates[curr_cu_info->best_candidate_idx[ref_list]].mv.ver_vector;
				}
				else
				{
#ifdef COMPUTE_AS_HM
					cost = hmr_motion_estimation_HM(et, ctu, curr_cu_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, curr_part_global_x, curr_part_global_y, 0, 0, curr_part_size, curr_part_size_shift, 64, 64, et->pict_width[Y_COMP], et->pict_height[Y_COMP], &mv);	
					mv_cost = select_mv_candidate(et, curr_cu_info, &amvp_candidate_list, &mv, &mvp_candidate_idx);//&curr_cu_info->best_candidate_idx[ref_list]);
//					curr_cu_info->best_dif_mv[ref_list].hor_vector = mv.hor_vector-et->amvp_candidates[ref_list].mv_candidates[curr_cu_info->best_candidate_idx[ref_list]].mv.hor_vector;
//					curr_cu_info->best_dif_mv[ref_list].ver_vector = mv.ver_vector-et->amvp_candidates[ref_list].mv_candidates[curr_cu_info->best_candidate_idx[ref_list]].mv.ver_vector;

					if(ref_list == REF_PIC_LIST_0)
					{
						mv_temp_l0[ref_idx] = mv;
						cost_temp_l0[ref_idx] = cost;
					}
				}
#else
					et->mv_search_candidates.num_mv_candidates = 0;
			//		get_mv_search_candidates(et, currslice, ctu, curr_cu_info, REF_PIC_LIST_0, ref_idx, part_size_type);//get candidates for motion search from the neigbour CUs
					for(i=0;i<amvp_candidate_list.num_mv_candidates;i++)//for(i=0;i<et->amvp_candidates[ref_list].num_mv_candidates;i++)
					{
						motion_vector_t mv = amvp_candidate_list.mv_candidates[i].mv; //et->amvp_candidates[ref_list].mv_candidates[i].mv;
						if(mv.hor_vector!=0 && mv.ver_vector!=0)
						{
							et->mv_search_candidates.mv_candidates[et->mv_search_candidates.num_mv_candidates++].mv = mv;	
						}
					}
					if(et->enc_engine->num_encoded_frames==2 && ctu->ctu_number==5 && curr_cu_info->abs_index==128 && curr_cu_info->depth==1)
						printf("\r\n ref_list:%d, parent_mv(%d, %d) \r\n", ref_list, curr_cu_info->parent->inter_mv[ref_list].hor_vector, curr_cu_info->parent->inter_mv[ref_list].ver_vector);
			//		if(curr_cu_info->parent && (curr_cu_info->parent->inter_mv[REF_PIC_LIST_0].hor_vector!=0 || curr_cu_info->parent->inter_mv[REF_PIC_LIST_0].ver_vector!=0))
					if(curr_cu_info->parent && (curr_cu_info->parent->inter_mv[ref_list].hor_vector!=0 && curr_cu_info->parent->inter_mv[ref_list].ver_vector!=0))
					{
						et->mv_search_candidates.mv_candidates[et->mv_search_candidates.num_mv_candidates++].mv = curr_cu_info->parent->inter_mv[ref_list];	
					}
		
					if((action & (MOTION_HALF_PEL_MASK|MOTION_QUARTER_PEL_MASK) && !(action & (MOTION_PEL_MASK))))
					{
						mv = curr_cu_info->inter_mv[ref_list];
						subpix_mv = curr_cu_info->subpix_mv[ref_list];		
					}


					cost = hmr_motion_estimation(et, ctu, curr_cu_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, curr_part_global_x, curr_part_global_y, 0, 0, curr_part_size, curr_part_size_shift, MOTION_SEARCH_RANGE_X, MOTION_SEARCH_RANGE_Y, et->pict_width[Y_COMP], et->pict_height[Y_COMP], &mv, &subpix_mv, &amvp_candidate_list, threshold, action);//|MOTION_HALF_PEL_MASK|MOTION_QUARTER_PEL_MASK);	

					//mv_cost = select_mv_candidate_fast(et, curr_cu_info, REF_PIC_LIST_0, &mv);
					mv_cost = select_mv_candidate_fast(et, curr_cu_info, &amvp_candidate_list, &mv, &mvp_candidate_idx);
//					select_mv_candidate_fast(henc_thread_t* et, cu_partition_info_t* curr_cu_info, mv_candiate_list_t	*search_candidate_list, motion_vector_t *mv, int *best_candidate)
					//mv_cost = select_mv_candidate(et, curr_cu_info, &amvp_candidate_list, &mv, &mvp_candidate_idx);
			//		mv_cost = select_mv_candidate(et, curr_cu_info, REF_PIC_LIST_0, &mv);
					curr_cu_info->subpix_mv[REF_PIC_LIST_0] = subpix_mv;

					if(ref_list == REF_PIC_LIST_0)
					{
						mv_temp_l0[ref_idx] = mv;
						cost_temp_l0[ref_idx] = cost;
					}
				}
#endif

				//-----------------Juan para deshabilitar la eleccion de la segunda referencia------------------------
/*				if(ref_idx>0)
				{
					cost = best_cost[ref_list]-1;
					mv_cost = mv_best_cost[ref_list]-1;
				}
*/				//-----------------Fin Juan para deshabilitar la eleccion de la segunda referencia------------------------
#ifdef COMPUTE_AS_HM
				if(cost<best_cost[ref_list])
#else
				if(cost+mv_cost<best_cost[ref_list]+mv_best_cost[ref_list])
#endif
				{
					best_cost[ref_list] = cost;
					mv_best_cost[ref_list] = mv_cost;
					//set mvs and ref_idx
					curr_cu_info->inter_mv[ref_list] = mv;
					curr_cu_info->inter_ref_index[ref_list] = ref_idx;					
					curr_cu_info->inter_mode = (1<<ref_list);

					//vector prediction
					et->amvp_candidates[ref_list] = amvp_candidate_list;
					curr_cu_info->best_candidate_idx[ref_list] = mvp_candidate_idx;
					curr_cu_info->best_dif_mv[ref_list].hor_vector = mv.hor_vector-et->amvp_candidates[ref_list].mv_candidates[mvp_candidate_idx].mv.hor_vector;
					curr_cu_info->best_dif_mv[ref_list].ver_vector = mv.ver_vector-et->amvp_candidates[ref_list].mv_candidates[mvp_candidate_idx].mv.ver_vector;
				}	
			}//for ( ref_idx = 0; ref_idx < currslice->num_ref_idx[ref_list]; ref_idx++)

		}//for (ref_list = 0; ref_list <= (int)REF_PIC_LIST_1; ref_list++ )

		//bi-prediction motion estimation
		if(currslice->slice_type == B_SLICE && curr_cu_info->size>=8)
		{
			int ref_idx = 0, num_ref_idx;
			uint cost, mv_cost;
			//curr_cu_info->best_candidate_idx[ref_pic_list] = best_idx;
			//curr_cu_info->best_dif_mv->hor_vector = mv->hor_vector-search_candidate_list->mv_candidates[best_idx].mv.hor_vector;
			//curr_cu_info->best_dif_mv->ver_vector = mv->ver_vector-search_candidate_list->mv_candidates[best_idx].mv.ver_vector;

			mv_bi[REF_PIC_LIST_0] = curr_cu_info->inter_mv[REF_PIC_LIST_0];
			mv_bi[REF_PIC_LIST_1] = curr_cu_info->inter_mv[REF_PIC_LIST_1];
			ref_index_bi[REF_PIC_LIST_0] = curr_cu_info->inter_ref_index[REF_PIC_LIST_0];
			ref_index_bi[REF_PIC_LIST_1] = curr_cu_info->inter_ref_index[REF_PIC_LIST_1];			
			mvp_idx[REF_PIC_LIST_0] = curr_cu_info->best_candidate_idx[0];
			mvp_idx[REF_PIC_LIST_1] = curr_cu_info->best_candidate_idx[1];

			mv_best_cand_bi[REF_PIC_LIST_0] = et->amvp_candidates[REF_PIC_LIST_0].mv_candidates[mvp_idx[REF_PIC_LIST_0]].mv;
			mv_best_cand_bi[REF_PIC_LIST_1] = et->amvp_candidates[REF_PIC_LIST_1].mv_candidates[mvp_idx[REF_PIC_LIST_1]].mv;

			mvp_bi[REF_PIC_LIST_0].hor_vector = mv_bi[REF_PIC_LIST_0].hor_vector - mv_best_cand_bi[REF_PIC_LIST_0].hor_vector;
			mvp_bi[REF_PIC_LIST_0].ver_vector = mv_bi[REF_PIC_LIST_0].ver_vector - mv_best_cand_bi[REF_PIC_LIST_0].ver_vector;
			mvp_bi[REF_PIC_LIST_1].hor_vector = mv_bi[REF_PIC_LIST_1].hor_vector - mv_best_cand_bi[REF_PIC_LIST_1].hor_vector;
			mvp_bi[REF_PIC_LIST_1].ver_vector = mv_bi[REF_PIC_LIST_1].ver_vector - mv_best_cand_bi[REF_PIC_LIST_1].ver_vector;

			if(currslice->mvd_l1_zero_flag)
			{
				int16_t *pred_other_y, *pred_other_u, *pred_other_v;
				int pred_other_y_stride, pred_other_ch_stride;
				
				motion_vector_t mv = mv_best_cand_bi[REF_PIC_LIST_1];

				ref_idx = curr_cu_info->inter_ref_index[REF_PIC_LIST_1];

				mv_bi[REF_PIC_LIST_1] = mv_best_cand_bi[REF_PIC_LIST_1];

				pred_other_y_stride = WND_STRIDE_2D(et->prediction_wnd[1+REF_PIC_LIST_1], Y_COMP);
				pred_other_y = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+REF_PIC_LIST_1], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				pred_other_ch_stride = WND_STRIDE_2D(et->prediction_wnd[1+REF_PIC_LIST_1], U_COMP);
				pred_other_u = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+REF_PIC_LIST_1], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
				pred_other_v = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+REF_PIC_LIST_1], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

				reference_wnd = &currslice->ref_pic_list[REF_PIC_LIST_1][ref_idx]->img;
				reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
				reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
				reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, U_COMP);
				reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
				reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);

				hmr_motion_compensation_luma(et, curr_cu_info, reference_buff_cu_position, reference_buff_stride, pred_other_y, pred_other_y_stride, curr_part_size, curr_part_size, curr_part_size_shift, &mv, 0);
				hmr_motion_compensation_chroma(et, reference_buff_cu_position_u, reference_buff_stride_chroma, pred_other_u, pred_other_ch_stride, curr_part_size_chroma, curr_part_size_shift_chroma, &mv, 0);
				hmr_motion_compensation_chroma(et, reference_buff_cu_position_v, reference_buff_stride_chroma, pred_other_v, pred_other_ch_stride, curr_part_size_chroma, curr_part_size_shift_chroma, &mv, 0);
			}

			if(best_cost[REF_PIC_LIST_0]>best_cost[REF_PIC_LIST_1] || currslice->mvd_l1_zero_flag)
				ref_list = REF_PIC_LIST_0;
			else
				ref_list = REF_PIC_LIST_1;

			if(!currslice->mvd_l1_zero_flag)
			{
				int16_t *pred_other_y, *pred_other_u, *pred_other_v;
				int pred_other_y_stride, pred_other_ch_stride;
				int ref_list_other = 1-ref_list;
				
				motion_vector_t mv = mv_bi[ref_list_other];

				ref_idx = curr_cu_info->inter_ref_index[ref_list_other];

//				mv_bi[ref_list_other] = mvp_bi[ref_list_other].mv;

				pred_other_y_stride = WND_STRIDE_2D(et->prediction_wnd[1+ref_list_other], Y_COMP);
				pred_other_y = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+ref_list_other], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				pred_other_ch_stride = WND_STRIDE_2D(et->prediction_wnd[1+ref_list_other], U_COMP);
				pred_other_u = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+ref_list_other], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
				pred_other_v = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+ref_list_other], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

				reference_wnd = &currslice->ref_pic_list[ref_list_other][ref_idx]->img;//[0] up to now we only use one reference	
				reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
				reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
				reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, U_COMP);
				reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
				reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);

				hmr_motion_compensation_luma(et, curr_cu_info, reference_buff_cu_position, reference_buff_stride, pred_other_y, pred_other_y_stride, curr_part_size, curr_part_size, curr_part_size_shift, &mv, 0);
				hmr_motion_compensation_chroma(et, reference_buff_cu_position_u, reference_buff_stride_chroma, pred_other_u, pred_other_ch_stride, curr_part_size_chroma, curr_part_size_shift_chroma, &mv, 0);
				hmr_motion_compensation_chroma(et, reference_buff_cu_position_v, reference_buff_stride_chroma, pred_other_v, pred_other_ch_stride, curr_part_size_chroma, curr_part_size_shift_chroma, &mv, 0);
			}


			num_ref_idx = currslice->num_ref_idx[ref_list];
			//for(iter=0;iter<num_iter;iter++)
			for (ref_idx = 0; ref_idx < num_ref_idx; ref_idx++)
			{
				int16_t *pred_other_y, *pred_other_u, *pred_other_v;
				int pred_other_y_stride, pred_other_ch_stride;
				int orig_buff_stride;
				int mvp_candidate_idx;
				mv_candiate_list_t	amvp_candidate_list;
				int16_t *aux_y, *aux_u, *aux_v;
				int aux_y_stride, aux_ch_stride;
				int list_other = 1-(int)ref_list;
//				int ref_idx = curr_cu_info->inter_ref_index[ref_list];
				motion_vector_t mv = curr_cu_info->inter_mv[ref_list];

				wnd_copy(et->funcs->sse_copy_16_16, &et->curr_mbs_wnd, &et->curr_mbs_aux_wnd);

				aux_y_stride = WND_STRIDE_2D(et->curr_mbs_aux_wnd, Y_COMP);
				aux_y = WND_POSITION_2D(int16_t *, et->curr_mbs_aux_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				aux_ch_stride = WND_STRIDE_2D(et->curr_mbs_aux_wnd, U_COMP);
				aux_u = WND_POSITION_2D(int16_t *, et->curr_mbs_aux_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
				aux_v = WND_POSITION_2D(int16_t *, et->curr_mbs_aux_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
				
				pred_other_y_stride = WND_STRIDE_2D(et->prediction_wnd[1+list_other], Y_COMP);
				pred_other_y = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+list_other], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				pred_other_ch_stride = WND_STRIDE_2D(et->prediction_wnd[1+list_other], U_COMP);
				pred_other_u = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+list_other], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
				pred_other_v = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+list_other], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

				remove_high_freq(pred_other_y, pred_other_y_stride, aux_y, aux_y_stride, curr_cu_info->size, curr_cu_info->size);
				remove_high_freq(pred_other_u, pred_other_ch_stride, aux_u, aux_ch_stride, curr_cu_info->size_chroma, curr_cu_info->size_chroma);
				remove_high_freq(pred_other_v, pred_other_ch_stride, aux_v, aux_ch_stride, curr_cu_info->size_chroma, curr_cu_info->size_chroma);			

				reference_wnd = &currslice->ref_pic_list[ref_list][ref_idx]->img;//[0] up to now we only use one reference	
				reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
				reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
				reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, U_COMP);
				reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
				reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);

				cost = hmr_bi_motion_estimation_HM(et, ctu, curr_cu_info, aux_y, aux_y_stride, reference_buff_cu_position, reference_buff_stride, curr_part_global_x, curr_part_global_y, 0, 0, curr_part_size, curr_part_size_shift, 4, 4, et->pict_width[Y_COMP], et->pict_height[Y_COMP], &mv);
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				//!!esto deberia guardarlo en un conjunto de listas para de candidatos y coger la lista de candidatos que queremos!!
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				get_amvp_candidates(et, currslice, ctu, curr_cu_info, &amvp_candidate_list, ref_list, ref_idx, part_size_type);//get candidates for advanced motion vector prediction from the neigbour CUs

				mv_cost = select_mv_candidate(et, curr_cu_info, &amvp_candidate_list, &mv, &mvp_candidate_idx);

				if (cost < best_cost_bi)
				{
					mv_bi[ref_list] = mv;
					ref_index_bi[ref_list] = ref_idx;
					mvp_bi[ref_list].hor_vector = mv.hor_vector-amvp_candidate_list.mv_candidates[mvp_candidate_idx].mv.hor_vector;
					mvp_bi[ref_list].ver_vector = mv.ver_vector-amvp_candidate_list.mv_candidates[mvp_candidate_idx].mv.ver_vector;
					mvp_idx[ref_list] = mvp_candidate_idx;
					best_cost_bi = cost;
					mv_best_cost_bi = mv_cost;
					//set mvs and ref_idx
				}
			}
		}

#ifdef COMPUTE_AS_HM
		if(best_cost_bi <= best_cost[REF_PIC_LIST_0] && best_cost_bi <= best_cost[REF_PIC_LIST_1])
#else
		if((best_cost_bi+mv_best_cost_bi) <= best_cost[REF_PIC_LIST_0]+mv_best_cost[REF_PIC_LIST_0] && best_cost_bi <= best_cost[REF_PIC_LIST_1]+mv_best_cost[REF_PIC_LIST_1])
#endif
		{
			sad+=best_cost_bi;
			mv_total_cost += mv_best_cost_bi;
			curr_cu_info->best_candidate_idx[REF_PIC_LIST_0] = mvp_idx[REF_PIC_LIST_0];
			curr_cu_info->best_dif_mv[REF_PIC_LIST_0] = mvp_bi[REF_PIC_LIST_0];
			curr_cu_info->best_candidate_idx[REF_PIC_LIST_1] = mvp_idx[REF_PIC_LIST_1];
			curr_cu_info->best_dif_mv[REF_PIC_LIST_1] = mvp_bi[REF_PIC_LIST_1];
			curr_cu_info->inter_mv[REF_PIC_LIST_0] = mv_bi[REF_PIC_LIST_0];
			curr_cu_info->inter_mv[REF_PIC_LIST_1] = mv_bi[REF_PIC_LIST_1];
//			curr_cu_info->inter_ref_list = //REF_PIC_LIST_0 | REF_PIC_LIST_1;
			curr_cu_info->inter_ref_index[REF_PIC_LIST_0] =  ref_index_bi[REF_PIC_LIST_0];
			curr_cu_info->inter_ref_index[REF_PIC_LIST_1] =  ref_index_bi[REF_PIC_LIST_1];
			curr_cu_info->inter_mode = (0x1<<REF_PIC_LIST_0) | (0x1<<REF_PIC_LIST_1);
		}
#ifdef COMPUTE_AS_HM
		else if(best_cost[REF_PIC_LIST_0]<=best_cost[REF_PIC_LIST_1])
#else
		else if(best_cost[REF_PIC_LIST_0]+mv_best_cost[REF_PIC_LIST_0]<=best_cost[REF_PIC_LIST_1]+mv_best_cost[REF_PIC_LIST_1])
#endif
		{
			sad+=best_cost[REF_PIC_LIST_0];
			mv_total_cost += mv_best_cost[REF_PIC_LIST_0];
			curr_cu_info->inter_ref_index[REF_PIC_LIST_1] = -1;					
//			curr_cu_info->inter_ref_list = REF_PIC_LIST_0;
			curr_cu_info->inter_mode = (0x1<<REF_PIC_LIST_0);
		}
		else
		{
			sad+=best_cost[REF_PIC_LIST_1];
			mv_total_cost += mv_best_cost[REF_PIC_LIST_1];
			curr_cu_info->inter_ref_index[REF_PIC_LIST_0] = -1;							
//			curr_cu_info->inter_ref_list = REF_PIC_LIST_1;
			curr_cu_info->inter_mode = (0x1<<REF_PIC_LIST_1);
		}

		if(part_size_type == SIZE_NxN)
		{
			SET_INTER_MV_BUFFS(et, ctu, curr_cu_info, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);
			memset(&ctu->mv_ref_idx[REF_PIC_LIST_0][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_0], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));	
			memset(&ctu->mv_ref_idx[REF_PIC_LIST_1][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_1], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));	
			memset(&ctu->inter_mode[curr_cu_info->abs_index], curr_cu_info->inter_mode, curr_cu_info->num_part_in_cu*sizeof(ctu->inter_mode[0]));	
			memset(&ctu->pred_mode[curr_cu_info->abs_index], INTER_MODE, curr_cu_info->num_part_in_cu*sizeof(ctu->pred_mode[0]));				
		}

//		SET_INTER_MV_BUFFS(et, ctu, curr_cu_info, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);
//		memset(&ctu->mv_ref_idx[REF_PIC_LIST_0][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_0], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));
		curr_cu_info++;
	}
	return sad + 1*mv_total_cost;//.5*mv_cost;
}


int check_unidirectional_motion(slice_t *currslice, cu_partition_info_t *curr_cu_info)
{
	if(currslice->slice_type != B_SLICE || curr_cu_info->inter_mode != 3)
		return TRUE;

	//if( pcCU->getSlice()->isInterB() && !pcCU->getSlice()->getPPS()->getWPBiPred() )
	if(currslice->slice_type == B_SLICE && !currslice->pps->weighted_bipred_flag)
	{
		if( curr_cu_info->inter_ref_index[REF_PIC_LIST_0] >= 0 && curr_cu_info->inter_ref_index[REF_PIC_LIST_1] >= 0)
		{
			//currslice->ref_pic_list[ref_pic_list][ref_idx]->temp_info.poc
			int RefPOCL0 = currslice->ref_pic_list[REF_PIC_LIST_0][curr_cu_info->inter_ref_index[REF_PIC_LIST_0]]->temp_info.poc; //pcCU->getSlice()->getRefPic(REF_PIC_LIST_0, pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(PartAddr))->getPOC();
			int RefPOCL1 = currslice->ref_pic_list[REF_PIC_LIST_1][curr_cu_info->inter_ref_index[REF_PIC_LIST_1]]->temp_info.poc; //pcCU->getSlice()->getRefPic(REF_PIC_LIST_1, pcCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx(PartAddr))->getPOC();
			//if(RefPOCL0 == RefPOCL1 && pcCU->getCUMvField(REF_PIC_LIST_0)->getMv(PartAddr) == pcCU->getCUMvField(REF_PIC_LIST_1)->getMv(PartAddr))
			if(RefPOCL0 == RefPOCL1 && curr_cu_info->inter_mv[REF_PIC_LIST_0].hor_vector == curr_cu_info->inter_mv[REF_PIC_LIST_1].hor_vector && curr_cu_info->inter_mv[REF_PIC_LIST_0].ver_vector == curr_cu_info->inter_mv[REF_PIC_LIST_1].ver_vector)
			{
				return TRUE;
			}
		}
	}
	return FALSE;
}

void weighted_average_motion(int16_t* src0, int src0_stride, int16_t* src1, int src1_stride, int16_t* dst, int dst_stride, int height, int width, int bit_depth)
{
  int x, y;
  int max_pix_val = (1<<bit_depth)-1;
  int shiftNum = IF_INTERNAL_PREC + 1 - bit_depth;
  int offset = ( 1 << ( shiftNum - 1 ) ) + 2 * IF_INTERNAL_OFFS;
  
  for ( y = 0; y < height; y++ )
  {
    for ( x = 0; x < width; x ++)
    {
      dst[x] = clip(( ( src0[x] + src1[x] + offset ) >> shiftNum ),0,max_pix_val);
    }
    src0 += src0_stride;
    src1 += src1_stride;
    dst +=dst_stride;
  }
}



int predict_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type)
{
	double distortion = 0.;

	slice_t *currslice = &et->enc_engine->current_pict.slice;
	ctu_info_t *ctu_rd = et->ctu_rd;
	int pred_buff_stride, orig_buff_stride, reference_buff_stride, residual_buff_stride;
	int pred_buff_stride_chroma, orig_buff_stride_chroma, reference_buff_stride_chroma, residual_buff_stride_chroma;
	int16_t  *orig_buff, *orig_buff_u, *orig_buff_v;
	int16_t  *reference_buff_cu_position, *reference_buff_cu_position_u, *reference_buff_cu_position_v;
	int16_t *pred_buff, *pred_buff_u, *pred_buff_v, *residual_buff, *residual_buff_u, *residual_buff_v;
	int16_t *pred_aux_buffs[NUM_REF_LISTS][NUM_PICT_COMPONENTS];
	wnd_t *reference_wnd=NULL;//, *resi_wnd = NULL;
	int curr_part_size, curr_part_size_shift;
	int curr_part_size_chroma, curr_part_size_shift_chroma;
	int curr_part_x, curr_part_y, curr_part_global_x, curr_part_global_y;
	int curr_part_x_chroma, curr_part_y_chroma, curr_part_global_x_chroma, curr_part_global_y_chroma;
	int curr_depth = depth;
	cu_partition_info_t *parent_part_info;
	cu_partition_info_t *curr_cu_info;
	int num_partitions, npart;
	int mv_cost = 0;
	
	curr_cu_info = &ctu->partition_list[et->partition_depth_start[curr_depth]]+part_position;

	if(part_size_type == SIZE_2Nx2N)
	{
		parent_part_info = curr_cu_info->parent;	
		num_partitions = 1;
	}
	else if(part_size_type == SIZE_NxN)
	{
		parent_part_info = curr_cu_info->parent;
		curr_cu_info = parent_part_info->children[0];
		num_partitions = 4;
	}

	for(npart=0;npart<num_partitions;npart++)
	{
		//set mvs and ref_idx
		motion_vector_t mv;

		int ref_list = REF_PIC_LIST_0, list_init = REF_PIC_LIST_0, list_last = REF_PIC_LIST_1;
		int uni_prediction = check_unidirectional_motion(currslice, curr_cu_info);

		if(uni_prediction)
		{
			list_init = list_last = ref_list = curr_cu_info->inter_mode==0x1?REF_PIC_LIST_0:REF_PIC_LIST_1;
		}

		curr_depth = curr_cu_info->depth;
		curr_part_x = curr_cu_info->x_position;
		curr_part_y = curr_cu_info->y_position;
		curr_part_global_x = ctu->x[Y_COMP]+curr_part_x;
		curr_part_global_y = ctu->y[Y_COMP]+curr_part_y;
		curr_part_size = curr_cu_info->size;
		curr_part_size_shift = et->max_cu_size_shift-curr_depth;
		curr_part_x_chroma = curr_cu_info->x_position_chroma;
		curr_part_y_chroma = curr_cu_info->y_position_chroma;
		curr_part_global_x_chroma = ctu->x[CHR_COMP]+curr_part_x_chroma;
		curr_part_global_y_chroma = ctu->y[CHR_COMP]+curr_part_y_chroma;
		curr_part_size_chroma = curr_cu_info->size_chroma;
		curr_part_size_shift_chroma = et->max_cu_size_shift-curr_depth-1;//420

		pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd[0], Y_COMP);
		pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
		pred_buff_stride_chroma = WND_STRIDE_2D(et->prediction_wnd[0], CHR_COMP);
		pred_buff_u = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
		pred_buff_v = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

		for (ref_list = list_init; ref_list <= list_last; ref_list++)
		{
			int ref_idx = curr_cu_info->inter_ref_index[ref_list];
			mv = curr_cu_info->inter_mv[ref_list];

			if(uni_prediction)
			{
				pred_aux_buffs[ref_list][Y_COMP] = pred_buff;
				pred_aux_buffs[ref_list][U_COMP] = pred_buff_u;
				pred_aux_buffs[ref_list][V_COMP] = pred_buff_v;
			}
			else
			{
				pred_aux_buffs[ref_list][Y_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+ref_list], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
				pred_aux_buffs[ref_list][U_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+ref_list], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
				pred_aux_buffs[ref_list][V_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[1+ref_list], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);			
			}

			orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
			orig_buff = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
			orig_buff_stride_chroma = WND_STRIDE_2D(et->curr_mbs_wnd, CHR_COMP);
			orig_buff_u = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
			orig_buff_v = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

			residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, Y_COMP);
			residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
			residual_buff_stride_chroma = WND_STRIDE_2D(et->residual_wnd, CHR_COMP);
			residual_buff_u = WND_POSITION_2D(int16_t *, et->residual_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
			residual_buff_v = WND_POSITION_2D(int16_t *, et->residual_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

			reference_wnd = &currslice->ref_pic_list[ref_list][ref_idx]->img;//[0] up to now we only use one reference	
			reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
			reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
			reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, CHR_COMP);
			reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
			reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);

			get_amvp_candidates(et, currslice, ctu, curr_cu_info, &et->amvp_candidates[ref_list], ref_list, ref_idx, part_size_type);//get candidates for motion search from the neigbour CUs
#ifdef COMPUTE_AS_HM
			mv_cost += select_mv_candidate(et, curr_cu_info, &et->amvp_candidates[ref_list], &mv, &curr_cu_info->best_candidate_idx[ref_list]);
#else
			mv_cost += select_mv_candidate(et, curr_cu_info, &et->amvp_candidates[ref_list], &mv, &curr_cu_info->best_candidate_idx[ref_list]);
#endif
			curr_cu_info->best_dif_mv[ref_list].hor_vector = mv.hor_vector-et->amvp_candidates[ref_list].mv_candidates[curr_cu_info->best_candidate_idx[ref_list]].mv.hor_vector;
			curr_cu_info->best_dif_mv[ref_list].ver_vector = mv.ver_vector-et->amvp_candidates[ref_list].mv_candidates[curr_cu_info->best_candidate_idx[ref_list]].mv.ver_vector;


			//##################################################### this should be moved somewhere in an upper herarchy ##################################################
			SET_INTER_MV_BUFFS(et, ctu, curr_cu_info, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);
			memset(&ctu->mv_ref_idx[REF_PIC_LIST_0][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_0], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));
			memset(&ctu->mv_ref_idx[REF_PIC_LIST_1][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_1], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));
			//##########################################################################################################################################################################

			hmr_motion_compensation_luma(et, curr_cu_info, reference_buff_cu_position, reference_buff_stride, pred_aux_buffs[ref_list][Y_COMP], pred_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, &mv, !uni_prediction);
			hmr_motion_compensation_chroma(et, reference_buff_cu_position_u, reference_buff_stride_chroma, pred_aux_buffs[ref_list][U_COMP], pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv, !uni_prediction);
			hmr_motion_compensation_chroma(et, reference_buff_cu_position_v, reference_buff_stride_chroma, pred_aux_buffs[ref_list][V_COMP], pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv, !uni_prediction);
		}

		if(!uni_prediction)
		{
			et->funcs->weighted_average_motion(pred_aux_buffs[REF_PIC_LIST_0][Y_COMP], pred_buff_stride, pred_aux_buffs[REF_PIC_LIST_1][Y_COMP], pred_buff_stride, pred_buff, pred_buff_stride, curr_cu_info->size, curr_cu_info->size, et->bit_depth);
			et->funcs->weighted_average_motion(pred_aux_buffs[REF_PIC_LIST_0][U_COMP], pred_buff_stride_chroma, pred_aux_buffs[REF_PIC_LIST_1][U_COMP], pred_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, curr_cu_info->size_chroma, curr_cu_info->size_chroma, et->bit_depth);
			et->funcs->weighted_average_motion(pred_aux_buffs[REF_PIC_LIST_0][V_COMP], pred_buff_stride_chroma, pred_aux_buffs[REF_PIC_LIST_1][V_COMP], pred_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, curr_cu_info->size_chroma, curr_cu_info->size_chroma, et->bit_depth);
		}

		et->funcs->predict(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride, residual_buff, residual_buff_stride, curr_part_size);
		et->funcs->predict(orig_buff_u, orig_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, residual_buff_u, residual_buff_stride_chroma, curr_part_size_chroma);
		et->funcs->predict(orig_buff_v, orig_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, residual_buff_v, residual_buff_stride_chroma, curr_part_size_chroma);

		curr_cu_info++;
	}

	return mv_cost;
}


//this function is referenced by the initial depth, not by the processing depth
int encode_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type)
{
	double cost;
	slice_t *currslice = &et->enc_engine->current_pict.slice;
	int curr_depth = depth;
	cu_partition_info_t *parent_part_info;
	cu_partition_info_t *curr_cu_info;
	int curr_sum_y = 0, curr_sum_u = 0, curr_sum_v = 0;
	int cu_min_tu_size_shift;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int max_tr_depth, max_tr_processing_depth;
	int initial_state, end_state;
//	int cbf_split[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
//	int acc_cost[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int log2cu_size;
	int ref_idx = 0;
	int tu_log_max_size = et->max_tu_size_shift;
	int qp;

//	curr_depth = curr_cu_info->depth;
	memset(depth_state, 0, sizeof(depth_state));


	if(depth==0 && et->max_cu_size == MAX_CU_SIZE)
	{
		parent_part_info = &ctu->partition_list[et->partition_depth_start[depth]];
		curr_cu_info = parent_part_info->children[0];
		parent_part_info->cost = INT_MAX;//the parent is consolidated if 2Nx2N
		initial_state = part_position & 0x3;//initial_state = 0;
		end_state = initial_state;//end_state = 1;
		qp = parent_part_info->qp;
	}
	else
	{
		curr_cu_info = &ctu->partition_list[et->partition_depth_start[curr_depth]]+part_position;
		parent_part_info = curr_cu_info->parent;
		initial_state = part_position & 0x3;
		end_state = initial_state+1;
		qp = curr_cu_info->qp;
	}
	
	curr_depth = curr_cu_info->depth;

	max_tr_depth = et->max_inter_tr_depth;

	log2cu_size = et->max_cu_size_shift-(depth);//-(part_size_type==SIZE_NxN));

	if(log2cu_size < et->min_tu_size_shift+et->max_inter_tr_depth-1 + (et->max_inter_tr_depth==1 && part_size_type!=SIZE_2Nx2N))
		cu_min_tu_size_shift = et->min_tu_size_shift;
	else
	{
		cu_min_tu_size_shift = log2cu_size - (max_tr_depth - 1 + (et->max_inter_tr_depth==1 && part_size_type!=SIZE_2Nx2N));
		if(cu_min_tu_size_shift>tu_log_max_size)
			cu_min_tu_size_shift=tu_log_max_size;
	}




#ifdef COMPUTE_AS_HM
	max_tr_processing_depth = et->max_cu_size_shift-cu_min_tu_size_shift;
#else
	max_tr_processing_depth = et->max_cu_size_shift-cu_min_tu_size_shift;
	if(et->performance_mode >= PERF_FAST_COMPUTATION)
	{
		int aux_proc_depth = depth==0?1:(depth+((part_size_type!=SIZE_2Nx2N)?1:0));
		max_tr_processing_depth = aux_proc_depth;//max_tr_processing_depth = (depth+2<=max_tr_processing_depth)?depth+2:((depth+1<=max_tr_processing_depth)?depth+1:max_tr_processing_depth);
	}
#endif

	//skip first level 
	if(et->max_inter_tr_depth==1 && part_size_type != SIZE_2Nx2N && curr_depth==depth && log2cu_size>max_tr_processing_depth)
	{
		parent_part_info = curr_cu_info;
		curr_cu_info = parent_part_info->children[0];
		parent_part_info->distortion = parent_part_info->cost = MAX_COST;
		initial_state = part_position & 0x3;
		end_state = initial_state;//+1;
	}

	memset(depth_state, 0, sizeof(depth_state));
	depth_state[curr_depth] = initial_state;

	curr_depth = curr_cu_info->depth;

	while(curr_depth!=depth || depth_state[curr_depth]!=end_state)//transform loop
	{
//		int fast_end_loop = FALSE;
		uint bit_cost = 0;
		uint dist_y, dist_u, dist_v;
		curr_cu_info = (parent_part_info==NULL)?curr_cu_info:parent_part_info->children[depth_state[curr_depth]];//if cu_size=64 we process 4 32x32 partitions, else just the curr_partition
		curr_cu_info->qp = qp;
		curr_depth = curr_cu_info->depth;

		dist_y = encode_inter_cu(et, ctu, curr_cu_info, depth, part_size_type, &curr_sum_y, gcnt);//depth = prediction depth

		if(curr_cu_info->size_chroma!=2 || (curr_cu_info->size_chroma==2 && depth_state[curr_depth]==0))
		{
			dist_u = encode_inter_cu_chroma(et, ctu, curr_cu_info, U_COMP, depth, part_size_type, &curr_sum_u, gcnt);//depth = prediction depth
			dist_v = encode_inter_cu_chroma(et, ctu, curr_cu_info, V_COMP, depth, part_size_type, &curr_sum_v, gcnt);//depth = prediction depth
		}
		else
		{
			dist_u = dist_v = 0;
			curr_cu_info->inter_cbf[U_COMP] = (curr_cu_info-1)->inter_cbf[U_COMP];
			curr_cu_info->inter_cbf[V_COMP] = (curr_cu_info-1)->inter_cbf[V_COMP];
		}
		curr_cu_info->distortion = dist_y+dist_u+dist_v;
		curr_cu_info->cost = curr_cu_info->distortion;
		curr_cu_info->sum = curr_sum_y+curr_sum_u+curr_sum_v;

		depth_state[curr_depth]++;

/*		if(curr_depth < max_tr_processing_depth && curr_sum_y==0 && curr_sum_u==0 && curr_sum_v==0)
		{
			fast_end_loop = TRUE;
//			if(curr_depth==depth || depth_state[curr_depth]==end_state)
			{
				SET_ENC_INFO_BUFFS(et, curr_cu_info, depth+(part_size_type!=SIZE_2Nx2N), curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);
				synchronize_motion_buffers_luma(et, curr_cu_info, et->transform_quant_wnd[curr_depth+(part_size_type!=SIZE_2Nx2N)], et->transform_quant_wnd[depth+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[depth+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
				synchronize_motion_buffers_chroma(et, curr_cu_info, et->transform_quant_wnd[curr_depth+(part_size_type!=SIZE_2Nx2N)], et->transform_quant_wnd[depth+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[depth+1+(part_size_type!=SIZE_2Nx2N)], gcnt);

			}
		}
*/
		//HM order
		if(curr_depth < max_tr_processing_depth)// && !fast_end_loop)//go one level down
		{
			curr_depth++;
			parent_part_info = curr_cu_info;
		}
		else if(depth_state[curr_depth]==4)//consolidate 
		{	
			while(depth_state[curr_depth]==4 && (curr_depth > (depth)))
			{
				double distortion =  parent_part_info->children[0]->distortion+parent_part_info->children[1]->distortion+parent_part_info->children[2]->distortion+parent_part_info->children[3]->distortion;
				uint sum =  parent_part_info->children[0]->sum+parent_part_info->children[1]->sum+parent_part_info->children[2]->sum+parent_part_info->children[3]->sum;
//				distortion = parent_part_info->children[0]->distortion+parent_part_info->children[1]->distortion+parent_part_info->children[2]->distortion+parent_part_info->children[3]->distortion;
				cost = distortion;

				depth_state[curr_depth] = 0;

				if(cost < parent_part_info->cost)// && (cbf_y||cbf_u||cbf_v))
				{
					parent_part_info->cost = (uint32_t) cost;
					parent_part_info->distortion = (uint32_t)distortion;
					parent_part_info->sum = sum;

					if(curr_depth == max_tr_processing_depth)	//create cbf and tr_idx buffs
					{
						int nchild;
//						uint cbf_split[NUM_PICT_COMPONENTS];
						uint cbf_split_y, cbf_split_u, cbf_split_v;
						int tr_mask = 0x1<<(curr_depth-depth);

						cbf_split_y = (parent_part_info->children[0]->inter_cbf[Y_COMP]&tr_mask)|(parent_part_info->children[1]->inter_cbf[Y_COMP]&tr_mask)|(parent_part_info->children[2]->inter_cbf[Y_COMP]&tr_mask)|(parent_part_info->children[3]->inter_cbf[Y_COMP]&tr_mask);
						cbf_split_u = (parent_part_info->children[0]->inter_cbf[U_COMP]&tr_mask)|(parent_part_info->children[1]->inter_cbf[U_COMP]&tr_mask)|(parent_part_info->children[2]->inter_cbf[U_COMP]&tr_mask)|(parent_part_info->children[3]->inter_cbf[U_COMP]&tr_mask);
						cbf_split_v = (parent_part_info->children[0]->inter_cbf[V_COMP]&tr_mask)|(parent_part_info->children[1]->inter_cbf[V_COMP]&tr_mask)|(parent_part_info->children[2]->inter_cbf[V_COMP]&tr_mask)|(parent_part_info->children[3]->inter_cbf[V_COMP]&tr_mask);

						cbf_split_y>>=1;
						cbf_split_u>>=1;
						cbf_split_v>>=1;
						for(nchild=0;nchild<4;nchild++)
						{
							cu_partition_info_t *cu_info = parent_part_info->children[nchild];
							cu_info->inter_cbf[Y_COMP]|=cbf_split_y;
							cu_info->inter_cbf[U_COMP]|=cbf_split_u;
							cu_info->inter_cbf[V_COMP]|=cbf_split_v;
							SET_ENC_INFO_BUFFS(et, cu_info, depth+(part_size_type!=SIZE_2Nx2N), cu_info->abs_index, cu_info->num_part_in_cu);//consolidate in prediction depth
						}
					}
					else
					{
						int ll;
						int buff_depth = depth+(part_size_type!=SIZE_2Nx2N);
						int abs_index = parent_part_info->abs_index;
						int num_part_in_cu = parent_part_info->num_part_in_cu;
//						uint cbf_split[NUM_PICT_COMPONENTS];
						uint cbf_y, cbf_u, cbf_v;
						int tr_mask = 0x1<<(curr_depth-depth);
						cbf_y = (et->cbf_buffs[Y_COMP][buff_depth][parent_part_info->children[0]->abs_index]&tr_mask)|(et->cbf_buffs[Y_COMP][buff_depth][parent_part_info->children[1]->abs_index]&tr_mask)|(et->cbf_buffs[Y_COMP][buff_depth][parent_part_info->children[2]->abs_index]&tr_mask)|(et->cbf_buffs[Y_COMP][buff_depth][parent_part_info->children[3]->abs_index]&tr_mask);
						cbf_u = (et->cbf_buffs[U_COMP][buff_depth][parent_part_info->children[0]->abs_index]&tr_mask)|(et->cbf_buffs[U_COMP][buff_depth][parent_part_info->children[1]->abs_index]&tr_mask)|(et->cbf_buffs[U_COMP][buff_depth][parent_part_info->children[2]->abs_index]&tr_mask)|(et->cbf_buffs[U_COMP][buff_depth][parent_part_info->children[3]->abs_index]&tr_mask);
						cbf_v = (et->cbf_buffs[V_COMP][buff_depth][parent_part_info->children[0]->abs_index]&tr_mask)|(et->cbf_buffs[V_COMP][buff_depth][parent_part_info->children[1]->abs_index]&tr_mask)|(et->cbf_buffs[V_COMP][buff_depth][parent_part_info->children[2]->abs_index]&tr_mask)|(et->cbf_buffs[V_COMP][buff_depth][parent_part_info->children[3]->abs_index]&tr_mask);

						cbf_y>>=1;//<<= (curr_depth-depth-1);
						cbf_u>>=1;//<<= (curr_depth-depth-1);
						cbf_v>>=1;//<<= (curr_depth-depth-1);
						cbf_y |= ((cbf_y&0xff)<<8) | ((cbf_y&0xff)<<16) | ((cbf_y&0xff)<<24);
						cbf_u |= ((cbf_u&0xff)<<8) | ((cbf_u&0xff)<<16) | ((cbf_u&0xff)<<24);
						cbf_v |= ((cbf_v&0xff)<<8) | ((cbf_v&0xff)<<16) | ((cbf_v&0xff)<<24);
						//consolidate cbf flags
						for(ll=(abs_index>>2);ll<((abs_index+num_part_in_cu)>>2);ll++)
						{
							((uint*)(et->cbf_buffs[Y_COMP][buff_depth]))[ll] |= cbf_y;
							((uint*)(et->cbf_buffs[U_COMP][buff_depth]))[ll] |= cbf_u;
							((uint*)(et->cbf_buffs[V_COMP][buff_depth]))[ll] |= cbf_v;
						}											
					}
					synchronize_motion_buffers_luma(et, parent_part_info, et->transform_quant_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], et->transform_quant_wnd[curr_depth-1+1+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[curr_depth-1+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
					synchronize_motion_buffers_chroma(et, parent_part_info, et->transform_quant_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], et->transform_quant_wnd[curr_depth-1+1+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[curr_depth-1+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
				}
				else 
				{
					SET_ENC_INFO_BUFFS(et, parent_part_info, depth+(part_size_type!=SIZE_2Nx2N), parent_part_info->abs_index, parent_part_info->num_part_in_cu);
//					synchronize_motion_buffers_luma(et, parent_part_info, et->transform_quant_wnd[curr_depth-1+(part_size_type!=SIZE_2Nx2N)], et->transform_quant_wnd[depth+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[curr_depth+1-1+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[depth+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
//					synchronize_motion_buffers_chroma(et, parent_part_info, et->transform_quant_wnd[curr_depth-1+(part_size_type!=SIZE_2Nx2N)], et->transform_quant_wnd[depth+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[curr_depth+1-1+(part_size_type!=SIZE_2Nx2N)], et->decoded_mbs_wnd[depth+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
				}

				curr_depth--;
				parent_part_info = parent_part_info->parent;
			}
		}
	}

	curr_cu_info = &ctu->partition_list[et->partition_depth_start[depth]]+part_position;
	if(depth==max_tr_processing_depth)
	{
		SET_ENC_INFO_BUFFS(et, curr_cu_info, depth+(part_size_type!=SIZE_2Nx2N), curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);	
	}

	memset(&ctu->mv_ref_idx[REF_PIC_LIST_0][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_0], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[REF_PIC_LIST_0][0]));

	return curr_cu_info->cost;
}



#define CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, depth, abs_idx, num_partitions)																	\
{																																			\
	memcpy(&ctu->cbf[Y_COMP][abs_idx], &et->cbf_buffs[Y_COMP][depth][abs_idx], num_partitions*sizeof(et->cbf_buffs[0][0][0]));				\
	memcpy(&ctu->cbf[U_COMP][abs_idx], &et->cbf_buffs[U_COMP][depth][abs_idx], num_partitions*sizeof(et->cbf_buffs[0][0][0]));				\
	memcpy(&ctu->cbf[V_COMP][abs_idx], &et->cbf_buffs[V_COMP][depth][abs_idx], num_partitions*sizeof(et->cbf_buffs[0][0][0]));				\
	memcpy(&ctu->tr_idx[abs_idx], &et->tr_idx_buffs[depth][abs_idx], num_partitions*sizeof(et->tr_idx_buffs[0][0]));						\
	memcpy(&ctu->intra_mode[Y_COMP][abs_idx], &et->intra_mode_buffs[Y_COMP][depth][abs_idx], num_partitions*sizeof(et->intra_mode_buffs[0][0][0]));		\
	memcpy(&ctu->intra_mode[CHR_COMP][abs_idx], &et->intra_mode_buffs[CHR_COMP][depth][abs_idx], num_partitions*sizeof(et->intra_mode_buffs[0][0][0]));	\
}


void SET_INTER_INFO_BUFFS(henc_thread_t *et, ctu_info_t *ctu, cu_partition_info_t *cu_info, int abs_idx, int num_partitions)//, int ref_list)
{																																							
	int i, ref_list;																																					
	if(cu_info->prediction_mode == INTER_MODE)																												
	{																																						
		memset(&ctu->inter_mode[abs_idx], cu_info->inter_mode, num_partitions*sizeof(ctu->inter_mode[0]));
		memset(&ctu->skipped[abs_idx], cu_info->skipped, num_partitions*sizeof(ctu->skipped[0]));															
		memset(&ctu->merge[abs_idx], cu_info->merge_flag, num_partitions*sizeof(ctu->merge[0]));															
		memset(&ctu->merge_idx[abs_idx], cu_info->merge_idx, num_partitions*sizeof(ctu->merge_idx[0]));													
		for(ref_list=REF_PIC_LIST_0;ref_list<=REF_PIC_LIST_1;ref_list++)
		{
			memset(&ctu->mv_ref_idx[ref_list][abs_idx], cu_info->inter_ref_index[ref_list], num_partitions*sizeof(ctu->mv_ref_idx[0][0]));
			if(cu_info->inter_mode & (0x1<<ref_list))
			{
				ctu->mv_diff_ref_idx[ref_list][abs_idx] = cu_info->best_candidate_idx[ref_list];
				ctu->mv_diff[ref_list][abs_idx] = cu_info->best_dif_mv[ref_list];														
				
				for(i=0;i<num_partitions;i++)																														
				{
					ctu->mv_ref[ref_list][abs_idx+i] = cu_info->inter_mv[ref_list];
				}
			}
		}
	}																																						
	else																																					
	{																																						
		memset(&ctu->mv_ref_idx[0][abs_idx], -1, num_partitions*sizeof(ctu->mv_ref_idx[0][0]));
		memset(&ctu->mv_ref_idx[1][abs_idx], -1, num_partitions*sizeof(ctu->mv_ref_idx[1][0]));
		memset(&ctu->skipped[abs_idx], 0, num_partitions*sizeof(ctu->skipped[0]));																			
		memset(&ctu->merge[abs_idx], 0, num_partitions*sizeof(ctu->merge[0]));															
	}																																						
	memset(&ctu->pred_mode[abs_idx], cu_info->prediction_mode, num_partitions*sizeof(ctu->pred_mode[0]));													
}


#define SET_INTER_INFO_BUFFS_(et, ctu, cu_info, abs_idx, num_partitions, ref_list)																			\
{																																							\
	int i;																																					\
	if(cu_info->prediction_mode == INTER_MODE)																												\
	{																																						\
		ctu->mv_diff_ref_idx[ref_list][abs_idx] = cu_info->best_candidate_idx[ref_list];																	\
		ctu->mv_diff[ref_list][abs_idx] = cu_info->best_dif_mv[ref_list];																					\
		memset(&ctu->inter_mode[abs_idx], 1<<ref_list, num_partitions*sizeof(ctu->inter_mode[0]));															\
		memset(&ctu->mv_ref_idx[ref_list][abs_idx], cu_info->inter_ref_index[ref_list], num_partitions*sizeof(ctu->mv_ref_idx[0][0]));					\
		memset(&ctu->skipped[abs_idx], cu_info->skipped, num_partitions*sizeof(ctu->skipped[0]));															\
		memset(&ctu->merge[abs_idx], cu_info->merge_flag, num_partitions*sizeof(ctu->merge[0]));															\
		memset(&ctu->merge_idx[abs_idx], cu_info->merge_idx, num_partitions*sizeof(ctu->merge_idx[0]));													\
		for(i=0;i<num_partitions;i++)																														\
		{																																					\
			ctu->mv_ref[ref_list][abs_idx+i] = cu_info->inter_mv[ref_list];																					\
		}																																					\
	}																																						\
	else																																					\
	{																																						\
		memset(&ctu->mv_ref_idx[ref_list][abs_idx], -1, num_partitions*sizeof(ctu->mv_ref_idx[0][0]));														\
		memset(&ctu->skipped[abs_idx], 0, num_partitions*sizeof(ctu->skipped[0]));																			\
		memset(&ctu->merge[abs_idx], 0, num_partitions*sizeof(ctu->merge[0]));															\
	}																																						\
	memset(&ctu->pred_mode[abs_idx], cu_info->prediction_mode, num_partitions*sizeof(ctu->pred_mode[0]));													\
}																																							



void consolidate_prediction_info(henc_thread_t *et, ctu_info_t *ctu, ctu_info_t *ctu_rd, cu_partition_info_t *parent_part_info, uint32_t parent_cost, uint32_t children_cost, int is_max_depth, uint32_t *cost_sum)
{
	int abs_index = parent_part_info->abs_index;
	int num_part_in_cu = parent_part_info->num_part_in_cu;
	int curr_depth = parent_part_info->depth + 1;
	int gcnt = 0;
	int parent_sum = parent_part_info->sum;
	int children_sum = parent_part_info->children[0]->sum+parent_part_info->children[1]->sum+parent_part_info->children[2]->sum+parent_part_info->children[3]->sum;

	//choose best
#ifdef COMPUTE_AS_HM
	if((children_cost<parent_cost || (/*et->enc_engine->current_pict.slice.slice_type != I_SLICE && */children_cost==parent_cost && is_max_depth && parent_part_info->prediction_mode == INTRA_MODE && parent_part_info->children[0]->prediction_mode == INTER_MODE)) || !(parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame))//if we get here, tl should be inside the frame
#else
	if(children_cost<parent_cost || !(parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame))//if we get here, tl should be inside the frame
//	if(((et->enc_engine->current_pict.slice.slice_type == I_SLICE && children_cost<parent_cost) || 
//		(et->enc_engine->current_pict.slice.slice_type != I_SLICE && (children_cost+clip(et->enc_engine->avg_dist/1.75,5.,20000.)*children_sum<parent_cost+clip(et->enc_engine->avg_dist/1.75,5.,20000.)*parent_sum))) || 
//		!(parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame))//if we get here, tl should be inside the frame
//	if((et->rd_mode != RD_FAST && (children_cost < parent_cost)) || (et->rd_mode == RD_FAST && ((children_cost+45*children_sum)<(parent_cost+45*parent_part_info->sum))) || !(parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame))//if we get here, tl should be inside the frame
#endif
	{
		//here we consolidate the bottom-up results being preferred to the top-down computation
		PartSize part_size_type2 = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//

		if(cost_sum!=NULL)
		{
			cost_sum[parent_part_info->depth] -= parent_cost;
			cost_sum[parent_part_info->depth] += children_cost;
		}
		parent_part_info->cost = children_cost;
		parent_part_info->distortion = parent_part_info->children[0]->distortion+parent_part_info->children[1]->distortion+parent_part_info->children[2]->distortion+parent_part_info->children[3]->distortion;
		parent_part_info->sum = children_sum;//parent_part_info->children[0]->sum+parent_part_info->children[1]->sum+parent_part_info->children[2]->sum+parent_part_info->children[3]->sum;

		if(is_max_depth)
		{
			int nchild;
			int num_part_in_sub_cu = parent_part_info->children[0]->num_part_in_cu;

			synchronize_motion_buffers_luma(et, parent_part_info, et->transform_quant_wnd[curr_depth+1], et->transform_quant_wnd[0], et->decoded_mbs_wnd[curr_depth+1], et->decoded_mbs_wnd[0], gcnt);
			synchronize_motion_buffers_chroma(et, parent_part_info, et->transform_quant_wnd[curr_depth+1], et->transform_quant_wnd[0], et->decoded_mbs_wnd[curr_depth+1], et->decoded_mbs_wnd[0], gcnt);

			CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu)

			for(nchild=0;nchild<4;nchild++)
			{
				cu_partition_info_t *cu_info = parent_part_info->children[nchild];
				SET_INTER_INFO_BUFFS(et, ctu, cu_info, cu_info->abs_index, cu_info->num_part_in_cu);//, REF_PIC_LIST_0);
				memset(&ctu->qp[cu_info->abs_index], cu_info->qp, cu_info->num_part_in_cu*sizeof(ctu->qp[0]));
			}

			//if we fill this in here we don't have to consolidate
			memset(&ctu->pred_depth[abs_index], curr_depth-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu->pred_depth[0]));
			memset(&ctu->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu->part_size_type[0]));
			if(et->rd_mode==RD_FULL)
			{
				memset(&ctu_rd->pred_depth[abs_index], curr_depth-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu_rd->pred_depth[0]));
				memset(&ctu_rd->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu_rd->part_size_type[0]));
			}
		}

	}
	else
	{
		//top-down computation results are prefered
		PartSize part_size_type2 = (parent_part_info->depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//
		int parent_depth = parent_part_info->depth;//curr_depth-1

		synchronize_motion_buffers_luma(et, parent_part_info, et->transform_quant_wnd[parent_depth+1], et->transform_quant_wnd[0], et->decoded_mbs_wnd[parent_depth+1], et->decoded_mbs_wnd[0], gcnt);
		synchronize_motion_buffers_chroma(et, parent_part_info, et->transform_quant_wnd[parent_depth+1], et->transform_quant_wnd[0], et->decoded_mbs_wnd[parent_depth+1], et->decoded_mbs_wnd[0], gcnt);

		CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, parent_depth, abs_index, num_part_in_cu);
		SET_INTER_INFO_BUFFS(et, ctu, parent_part_info, abs_index, num_part_in_cu);//, REF_PIC_LIST_0);

		memset(&ctu->qp[abs_index], parent_part_info->qp, parent_part_info->num_part_in_cu*sizeof(ctu->qp[0]));
		//if we fill this in here we don't have to consolidate
		memset(&ctu->pred_depth[abs_index], parent_part_info->depth-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu->pred_depth[0]));
		memset(&ctu->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu->part_size_type[0]));
		if(et->rd_mode==RD_FULL)//rd
		{
			memset(&ctu_rd->pred_depth[abs_index], parent_depth -(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu_rd->pred_depth[0]));
			memset(&ctu_rd->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu_rd->part_size_type[0]));
		}
	}
}

void get_back_consolidated_info(henc_thread_t *et, ctu_info_t *ctu, cu_partition_info_t *curr_cu_info, int curr_depth)
{
	int abs_index = curr_cu_info->abs_index;
	int num_part_in_cu = curr_cu_info->num_part_in_cu;

	memcpy(&et->cbf_buffs[Y_COMP][curr_depth][abs_index], &ctu->cbf[Y_COMP][abs_index], num_part_in_cu*sizeof(et->cbf_buffs[0][0][0]));
	memcpy(&et->cbf_buffs[U_COMP][curr_depth][abs_index], &ctu->cbf[U_COMP][abs_index], num_part_in_cu*sizeof(et->cbf_buffs[0][0][0]));
	memcpy(&et->cbf_buffs[V_COMP][curr_depth][abs_index], &ctu->cbf[V_COMP][abs_index], num_part_in_cu*sizeof(et->cbf_buffs[0][0][0]));
	memcpy(&et->tr_idx_buffs[curr_depth][abs_index], &ctu->tr_idx[abs_index], num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));
	memcpy(&et->intra_mode_buffs[Y_COMP][curr_depth][abs_index], &ctu->intra_mode[Y_COMP][abs_index], num_part_in_cu*sizeof(et->intra_mode_buffs[0][0][0]));
	memcpy(&et->intra_mode_buffs[CHR_COMP][curr_depth][abs_index], &ctu->intra_mode[CHR_COMP][abs_index], num_part_in_cu*sizeof(et->intra_mode_buffs[0][0][0]));

	synchronize_motion_buffers_luma(et, curr_cu_info, et->transform_quant_wnd[0], et->transform_quant_wnd[curr_depth+1], et->decoded_mbs_wnd[0], et->decoded_mbs_wnd[curr_depth+1], 0);
	synchronize_motion_buffers_chroma(et, curr_cu_info, et->transform_quant_wnd[0], et->transform_quant_wnd[curr_depth+1], et->decoded_mbs_wnd[0], et->decoded_mbs_wnd[curr_depth+1], 0);
}

void put_consolidated_info(henc_thread_t *et, ctu_info_t *ctu, cu_partition_info_t *curr_cu_info, int curr_depth)
{
//	int curr_depth = curr_cu_info->depth;
	int abs_index = curr_cu_info->abs_index;
	int num_part_in_cu = curr_cu_info->num_part_in_cu;

	memcpy(&ctu->cbf[Y_COMP][abs_index], &et->cbf_buffs[Y_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(et->cbf_buffs[0][0][0]));
	memcpy(&ctu->cbf[U_COMP][abs_index], &et->cbf_buffs[U_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(et->cbf_buffs[0][0][0]));
	memcpy(&ctu->cbf[V_COMP][abs_index], &et->cbf_buffs[V_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(et->cbf_buffs[0][0][0]));
	memcpy(&ctu->tr_idx[abs_index], &et->tr_idx_buffs[curr_depth][abs_index], num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));
	memcpy(&ctu->intra_mode[Y_COMP][abs_index], &et->intra_mode_buffs[Y_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(et->intra_mode_buffs[0][0][0]));
	memcpy(&ctu->intra_mode[CHR_COMP][abs_index], &et->intra_mode_buffs[CHR_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(et->intra_mode_buffs[0][0][0]));

	synchronize_motion_buffers_luma(et, curr_cu_info, et->transform_quant_wnd[curr_depth+1], et->transform_quant_wnd[0], et->decoded_mbs_wnd[curr_depth+1], et->decoded_mbs_wnd[0], 0);
	synchronize_motion_buffers_chroma(et, curr_cu_info, et->transform_quant_wnd[curr_depth+1], et->transform_quant_wnd[0], et->decoded_mbs_wnd[curr_depth+1], et->decoded_mbs_wnd[0], 0);
}


#define EQUAL_MOTION(mv_cand0, poc_cand0, mv_cand1, poc_cand1) (poc_cand0 == poc_cand1 && mv_cand0->ref_idx == mv_cand1->ref_idx && mv_cand0->mv.hor_vector==mv_cand1->mv.hor_vector && mv_cand0->mv.ver_vector==mv_cand1->mv.ver_vector)

//xCheckRDCostMerge2Nx2N
uint32_t check_rd_cost_merge_2nx2n(henc_thread_t* et, ctu_info_t* ctu, int depth, int position)
{
	int gcnt = 0;
	picture_t *currpict = &et->enc_engine->current_pict;
	slice_t *currslice = &currpict->slice;
	int ref_idx = 0;
	int uiNoResidual, uiMergeCand;
	int mergeCandBuffer[MERGE_MVP_MAX_NUM_CANDS] = {0,0,0,0,0};
	int numValidMergeCand = et->enc_engine->num_merge_mvp_candidates;//MERGE_MVP_MAX_NUM_CANDS;
	int bestIsSkip = FALSE;
	int mv_cost;
	cu_partition_info_t	*curr_cu_info = &ctu->partition_list[et->partition_depth_start[depth]]+position;
	int abs_index = curr_cu_info->abs_index;
	PartSize part_size_type = SIZE_2Nx2N;
	int pred_buff_stride, orig_buff_stride, reference_buff_stride, residual_buff_stride;
	int pred_buff_stride_chroma, orig_buff_stride_chroma, reference_buff_stride_chroma, residual_buff_stride_chroma;
	int16_t  *orig_buff, *orig_buff_u, *orig_buff_v;
	int16_t  *reference_buff_cu_position, *reference_buff_cu_position_u, *reference_buff_cu_position_v;
	int16_t *pred_buff, *pred_buff_u, *pred_buff_v, *residual_buff, *residual_buff_u, *residual_buff_v;
	wnd_t *reference_wnd=NULL;//, *resi_wnd = NULL;
//	uint8_t *cbf_buff = NULL;
	int curr_part_size, curr_part_size_shift;
	int curr_part_size_chroma, curr_part_size_shift_chroma;
	int curr_part_x, curr_part_y, curr_part_global_x, curr_part_global_y;
	int curr_part_x_chroma, curr_part_y_chroma, curr_part_global_x_chroma, curr_part_global_y_chroma;
	int curr_depth = depth;
	uint32_t dist, best_dist = MAX_COST;
	uint32_t cost, best_cost = MAX_COST;
	uint32_t best_sum = 0;
	mv_candidate_t best_mv0, best_mv1, mv;
	int best_candidate = 0, is_skipped;
	int chr_qp_offset = et->enc_engine->chroma_qp_offset;
	double weight = pow( 2.0, (currslice->qp-chroma_scale_conversion_table[clip(currslice->qp+chr_qp_offset,0,57)])/3.0 );
	int motion_compensation_done = FALSE;
	byte inter_mode_neighbours[MERGE_MVP_MAX_NUM_CANDS] = {-1,-1,-1,-1,-1};

	get_merge_mvp_candidates(et, currslice, ctu, curr_cu_info, part_size_type, inter_mode_neighbours);//get candidates for merge motion motion vector prediction from the neigbour CUs	

	curr_depth = curr_cu_info->depth;
	curr_part_x = curr_cu_info->x_position;
	curr_part_y = curr_cu_info->y_position;
	curr_part_global_x = ctu->x[Y_COMP]+curr_part_x;
	curr_part_global_y = ctu->y[Y_COMP]+curr_part_y;
	curr_part_size = curr_cu_info->size;
	curr_part_size_shift = et->max_cu_size_shift-curr_depth;
	curr_part_x_chroma = curr_cu_info->x_position_chroma;
	curr_part_y_chroma = curr_cu_info->y_position_chroma;
	curr_part_global_x_chroma = ctu->x[CHR_COMP]+curr_part_x_chroma;
	curr_part_global_y_chroma = ctu->y[CHR_COMP]+curr_part_y_chroma;
	curr_part_size_chroma = curr_cu_info->size_chroma;
	curr_part_size_shift_chroma = et->max_cu_size_shift-curr_depth-1;//420
//	curr_adi_size = 2*2*curr_part_size+1;

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd[0], Y_COMP);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	pred_buff_stride_chroma = WND_STRIDE_2D(et->prediction_wnd[0], CHR_COMP);
	pred_buff_u = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	pred_buff_v = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
	orig_buff = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride_chroma = WND_STRIDE_2D(et->curr_mbs_wnd, CHR_COMP);
	orig_buff_u = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	orig_buff_v = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

	residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, Y_COMP);
	residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_buff_stride_chroma = WND_STRIDE_2D(et->residual_wnd, CHR_COMP);
	residual_buff_u = WND_POSITION_2D(int16_t *, et->residual_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	residual_buff_v = WND_POSITION_2D(int16_t *, et->residual_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

	for(uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand)
	{
		int data_size_x, data_size_y, data_padding_x, data_padding_y;
		int search_point_x, search_point_y;
		int xlow, xhigh, ylow, yhigh;//limits
		mv_candidate_t *mv_cand0, *mv_cand1;

		motion_compensation_done = FALSE;

		for(uiNoResidual = 0; uiNoResidual < 2; ++uiNoResidual )
		{
			if(!(uiNoResidual==1 && mergeCandBuffer[uiMergeCand]==1))
			{
				if( !(bestIsSkip && uiNoResidual == 0) )
				{
					int16_t *pred_aux_buffs[2][NUM_PICT_COMPONENTS];

					if(!motion_compensation_done)
					{
						uint64_t ref_poc0, ref_poc1;
						int ref_list, ref_list_init, ref_list_max;
						int is_bi_predict = 0;	
						mv_cand0 = &et->merge_mvp_candidates[REF_PIC_LIST_0].mv_candidates[uiMergeCand];
						mv_cand1 = &et->merge_mvp_candidates[REF_PIC_LIST_1].mv_candidates[uiMergeCand];

						ref_poc0 = currslice->ref_poc_list[REF_PIC_LIST_0][mv_cand0->ref_idx];
						ref_poc1 = currslice->ref_poc_list[REF_PIC_LIST_1][mv_cand1->ref_idx];

						if(currslice->slice_type == P_SLICE || (currslice->slice_type == B_SLICE && ((mv_cand0->ref_idx>=0 && mv_cand1->ref_idx<0) || EQUAL_MOTION(mv_cand0, ref_poc0, mv_cand1, ref_poc1))))
						{
							//just 1 reference using list 0
							ref_list_max = ref_list_init = REF_PIC_LIST_0;
							pred_aux_buffs[REF_PIC_LIST_0][Y_COMP] = pred_buff;
							pred_aux_buffs[REF_PIC_LIST_0][U_COMP] = pred_buff_u;
							pred_aux_buffs[REF_PIC_LIST_0][V_COMP] = pred_buff_v;
						}
						else if(currslice->slice_type == B_SLICE && (mv_cand0->ref_idx<0 && mv_cand1->ref_idx>=0))
						{
							//just 1 reference using list 1
							ref_list_max = ref_list_init = REF_PIC_LIST_1;
							pred_aux_buffs[REF_PIC_LIST_1][Y_COMP] = pred_buff;
							pred_aux_buffs[REF_PIC_LIST_1][U_COMP] = pred_buff_u;
							pred_aux_buffs[REF_PIC_LIST_1][V_COMP] = pred_buff_v;						
						}
						else if(currslice->slice_type == B_SLICE && (mv_cand0->ref_idx>=0 && mv_cand1->ref_idx>=0))
						{
							//Two references using weighted prediction
							is_bi_predict = 1;
							ref_list_init = REF_PIC_LIST_0;
							ref_list_max =  REF_PIC_LIST_1;

							pred_aux_buffs[REF_PIC_LIST_0][Y_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[1], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
							pred_aux_buffs[REF_PIC_LIST_0][U_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[1], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
							pred_aux_buffs[REF_PIC_LIST_0][V_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[1], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

							pred_aux_buffs[REF_PIC_LIST_1][Y_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[2], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
							pred_aux_buffs[REF_PIC_LIST_1][U_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[2], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
							pred_aux_buffs[REF_PIC_LIST_1][V_COMP] = WND_POSITION_2D(int16_t *, et->prediction_wnd[2], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
						}
						else
						{
							int iiiii=0;
						}

						for (ref_list = ref_list_init; ref_list <= ref_list_max; ref_list++)
						{
							mv = et->merge_mvp_candidates[ref_list].mv_candidates[uiMergeCand];
							reference_wnd = &currslice->ref_pic_list[ref_list][mv.ref_idx]->img;//[0] up to now we only use one reference, for more than one reference maybe this should be moved down
							reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
							reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
							reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, CHR_COMP);
							reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
							reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
#ifndef COMPUTE_AS_HM
							data_size_x = WND_WIDTH_2D(*reference_wnd, Y_COMP);
							data_size_y = WND_HEIGHT_2D(*reference_wnd, Y_COMP);
							data_padding_x = reference_wnd->data_padding_x[Y_COMP];
							data_padding_y = reference_wnd->data_padding_y[Y_COMP];
	
							xlow= 0 - data_padding_x;
							xhigh= data_size_x + data_padding_x;
							ylow= 0 - data_padding_y;
							yhigh= data_size_y + data_padding_y;

							search_point_x = curr_part_global_x + mv.mv.hor_vector/4;
							search_point_y = curr_part_global_y + mv.mv.ver_vector/4;

							if (search_point_x<xlow || search_point_x+curr_part_size>xhigh || search_point_y<ylow || search_point_y+curr_part_size>yhigh)
								continue;
#endif
							//create prediction if needed
							hmr_motion_compensation_luma(et, curr_cu_info, reference_buff_cu_position, reference_buff_stride, pred_aux_buffs[ref_list][Y_COMP], pred_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, &mv.mv, is_bi_predict);
							hmr_motion_compensation_chroma(et, reference_buff_cu_position_u, reference_buff_stride_chroma, pred_aux_buffs[ref_list][U_COMP], pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv.mv, is_bi_predict);
							hmr_motion_compensation_chroma(et, reference_buff_cu_position_v, reference_buff_stride_chroma, pred_aux_buffs[ref_list][V_COMP], pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv.mv, is_bi_predict);				
						}

						if(is_bi_predict)//currslice->slice_type == B_SLICE && (mv_cand0->ref_idx>=0 && mv_cand1->ref_idx>=0) && !EQUAL_MOTION(mv_cand0, mv_cand1))
						{
							et->funcs->weighted_average_motion(pred_aux_buffs[REF_PIC_LIST_0][Y_COMP], pred_buff_stride, pred_aux_buffs[REF_PIC_LIST_1][Y_COMP], pred_buff_stride, pred_buff, pred_buff_stride, curr_cu_info->size, curr_cu_info->size, et->bit_depth);
							et->funcs->weighted_average_motion(pred_aux_buffs[REF_PIC_LIST_0][U_COMP], pred_buff_stride_chroma, pred_aux_buffs[REF_PIC_LIST_1][U_COMP], pred_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, curr_cu_info->size_chroma, curr_cu_info->size_chroma, et->bit_depth);
							et->funcs->weighted_average_motion(pred_aux_buffs[REF_PIC_LIST_0][V_COMP], pred_buff_stride_chroma, pred_aux_buffs[REF_PIC_LIST_1][V_COMP], pred_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, curr_cu_info->size_chroma, curr_cu_info->size_chroma, et->bit_depth);
						}
						motion_compensation_done = TRUE;
					}

					// set MC parameters
					if(uiNoResidual==0)
					{
						et->funcs->predict(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride, residual_buff, residual_buff_stride, curr_part_size);
						et->funcs->predict(orig_buff_u, orig_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, residual_buff_u, residual_buff_stride_chroma, curr_part_size_chroma);
						et->funcs->predict(orig_buff_v, orig_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, residual_buff_v, residual_buff_stride_chroma, curr_part_size_chroma);
						cost = dist = encode_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
#ifndef COMPUTE_AS_HM
						cost += cost_rd(et->enc_engine->avg_dist, curr_cu_info->sum);
#endif
					}
					else
					{
						dist = (uint32_t) et->funcs->ssd16b(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride, curr_part_size);
						dist += (uint32_t) (weight*et->funcs->ssd16b(orig_buff_u, orig_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, curr_part_size_chroma));
						dist += (uint32_t) (weight*et->funcs->ssd16b(orig_buff_v, orig_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, curr_part_size_chroma));

//						dist = MAX_COST;
						curr_cu_info->inter_cbf[Y_COMP] = curr_cu_info->inter_cbf[U_COMP] = curr_cu_info->inter_cbf[V_COMP] = 0;
						curr_cu_info->inter_tr_idx = 0;		
						curr_cu_info->sum = 0;		
						cost = dist;
					}

					if(cost<best_cost)
					{
						best_mv0 = *mv_cand0;
						best_mv1 = *mv_cand1;
						best_candidate = uiMergeCand;
						best_dist = dist;
						best_cost = cost;
						best_sum = curr_cu_info->sum;
						if(uiNoResidual==1)//if it is skipped write 0 in the residual and copy the prediction wnd to the decoder wnd
						{
							wnd_copy_cu_2D(et->funcs->sse_copy_16_16, curr_cu_info, &et->prediction_wnd[0], et->decoded_mbs_wnd[curr_depth+1], ALL_COMP);
							wnd_zero_cu_1D(et, curr_cu_info, et->transform_quant_wnd[curr_depth+1]);
							SET_ENC_INFO_BUFFS(et, curr_cu_info, curr_depth/*+(part_size_type!=SIZE_2Nx2N)*/, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);//consolidate in prediction depth
						}
						put_consolidated_info(et, ctu, curr_cu_info, curr_depth);
						bestIsSkip = CBF_ALL(ctu, abs_index, 0) == 0;
						//exchange we keep in currdepth and the aux_depth, where we keep the best buffers
						//we keep best windows in the last depth as we are only processing 2nx2n
//						ptrswap(wnd_t*, et->transform_quant_wnd[curr_depth+1],et->transform_quant_wnd[NUM_QUANT_WNDS-1]);
//						ptrswap(wnd_t*, et->decoded_mbs_wnd[curr_depth+1],et->decoded_mbs_wnd[NUM_DECODED_WNDS-1]);
					}

					if ( uiNoResidual == 0 && CBF_ALL(ctu, abs_index, curr_depth) == 0)
					{
						// If no residual when allowing for one, then set mark to not try case where residual is forced to 0
						mergeCandBuffer[uiMergeCand] = 1;
					}
				}
			}
		}
	}


	//put the buffers with the best option in the correct depth
//	ptrswap(wnd_t*, et->transform_quant_wnd[NUM_QUANT_WNDS-1], et->transform_quant_wnd[curr_depth+1]);
//	ptrswap(wnd_t*, et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], et->decoded_mbs_wnd[curr_depth+1]);

	//save the information
	curr_cu_info->skipped = bestIsSkip;
	curr_cu_info->inter_mv[REF_PIC_LIST_0] = best_mv0.mv;
	curr_cu_info->inter_ref_index[REF_PIC_LIST_0] = best_mv0.ref_idx;
	curr_cu_info->inter_mv[REF_PIC_LIST_1] = best_mv1.mv;
	curr_cu_info->inter_ref_index[REF_PIC_LIST_1] = best_mv1.ref_idx;
	curr_cu_info->cost = curr_cu_info->distortion = best_dist;
	curr_cu_info->merge_flag = TRUE;
	curr_cu_info->merge_idx = best_candidate;
	curr_cu_info->inter_mode = inter_mode_neighbours[best_candidate];
	curr_cu_info->sum = best_sum;
	return best_dist;
}


extern int aux_dbg;
uint32_t motion_inter_full(henc_thread_t* et, ctu_info_t* ctu)
{
	int gcnt = 0;
	picture_t *currpict = &et->enc_engine->current_pict;
	slice_t *currslice = &currpict->slice;
	double dist, best_cost;//, cost_luma, cost_chroma;
	int position = 0;
	int curr_depth = 0;
	ctu_info_t *ctu_rd = et->ctu_rd;
	cu_partition_info_t	*parent_part_info = NULL;
	cu_partition_info_t	*curr_cu_info = ctu->partition_list;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	uint cost_sum[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int abs_index;
	int num_part_in_cu;
	double consumed_distortion = 0, avg_distortion = et->enc_engine->avg_dist;
	int perf_min_depth = et->enc_engine->performance_min_depth;
	int perf_fast_skip = et->enc_engine->performance_fast_skip_loop;

#ifndef COMPUTE_AS_HM
	int ithreads;
	uint total_intra_partitions = 0, total_partitions = 0;

	for(ithreads=0;ithreads<et->wfpp_num_threads;ithreads++)
	{
		henc_thread_t* henc_th = et->enc_engine->thread[ithreads];
		
		consumed_distortion += henc_th->acc_dist;
		total_intra_partitions += henc_th->num_intra_partitions;
		total_partitions += henc_th->num_total_partitions;
	}
//	total_partitions = 0;

	if(total_partitions==0)
		total_partitions = 1;

	if(total_partitions>10*et->num_partitions_in_cu || total_partitions>et->enc_engine->pict_total_ctu*et->num_partitions_in_cu/15)
	{
		avg_distortion = consumed_distortion/(total_partitions);
		if(et->enc_engine->is_scene_change)
			et->enc_engine->avg_dist = avg_distortion;//update avg_dist as it evolves
	}
//	else
		avg_distortion = et->enc_engine->avg_dist;

	if(et->index==0 && currslice->slice_type == P_SLICE && et->enc_engine->num_encoded_frames >1 && et->enc_engine->is_scene_change == 0 && 20<currslice->poc-et->enc_engine->hvenc->last_gop_reinit && total_partitions>et->enc_engine->pict_total_ctu*et->num_partitions_in_cu/10)
	{
		if(total_intra_partitions > (total_partitions*.7))
		{
			hmr_rc_change_pic_mode(et, currslice);
			et->enc_engine->is_scene_change = 1;
			if(et->enc_engine->gop_reinit_on_scene_change)
				et->enc_engine->hvenc->last_intra = et->enc_engine->last_intra = currslice->poc;
			et->enc_engine->last_gop_reinit = currslice->poc;
			et->enc_engine->hvenc->last_gop_reinit = currslice->poc;
#ifdef DBG_TRACE
			printf("\r\n---------------------scene change detected. frame:%d, total_intra_partitions:%d, total_partitions:%d , enc_engine->avg_dist:%.2f, avg_distortion:%.2f, ----------------------\r\n", et->enc_engine->num_encoded_frames, total_intra_partitions, total_partitions, et->enc_engine->avg_dist, avg_distortion);
#endif

		}
	}

//	avg_distortion = et->enc_engine->avg_dist;
#endif

	//init rd auxiliar ctu
	if(et->rd_mode != RD_DIST_ONLY)
	{
		copy_ctu(ctu, ctu_rd);
	}

	while(curr_depth!=0 || depth_state[curr_depth]!=1)
	{
		double cost = 0, intra_cost = 0;
		int stop_recursion = FALSE;
		int is_skipped = FALSE;
		PartSize part_size_type = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;
		int ref_list;
		int ref_idx;
		curr_depth = curr_cu_info->depth;
		num_part_in_cu = curr_cu_info->num_part_in_cu;
		abs_index = curr_cu_info->abs_index;

		position = curr_cu_info->list_index - et->partition_depth_start[curr_depth];

		//rc
		if(currslice->slice_type != I_SLICE && curr_depth<=et->enc_engine->qp_depth)
		{
//			curr_cu_info->variance = calc_variance_cu(et, curr_cu_info);
		}

		curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);

		//  Uni-directional prediction
		if(curr_cu_info->is_b_inside_frame && curr_cu_info->is_r_inside_frame)//if br (and tl) are inside the frame, process
		{
			int mv_cost = 0;

			if(part_size_type == SIZE_2Nx2N)
			{
				uint sad, merge_dist = MAX_COST, merge_cost = MAX_COST, merge_sum;
				motion_vector_t merge_mv0, merge_mv1;
				int merge_ref_idx0, merge_ref_idx1;
				uint merge_inter_mode;
				uint motion_estimation_precision = (et->enc_engine->motion_estimation_precision*2-1);//compute all precisions below the configured

				//encode inter
				curr_cu_info->prediction_mode = INTER_MODE;
				if(curr_cu_info->depth >= perf_min_depth)//-1)
				{
					if(et->enc_engine->num_encoded_frames==2 && ctu->ctu_number==5 && curr_cu_info->abs_index>=144 && curr_cu_info->depth==2)// && curr_cu_info->abs_index>=32)// && curr_partition_info->abs_index>=192)// && curr_depth==2)
					{
						int iiiiii=0;
					}

					merge_dist = check_rd_cost_merge_2nx2n(et, ctu, curr_depth, position);//after this call we just need to keep the merge_idx, the motion vector is no longer needed for merge types
					merge_mv0 = curr_cu_info->inter_mv[REF_PIC_LIST_0];
					merge_ref_idx0 = curr_cu_info->inter_ref_index[REF_PIC_LIST_0];
					merge_mv1 = curr_cu_info->inter_mv[REF_PIC_LIST_1];
					merge_ref_idx1 = curr_cu_info->inter_ref_index[REF_PIC_LIST_1];
					merge_inter_mode = curr_cu_info->inter_mode;
					merge_sum = curr_cu_info->sum;
					merge_cost = merge_dist;

#ifndef COMPUTE_AS_HM
					merge_cost +=curr_cu_info->merge_idx;
					merge_cost=calc_cost_full(merge_cost, curr_depth, avg_distortion);
					merge_cost+=cost_rd(et->enc_engine->avg_dist, curr_cu_info->sum);
					cost = merge_cost;
					dist = merge_dist;
	//				merge_cost += clip(et->enc_engine->avg_dist/1.75,5.,20000.)*curr_cu_info->sum;

					if(curr_cu_info->skipped)
					{
						merge_cost*=.95;
					}
					if((curr_cu_info->skipped && merge_cost < avg_distortion*ctu->num_part_in_ctu/2.5))
						is_skipped = TRUE;
#endif
//					if(curr_cu_info->depth == perf_min_depth-1)
//						get_back_consolidated_info(et, ctu, curr_cu_info, curr_depth);
				}
				else
				{
					memset(curr_cu_info->inter_mv,0,sizeof(curr_cu_info->inter_mv));
					curr_cu_info->merge_flag = FALSE;
					curr_cu_info->skipped = FALSE;
					cost = curr_cu_info->cost = merge_cost = (uint32_t)MAX_COST;
					dist = curr_cu_info->distortion = merge_dist = (uint32_t)MAX_COST;				
					is_skipped = FALSE;
				}
//				merge_cost = MAX_COST;
//				put_consolidated_info(et, ctu, curr_cu_info, curr_depth);

				if(curr_cu_info->depth >= perf_min_depth)
				{
					if(et->enc_engine->num_encoded_frames==2 && ctu->ctu_number==2 && curr_cu_info->depth==1)// && curr_cu_info->abs_index>=32)// && curr_partition_info->abs_index>=192)// && curr_depth==2)
					{
						int iiiiii=0;
					}


					if(!is_skipped)
						sad = hmr_cu_motion_estimation(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N, 2*curr_cu_info->size*curr_cu_info->size, motion_estimation_precision);//(MOTION_PEL_MASK|MOTION_HALF_PEL_MASK|MOTION_QUARTER_PEL_MASK));//.25*avg_distortion*curr_cu_info->num_part_in_cu);
#ifdef COMPUTE_AS_HM

					mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
/*					if(check_unidirectional_motion(currslice, curr_cu_info))
					{
						mv_cost = predict_inter_uni(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					}
					else
					{
						mv_cost = predict_inter_bi(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);//.25*avg_distortion*curr_cu_info->num_part_in_cu);					
					}
*/
					dist = encode_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
#else
					if(!is_skipped && (curr_cu_info->size < 64 || (sad<100*curr_cu_info->num_part_in_cu)))// && curr_cu_info->size == 64))
					{
						mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
						dist = encode_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
					}
					else
					{
	//					printf("64x64 inter computation skipped");
						mv_cost = 0;
						dist = MAX_COST;
						curr_cu_info->sum = 0;
						curr_cu_info->inter_mv[REF_PIC_LIST_0].hor_vector = curr_cu_info->inter_mv[REF_PIC_LIST_0].ver_vector = 0;
						curr_cu_info->inter_mv[REF_PIC_LIST_1].hor_vector = curr_cu_info->inter_mv[REF_PIC_LIST_1].ver_vector = 0;
					}
#endif
					cost = dist;

					cost += 2*mv_cost;
					//cost=calc_cost(cost, curr_depth);
					cost=calc_cost_full(cost, curr_depth, avg_distortion);
					cost+=cost_rd(et->enc_engine->avg_dist, curr_cu_info->sum);
	//				cost += clip(et->enc_engine->avg_dist/1.75,5.,20000.)*curr_cu_info->sum;
	//				cost=cost*DEPHT_SCALE+DEPHT_ADD*curr_depth;
#ifdef COMPUTE_AS_HM
					cost=dist+5*curr_depth;
#endif
	//				if(cost+clip(avg_distortion/1.75,5.,20000.)*(curr_cu_info->sum+2)<merge_cost+clip(avg_distortion/1.75,5.,20000.)*merge_sum)// || dist/curr_cu_info->num_part_in_cu < 10000)
					if(cost<merge_cost)// || dist/curr_cu_info->num_part_in_cu < 10000)
					{
						curr_cu_info->merge_flag = FALSE;
						curr_cu_info->skipped = FALSE;
//						curr_cu_info->inter_mode = (1<<REF_PIC_LIST_0);
	//					is_skipped = FALSE;
	//					put_consolidated_info(et, ctu, curr_cu_info, curr_depth);
					}
					else
					{
						curr_cu_info->inter_mv[REF_PIC_LIST_0] = merge_mv0;
						curr_cu_info->inter_ref_index[REF_PIC_LIST_0] = merge_ref_idx0;
						curr_cu_info->inter_mv[REF_PIC_LIST_1] = merge_mv1;
						curr_cu_info->inter_ref_index[REF_PIC_LIST_1] = merge_ref_idx1;
						curr_cu_info->inter_mode = merge_inter_mode;
						curr_cu_info->sum = merge_sum;
						cost = merge_cost;
						dist = merge_dist;
					}
				}

				curr_cu_info->cost = (uint32_t)cost;
				curr_cu_info->distortion = (uint32_t)dist;

#ifndef COMPUTE_AS_HM
				//fast skip
				if(perf_fast_skip && (dist == 0 || (curr_cu_info->sum < curr_cu_info->num_part_in_cu && dist<.25*avg_distortion*curr_cu_info->num_part_in_cu) || (curr_cu_info->sum == 0 && dist<avg_distortion*curr_cu_info->num_part_in_cu)))
				{
					int max_processing_depth;
					if(curr_cu_info->merge_flag)
						get_back_consolidated_info(et, ctu, curr_cu_info, curr_depth);
					
					stop_recursion = TRUE;
					consolidate_prediction_info(et, ctu, ctu_rd, curr_cu_info, curr_cu_info->cost, MAX_COST, FALSE, NULL);
					max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

					if(curr_depth <= max_processing_depth)//
					{
						int aux_depth;
						cu_partition_info_t* aux_partition_info = curr_cu_info;//(parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth])]:&ctu->partition_list[0];
						abs_index = aux_partition_info->abs_index;
						num_part_in_cu  = aux_partition_info->num_part_in_cu;

						for(aux_depth=curr_depth;aux_depth<=max_processing_depth;aux_depth++)
						{
							synchronize_reference_buffs(et, aux_partition_info, et->decoded_mbs_wnd[0], et->decoded_mbs_wnd[aux_depth+1], gcnt);	
							//for rd
							if(et->rd_mode!=RD_DIST_ONLY)
								consolidate_info_buffers_for_rd(et, ctu, aux_depth, abs_index, num_part_in_cu);
						}
						synchronize_reference_buffs_chroma(et, aux_partition_info, et->decoded_mbs_wnd[0], et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
					}
				}

				if(perf_fast_skip && ((curr_cu_info->depth >= perf_min_depth) && !stop_recursion && !is_skipped && (curr_cu_info->size<32 || sad>400*curr_cu_info->num_part_in_cu)))//if(0)
#endif
				{
					//encode intra
					uint inter_sum = curr_cu_info->sum;
					uint intra_dist;
					if(!curr_cu_info->merge_flag)
						put_consolidated_info(et, ctu, curr_cu_info, curr_depth);
					intra_dist = encode_intra(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
#ifdef COMPUTE_AS_HM
					intra_cost = intra_dist+5*curr_depth;
					if(intra_cost < cost)
#else
					intra_cost = intra_dist*(1.275-clip(((double)total_intra_partitions/(double)total_partitions), .0, .15))+clip(avg_distortion-400,40,avg_distortion)/1.75*curr_depth;
					intra_cost+=cost_rd(et->enc_engine->avg_dist, curr_cu_info->sum);
					//if(intra_cost*(1.+curr_cu_info[0].qp/50.0)<cost)
					if(intra_cost<cost)
#endif
					{	//we prefer intra and it is already in its buffer
						curr_cu_info->cost = (uint32_t) intra_cost;
						curr_cu_info->distortion = (uint32_t) intra_dist;
						//curr_cu_info->sum = curr_cu_info->sum;
						curr_cu_info->prediction_mode = INTRA_MODE;
						curr_cu_info->merge_flag = FALSE;
						curr_cu_info->skipped = FALSE;
					}
					else
					{	//we prefer inter, bring it back
						get_back_consolidated_info(et, ctu, curr_cu_info, curr_depth);
						curr_cu_info->cost = (uint32_t) cost;
						curr_cu_info->distortion = (uint32_t) dist;
						curr_cu_info->sum = inter_sum;
						curr_cu_info->prediction_mode = INTER_MODE;
					}
				}
#ifndef COMPUTE_AS_HM
				else if((curr_cu_info->depth >= perf_min_depth) && !stop_recursion)
				{
					if(curr_cu_info->merge_flag)
						get_back_consolidated_info(et, ctu, curr_cu_info, curr_depth);
				}
#endif
#ifndef COMPUTE_AS_HM
			//if this matches, it is useless to continue the recursion. the case where curr_depth!=et->max_pred_partition_depth is checked at the end of the consolidation loop)
/*				if(curr_depth>0 && depth_state[curr_depth]!=3 && curr_depth == et->max_pred_partition_depth && cost_sum[curr_depth]+curr_cu_info->cost>parent_part_info->cost && parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
				{
					depth_state[curr_depth] = 3;
				}
*/
#endif
			}
			else if(part_size_type == SIZE_NxN)//intra NxN is processed in its current depth, while inter NxN is processed in its father�s depth. So, intra NxN does not have to be compaired
			{
				cost = dist = curr_cu_info->cost = curr_cu_info->distortion = MAX_COST;
				if((curr_depth-1) == (et->max_cu_depth - et->mincu_mintr_shift_diff) && curr_cu_info->parent->size>8)	//SIZE_NxN
				{
					uint sad;
					uint motion_estimation_precision = (et->enc_engine->motion_estimation_precision*2-1);//compute all precisions below the configured
					int position_aux = curr_cu_info->parent->list_index - et->partition_depth_start[curr_depth-1];
					uint aux_cost = curr_cu_info->parent->cost;
					uint aux_dist = curr_cu_info->parent->distortion;
					uint aux_sum = curr_cu_info->parent->sum;

					if(et->enc_engine->num_encoded_frames==1 && curr_cu_info->abs_index>=32)// && ctu->ctu_number==5 )// && curr_depth==3)
					{
						int iiiiii=0;
					}

//					mv_cost = predict_inter_uni(et, ctu, gcnt, curr_depth, position, SIZE_NxN, 0);//.25*avg_distortion*4*curr_cu_info->num_part_in_cu);
					sad = hmr_cu_motion_estimation(et, ctu, gcnt, curr_depth, position, SIZE_NxN, 0, motion_estimation_precision);//(MOTION_PEL_MASK|MOTION_HALF_PEL_MASK|MOTION_QUARTER_PEL_MASK));//.25*avg_distortion*curr_cu_info->num_part_in_cu);


					if(et->enc_engine->num_encoded_frames == 1 && ctu->ctu_number == 5 && curr_cu_info->abs_index == 176)
					{
						int iiiii=0;
					}

					mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_NxN);//.25*avg_distortion*curr_cu_info->num_part_in_cu);

/*					if(check_unidirectional_motion(currslice, curr_cu_info))
					{
						mv_cost = predict_inter_uni(et, ctu, gcnt, curr_depth, position, SIZE_NxN);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					}
					else
					{
						mv_cost = predict_inter_bi(et, ctu, gcnt, curr_depth, position, SIZE_NxN);//.25*avg_distortion*curr_cu_info->num_part_in_cu);					
					}
*/
//					mv_cost = predict_inter_uni(et, ctu, gcnt, curr_depth, position, SIZE_NxN);//.25*avg_distortion*curr_cu_info->num_part_in_cu);

					dist = encode_inter(et, ctu, gcnt, curr_depth-1, position_aux, SIZE_NxN);//this function is referenced by the initial depth, not by the processing depth
#ifdef COMPUTE_AS_HM
					cost=dist+5*curr_depth;
#else
//					cost = calc_cost(dist, curr_depth);
					cost=calc_cost_full(dist, curr_depth, avg_distortion);
//					cost=dist*DEPHT_SCALE+DEPHT_ADD*curr_depth;
					cost += mv_cost;
#endif
					curr_cu_info[0].cost = curr_cu_info[1].cost = curr_cu_info[2].cost = curr_cu_info[3].cost = 0;
					curr_cu_info[0].sum = curr_cu_info[1].sum = curr_cu_info[2].sum = curr_cu_info[3].sum = 0;
					curr_cu_info[0].distortion = curr_cu_info[1].distortion = curr_cu_info[2].distortion = curr_cu_info[3].distortion = 0;

					curr_cu_info->cost = (uint32_t)cost;
					curr_cu_info->distortion = (uint32_t)dist;
					curr_cu_info->sum = curr_cu_info->parent->sum;
					curr_cu_info[0].prediction_mode = curr_cu_info[1].prediction_mode = curr_cu_info[2].prediction_mode = curr_cu_info[3].prediction_mode = INTER_MODE;

					//restore parent values
					curr_cu_info->parent->cost = aux_cost;//este valor se sobreescribe en encode_inter. Deberia intentar hacer que el intra funcionase igual, con 1 sola llamada que hiciese las 4 particiones NxN. Asi se consolidaria el cbf totalmente
					curr_cu_info->parent->distortion = aux_dist;
					curr_cu_info->parent->sum = aux_sum;
#ifdef COMPUTE_AS_HM
				}
#else				
					if(curr_cu_info->size<64 && dist>(1.5*avg_distortion*(4*curr_cu_info->num_part_in_cu+4*curr_depth)/*+4000*curr_cu_info->num_part_in_cu*/))
//					if((dist>(1.5*avg_distortion*4*curr_cu_info->num_part_in_cu+4000*4*curr_cu_info->num_part_in_cu)) &&  curr_cu_info->cost>=curr_cu_info->parent->cost/5)
#endif	
					{
						//intra
						uint previous_sum = curr_cu_info->sum;
						uint intra_dist, intra_sum;

						put_consolidated_info(et, ctu, curr_cu_info->parent, curr_depth);
						intra_cost = 0;
						intra_dist = encode_intra(et, ctu, gcnt, curr_depth, position, part_size_type);
						intra_cost = calc_cost_full(intra_dist, curr_depth, avg_distortion);
						intra_sum = curr_cu_info[0].sum + curr_cu_info[1].sum + curr_cu_info[2].sum + curr_cu_info[3].sum;

						curr_cu_info[0].cost = curr_cu_info[1].cost = curr_cu_info[2].cost = curr_cu_info[3].cost = 0;
						curr_cu_info[0].sum = curr_cu_info[1].sum = curr_cu_info[2].sum = curr_cu_info[3].sum = 0;
						curr_cu_info[0].distortion = curr_cu_info[1].distortion = curr_cu_info[2].distortion = curr_cu_info[3].distortion = 0;
						curr_cu_info[0].skipped = curr_cu_info[1].skipped = curr_cu_info[2].skipped = curr_cu_info[3].skipped = 0;
#ifdef COMPUTE_AS_HM
						intra_cost=intra_dist+5*curr_depth;

/*						if(et->enc_engine->num_encoded_frames==10 && ctu->ctu_number == 72 && curr_cu_info->abs_index>=240)
						{
							intra_cost = cost+1;
						}
*/
						if(intra_cost < cost)
#else
						//if((1.25+curr_cu_info[0].qp/25.0)*intra_cost+45*intra_sum<cost+45*previous_sum)// && intra_cost<64*curr_cu_info->variance)
						if((1.25)*intra_cost+45*intra_sum<cost+45*previous_sum)// && intra_cost<64*curr_cu_info->variance)
#endif
						{	//we prefer intra and it is already in its buffer
							curr_cu_info->cost = (uint32_t)intra_cost;
							curr_cu_info->sum = intra_sum;
							curr_cu_info->distortion = (uint32_t)intra_dist;
							curr_cu_info[0].prediction_mode = curr_cu_info[1].prediction_mode = curr_cu_info[2].prediction_mode = curr_cu_info[3].prediction_mode = INTRA_MODE;
						}
						else
						{	//we prefer inter, bring it back
							get_back_consolidated_info(et, ctu, curr_cu_info->parent, curr_depth);
							curr_cu_info->cost = (uint32_t)cost;
							curr_cu_info->sum = previous_sum;
							curr_cu_info->distortion = (uint32_t)dist;
							curr_cu_info[0].prediction_mode = curr_cu_info[1].prediction_mode = curr_cu_info[2].prediction_mode = curr_cu_info[3].prediction_mode = INTER_MODE;
						}
					}
#ifndef COMPUTE_AS_HM
				}
#endif
				depth_state[curr_depth] = 3;
			}
		}
		else
		{
			curr_cu_info->cost = MAX_COST;
		}

		cost_sum[curr_depth] += curr_cu_info->cost;

		depth_state[curr_depth]++;

		if((curr_depth)<et->max_pred_partition_depth && curr_cu_info->is_tl_inside_frame && !stop_recursion)//depth_state[curr_depth]!=4 is for fast skip//if tl is not inside the frame don't process the next depths
		{
			curr_depth++;
			parent_part_info = curr_cu_info;
		}
		else if(depth_state[curr_depth]==4)//depth =1 already consolidated
		{
			int max_processing_depth;

			while(depth_state[curr_depth]==4 && curr_depth>0)//>0 pq consolidamos sobre el padre, 
			{
				int is_max_depth = (curr_depth==et->max_pred_partition_depth);
				cost = parent_part_info->children[0]->cost + parent_part_info->children[1]->cost +parent_part_info->children[2]->cost+parent_part_info->children[3]->cost;

				depth_state[curr_depth] = 0;

				best_cost = parent_part_info->cost;

				consolidate_prediction_info(et, ctu, ctu_rd, parent_part_info, (uint32_t)best_cost, (uint32_t)cost, is_max_depth, cost_sum);

				depth_state[curr_depth] = 0;
				cost_sum[curr_depth] = 0;

				curr_depth--;
				parent_part_info = parent_part_info->parent;

#ifndef COMPUTE_AS_HM
				if(curr_depth>perf_min_depth/*-1*/ && cost_sum[curr_depth] > parent_part_info->cost && depth_state[curr_depth]<4 && ctu->partition_list[0].is_b_inside_frame && ctu->partition_list[0].is_r_inside_frame)
				{
					depth_state[curr_depth] = 4;
				}
#endif
			}

			max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

			if(curr_depth <= max_processing_depth)//= is needed when et->max_intra_tr_depth!=4
			{
				int aux_depth;
				cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
				abs_index = aux_partition_info->abs_index;
				num_part_in_cu  = aux_partition_info->num_part_in_cu;

				for(aux_depth=curr_depth;aux_depth<=max_processing_depth;aux_depth++)
				{
					synchronize_reference_buffs(et, aux_partition_info, et->decoded_mbs_wnd[0], et->decoded_mbs_wnd[aux_depth+1], gcnt);	
					//for rd
//					if(et->rd_mode!=RD_DIST_ONLY)
//						CONSOLIDATE_INTRA_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu);
				}
				synchronize_reference_buffs_chroma(et, aux_partition_info, et->decoded_mbs_wnd[0], et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
			}
		}

		if(parent_part_info!=NULL)
			curr_cu_info = parent_part_info->children[depth_state[curr_depth]];
	}
	
	curr_cu_info = &ctu->partition_list[0];
	abs_index = curr_cu_info->abs_index;
	curr_depth = curr_cu_info->depth;
	num_part_in_cu = curr_cu_info->num_part_in_cu;

	//if pred_depth==0 there is no NxN subdivision. we need to collect the information of the ctu
	if(et->max_pred_partition_depth==0)
	{
		CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu);
		SET_INTER_INFO_BUFFS(et, ctu, curr_cu_info, abs_index, num_part_in_cu);//, REF_PIC_LIST_0);
	}


//	memset(&ctu->pred_mode[abs_index], INTER_MODE, num_part_in_cu*sizeof(ctu->pred_mode[0]));//signal all partitions as inter
//	memset(&ctu->skipped[abs_index], FALSE, num_part_in_cu*sizeof(ctu->skipped[0]));//signal all partitions as non skipped
	return curr_cu_info->cost;
}



//xCheckRDCostMerge2Nx2N
uint32_t check_rd_cost_merge_2nx2n_fast(henc_thread_t* et, ctu_info_t* ctu, int depth, int position)
{
	int gcnt = 0;
	picture_t *currpict = &et->enc_engine->current_pict;
	slice_t *currslice = &currpict->slice;
	int ref_idx = 0;
	int uiNoResidual, uiMergeCand;
	int mergeCandBuffer[MERGE_MVP_MAX_NUM_CANDS] = {0,0,0,0,0};
	int numValidMergeCand = et->enc_engine->num_merge_mvp_candidates;//MERGE_MVP_MAX_NUM_CANDS;
	int bestIsSkip = FALSE;
	int mv_cost;
	cu_partition_info_t	*curr_cu_info = &ctu->partition_list[et->partition_depth_start[depth]]+position;
	int abs_index = curr_cu_info->abs_index;
	PartSize part_size_type = SIZE_2Nx2N;
	int pred_buff_stride, orig_buff_stride, reference_buff_stride, residual_buff_stride;
	int pred_buff_stride_chroma, orig_buff_stride_chroma, reference_buff_stride_chroma, residual_buff_stride_chroma;
	int16_t  *orig_buff, *orig_buff_u, *orig_buff_v;
	int16_t  *reference_buff_cu_position, *reference_buff_cu_position_u, *reference_buff_cu_position_v;
	int16_t *pred_buff, *pred_buff_u, *pred_buff_v, *residual_buff, *residual_buff_u, *residual_buff_v;
	wnd_t *reference_wnd=NULL;//, *resi_wnd = NULL;
//	uint8_t *cbf_buff = NULL;
	int curr_part_size, curr_part_size_shift;
	int curr_part_size_chroma, curr_part_size_shift_chroma;
	int curr_part_x, curr_part_y, curr_part_global_x, curr_part_global_y;
	int curr_part_x_chroma, curr_part_y_chroma, curr_part_global_x_chroma, curr_part_global_y_chroma;
	int curr_depth = depth;
	uint32_t dist, best_dist = MAX_COST;
	uint32_t cost, best_cost = MAX_COST;
	uint32_t best_sum = 0;
	motion_vector_t mv, best_mv;
	int best_candidate = 0, is_skipped;
	int chr_qp_offset = et->enc_engine->chroma_qp_offset;
	double weight = pow( 2.0, (currslice->qp-chroma_scale_conversion_table[clip(currslice->qp+chr_qp_offset,0,57)])/3.0 );
	int motion_compensation_done = FALSE;
	uint8_t inter_mode_neighbours[MERGE_MVP_MAX_NUM_CANDS] = {-1,-1,-1,-1,-1};

	get_merge_mvp_candidates(et, currslice, ctu, curr_cu_info, part_size_type, inter_mode_neighbours);//get candidates for merge motion motion vector prediction from the neigbour CUs	

	curr_depth = curr_cu_info->depth;
	curr_part_x = curr_cu_info->x_position;
	curr_part_y = curr_cu_info->y_position;
	curr_part_global_x = ctu->x[Y_COMP]+curr_part_x;
	curr_part_global_y = ctu->y[Y_COMP]+curr_part_y;
	curr_part_size = curr_cu_info->size;
	curr_part_size_shift = et->max_cu_size_shift-curr_depth;
	curr_part_x_chroma = curr_cu_info->x_position_chroma;
	curr_part_y_chroma = curr_cu_info->y_position_chroma;
	curr_part_global_x_chroma = ctu->x[CHR_COMP]+curr_part_x_chroma;
	curr_part_global_y_chroma = ctu->y[CHR_COMP]+curr_part_y_chroma;
	curr_part_size_chroma = curr_cu_info->size_chroma;
	curr_part_size_shift_chroma = et->max_cu_size_shift-curr_depth-1;//420
//	curr_adi_size = 2*2*curr_part_size+1;

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd[0], Y_COMP);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	pred_buff_stride_chroma = WND_STRIDE_2D(et->prediction_wnd[0], CHR_COMP);
	pred_buff_u = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	pred_buff_v = WND_POSITION_2D(int16_t *, et->prediction_wnd[0], V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
	orig_buff = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride_chroma = WND_STRIDE_2D(et->curr_mbs_wnd, CHR_COMP);
	orig_buff_u = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	orig_buff_v = WND_POSITION_2D(int16_t *, et->curr_mbs_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

	residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, Y_COMP);
	residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_buff_stride_chroma = WND_STRIDE_2D(et->residual_wnd, CHR_COMP);
	residual_buff_u = WND_POSITION_2D(int16_t *, et->residual_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	residual_buff_v = WND_POSITION_2D(int16_t *, et->residual_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

	reference_wnd = &currslice->ref_pic_list[REF_PIC_LIST_0][ref_idx]->img;//[0] up to now we only use one reference	
	reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
	reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
	reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, CHR_COMP);
	reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
	reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);


	for(uiMergeCand=0; uiMergeCand<numValidMergeCand; ++uiMergeCand )
	{
		mv = et->merge_mvp_candidates[REF_PIC_LIST_0].mv_candidates[uiMergeCand].mv;
		motion_compensation_done = FALSE;
		//prune duplicated candidates
		if(uiMergeCand>0 && mv.hor_vector == et->merge_mvp_candidates[REF_PIC_LIST_0].mv_candidates[uiMergeCand-1].mv.hor_vector && mv.ver_vector == et->merge_mvp_candidates[REF_PIC_LIST_0].mv_candidates[uiMergeCand-1].mv.ver_vector)
			continue;

		for(uiNoResidual = 0; uiNoResidual < 2; ++uiNoResidual )
		{
			if(!(uiNoResidual==0 && mergeCandBuffer[uiMergeCand]==1))
			{
				if( !(bestIsSkip && uiNoResidual == 0) )
				{
					if(!motion_compensation_done)
					{
						//create prediction if needed
						hmr_motion_compensation_luma(et, curr_cu_info, reference_buff_cu_position, reference_buff_stride, pred_buff, pred_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, &mv, 0);
						hmr_motion_compensation_chroma(et, reference_buff_cu_position_u, reference_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv, 0);
						hmr_motion_compensation_chroma(et, reference_buff_cu_position_v, reference_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv, 0);	
					
						motion_compensation_done = TRUE;
					}

					// set MC parameters
					if(uiNoResidual==0)
					{
						et->funcs->predict(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride, residual_buff, residual_buff_stride, curr_part_size);
						et->funcs->predict(orig_buff_u, orig_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, residual_buff_u, residual_buff_stride_chroma, curr_part_size_chroma);
						et->funcs->predict(orig_buff_v, orig_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, residual_buff_v, residual_buff_stride_chroma, curr_part_size_chroma);
						cost = dist = encode_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
#ifndef COMPUTE_AS_HM
						cost += cost_rd(et->enc_engine->avg_dist, curr_cu_info->sum);
#endif
					}
					else
					{
						dist = (uint32_t) et->funcs->ssd16b(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride, curr_part_size);
						dist += (uint32_t) (weight*et->funcs->ssd16b(orig_buff_u, orig_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, curr_part_size_chroma));
						dist += (uint32_t) (weight*et->funcs->ssd16b(orig_buff_v, orig_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, curr_part_size_chroma));

//						dist = MAX_COST;
						curr_cu_info->inter_cbf[Y_COMP] = curr_cu_info->inter_cbf[U_COMP] = curr_cu_info->inter_cbf[V_COMP] = 0;
						curr_cu_info->inter_tr_idx = 0;		
						curr_cu_info->sum = 0;		
						cost = dist;
					}

					if(cost<best_cost)
					{
						best_mv = mv;
						best_candidate = uiMergeCand;
						best_dist = dist;
						best_cost = cost;
						best_sum = curr_cu_info->sum;
						if(uiNoResidual==1)//if it is skipped write 0 in the residual and copy the prediction wnd to the decoder wnd
						{
							wnd_copy_cu_2D(et->funcs->sse_copy_16_16, curr_cu_info, &et->prediction_wnd[0], et->decoded_mbs_wnd[curr_depth+1], ALL_COMP);
							wnd_zero_cu_1D(et, curr_cu_info, et->transform_quant_wnd[curr_depth+1]);
							SET_ENC_INFO_BUFFS(et, curr_cu_info, curr_depth/*+(part_size_type!=SIZE_2Nx2N)*/, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);//consolidate in prediction depth
						}
						put_consolidated_info(et, ctu, curr_cu_info, curr_depth);
						bestIsSkip = CBF_ALL(ctu, abs_index, 0) == 0;
						//exchange we keep in currdepth and the aux_depth, where we keep the best buffers
						//we keep best windows in the last depth as we are only processing 2nx2n
//						ptrswap(wnd_t*, et->transform_quant_wnd[curr_depth+1],et->transform_quant_wnd[NUM_QUANT_WNDS-1]);
//						ptrswap(wnd_t*, et->decoded_mbs_wnd[curr_depth+1],et->decoded_mbs_wnd[NUM_DECODED_WNDS-1]);
					}

					if ( uiNoResidual == 0 && CBF_ALL(ctu, abs_index, curr_depth) == 0)
					{
						// If no residual when allowing for one, then set mark to not try case where residual is forced to 0
						mergeCandBuffer[uiMergeCand] = 1;
					}
				}
			}
		}
	}
	get_back_consolidated_info(et, ctu, curr_cu_info, curr_depth);

	//save the information
	curr_cu_info->skipped = bestIsSkip;
	curr_cu_info->inter_mv[REF_PIC_LIST_0] = best_mv;
	curr_cu_info->cost = curr_cu_info->distortion = best_dist;
	curr_cu_info->merge_flag = TRUE;
	curr_cu_info->merge_idx = best_candidate;
	curr_cu_info->sum = best_sum;

	consolidate_prediction_info(et, ctu, NULL, curr_cu_info, curr_cu_info->cost, MAX_COST, FALSE, NULL);
	return best_dist;
}



#define format_vector(buff, ctu_index, part, name, vector_x, vector_y)	sprintf(buff, "\r\nctu:%d, %s, size:%d, abs_index:%d, v:(%d,%d)", ctu_index, name, part->size, part->abs_index, vector_x, vector_y);
#define format_candidates(buff, cand_list)  sprintf(buff, ", candidates:(%d,%d),(%d,%d)", cand_list.mv_candidates[0].mv.hor_vector, cand_list.mv_candidates[0].mv.ver_vector, cand_list.mv_candidates[1].mv.hor_vector, cand_list.mv_candidates[1].mv.ver_vector);
void print_vector(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t *curr_cu_info)
{
	char str[256];
	int n;
	FILE *vector_file = fopen("C:\\Patrones\\vector_file.txt", "a");
	if(curr_cu_info->prediction_mode == INTRA_MODE)
	{
		n = format_vector(str, ctu->ctu_number, curr_cu_info, "intra", 0, 0);
		fwrite(str, sizeof(char), n, vector_file);
	}
	else
	{
		if(curr_cu_info->merge_flag)
		{
			n = format_vector(str, ctu->ctu_number, curr_cu_info, "merge", curr_cu_info->inter_mv[REF_PIC_LIST_0].hor_vector, curr_cu_info->inter_mv[REF_PIC_LIST_0].ver_vector);
			fwrite(str, sizeof(char), n, vector_file);
			str[0] = 0;
			n = format_candidates(str, et->merge_mvp_candidates[REF_PIC_LIST_0]);
			fwrite(str, sizeof(char), n, vector_file);
		}
		else
		{
			n = format_vector(str, ctu->ctu_number, curr_cu_info, "inter", curr_cu_info->inter_mv[REF_PIC_LIST_0].hor_vector, curr_cu_info->inter_mv[REF_PIC_LIST_0].ver_vector);
			fwrite(str, sizeof(char), n, vector_file);
			n = format_candidates(str, et->amvp_candidates[REF_PIC_LIST_0]);
			fwrite(str, sizeof(char), n, vector_file);
		}
	}
	fclose(vector_file);
}


uint32_t motion_inter(henc_thread_t* et, ctu_info_t* ctu)
{
	return motion_inter_full(et, ctu);
}
