/*****************************************************************************
 * hmr_motion_inter.c : homerHEVC encoding library
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
//#include "hmr_profiler.h"

#include "hmr_sse42_functions.h"


int encode_inter_cu(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_cu_info, int depth, PartSize part_size_type, int *curr_sum, int gcnt)//depth = prediction depth
{		
	int ssd_;
	int pred_buff_stride, orig_buff_stride, residual_buff_stride, residual_dec_buff_stride, decoded_buff_stride;
	uint8_t *orig_buff;
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

	quant_wnd = &et->transform_quant_wnd[curr_depth + 1 + (part_size_type!=SIZE_2Nx2N)];
	decoded_wnd = &et->decoded_mbs_wnd[curr_depth + 1  + (part_size_type!=SIZE_2Nx2N)];
//	cbf_buff = et->cbf_buffs[Y_COMP][curr_depth];

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, Y_COMP);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
	orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
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
	et->funcs->transform(et->bit_depth, residual_buff, et->pred_aux_buff, residual_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, curr_part_size_shift, cu_mode, quant_buff);//usamos quant buff como auxiliar

	et->funcs->quant(et, et->pred_aux_buff, quant_buff, curr_scan_mode, curr_depth, Y_COMP, cu_mode, 0, curr_sum, curr_part_size, per, rem);//Si queremos quitar el bit de signo necesitamos hacerlo en dos arrays distintos

	curr_cu_info->inter_cbf[Y_COMP] = (( *curr_sum ? 1 : 0 ) << (curr_depth - depth));// + (part_size_type == SIZE_NxN)));
	curr_cu_info->inter_tr_idx = (curr_depth - depth);// + (part_size_type == SIZE_NxN));


	if(*curr_sum>0)//curr_cu_info->size)
	{
		uint32_t ssd_zero;
		int16_t zero_buff[256];
		memset(zero_buff,0,sizeof(zero_buff));
		ssd_zero = ssd_16(residual_buff, residual_buff_stride, zero_buff, 0, curr_part_size);

		et->funcs->inv_quant(et, quant_buff, iquant_buff, curr_depth, Y_COMP, 0, curr_part_size, per, rem);

		//1D ->2D buffer
		et->funcs->itransform(et->bit_depth, residual_dec_buff, iquant_buff, residual_buff_stride, curr_part_size, curr_part_size, cu_mode, et->pred_aux_buff);

		ssd_ = ssd_16(residual_buff, residual_buff_stride, residual_dec_buff, residual_buff_stride, curr_part_size);

		if(ssd_zero < clip((200./et->ed->avg_dist),1.01,1.25)*ssd_)
		{
			memset(quant_buff, 0, curr_part_size*curr_part_size*sizeof(quant_buff[0]));
			*curr_sum = 0;
			curr_cu_info->inter_cbf[Y_COMP] = 0;//
			et->funcs->reconst(pred_buff, pred_buff_stride, quant_buff, 0, decoded_buff, decoded_buff_stride, curr_part_size);//quant buff is full of zeros - a memcpy could do
		}
		else
			et->funcs->reconst(pred_buff, pred_buff_stride, residual_dec_buff, residual_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	}
	else
	{
		ssd_ = ssd_16(residual_buff, residual_buff_stride, quant_buff, 0, curr_part_size);

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
	picture_t *currpict = &et->ed->current_pict;
	slice_t *currslice = &currpict->slice;
	int pred_buff_stride, orig_buff_stride, residual_buff_stride, residual_dec_buff_stride, decoded_buff_stride;
	uint8_t *orig_buff;
	int16_t *pred_buff, *residual_buff, *residual_dec_buff, *quant_buff, *iquant_buff, *decoded_buff;
//	uint8_t *cbf_buff = NULL;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL;
	int inv_depth, diff;//, is_filtered;
	int cu_mode = REG_DCT;//if !IsIntra(abs_index)cu_mode = REG_DCT;
	int original_depth = curr_cu_info->depth;
	cu_partition_info_t* processing_partition_info = (curr_cu_info->size_chroma!=2)?curr_cu_info:curr_cu_info->parent;
	int curr_depth = processing_partition_info->depth;
	int curr_part_x = processing_partition_info->x_position_chroma;
	int curr_part_y = processing_partition_info->y_position_chroma;
	int curr_part_size = processing_partition_info->size_chroma;
	int curr_part_size_shift = et->max_cu_size_shift-curr_depth-1;//420
	int curr_scan_mode = find_scan_mode(TRUE, TRUE, curr_part_size, cu_mode, 0);
	int chr_qp_offset = et->ed->chroma_qp_offset;
	int qp_chroma = chroma_scale_conversion_table[clip(curr_cu_info->qp+chr_qp_offset,0,57)];
	double weight = pow( 2.0, (currslice->qp-chroma_scale_conversion_table[clip(currslice->qp+chr_qp_offset,0,57)])/3.0 );

	int per = qp_chroma/6;
	int rem = qp_chroma%6;

	quant_wnd = &et->transform_quant_wnd[original_depth+1+(part_size_type!=SIZE_2Nx2N)];
	decoded_wnd = &et->decoded_mbs_wnd[original_depth+1+(part_size_type!=SIZE_2Nx2N)];
//	cbf_buff = et->cbf_buffs[component][original_depth];

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, component);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, component);
	orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
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
	et->funcs->transform(et->bit_depth, residual_buff, et->pred_aux_buff, residual_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, curr_part_size_shift, cu_mode, quant_buff);//usamos quant buff como auxiliar

	et->funcs->quant(et, et->pred_aux_buff, quant_buff, curr_scan_mode, curr_depth, component, cu_mode, 0, curr_sum, curr_part_size, per, rem);//Si queremos quitar el bit de signo necesitamos hacerlo en dos arrays distintos


	curr_cu_info->inter_cbf[component] = (( *curr_sum ? 1 : 0 ) << (original_depth-depth));//+(part_size_type==SIZE_NxN)));

	if(*curr_sum>0)//curr_cu_info->size_chroma)
	{
		uint32_t ssd_zero;
		int16_t zero_buff[256];
		memset(zero_buff,0,sizeof(zero_buff));
		ssd_zero = weight*ssd_16(residual_buff, residual_buff_stride, zero_buff, 0, curr_part_size);

		et->funcs->inv_quant(et, quant_buff, iquant_buff, curr_depth, component, 0, curr_part_size, per, rem);

		//1D ->2D buffer
		et->funcs->itransform(et->bit_depth, residual_dec_buff, iquant_buff, residual_buff_stride, curr_part_size, curr_part_size, cu_mode, et->pred_aux_buff);
		ssd_ = weight*ssd_16(residual_buff, residual_buff_stride, residual_dec_buff, residual_buff_stride, curr_part_size);

		if(ssd_zero < clip((200./et->ed->avg_dist),1.01,1.25)*ssd_)
		{
			memset(quant_buff, 0, curr_part_size*curr_part_size*sizeof(quant_buff[0]));
			*curr_sum = 0;
			curr_cu_info->inter_cbf[component] = 0;//
			et->funcs->reconst(pred_buff, pred_buff_stride, quant_buff, 0, decoded_buff, decoded_buff_stride, curr_part_size);//quant buff is full of zeros - a memcpy could do
		}
		else
			et->funcs->reconst(pred_buff, pred_buff_stride, residual_dec_buff, residual_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	}
	else
	{
		ssd_ = weight*ssd_16(residual_buff, residual_buff_stride, quant_buff, 0, curr_part_size);

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


//TComInterpolationFilter::filterCopy(Int bitDepth, const Pel *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast)
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
	int x_fraction;// = mv->hor_vector&0x3;
	int y_fraction;// = mv->ver_vector&0x3;

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
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width+1, height+filter_size, INTERPOLATE_HOR, TRUE, FALSE);

	//half pixel horizontal interpolation 
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width+1, height+filter_size, INTERPOLATE_HOR, TRUE, FALSE);

	//copy original samples
	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], Y_COMP);
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + half_filter_size*src_stride+1;
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[0][0], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);

	//half pixel vertical interpolation of original pixels
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride+1;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height+1, INTERPOLATE_VERT, FALSE, TRUE);
  
	//copy half pixel horizontal interpolation
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width+1, height, INTERPOLATE_VERT, FALSE, TRUE);
	
	//half pixel vertical interpolation of horizontaly interpolated pixels
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width+1, height+1, INTERPOLATE_VERT, FALSE, TRUE);
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
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);

	src = reference_buff - half_filter_size*src_stride - 1;
	//horizontal filter 3,0
	if(half_pel_mv->ver_vector>0)
		src += src_stride;
	if(half_pel_mv->hor_vector>0)
		src += 1;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);

	//vertical filter 1,1
	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[1][1], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	if(half_pel_mv->ver_vector==0)
		src += src_stride;
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

	//vertical filter 1,3
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);

	if(half_pel_mv->ver_vector != 0)
	{
		//vertical filter 1,2
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector==0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//vertical filter 3,2
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector==0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
	}
	else
	{
		//vertical filter 1,0
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//vertical filter 3,0
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);    		
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
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//filter 2,3
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->hor_vector>0)
			src += 1;
		if(half_pel_mv->ver_vector>0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
	}
	else
	{
		//filter 0,1
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector>=0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//filter 0,3
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector>0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
	}

	//filter 3,1
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	if(half_pel_mv->ver_vector==0)
		src += src_stride;
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

	//filter 3,3
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
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
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);

	src = reference_buff - half_filter_size*src_stride - 1;
	//horizontal filter 3,0
	if(half_pel_mv->ver_vector>0)
		src += src_stride;
	if(half_pel_mv->hor_vector>0)
		src += 1;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);

	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[1][1], Y_COMP);

	if(zero_curr_best)
	{
		if(dir_x<0)
		{
			//vertical filter 3,0
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);    					
			processed++;
		}
		else if(dir_x>0)
		{
			//vertical filter 1,0
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);			
			processed++;		
		}

		if(dir_y>0)
		{
			//filter 0,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector>=0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}
		else if(dir_y<0)
		{
			//filter 0,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector>0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x<=0 && dir_y<=0)
		{
			//filter 3,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;		
		}

		if(dir_x<=0 && dir_y>=0)
		{
			//filter 3,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x>=0 && dir_y<=0)
		{
			//vertical filter 1,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x>=0 && dir_y>=0)
		{
			//vertical filter 1,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
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
				et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);			
				processed++;
			}
			else if(dir_x<0)
			{
				//vertical filter 3,0
				src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
				dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
				et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);    					
				processed++;
			}
			else
			{
				int iiiii=0;			
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
				et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
				processed++;
			}
			else if(dir_y<0)
			{
				//filter 0,3
				src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
				dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
				if(half_pel_mv->ver_vector>0)
					src += src_stride;
				et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
				processed++;
			}
			else
			{
				int iiiii=0;			
			}
		}

		if(dir_x>0 && dir_y>0)
		{
			//vertical filter 1,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x>0 && dir_y<0)
		{
			//vertical filter 1,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x<0 && dir_y>0)
		{
			//filter 3,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}

		if(dir_x<0 && dir_y<0)
		{
			//filter 3,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}


		if(dir_x>=0 && dir_y!=0)
		{
			//vertical filter 1,2
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			processed++;
		}
		
		if(dir_x<=0 && dir_y!=0)
		{
			//vertical filter 3,2
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
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
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
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
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
			processed++;
		}
	}
	if(processed!=3)
		printf("processed=%d", processed);

}


void hmr_half_pixel_estimation_luma_fast(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, cu_partition_info_t* curr_cu_info, int width, int height, int curr_part_size_shift, motion_vector_t *mv, int dir_x, int dir_y)
{
	int is_bi_predict = 0;
	int x_fraction;// = mv->hor_vector&0x3;
	int y_fraction;// = mv->ver_vector&0x3;

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
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width+1, height+filter_size, INTERPOLATE_HOR, TRUE, FALSE);

		//half pixel vertical interpolation of original pixels
		src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], Y_COMP);
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride+1;
		dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[2][0], Y_COMP);
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height+1, INTERPOLATE_VERT, FALSE, TRUE);
/*
		src_stride = reference_buff_stride;
		src = reference_buff - half_filter_size*src_stride - 1  +  (half_filter_size-1)*src_stride+1;
		dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[2][0], Y_COMP);
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height+1, INTERPOLATE_VERT, TRUE, TRUE);
*/	}

	//half pixel horizontal interpolation 
	src_stride = reference_buff_stride;
	src = reference_buff - half_filter_size*src_stride - 1;
	dst_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[2], Y_COMP);
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width+1, height+filter_size, INTERPOLATE_HOR, TRUE, FALSE);

	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[2], Y_COMP);
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[0][2], Y_COMP);

	//copy half pixel horizontal interpolation
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width+1, height, INTERPOLATE_HOR, FALSE, TRUE);

	//half pixel vertical interpolation of horizontaly interpolated pixels
	src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
	dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
	et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width+1, height+1, INTERPOLATE_VERT, FALSE, TRUE);
}

void hmr_quarter_pixel_estimation_luma_fast_(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, cu_partition_info_t* curr_cu_info, int width, int height, int curr_part_size_shift, motion_vector_t *half_pel_mv, int dir_x, int dir_y)
{
	int filter_size = NTAPS_LUMA;
	int half_filter_size = (filter_size>>1);
	int src_stride = reference_buff_stride;
	int16_t *src;
	int16_t *dst;
	int dst_stride;
	int curr_part_x = curr_cu_info->x_position, curr_part_y = curr_cu_info->y_position;
	int ext_height = (half_pel_mv->ver_vector == 0) ? height + filter_size : height + filter_size-1;

	dst_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);

	if(dir_x<=0)
	{
		src = reference_buff - half_filter_size*src_stride - 1;
		//horizontal filter 3,0
		if(half_pel_mv->ver_vector>0)
			src += src_stride;
		if(half_pel_mv->hor_vector>0)
			src += 1;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);
	}

	if(dir_x>=0)
	{
		//horizontal filter 1,0
		src = reference_buff - half_filter_size*src_stride - 1;
		if(half_pel_mv->ver_vector>0)
			src += src_stride;
		if(half_pel_mv->hor_vector>=0)
			src += 1;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, ext_height, INTERPOLATE_HOR, TRUE, FALSE);
	}

	src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);
	dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[3][2], Y_COMP);


	if(dir_y==0)
	{
		//2,3 final
		//filter 2,3
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->hor_vector>0)
			src += 1;
		if(half_pel_mv->ver_vector>0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);    

		//2,1 final		
		//filter 2,1
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->hor_vector>0)
			src += 1;
		if(half_pel_mv->ver_vector>=0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		if(dir_x<0)
		{			
			//3,0 final
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);    		
		}
		else if(dir_x>0)
		{
			//1,0 final
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[0][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 0, width, height, INTERPOLATE_VERT, FALSE, TRUE);
		}
	}
	else if(dir_x==0)
	{
		//1,2 final
		//vertical filter 1,2
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector==0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		//3,2 final
		//vertical filter 3,2
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector==0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);    

		if(dir_y<0)
		{
			//0,3 final
			//filter 0,3
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector>0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
		}
		else
		{
			//0,1 final
			//filter 0,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride + 1;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][0], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector>=0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
		}
	}
	else if(dir_x<0)
	{
		//3,2 final
		//vertical filter 3,2
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector==0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);    

		if(dir_y<0)
		{
			//2,3 final
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->hor_vector>0)
				src += 1;
			if(half_pel_mv->ver_vector>0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
			//3,3 final
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);
		}
		else
		{
			//2,1 final
			//filter 2,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->hor_vector>0)
				src += 1;
			if(half_pel_mv->ver_vector>=0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
			//3,1 final
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][3], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
		}
	}
	else //if(dir_x>0)
	{
		//1,2 final
		src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
		dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[2][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
		if(half_pel_mv->ver_vector==0)
			src += src_stride;
		et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 2, width, height, INTERPOLATE_VERT, FALSE, TRUE);

		if(dir_y<0)
		{
			//1,3 final
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);

			//2,3 final	
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[3][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->hor_vector>0)
				src += 1;
			if(half_pel_mv->ver_vector>0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 3, width, height, INTERPOLATE_VERT, FALSE, TRUE);    
		}
		else //if(dir_y<0)
		{
			//1,1 final
			src_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[1], Y_COMP);
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst_stride = WND_STRIDE_2D(et->filtered_block_wnd[1][1], Y_COMP);
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][1], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->ver_vector==0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);

			//2,1 final
			//filter 2,1
			src = WND_POSITION_2D(int16_t *, et->filtered_block_temp_wnd[2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width) + (half_filter_size-1)*src_stride;
			dst = WND_POSITION_2D(int16_t *, et->filtered_block_wnd[1][2], Y_COMP, curr_part_x, curr_part_y, 0, et->ctu_width);
			if(half_pel_mv->hor_vector>0)
				src += 1;
			if(half_pel_mv->ver_vector>=0)
				src += src_stride;
			et->funcs->interpolate_luma(src, src_stride, dst, dst_stride, 1, width, height, INTERPOLATE_VERT, FALSE, TRUE);
		}	
	}
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


static const int diamond_small[][2] = {{0,-1},{-1,0},{1,0},{0,1}};
static const int diamond_big[][2] = {{0,-2},{-1,-1},{1,-1},{-2,0},{2,0},{-1,1},{1,1},{0,2}};
static const int square[][2] = {{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0}};

uint32_t hmr_motion_estimation_HM(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_cu_info, uint8_t *orig_buff, int orig_buff_stride, int16_t *reference_buff, int reference_buff_stride, int curr_part_global_x, 
						int curr_part_global_y, int init_x, int init_y, int curr_part_size, int curr_part_size_shift, int search_range_x, int search_range_y, int frame_size_x, int frame_size_y, motion_vector_t *mv)
{
	int i,j; 
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
//	curr_best_x = 0;
//	curr_best_y = 0;

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
//	best_x = (curr_best_x & 0xfffffffe)+1;
//	best_y = (curr_best_y & 0xfffffffe)+1;
//	curr_best_x = 0;
//	curr_best_y = 0;

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



static const int quarter_search_points_top_left[3][2] =
{
	{ -1, -2 }, // 0
	{ -1, -1 }, // 1
	{ -2, -1 }, // 2
};

static const int quarter_search_points_top[3][2] =
{
	{ -1, -2 }, // 0
	{  0, -1 }, // 1
	{  1, -2 }, // 2
};

static const int quarter_search_points_top_right[3][2] =
{
	{  1, -2 }, // 0
	{  1, -1 }, // 1
	{  2, -1 }, // 2
};

static const int quarter_search_points_right[3][2] =
{
	{  2, -1 }, // 0
	{  1,  0 }, // 1
	{  2,  1 }, // 2
};

static const int quarter_search_points_bottom_right[3][2] =
{
	{  2,  1 }, // 0
	{  1,  1 }, // 1
	{  1,  2 }, // 2
};

static const int quarter_search_points_bottom[3][2] =
{
	{  1,  2 }, // 0
	{  0,  1 }, // 1
	{ -1,  2 }, // 2
};

static const int quarter_search_points_bottom_left[3][2] =
{
	{  2, -1 }, // 0
	{ -1,  1 }, // 1
	{ -2,  1 }, // 2
};

static const int quarter_search_points_left[3][2] =
{
	{ -2,  1 }, // 0
	{ -1,  0 }, // 1
	{ -2, -1 }, // 2
};


uint32_t hmr_motion_estimation(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_cu_info, uint8_t *orig_buff, int orig_buff_stride, int16_t *reference_buff, int reference_buff_stride, int curr_part_global_x, 
						int curr_part_global_y, int init_x, int init_y, int curr_part_size, int curr_part_size_shift, int search_range_x, int search_range_y, int frame_size_x, int frame_size_y, motion_vector_t *mv, motion_vector_t *subpix_mv, uint32_t threshold)
{
	int i,j, l; 
	int xlow, xhigh, ylow, yhigh;
	uint32_t prev_sad, prev_best_sad, prev_best_rd, curr_best_sad, curr_best_rd, best_sad, best_rd;
	int prev_best_x, prev_best_y, curr_best_x, curr_best_y, best_x, best_y;
	int dist = 1;
	int end;
	mv_candiate_list_t	*mv_candidate_list = &et->mv_search_candidates;//&et->mv_candidates[REF_PIC_LIST_0];
	int mv_cost = 0;

	subpix_mv->hor_vector = 0;
	subpix_mv->ver_vector = 0;

	threshold = 0;


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

	mv->hor_vector = curr_best_x<<2;
	mv->ver_vector = curr_best_y<<2;
	mv_cost = select_mv_candidate_fast(et, curr_cu_info, REF_PIC_LIST_0, mv);
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

	for(i=0;i<mv_candidate_list->num_mv_candidates;i++)
	{
		uint32_t curr_sad, curr_rd;
		int curr_x = mv_candidate_list->mv_candidates[i].hor_vector>>2;
		int curr_y = mv_candidate_list->mv_candidates[i].ver_vector>>2;

		if(curr_x == 0 && curr_y == 0)
		{
			continue;
		}
		if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
		{
			curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);
			mv->hor_vector = curr_x<<2;
			mv->ver_vector = curr_y<<2;
			mv_cost = select_mv_candidate_fast(et, curr_cu_info, REF_PIC_LIST_0, mv);
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
			mv_cost = select_mv_candidate_fast(et, curr_cu_info, REF_PIC_LIST_0, mv);
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

	for(l = 0; l < 1; l++)
	{
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
			for(i=0;i<sizeof(diamond_big)/sizeof(diamond_big[0]);i++)
			{
				int curr_x = best_x+diamond_big[i][0]*dist;
				int curr_y = best_y+diamond_big[i][1]*dist;
				uint32_t curr_sad, curr_rd;

				if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
				{
					curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);

					mv->hor_vector = curr_x<<2;
					mv->ver_vector = curr_y<<2;
					mv_cost = select_mv_candidate_fast(et, curr_cu_info, REF_PIC_LIST_0, mv);
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


	while(1)
	{
		for(i=0;i<sizeof(square)/sizeof(square[0]);i++)
		{
			int curr_x = best_x+square[i][0];
			int curr_y = best_y+square[i][1];
			uint32_t curr_sad, curr_rd;

			if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
			{
				curr_sad = et->funcs->sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_cu_info->size);
				mv->hor_vector = curr_x<<2;
				mv->ver_vector = curr_y<<2;
				mv_cost = select_mv_candidate_fast(et, curr_cu_info, REF_PIC_LIST_0, mv);
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



void hmr_motion_compensation_luma(henc_thread_t *et, ctu_info_t *ctu, cu_partition_info_t* curr_cu_info, int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int width, int height, int curr_part_size_shift, motion_vector_t *mv)
{
	int is_bi_predict = 0;
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
		et->funcs->interpolate_luma(src, src_stride, pred_buff, pred_buff_stride, y_fraction, width, height, INTERPOLATE_VERT, TRUE, !is_bi_predict);
	}
	else if (y_fraction  == 0)
	{
		et->funcs->interpolate_luma(src, src_stride, pred_buff, pred_buff_stride, x_fraction, width, height, INTERPOLATE_HOR, TRUE, !is_bi_predict);
	}
	else
	{
		int filter_size = NTAPS_LUMA;
		int half_filter_size = filter_size>>1;
		int16_t *temp_buff = WND_DATA_PTR(int16_t *, et->filtered_block_temp_wnd[0], Y_COMP);
		int temp_buff_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], Y_COMP);

		//horizontal
		et->funcs->interpolate_luma(src - (half_filter_size-1)*src_stride, src_stride, temp_buff, temp_buff_stride, x_fraction, width, height+filter_size-1, INTERPOLATE_HOR, TRUE, FALSE);
		//vertical filter 
		et->funcs->interpolate_luma(temp_buff + (half_filter_size-1)*temp_buff_stride, temp_buff_stride, pred_buff, pred_buff_stride, y_fraction, width, height, INTERPOLATE_VERT, FALSE, !is_bi_predict);
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





void hmr_motion_compensation_chroma(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int curr_part_size, int curr_part_size_shift, motion_vector_t *mv)
{
	int is_bi_predict = 0;
	int x_fraction = mv->hor_vector&0x7;
	int y_fraction = mv->ver_vector&0x7;

	int x_vect = mv->hor_vector>>3;
	int y_vect = mv->ver_vector>>3;
	reference_buff+= y_vect*reference_buff_stride+x_vect;


	if(x_fraction==0 && y_fraction==0)
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
	else if(x_fraction == 0)
	{
		//vertical filter 
		et->funcs->interpolate_chroma(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, y_fraction, curr_part_size, curr_part_size, INTERPOLATE_VERT, TRUE, !is_bi_predict);
	}
	else if(y_fraction == 0)
	{
		//horizontal filter 
		et->funcs->interpolate_chroma(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, x_fraction, curr_part_size, curr_part_size, INTERPOLATE_HOR, TRUE, !is_bi_predict);	
	}
	else //if(x_fraction!=0 && y_fraction!=0)
	{
		int filter_size = NTAPS_CHROMA;
		int half_filter_size = filter_size>>1;
		int16_t *temp_buff = WND_DATA_PTR(int16_t *, et->filtered_block_temp_wnd[0], U_COMP);//this is only for temporal results so the component does not care
		int temp_buff_stride = WND_STRIDE_2D(et->filtered_block_temp_wnd[0], U_COMP);
		//horizontal
		et->funcs->interpolate_chroma(reference_buff - (half_filter_size-1)*reference_buff_stride, reference_buff_stride, temp_buff, temp_buff_stride, x_fraction, curr_part_size, curr_part_size+filter_size+1, INTERPOLATE_HOR, TRUE, FALSE);
		//vertical filter 
		et->funcs->interpolate_chroma(temp_buff + (half_filter_size-1)*temp_buff_stride, temp_buff_stride, pred_buff, pred_buff_stride, y_fraction, curr_part_size, curr_part_size, INTERPOLATE_VERT, FALSE, !is_bi_predict);
	}
}


//fillMvpCand
void get_mv_candidates(henc_thread_t* et, slice_t *currslice, ctu_info_t* ctu, cu_partition_info_t *curr_cu_info, int ref_pic_list, int ref_idx, PartSize part_size_type)//get candidates for motion search from the neigbour CUs
{
	mv_candiate_list_t	*mv_candidate_list = &et->mv_candidates[ref_pic_list];
	ctu_info_t	*ctu_left=NULL, *ctu_left_bottom=NULL, *ctu_top=NULL, *ctu_top_right=NULL, *ctu_top_left=NULL;
	uint	part_idx_lb, part_idx_l, aux_part_idx;
	int	added = FALSE;
	int	num_partitions_height = curr_cu_info->size>>2;
	int	num_partitions_width = curr_cu_info->size>>2;
	int	max_partitions_width = et->max_cu_size>>2;
	uint num_partitions_mask = curr_cu_info->num_part_in_cu-1;
	int abs_index_lb = et->ed->raster2abs_table[curr_cu_info->raster_index + max_partitions_width*(num_partitions_height-1)];
	int abs_index_tl = et->ed->raster2abs_table[curr_cu_info->raster_index];
	int abs_index_tr = et->ed->raster2abs_table[curr_cu_info->raster_index + num_partitions_width-1];
	cu_partition_info_t *partition_info_lb = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_lb;
	cu_partition_info_t *partition_tl = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_tl;
	cu_partition_info_t *partition_tr = &ctu->partition_list[et->partition_depth_start[et->max_cu_depth]]+abs_index_tr;
	int left_pixel_x = ctu->x[Y_COMP]+curr_cu_info->x_position+curr_cu_info->size;

	mv_candidate_list->num_mv_candidates = 0;

	//get spatial candidates
	partition_info_lb->left_bottom_neighbour = curr_cu_info->left_bottom_neighbour;
	ctu_left_bottom = get_pu_left_bottom(et, ctu, partition_info_lb, &part_idx_lb);

	if(ctu_left_bottom!=NULL && ctu_left_bottom->mv_ref_idx[ref_pic_list][part_idx_lb]>=0)// && currslice->ref_pic_list[ref_pic_list][ref_idx]->temp_info.poc == currslice->re)
	{
		mv_candidate_list->mv_candidates[mv_candidate_list->num_mv_candidates++] = ctu_left_bottom->mv_ref[ref_pic_list][part_idx_lb];
		added = TRUE;
	}
	else
	{
		ctu_left = get_pu_left(ctu, partition_info_lb, &part_idx_l);//ctu->ctu_left;
		if(ctu_left!=NULL && ctu_left->mv_ref_idx[ref_pic_list][part_idx_l]>=0)
		{
			mv_candidate_list->mv_candidates[mv_candidate_list->num_mv_candidates++] = ctu_left->mv_ref[ref_pic_list][part_idx_l];
		}
	}

	added = FALSE;
	partition_tr->top_right_neighbour = curr_cu_info->top_right_neighbour;
	ctu_top_right = get_pu_top_right(ctu, partition_tr, &aux_part_idx);
//	ctu_top_right = get_pu_top_right(ctu, curr_cu_info, &aux_part_idx);
	if(ctu_top_right!=NULL && ctu_top_right->mv_ref_idx[ref_pic_list][aux_part_idx]>=0)
	{
		mv_candidate_list->mv_candidates[mv_candidate_list->num_mv_candidates++] = ctu_top_right->mv_ref[ref_pic_list][aux_part_idx];
		added = TRUE;		
	}

	if(!added)
	{
		ctu_top = get_pu_top(ctu, partition_tr, &aux_part_idx, 0);
		if(ctu_top!=NULL && ctu_top->mv_ref_idx[ref_pic_list][aux_part_idx]>=0)
		{
			mv_candidate_list->mv_candidates[mv_candidate_list->num_mv_candidates++] = ctu_top->mv_ref[ref_pic_list][aux_part_idx];
			added = TRUE;
		}
		else
			added = FALSE;
	}

	if(!added)
	{
		ctu_top_left = get_pu_top_left(ctu, partition_tl, &aux_part_idx);
		if(ctu_top_left!=NULL && ctu_top_left->mv_ref_idx[ref_pic_list][aux_part_idx]>=0)
		{
			mv_candidate_list->mv_candidates[mv_candidate_list->num_mv_candidates++] = ctu_top_left->mv_ref[ref_pic_list][aux_part_idx];
			added = TRUE;
		}
		else
			added = FALSE;
	}

	if((ctu_left_bottom!=NULL && ctu_left_bottom->pred_mode[part_idx_lb] != INTRA_MODE)|| (ctu_left!=NULL && ctu_left->pred_mode[part_idx_l] != INTRA_MODE))
	{
		//reorder		¿?
	}

	if(mv_candidate_list->num_mv_candidates==2 && mv_candidate_list->mv_candidates[0].hor_vector == mv_candidate_list->mv_candidates[1].hor_vector && mv_candidate_list->mv_candidates[0].ver_vector == mv_candidate_list->mv_candidates[1].ver_vector)
	{
		mv_candidate_list->num_mv_candidates = 1;
	}

	//get temporal candidates

	//.....
	if(mv_candidate_list->num_mv_candidates>AMVP_MAX_NUM_CANDS)
		mv_candidate_list->num_mv_candidates=AMVP_MAX_NUM_CANDS;

	while(mv_candidate_list->num_mv_candidates<AMVP_MAX_NUM_CANDS)
	{
		mv_candidate_list->mv_candidates[mv_candidate_list->num_mv_candidates].hor_vector = 0;
		mv_candidate_list->mv_candidates[mv_candidate_list->num_mv_candidates++].ver_vector = 0;
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

//xCheckBestMVP
int select_mv_candidate(henc_thread_t* et, cu_partition_info_t* curr_cu_info, int ref_pic_list, motion_vector_t *mv)
{
	mv_candiate_list_t	*mv_candidate_list = &et->mv_candidates[ref_pic_list];
	int idx;
	int best_cost = INT_MAX, best_idx = 0;

	for (idx = 0; idx < mv_candidate_list->num_mv_candidates; idx++)
	{
#ifdef COMPUTE_AS_HM
		int cost_mvx = 30*sqrt((float)abs(mv_candidate_list->mv_candidates[idx].hor_vector - mv->hor_vector));
		int cost_mvy = 30*sqrt((float)abs(mv_candidate_list->mv_candidates[idx].ver_vector - mv->ver_vector));
#else
		int cost_mvx = curr_cu_info->qp*squareRoot((float)abs(mv_candidate_list->mv_candidates[idx].hor_vector - mv->hor_vector));
		int cost_mvy = curr_cu_info->qp*squareRoot((float)abs(mv_candidate_list->mv_candidates[idx].ver_vector - mv->ver_vector));
#endif
		int cost = 3+cost_mvx+cost_mvy;

		if(best_cost>cost)
		{
			best_cost = cost;
			best_idx = idx;		
		}
	}
	curr_cu_info->best_candidate_idx[ref_pic_list] = best_idx;
	curr_cu_info->best_dif_mv->hor_vector = mv->hor_vector-mv_candidate_list->mv_candidates[best_idx].hor_vector;
	curr_cu_info->best_dif_mv->ver_vector = mv->ver_vector-mv_candidate_list->mv_candidates[best_idx].ver_vector;
	return best_cost;
}

int select_mv_candidate_fast(henc_thread_t* et, cu_partition_info_t* curr_cu_info, int ref_pic_list, motion_vector_t *mv)
{
#ifdef COMPUTE_AS_HM
	return 0;
#else
	mv_candiate_list_t	*mv_candidate_list = &et->mv_candidates[ref_pic_list];
	int idx;
	int best_cost = INT_MAX, best_idx = 0;

	for (idx = 0; idx < mv_candidate_list->num_mv_candidates; idx++)
	{
		double correction = calc_mv_correction(curr_cu_info->qp, et->ed->avg_dist);//.25+et->ed->avg_dist*et->ed->avg_dist/5000000.;
//		int cost_mvx = 10*abs(mv_candidate_list->mv_candidates[idx].hor_vector - mv->hor_vector);
//		int cost_mvy = 10*abs(mv_candidate_list->mv_candidates[idx].ver_vector - mv->ver_vector);
//		int cost_mvx = curr_cu_info->qp/clip((3500000/(et->ed->avg_dist*et->ed->avg_dist)),.35,4.)*squareRoot((float)abs(mv_candidate_list->mv_candidates[idx].hor_vector - mv->hor_vector));
//		int cost_mvy = curr_cu_info->qp/clip((3500000/(et->ed->avg_dist*et->ed->avg_dist)),.35,4.)*squareRoot((float)abs(mv_candidate_list->mv_candidates[idx].ver_vector - mv->ver_vector));
		int cost_mvx = correction*((float)abs(mv_candidate_list->mv_candidates[idx].hor_vector - mv->hor_vector));
		int cost_mvy = correction*((float)abs(mv_candidate_list->mv_candidates[idx].ver_vector - mv->ver_vector));


		int cost = 3+cost_mvx+cost_mvy;

		if(best_cost>cost)
		{
			best_cost = cost;
//			best_idx = idx;		
		}
	}
	return best_cost;
#endif
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
	}																																						\
}


int hmr_cu_motion_estimation(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type, uint threshold)
{
	int k;
	int cu_mode;
	double distortion = 0.;

	uint32_t sad, cost, best_cost;
	slice_t *currslice = &et->ed->current_pict.slice;
	ctu_info_t *ctu_rd = et->ctu_rd;
	int orig_buff_stride, reference_buff_stride;
	int orig_buff_stride_chroma, reference_buff_stride_chroma;
	uint8_t  *orig_buff, *orig_buff_u, *orig_buff_v;
	int16_t  *reference_buff_cu_position, *reference_buff_cu_position_u, *reference_buff_cu_position_v;
	wnd_t *reference_wnd=NULL;//, *resi_wnd = NULL;
	int curr_part_size, curr_part_size_shift;
	int curr_part_size_chroma, curr_part_size_shift_chroma;
	int curr_part_x, curr_part_y, curr_part_global_x, curr_part_global_y;
	int curr_part_x_chroma, curr_part_y_chroma, curr_part_global_x_chroma, curr_part_global_y_chroma;
	int curr_depth = depth;
	cu_partition_info_t *parent_part_info;
	cu_partition_info_t *curr_cu_info;
	int curr_sum = 0, best_sum;
	int num_part_in_cu;
	int partition_cost;
	int cu_min_tu_size_shift;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int max_tr_depth, max_tr_processing_depth;
	int initial_state, end_state;
//	int cbf_split[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
//	int acc_cost[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int bitcost_cu_mode;
	int log2cu_size;
	int cu_x_position;
	motion_vector_t mv, subpix_mv;
	int ref_idx = 0;
	int num_partitions, npart, part_incr = 1;
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
		part_incr = 1;
	}

	sad = 0;
	for(npart=0;npart<num_partitions;npart+=part_incr)
	{
		int i;
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
		orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
		orig_buff_stride_chroma = WND_STRIDE_2D(et->curr_mbs_wnd, CHR_COMP);
		orig_buff_u = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
		orig_buff_v = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

		reference_wnd = &currslice->ref_pic_list[REF_PIC_LIST_0][ref_idx]->img;//[0] up to now we only use one reference	
		reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
		reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
		reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, CHR_COMP);
		reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
		reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);

		get_mv_candidates(et, currslice, ctu, curr_cu_info, REF_PIC_LIST_0, ref_idx, part_size_type);//get candidates for motion search from the neigbour CUs

#ifdef COMPUTE_AS_HM
		sad += hmr_motion_estimation_HM(et, ctu, curr_cu_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, curr_part_global_x, curr_part_global_y, 0, 0, curr_part_size, curr_part_size_shift, 64, 64, et->pict_width[Y_COMP], et->pict_height[Y_COMP], &mv);	
		select_mv_candidate(et, curr_cu_info, REF_PIC_LIST_0, &mv);
#else
		et->mv_search_candidates.num_mv_candidates = 0;
//		get_mv_search_candidates(et, currslice, ctu, curr_cu_info, REF_PIC_LIST_0, ref_idx, part_size_type);//get candidates for motion search from the neigbour CUs
		for(i=0;i<et->mv_candidates[REF_PIC_LIST_0].num_mv_candidates;i++)
		{
			motion_vector_t mv = et->mv_candidates[REF_PIC_LIST_0].mv_candidates[i];
			if(mv.hor_vector!=0 && mv.ver_vector!=0)
			{
				et->mv_search_candidates.mv_candidates[et->mv_search_candidates.num_mv_candidates++] = mv;	
			}
		}
		if(curr_cu_info->parent && curr_cu_info->parent->inter_mv[REF_PIC_LIST_0].hor_vector!=0 && curr_cu_info->parent->inter_mv[REF_PIC_LIST_0].ver_vector!=0)
		{
			et->mv_search_candidates.mv_candidates[et->mv_search_candidates.num_mv_candidates++] = curr_cu_info->parent->inter_mv[REF_PIC_LIST_0];	
		}
		
		sad = hmr_motion_estimation(et, ctu, curr_cu_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, curr_part_global_x, curr_part_global_y, 0, 0, curr_part_size, curr_part_size_shift, 64, 64, et->pict_width[Y_COMP], et->pict_height[Y_COMP], &mv, &subpix_mv, threshold);	
		mv_cost = select_mv_candidate_fast(et, curr_cu_info, REF_PIC_LIST_0, &mv);
//		mv_cost = select_mv_candidate(et, curr_cu_info, REF_PIC_LIST_0, &mv);
		curr_cu_info->subpix_mv[REF_PIC_LIST_0] = subpix_mv;
#endif
		//set mvs and ref_idx
		curr_cu_info->inter_mv[REF_PIC_LIST_0] = mv;
		curr_cu_info->inter_ref_index[REF_PIC_LIST_0] = ref_idx;

//		SET_INTER_MV_BUFFS(et, ctu, curr_cu_info, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);
//		memset(&ctu->mv_ref_idx[REF_PIC_LIST_0][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_0], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));
		curr_cu_info++;
	}
	return sad + 1*mv_cost;
}

int predict_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type)
{
	int k;
	int cu_mode;
	double distortion = 0.;

	uint32_t sad, cost, best_cost;
	slice_t *currslice = &et->ed->current_pict.slice;
	ctu_info_t *ctu_rd = et->ctu_rd;
	int pred_buff_stride, orig_buff_stride, reference_buff_stride, residual_buff_stride;
	int pred_buff_stride_chroma, orig_buff_stride_chroma, reference_buff_stride_chroma, residual_buff_stride_chroma;
	uint8_t  *orig_buff, *orig_buff_u, *orig_buff_v;
	int16_t  *reference_buff_cu_position, *reference_buff_cu_position_u, *reference_buff_cu_position_v;
	int16_t *pred_buff, *pred_buff_u, *pred_buff_v, *residual_buff, *residual_buff_u, *residual_buff_v;
	wnd_t *reference_wnd=NULL;//, *resi_wnd = NULL;
//	uint8_t *cbf_buff = NULL;
	int curr_part_size, curr_part_size_shift;
	int curr_part_size_chroma, curr_part_size_shift_chroma;
	int curr_part_x, curr_part_y, curr_part_global_x, curr_part_global_y;
	int curr_part_x_chroma, curr_part_y_chroma, curr_part_global_x_chroma, curr_part_global_y_chroma;
	int curr_depth = depth;
	cu_partition_info_t *parent_part_info;
	cu_partition_info_t *curr_cu_info;
	int curr_sum = 0, best_sum;
	int num_part_in_cu;
	int partition_cost;
	int cu_min_tu_size_shift;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int max_tr_depth, max_tr_processing_depth;
	int initial_state, end_state;
//	int cbf_split[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
//	int acc_cost[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int bitcost_cu_mode;
	int log2cu_size;
	int cu_x_position;
	motion_vector_t mv;
	int ref_idx = 0;
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
		mv = curr_cu_info->inter_mv[REF_PIC_LIST_0];
		ref_idx = curr_cu_info->inter_ref_index[REF_PIC_LIST_0];


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

		pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, Y_COMP);
		pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
		pred_buff_stride_chroma = WND_STRIDE_2D(et->prediction_wnd, CHR_COMP);
		pred_buff_u = WND_POSITION_2D(int16_t *, et->prediction_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
		pred_buff_v = WND_POSITION_2D(int16_t *, et->prediction_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

		orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
		orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
		orig_buff_stride_chroma = WND_STRIDE_2D(et->curr_mbs_wnd, CHR_COMP);
		orig_buff_u = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
		orig_buff_v = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

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

		get_mv_candidates(et, currslice, ctu, curr_cu_info, REF_PIC_LIST_0, ref_idx, part_size_type);//get candidates for motion search from the neigbour CUs
#ifdef COMPUTE_AS_HM
//		sad = hmr_motion_estimation_HM(et, ctu, curr_cu_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, curr_part_global_x, curr_part_global_y, 0, 0, curr_part_size, curr_part_size_shift, 64, 64, et->pict_width[Y_COMP], et->pict_height[Y_COMP], &mv);	
		mv_cost += select_mv_candidate(et, curr_cu_info, REF_PIC_LIST_0, &mv);
#else
//		sad = hmr_motion_estimation(et, ctu, curr_cu_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, curr_part_global_x, curr_part_global_y, 0, 0, curr_part_size, curr_part_size_shift, 64, 64, et->pict_width[Y_COMP], et->pict_height[Y_COMP], &mv);	
		mv_cost += select_mv_candidate(et, curr_cu_info, REF_PIC_LIST_0, &mv);
#endif
		SET_INTER_MV_BUFFS(et, ctu, curr_cu_info, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu);
		memset(&ctu->mv_ref_idx[REF_PIC_LIST_0][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_0], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));
		hmr_motion_compensation_luma(et, ctu, curr_cu_info, reference_buff_cu_position, reference_buff_stride, pred_buff, pred_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, &mv);
		et->funcs->predict(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride, residual_buff, residual_buff_stride, curr_part_size);

		hmr_motion_compensation_chroma(et, reference_buff_cu_position_u, reference_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv);
		hmr_motion_compensation_chroma(et, reference_buff_cu_position_v, reference_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv);	
		et->funcs->predict(orig_buff_u, orig_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, residual_buff_u, residual_buff_stride_chroma, curr_part_size_chroma);
		et->funcs->predict(orig_buff_v, orig_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, residual_buff_v, residual_buff_stride_chroma, curr_part_size_chroma);
		curr_cu_info++;
	}
	return mv_cost;
}



//this function is referenced by the initial depth, not by the processing depth
int encode_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type)
{
	int k;
	int cu_mode;
	uint32_t sad, cost, best_cost;
	slice_t *currslice = &et->ed->current_pict.slice;
	int curr_depth = depth;
	cu_partition_info_t *parent_part_info;
	cu_partition_info_t *curr_cu_info;
	int curr_sum_y = 0, curr_sum_u = 0, curr_sum_v = 0, best_sum;
	int num_part_in_cu;
	int partition_cost;
	int cu_min_tu_size_shift;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int max_tr_depth, max_tr_processing_depth;
	int initial_state, end_state;
//	int cbf_split[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
//	int acc_cost[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int bitcost_cu_mode;
	int log2cu_size;
	int cu_x_position;
	motion_vector_t mv;
	int ref_idx = 0;
	int processing_buff_depth;
	int tu_log_max_size = et->max_tu_size_shift;
	int qp, rem, per;

	best_cost = INT_MAX;
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
	if(et->performance_mode == 2)
		max_tr_processing_depth = depth==0?1:depth;//max_tr_processing_depth = (depth+2<=max_tr_processing_depth)?depth+2:((depth+1<=max_tr_processing_depth)?depth+1:max_tr_processing_depth);
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

//	if(et->performance_mode != 0 && curr_cu_info->recursive_split && curr_cu_info->children[0]!=NULL && curr_cu_info->children[0]->recursive_split==0)
//		max_tr_processing_depth = (depth+2<=max_tr_processing_depth)?depth+2:((depth+1<=max_tr_processing_depth)?depth+1:max_tr_processing_depth);


	memset(depth_state, 0, sizeof(depth_state));
	depth_state[curr_depth] = initial_state;

	curr_depth = curr_cu_info->depth;

//	processing_buff_depth = depth+(part_size_type==SIZE_NxN);

	if(et->ed->num_encoded_frames == 7 && ctu->ctu_number == 1 && curr_cu_info->abs_index == 64)
	{
		int iiii=0;
	}

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
				synchronize_motion_buffers_luma(et, curr_cu_info, &et->transform_quant_wnd[curr_depth+(part_size_type!=SIZE_2Nx2N)], &et->transform_quant_wnd[depth+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[depth+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
				synchronize_motion_buffers_chroma(et, curr_cu_info, &et->transform_quant_wnd[curr_depth+(part_size_type!=SIZE_2Nx2N)], &et->transform_quant_wnd[depth+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[depth+1+(part_size_type!=SIZE_2Nx2N)], gcnt);

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
					parent_part_info->cost = cost;
					parent_part_info->distortion = distortion;
					parent_part_info->sum = sum;

					if(parent_part_info->abs_index == 96)
					{
						int iiiii=0;
					}

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
					synchronize_motion_buffers_luma(et, parent_part_info, &et->transform_quant_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], &et->transform_quant_wnd[curr_depth-1+1+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[curr_depth-1+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
					synchronize_motion_buffers_chroma(et, parent_part_info, &et->transform_quant_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], &et->transform_quant_wnd[curr_depth-1+1+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[curr_depth+1+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[curr_depth-1+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
				}
				else 
				{
					SET_ENC_INFO_BUFFS(et, parent_part_info, depth+(part_size_type!=SIZE_2Nx2N), parent_part_info->abs_index, parent_part_info->num_part_in_cu);
//					synchronize_motion_buffers_luma(et, parent_part_info, &et->transform_quant_wnd[curr_depth-1+(part_size_type!=SIZE_2Nx2N)], &et->transform_quant_wnd[depth+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[curr_depth+1-1+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[depth+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
//					synchronize_motion_buffers_chroma(et, parent_part_info, &et->transform_quant_wnd[curr_depth-1+(part_size_type!=SIZE_2Nx2N)], &et->transform_quant_wnd[depth+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[curr_depth+1-1+(part_size_type!=SIZE_2Nx2N)], &et->decoded_mbs_wnd[depth+1+(part_size_type!=SIZE_2Nx2N)], gcnt);
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
	memcpy(&ctu->intra_mode[Y_COMP][abs_index], &et->intra_mode_buffs[Y_COMP][depth][abs_idx], num_part_in_cu*sizeof(et->intra_mode_buffs[0][0][0]));		\
	memcpy(&ctu->intra_mode[CHR_COMP][abs_index], &et->intra_mode_buffs[CHR_COMP][depth][abs_idx], num_part_in_cu*sizeof(et->intra_mode_buffs[0][0][0]));	\
}

#define SET_INTER_INFO_BUFFS(et, ctu, cu_info, abs_idx, num_partitions, ref_list)																			\
{																																							\
	int i;																																					\
	if(cu_info->prediction_mode == INTER_MODE)																												\
	{																																						\
		ctu->mv_diff_ref_idx[ref_list][abs_idx] = cu_info->best_candidate_idx[ref_list];																	\
		ctu->mv_diff[ref_list][abs_idx] = cu_info->best_dif_mv[ref_list];																					\
		memset(&ctu->inter_mode[abs_idx], 1<<ref_list, num_partitions*sizeof(ctu->inter_mode[0]));															\
		memset(&ctu->mv_ref_idx[ref_list][abs_index], cu_info->inter_ref_index[ref_list], num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));					\
		for(i=0;i<num_partitions;i++)																														\
		{																																					\
			ctu->mv_ref[ref_list][abs_idx+i] = cu_info->inter_mv[ref_list];																					\
		}																																					\
	}																																						\
	else																																					\
	{																																						\
		memset(&ctu->mv_ref_idx[ref_list][abs_index], -1, num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));													\
	}																																						\
	memset(&ctu->pred_mode[abs_index], cu_info->prediction_mode, num_part_in_cu*sizeof(ctu->pred_mode[0]));													\
}																																							


void consolidate_prediction_info(henc_thread_t *et, ctu_info_t *ctu, ctu_info_t *ctu_rd, cu_partition_info_t *parent_part_info, uint parent_cost, uint children_cost, int is_max_depth, uint *cost_sum)
{
	int abs_index = parent_part_info->abs_index;
	int num_part_in_cu = parent_part_info->num_part_in_cu;
	int curr_depth = parent_part_info->depth + 1;
	int gcnt = 0;
	int children_sum = parent_part_info->children[0]->sum+parent_part_info->children[1]->sum+parent_part_info->children[2]->sum+parent_part_info->children[3]->sum;

	//choose best
#ifdef COMPUTE_AS_HM
	if((children_cost<parent_cost || (/*et->ed->current_pict.slice.slice_type != I_SLICE && */children_cost==parent_cost && is_max_depth && parent_part_info->prediction_mode == INTRA_MODE && parent_part_info->children[0]->prediction_mode == INTER_MODE)) || !(parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame))//if we get here, tl should be inside the frame
#else
	if(children_cost<parent_cost || !(parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame))//if we get here, tl should be inside the frame
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
		parent_part_info->sum = parent_part_info->children[0]->sum+parent_part_info->children[1]->sum+parent_part_info->children[2]->sum+parent_part_info->children[3]->sum;

		if(is_max_depth)
		{
			int ll;
			int nchild;
			int num_part_in_sub_cu = parent_part_info->children[0]->num_part_in_cu;

			synchronize_motion_buffers_luma(et, parent_part_info, &et->transform_quant_wnd[curr_depth+1], &et->transform_quant_wnd[0], &et->decoded_mbs_wnd[curr_depth+1], &et->decoded_mbs_wnd[0], gcnt);
			synchronize_motion_buffers_chroma(et, parent_part_info, &et->transform_quant_wnd[curr_depth+1], &et->transform_quant_wnd[0], &et->decoded_mbs_wnd[curr_depth+1], &et->decoded_mbs_wnd[0], gcnt);

			CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu)

			for(nchild=0;nchild<4;nchild++)
			{
				cu_partition_info_t *cu_info = parent_part_info->children[nchild];
				SET_INTER_INFO_BUFFS(et, ctu, cu_info, cu_info->abs_index, cu_info->num_part_in_cu, REF_PIC_LIST_0)
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
		int part_size_type2 = (parent_part_info->depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//
		int parent_depth = parent_part_info->depth;//curr_depth-1

		synchronize_motion_buffers_luma(et, parent_part_info, &et->transform_quant_wnd[parent_depth+1], &et->transform_quant_wnd[0], &et->decoded_mbs_wnd[parent_depth+1], &et->decoded_mbs_wnd[0], gcnt);
		synchronize_motion_buffers_chroma(et, parent_part_info, &et->transform_quant_wnd[parent_depth+1], &et->transform_quant_wnd[0], &et->decoded_mbs_wnd[parent_depth+1], &et->decoded_mbs_wnd[0], gcnt);

		CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, parent_depth, abs_index, num_part_in_cu)
		SET_INTER_INFO_BUFFS(et, ctu, parent_part_info, abs_index, num_part_in_cu, REF_PIC_LIST_0)

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

	synchronize_motion_buffers_luma(et, curr_cu_info, &et->transform_quant_wnd[0], &et->transform_quant_wnd[curr_depth+1], &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[curr_depth+1], 0);
	synchronize_motion_buffers_chroma(et, curr_cu_info, &et->transform_quant_wnd[0], &et->transform_quant_wnd[curr_depth+1], &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[curr_depth+1], 0);
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

	synchronize_motion_buffers_luma(et, curr_cu_info, &et->transform_quant_wnd[curr_depth+1], &et->transform_quant_wnd[0], &et->decoded_mbs_wnd[curr_depth+1], &et->decoded_mbs_wnd[0], 0);
	synchronize_motion_buffers_chroma(et, curr_cu_info, &et->transform_quant_wnd[curr_depth+1], &et->transform_quant_wnd[0], &et->decoded_mbs_wnd[curr_depth+1], &et->decoded_mbs_wnd[0], 0);
}


uint motion_inter_(henc_thread_t* et, ctu_info_t* ctu, int gcnt)
{
	picture_t *currpict = &et->ed->current_pict;
	slice_t *currslice = &currpict->slice;
	double dist, dist_aux, best_cost;//, cost_luma, cost_chroma;
	int position = 0;
	int curr_depth = 0;
	ctu_info_t *ctu_rd = et->ctu_rd;
	cu_partition_info_t	*parent_part_info = NULL;
	cu_partition_info_t	*curr_cu_info = ctu->partition_list;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	uint cost_sum[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int abs_index;
	int num_part_in_cu;
	int ll, ithreads;
	double consumed_distortion = 0, avg_distortion = 0;
	uint consumed_ctus = 0;
	uint total_intra_partitions = 0, total_partitions;
//	int cbf_split[NUM_PICT_COMPONENTS] = {0,0,0};

#ifndef COMPUTE_AS_HM
	for(ithreads=0;ithreads<et->wfpp_num_threads;ithreads++)
	{
		henc_thread_t* henc_th = et->ed->thread[ithreads];
		
		consumed_distortion += henc_th->acc_dist;
		total_intra_partitions += henc_th->num_intra_partitions;
		consumed_ctus += henc_th->num_encoded_ctus;
	}
	total_partitions = consumed_ctus*et->num_partitions_in_cu;

	if(consumed_ctus>10 || consumed_ctus>et->ed->pict_total_ctu/15)
		avg_distortion = consumed_distortion/(consumed_ctus*ctu->num_part_in_ctu);		
	else
		avg_distortion = et->ed->avg_dist;

	if(et->index==0 && et->ed->num_encoded_frames >1 && et->ed->is_scene_change == 0 && consumed_ctus>et->ed->pict_total_ctu/15)
	{
		if(total_intra_partitions > (total_partitions/2.5))
		{
			et->ed->is_scene_change = 1;
			if(et->ed->gop_reinit_on_scene_change)
				et->ed->last_intra = currslice->poc;
			printf("\r\n---------------------scene change detected. total_intra_partitions:%d, total_partitions:%d , ed->avg_dist:%.2f, avg_distortion:%.2f, ----------------------\r\n", total_intra_partitions, total_partitions, et->ed->avg_dist, avg_distortion);
			hmr_rc_change_pic_mode(et, currslice);
//			int iiii=0;
		}
	}

//	avg_distortion = et->ed->avg_dist;
#endif

	//init rd auxiliar ctu
	if(et->rd_mode != RD_DIST_ONLY)
	{
		copy_ctu(ctu, ctu_rd);
	}

	if(ctu->y[Y_COMP] != 0)
	{
		int iiii=0;
	}


	if(et->ed->num_encoded_frames == 70 && ctu->ctu_number == 123)
	{
		int iiii=0;
	}

	while(curr_depth!=0 || depth_state[curr_depth]!=1)
	{
		double cost = 0, intra_cost = 0;
		int stop_recursion = FALSE;
		PartSize part_size_type = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;
		curr_depth = curr_cu_info->depth;
		num_part_in_cu = curr_cu_info->num_part_in_cu;
		abs_index = curr_cu_info->abs_index;
		
		position = curr_cu_info->list_index - et->partition_depth_start[curr_depth];

		//rc
//		if(currslice->slice_type != I_SLICE && curr_depth<=et->ed->qp_depth)
		{
			int orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
			uint8_t *orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, Y_COMP, curr_cu_info->x_position, curr_cu_info->y_position, 0, et->ctu_width);
			curr_cu_info->variance_luma = et->funcs->modified_variance(orig_buff, curr_cu_info->size, orig_buff_stride, 1)/(curr_cu_info->size*curr_cu_info->size);//for intra imgs this is done in analyse_intra_recursive_info
			orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, U_COMP);
			orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, U_COMP, curr_cu_info->x_position_chroma, curr_cu_info->y_position_chroma, gcnt, et->ctu_width);
			curr_cu_info->variance_chroma = 1.25*et->funcs->modified_variance(orig_buff, curr_cu_info->size_chroma, orig_buff_stride, 2)/(curr_cu_info->size_chroma*curr_cu_info->size_chroma);
			orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, V_COMP, curr_cu_info->x_position_chroma, curr_cu_info->y_position_chroma, gcnt, et->ctu_width);
			curr_cu_info->variance_chroma += 1.25*et->funcs->modified_variance(orig_buff, curr_cu_info->size_chroma, orig_buff_stride, 2)/(curr_cu_info->size_chroma*curr_cu_info->size_chroma);
			curr_cu_info->variance = curr_cu_info->variance_luma + curr_cu_info->variance_chroma;
		}

		curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);

#define DEPHT_ADD	40
#define DEPHT_MULT	1.1
		//if(ctu->ctu_number == 0 && abs_index==64)// && curr_depth==1)//ctu->ctu_number == 97 && et->ed->num_encoded_frames == 10 && && curr_depth==2  && abs_index == 64)
		if(curr_cu_info->is_b_inside_frame && curr_cu_info->is_r_inside_frame)//if br (and tl) are inside the frame, process
		{
			int mv_cost;

			if(part_size_type == SIZE_2Nx2N)
			{
				uint sad;
				//encode inter
				sad = hmr_cu_motion_estimation(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N, 	2.*curr_cu_info->size*curr_cu_info->size);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
#ifdef COMPUTE_AS_HM
				mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
				dist = encode_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
#else
				if(curr_cu_info->size < 64 || (sad<50000))// && curr_cu_info->size == 64))
				{
					mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					dist = encode_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
				}
				else
				{
//					printf("64x64 inter skipped");
					dist = MAX_COST;
					curr_cu_info->sum = MAX_COST;
					curr_cu_info->inter_mv[REF_PIC_LIST_0].hor_vector = curr_cu_info->inter_mv[REF_PIC_LIST_0].hor_vector = 0;
				}
#endif
				curr_cu_info->distortion = dist;
				cost = curr_cu_info->distortion;

#ifdef COMPUTE_AS_HM
				cost+=5*curr_depth;
#else
				cost += 2*mv_cost;
				cost=cost*DEPHT_MULT+DEPHT_ADD*curr_depth;
#endif
				curr_cu_info->cost = cost;
				curr_cu_info->prediction_mode = INTER_MODE;

//#ifndef COMPUTE_AS_HM
//				if((dist<.25*avg_distortion*curr_cu_info->num_part_in_cu ||/* (dist<1.5*avg_distortion*curr_cu_info->num_part_in_cu && curr_cu_info->sum <= curr_cu_info->size*.2) || */curr_cu_info->sum <= curr_cu_info->size*.05  ||
				//if((dist<.25*avg_distortion*curr_cu_info->num_part_in_cu || (dist<=1.5*avg_distortion*curr_cu_info->num_part_in_cu && curr_cu_info->sum <= (curr_cu_info->size*((float)MAX_QP-curr_cu_info->qp)*.025)) || curr_cu_info->sum <= (curr_cu_info->size*((float)MAX_QP-curr_cu_info->qp)*.005) ||
//				if((dist<.25*avg_distortion*curr_cu_info->num_part_in_cu || (dist<1.5*avg_distortion*curr_cu_info->num_part_in_cu && curr_cu_info->sum <= curr_cu_info->size*10./(double)curr_cu_info->qp) || curr_cu_info->sum <= curr_cu_info->size*2./(double)curr_cu_info->qp  ||
//					(curr_cu_info->parent!=NULL && curr_cu_info->parent->distortion!=MAX_COST && curr_cu_info->cost<curr_cu_info->parent->cost/8)) && (curr_depth+1)<et->max_pred_partition_depth && curr_cu_info->is_b_inside_frame && curr_cu_info->is_r_inside_frame)//stop recursion calls

				//if(((2*dist<curr_cu_info->variance) && (curr_cu_info->parent!=NULL && curr_cu_info->parent->distortion!=MAX_COST && curr_cu_info->cost<curr_cu_info->parent->cost/8)) || curr_cu_info->sum==0)
				if((curr_cu_info->sum==0 || (2*dist<curr_cu_info->variance)) && (curr_cu_info->parent!=NULL && curr_cu_info->parent->distortion!=MAX_COST && curr_cu_info->cost<curr_cu_info->parent->cost/6))
				{
					int max_processing_depth;
					consolidate_prediction_info(et, ctu, ctu_rd, curr_cu_info, curr_cu_info->cost, MAX_COST, FALSE, NULL);
					stop_recursion = TRUE;

					max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

					if(curr_depth <= max_processing_depth)//
					{
						int aux_depth;
						cu_partition_info_t* aux_partition_info = curr_cu_info;//(parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth])]:&ctu->partition_list[0];
						abs_index = aux_partition_info->abs_index;
						num_part_in_cu  = aux_partition_info->num_part_in_cu;

						for(aux_depth=curr_depth;aux_depth<=max_processing_depth;aux_depth++)
						{
							synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[aux_depth+1], gcnt);	
							//for rd
							if(et->rd_mode!=RD_DIST_ONLY)
								consolidate_info_buffers_for_rd(et, ctu, aux_depth, abs_index, num_part_in_cu);
						}
						synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
					}
				}

//				if(!stop_recursion && ((curr_cu_info->size<64 && dist>(1.5*avg_distortion*(curr_cu_info->num_part_in_cu+curr_depth)/*+4000*curr_cu_info->num_part_in_cu*/))))// || (curr_cu_info->size>16 && curr_cu_info->variance<curr_cu_info->size/4)))
				if(!stop_recursion && curr_cu_info->size<32 && curr_cu_info->variance<2*dist)// && dist>(1*avg_distortion*(curr_cu_info->num_part_in_cu)/*+4000*curr_cu_info->num_part_in_cu*/))))// || (curr_cu_info->size>16 && curr_cu_info->variance<curr_cu_info->size/4)))
//#endif
				{
					//encode intra
					uint inter_sum = curr_cu_info->sum;
					uint intra_dist;
					put_consolidated_info(et, ctu, curr_cu_info, curr_depth);
					intra_dist = encode_intra(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
#ifdef COMPUTE_AS_HM
					intra_cost = intra_dist+5*curr_depth;
					if(intra_cost < cost)
#else
					intra_cost = intra_dist*DEPHT_MULT+DEPHT_ADD*curr_depth;
					if(intra_cost+90*curr_cu_info->sum<cost+90*inter_sum)// && intra_cost<64*curr_cu_info->variance)
#endif
					{	//we prefer intra and it is already in its buffer
						curr_cu_info->cost = intra_cost;
						curr_cu_info->distortion = intra_dist;
						curr_cu_info->sum = curr_cu_info->sum;
						curr_cu_info->prediction_mode = INTRA_MODE;
					}
					else
					{	//we prefer inter, bring it back
						get_back_consolidated_info(et, ctu, curr_cu_info, curr_depth);
						curr_cu_info->cost = cost;
						curr_cu_info->distortion = dist;
						curr_cu_info->sum = inter_sum;
						curr_cu_info->prediction_mode = INTER_MODE;
					}
				}
//#ifndef COMPUTE_AS_HM

				dist = curr_cu_info->distortion;

				if((curr_cu_info->sum==0 || (2*dist<curr_cu_info->variance)) && (curr_cu_info->parent!=NULL && curr_cu_info->parent->distortion!=MAX_COST && curr_cu_info->cost<curr_cu_info->parent->cost/6))
				{
					int max_processing_depth;
					consolidate_prediction_info(et, ctu, ctu_rd, curr_cu_info, curr_cu_info->cost, MAX_COST, FALSE, NULL);
					stop_recursion = TRUE;

					max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

					if(curr_depth <= max_processing_depth)//el = es para cuando et->max_intra_tr_depth!=4
					{
						int aux_depth;
						cu_partition_info_t* aux_partition_info = curr_cu_info;//(parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth])]:&ctu->partition_list[0];
						abs_index = aux_partition_info->abs_index;
						num_part_in_cu  = aux_partition_info->num_part_in_cu;

						for(aux_depth=curr_depth;aux_depth<=max_processing_depth;aux_depth++)
						{
							synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[aux_depth+1], gcnt);	
							//for rd
							if(et->rd_mode!=RD_DIST_ONLY)
								consolidate_info_buffers_for_rd(et, ctu, aux_depth, abs_index, num_part_in_cu);
						}
						synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
					}
				}

				//if this matches, it is useless to continue the recursion. the case where curr_depth!=et->max_pred_partition_depth is checked at the end of the consolidation loop)
/*				if(curr_depth>0 && depth_state[curr_depth]!=3 && curr_depth == et->max_pred_partition_depth && cost_sum[curr_depth]+curr_cu_info->cost>parent_part_info->cost && parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
				{
					depth_state[curr_depth] = 3;
				}
*/
//#endif
			}
			else if(part_size_type == SIZE_NxN)//intra NxN is processed in its current depth, while inter NxN is processed in its father´s depth. So, intra NxN does not have to be compaired
			{
				if(curr_cu_info->abs_index==48)
				{
					int iiii=0;
				}

				curr_cu_info->cost = curr_cu_info->distortion = cost = dist = MAX_COST;
				if((curr_depth-1) == (et->max_cu_depth - et->mincu_mintr_shift_diff) && curr_cu_info->parent->size>8)	//SIZE_NxN
				{
					uint sad;
					int position_aux = curr_cu_info->parent->list_index - et->partition_depth_start[curr_depth-1];
					uint aux_cost = curr_cu_info->parent->cost;
					uint aux_dist = curr_cu_info->parent->distortion;
					uint aux_sum = curr_cu_info->parent->sum;

//					mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_NxN, 0);//.25*avg_distortion*4*curr_cu_info->num_part_in_cu);
					sad = hmr_cu_motion_estimation(et, ctu, gcnt, curr_depth, position, SIZE_NxN, 0);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_NxN);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					dist = encode_inter(et, ctu, gcnt, curr_depth-1, position_aux, SIZE_NxN);//this function is referenced by the initial depth, not by the processing depth
#ifdef COMPUTE_AS_HM
					cost=dist+5*curr_depth;
#else
					cost=dist*DEPHT_MULT+DEPHT_ADD*curr_depth;
					cost += mv_cost;
#endif
					curr_cu_info[0].cost = curr_cu_info[1].cost = curr_cu_info[2].cost = curr_cu_info[3].cost = 0;
					curr_cu_info[0].sum = curr_cu_info[1].sum = curr_cu_info[2].sum = curr_cu_info[3].sum = 0;
					curr_cu_info[0].distortion = curr_cu_info[1].distortion = curr_cu_info[2].distortion = curr_cu_info[3].distortion = 0;

					curr_cu_info->cost = cost;
					curr_cu_info->distortion = dist;
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
						intra_cost = intra_dist*DEPHT_MULT+DEPHT_ADD*curr_depth;
						intra_sum = curr_cu_info[0].sum + curr_cu_info[1].sum + curr_cu_info[2].sum + curr_cu_info[3].sum;

						curr_cu_info[0].cost = curr_cu_info[1].cost = curr_cu_info[2].cost = curr_cu_info[3].cost = 0;
						curr_cu_info[0].sum = curr_cu_info[1].sum = curr_cu_info[2].sum = curr_cu_info[3].sum = 0;
						curr_cu_info[0].distortion = curr_cu_info[1].distortion = curr_cu_info[2].distortion = curr_cu_info[3].distortion = 0;
#ifdef COMPUTE_AS_HM
						intra_cost=intra_dist+5*curr_depth;
						if(intra_cost < cost)
#else
						if(1.25*intra_cost+45*intra_sum<cost+45*previous_sum)// && intra_cost<64*curr_cu_info->variance)
#endif
						{	//we prefer intra and it is already in its buffer
							curr_cu_info->cost = intra_cost;
							curr_cu_info->sum = intra_sum;
							curr_cu_info->distortion = intra_dist;
							curr_cu_info[0].prediction_mode = curr_cu_info[1].prediction_mode = curr_cu_info[2].prediction_mode = curr_cu_info[3].prediction_mode = INTRA_MODE;
						}
						else
						{	//we prefer inter, bring it back
							get_back_consolidated_info(et, ctu, curr_cu_info->parent, curr_depth);
							curr_cu_info->cost = cost;
							curr_cu_info->sum = previous_sum;
							curr_cu_info->distortion = dist;
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
		else if(depth_state[curr_depth]==4)//la depth =1 lo hemos consolidado antes del bucle
		{
			int max_processing_depth;

			while(depth_state[curr_depth]==4 && curr_depth>0)//>0 pq consolidamos sobre el padre, 
			{
				int is_max_depth = (curr_depth==et->max_pred_partition_depth);
				cost = parent_part_info->children[0]->cost + parent_part_info->children[1]->cost +parent_part_info->children[2]->cost+parent_part_info->children[3]->cost;

				depth_state[curr_depth] = 0;

				best_cost = parent_part_info->cost;

				if(curr_depth == 1 && et->ed->num_encoded_frames == 25 && ctu->ctu_number == 123)
				{
					int iiiiii=0;
				}

				consolidate_prediction_info(et, ctu, ctu_rd, parent_part_info, best_cost, cost, is_max_depth, cost_sum);

				depth_state[curr_depth] = 0;
				cost_sum[curr_depth] = 0;

				curr_depth--;
				parent_part_info = parent_part_info->parent;

#ifndef COMPUTE_AS_HM
				if(curr_depth>0 && cost_sum[curr_depth] > parent_part_info->cost && depth_state[curr_depth]<4 && ctu->partition_list[0].is_b_inside_frame && ctu->partition_list[0].is_r_inside_frame)
				{
					depth_state[curr_depth] = 4;
				}
#endif
			}

			max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

			if(curr_depth <= max_processing_depth)//el = es para cuando et->max_intra_tr_depth!=4
			{
				int aux_depth;
				cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
				abs_index = aux_partition_info->abs_index;
				num_part_in_cu  = aux_partition_info->num_part_in_cu;

				for(aux_depth=curr_depth;aux_depth<=max_processing_depth;aux_depth++)
				{
					synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[aux_depth+1], gcnt);	
					//for rd
//					if(et->rd_mode!=RD_DIST_ONLY)
//						CONSOLIDATE_INTRA_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu);
				}
				synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
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
		CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu)
		SET_INTER_INFO_BUFFS(et, ctu, curr_cu_info, abs_index, num_part_in_cu, REF_PIC_LIST_0)	
	}


//	memset(&ctu->pred_mode[abs_index], INTER_MODE, num_part_in_cu*sizeof(ctu->pred_mode[0]));//signal all partitions as inter
	memset(&ctu->skipped[abs_index], FALSE, num_part_in_cu*sizeof(ctu->skipped[0]));//signal all partitions as non skipped
	return curr_cu_info->cost;
}


uint motion_inter__(henc_thread_t* et, ctu_info_t* ctu, int gcnt)
{
	picture_t *currpict = &et->ed->current_pict;
	slice_t *currslice = &currpict->slice;
	double dist, dist_aux, best_cost;//, cost_luma, cost_chroma;
	int position = 0;
	int curr_depth = 0;
	ctu_info_t *ctu_rd = et->ctu_rd;
	cu_partition_info_t	*parent_part_info = NULL;
	cu_partition_info_t	*curr_cu_info = ctu->partition_list;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	uint cost_sum[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int abs_index;
	int num_part_in_cu;
	int ll, ithreads;
	double consumed_distortion = 0, avg_distortion = 0;
	uint consumed_ctus = 0;
	uint total_intra_partitions = 0, total_partitions;
//	int cbf_split[NUM_PICT_COMPONENTS] = {0,0,0};

#ifndef COMPUTE_AS_HM
	for(ithreads=0;ithreads<et->wfpp_num_threads;ithreads++)
	{
		henc_thread_t* henc_th = et->ed->thread[ithreads];
		
		consumed_distortion += henc_th->acc_dist;
		total_intra_partitions += henc_th->num_intra_partitions;
		consumed_ctus += henc_th->num_encoded_ctus;
	}
	total_partitions = consumed_ctus*et->num_partitions_in_cu;

	if(consumed_ctus>10 || consumed_ctus>et->ed->pict_total_ctu/15)
		avg_distortion = consumed_distortion/(consumed_ctus*ctu->num_part_in_ctu);		
	else
		avg_distortion = et->ed->avg_dist;

	if(et->index==0 && et->ed->num_encoded_frames >1 && et->ed->is_scene_change == 0 && consumed_ctus>et->ed->pict_total_ctu/15)
	{
		if(total_intra_partitions > (total_partitions/2.5))
		{
			et->ed->is_scene_change = 1;
			if(et->ed->gop_reinit_on_scene_change)
				et->ed->last_intra = currslice->poc;
			printf("\r\n---------------------scene change detected. total_intra_partitions:%d, total_partitions:%d , ed->avg_dist:%.2f, avg_distortion:%.2f, ----------------------\r\n", total_intra_partitions, total_partitions, et->ed->avg_dist, avg_distortion);
			hmr_rc_change_pic_mode(et, currslice);
//			int iiii=0;
		}
	}

//	avg_distortion = et->ed->avg_dist;
#endif

	//init rd auxiliar ctu
	if(et->rd_mode == RD_DIST_ONLY)
	{
		copy_ctu(ctu, ctu_rd);
	}

	if(ctu->y[Y_COMP] != 0)
	{
		int iiii=0;
	}


	curr_depth = 0;
	position = 0;
	if(curr_cu_info->is_b_inside_frame && curr_cu_info->is_r_inside_frame)
		curr_cu_info->sad = hmr_cu_motion_estimation(et, ctu, gcnt, 0, 0, SIZE_2Nx2N, 2.*curr_cu_info->size*curr_cu_info->size);
	else
		curr_cu_info->sad = MAX_COST;


	while(curr_depth!=0 || depth_state[curr_depth]!=1)
	{
		uint child_sad, child_sum, curr_sad, curr_sum;
		uint child_depth;
		uint child_position;
		double cost = 0, intra_cost = 0;
		int stop_recursion = FALSE;
		PartSize part_size_type = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;
		curr_depth = curr_cu_info->depth;
		child_depth = curr_cu_info->children[0]->depth;

		num_part_in_cu = curr_cu_info->num_part_in_cu;
		abs_index = curr_cu_info->abs_index;
		
		position = curr_cu_info->list_index - et->partition_depth_start[curr_depth];
		child_position = curr_cu_info->children[0]->list_index - et->partition_depth_start[child_depth];

		//rc
		if(currslice->slice_type != I_SLICE && curr_depth<=et->ed->qp_depth)
		{
			int orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
			uint8_t *orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, Y_COMP, curr_cu_info->x_position, curr_cu_info->y_position, 0, et->ctu_width);
			curr_cu_info->variance_luma = et->funcs->modified_variance(orig_buff, curr_cu_info->size, orig_buff_stride, 1)/(curr_cu_info->size*curr_cu_info->size);//for intra imgs this is done in analyse_intra_recursive_info
			orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, U_COMP);
			orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, U_COMP, curr_cu_info->x_position_chroma, curr_cu_info->y_position_chroma, gcnt, et->ctu_width);
			curr_cu_info->variance_chroma = 1.25*et->funcs->modified_variance(orig_buff, curr_cu_info->size_chroma, orig_buff_stride, 2)/(curr_cu_info->size_chroma*curr_cu_info->size_chroma);
			orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, V_COMP, curr_cu_info->x_position_chroma, curr_cu_info->y_position_chroma, gcnt, et->ctu_width);
			curr_cu_info->variance_chroma += 1.25*et->funcs->modified_variance(orig_buff, curr_cu_info->size_chroma, orig_buff_stride, 2)/(curr_cu_info->size_chroma*curr_cu_info->size_chroma);
			curr_cu_info->variance = curr_cu_info->variance_luma + curr_cu_info->variance_chroma;
		}

		curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);

#define DEPHT_ADD	40
#define DEPHT_MULT	1.1
		//if(ctu->ctu_number == 0 && abs_index==64)// && curr_depth==1)//ctu->ctu_number == 97 && et->ed->num_encoded_frames == 10 && && curr_depth==2  && abs_index == 64)
		if(curr_cu_info->is_tl_inside_frame)//is_b_inside_frame && curr_cu_info->is_r_inside_frame)//if br (and tl) are inside the frame, process
		{
			int mv_cost;

			if(part_size_type == SIZE_2Nx2N)
			{
				uint child_size = curr_cu_info->children[0]->size;

				if(curr_depth!=3)
				{
					int i;
					for(i=0;i<4;i++)
					{
						curr_cu_info->children[i]->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info->children[i], currslice);
						if(curr_cu_info->children[i]->is_b_inside_frame && curr_cu_info->children[i]->is_r_inside_frame)
						{
							curr_cu_info->children[i]->sad = hmr_cu_motion_estimation(et, ctu, gcnt, child_depth, child_position+i, SIZE_2Nx2N, 2.*curr_cu_info->size*curr_cu_info->size);//.25*avg_distortion*curr_cu_info->num_part_in_cu);		
							SET_INTER_MV_BUFFS(et, ctu, curr_cu_info, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu)				
							memset(&ctu->mv_ref_idx[REF_PIC_LIST_0][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_0], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));
						}
						else
						{
							if(curr_cu_info->children[i]->is_tl_inside_frame)
								curr_cu_info->children[i]->sad = MAX_COST;
							else
								curr_cu_info->children[i]->sad = 0;
						}
					}
				}
				else
					depth_state[curr_depth]=3;
			}
		}

		depth_state[curr_depth]++;

		if((curr_depth)<et->max_pred_partition_depth-1 && curr_cu_info->is_tl_inside_frame && !stop_recursion)//depth_state[curr_depth]!=4 is for fast skip//if tl is not inside the frame don't process the next depths
		{
			curr_depth++;
			parent_part_info = curr_cu_info;
		}
		else if(depth_state[curr_depth]==4)//la depth =1 lo hemos consolidado antes del bucle
		{
			int max_processing_depth;

			while(depth_state[curr_depth]==4 && curr_depth>0)//>0 pq consolidamos sobre el padre, 
			{
				uint abs_index = parent_part_info->abs_index;
				uint parent_sad = parent_part_info->sad;
				uint child_sad = parent_part_info->children[0]->sad+parent_part_info->children[1]->sad+parent_part_info->children[2]->sad+parent_part_info->children[3]->sad;

				if(parent_sad<child_sad+parent_part_info->size && parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
				{
					memset(&ctu->pred_depth[abs_index], parent_part_info->depth, parent_part_info->num_part_in_cu*sizeof(ctu->pred_depth[0]));
				}
				else
				{
					if(curr_depth == et->max_pred_partition_depth-1)// ||  parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
						memset(&ctu->pred_depth[abs_index], parent_part_info->children[0]->depth, parent_part_info->num_part_in_cu*sizeof(ctu->pred_depth[0]));
					parent_part_info->sad = child_sad;
				}
				depth_state[curr_depth] = 0;
				cost_sum[curr_depth] = 0;

				if(ctu->ctu_number >= 220)// && curr_cu_info->abs_index==128)
				{
					int iiiii=0;
				}

				curr_depth--;
				parent_part_info = parent_part_info->parent;
			}

		}

		if(parent_part_info!=NULL)
			curr_cu_info = parent_part_info->children[depth_state[curr_depth]];
	}
	

	parent_part_info = curr_cu_info = &ctu->partition_list[0];
	curr_depth = 0;
	memset(depth_state,0,sizeof(depth_state));

	while(curr_depth!=0 || depth_state[curr_depth]!=1)
	{
		uint child_sad, child_sum, curr_sad, curr_sum;
		uint child_depth;
		uint child_position;
		double cost = 0, intra_cost = 0;
		int stop_recursion = FALSE;
		PartSize part_size_type = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;
		int pred_depth = ctu->pred_depth[curr_cu_info->abs_index];

		curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);

		if(curr_cu_info->is_tl_inside_frame)
		{
			while(curr_depth<pred_depth)
			{
				depth_state[curr_depth]++;
				curr_depth++;
				curr_cu_info = curr_cu_info->children[depth_state[curr_depth]];
				curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);
			}
			curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);
			curr_depth = curr_cu_info->depth;
			parent_part_info = curr_cu_info->parent;

			if(curr_depth!=pred_depth)
			{
				int iiiii=0;
			}

			if(ctu->ctu_number==48)
			{
				int iiiii=0;
			}

			num_part_in_cu = curr_cu_info->num_part_in_cu;
			abs_index = curr_cu_info->abs_index;
		
			position = curr_cu_info->list_index - et->partition_depth_start[curr_depth];


#define DEPHT_ADD	40
#define DEPHT_MULT	1.1
			//if(ctu->ctu_number == 0 && abs_index==64)// && curr_depth==1)//ctu->ctu_number == 97 && et->ed->num_encoded_frames == 10 && && curr_depth==2  && abs_index == 64)
			if(curr_cu_info->is_b_inside_frame && curr_cu_info->is_r_inside_frame)//if br (and tl) are inside the frame, process
			{
				int mv_cost;
				int is_max_depth = (curr_depth==et->max_pred_partition_depth-1);
	//			if(part_size_type == SIZE_2Nx2N)
				{
	//				uint child_size = curr_cu_info->children[0]->size;
//					curr_cu_info->sad = hmr_cu_motion_estimation(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N, 2.*curr_cu_info->size*curr_cu_info->size);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					curr_cu_info->distortion = dist = encode_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);

					cost = curr_cu_info->distortion + 2*mv_cost;
					cost=cost*DEPHT_MULT+DEPHT_ADD*curr_depth;
					curr_cu_info->cost = cost;
					curr_cu_info->prediction_mode = INTER_MODE;
					{
						//encode intra
						uint inter_sum = curr_cu_info->sum;
						uint intra_dist;
						put_consolidated_info(et, ctu, curr_cu_info, curr_depth);
						intra_dist = encode_intra(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
#ifdef COMPUTE_AS_HM
						intra_cost = intra_dist+5*curr_depth;
						if(intra_cost < cost)
#else
						intra_cost = intra_dist*DEPHT_MULT+DEPHT_ADD*curr_depth;
						if(intra_cost+90*curr_cu_info->sum<cost+90*inter_sum)// && intra_cost<64*curr_cu_info->variance)
#endif
						{	//we prefer intra and it is already in its buffer
							curr_cu_info->cost = intra_cost;
							curr_cu_info->distortion = intra_dist;
							curr_cu_info->sum = curr_cu_info->sum;
							curr_cu_info->prediction_mode = INTRA_MODE;
						}
						else
						{	//we prefer inter, bring it back
							get_back_consolidated_info(et, ctu, curr_cu_info, curr_depth);
							curr_cu_info->cost = cost;
							curr_cu_info->distortion = dist;
							curr_cu_info->sum = inter_sum;
							curr_cu_info->prediction_mode = INTER_MODE;
						}
						consolidate_prediction_info(et, ctu, ctu_rd, curr_cu_info, curr_cu_info->cost, MAX_COST, FALSE, NULL);
					}
				}
			}
			else
			{
				curr_cu_info->cost = MAX_COST;
			}
		}

		depth_state[curr_depth]++;

		if(depth_state[curr_depth]==4)//la depth =1 lo hemos consolidado antes del bucle
		{
			int max_processing_depth;

			while(depth_state[curr_depth]==4 && curr_depth>0)//>0 pq consolidamos sobre el padre, 
			{

				depth_state[curr_depth] = 0;
				cost_sum[curr_depth] = 0;

				curr_depth--;
				curr_cu_info = curr_cu_info->parent;
			}
		}
		curr_cu_info++;// = parent_part_info->children[depth_state[curr_depth]];
	}


	curr_cu_info = &ctu->partition_list[0];
	abs_index = curr_cu_info->abs_index;
	curr_depth = curr_cu_info->depth;
	num_part_in_cu = curr_cu_info->num_part_in_cu;

	//if pred_depth==0 there is no NxN subdivision. we need to collect the information of the ctu
	if(et->max_pred_partition_depth==0)
	{
		CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu)
		SET_INTER_INFO_BUFFS(et, ctu, curr_cu_info, abs_index, num_part_in_cu, REF_PIC_LIST_0)	
	}


//	memset(&ctu->pred_mode[abs_index], INTER_MODE, num_part_in_cu*sizeof(ctu->pred_mode[0]));//signal all partitions as inter
	memset(&ctu->skipped[abs_index], FALSE, num_part_in_cu*sizeof(ctu->skipped[0]));//signal all partitions as non skipped
	return curr_cu_info->cost;
}



uint motion_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt)
{
	picture_t *currpict = &et->ed->current_pict;
	slice_t *currslice = &currpict->slice;
	double dist, dist_aux, best_cost;//, cost_luma, cost_chroma;
	int position = 0;
	int curr_depth = 0;
	ctu_info_t *ctu_rd = et->ctu_rd;
	cu_partition_info_t	*parent_part_info = NULL;
	cu_partition_info_t	*curr_cu_info = ctu->partition_list;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	uint cost_sum[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int abs_index;
	int num_part_in_cu;
	int ll, ithreads;
	double consumed_distortion = 0, avg_distortion = 0;
	uint consumed_ctus = 0;
	uint total_intra_partitions = 0, total_partitions;
//	int cbf_split[NUM_PICT_COMPONENTS] = {0,0,0};

#ifndef COMPUTE_AS_HM
	for(ithreads=0;ithreads<et->wfpp_num_threads;ithreads++)
	{
		henc_thread_t* henc_th = et->ed->thread[ithreads];
		
		consumed_distortion += henc_th->acc_dist;
		total_intra_partitions += henc_th->num_intra_partitions;
		consumed_ctus += henc_th->num_encoded_ctus;
	}
	total_partitions = consumed_ctus*et->num_partitions_in_cu;

	if(consumed_ctus>10 || consumed_ctus>et->ed->pict_total_ctu/15)
	{
		avg_distortion = consumed_distortion/(consumed_ctus*ctu->num_part_in_ctu);		
		et->ed->avg_dist = avg_distortion;//update avg_dist as it evolves
	}
	else
		avg_distortion = et->ed->avg_dist;

	if(et->index==0 && et->ed->num_encoded_frames >1 && et->ed->is_scene_change == 0 && consumed_ctus>et->ed->pict_total_ctu/16)
	{
		if(total_intra_partitions > (total_partitions*.8))
		{
			et->ed->is_scene_change = 1;
			if(et->ed->gop_reinit_on_scene_change)
				et->ed->last_intra = currslice->poc;
			printf("\r\n---------------------scene change detected. total_intra_partitions:%d, total_partitions:%d , ed->avg_dist:%.2f, avg_distortion:%.2f, ----------------------\r\n", total_intra_partitions, total_partitions, et->ed->avg_dist, avg_distortion);
			hmr_rc_change_pic_mode(et, currslice);
//			int iiii=0;
		}
	}

//	avg_distortion = et->ed->avg_dist;
#endif

	//init rd auxiliar ctu
	if(et->rd_mode == RD_DIST_ONLY)
	{
		copy_ctu(ctu, ctu_rd);
	}

	if(ctu->y[Y_COMP] != 0)
	{
		int iiii=0;
	}

	curr_depth = 0;
	position = 0;
	if(curr_cu_info->is_b_inside_frame && curr_cu_info->is_r_inside_frame)
		curr_cu_info->sad = hmr_cu_motion_estimation(et, ctu, gcnt, 0, 0, SIZE_2Nx2N, 2.*curr_cu_info->size*curr_cu_info->size);
	else
		curr_cu_info->sad = MAX_COST;


	while(curr_depth!=0 || depth_state[curr_depth]!=1)
	{
		uint child_sad, child_sum, curr_sad, curr_sum;
		uint child_depth;
		uint child_position;
		double cost = 0, intra_cost = 0;
		int stop_recursion = FALSE;
		PartSize part_size_type = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;
		curr_depth = curr_cu_info->depth;
		child_depth = curr_cu_info->children[0]->depth;

		num_part_in_cu = curr_cu_info->num_part_in_cu;
		abs_index = curr_cu_info->abs_index;
		
		position = curr_cu_info->list_index - et->partition_depth_start[curr_depth];
		child_position = curr_cu_info->children[0]->list_index - et->partition_depth_start[child_depth];

		//rc
		if(currslice->slice_type != I_SLICE && curr_depth<=et->ed->qp_depth)
		{
			int orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
			uint8_t *orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, Y_COMP, curr_cu_info->x_position, curr_cu_info->y_position, 0, et->ctu_width);
			curr_cu_info->variance_luma = et->funcs->modified_variance(orig_buff, curr_cu_info->size, orig_buff_stride, 1)/(curr_cu_info->size*curr_cu_info->size);//for intra imgs this is done in analyse_intra_recursive_info
			orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, U_COMP);
			orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, U_COMP, curr_cu_info->x_position_chroma, curr_cu_info->y_position_chroma, gcnt, et->ctu_width);
			curr_cu_info->variance_chroma = 1.25*et->funcs->modified_variance(orig_buff, curr_cu_info->size_chroma, orig_buff_stride, 2)/(curr_cu_info->size_chroma*curr_cu_info->size_chroma);
			orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, V_COMP, curr_cu_info->x_position_chroma, curr_cu_info->y_position_chroma, gcnt, et->ctu_width);
			curr_cu_info->variance_chroma += 1.25*et->funcs->modified_variance(orig_buff, curr_cu_info->size_chroma, orig_buff_stride, 2)/(curr_cu_info->size_chroma*curr_cu_info->size_chroma);
			curr_cu_info->variance = curr_cu_info->variance_luma + curr_cu_info->variance_chroma;
		}

		curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);

#define DEPHT_ADD	40
#define DEPHT_MULT	1.1
		//if(ctu->ctu_number == 0 && abs_index==64)// && curr_depth==1)//ctu->ctu_number == 97 && et->ed->num_encoded_frames == 10 && && curr_depth==2  && abs_index == 64)
		if(curr_cu_info->is_tl_inside_frame)//is_b_inside_frame && curr_cu_info->is_r_inside_frame)//if br (and tl) are inside the frame, process
		{
			int mv_cost;

			if(part_size_type == SIZE_2Nx2N)
			{
				uint child_size = curr_cu_info->children[0]->size;

				if(curr_depth!=(et->max_pred_partition_depth-1))
				{
					int i;
					for(i=0;i<4;i++)
					{
						curr_cu_info->children[i]->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info->children[i], currslice);
						if(curr_cu_info->children[i]->is_b_inside_frame && curr_cu_info->children[i]->is_r_inside_frame)
						{
							curr_cu_info->children[i]->sad = hmr_cu_motion_estimation(et, ctu, gcnt, child_depth, child_position+i, SIZE_2Nx2N, 2.*curr_cu_info->size*curr_cu_info->size);//.25*avg_distortion*curr_cu_info->num_part_in_cu);		
//							SET_INTER_MV_BUFFS(et, ctu, curr_cu_info, curr_cu_info->abs_index, curr_cu_info->num_part_in_cu)				
//							memset(&ctu->mv_ref_idx[REF_PIC_LIST_0][curr_cu_info->abs_index], curr_cu_info->inter_ref_index[REF_PIC_LIST_0], curr_cu_info->num_part_in_cu*sizeof(ctu->mv_ref_idx[0][0]));
						}
						else
						{
							if(curr_cu_info->children[i]->is_tl_inside_frame)
								curr_cu_info->children[i]->sad = MAX_COST;
							else
								curr_cu_info->children[i]->sad = 0;
						}
					}
				}
				else
					depth_state[curr_depth]=3;
			}
		}

		depth_state[curr_depth]++;

		if((curr_depth)<et->max_pred_partition_depth-1 && curr_cu_info->is_tl_inside_frame && !stop_recursion)//depth_state[curr_depth]!=4 is for fast skip//if tl is not inside the frame don't process the next depths
		{
			curr_depth++;
			parent_part_info = curr_cu_info;
		}
		else if(depth_state[curr_depth]==4)//la depth =1 lo hemos consolidado antes del bucle
		{
			int max_processing_depth;

			while(depth_state[curr_depth]==4 && curr_depth>0)//>0 pq consolidamos sobre el padre, 
			{
				uint abs_index = parent_part_info->abs_index;
				uint parent_sad = parent_part_info->sad;
				uint child_sad = parent_part_info->children[0]->sad+parent_part_info->children[1]->sad+parent_part_info->children[2]->sad+parent_part_info->children[3]->sad;

				if(parent_sad<child_sad+parent_part_info->size_chroma && parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
				{
					memset(&ctu->pred_depth[abs_index], parent_part_info->depth, parent_part_info->num_part_in_cu*sizeof(ctu->pred_depth[0]));
				}
				else
				{
					if(curr_depth == et->max_pred_partition_depth-1)// ||  parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
						memset(&ctu->pred_depth[abs_index], parent_part_info->children[0]->depth, parent_part_info->num_part_in_cu*sizeof(ctu->pred_depth[0]));
					parent_part_info->sad = child_sad;
				}
				depth_state[curr_depth] = 0;
				cost_sum[curr_depth] = 0;

				curr_depth--;
				parent_part_info = parent_part_info->parent;
			}

		}

		if(parent_part_info!=NULL)
			curr_cu_info = parent_part_info->children[depth_state[curr_depth]];
	}
	

	parent_part_info = curr_cu_info = &ctu->partition_list[0];
	curr_depth = 0;
	memset(depth_state,0,sizeof(depth_state));

	while(curr_depth!=0 || depth_state[curr_depth]!=1)
	{
		uint child_sad, child_sum, curr_sad, curr_sum;
		uint child_depth;
		uint child_position;
		double cost = 0, intra_cost = 0;
		int stop_recursion = FALSE;
		PartSize part_size_type = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;
		int pred_depth = ctu->pred_depth[curr_cu_info->abs_index];

		curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);

		if(curr_cu_info->is_tl_inside_frame)
		{
			while(curr_depth<pred_depth)
			{
				depth_state[curr_depth]++;
				curr_depth++;
				curr_cu_info = curr_cu_info->children[depth_state[curr_depth]];
				curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);
			}
			curr_cu_info->qp = hmr_rc_get_cu_qp(et, ctu, curr_cu_info, currslice);
			curr_depth = curr_cu_info->depth;
			parent_part_info = curr_cu_info->parent;

			if(curr_depth!=pred_depth)
			{
				int iiiii=0;
			}

			num_part_in_cu = curr_cu_info->num_part_in_cu;
			abs_index = curr_cu_info->abs_index;
		
			position = curr_cu_info->list_index - et->partition_depth_start[curr_depth];


#define DEPHT_ADD	40
#define DEPHT_MULT	1.1
			//if(ctu->ctu_number == 0 && abs_index==64)// && curr_depth==1)//ctu->ctu_number == 97 && et->ed->num_encoded_frames == 10 && && curr_depth==2  && abs_index == 64)
			if(curr_cu_info->is_b_inside_frame && curr_cu_info->is_r_inside_frame)//if br (and tl) are inside the frame, process
			{
				int mv_cost;
				int is_max_depth = (curr_depth==et->max_pred_partition_depth-1);
	//			if(part_size_type == SIZE_2Nx2N)
				{
	//				uint child_size = curr_cu_info->children[0]->size;
//					curr_cu_info->sad = hmr_cu_motion_estimation(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N, 2.*curr_cu_info->size*curr_cu_info->size);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					mv_cost = predict_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);//.25*avg_distortion*curr_cu_info->num_part_in_cu);
					curr_cu_info->distortion = dist = encode_inter(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);

					if(ctu->ctu_number == 194 && curr_cu_info->abs_index==48)
					{
						int iiii=0;
					}


					cost = curr_cu_info->distortion + 2*mv_cost;
					cost=cost*DEPHT_MULT+DEPHT_ADD*curr_depth;
					curr_cu_info->cost = cost;
					curr_cu_info->prediction_mode = INTER_MODE;
					consolidate_prediction_info(et, ctu, ctu_rd, curr_cu_info, curr_cu_info->cost, MAX_COST, FALSE, NULL);
					if(curr_cu_info->size < 64)
					{
						//encode intra
						uint inter_sum = curr_cu_info->sum;
						uint intra_dist;
						int max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);
//						put_consolidated_info(et, ctu, curr_cu_info, curr_depth);

						intra_dist = encode_intra(et, ctu, gcnt, curr_depth, position, SIZE_2Nx2N);
#ifdef COMPUTE_AS_HM
						intra_cost = intra_dist+5*curr_depth;
						if(intra_cost < cost)
#else
//						intra_cost = intra_dist*1.05+DEPHT_ADD*curr_depth;
						intra_cost = intra_dist*(1.25-clip(((double)total_intra_partitions/(double)total_partitions), .0, .25))+DEPHT_ADD*curr_depth;
						
//						if(intra_cost+2000*curr_cu_info->sum<cost+2000*inter_sum)// && intra_cost<64*curr_cu_info->variance)
						if(intra_cost+clip(avg_distortion,100.,2000.)*curr_cu_info->sum<cost+clip(avg_distortion,100.,2000.)*inter_sum)// && intra_cost<64*curr_cu_info->variance)
//						if(intra_cost<cost)
#endif
						{	//we prefer intra and it is already in its buffer
							curr_cu_info->cost = intra_cost;
							curr_cu_info->distortion = intra_dist;
							curr_cu_info->sum = curr_cu_info->sum;
							curr_cu_info->prediction_mode = INTRA_MODE;
							consolidate_prediction_info(et, ctu, ctu_rd, curr_cu_info, curr_cu_info->cost, MAX_COST, FALSE, NULL);
						}
						else
						{	//we prefer inter, bring it back
//							get_back_consolidated_info(et, ctu, curr_cu_info, curr_depth);
							curr_cu_info->cost = cost;
							curr_cu_info->distortion = dist;
							curr_cu_info->sum = inter_sum;
							curr_cu_info->prediction_mode = INTER_MODE;
						}
						if(curr_depth <= max_processing_depth)//el = es para cuando et->max_intra_tr_depth!=4
						{
							int aux_depth;
//							cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
							abs_index = curr_cu_info->abs_index;
							num_part_in_cu  = curr_cu_info->num_part_in_cu;

							for(aux_depth=1;aux_depth<=max_processing_depth;aux_depth++)
							{
								synchronize_reference_buffs(et, curr_cu_info, &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[aux_depth+1], gcnt);	
								//for rd
			//					if(et->rd_mode!=RD_DIST_ONLY)
			//						CONSOLIDATE_INTRA_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu);
							}
							synchronize_reference_buffs_chroma(et, curr_cu_info, &et->decoded_mbs_wnd[0], &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
						}


//						consolidate_prediction_info(et, ctu, ctu_rd, curr_cu_info, curr_cu_info->cost, MAX_COST, FALSE, NULL);
					}

				}
			}
			else
			{
				curr_cu_info->cost = MAX_COST;
			}
		}

		depth_state[curr_depth]++;

		if(depth_state[curr_depth]==4)//la depth =1 lo hemos consolidado antes del bucle
		{
			int max_processing_depth;

			while(depth_state[curr_depth]==4 && curr_depth>0)//>0 pq consolidamos sobre el padre, 
			{

				depth_state[curr_depth] = 0;
				cost_sum[curr_depth] = 0;

				curr_depth--;
				curr_cu_info = curr_cu_info->parent;
			}
		}
		curr_cu_info++;// = parent_part_info->children[depth_state[curr_depth]];
	}


	curr_cu_info = &ctu->partition_list[0];
	abs_index = curr_cu_info->abs_index;
	curr_depth = curr_cu_info->depth;
	num_part_in_cu = curr_cu_info->num_part_in_cu;

	//if pred_depth==0 there is no NxN subdivision. we need to collect the information of the ctu
	if(et->max_pred_partition_depth==0)
	{
		CONSOLIDATE_ENC_INFO_BUFFS(et, ctu, curr_depth, abs_index, num_part_in_cu)
		SET_INTER_INFO_BUFFS(et, ctu, curr_cu_info, abs_index, num_part_in_cu, REF_PIC_LIST_0)	
	}


//	memset(&ctu->pred_mode[abs_index], INTER_MODE, num_part_in_cu*sizeof(ctu->pred_mode[0]));//signal all partitions as inter
	memset(&ctu->skipped[abs_index], FALSE, num_part_in_cu*sizeof(ctu->skipped[0]));//signal all partitions as non skipped
	return curr_cu_info->cost;
}




