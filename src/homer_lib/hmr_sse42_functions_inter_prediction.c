/*****************************************************************************
* hmr_sse42_functions_inter_prediction.c : homerHEVC encoding library
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



#include <string.h>
#include <math.h>

#include "hmr_sse42_primitives.h"
#include "hmr_sse42_macros.h"
#include "hmr_os_primitives.h"
#include "hmr_sse42_functions.h"


#define NTAPS_LUMA        8 ///< Number of taps for luma
#define NTAPS_CHROMA      4 ///< Number of taps for chroma
#define IF_INTERNAL_PREC 14 ///< Number of bits for internal precision
#define IF_FILTER_PREC    6 ///< Log2 of sum of filter taps
#define IF_INTERNAL_OFFS (1<<(IF_INTERNAL_PREC-1)) ///< Offset used internally


ALIGN(16) const int16_t sse_luma_filter_coeffs[4][NTAPS_LUMA] =
{
  {  0, 0,   0, 64,  0,   0, 0,  0 },
  { -1, 4, -10, 58, 17,  -5, 1,  0 },
  { -1, 4, -11, 40, 40, -11, 4, -1 },
  {  0, 1,  -5, 17, 58, -10, 4, -1 }
};

void sse_filter_copy_4xn(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int row, col;
  
	if ( is_first == is_last )
	{
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col+=4)
			{
				sse_64_storel_vector_u(dst+col, sse_128_loadlo_vector64(src+col));
			}
			src += src_stride;
			dst += dst_stride;
		}              
	}
	else if ( is_first )
	{
		int shift = IF_INTERNAL_PREC - bit_depth;
		__m128_i16 _128_offset = sse_128_vector_i16(IF_INTERNAL_OFFS);
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col+=8)
			{
				__m128_i16 val = sse_128_shift_l_i16(sse_128_loadlo_vector64(src+col), shift);
				sse_64_storel_vector_u(dst+col, sse_128_sub_i16(val,_128_offset));
			}
			src += src_stride;
			dst += dst_stride;
		}          
	}
	else
	{
		int shift = IF_INTERNAL_PREC - bit_depth;
		__m128_i16 _128_offset = sse_128_vector_i16(IF_INTERNAL_OFFS + (shift?(1 << (shift - 1)):0));
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col+=8)
			{
				__m128_i16 val = sse_128_add_i16(sse_128_loadlo_vector64(src+col),_128_offset);
				val =sse_128_shift_r_i16(val, shift);
				sse_64_storel_vector_u(dst+col, sse_128_convert_u8_i16(sse128_packs_i16_u8(val,val)));
			}
			src += src_stride;
			dst += dst_stride;
		}              
	}
}

void sse_filter_copy_8nxm(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int row, col;
  
	if ( is_first == is_last )
	{
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col+=8)
			{
				sse_128_store_vector_a(dst+col, sse_128_load_vector_u(src+col));
			}
			src += src_stride;
			dst += dst_stride;
		}              
	}
	else if ( is_first )
	{
		int shift = IF_INTERNAL_PREC - bit_depth;
		__m128_i16 _128_offset = sse_128_vector_i16(IF_INTERNAL_OFFS);
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col+=8)
			{
				__m128_i16 val = sse_128_shift_l_i16(sse_128_load_vector_u(src+col), shift);
				sse_128_store_vector_a(dst+col, sse_128_sub_i16(val,_128_offset));
			}
			src += src_stride;
			dst += dst_stride;
		}          
	}
	else
	{
		int shift = IF_INTERNAL_PREC - bit_depth;
		__m128_i16 _128_offset = sse_128_vector_i16(IF_INTERNAL_OFFS + (shift?(1 << (shift - 1)):0));
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col+=8)
			{
				__m128_i16 val = sse_128_add_i16(sse_128_load_vector_u(src+col),_128_offset);
				val =sse_128_shift_r_i16(val, shift);
				sse_128_store_vector_a(dst+col, sse_128_convert_u8_i16(sse128_packs_i16_u8(val,val)));
			}
			src += src_stride;
			dst += dst_stride;
		}              
	}
}



void sse_hmr_interpolation_filter_luma_8xn(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int num_taps = NTAPS_LUMA;//argument
	int y;
	int offset;
	int headRoom = IF_INTERNAL_PREC - bit_depth;
	int shift = IF_FILTER_PREC;
	__m128_i16 max_limit;
	__m128_i16 min_limit;

	if ( is_last )
	{
		shift += (is_first) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (is_first) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
		max_limit = sse_128_vector_i16((1<<bit_depth)-1);
		min_limit = sse_128_zero_vector();
	}
	else
	{
		shift -= (is_first) ? headRoom : 0;
		offset = (is_first) ? -IF_INTERNAL_OFFS << shift : 0;
	}

	if(!is_vertical)
	{
		//horizontal
		int16_t *src_buff_init = src - ( num_taps/2 - 1 );
		__m128_i16 coeff = sse_128_load_vector_a(sse_luma_filter_coeffs[fraction]);
		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)
		{
			int16_t *src_buff = src_buff_init;
			int16_t *dst_buff = dst;
			__m128_i16 pixels0 = sse_128_load_vector_u(src_buff);
			__m128_i16 pixels1 = sse_128_load_vector_u(src_buff+1);
			__m128_i16 pixels2 = sse_128_load_vector_u(src_buff+2);
			__m128_i16 pixels3 = sse_128_load_vector_u(src_buff+3);
			__m128_i16 pixels4 = sse_128_load_vector_u(src_buff+4);
			__m128_i16 pixels5 = sse_128_load_vector_u(src_buff+5);
			__m128_i16 pixels6 = sse_128_load_vector_u(src_buff+6);
			__m128_i16 pixels7 = sse_128_load_vector_u(src_buff+7);

			__m128_i32 sum0 = sse_128_madd_i16_i32(pixels0, coeff);
			__m128_i32 sum1 = sse_128_madd_i16_i32(pixels1, coeff);
			__m128_i32 sum2 = sse_128_madd_i16_i32(pixels2, coeff);
			__m128_i32 sum3 = sse_128_madd_i16_i32(pixels3, coeff);
			__m128_i32 sum4 = sse_128_madd_i16_i32(pixels4, coeff);
			__m128_i32 sum5 = sse_128_madd_i16_i32(pixels5, coeff);
			__m128_i32 sum6 = sse_128_madd_i16_i32(pixels6, coeff);
			__m128_i32 sum7 = sse_128_madd_i16_i32(pixels7, coeff);

			__m128_i32 sum01 = sse_128_hadd_i32(sum0, sum1);
			__m128_i32 sum23 = sse_128_hadd_i32(sum2, sum3);
			__m128_i32 sum0123 = sse_128_hadd_i32(sum01, sum23);
			__m128_i32 sum45 = sse_128_hadd_i32(sum4, sum5);
			__m128_i32 sum67 = sse_128_hadd_i32(sum6, sum7);
			__m128_i32 sum4567 = sse_128_hadd_i32(sum45, sum67);
			__m128_i16 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum0123, _m128_offset), shift), sse_128_shift_r_i32(sse_128_add_i32(sum4567, _m128_offset), shift));// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate
			if ( is_last)
			{
				val = sse_128_clip_16(val, min_limit, max_limit);
			}
			sse_128_store_vector_u(dst_buff, val);
			src_buff_init += src_stride;
			dst += dst_stride;
		}
	}
	else
	{
		//vertical
		int16_t *src_buff_init = src - (( num_taps/2 - 1 )*src_stride);
		__m128_i16 coeff = sse_128_load_vector_a(sse_luma_filter_coeffs[fraction]);
		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)
		{
			int16_t *src_buff = src_buff_init;
			int16_t *dst_buff = dst;

			__m128_i16 l0 = sse_128_load_vector_u(src_buff);
			__m128_i16 l1 = sse_128_load_vector_u(src_buff+src_stride);
			__m128_i16 l2 = sse_128_load_vector_u(src_buff+2*src_stride);
			__m128_i16 l3 = sse_128_load_vector_u(src_buff+3*src_stride);
			__m128_i16 l4 = sse_128_load_vector_u(src_buff+4*src_stride);
			__m128_i16 l5 = sse_128_load_vector_u(src_buff+5*src_stride);
			__m128_i16 l6 = sse_128_load_vector_u(src_buff+6*src_stride);
			__m128_i16 l7 = sse_128_load_vector_u(src_buff+7*src_stride);

			//TRANSPOSE_MATRIX_8x8_16BITS
			__m128_i16 l0l1_l = sse128_unpacklo_u16(l0,l1);		
			__m128_i16 l0l1_h = sse128_unpackhi_u16(l0,l1);		
			__m128_i16 l2l3_l = sse128_unpacklo_u16(l2,l3);		
			__m128_i16 l2l3_h = sse128_unpackhi_u16(l2,l3);		
			__m128_i16 l4l5_l = sse128_unpacklo_u16(l4,l5);		
			__m128_i16 l4l5_h = sse128_unpackhi_u16(l4,l5);		
			__m128_i16 l6l7_l = sse128_unpacklo_u16(l6,l7);		
			__m128_i16 l6l7_h = sse128_unpackhi_u16(l6,l7);		
																	
			__m128_i16 l0l1l2l3_ll = sse128_unpacklo_u32(l0l1_l,l2l3_l);	
			__m128_i16 l0l1l2l3_lh = sse128_unpackhi_u32(l0l1_l,l2l3_l);	
			__m128_i16 l0l1l2l3_hl = sse128_unpacklo_u32(l0l1_h,l2l3_h);	
			__m128_i16 l0l1l2l3_hh = sse128_unpackhi_u32(l0l1_h,l2l3_h);	
			__m128_i16 l4l5l6l7_ll = sse128_unpacklo_u32(l4l5_l,l6l7_l);	
			__m128_i16 l4l5l6l7_lh = sse128_unpackhi_u32(l4l5_l,l6l7_l);	
			__m128_i16 l4l5l6l7_hl = sse128_unpacklo_u32(l4l5_h,l6l7_h);	
			__m128_i16 l4l5l6l7_hh = sse128_unpackhi_u32(l4l5_h,l6l7_h);	
																	
			__m128_i16 col0 = sse128_unpacklo_u64(l0l1l2l3_ll,l4l5l6l7_ll);	
			__m128_i16 col1 = sse128_unpackhi_u64(l0l1l2l3_ll,l4l5l6l7_ll);	
			__m128_i16 col2 = sse128_unpacklo_u64(l0l1l2l3_lh,l4l5l6l7_lh);	
			__m128_i16 col3 = sse128_unpackhi_u64(l0l1l2l3_lh,l4l5l6l7_lh);	
			__m128_i16 col4 = sse128_unpacklo_u64(l0l1l2l3_hl,l4l5l6l7_hl);	
			__m128_i16 col5 = sse128_unpackhi_u64(l0l1l2l3_hl,l4l5l6l7_hl);	
			__m128_i16 col6 = sse128_unpacklo_u64(l0l1l2l3_hh,l4l5l6l7_hh);	
			__m128_i16 col7 = sse128_unpackhi_u64(l0l1l2l3_hh,l4l5l6l7_hh);	
			//----------------TRANSPOSE_MATRIX_8x8_16BITS---------------

			__m128_i32 sum0 = sse_128_madd_i16_i32(col0, coeff);
			__m128_i32 sum1 = sse_128_madd_i16_i32(col1, coeff);
			__m128_i32 sum2 = sse_128_madd_i16_i32(col2, coeff);
			__m128_i32 sum3 = sse_128_madd_i16_i32(col3, coeff);
			__m128_i32 sum4 = sse_128_madd_i16_i32(col4, coeff);
			__m128_i32 sum5 = sse_128_madd_i16_i32(col5, coeff);
			__m128_i32 sum6 = sse_128_madd_i16_i32(col6, coeff);
			__m128_i32 sum7 = sse_128_madd_i16_i32(col7, coeff);

			__m128_i32 sum01 = sse_128_hadd_i32(sum0, sum1);
			__m128_i32 sum23 = sse_128_hadd_i32(sum2, sum3);
			__m128_i32 sum0123 = sse_128_hadd_i32(sum01, sum23);
			__m128_i32 sum45 = sse_128_hadd_i32(sum4, sum5);
			__m128_i32 sum67 = sse_128_hadd_i32(sum6, sum7);
			__m128_i32 sum4567 = sse_128_hadd_i32(sum45, sum67);
			__m128_i16 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum0123, _m128_offset), shift), sse_128_shift_r_i32(sse_128_add_i32(sum4567, _m128_offset), shift));// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate
			if ( is_last)
			{
				val = sse_128_clip_16(val, min_limit, max_limit);
			}
			sse_128_store_vector_u(dst_buff, val);
			src_buff_init += src_stride;
			dst += dst_stride;
		}
	}
}


void sse_hmr_interpolation_filter_luma_8nxm(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int num_taps = NTAPS_LUMA;//argument
	int y, x;
	int offset;
	int headRoom = IF_INTERNAL_PREC - bit_depth;
	int shift = IF_FILTER_PREC;
	__m128_i16 max_limit;
	__m128_i16 min_limit;

	if ( is_last )
	{
		shift += (is_first) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (is_first) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
		max_limit = sse_128_vector_i16((1<<bit_depth)-1);
		min_limit = sse_128_zero_vector();
	}
	else
	{
		shift -= (is_first) ? headRoom : 0;
		offset = (is_first) ? -IF_INTERNAL_OFFS << shift : 0;
	}

	if(!is_vertical)
	{
		//horizontal
		int16_t *src_buff_init = src - ( num_taps/2 - 1 );
		__m128_i16 coeff = sse_128_load_vector_a(sse_luma_filter_coeffs[fraction]);
		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)
		{
			for (x = 0; x < width; x+=8)
			{
				int16_t *src_buff = src_buff_init+x;
				int16_t *dst_buff = dst+x;
				__m128_i16 pixels0 = sse_128_load_vector_u(src_buff);
				__m128_i16 pixels1 = sse_128_load_vector_u(src_buff+1);
				__m128_i16 pixels2 = sse_128_load_vector_u(src_buff+2);
				__m128_i16 pixels3 = sse_128_load_vector_u(src_buff+3);
				__m128_i16 pixels4 = sse_128_load_vector_u(src_buff+4);
				__m128_i16 pixels5 = sse_128_load_vector_u(src_buff+5);
				__m128_i16 pixels6 = sse_128_load_vector_u(src_buff+6);
				__m128_i16 pixels7 = sse_128_load_vector_u(src_buff+7);

				__m128_i32 sum0 = sse_128_madd_i16_i32(pixels0, coeff);
				__m128_i32 sum1 = sse_128_madd_i16_i32(pixels1, coeff);
				__m128_i32 sum2 = sse_128_madd_i16_i32(pixels2, coeff);
				__m128_i32 sum3 = sse_128_madd_i16_i32(pixels3, coeff);
				__m128_i32 sum4 = sse_128_madd_i16_i32(pixels4, coeff);
				__m128_i32 sum5 = sse_128_madd_i16_i32(pixels5, coeff);
				__m128_i32 sum6 = sse_128_madd_i16_i32(pixels6, coeff);
				__m128_i32 sum7 = sse_128_madd_i16_i32(pixels7, coeff);

				__m128_i32 sum01 = sse_128_hadd_i32(sum0, sum1);
				__m128_i32 sum23 = sse_128_hadd_i32(sum2, sum3);
				__m128_i32 sum0123 = sse_128_hadd_i32(sum01, sum23);
				__m128_i32 sum45 = sse_128_hadd_i32(sum4, sum5);
				__m128_i32 sum67 = sse_128_hadd_i32(sum6, sum7);
				__m128_i32 sum4567 = sse_128_hadd_i32(sum45, sum67);
				__m128_i16 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum0123, _m128_offset), shift), sse_128_shift_r_i32(sse_128_add_i32(sum4567, _m128_offset), shift));// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate
				if ( is_last)
				{
					val = sse_128_clip_16(val, min_limit, max_limit);
				}
				sse_128_store_vector_u(dst_buff, val);
			}
			src_buff_init += src_stride;
			dst += dst_stride;
		}
	}
	else
	{
		//vertical
		int16_t *src_buff_init = src - (( num_taps/2 - 1 )*src_stride);
		__m128_i16 coeff = sse_128_load_vector_a(sse_luma_filter_coeffs[fraction]);
		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)
		{
			for (x = 0; x < width; x+=8)
			{
				int16_t *src_buff = src_buff_init+x;
				int16_t *dst_buff = dst+x;

				__m128_i16 l0 = sse_128_load_vector_u(src_buff);
				__m128_i16 l1 = sse_128_load_vector_u(src_buff+src_stride);
				__m128_i16 l2 = sse_128_load_vector_u(src_buff+2*src_stride);
				__m128_i16 l3 = sse_128_load_vector_u(src_buff+3*src_stride);
				__m128_i16 l4 = sse_128_load_vector_u(src_buff+4*src_stride);
				__m128_i16 l5 = sse_128_load_vector_u(src_buff+5*src_stride);
				__m128_i16 l6 = sse_128_load_vector_u(src_buff+6*src_stride);
				__m128_i16 l7 = sse_128_load_vector_u(src_buff+7*src_stride);

				//TRANSPOSE_MATRIX_8x8_16BITS
				__m128_i16 l0l1_l = sse128_unpacklo_u16(l0,l1);		
				__m128_i16 l0l1_h = sse128_unpackhi_u16(l0,l1);		
				__m128_i16 l2l3_l = sse128_unpacklo_u16(l2,l3);		
				__m128_i16 l2l3_h = sse128_unpackhi_u16(l2,l3);		
				__m128_i16 l4l5_l = sse128_unpacklo_u16(l4,l5);		
				__m128_i16 l4l5_h = sse128_unpackhi_u16(l4,l5);		
				__m128_i16 l6l7_l = sse128_unpacklo_u16(l6,l7);		
				__m128_i16 l6l7_h = sse128_unpackhi_u16(l6,l7);		
																	
				__m128_i16 l0l1l2l3_ll = sse128_unpacklo_u32(l0l1_l,l2l3_l);	
				__m128_i16 l0l1l2l3_lh = sse128_unpackhi_u32(l0l1_l,l2l3_l);	
				__m128_i16 l0l1l2l3_hl = sse128_unpacklo_u32(l0l1_h,l2l3_h);	
				__m128_i16 l0l1l2l3_hh = sse128_unpackhi_u32(l0l1_h,l2l3_h);	
				__m128_i16 l4l5l6l7_ll = sse128_unpacklo_u32(l4l5_l,l6l7_l);	
				__m128_i16 l4l5l6l7_lh = sse128_unpackhi_u32(l4l5_l,l6l7_l);	
				__m128_i16 l4l5l6l7_hl = sse128_unpacklo_u32(l4l5_h,l6l7_h);	
				__m128_i16 l4l5l6l7_hh = sse128_unpackhi_u32(l4l5_h,l6l7_h);	
																	
				__m128_i16 col0 = sse128_unpacklo_u64(l0l1l2l3_ll,l4l5l6l7_ll);	
				__m128_i16 col1 = sse128_unpackhi_u64(l0l1l2l3_ll,l4l5l6l7_ll);	
				__m128_i16 col2 = sse128_unpacklo_u64(l0l1l2l3_lh,l4l5l6l7_lh);	
				__m128_i16 col3 = sse128_unpackhi_u64(l0l1l2l3_lh,l4l5l6l7_lh);	
				__m128_i16 col4 = sse128_unpacklo_u64(l0l1l2l3_hl,l4l5l6l7_hl);	
				__m128_i16 col5 = sse128_unpackhi_u64(l0l1l2l3_hl,l4l5l6l7_hl);	
				__m128_i16 col6 = sse128_unpacklo_u64(l0l1l2l3_hh,l4l5l6l7_hh);	
				__m128_i16 col7 = sse128_unpackhi_u64(l0l1l2l3_hh,l4l5l6l7_hh);	
				//----------------TRANSPOSE_MATRIX_8x8_16BITS---------------

				__m128_i32 sum0 = sse_128_madd_i16_i32(col0, coeff);
				__m128_i32 sum1 = sse_128_madd_i16_i32(col1, coeff);
				__m128_i32 sum2 = sse_128_madd_i16_i32(col2, coeff);
				__m128_i32 sum3 = sse_128_madd_i16_i32(col3, coeff);
				__m128_i32 sum4 = sse_128_madd_i16_i32(col4, coeff);
				__m128_i32 sum5 = sse_128_madd_i16_i32(col5, coeff);
				__m128_i32 sum6 = sse_128_madd_i16_i32(col6, coeff);
				__m128_i32 sum7 = sse_128_madd_i16_i32(col7, coeff);

				__m128_i32 sum01 = sse_128_hadd_i32(sum0, sum1);
				__m128_i32 sum23 = sse_128_hadd_i32(sum2, sum3);
				__m128_i32 sum0123 = sse_128_hadd_i32(sum01, sum23);
				__m128_i32 sum45 = sse_128_hadd_i32(sum4, sum5);
				__m128_i32 sum67 = sse_128_hadd_i32(sum6, sum7);
				__m128_i32 sum4567 = sse_128_hadd_i32(sum45, sum67);
				__m128_i16 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum0123, _m128_offset), shift), sse_128_shift_r_i32(sse_128_add_i32(sum4567, _m128_offset), shift));// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate
				if ( is_last)
				{
					val = sse_128_clip_16(val, min_limit, max_limit);
				}
				sse_128_store_vector_u(dst_buff, val);
			}
			src_buff_init += src_stride;
			dst += dst_stride;
		}
	}
}


/*void sse_interpolate_luma(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	if(fraction==0)
	{
		if(width<8)
		{
			sse_filter_copy_4xn(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
		}
		else
		{
			sse_filter_copy_8nxm(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
		}
	}
	else
	{
		if(width==8)
			sse_hmr_interpolation_filter_luma_8xn(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
		else
			sse_hmr_interpolation_filter_luma_8nxm(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
	}
}
*/


ALIGN(16) int16_t sse_chroma_filter_coeffs_horizontal[8][2*NTAPS_CHROMA] =
{
  {  0, 64,  0,  0,  0, 64,  0,  0 },
  { -2, 58, 10, -2, -2, 58, 10, -2 },
  { -4, 54, 16, -2, -4, 54, 16, -2 },
  { -6, 46, 28, -4, -6, 46, 28, -4 },
  { -4, 36, 36, -4, -4, 36, 36, -4 },
  { -4, 28, 46, -6, -4, 28, 46, -6 },
  { -2, 16, 54, -4, -2, 16, 54, -4 },
  { -2, 10, 58, -2, -2, 10, 58, -2 }
};


//this function differs to hmr_interpolate_chroma and to HM in higher range values because val here is int32 and in HM is int16 and so it may overflow with high values - I think is an error in HM as it is saturated later
void sse_hmr_interpolate_chroma_4xn(int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int num_taps = NTAPS_CHROMA;//argument

	int offset = 0;
	int headRoom = IF_INTERNAL_PREC - bit_depth;
	int shift = IF_FILTER_PREC;
	int y = 0;
	__m128_i16 max_limit;
	__m128_i16 min_limit;

	if ( is_last )
	{
		shift += (is_first) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (is_first) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
		max_limit = sse_128_vector_i16((1<<bit_depth)-1);
		min_limit = sse_128_zero_vector();
	}
	else
	{
		shift -= (is_first) ? headRoom : 0;
		offset = (is_first) ? -IF_INTERNAL_OFFS << shift : 0;
	}



	if(!is_vertical)
	{
		//horizontal
		int16_t *ref_buff = reference_buff - (( num_taps/2 - 1 ));
		__m128_i16 coeff = sse_128_load_vector_u(sse_chroma_filter_coeffs_horizontal[fraction]);
		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)//row
		{
			__m128_i16 pixels01 = sse128_unpacklo_u64(sse_128_load_vector_u(ref_buff), sse_128_load_vector_u(ref_buff+1));
			__m128_i16 pixels23 = sse128_unpacklo_u64(sse_128_load_vector_u(ref_buff+2), sse_128_load_vector_u(ref_buff+3));
			__m128_i32 sum12 = sse_128_madd_i16_i32(pixels01, coeff);
			__m128_i32 sum23 = sse_128_madd_i16_i32(pixels23, coeff);
			__m128_i32 sum = sse_128_hadd_i32(sum12, sum23);
			__m128_i32 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum, _m128_offset), shift), _m128_offset);// in HM val is type short. 
			if ( is_last)
			{
				val = sse_128_clip_16(val, min_limit, max_limit);
				//val = sse_128_convert_u8_i16(sse128_packs_i16_u8(val,val));
			}
			sse_64_storel_vector_u(pred_buff, val);

			ref_buff += reference_buff_stride;
			pred_buff += pred_buff_stride;
		}
	}
	else
	{
		//vertical
		int16_t *ref_buff = reference_buff - (( num_taps/2 - 1 )*reference_buff_stride);
		__m128_i16 coeff = sse_128_load_vector_u(sse_chroma_filter_coeffs_horizontal[fraction]);

		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)//row
		{
			__m128_i16 l0l1 = sse128_unpacklo_u16(sse_128_load_vector_u(ref_buff), sse_128_load_vector_u(ref_buff+reference_buff_stride));
			__m128_i16 l2l3 = sse128_unpacklo_u16(sse_128_load_vector_u(ref_buff+2*reference_buff_stride), sse_128_load_vector_u(ref_buff+3*reference_buff_stride));
			__m128_i16 c0c1 = sse128_unpacklo_u32(l0l1,l2l3);
			__m128_i16 c2c3 = sse128_unpackhi_u32(l0l1,l2l3);
			__m128_i32 sum12 = sse_128_madd_i16_i32(c0c1, coeff);
			__m128_i32 sum23 = sse_128_madd_i16_i32(c2c3, coeff);
			__m128_i32 sum = sse_128_hadd_i32(sum12, sum23);
			__m128_i32 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum, _m128_offset), shift), _m128_offset);//_m128_offset is whatever
			if ( is_last)
			{
				val = sse_128_clip_16(val, min_limit, max_limit);
				//val = sse_128_convert_u8_i16(sse128_packs_i16_u8(val,val));
			}
			sse_64_storel_vector_u(pred_buff, val);

			ref_buff += reference_buff_stride;
			pred_buff += pred_buff_stride;
		}	
	}
}


//this function differs to hmr_interpolate_chroma and to HM in higher range values because val here is int32 and in HM is int16 and so it may overflow with high values - I think is an error in HM as it is saturated later
void sse_hmr_interpolate_chroma_8xn(int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int num_taps = NTAPS_CHROMA;//argument

	int offset = 0;
	int headRoom = IF_INTERNAL_PREC - bit_depth;
	int shift = IF_FILTER_PREC;
	int y;
	__m128_i16 max_limit;
	__m128_i16 min_limit;

	if ( is_last )
	{
		shift += (is_first) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (is_first) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
		max_limit = sse_128_vector_i16((1<<bit_depth)-1);
		min_limit = sse_128_zero_vector();
	}
	else
	{
		shift -= (is_first) ? headRoom : 0;
		offset = (is_first) ? -IF_INTERNAL_OFFS << shift : 0;
	}

	if(!is_vertical)
	{
		//horizontal
		int16_t *ref_buff = reference_buff - (( num_taps/2 - 1 ));
		__m128_i16 coeff = sse_128_load_vector_u(sse_chroma_filter_coeffs_horizontal[fraction]);
		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)
		{
			__m128_i16 pixels04 = sse_128_load_vector_u(ref_buff);
			__m128_i16 pixels15 = sse_128_load_vector_u(ref_buff+1);
			__m128_i16 pixels26 = sse_128_load_vector_u(ref_buff+2);
			__m128_i16 pixels37 = sse_128_load_vector_u(ref_buff+3);
			__m128_i32 sum04 = sse_128_madd_i16_i32(pixels04, coeff);
			__m128_i32 sum15 = sse_128_madd_i16_i32(pixels15, coeff);
			__m128_i32 sum26 = sse_128_madd_i16_i32(pixels26, coeff);
			__m128_i32 sum37 = sse_128_madd_i16_i32(pixels37, coeff);
			__m128_i32 sum0426 = sse_128_hadd_i32(sum04, sum26);
			__m128_i32 sum1537 = sse_128_hadd_i32(sum15, sum37);
			__m128_i32 sum0145 = sse128_unpacklo_u32(sum0426, sum1537);
			__m128_i32 sum2367 = sse128_unpackhi_u32(sum0426, sum1537);
			__m128_i32 sum0123 = sse128_unpacklo_u64(sum0145, sum2367);
			__m128_i32 sum4567 = sse128_unpackhi_u64(sum0145, sum2367);
			__m128_i16 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum0123, _m128_offset), shift), sse_128_shift_r_i32(sse_128_add_i32(sum4567, _m128_offset), shift));// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate
			if ( is_last)
			{
				val = sse_128_clip_16(val, min_limit, max_limit);
//				val = sse_128_convert_u8_i16(sse128_packs_i16_u8(val,val));
			}
			sse_128_store_vector_u(pred_buff, val);

			ref_buff += reference_buff_stride;
			pred_buff += pred_buff_stride;
		}
	}
	else
	{
		//vertical
		int16_t *ref_buff = reference_buff - (( num_taps/2 - 1 )*reference_buff_stride);
		__m128_i16 coeff = sse_128_load_vector_u(sse_chroma_filter_coeffs_horizontal[fraction]);

		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)//row
		{
			__m128_i16 l0 = sse_128_load_vector_u(ref_buff);
			__m128_i16 l1 = sse_128_load_vector_u(ref_buff+reference_buff_stride);
			__m128_i16 l2 = sse_128_load_vector_u(ref_buff+2*reference_buff_stride);
			__m128_i16 l3 = sse_128_load_vector_u(ref_buff+3*reference_buff_stride);
			__m128_i16 l0l1_l = sse128_unpacklo_u16(l0, l1);
			__m128_i16 l0l1_h = sse128_unpackhi_u16(l0, l1);
			__m128_i16 l2l3_l = sse128_unpacklo_u16(l2, l3);
			__m128_i16 l2l3_h = sse128_unpackhi_u16(l2, l3);
			__m128_i16 pixels01 = sse128_unpacklo_u32(l0l1_l, l2l3_l);
			__m128_i16 pixels23 = sse128_unpackhi_u32(l0l1_l, l2l3_l);
			__m128_i16 pixels45 = sse128_unpacklo_u32(l0l1_h, l2l3_h);
			__m128_i16 pixels67 = sse128_unpackhi_u32(l0l1_h, l2l3_h);
			__m128_i32 sum01 = sse_128_madd_i16_i32(pixels01, coeff);
			__m128_i32 sum23 = sse_128_madd_i16_i32(pixels23, coeff);
			__m128_i32 sum45 = sse_128_madd_i16_i32(pixels45, coeff);
			__m128_i32 sum67 = sse_128_madd_i16_i32(pixels67, coeff);
			__m128_i32 sum0123 = sse_128_hadd_i32(sum01, sum23);
			__m128_i32 sum4567 = sse_128_hadd_i32(sum45, sum67);
			__m128_i16 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum0123, _m128_offset), shift), sse_128_shift_r_i32(sse_128_add_i32(sum4567, _m128_offset), shift));// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate
			if ( is_last)
			{
				val = sse_128_clip_16(val, min_limit, max_limit);
//				val = sse_128_convert_u8_i16(sse128_packs_i16_u8(val,val));
			}
			sse_128_store_vector_u(pred_buff, val);

			ref_buff += reference_buff_stride;
			pred_buff += pred_buff_stride;
		}	
	}
}


//this function differs to hmr_interpolate_chroma and to HM in higher range values because val here is int32 and in HM is int16 and so it may overflow with high values - I think there is an error in HM as it is saturated later
void sse_hmr_interpolate_chroma_nxn(int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int num_taps = NTAPS_CHROMA;//argument

	int offset = 0;
	int headRoom = IF_INTERNAL_PREC - bit_depth;
	int shift = IF_FILTER_PREC;
	int y, x;
	__m128_i16 max_limit;
	__m128_i16 min_limit;

	if ( is_last )
	{
		shift += (is_first) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (is_first) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
		max_limit = sse_128_vector_i16((1<<bit_depth)-1);
		min_limit = sse_128_zero_vector();
	}
	else
	{
		shift -= (is_first) ? headRoom : 0;
		offset = (is_first) ? -IF_INTERNAL_OFFS << shift : 0;
	}

	if(!is_vertical)
	{
		//horizontal
		int16_t *ref_buff_init = reference_buff - (( num_taps/2 - 1 ));
		__m128_i16 coeff = sse_128_load_vector_u(sse_chroma_filter_coeffs_horizontal[fraction]);
		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)
		{
			for (x = 0; x < width; x+=8)
			{
				int16_t *ref_buff = ref_buff_init+x;
				int16_t *predict_buff = pred_buff+x;
				__m128_i16 pixels04 = sse_128_load_vector_u(ref_buff);
				__m128_i16 pixels15 = sse_128_load_vector_u(ref_buff+1);
				__m128_i16 pixels26 = sse_128_load_vector_u(ref_buff+2);
				__m128_i16 pixels37 = sse_128_load_vector_u(ref_buff+3);
				__m128_i32 sum04 = sse_128_madd_i16_i32(pixels04, coeff);
				__m128_i32 sum15 = sse_128_madd_i16_i32(pixels15, coeff);
				__m128_i32 sum26 = sse_128_madd_i16_i32(pixels26, coeff);
				__m128_i32 sum37 = sse_128_madd_i16_i32(pixels37, coeff);
				__m128_i32 sum0426 = sse_128_hadd_i32(sum04, sum26);
				__m128_i32 sum1537 = sse_128_hadd_i32(sum15, sum37);
				__m128_i32 sum0145 = sse128_unpacklo_u32(sum0426, sum1537);
				__m128_i32 sum2367 = sse128_unpackhi_u32(sum0426, sum1537);
				__m128_i32 sum0123 = sse128_unpacklo_u64(sum0145, sum2367);
				__m128_i32 sum4567 = sse128_unpackhi_u64(sum0145, sum2367);
				__m128_i16 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum0123, _m128_offset), shift), sse_128_shift_r_i32(sse_128_add_i32(sum4567, _m128_offset), shift));// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate
				if ( is_last)
				{
					val = sse_128_clip_16(val, min_limit, max_limit);
//					val = sse_128_convert_u8_i16(sse128_packs_i16_u8(val,val));
				}
				sse_128_store_vector_u(predict_buff, val);
			}
			ref_buff_init += reference_buff_stride;
			pred_buff += pred_buff_stride;
		}
	}
	else
	{
		//vertical
		int16_t *ref_buff_init = reference_buff - (( num_taps/2 - 1 )*reference_buff_stride);
		__m128_i16 coeff = sse_128_load_vector_u(sse_chroma_filter_coeffs_horizontal[fraction]);

		__m128_i32 _m128_offset = sse_128_vector_i32(offset);

		for (y = 0; y < height; y++)
		{
			for (x = 0; x < width; x+=8)
			{
				int16_t *ref_buff = ref_buff_init+x;
				int16_t *predict_buff = pred_buff+x;
				__m128_i16 l0 = sse_128_load_vector_u(ref_buff);
				__m128_i16 l1 = sse_128_load_vector_u(ref_buff+reference_buff_stride);
				__m128_i16 l0l1_l = sse128_unpacklo_u16(l0, l1);
				__m128_i16 l0l1_h = sse128_unpackhi_u16(l0, l1);
				__m128_i16 l2 = sse_128_load_vector_u(ref_buff+2*reference_buff_stride);
				__m128_i16 l3 = sse_128_load_vector_u(ref_buff+3*reference_buff_stride);
				__m128_i16 l2l3_l = sse128_unpacklo_u16(l2, l3);
				__m128_i16 l2l3_h = sse128_unpackhi_u16(l2, l3);
				__m128_i16 pixels01 = sse128_unpacklo_u32(l0l1_l, l2l3_l);
				__m128_i16 pixels23 = sse128_unpackhi_u32(l0l1_l, l2l3_l);
				__m128_i16 pixels45 = sse128_unpacklo_u32(l0l1_h, l2l3_h);
				__m128_i16 pixels67 = sse128_unpackhi_u32(l0l1_h, l2l3_h);
				__m128_i32 sum01 = sse_128_madd_i16_i32(pixels01, coeff);
				__m128_i32 sum23 = sse_128_madd_i16_i32(pixels23, coeff);
				__m128_i32 sum45 = sse_128_madd_i16_i32(pixels45, coeff);
				__m128_i32 sum67 = sse_128_madd_i16_i32(pixels67, coeff);
				__m128_i32 sum0123 = sse_128_hadd_i32(sum01, sum23);
				__m128_i32 sum4567 = sse_128_hadd_i32(sum45, sum67);
				__m128_i16 val = sse128_packs_i32_i16(sse_128_shift_r_i32(sse_128_add_i32(sum0123, _m128_offset), shift), sse_128_shift_r_i32(sse_128_add_i32(sum4567, _m128_offset), shift));// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate
				if ( is_last)
				{
					val = sse_128_clip_16(val, min_limit, max_limit);
//					val = sse_128_convert_u8_i16(sse128_packs_i16_u8(val,val));
				}
				sse_128_store_vector_u(predict_buff, val);
			}
			ref_buff_init += reference_buff_stride;
			pred_buff += pred_buff_stride;
		}	
	}
}



void sse_interpolate_luma(int16_t *src, int src_stride, int16_t *dst, int dst_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	if(fraction==0)
	{
		if(width<8)
		{
			sse_filter_copy_4xn(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
		}
		else
		{
			sse_filter_copy_8nxm(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
		}
	}
	else
	{
		if(width==8)
			sse_hmr_interpolation_filter_luma_8xn(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
		else
			sse_hmr_interpolation_filter_luma_8nxm(src, src_stride, dst, dst_stride, fraction, width, height, is_vertical, is_first, is_last);
	}
}

void sse_interpolate_chroma(int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last)
{
	if(fraction==0)
	{
		if(width<4)
		{
			int iiiii=0;
		}
		else if(width<8)
		{
			sse_filter_copy_4xn(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, fraction, width, height, is_vertical, is_first, is_last);
		}
		else
		{
			sse_filter_copy_8nxm(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, fraction, width, height, is_vertical, is_first, is_last);
		}
	}
	else
	{
		if(width==4)
			sse_hmr_interpolate_chroma_4xn(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, fraction, width, height, is_vertical, is_first, is_last);
		else if(width==8)
			sse_hmr_interpolate_chroma_8xn(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, fraction, width, height, is_vertical, is_first, is_last);
		else
			sse_hmr_interpolate_chroma_nxn(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, fraction, width, height, is_vertical, is_first, is_last);
	}
}



void sse_weighted_average_motion_4xm(int16_t* src0, int src0_stride, int16_t* src1, int src1_stride, int16_t* dst, int dst_stride, int height, int width, int bit_depth)
{
	int x, y;
	int max_pix_val = (1<<bit_depth)-1;
	int shift = IF_INTERNAL_PREC + 1 - bit_depth;
	int offset = ( 1 << ( shift - 1 ) ) + 2 * IF_INTERNAL_OFFS;

	__m128_i32 _m128_offset = sse_128_vector_i32(offset);
	__m128_i16 _m128_min = sse_128_vector_i16(0);
	__m128_i16 _m128_max = sse_128_vector_i16(max_pix_val);
	for ( y = 0; y < height; y++ )
	{
//		for ( x = 0; x < width; x+=8)
		{
			__m128_i32 s0_0 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src0));
			__m128_i32 s1_0 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src1));

			__m128_i32 sum0 = sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(s0_0, s1_0), _m128_offset), shift);

			__m128_i16 val = sse128_packs_i32_i16(sum0, sum0);// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate - if computed in 16 bits this overflows the 16 bit range
			val = sse_128_clip_16(val, _m128_min, _m128_max);
			sse_64_storel_vector_u(dst, val);
		}
		src0 += src0_stride;
		src1 += src1_stride;
		dst +=dst_stride;
	}
}



void sse_weighted_average_motion_8xm(int16_t* src0, int src0_stride, int16_t* src1, int src1_stride, int16_t* dst, int dst_stride, int height, int width, int bit_depth)
{
	int x, y;
	int max_pix_val = (1<<bit_depth)-1;
	int shift = IF_INTERNAL_PREC + 1 - bit_depth;
	int offset = ( 1 << ( shift - 1 ) ) + 2 * IF_INTERNAL_OFFS;

	__m128_i32 _m128_offset = sse_128_vector_i32(offset);
	__m128_i16 _m128_min = sse_128_vector_i16(0);
	__m128_i16 _m128_max = sse_128_vector_i16(max_pix_val);
	for ( y = 0; y < height; y++ )
	{
//		for ( x = 0; x < width; x+=8)
		{
			__m128_i32 s0_0 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src0));
			__m128_i32 s1_0 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src1));
			__m128_i32 s0_1 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src0+4));
			__m128_i32 s1_1 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src1+4));

			__m128_i32 sum0 = sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(s0_0, s1_0), _m128_offset), shift);
			__m128_i32 sum1 = sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(s0_1, s1_1), _m128_offset), shift);
			__m128_i16 val = sse128_packs_i32_i16(sum0, sum1);// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate - if computed in 16 bits this overflows the 16 bit range
			val = sse_128_clip_16(val, _m128_min, _m128_max);
			sse_128_store_vector_u(dst, val);
			//		  dst[x] = clip(( ( src0[x] + src1[x] + offset ) >> shiftNum ),0,max_pix_val);
		}
		src0 += src0_stride;
		src1 += src1_stride;
		dst +=dst_stride;
	}
}

void sse_weighted_average_motion_8nxm(int16_t* src0, int src0_stride, int16_t* src1, int src1_stride, int16_t* dst, int dst_stride, int height, int width, int bit_depth)
{
	int x, y;
	int max_pix_val = (1<<bit_depth)-1;
	int shift = IF_INTERNAL_PREC + 1 - bit_depth;
	int offset = ( 1 << ( shift - 1 ) ) + 2 * IF_INTERNAL_OFFS;

	__m128_i32 _m128_offset = sse_128_vector_i32(offset);
	__m128_i16 _m128_min = sse_128_vector_i16(0);
	__m128_i16 _m128_max = sse_128_vector_i16(max_pix_val);
	for ( y = 0; y < height; y++ )
	{
		for ( x = 0; x < width; x+=8)
		{
			__m128_i32 s0_0 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src0+x));
			__m128_i32 s1_0 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src1+x));
			__m128_i32 s0_1 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src0+x+4));
			__m128_i32 s1_1 = sse_128_convert_i16_i32(sse_128_loadlo_vector64(src1+x+4));

			__m128_i32 sum0 = sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(s0_0, s1_0), _m128_offset), shift);
			__m128_i32 sum1 = sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(s0_1, s1_1), _m128_offset), shift);
			__m128_i16 val = sse128_packs_i32_i16(sum0, sum1);// in HM val is type short. it is equivalent to packing whithout saturation, and then saturate - if computed in 16 bits this overflows the 16 bit range because of the "+ 2 * IF_INTERNAL_OFFS" of the offset
			val = sse_128_clip_16(val, _m128_min, _m128_max);
			sse_128_store_vector_u(dst+x, val);
		}
		src0 += src0_stride;
		src1 += src1_stride;
		dst +=dst_stride;
	}
}



void sse_weighted_average_motion(int16_t* src0, int src0_stride, int16_t* src1, int src1_stride, int16_t* dst, int dst_stride, int height, int width, int bit_depth)
{
	if(width == 4)
		sse_weighted_average_motion_4xm(src0, src0_stride, src1, src1_stride, dst, dst_stride, height, width, bit_depth);
	else if(width == 8)
		sse_weighted_average_motion_8xm(src0, src0_stride, src1, src1_stride, dst, dst_stride, height, width, bit_depth);
	else
		sse_weighted_average_motion_8nxm(src0, src0_stride, src1, src1_stride, dst, dst_stride, height, width, bit_depth);
}