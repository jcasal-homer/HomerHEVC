/*****************************************************************************
* hmr_sse42_functions_quant.c : homerHEVC encoding library
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

#include "hmr_os_primitives.h"
#include "hmr_sse42_primitives.h"
#include "hmr_sse42_macros.h"
#include "hmr_sse42_functions.h"





void sse_aligned_quant(henc_thread_t* et, ctu_info_t *ctu, int16_t* src, int16_t* dst, int scan_mode, int depth, int comp, int cu_mode, int is_intra, int *ac_sum, int cu_size)
{
	int16_t *psrc = src, *pdst = dst;
	int iLevel, auxLevel;
	int  iSign;
	picture_t *currpict = &et->ed->current_pict;	//esto y
	slice_t *currslice = &currpict->slice;		//esto deberia pasarse como parametro
	int inv_depth = (et->max_cu_size_shift - (depth + (comp!=Y_COMP)));//ed->max_cu_size_shift
	uint32_t *scan = et->ed->scan_pyramid[scan_mode][inv_depth-1];
	int scan_type = (is_intra?0:3) + comp;
	int32_t *quant_coeff = et->ed->quant_pyramid[inv_depth-2][scan_type][ctu->rem];
	uint bit_depth = et->bit_depth;
	uint transform_shift = MAX_TR_DYNAMIC_RANGE - et->bit_depth - inv_depth;
	int qbits = QUANT_SHIFT + ctu->per + transform_shift;                // Right shift of non-RDOQ quantizer;  level = (coeff*uiQ + offset)>>q_bits
	int qbits8 = qbits-8;
	int add = (currslice->slice_type==I_SLICE ? 171 : 85) << (qbits-9);
	short *deltaU = et->aux_buff;//[32*32];//hay que usar el buffer auxiliar
	int n;

	__m128i isum = sse_128_vector_i16(0);

	__m128i _128add = sse_128_vector_i32(add);
	__m128i _128one = sse_128_vector_i16(1);
	__m128i _128zero = sse_128_vector_i16(0);

	for(n=0;n<(cu_size*cu_size);n+=16)
	{
		__m128i q0 = sse_128_load_vector_a(quant_coeff);

		__m128i l0l1 = sse_128_load_vector_a(psrc);

		__m128i l0l1sign = sse_128_sign_16(_128one, l0l1);

		__m128i l0l1abs = sse_128_abs_i16(l0l1);

		__m128i l0abs = sse128_unpacklo_u16(l0l1abs,_128zero);

		__m128i aux0 = sse_128_mul_i32(l0abs,q0);
		__m128i c0 = sse_128_shift_r_i32(sse_128_add_i32(aux0,_128add), qbits); 
		__m128i delta0 = sse_128_shift_r_i32(sse_128_sub_i32(aux0, sse_128_shift_l_i32(c0, qbits)),qbits8);
		__m128i l1abs = sse128_unpackhi_u16(l0l1abs,_128zero);
		__m128i q1 = sse_128_load_vector_a(quant_coeff+4);

		__m128i aux1 = sse_128_mul_i32(l1abs,q1);
		__m128i c1 = sse_128_shift_r_i32(sse_128_add_i32(aux1,_128add), qbits); 
		__m128i delta1 = sse_128_shift_r_i32(sse_128_sub_i32(aux1, sse_128_shift_l_i32(c1, qbits)),qbits8);
		
		isum = sse_128_hadd_i32(isum, sse_128_hadd_i32(c0,c1));

		sse_128_store_vector_a(pdst, sse_128_mul_i16(l0l1sign, sse128_packs_u32_u16(c0,c1)));
		sse_128_store_vector_a(deltaU, sse128_packs_u32_u16(delta0,delta1));

		quant_coeff+=8;
		psrc+=8;
		pdst+=8;
		deltaU+=8;

		q0 = sse_128_load_vector_a(quant_coeff);

		l0l1 = sse_128_load_vector_a(psrc);

		l0l1sign = sse_128_sign_16(_128one, l0l1);

		l0l1abs = sse_128_abs_i16(l0l1);

		l0abs = sse128_unpacklo_u16(l0l1abs,_128zero);

		aux0 = sse_128_mul_i32(l0abs,q0);
		c0 = sse_128_shift_r_i32(sse_128_add_i32(aux0,_128add), qbits); 
		delta0 = sse_128_shift_r_i32(sse_128_sub_i32(aux0, sse_128_shift_l_i32(c0, qbits)),qbits8);

		isum = sse_128_hadd_i32(isum,c0);

		l1abs = sse128_unpackhi_u16(l0l1abs,_128zero);
		q1 = sse_128_load_vector_a(quant_coeff+4);

		aux1 = sse_128_mul_i32(l1abs,q1);
		c1 = sse_128_shift_r_i32(sse_128_add_i32(aux1,_128add), qbits); 
		delta1 = sse_128_shift_r_i32(sse_128_sub_i32(aux1, sse_128_shift_l_i32(c1, qbits)),qbits8);

		isum = sse_128_hadd_i32(isum,c1);

		sse_128_store_vector_a(pdst, sse_128_mul_i16(l0l1sign, sse128_packs_u32_u16(c0,c1)));
		sse_128_store_vector_a(deltaU, sse128_packs_u32_u16(delta0,delta1));

		quant_coeff+=8;
		psrc+=8;
		pdst+=8;
		deltaU+=8;
	}
	*ac_sum = sse_128_get_data_u32(isum,0)+sse_128_get_data_u32(isum,1)+sse_128_get_data_u32(isum,2)+sse_128_get_data_u32(isum,3);

	deltaU = et->aux_buff;
	if(et->pps->sign_data_hiding_flag)
	{
		if(*ac_sum>=2)
		{
			sign_bit_hidding( dst, src, scan, deltaU, cu_size, cu_size ) ;
		}
	}
}


void sse_aligned_inv_quant(henc_thread_t* et, ctu_info_t *ctu, short *src, short *dst, int depth, int comp, int is_intra, int cu_size)
{
	int iLevel;
	int inv_depth = (et->max_cu_size_shift - (depth+(comp!=Y_COMP)));//ed->max_cu_size_shift
	int scan_type = is_intra?0:3 + comp;
	int32_t *dequant_coeff = et->ed->dequant_pyramid[inv_depth-2][scan_type][ctu->rem];
	uint bit_depth = et->bit_depth;
	uint transform_shift = MAX_TR_DYNAMIC_RANGE - et->bit_depth - inv_depth;
	int iq_shift = QUANT_IQUANT_SHIFT - QUANT_SHIFT - transform_shift + 4;
	int iq_add = (iq_shift>ctu->per)? 1 << (iq_shift - ctu->per - 1): 0;
	int n;

	__m128i _128one = sse_128_vector_i16(1);
	__m128i _128zero = sse_128_vector_i16(0);

	if(iq_shift>ctu->per)
	{
		__m128i _128add = sse_128_vector_i32(iq_add);
		iq_shift=iq_shift-ctu->per;

		for(n=0;n<(cu_size*cu_size);n+=16)
		{
			__m128i q0 = sse_128_load_vector_a(dequant_coeff);

			__m128i l0l1 = sse_128_load_vector_a(src);

			__m128i l0 = sse_128_convert_i16_i32(l0l1);

			__m128i aux0 = sse_128_mul_i32(l0,q0);
			__m128i c0 = sse_128_shift_r_i32(sse_128_add_i32(aux0,_128add), iq_shift); 

			__m128i l1 = sse_128_convert_i16_i32(sse128_unpackhi_u64(l0l1,l0l1));
			__m128i q1 = sse_128_load_vector_a(dequant_coeff+4);

			__m128i aux1 = sse_128_mul_i32(l1,q1);
			__m128i c1 = sse_128_shift_r_i32(sse_128_add_i32(aux1,_128add), iq_shift); 

			sse_128_store_vector_a(dst, sse128_packs_u32_u16(c0,c1));
			src+=8;
			dst+=8;
			dequant_coeff+=8;

			q0 = sse_128_load_vector_a(dequant_coeff);

			l0l1 = sse_128_load_vector_a(src);

			l0 = sse_128_convert_i16_i32(l0l1);

			aux0 = sse_128_mul_i32(l0,q0);
			c0 = sse_128_shift_r_i32(sse_128_add_i32(aux0,_128add), iq_shift); 

			l1 = sse_128_convert_i16_i32(sse128_unpackhi_u64(l0l1,l0l1));
			q1 = sse_128_load_vector_a(dequant_coeff+4);

			aux1 = sse_128_mul_i32(l1,q1);
			c1 = sse_128_shift_r_i32(sse_128_add_i32(aux1,_128add), iq_shift); 

			sse_128_store_vector_a(dst, sse128_packs_u32_u16(c0,c1));
			src+=8;
			dst+=8;
			dequant_coeff+=8;
		}
	}
	else 
	{
		iq_shift=(ctu->per-iq_shift);

		for(n=0;n<(cu_size*cu_size);n+=16)
		{
			__m128i q0 = sse_128_load_vector_a(dequant_coeff);

			__m128i l0l1 = sse_128_load_vector_a(src);

			__m128i l0 = sse_128_convert_i16_i32(l0l1);

			__m128i aux0 = sse_128_mul_i32(l0,q0);
			__m128i c0 = sse_128_shift_l_i32(aux0, iq_shift); 

			__m128i l1 = sse_128_convert_i16_i32(sse128_unpackhi_u64(l0l1,l0l1));
			__m128i q1 = sse_128_load_vector_a(dequant_coeff+4);

			__m128i aux1 = sse_128_mul_i32(l1,q1);
			__m128i c1 = sse_128_shift_l_i32(aux1, iq_shift); 

			sse_128_store_vector_a(dst, sse128_packs_u32_u16(c0,c1));
			src+=8;
			dst+=8;
			dequant_coeff+=8;

			q0 = sse_128_load_vector_a(dequant_coeff);

			l0l1 = sse_128_load_vector_a(src);

			l0 = sse_128_convert_i16_i32(l0l1);

			aux0 = sse_128_mul_i32(l0,q0);
			c0 = sse_128_shift_l_i32(aux0, iq_shift); 

			l1 = sse_128_convert_i16_i32(sse128_unpackhi_u64(l0l1,l0l1));
			q1 = sse_128_load_vector_a(dequant_coeff+4);

			aux1 = sse_128_mul_i32(l1,q1);
			c1 = sse_128_shift_l_i32(aux1, iq_shift); 

			sse_128_store_vector_a(dst, sse128_packs_u32_u16(c0,c1));
			src+=8;
			dst+=8;
			dequant_coeff+=8;
		}
	}


}

