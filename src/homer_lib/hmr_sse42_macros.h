/*****************************************************************************
* hmr_sse42_functions_macros.c : homerHEVC encoding library
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


#ifndef __HOMER_HEVC_SSE42_MACROS__
#define __HOMER_HEVC_SSE42_MACROS__


#include "hmr_sse42_primitives.h"

#define CALC_ALIGNED_SAD_2x4(result, src_ln1, src_ln2, pred_ln1, pred_ln2)	\
	result = sse_128_add_i64(result, sse_128_sad_u8(sse_128_loadlo_vector64(&sse128_unpacklo_u8(sse_128_loadlo_vector64(src_ln1),sse_128_loadlo_vector64(src_ln2))),sse_128_loadlo_vector64(&sse128_unpacklo_u8(sse_128_loadlo_vector64(pred_ln1),sse_128_loadlo_vector64(pred_ln2)))));	
	//result = sse_128_add_i64(result, sse_128_sad_u8(sse_128_loadlo_vector64(&sse128_unpacklo_u8(sse_128_loadlo_vector64(src_ln1),sse_128_loadlo_vector64(src_ln2))),sse_128_loadlo_vector64(&sse128_unpacklo_u8(sse_128_loadlo_vector64(pred_ln1),sse_128_loadlo_vector64(pred_ln2)))));	

#define CALC_ALIGNED_SAD_2x8(result, src_ln1, src_ln2, pred_ln1, pred_ln2)	\
	result = sse_128_add_i64(result,sse_128_sad_u8(sse128_unpacklo_u64(sse_128_load_vector_u(src_ln1),sse_128_load_vector_u(src_ln2)),sse128_packs_i16_u8(sse_128_load_vector_u(pred_ln1),sse_128_load_vector_u(pred_ln2))));
//	result = sse_128_add_i64(result, sse_128_sad_u8(sse128_unpacklo_u8(sse_128_load_vector_u(src_ln1),sse_128_load_vector_u(src_ln2)),sse128_unpacklo_u8(sse_128_load_vector_u(pred_ln1),sse_128_load_vector_u(pred_ln2))));

#define CALC_ALIGNED_SAD_16(result, src, pred)																								\
	result = sse_128_add_i64(result, sse_128_sad_u8(sse_128_load_vector_a(src),sse128_packs_i16_u8(sse_128_load_vector_u(pred),sse_128_load_vector_u(pred+8))));
//	result = sse_128_add_i64(result, sse_128_sad_u8(sse_128_load_vector_a(src),sse_128_load_vector_a(pred)));						

#define CALC_ALIGNED_SAD_32(result, src, pred)																								\
	CALC_ALIGNED_SAD_16(result, src, pred)																									\
	CALC_ALIGNED_SAD_16(result, src+16, pred+16)																							

#define CALC_ALIGNED_SAD_64(result, src, pred)																								\
	CALC_ALIGNED_SAD_16(result, src, pred)																									\
	CALC_ALIGNED_SAD_16(result, src+16, pred+16)																							\
	CALC_ALIGNED_SAD_16(result, src+32, pred+32)																							\
	CALC_ALIGNED_SAD_16(result, src+48, pred+48)																							


//--------------------------------------- SSD -----------------------------------------------------------------------------------------------------------
#define CALC_ALIGNED_SSD_2x4(result, src_ln1, src_ln2, pred_ln1, pred_ln2, aux)																				\
	aux = sse_128_sub_i16(sse_128_convert_u8_i16(sse128_unpacklo_u8(sse_128_load_vector_u(src_ln1),sse_128_load_vector_u(src_ln2))),sse128_unpacklo_u16(sse_128_load_vector_u(pred_ln1),sse_128_load_vector_u(pred_ln2)));									\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));

#define CALC_ALIGNED_SSD_2x8(result, src_ln1, src_ln2, pred_ln1, pred_ln2, aux)																						\
	aux = sse_128_sub_i16(sse_128_convert_u8_i16(sse_128_load_vector_u(src_ln1)),sse_128_load_vector_u(pred_ln1));													\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));																								\
	aux = sse_128_sub_i16(sse_128_convert_u8_i16(sse_128_load_vector_u(src_ln2)),sse_128_load_vector_u(pred_ln2));													\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));																														

#define CALC_ALIGNED_SSD_16(result, src, pred, aux)																													\
	aux = sse_128_sub_i16(sse_128_convert_u8_i16(sse_128_load_vector_a(src)),sse_128_load_vector_a(pred));															\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));																								\
	aux = sse_128_sub_i16(sse_128_convert_u8_i16(sse_128_load_vector_u(src+8)),sse_128_load_vector_a(pred+8));														\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));

#define CALC_ALIGNED_SSD_32(result, src, pred, aux)																													\
	CALC_ALIGNED_SSD_16(result, src, pred,  aux)																													\
	CALC_ALIGNED_SSD_16(result, (src+16), (pred+16), aux)										


#define CALC_ALIGNED_SSD_64(result, src, pred, aux)																													\
	CALC_ALIGNED_SSD_16(result, src, pred, aux)																														\
	CALC_ALIGNED_SSD_16(result, src+16, pred+16, aux)																												\
	CALC_ALIGNED_SSD_16(result, src+32, pred+32, aux)																												\
	CALC_ALIGNED_SSD_16(result, src+48, pred+48, aux)																				


//--------------------------------------- SSD16b -----------------------------------------------------------------------------------------------------------

#define CALC_ALIGNED_SSD16b_2x4(result, src_ln1, src_ln2, pred_ln1, pred_ln2, aux)																				\
	aux = sse_128_sub_i16(sse128_unpacklo_u16(sse_128_load_vector_u(src_ln1),sse_128_load_vector_u(src_ln2)),sse128_unpacklo_u16(sse_128_load_vector_u(pred_ln1),sse_128_load_vector_u(pred_ln2)));									\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));


#define CALC_ALIGNED_SSD16b_8(result, src, pred, aux)																												\
	aux = sse_128_sub_i16(sse_128_load_vector_u(src),sse_128_load_vector_u(pred));																					\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));

#define CALC_ALIGNED_SSD16b_16(result, src, pred, aux)																												\
	CALC_ALIGNED_SSD16b_8(result, src, pred,  aux)																													\
	CALC_ALIGNED_SSD16b_8(result, (src+8), (pred+8), aux)

#define CALC_ALIGNED_SSD16b_32(result, src, pred, aux)																												\
	CALC_ALIGNED_SSD16b_16(result, src, pred,  aux)																													\
	CALC_ALIGNED_SSD16b_16(result, (src+16), (pred+16), aux)										

#define CALC_ALIGNED_SSD16b_64(result, src, pred, aux)																												\
	CALC_ALIGNED_SSD16b_16(result, src, pred, aux)																													\
	CALC_ALIGNED_SSD16b_16(result, src+16, pred+16, aux)																											\
	CALC_ALIGNED_SSD16b_16(result, src+32, pred+32, aux)																											\
	CALC_ALIGNED_SSD16b_16(result, src+48, pred+48, aux)																				



//--------------------------------------- PREDICT -----------------------------------------------------------------------------------------------------------

#define CALC_ALIGNED_PREDICT_4(src, pred, dst, zero)																														\
	sse_64_storel_vector_u((dst), sse_128_sub_i16(sse128_unpacklo_u8(sse_128_load_vector_u(src),(zero)),sse_128_load_vector_u(pred)));	


#define CALC_ALIGNED_PREDICT_8(src, pred, dst, zero)																														\
	sse_128_store_vector_u((dst), sse_128_sub_i16(sse128_unpacklo_u8(sse_128_load_vector_u(src),zero),sse_128_load_vector_u(pred)));					

#define CALC_ALIGNED_PREDICT_16(src, pred, dst, zero)																														\
	sse_128_store_vector_a((dst), sse_128_sub_i16(sse128_unpacklo_u8(sse_128_load_vector_a(src),zero),sse_128_load_vector_a(pred)));										\
	sse_128_store_vector_a((dst+8), sse_128_sub_i16(sse128_unpackhi_u8(sse_128_load_vector_a(src),zero),sse_128_load_vector_a(pred+8)));	


#define CALC_ALIGNED_PREDICT_32(src, pred, dst, zero)													\
	CALC_ALIGNED_PREDICT_16(src, pred, dst, zero)														\
	CALC_ALIGNED_PREDICT_16(src+16, pred+16, dst+16, zero)										


#define CALC_ALIGNED_PREDICT_64(src, pred, dst, zero)															\
	CALC_ALIGNED_PREDICT_16(src, pred, dst, zero)															\
	CALC_ALIGNED_PREDICT_16(src+16, pred+16, dst+16, zero)													\
	CALC_ALIGNED_PREDICT_16(src+32, pred+32, dst+32, zero)													\
	CALC_ALIGNED_PREDICT_16(src+48, pred+48, dst+48, zero)


//--------------------------------------- RECONST -----------------------------------------------------------------------------------------------------------

#define CALC_ALIGNED_RECONST_4(pred, resi, deco, zero)																															\
	sse_32_store_vector0_u(deco, sse128_packs_i16_u8(sse_128_adds_i16(sse128_unpacklo_u8(sse_128_load_vector_u(pred),zero), sse_128_load_vector_u(resi)),zero));

#define CALC_ALIGNED_RECONST_8(pred, resi, deco, zero)																															\
	sse_128_store_vector_u(deco, sse_128_convert_u8_i16(sse128_packs_i16_u8(sse_128_adds_i16(sse_128_load_vector_u(pred), sse_128_load_vector_u(resi)),zero)));
	//sse_64_storel_vector_u(deco, sse128_packs_i16_u8(sse_128_adds_i16(sse128_unpacklo_u8(sse_128_load_vector_u(pred),zero), sse_128_load_vector_a(resi)),zero));
	//sse_64_storel_vector_u(deco, _mm_adds_epu8(sse_128_load_vector_u(pred), sse128_packs_i16_u8(sse_128_load_vector_a(resi),sse_128_load_vector_a(resi)))); 						

#define CALC_ALIGNED_RECONST_16(pred, resi, deco, zero)																															\
	CALC_ALIGNED_RECONST_8(pred, resi, deco, zero)																																\
	CALC_ALIGNED_RECONST_8(pred+8, resi+8, deco+8, zero)
	//sse_128_store_vector_a(deco, sse128_packs_i16_u8(sse_128_adds_i16(sse_128_load_vector_a(resi), sse128_unpacklo_u8(sse_128_load_vector_a(pred), zero)),sse_128_adds_i16(sse_128_load_vector_a(resi+8), sse128_unpackhi_u8(sse_128_load_vector_a(pred), zero))));

//		sse_128_store_vector_a(deco, _mm_adds_epu8(sse_128_load_vector_a(pred), sse128_packs_i16_u8(sse_128_load_vector_a(resi),sse_128_load_vector_a(resi+8)))); 
//		sse_128_store_vector_a(deco, _mm_adds_epu8(sse_128_load_vector_a(pred), sse128_packs_i16_u8(sse_128_load_vector_a(resi),sse_128_load_vector_a(resi+8)))); 

#define CALC_ALIGNED_RECONST_32(pred, resi, deco, zero)													\
	CALC_ALIGNED_RECONST_16(pred, resi, deco, zero)														\
	CALC_ALIGNED_RECONST_16(pred+16, resi+16, deco+16, zero)										


#define CALC_ALIGNED_RECONST_64(pred, resi, deco, zero)														\
	CALC_ALIGNED_RECONST_16(pred, resi, deco, zero)															\
	CALC_ALIGNED_RECONST_16(pred+16, resi+16, deco+16, zero)													\
	CALC_ALIGNED_RECONST_16(pred+32, resi+32, deco+32, zero)													\
	CALC_ALIGNED_RECONST_16(pred+48, resi+48, deco+48, zero)													


#define TRANSPOSE_MATRIX_4x4_8BITS(regs, matrix, stride)												\
{																										\
	__m128_i16 l0l1_l = sse128_unpacklo_u8(regs[0],regs[1]);											\
	__m128_i16 l2l3_l = sse128_unpacklo_u8(regs[2],regs[3]);											\
	__m128_i16 c0c1c2c3 = sse128_unpacklo_u16(l0l1_l,l2l3_l);											\
	sse_64_storel_vector_u(matrix, c0c1c2c3);															\
	sse_64_storeh_vector_u(matrix+2*stride, c0c1c2c3);													\
	c0c1c2c3 =  sse_128_shift_r_u32(c0c1c2c3,32);														\
	sse_64_storel_vector_u(matrix+1*stride, c0c1c2c3);													\
	sse_64_storeh_vector_u(matrix+3*stride, c0c1c2c3);													\
}


#define TRANSPOSE_MATRIX_8x8_8BITS(regs, matrix, stride)												\
{																										\
	__m128_i16 l0l1_l = sse128_unpacklo_u8(regs[0],regs[1]);											\
	__m128_i16 l2l3_l = sse128_unpacklo_u8(regs[2],regs[3]);											\
	__m128_i16 l4l5_l = sse128_unpacklo_u8(regs[4],regs[5]);											\
	__m128_i16 l6l7_l = sse128_unpacklo_u8(regs[6],regs[7]);											\
																										\
	__m128_i16 l0l1l2l3_ll = sse128_unpacklo_u16(l0l1_l,l2l3_l);										\
	__m128_i16 l0l1l2l3_lh = sse128_unpackhi_u16(l0l1_l,l2l3_l);										\
	__m128_i16 l4l5l6l7_ll = sse128_unpacklo_u16(l4l5_l,l6l7_l);										\
	__m128_i16 l4l5l6l7_lh = sse128_unpackhi_u16(l4l5_l,l6l7_l);										\
																										\
	__m128_i16 c0c1_l = sse128_unpacklo_u32(l0l1l2l3_ll,l4l5l6l7_ll);									\
	__m128_i16 c2c3_l = sse128_unpackhi_u32(l0l1l2l3_ll,l4l5l6l7_ll);									\
	__m128_i16 c4c5_l = sse128_unpacklo_u32(l0l1l2l3_lh,l4l5l6l7_lh);									\
	__m128_i16 c6c7_l = sse128_unpackhi_u32(l0l1l2l3_lh,l4l5l6l7_lh);									\
																										\
	sse_64_storel_vector_u(matrix, c0c1_l);																\
	sse_64_storeh_vector_u(matrix+1*stride, c0c1_l);													\
	sse_64_storel_vector_u(matrix+2*stride, c2c3_l);													\
	sse_64_storeh_vector_u(matrix+3*stride, c2c3_l);													\
	sse_64_storel_vector_u(matrix+4*stride, c4c5_l);													\
	sse_64_storeh_vector_u(matrix+5*stride, c4c5_l);													\
	sse_64_storel_vector_u(matrix+6*stride, c6c7_l);													\
	sse_64_storeh_vector_u(matrix+7*stride, c6c7_l);													\
}

//if the matrix is bigger than 8x8, data is transposed inside the 16x16 blocks, but these 16x16 blocks have to be already transposed
#define TRANSPOSE_MATRIX_8BIT(matrix, matrix_size, stride)															\
{																													\
	int i, j;																										\
	if(matrix_size>8)																								\
	{																												\
		for (j=0;j<matrix_size;j+=16)																				\
		{																											\
			for (i=0;i<matrix_size;i+=16)																			\
			{																										\
				__m128_i16 c0c1_h, c2c3_h, c4c5_h, c6c7_h, c8c9_h, c10c11_h, c12c13_h, c14c15_h;					\
																													\
				__m128_i16 l0 = sse_128_load_vector_u(matrix+j*stride+i);											\
				__m128_i16 l1 = sse_128_load_vector_u(matrix+(j+1)*stride+i);										\
				__m128_i16 l2 = sse_128_load_vector_u(matrix+(j+2)*stride+i);										\
				__m128_i16 l3 = sse_128_load_vector_u(matrix+(j+3)*stride+i);										\
				__m128_i16 l4 = sse_128_load_vector_u(matrix+(j+4)*stride+i);										\
				__m128_i16 l5 = sse_128_load_vector_u(matrix+(j+5)*stride+i);										\
				__m128_i16 l6 = sse_128_load_vector_u(matrix+(j+6)*stride+i);										\
				__m128_i16 l7 = sse_128_load_vector_u(matrix+(j+7)*stride+i);										\
																													\
				__m128_i16 l0l1_l = sse128_unpacklo_u8(l0,l1);														\
				__m128_i16 l0l1_h = sse128_unpackhi_u8(l0,l1);														\
				__m128_i16 l2l3_l = sse128_unpacklo_u8(l2,l3);														\
				__m128_i16 l2l3_h = sse128_unpackhi_u8(l2,l3);														\
				__m128_i16 l4l5_l = sse128_unpacklo_u8(l4,l5);														\
				__m128_i16 l4l5_h = sse128_unpackhi_u8(l4,l5);														\
				__m128_i16 l6l7_l = sse128_unpacklo_u8(l6,l7);														\
				__m128_i16 l6l7_h = sse128_unpackhi_u8(l6,l7);														\
																													\
				__m128_i16 l0l1l2l3_ll = sse128_unpacklo_u16(l0l1_l,l2l3_l);										\
				__m128_i16 l0l1l2l3_lh = sse128_unpackhi_u16(l0l1_l,l2l3_l);										\
				__m128_i16 l0l1l2l3_hl = sse128_unpacklo_u16(l0l1_h,l2l3_h);										\
				__m128_i16 l0l1l2l3_hh = sse128_unpackhi_u16(l0l1_h,l2l3_h);										\
				__m128_i16 l4l5l6l7_ll = sse128_unpacklo_u16(l4l5_l,l6l7_l);										\
				__m128_i16 l4l5l6l7_lh = sse128_unpackhi_u16(l4l5_l,l6l7_l);										\
				__m128_i16 l4l5l6l7_hl = sse128_unpacklo_u16(l4l5_h,l6l7_h);										\
				__m128_i16 l4l5l6l7_hh = sse128_unpackhi_u16(l4l5_h,l6l7_h);										\
																													\
				__m128_i16 c0c1_l = sse128_unpacklo_u32(l0l1l2l3_ll,l4l5l6l7_ll);									\
				__m128_i16 c2c3_l = sse128_unpackhi_u32(l0l1l2l3_ll,l4l5l6l7_ll);									\
				__m128_i16 c4c5_l = sse128_unpacklo_u32(l0l1l2l3_lh,l4l5l6l7_lh);									\
				__m128_i16 c6c7_l = sse128_unpackhi_u32(l0l1l2l3_lh,l4l5l6l7_lh);									\
				__m128_i16 c8c9_l = sse128_unpacklo_u32(l0l1l2l3_hl,l4l5l6l7_hl);									\
				__m128_i16 c10c11_l = sse128_unpackhi_u32(l0l1l2l3_hl,l4l5l6l7_hl);									\
				__m128_i16 c12c13_l = sse128_unpacklo_u32(l0l1l2l3_hh,l4l5l6l7_hh);									\
				__m128_i16 c14c15_l = sse128_unpackhi_u32(l0l1l2l3_hh,l4l5l6l7_hh);									\
																													\
				l0 = sse_128_load_vector_u(matrix+(j+8)*stride+i);													\
				l1 = sse_128_load_vector_u(matrix+(j+9)*stride+i);													\
				l2 = sse_128_load_vector_u(matrix+(j+10)*stride+i);													\
				l3 = sse_128_load_vector_u(matrix+(j+11)*stride+i);													\
				l4 = sse_128_load_vector_u(matrix+(j+12)*stride+i);													\
				l5 = sse_128_load_vector_u(matrix+(j+13)*stride+i);													\
				l6 = sse_128_load_vector_u(matrix+(j+14)*stride+i);													\
				l7 = sse_128_load_vector_u(matrix+(j+15)*stride+i);													\
																													\
				l0l1_l = sse128_unpacklo_u8(l0,l1);																	\
				l0l1_h = sse128_unpackhi_u8(l0,l1);																	\
				l2l3_l = sse128_unpacklo_u8(l2,l3);																	\
				l2l3_h = sse128_unpackhi_u8(l2,l3);																	\
				l4l5_l = sse128_unpacklo_u8(l4,l5);																	\
				l4l5_h = sse128_unpackhi_u8(l4,l5);																	\
				l6l7_l = sse128_unpacklo_u8(l6,l7);																	\
				l6l7_h = sse128_unpackhi_u8(l6,l7);																	\
																													\
				l0l1l2l3_ll = sse128_unpacklo_u16(l0l1_l,l2l3_l);													\
				l0l1l2l3_lh = sse128_unpackhi_u16(l0l1_l,l2l3_l);													\
				l0l1l2l3_hl = sse128_unpacklo_u16(l0l1_h,l2l3_h);													\
				l0l1l2l3_hh = sse128_unpackhi_u16(l0l1_h,l2l3_h);													\
				l4l5l6l7_ll = sse128_unpacklo_u16(l4l5_l,l6l7_l);													\
				l4l5l6l7_lh = sse128_unpackhi_u16(l4l5_l,l6l7_l);													\
				l4l5l6l7_hl = sse128_unpacklo_u16(l4l5_h,l6l7_h);													\
				l4l5l6l7_hh = sse128_unpackhi_u16(l4l5_h,l6l7_h);													\
																													\
				c0c1_h = sse128_unpacklo_u32(l0l1l2l3_ll,l4l5l6l7_ll);												\
				c2c3_h = sse128_unpackhi_u32(l0l1l2l3_ll,l4l5l6l7_ll);												\
				c4c5_h = sse128_unpacklo_u32(l0l1l2l3_lh,l4l5l6l7_lh);												\
				c6c7_h = sse128_unpackhi_u32(l0l1l2l3_lh,l4l5l6l7_lh);												\
				c8c9_h = sse128_unpacklo_u32(l0l1l2l3_hl,l4l5l6l7_hl);												\
				c10c11_h = sse128_unpackhi_u32(l0l1l2l3_hl,l4l5l6l7_hl);											\
				c12c13_h = sse128_unpacklo_u32(l0l1l2l3_hh,l4l5l6l7_hh);											\
				c14c15_h = sse128_unpackhi_u32(l0l1l2l3_hh,l4l5l6l7_hh);											\
																													\
				sse_128_store_vector_u(matrix+(j)*stride+i, sse128_unpacklo_u64(c0c1_l,c0c1_h));					\
				sse_128_store_vector_u(matrix+(j+1)*stride+i, sse128_unpackhi_u64(c0c1_l,c0c1_h));					\
				sse_128_store_vector_u(matrix+(j+2)*stride+i, sse128_unpacklo_u64(c2c3_l,c2c3_h));					\
				sse_128_store_vector_u(matrix+(j+3)*stride+i, sse128_unpackhi_u64(c2c3_l,c2c3_h));					\
				sse_128_store_vector_u(matrix+(j+4)*stride+i, sse128_unpacklo_u64(c4c5_l,c4c5_h));					\
				sse_128_store_vector_u(matrix+(j+5)*stride+i, sse128_unpackhi_u64(c4c5_l,c4c5_h));					\
				sse_128_store_vector_u(matrix+(j+6)*stride+i, sse128_unpacklo_u64(c6c7_l,c6c7_h));					\
				sse_128_store_vector_u(matrix+(j+7)*stride+i, sse128_unpackhi_u64(c6c7_l,c6c7_h));					\
				sse_128_store_vector_u(matrix+(j+8)*stride+i, sse128_unpacklo_u64(c8c9_l,c8c9_h));					\
				sse_128_store_vector_u(matrix+(j+9)*stride+i, sse128_unpackhi_u64(c8c9_l,c8c9_h));					\
				sse_128_store_vector_u(matrix+(j+10)*stride+i, sse128_unpacklo_u64(c10c11_l,c10c11_h));				\
				sse_128_store_vector_u(matrix+(j+11)*stride+i, sse128_unpackhi_u64(c10c11_l,c10c11_h));				\
				sse_128_store_vector_u(matrix+(j+12)*stride+i, sse128_unpacklo_u64(c12c13_l,c12c13_h));				\
				sse_128_store_vector_u(matrix+(j+13)*stride+i, sse128_unpackhi_u64(c12c13_l,c12c13_h));				\
				sse_128_store_vector_u(matrix+(j+14)*stride+i, sse128_unpacklo_u64(c14c15_l,c14c15_h));				\
				sse_128_store_vector_u(matrix+(j+15)*stride+i, sse128_unpackhi_u64(c14c15_l,c14c15_h));				\
			}																										\
		}																											\
	}																												\
	else if(matrix_size==8)																							\
	{																												\
		__m128_i16 l[8];																							\
		l[0] = sse_128_load_vector_u(matrix);																		\
		l[1] = sse_128_load_vector_u(matrix+1*stride);																\
		l[2] = sse_128_load_vector_u(matrix+2*stride);																\
		l[3] = sse_128_load_vector_u(matrix+3*stride);																\
		l[4] = sse_128_load_vector_u(matrix+4*stride);																\
		l[5] = sse_128_load_vector_u(matrix+5*stride);																\
		l[6] = sse_128_load_vector_u(matrix+6*stride);																\
		l[7] = sse_128_load_vector_u(matrix+7*stride);																\
		TRANSPOSE_MATRIX_8x8_8BITS(l,matrix,stride)																	\
	}																												\
	else if(matrix_size==4)																							\
	{																												\
		__m128_i16 l[4];																							\
		l[0] = sse_128_load_vector_u(matrix);																		\
		l[1] = sse_128_load_vector_u(matrix+1*stride);																\
		l[2] = sse_128_load_vector_u(matrix+2*stride);																\
		l[3] = sse_128_load_vector_u(matrix+3*stride);																\
		TRANSPOSE_MATRIX_4x4_8BITS(l, matrix, stride)																\
	}																												\
}



#define TRANSPOSE_MATRIX_4x4_16BITS(regs, matrix, stride)															\
{																													\
	__m128_i16 l0l1_l = sse128_unpacklo_u16(regs[0],regs[1]);														\
	__m128_i16 l2l3_l = sse128_unpacklo_u16(regs[2],regs[3]);														\
	__m128_i16 c0c1c2c3 = sse128_unpacklo_u32(l0l1_l,l2l3_l);														\
	sse_64_storel_vector_u(matrix, c0c1c2c3);																		\
	sse_64_storeh_vector_u(matrix+1*stride, c0c1c2c3);																\
	c0c1c2c3 = sse128_unpackhi_u32(l0l1_l,l2l3_l);																	\
	sse_64_storel_vector_u(matrix+2*stride, c0c1c2c3);																\
	sse_64_storeh_vector_u(matrix+3*stride, c0c1c2c3);																\
}


#define TRANSPOSE_MATRIX_8x8_16BITS(regs, matrix, stride)												\
{																										\
	__m128_i16 l0l1_l = sse128_unpacklo_u16(regs[0],regs[1]);											\
	__m128_i16 l0l1_h = sse128_unpackhi_u16(regs[0],regs[1]);											\
	__m128_i16 l2l3_l = sse128_unpacklo_u16(regs[2],regs[3]);											\
	__m128_i16 l2l3_h = sse128_unpackhi_u16(regs[2],regs[3]);											\
	__m128_i16 l4l5_l = sse128_unpacklo_u16(regs[4],regs[5]);											\
	__m128_i16 l4l5_h = sse128_unpackhi_u16(regs[4],regs[5]);											\
	__m128_i16 l6l7_l = sse128_unpacklo_u16(regs[6],regs[7]);											\
	__m128_i16 l6l7_h = sse128_unpackhi_u16(regs[6],regs[7]);											\
																										\
	__m128_i16 l0l1l2l3_ll = sse128_unpacklo_u32(l0l1_l,l2l3_l);										\
	__m128_i16 l0l1l2l3_lh = sse128_unpackhi_u32(l0l1_l,l2l3_l);										\
	__m128_i16 l0l1l2l3_hl = sse128_unpacklo_u32(l0l1_h,l2l3_h);										\
	__m128_i16 l0l1l2l3_hh = sse128_unpackhi_u32(l0l1_h,l2l3_h);										\
	__m128_i16 l4l5l6l7_ll = sse128_unpacklo_u32(l4l5_l,l6l7_l);										\
	__m128_i16 l4l5l6l7_lh = sse128_unpackhi_u32(l4l5_l,l6l7_l);										\
	__m128_i16 l4l5l6l7_hl = sse128_unpacklo_u32(l4l5_h,l6l7_h);										\
	__m128_i16 l4l5l6l7_hh = sse128_unpackhi_u32(l4l5_h,l6l7_h);										\
																										\
	__m128_i16 c0 = sse128_unpacklo_u64(l0l1l2l3_ll,l4l5l6l7_ll);										\
	__m128_i16 c1 = sse128_unpackhi_u64(l0l1l2l3_ll,l4l5l6l7_ll);										\
	__m128_i16 c2 = sse128_unpacklo_u64(l0l1l2l3_lh,l4l5l6l7_lh);										\
	__m128_i16 c3 = sse128_unpackhi_u64(l0l1l2l3_lh,l4l5l6l7_lh);										\
	__m128_i16 c4 = sse128_unpacklo_u64(l0l1l2l3_hl,l4l5l6l7_hl);										\
	__m128_i16 c5 = sse128_unpackhi_u64(l0l1l2l3_hl,l4l5l6l7_hl);										\
	__m128_i16 c6 = sse128_unpacklo_u64(l0l1l2l3_hh,l4l5l6l7_hh);										\
	__m128_i16 c7 = sse128_unpackhi_u64(l0l1l2l3_hh,l4l5l6l7_hh);										\
																										\
	sse_128_store_vector_u(matrix, c0);																	\
	sse_128_store_vector_u(matrix+1*stride, c1);														\
	sse_128_store_vector_u(matrix+2*stride, c2);														\
	sse_128_store_vector_u(matrix+3*stride, c3);														\
	sse_128_store_vector_u(matrix+4*stride, c4);														\
	sse_128_store_vector_u(matrix+5*stride, c5);														\
	sse_128_store_vector_u(matrix+6*stride, c6);														\
	sse_128_store_vector_u(matrix+7*stride, c7);														\
}


//Traspositions
//Matrices bigger than 8x8 are transposed in 8x8 blocks. First the diagonal, then the rest in groups of two 8x8 blocks
#define TRANSPOSE_MATRIX_16BIT(matrix, matrix_size, stride)												\
{																										\
	if(matrix_size==4)																					\
	{																									\
		__m128_i16 l[4];																				\
		l[0] = sse_128_load_vector_u(matrix);															\
		l[1] = sse_128_load_vector_u(matrix+1*stride);													\
		l[2] = sse_128_load_vector_u(matrix+2*stride);													\
		l[3] = sse_128_load_vector_u(matrix+3*stride);													\
																										\
		TRANSPOSE_MATRIX_4x4_16BITS(l, matrix, stride)													\
	}																									\
	else if(matrix_size==8)																				\
	{																									\
		int i, j;																						\
		for (j=0;j<matrix_size;j+=8)																	\
		{																								\
			for (i=0;i<matrix_size;i+=8)																\
			{																							\
				__m128_i16 l[8];																		\
				l[0] = sse_128_load_vector_u(matrix+(j+0)*stride+i);									\
				l[1] = sse_128_load_vector_u(matrix+(j+1)*stride+i);									\
				l[2] = sse_128_load_vector_u(matrix+(j+2)*stride+i);									\
				l[3] = sse_128_load_vector_u(matrix+(j+3)*stride+i);									\
				l[4] = sse_128_load_vector_u(matrix+(j+4)*stride+i);									\
				l[5] = sse_128_load_vector_u(matrix+(j+5)*stride+i);									\
				l[6] = sse_128_load_vector_u(matrix+(j+6)*stride+i);									\
				l[7] = sse_128_load_vector_u(matrix+(j+7)*stride+i);									\
				TRANSPOSE_MATRIX_8x8_16BITS(l, matrix+(j+0)*stride+i, stride)							\
			}																							\
		}																								\
	}																									\
	else																								\
	{																									\
		int i, j;																				\
		for(i=0;i<matrix_size;i+=8)																	\
		{																							\
			__m128_i16 block8x8a[8];																\
			block8x8a[0] = sse_128_load_vector_u(matrix+(i+0)*stride+i);							\
			block8x8a[1] = sse_128_load_vector_u(matrix+(i+1)*stride+i);							\
			block8x8a[2] = sse_128_load_vector_u(matrix+(i+2)*stride+i);							\
			block8x8a[3] = sse_128_load_vector_u(matrix+(i+3)*stride+i);							\
			block8x8a[4] = sse_128_load_vector_u(matrix+(i+4)*stride+i);							\
			block8x8a[5] = sse_128_load_vector_u(matrix+(i+5)*stride+i);							\
			block8x8a[6] = sse_128_load_vector_u(matrix+(i+6)*stride+i);							\
			block8x8a[7] = sse_128_load_vector_u(matrix+(i+7)*stride+i);							\
			TRANSPOSE_MATRIX_8x8_16BITS(block8x8a, matrix+(i+0)*stride+i, stride)					\
		}																							\
		for (j=0;j<matrix_size;j+=8)																\
		{																							\
			for(i=j+8;i<matrix_size;i+=8)															\
			{																						\
				__m128_i16 block8x8a[8];															\
				__m128_i16 block8x8b[8];															\
				block8x8a[0] = sse_128_load_vector_u(matrix+(j+0)*stride+i);						\
				block8x8a[1] = sse_128_load_vector_u(matrix+(j+1)*stride+i);						\
				block8x8a[2] = sse_128_load_vector_u(matrix+(j+2)*stride+i);						\
				block8x8a[3] = sse_128_load_vector_u(matrix+(j+3)*stride+i);						\
				block8x8a[4] = sse_128_load_vector_u(matrix+(j+4)*stride+i);						\
				block8x8a[5] = sse_128_load_vector_u(matrix+(j+5)*stride+i);						\
				block8x8a[6] = sse_128_load_vector_u(matrix+(j+6)*stride+i);						\
				block8x8a[7] = sse_128_load_vector_u(matrix+(j+7)*stride+i);						\
																									\
				block8x8b[0] = sse_128_load_vector_u(matrix+(i+0)*stride+j);						\
				block8x8b[1] = sse_128_load_vector_u(matrix+(i+1)*stride+j);						\
				block8x8b[2] = sse_128_load_vector_u(matrix+(i+2)*stride+j);						\
				block8x8b[3] = sse_128_load_vector_u(matrix+(i+3)*stride+j);						\
				block8x8b[4] = sse_128_load_vector_u(matrix+(i+4)*stride+j);						\
				block8x8b[5] = sse_128_load_vector_u(matrix+(i+5)*stride+j);						\
				block8x8b[6] = sse_128_load_vector_u(matrix+(i+6)*stride+j);						\
				block8x8b[7] = sse_128_load_vector_u(matrix+(i+7)*stride+j);						\
																									\
				TRANSPOSE_MATRIX_8x8_16BITS(block8x8a, matrix+(i+0)*stride+j, stride)				\
				TRANSPOSE_MATRIX_8x8_16BITS(block8x8b, matrix+(j+0)*stride+i, stride)				\
			}																						\
		}																							\
	}																								\
}

#endif /*__HOMER_HEVC_SSE42_MACROS__*/