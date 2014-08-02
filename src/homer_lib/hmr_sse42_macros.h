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
#define CALC_ALIGNED_SSD_2x4(result, src_ln1, src_ln2, pred_ln1, pred_ln2, zero, aux)																				\
	aux = sse_128_sub_i16(sse128_unpacklo_u8(sse128_unpacklo_u8(sse_128_load_vector_u(src_ln1),sse_128_load_vector_u(src_ln2)),zero),sse128_unpacklo_u8(sse128_unpacklo_u8(sse_128_load_vector_u(pred_ln1),sse_128_load_vector_u(pred_ln2)),zero));									\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));


#define CALC_ALIGNED_SSD_2x8(result, src_ln1, src_ln2, pred_ln1, pred_ln2, zero, aux)																				\
	aux = sse_128_sub_i16(sse128_unpacklo_u8(sse_128_load_vector_u(src_ln1),zero),sse128_unpacklo_u8(sse_128_load_vector_u(pred_ln1),zero));					\
	aux = sse_128_madd_i16_i32(aux, aux);																														\
	result = sse_128_add_i32(result,aux);																													\
	aux = sse_128_sub_i16(sse128_unpacklo_u8(sse_128_load_vector_u(src_ln2),zero),sse128_unpacklo_u8(sse_128_load_vector_u(pred_ln2),zero));					\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));																														

#define CALC_ALIGNED_SSD_16(result, src, pred, zero, aux)																											\
	aux = sse_128_sub_i16(sse128_unpacklo_u8(sse_128_load_vector_a(src),zero),sse128_unpacklo_u8(sse_128_load_vector_a(pred),zero));							\
	aux = sse_128_madd_i16_i32(aux, aux);																														\
	result = sse_128_add_i32(result,aux);																													\
	aux = sse_128_sub_i16(sse128_unpackhi_u8(sse_128_load_vector_a(src),zero),sse128_unpackhi_u8(sse_128_load_vector_a(pred),zero));							\
	result = sse_128_add_i32(result,sse_128_madd_i16_i32(aux, aux));

#define CALC_ALIGNED_SSD_32(result, src, pred, zero, aux)											\
	CALC_ALIGNED_SSD_16(result, src, pred, zero, aux)												\
	CALC_ALIGNED_SSD_16(result, (src+16), (pred+16), zero, aux)										


#define CALC_ALIGNED_SSD_64(result, src, pred, zero, aux)								\
	CALC_ALIGNED_SSD_16(result, src, pred, zero, aux)									\
	CALC_ALIGNED_SSD_16(result, src+16, pred+16, zero, aux)								\
	CALC_ALIGNED_SSD_16(result, src+32, pred+32, zero, aux)								\
	CALC_ALIGNED_SSD_16(result, src+48, pred+48, zero, aux)																				


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


#endif /*__HOMER_HEVC_SSE42_MACROS__*/