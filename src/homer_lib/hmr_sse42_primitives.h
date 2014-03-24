/*****************************************************************************
* hmr_sse42_functions_primitives.c : homerHEVC encoding library
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


#ifndef __HOMER_HEVC_SSE42_PRIMITIVES__
#define __HOMER_HEVC_SSE42_PRIMITIVES__

#include "hmr_os_primitives.h"
#include <smmintrin.h>

typedef __m128i	__m128_;				//generico
typedef __m128i	__m128_u64;				//array de 128 bits distribuido en datos unsigned de 64 bits (2 datos)
typedef __m128i	__m128_i64;				
typedef __m128i	__m128_u32;				//array de 128 bits distribuido en datos unsigned de 32 bits (4 datos)
typedef __m128i	__m128_i32;				
typedef __m128i	__m128_u16;				//array de 128 bits distribuido en datos unsigned de 16 bits (8 datos)
typedef __m128i	__m128_i16;				
typedef __m128i	__m128_u8;				//array de 128 bits distribuido en datos unsigned de 8 bits (16 datos)
typedef __m128i	__m128_i8;				


#define sse_128_load_vector_a(p)		_mm_load_si128((__m128i const*)(p))
#define sse_128_load_vector_u(p)		_mm_loadu_si128((__m128i const*)(p))
#define sse_128_loadlo_vector64(p)		_mm_loadl_epi64((__m128i const*)(p))
//#define sse_128_loadhi_vector64(p)		_mm_loadh_epi64((__m128i const*)p)

#define sse_128_get_data_u64(a,index)	_mm_extract_epi64(a,index)//(a.m128i_u64[index])
#define sse_128_get_data_u32(a,index)	_mm_extract_epi32(a,index)//(a.m128i_u32[index])
#define sse_128_get_data_u16(a,index)	_mm_extract_epi16(a,index)//(a.m128i_u16[index])

#define sse_128_store_vector_a(p, val)	_mm_store_si128((__m128i*)(p),val)
#define sse_128_store_vector_u(p, val)	_mm_storeu_si128((__m128i*)(p),val)
#define sse_64_storel_vector_u(p, val)	_mm_storel_epi64((__m128i*)(p),val)
#define sse_64_storeh_vector_u(p, val)	sse_64_storel_vector_u((p), sse128_unpackhi_u64(val, val))//*((uint64_t*)(p)) = sse_128_get_data_u64(val,1)
#define sse_32_store_vector0_u(p,val)	*((uint32_t*)(p)) = sse_128_get_data_u32(val,0)
#define sse_32_store_vector1_u(p,val)	*((uint32_t*)(p)) = sse_128_get_data_u32(val,1)
#define sse_32_store_vector2_u(p,val)	*((uint32_t*)(p)) = sse_128_get_data_u32(val,2)
#define sse_32_store_vector3_u(p,val)	*((uint32_t*)(p)) = sse_128_get_data_u32(val,3)

#define sse_128_zero_vector()			_mm_setzero_si128()
#define sse_128_vector_i32(a)			_mm_set1_epi32(a)
#define sse_128_vector_i16(a)			_mm_set1_epi16(a)
#define sse_128_vector_i8(a)			_mm_set1_epi8(a)
#define sse_128_set_vector_u8(i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) _mm_setr_epi8(i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15)


#define sse128_packs_i16_u8(a,b)			_mm_packus_epi16(a, b)
//#define sse128_packs_i16_i8(a,b)			_mm_packs_epi16(a, b)
#define sse128_packs_u32_u16(a,b)			_mm_packs_epi32(a, b)//this is not correct but widely used. It should be _mm_packus_epi32. It must be checked before being changed


#define sse128_unpacklo_u8(val1,val2)		_mm_unpacklo_epi8(val1,val2)
#define sse128_unpackhi_u8(val1,val2)		_mm_unpackhi_epi8(val1,val2)
#define sse128_unpacklo_u16(val1,val2)		_mm_unpacklo_epi16(val1,val2)
#define sse128_unpackhi_u16(val1,val2)		_mm_unpackhi_epi16(val1,val2)
#define sse128_unpacklo_u32(val1,val2)		_mm_unpacklo_epi32(val1,val2)
#define sse128_unpackhi_u32(val1,val2)		_mm_unpackhi_epi32(val1,val2)
#define sse128_unpacklo_u64(val1,val2)		_mm_unpacklo_epi64(val1,val2)
#define sse128_unpackhi_u64(val1,val2)		_mm_unpackhi_epi64(val1,val2)

#define sse_128_convert_u8_i16(a)			_mm_cvtepu8_epi16(a)
#define sse_128_convert_i16_i32(a)			_mm_cvtepi16_epi32(a)


#define sse_128_and(a,b)				_mm_and_si128(a,b)

#define sse_128_sign_16(one, a)				_mm_sign_epi16(one,a)



#define sse_128_sub_i16(a,b)				_mm_sub_epi16(a,b)
#define sse_128_sub_i32(a,b)				_mm_sub_epi32(a,b)

#define sse_128_add_i8(a,b)					_mm_add_epi8(a,b)
#define sse_128_add_i16(a,b)				_mm_add_epi16(a,b)
#define sse_128_add_i32(a,b)				_mm_add_epi32(a,b)
#define sse_128_add_i64(a,b)				_mm_add_epi64(a,b)

#define sse_128_adds_i16(a,b)				_mm_adds_epi16(a,b)

#define sse_128_hsub_i16(a,b)				_mm_hsub_epi16(a,b)
#define sse_128_hsub_i32(a,b)				_mm_hsub_epi32(a,b)

#define sse_128_hadd_i16(a,b)				_mm_hadd_epi16(a,b)
#define sse_128_hadd_i32(a,b)				_mm_hadd_epi32(a,b)
#define sse_128_hadd_i64(a,b)				_mm_hadd_epi64(a,b)

#define sse_128_madd_i16_i32(a,b)			_mm_madd_epi16(a,b)

#define sse_128_mul_i16(a,b)				_mm_mullo_epi16(a,b)
#define sse_128_mul_i32(a,b)				_mm_mullo_epi32(a,b)

#define sse_128_abs_i16(a)					_mm_abs_epi16(a)
#define sse_128_sad_u8(a,b)					_mm_sad_epu8(a,b)

#define sse_128_shift_r_u32(a,shift)		_mm_srli_epi64(a,shift)
#define sse_128_shift_r_i32(a,shift)		_mm_srai_epi32(a,shift)
#define sse_128_shift_r_i16(a,shift)		_mm_srai_epi16(a,shift)
#define sse_128_shift_l_i32(a,shift)		_mm_slli_epi32(a,shift)
#define sse_128_shift_l_i16(a,shift)		_mm_slli_epi16(a,shift)

#define sse_128_shuffle_32(a,b)				_mm_shuffle_epi32(a,b)
#define sse_128_shuffle_8(a,b)				_mm_shuffle_epi8(a,b)


#endif	/*__HOMER_HEVC_SSE42_PRIMITIVES__*/
