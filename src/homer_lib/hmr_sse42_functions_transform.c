/*****************************************************************************
* hmr_sse42_functions_transform.c : homerHEVC encoding library
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
#include "hmr_common.h"

#include "hmr_sse42_primitives.h"
#include "hmr_sse42_macros.h"
#include "hmr_sse42_functions.h"


//----------------------16 bits shuffle masks-----------------------------------
ALIGN(16) static const int8_t shuffle_mask_dct_16_0[16] ={ 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1};//0,1,2,3,4,5,6,7 -> 7,6,5,4,3,2,1,0
ALIGN(16) static const int8_t shuffle_mask_dct_16_1[16] ={ 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9};//0,1,2,3,4,5,6,7 -> 3,2,1,0,7,6,5,4
ALIGN(16) static const int8_t shuffle_mask_dct_16_2[16] ={ 2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13};//0,1,2,3,4,5,6,7 -> 1,0,3,2,5,4,7,6
ALIGN(16) static const int8_t shuffle_mask_dct_16_3[16] ={ 0, 1, 14, 15, 2, 3, 12, 13, 4, 5, 10, 11, 6, 7, 8, 9};//0,1,2,3,4,5,6,7 -> 0,7,1,6,2,5,3,4



ALIGN(16) static const int8_t shuffle_mask_idct16_0[16] ={ 0, 1, 2, 3, 4, 5, 6, 7, 14, 15, 12, 13, 10, 11, 8, 9};//0,1,2,3,4,5,6,7 -> //0,1,2,3,7,6,5,4

//----------------------32 bits shuffle masks-----------------------------------
ALIGN(16) static const int8_t shuffle_mask_dct_32_0[16] ={ 12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3};//0,1,2,3 -> 3,2,1,0
ALIGN(16) static const int8_t shuffle_mask_dct_32_1[16] ={ 4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11};//0,1,2,3 -> 1,0,3,2



//---------------------------------------DST-4x4----------------------------------------------



ALIGN(16) static const int16_t dst4x4_coeffs[4][8] =
{
    { 29, 55, 74, 84, 29, 55, 74, 84 },
    { 74, 74,  0,-74, 74, 74, 0, -74 },
    { 84,-29,-74, 55, 84,-29,-74, 55 },
    { 55,-84, 74,-29, 55,-84, 74,-29 }
};

ALIGN(16) static const int16_t inv_dst4x4_coeffs[4][8] =
{
    { 29, 74, 84, 55, 29, 74, 84, 55 },
    { 55, 74,-29,-84, 55, 74,-29,-84 },
    { 74, 0, -74, 74, 74, 0, -74, 74 },
    { 84,-74, 55,-29, 84,-74, 55,-29 }
};



#define DST4x4(result, round, shift, coeffs)																	\
		__m128i c0 = sse_128_madd_i16_i32(result[0], sse_128_load_vector_a(coeffs[0]));							\
		__m128i c1 = sse_128_madd_i16_i32(result[1], sse_128_load_vector_a(coeffs[0]));							\
		__m128i c2 = sse_128_madd_i16_i32(result[0], sse_128_load_vector_a(coeffs[1]));							\
		__m128i c3 = sse_128_madd_i16_i32(result[1], sse_128_load_vector_a(coeffs[1]));							\
		__m128i c4 = sse_128_madd_i16_i32(result[0], sse_128_load_vector_a(coeffs[2]));							\
		__m128i c5 = sse_128_madd_i16_i32(result[1], sse_128_load_vector_a(coeffs[2]));							\
		__m128i c6 = sse_128_madd_i16_i32(result[0], sse_128_load_vector_a(coeffs[3]));							\
		__m128i c7 = sse_128_madd_i16_i32(result[1], sse_128_load_vector_a(coeffs[3]));							\
		__m128i sum01 = sse_128_hadd_i32(c0, c1);																\
		__m128i sum23 = sse_128_hadd_i32(c2, c3);																\
		__m128i sum45 = sse_128_hadd_i32(c4, c5);																\
		__m128i sum67 = sse_128_hadd_i32(c6, c7);																\
		sum01 = sse_128_shift_r_i32(sse_128_add_i32(sum01, round), shift);										\
		sum23 = sse_128_shift_r_i32(sse_128_add_i32(sum23, round), shift);										\
		sum45 = sse_128_shift_r_i32(sse_128_add_i32(sum45, round), shift);										\
		sum67 = sse_128_shift_r_i32(sse_128_add_i32(sum67, round), shift);										\
		result[0] = sse128_packs_u32_u16(sum01, sum23);															\
		result[1] = sse128_packs_u32_u16(sum45, sum67);




void sse_aligned_dst_4x4(int16_t *src, int16_t *dst, int stride) 
{
	int i, shift[2] = {1,8};
	__m128i result[2];
	__m128i round[2];// = {sse_128_vector_i32(1),sse_128_vector_i32(128)};
//	__m128i final;
	__m128i l0  = sse_128_loadlo_vector64((src));
	__m128i l1  = sse_128_loadlo_vector64((src+stride));
	__m128i l2  = sse_128_loadlo_vector64((src+2*stride));
	__m128i l3  = sse_128_loadlo_vector64((src+3*stride));

	round[0] = sse_128_vector_i32(1);
	round[1] = sse_128_vector_i32(128);

	result[0]  = sse128_unpacklo_u64(l0, l1);
	result[1]  = sse128_unpacklo_u64(l2, l3);

	for (i = 0; i < 2; i++)
	{
		DST4x4(result, round[i], shift[i], dst4x4_coeffs)
	}

	sse_128_store_vector_a(dst, result[0]);
	sse_128_store_vector_a((dst+8), result[1]);
}



void sse_aligned_inv_dst_4x4(int16_t *src, int16_t *dst, int stride) 
{
	int i, shift[2] = {7,12};
	__m128i round[2];// = {sse_128_vector_i32(64),sse_128_vector_i32(2048)};
	__m128i result[2];

	__m128i l0l1 = sse128_unpacklo_u16(sse_128_loadlo_vector64((src)), sse_128_loadlo_vector64((src+4)));
	__m128i l2l3 = sse128_unpacklo_u16(sse_128_loadlo_vector64((src+8)), sse_128_loadlo_vector64((src+12)));

	round[0] = sse_128_vector_i32(64);
	round[1] = sse_128_vector_i32(2048);

	result[0]	= sse128_unpacklo_u32(l0l1, l2l3);
	result[1]	= sse128_unpackhi_u32(l0l1, l2l3);

	for (i = 0; i < 2; i++)
	{
		DST4x4(result, round[i], shift[i], inv_dst4x4_coeffs)
	}

	l0l1 = sse128_unpacklo_u16(result[0], result[1]);
	l2l3 = sse128_unpackhi_u16(result[0], result[1]);

	result[0] = sse128_unpacklo_u16(l0l1,l2l3);
	result[1] = sse128_unpackhi_u16(l0l1,l2l3);

	sse_64_storel_vector_u(dst, result[0]);
	sse_64_storeh_vector_u((dst+stride), result[0]);
	sse_64_storel_vector_u((dst+2*stride), result[1]);
	sse_64_storeh_vector_u((dst+3*stride), result[1]);
}


//---------------------------------------DCT-4x4 ----------------------------------------------


ALIGN(16) static const int16_t dct4x4_coeffs[4][8] =
{
    { 64, 64, 64, 64, 64, 64, 64, 64 },
    { 83, 36,-36,-83, 83, 36,-36,-83 },
    { 64,-64,-64, 64, 64,-64,-64, 64 },
    { 36,-83, 83,-36, 36,-83, 83,-36 }
};

ALIGN(16) static const int16_t inv_dct4x4_coeffs[4][8] =
{
    { 64, 83, 64, 36, 64, 83, 64, 36 },
    { 64, 36,-64,-83, 64, 36,-64,-83 },
    { 64,-36,-64, 83, 64,-36,-64, 83 },
    { 64,-83, 64,-36, 64,-83, 64,-36 }
};

#define DCT4x4(result, round, shift, coeffs)		DST4x4(result, round, shift, coeffs)

void sse_aligned_dct_4x4(int16_t *src, int16_t *dst, int stride) 
{
	int i, shift[2] = {1,8};
	__m128i round[2];// = {sse_128_vector_i32(1),sse_128_vector_i32(128)};
	__m128i result[2];
//	__m128i final;
	__m128i l0  = sse_128_loadlo_vector64((src));
	__m128i l1  = sse_128_loadlo_vector64((src+stride));
	__m128i l2  = sse_128_loadlo_vector64((src+2*stride));
	__m128i l3  = sse_128_loadlo_vector64((src+3*stride));

	round[0] = sse_128_vector_i32(1);
	round[1] = sse_128_vector_i32(128);

	result[0]  = sse128_unpacklo_u64(l0, l1);
	result[1]  = sse128_unpacklo_u64(l2, l3);

	for (i = 0; i < 2; i++)
	{
		DCT4x4(result, round[i], shift[i], dct4x4_coeffs)
	}

	sse_128_store_vector_a(dst, result[0]);
	sse_128_store_vector_a((dst+8), result[1]);
}



void sse_aligned_inv_dct_4x4(int16_t *src, int16_t *dst, int stride) 
{
	int i, shift[2] = {7,12};
	__m128i round[2];// = {sse_128_vector_i32(64),sse_128_vector_i32(2048)};
	__m128i result[2];

	__m128i l0l1 = sse128_unpacklo_u16(sse_128_loadlo_vector64((src)), sse_128_loadlo_vector64((src+4)));
	__m128i l2l3 = sse128_unpacklo_u16(sse_128_loadlo_vector64((src+8)), sse_128_loadlo_vector64((src+12)));

	round[0] = sse_128_vector_i32(64);
	round[1] = sse_128_vector_i32(2048);

	result[0]	= sse128_unpacklo_u32(l0l1, l2l3);
	result[1]	= sse128_unpackhi_u32(l0l1, l2l3);

	for (i = 0; i < 2; i++)
	{
		DCT4x4(result, round[i], shift[i], inv_dct4x4_coeffs)
	}

	l0l1 = sse128_unpacklo_u16(result[0], result[1]);
	l2l3 = sse128_unpackhi_u16(result[0], result[1]);

	result[0] = sse128_unpacklo_u16(l0l1,l2l3);
	result[1] = sse128_unpackhi_u16(l0l1,l2l3);

	sse_64_storel_vector_u(dst, result[0]);
	sse_64_storeh_vector_u((dst+stride), result[0]);
	sse_64_storel_vector_u((dst+2*stride), result[1]);
	sse_64_storeh_vector_u((dst+3*stride), result[1]);
}

//---------------------------------------DCT-8x8 ----------------------------------------------



ALIGN(16) const int16_t dct8x8_buterfly_coeffs[8][8] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64},	//EE
  { 89, 75, 50, 18, 89, 75, 50, 18},	
  { 83, 36,-36,-83, 83, 36,-36,-83},	//EO
  { 75,-18,-89,-50, 75,-18,-89,-50},
  { 64,-64,-64, 64, 64,-64,-64, 64},	//EE
  { 50,-89, 18, 75, 50,-89, 18, 75},
  { 36,-83, 83,-36, 36,-83, 83,-36},	//EO
  { 18,-50, 75,-89, 18,-50, 75,-89}
};


ALIGN(16) const int16_t dct8x8_coeffs[8][8] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64},
  { 89, 75, 50, 18,-18,-50,-75,-89},
  { 83, 36,-36,-83,-83,-36, 36, 83},
  { 75,-18,-89,-50, 50, 89, 18,-75},
  { 64,-64,-64, 64, 64,-64,-64, 64},
  { 50,-89, 18, 75,-75,-18, 89,-50},
  { 36,-83, 83,-36,-36, 83,-83, 36},
  { 18,-50, 75,-89, 89,-75, 50,-18}
};


#define DCT8x8_1LINE(src, coeffs, round, shift, dst)																															\
		__m128i l0l1 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[0], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[1], sse_128_load_vector_a(coeffs)));				\
		__m128i l2l3 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[2], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[3], sse_128_load_vector_a(coeffs)));				\
		__m128i l0l1l2l3 = sse_128_hadd_i32(l0l1, l2l3);																														\
		__m128i l4l5 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[4], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[5], sse_128_load_vector_a(coeffs)));				\
		__m128i l6l7 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[6], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[7], sse_128_load_vector_a(coeffs)));				\
		__m128i l4l5l6l7 = sse_128_hadd_i32(l4l5, l6l7);																														\
		dst = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(l0l1l2l3, round), shift), sse_128_shift_r_i32(sse_128_add_i32(l4l5l6l7, round), shift));


void sse_aligned_dct_8x8_1(int16_t *src, int16_t *dst, int stride) 
{
	int i;
	__m128i round;// = {sse_128_vector_i32(2),sse_128_vector_i32(256)};

	__m128i aux[8], result[8];

	result[0] = sse_128_load_vector_a(src);
	result[1] = sse_128_load_vector_a(src+stride);
	result[2] = sse_128_load_vector_a(src+2*stride);
	result[3] = sse_128_load_vector_a(src+3*stride);
	result[4] = sse_128_load_vector_a(src+4*stride);
	result[5] = sse_128_load_vector_a(src+5*stride);
	result[6] = sse_128_load_vector_a(src+6*stride);
	result[7] = sse_128_load_vector_a(src+7*stride);

	round = sse_128_vector_i32(2);
	for (i = 0; i < 8; i++)
	{
		DCT8x8_1LINE(result, dct8x8_coeffs[i], round, 2, aux[i])
	}

	round = sse_128_vector_i32(256);
	for (i = 0; i < 8; i++)
	{
		DCT8x8_1LINE(aux, dct8x8_coeffs[i], round, 9, result[i])
	}

	sse_128_store_vector_a(dst, result[0]);
	sse_128_store_vector_a(dst+8, result[1]);
	sse_128_store_vector_a(dst+16, result[2]);
	sse_128_store_vector_a(dst+24, result[3]);
	sse_128_store_vector_a(dst+32, result[4]);
	sse_128_store_vector_a(dst+40, result[5]);
	sse_128_store_vector_a(dst+48, result[6]);
	sse_128_store_vector_a(dst+56, result[7]);
}




#define DCT8x8_PASS1(E_l0l1,E_l2l3,E_l4l5,E_l6l7, coeffs, round, shift, result)																							\
		__m128i  dst0l = sse_128_hadd_i32(sse_128_madd_i16_i32(E_l0l1, sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(E_l2l3, sse_128_load_vector_a(coeffs)));		\
		__m128i  dst0h = sse_128_hadd_i32(sse_128_madd_i16_i32(E_l4l5, sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(E_l6l7, sse_128_load_vector_a(coeffs)));		\
		dst0l = sse_128_shift_r_i32(sse_128_add_i32(dst0l, round), shift);																									\
		dst0h = sse_128_shift_r_i32(sse_128_add_i32(dst0h, round), shift);																									\
		result = sse128_packs_u32_u16(dst0l, dst0h);													

//just for 8bit depth
void sse_aligned_dct_8x8(int16_t *src, int16_t *dst, int stride) 
{
	int i, shift[2] = {2,9};
	__m128i round[2];// = {sse_128_vector_i32(2),sse_128_vector_i32(256)};
	__m128i aux[8];
	__m128i result[8];
	__m128i shuff_mask;

	round[0] = sse_128_vector_i32(2);
	round[1] = sse_128_vector_i32(256);

	result[0] = sse_128_load_vector_a(src);
	result[1] = sse_128_load_vector_a(src+stride);
	result[2] = sse_128_load_vector_a(src+2*stride);
	result[3] = sse_128_load_vector_a(src+3*stride);
	result[4] = sse_128_load_vector_a(src+4*stride);
	result[5] = sse_128_load_vector_a(src+5*stride);
	result[6] = sse_128_load_vector_a(src+6*stride);
	result[7] = sse_128_load_vector_a(src+7*stride);

	shuff_mask = sse_128_load_vector_a(shuffle_mask_dct_16_3);//0,1,2,3,4,5,6,7 -> 0,7,1,6,2,5,3,4

	{
		__m128i l0 = sse_128_shuffle_8(result[0], shuff_mask);
		__m128i l1 = sse_128_shuffle_8(result[1], shuff_mask);
		__m128i l2 = sse_128_shuffle_8(result[2], shuff_mask);
		__m128i l3 = sse_128_shuffle_8(result[3], shuff_mask);
		__m128i l4 = sse_128_shuffle_8(result[4], shuff_mask);
		__m128i l5 = sse_128_shuffle_8(result[5], shuff_mask);
		__m128i l6 = sse_128_shuffle_8(result[6], shuff_mask);
		__m128i l7 = sse_128_shuffle_8(result[7], shuff_mask);

		//esto satura en la segunda pasada - Si se mantuviese la precision se podrian hacer asi las dos pasadas
		__m128i  E_l0l1= sse_128_hadd_i16(l0, l1);
		__m128i  O_l0l1= sse_128_hsub_i16(l0, l1);
		__m128i  E_l2l3= sse_128_hadd_i16(l2, l3);
		__m128i  O_l2l3= sse_128_hsub_i16(l2, l3);
		__m128i  E_l4l5= sse_128_hadd_i16(l4, l5);
		__m128i  O_l4l5= sse_128_hsub_i16(l4, l5);
		__m128i  E_l6l7= sse_128_hadd_i16(l6, l7);
		__m128i  O_l6l7= sse_128_hsub_i16(l6, l7);

		{DCT8x8_PASS1(E_l0l1,E_l2l3,E_l4l5,E_l6l7, dct8x8_buterfly_coeffs[0], round[0], shift[0], aux[0])}
		{DCT8x8_PASS1(E_l0l1,E_l2l3,E_l4l5,E_l6l7, dct8x8_buterfly_coeffs[2], round[0], shift[0], aux[2])}
		{DCT8x8_PASS1(E_l0l1,E_l2l3,E_l4l5,E_l6l7, dct8x8_buterfly_coeffs[4], round[0], shift[0], aux[4])}
		{DCT8x8_PASS1(E_l0l1,E_l2l3,E_l4l5,E_l6l7, dct8x8_buterfly_coeffs[6], round[0], shift[0], aux[6])}


	/*	__m128i  dst0l = sse_128_hadd_i32(sse_128_madd_i16_i32(E_l0l1, sse_128_load_vector_a(dct8x8_buterfly_coeffs[6])), sse_128_madd_i16_i32(E_l2l3, sse_128_load_vector_a(dct8x8_buterfly_coeffs[6])));
		__m128i  dst0h = sse_128_hadd_i32(sse_128_madd_i16_i32(E_l4l5, sse_128_load_vector_a(dct8x8_buterfly_coeffs[6])), sse_128_madd_i16_i32(E_l6l7, sse_128_load_vector_a(dct8x8_buterfly_coeffs[6])));
		dst0l = sse_128_shift_r_i32(sse_128_add_i32(dst0l, round[0]), shift[0]);																		
		dst0h = sse_128_shift_r_i32(sse_128_add_i32(dst0h, round[0]), shift[0]);																		
		aux[6] = sse128_packs_u32_u16(dst0l, dst0h);													
	*/



		{DCT8x8_PASS1(O_l0l1,O_l2l3,O_l4l5,O_l6l7,dct8x8_buterfly_coeffs[1], round[0], shift[0], aux[1])}
		{DCT8x8_PASS1(O_l0l1,O_l2l3,O_l4l5,O_l6l7,dct8x8_buterfly_coeffs[3], round[0], shift[0], aux[3])}
		{DCT8x8_PASS1(O_l0l1,O_l2l3,O_l4l5,O_l6l7,dct8x8_buterfly_coeffs[5], round[0], shift[0], aux[5])}
		{DCT8x8_PASS1(O_l0l1,O_l2l3,O_l4l5,O_l6l7,dct8x8_buterfly_coeffs[7], round[0], shift[0], aux[7])}
	}

	for (i = 0; i < 8; i++)
	{
		DCT8x8_1LINE(aux, dct8x8_coeffs[i], round[1], 9, result[i])
	}

	sse_128_store_vector_a(dst, result[0]);
	sse_128_store_vector_a(dst+8, result[1]);
	sse_128_store_vector_a(dst+16, result[2]);
	sse_128_store_vector_a(dst+24, result[3]);
	sse_128_store_vector_a(dst+32, result[4]);
	sse_128_store_vector_a(dst+40, result[5]);
	sse_128_store_vector_a(dst+48, result[6]);
	sse_128_store_vector_a(dst+56, result[7]);
}






ALIGN(16) const int16_t inv_dct8x8_odd[2][8] =
{
  { 89, 75, 50, 18, 75,-18,-89,-50},
  { 50,-89, 18, 75, 18,-50, 75,-89},
};

ALIGN(16) const int16_t inv_dct8x8_even[2][8] =
{
  { 64, 83, 64, 36, 64, 36,-64,-83},
  { 64,-36,-64, 83, 64,-83, 64,-36},
};


void sse_aligned_inv_dct_8x8(int16_t *src, int16_t *dst, int stride) 
{
	int i, shift;
	int iround[2] = {64, 2048};
	int ishift[2] = {7, 12};
	__m128i round;

	__m128i result[8];

	result[0] = sse_128_load_vector_a(src);
	result[1] = sse_128_load_vector_a(src+8);
	result[2] = sse_128_load_vector_a(src+2*8);
	result[3] = sse_128_load_vector_a(src+3*8);
	result[4] = sse_128_load_vector_a(src+4*8);
	result[5] = sse_128_load_vector_a(src+5*8);
	result[6] = sse_128_load_vector_a(src+6*8);
	result[7] = sse_128_load_vector_a(src+7*8);

	for (i = 0; i < 2; i++)
	{
		__m128i l1l3low = sse128_unpacklo_u16(result[1],result[3]);
		__m128i l1l3high = sse128_unpackhi_u16(result[1],result[3]);
		__m128i l5l7low = sse128_unpacklo_u16(result[5],result[7]);
		__m128i l5l7high = sse128_unpackhi_u16(result[5],result[7]);
		__m128i odd_c0c1 = sse128_unpacklo_u32(l1l3low,l5l7low);		//odd members of colums 0 and 1 in a row (10,30,50,70,11,31,51,71)
		__m128i odd_c2c3 = sse128_unpackhi_u32(l1l3low,l5l7low);		//odd members of colums 2 and 3 in a row (12,32,52,72,13,33,53,73)
		__m128i odd_c4c5 = sse128_unpacklo_u32(l1l3high,l5l7high);		//odd members of colums 4 and 5 in a row (14,34,54,74,15,35,55,75)
		__m128i odd_c6c7 = sse128_unpackhi_u32(l1l3high,l5l7high);		//odd members of colums 6 and 7 in a row (16,36,56,76,17,37,57,77)

		__m128i oddc0_x2  = sse128_unpacklo_u64(odd_c0c1,odd_c0c1);
		__m128i oddc1_x2  = sse128_unpackhi_u64(odd_c0c1,odd_c0c1);
		__m128i oddc2_x2  = sse128_unpacklo_u64(odd_c2c3,odd_c2c3);
		__m128i oddc3_x2  = sse128_unpackhi_u64(odd_c2c3,odd_c2c3);
		__m128i oddc4_x2  = sse128_unpacklo_u64(odd_c4c5,odd_c4c5);
		__m128i oddc5_x2  = sse128_unpackhi_u64(odd_c4c5,odd_c4c5);
		__m128i oddc6_x2  = sse128_unpacklo_u64(odd_c6c7,odd_c6c7);
		__m128i oddc7_x2  = sse128_unpackhi_u64(odd_c6c7,odd_c6c7);
		__m128i O_c0 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc0_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc0_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c1 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc1_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc1_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c2 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc2_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc2_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c3 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc3_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc3_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c4 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc4_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc4_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c5 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc5_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc5_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c6 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc6_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc6_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c7 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc7_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc7_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));


		__m128i l0l2low = sse128_unpacklo_u16(result[0],result[2]);
		__m128i l0l2high = sse128_unpackhi_u16(result[0],result[2]);
		__m128i l4l6low = sse128_unpacklo_u16(result[4],result[6]);
		__m128i l4l6high = sse128_unpackhi_u16(result[4],result[6]);
		__m128i even_c0c1 = sse128_unpacklo_u32(l0l2low,l4l6low);		//even members of colums 0 and 1 in a row (00,20,40,60,01,21,41,61)
		__m128i even_c2c3 = sse128_unpackhi_u32(l0l2low,l4l6low);		//even members of colums 2 and 3 in a row (02,22,42,62,03,23,43,63)
		__m128i even_c4c5 = sse128_unpacklo_u32(l0l2high,l4l6high);		//even members of colums 4 and 5 in a row (04,24,44,64,05,25,45,65)
		__m128i even_c6c7 = sse128_unpackhi_u32(l0l2high,l4l6high);		//even members of colums 6 and 7 in a row (06,26,46,66,07,27,47,67)

		__m128i evenc0_x2  = sse128_unpacklo_u64(even_c0c1,even_c0c1);
		__m128i evenc1_x2  = sse128_unpackhi_u64(even_c0c1,even_c0c1);
		__m128i evenc2_x2  = sse128_unpacklo_u64(even_c2c3,even_c2c3);
		__m128i evenc3_x2  = sse128_unpackhi_u64(even_c2c3,even_c2c3);
		__m128i evenc4_x2  = sse128_unpacklo_u64(even_c4c5,even_c4c5);
		__m128i evenc5_x2  = sse128_unpackhi_u64(even_c4c5,even_c4c5);
		__m128i evenc6_x2  = sse128_unpacklo_u64(even_c6c7,even_c6c7);
		__m128i evenc7_x2  = sse128_unpackhi_u64(even_c6c7,even_c6c7);
		__m128i E_c0 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc0_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc0_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c1 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc1_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc1_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c2 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc2_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc2_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c3 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc3_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc3_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c4 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc4_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc4_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c5 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc5_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc5_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c6 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc6_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc6_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c7 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc7_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc7_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));

		round = sse_128_vector_i32(iround[i]);
		shift = ishift[i];

		result[0] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c0, O_c0), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c0, O_c0), round), shift), 0x1b));	
		result[1] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c1, O_c1), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c1, O_c1), round), shift), 0x1b));	
		result[2] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c2, O_c2), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c2, O_c2), round), shift), 0x1b));	
		result[3] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c3, O_c3), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c3, O_c3), round), shift), 0x1b));	
		result[4] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c4, O_c4), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c4, O_c4), round), shift), 0x1b));	
		result[5] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c5, O_c5), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c5, O_c5), round), shift), 0x1b));	
		result[6] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c6, O_c6), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c6, O_c6), round), shift), 0x1b));	
		result[7] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c7, O_c7), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c7, O_c7), round), shift), 0x1b));	
	}

	sse_128_store_vector_a(dst, result[0]);
	sse_128_store_vector_a(dst+stride, result[1]);
	sse_128_store_vector_a(dst+2*stride, result[2]);
	sse_128_store_vector_a(dst+3*stride, result[3]);
	sse_128_store_vector_a(dst+4*stride, result[4]);
	sse_128_store_vector_a(dst+5*stride, result[5]);
	sse_128_store_vector_a(dst+6*stride, result[6]);
	sse_128_store_vector_a(dst+7*stride, result[7]);

}

/*
void sse_aligned_inv_dct_8x8_v2(int16_t *src, int16_t *dst, int dst_stride, int16_t *aux) 
{
	int i, shift;
	int iround[2] = {64, 2048};
	int ishift[2] = {7, 12};
	int astride[2] = {8,dst_stride};
	int stride;
	__m128i round;

//	__m128i aux[8], result[8];
	__m128i _128aux;
	int16_t *psrc, *pdst;
	int16_t *asrc[2] = {src,aux};
	int16_t *adst[2] = {aux,dst};

	for (i = 0; i < 2; i++)
	{
		psrc = asrc[i];
		pdst = adst[i];
		stride = astride[i];
		__m128i l1l3low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+1*8),sse_128_load_vector_a(psrc+3*8));
		__m128i l1l3high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+1*8),sse_128_load_vector_a(psrc+3*8));
		__m128i l5l7low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+5*8),sse_128_load_vector_a(psrc+7*8));
		__m128i l5l7high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+5*8),sse_128_load_vector_a(psrc+7*8));
		__m128i odd_c0c1 = sse128_unpacklo_u32(l1l3low,l5l7low);		//odd members of colums 0 and 1 in a row (10,30,50,70,11,31,51,71)
		__m128i odd_c2c3 = sse128_unpackhi_u32(l1l3low,l5l7low);		//odd members of colums 2 and 3 in a row (12,32,52,72,13,33,53,73)
		__m128i odd_c4c5 = sse128_unpacklo_u32(l1l3high,l5l7high);		//odd members of colums 4 and 5 in a row (14,34,54,74,15,35,55,75)
		__m128i odd_c6c7 = sse128_unpackhi_u32(l1l3high,l5l7high);		//odd members of colums 6 and 7 in a row (16,36,56,76,17,37,57,77)

		__m128i oddc0_x2  = sse128_unpacklo_u64(odd_c0c1,odd_c0c1);
		__m128i oddc1_x2  = sse128_unpackhi_u64(odd_c0c1,odd_c0c1);
		__m128i oddc2_x2  = sse128_unpacklo_u64(odd_c2c3,odd_c2c3);
		__m128i oddc3_x2  = sse128_unpackhi_u64(odd_c2c3,odd_c2c3);
		__m128i oddc4_x2  = sse128_unpacklo_u64(odd_c4c5,odd_c4c5);
		__m128i oddc5_x2  = sse128_unpackhi_u64(odd_c4c5,odd_c4c5);
		__m128i oddc6_x2  = sse128_unpacklo_u64(odd_c6c7,odd_c6c7);
		__m128i oddc7_x2  = sse128_unpackhi_u64(odd_c6c7,odd_c6c7);
		__m128i O_c0 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc0_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc0_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c1 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc1_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc1_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c2 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc2_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc2_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c3 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc3_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc3_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c4 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc4_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc4_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c5 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc5_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc5_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c6 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc6_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc6_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));
		__m128i O_c7 = sse_128_hadd_i32(sse_128_madd_i16_i32(oddc7_x2,sse_128_load_vector_a(inv_dct8x8_odd[0])),  sse_128_madd_i16_i32(oddc7_x2,sse_128_load_vector_a(inv_dct8x8_odd[1])));


		__m128i l0l2low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+0*8),sse_128_load_vector_a(psrc+2*8));
		__m128i l0l2high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+0*8),sse_128_load_vector_a(psrc+2*8));
		__m128i l4l6low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+4*8),sse_128_load_vector_a(psrc+6*8));
		__m128i l4l6high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+4*8),sse_128_load_vector_a(psrc+6*8));
		__m128i even_c0c1 = sse128_unpacklo_u32(l0l2low,l4l6low);		//even members of colums 0 and 1 in a row (00,20,40,60,01,21,41,61)
		__m128i even_c2c3 = sse128_unpackhi_u32(l0l2low,l4l6low);		//even members of colums 2 and 3 in a row (02,22,42,62,03,23,43,63)
		__m128i even_c4c5 = sse128_unpacklo_u32(l0l2high,l4l6high);		//even members of colums 4 and 5 in a row (04,24,44,64,05,25,45,65)
		__m128i even_c6c7 = sse128_unpackhi_u32(l0l2high,l4l6high);		//even members of colums 6 and 7 in a row (06,26,46,66,07,27,47,67)

		__m128i evenc0_x2  = sse128_unpacklo_u64(even_c0c1,even_c0c1);
		__m128i evenc1_x2  = sse128_unpackhi_u64(even_c0c1,even_c0c1);
		__m128i evenc2_x2  = sse128_unpacklo_u64(even_c2c3,even_c2c3);
		__m128i evenc3_x2  = sse128_unpackhi_u64(even_c2c3,even_c2c3);
		__m128i evenc4_x2  = sse128_unpacklo_u64(even_c4c5,even_c4c5);
		__m128i evenc5_x2  = sse128_unpackhi_u64(even_c4c5,even_c4c5);
		__m128i evenc6_x2  = sse128_unpacklo_u64(even_c6c7,even_c6c7);
		__m128i evenc7_x2  = sse128_unpackhi_u64(even_c6c7,even_c6c7);
		__m128i E_c0 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc0_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc0_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c1 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc1_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc1_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c2 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc2_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc2_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c3 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc3_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc3_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c4 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc4_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc4_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c5 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc5_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc5_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c6 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc6_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc6_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));
		__m128i E_c7 = sse_128_hadd_i32(sse_128_madd_i16_i32(evenc7_x2,sse_128_load_vector_a(inv_dct8x8_even[0])),  sse_128_madd_i16_i32(evenc7_x2,sse_128_load_vector_a(inv_dct8x8_even[1])));

		round = sse_128_vector_i32(iround[i]);
		shift = ishift[i];

		_128aux = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c0, O_c0), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c0, O_c0), round), shift), 0x1b));	
		sse_128_store_vector_a(pdst, _128aux);
		_128aux = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c1, O_c1), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c1, O_c1), round), shift), 0x1b));	
		sse_128_store_vector_a(pdst+stride, _128aux);
		_128aux = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c2, O_c2), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c2, O_c2), round), shift), 0x1b));	
		sse_128_store_vector_a(pdst+2*stride, _128aux);
		_128aux = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c3, O_c3), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c3, O_c3), round), shift), 0x1b));	
		sse_128_store_vector_a(pdst+3*stride, _128aux);
		_128aux = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c4, O_c4), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c4, O_c4), round), shift), 0x1b));	
		sse_128_store_vector_a(pdst+4*stride, _128aux);
		_128aux = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c5, O_c5), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c5, O_c5), round), shift), 0x1b));	
		sse_128_store_vector_a(pdst+5*stride, _128aux);
		_128aux = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c6, O_c6), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c6, O_c6), round), shift), 0x1b));	
		sse_128_store_vector_a(pdst+6*stride, _128aux);
		_128aux = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E_c7, O_c7), round), shift), sse_128_shuffle_32(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E_c7, O_c7), round), shift), 0x1b));	
		sse_128_store_vector_a(pdst+7*stride, _128aux);
	}
}
*/


//---------------------------------------DCT-16x16 ----------------------------------------------



ALIGN(16) const int16_t dct16x16_coeffs[16][16] =	//same as g_aiT16
{

  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90},
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87},
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80},
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70},
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57},
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43},
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25},
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9}
};

#define DCT16x16_BUTTERFLY_EVEN(src, coeffs, round, shift, result)																							\
		__m128i  dst0l = sse_128_hadd_i32(sse_128_madd_i16_i32(src[0], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[1], sse_128_load_vector_a(coeffs)));		\
		__m128i  dst0h = sse_128_hadd_i32(sse_128_madd_i16_i32(src[2], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[3], sse_128_load_vector_a(coeffs)));		\
		dst0l = sse_128_shift_r_i32(sse_128_add_i32(dst0l, round), shift);																									\
		dst0h = sse_128_shift_r_i32(sse_128_add_i32(dst0h, round), shift);																									\
		result = sse128_packs_u32_u16(dst0l, dst0h);													

#define DCT16x16_BUTTERFLY_ODD(src, coeffs, round, shift, dst)																															\
		__m128i l0l1 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[0], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[1], sse_128_load_vector_a(coeffs)));				\
		__m128i l2l3 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[2], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[3], sse_128_load_vector_a(coeffs)));				\
		__m128i l0l1l2l3 = sse_128_hadd_i32(l0l1, l2l3);																														\
		__m128i l4l5 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[4], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[5], sse_128_load_vector_a(coeffs)));				\
		__m128i l6l7 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[6], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[7], sse_128_load_vector_a(coeffs)));				\
		__m128i l4l5l6l7 = sse_128_hadd_i32(l4l5, l6l7);																														\
		dst = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(l0l1l2l3, round), shift), sse_128_shift_r_i32(sse_128_add_i32(l4l5l6l7, round), shift));		

#define DCT16x16_1LINE(src, coeffs, round, shift, dst)																															\
		__m128i l0l1 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[0], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[1], sse_128_load_vector_a(coeffs)));				\
		__m128i l2l3 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[2], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[3], sse_128_load_vector_a(coeffs)));				\
		__m128i l0l1l2l3 = sse_128_hadd_i32(l0l1, l2l3);																														\
		__m128i l4l5 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[4], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[5], sse_128_load_vector_a(coeffs)));				\
		__m128i l6l7 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[6], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[7], sse_128_load_vector_a(coeffs)));				\
		__m128i l4l5l6l7 = sse_128_hadd_i32(l4l5, l6l7);																														\
		dst[0] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(l0l1l2l3, round), shift), sse_128_shift_r_i32(sse_128_add_i32(l4l5l6l7, round), shift));				\
		l0l1 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[8], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[9], sse_128_load_vector_a(coeffs)));						\
		l2l3 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[10], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[11], sse_128_load_vector_a(coeffs)));					\
		l0l1l2l3 = sse_128_hadd_i32(l0l1, l2l3);																																\
		l4l5 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[12], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[13], sse_128_load_vector_a(coeffs)));					\
		l6l7 = sse_128_hadd_i32(sse_128_madd_i16_i32(src[14], sse_128_load_vector_a(coeffs)), sse_128_madd_i16_i32(src[15], sse_128_load_vector_a(coeffs)));					\
		l4l5l6l7 = sse_128_hadd_i32(l4l5, l6l7);																																\
		dst[1] = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(l0l1l2l3, round), shift), sse_128_shift_r_i32(sse_128_add_i32(l4l5l6l7, round), shift));		


ALIGN(16) const int16_t dct16x16_butterfly_coeffs[16][8] =
{
	{64, 64, 64, 64, 64, 64, 64, 64},	//ln 0 -> coeff00, coeff01, coeff01, coeff00, coeff00, coeff01, coeff01, coeff00					- EE
	{90, 87, 80, 70, 57, 43, 25,  9},	//ln 1 -> coeff10, coeff11, coeff12, coeff13, coeff14, coeff15, coeff16, coeff17					- O
	{89, 75, 50, 18, 89, 75, 50, 18},	//ln 2 -> coeff20, coeff21, coeff22, coeff23, coeff20, coeff21, coeff22, coeff23					- EO
	{87, 57,  9,-43,-80,-90,-70,-25},	//ln 3 -> coeff30, coeff31, coeff32, coeff33, coeff34, coeff35, coeff36, coeff37					- O
	{83, 36,-36,-83, 83, 36,-36,-83},	//ln 4 -> coeff40, coeff41,-coeff41,-coeff40, coeff40, coeff41,-coeff41,-coeff40					- EE
	{80,  9,-70,-87,-25, 57, 90, 43},	//ln 5 -> coeff50, coeff51, coeff52, coeff53, coeff54, coeff55, coeff56, coeff57					- O
	{75,-18,-89,-50, 75,-18,-89,-50},	//ln 6 -> coeff60, coeff61, coeff62, coeff63, coeff60, coeff61, coeff62, coeff63					- EO
	{70,-43,-87,  9, 90, 25,-80,-57},	//ln 7 -> coeff70, coeff71, coeff72, coeff73, coeff74, coeff75, coeff76, coeff77					- O
	{64,-64,-64, 64, 64,-64,-64, 64},	//ln 8 -> coeff80, coeff81, coeff81, coeff80, coeff80, coeff81, coeff81, coeff80					- EE
	{57,-80,-25, 90, -9,-87, 43, 70},	//ln 9 -> coeff90, coeff91, coeff92, coeff93, coeff94, coeff95, coeff96, coeff97					- O
	{50,-89, 18, 75, 50,-89, 18, 75},	//ln 10 -> coeff100, coeff101, coeff102, coeff103, coeff100, coeff101, coeff102, coeff103			- EO
	{43,-90, 57, 25,-87, 70,  9,-80},	//ln 11 -> coeff110, coeff111, coeff112, coeff113, coeff114, coeff115, coeff116, coeff117			- O
	{36,-83, 83,-36, 36,-83, 83,-36},	//ln 12 -> coeff120, coeff121,-coeff121,-coeff120, coeff120, coeff121,-coeff121,-coeff120			- EE
	{25,-70, 90,-80, 43,  9,-57, 87},	//ln 13 -> coeff130, coeff131, coeff132, coeff133, coeff134, coeff135, coeff136, coeff137			- O
	{18,-50, 75,-89, 18,-50, 75,-89},	//ln 14 -> coeff140, coeff141, coeff142, coeff143, coeff140, coeff141, coeff142, coeff143			- EO
	{ 9,-25, 43,-57, 70,-80, 87,-90}	//ln 15 -> coeff150, coeff151, coeff152, coeff153, coeff154, coeff155, coeff156, coeff157			- O
};



void sse_aligned_dct_16x16(int16_t *src, int16_t *dst, int stride) 
{
	int i;
	__m128i round[2];

	__m128i O[16];
	__m128i EE[8], EO[8];
	__m128i /*result[32], */aux[32];
	__m128i  aux1, aux0;

	__m128i shuff_mask = sse_128_load_vector_a(shuffle_mask_dct_16_0);//0,1,2,3,4,5,6,7 -> 7,6,5,4,3,2,1,0
	__m128i E_l[16];

	round[0] = sse_128_vector_i32(4);
	round[1] = sse_128_vector_i32(512);

	for (i = 0; i < 16; i++)
	{
		aux0 = sse_128_load_vector_a(src+i*stride);
		aux1 = sse_128_shuffle_8(sse_128_load_vector_a(src+i*stride+8), shuff_mask);
		E_l[i] = sse_128_add_i16(aux0, aux1);
		O[i] = sse_128_sub_i16(aux0, aux1);
	}

	shuff_mask = sse_128_load_vector_a(shuffle_mask_dct_16_1);//0,1,2,3,4,5,6,7 -> 3,2,1,0,7,6,5,4

	for (i = 0; i < 8; i++)
	{
		aux0 = sse128_unpacklo_u64(E_l[2*i],E_l[2*i+1]);
		aux1 = sse_128_shuffle_8(sse128_unpackhi_u64(E_l[2*i],E_l[2*i+1]), shuff_mask);
		EE[i] =  sse_128_add_i16(aux0, aux1);
		EO[i] =  sse_128_sub_i16(aux0, aux1);
	}

//	int write_result = 0;
	//first pass (butterfly)
	for (i = 0; i < 4; i++)
	{
		{DCT16x16_BUTTERFLY_EVEN(EE, dct16x16_butterfly_coeffs[4*i+0], round[0], 3, aux[8*i+0])}
		{DCT16x16_BUTTERFLY_EVEN((&EE[4]), dct16x16_butterfly_coeffs[4*i+0], round[0], 3, aux[8*i+1])}
		{DCT16x16_BUTTERFLY_ODD(O, dct16x16_butterfly_coeffs[4*i+1], round[0], 3, aux[8*i+2])}
		{DCT16x16_BUTTERFLY_ODD((&O[8]), dct16x16_butterfly_coeffs[4*i+1], round[0], 3, aux[8*i+3])}
		{DCT16x16_BUTTERFLY_EVEN(EO, dct16x16_butterfly_coeffs[4*i+2], round[0], 3, aux[8*i+4])}
		{DCT16x16_BUTTERFLY_EVEN((&EO[4]), dct16x16_butterfly_coeffs[4*i+2], round[0], 3, aux[8*i+5])}
		{DCT16x16_BUTTERFLY_ODD(O, dct16x16_butterfly_coeffs[4*i+3], round[0], 3, aux[8*i+6])}
		{DCT16x16_BUTTERFLY_ODD((&O[8]), dct16x16_butterfly_coeffs[4*i+3], round[0], 3, aux[8*i+7])}
/*		if(write_result)
		{
			sse_128_store_vector_a(dst, aux[8*i]);
			sse_128_store_vector_a(dst+8, aux[8*i+1]);	
			sse_128_store_vector_a(dst+16, aux[8*i+2]);	
			sse_128_store_vector_a(dst+24, aux[8*i+3]);	
			sse_128_store_vector_a(dst+32, aux[8*i+4]);	
			sse_128_store_vector_a(dst+40, aux[8*i+5]);	
			sse_128_store_vector_a(dst+48, aux[8*i+6]);	
			sse_128_store_vector_a(dst+56, aux[8*i+7]);	
			dst+=64;
		}
*/	}

	//second pass (no butterfly)
	for (i = 0; i < 16; i++)
	{
//		DCT16x16_1LINE(aux, dct16x16_coeffs[i], round[1], 10, (&(result[2*i])))
		__m128i l0 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[0], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[1], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		__m128i l1 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[2], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[3], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		__m128i l0l1 = sse_128_hadd_i32(l0, l1);					
		__m128i l2 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[4], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[5], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		__m128i l3 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[6], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[7], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		__m128i l2l3 = sse_128_hadd_i32(l2, l3);
		__m128i l0l1l2l3 = sse_128_hadd_i32(l0l1, l2l3);
		__m128i l4 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[8], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[9], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		__m128i l5 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[10], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[11], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		__m128i l4l5 = sse_128_hadd_i32(l4, l5);					
		__m128i l6 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[12], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[13], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		__m128i l7 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[14], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[15], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		__m128i l6l7 = sse_128_hadd_i32(l6, l7);
		__m128i l4l5l6l7 = sse_128_hadd_i32(l4l5, l6l7);
		aux0 = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(l0l1l2l3, round[1]), 10), sse_128_shift_r_i32(sse_128_add_i32(l4l5l6l7, round[1]), 10));	
		sse_128_store_vector_a(dst, aux0);

		l0 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[16], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[17], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		l1 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[18], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[19], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		l0l1 = sse_128_hadd_i32(l0, l1);					
		l2 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[20], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[21], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		l3 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[22], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[23], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		l2l3 = sse_128_hadd_i32(l2, l3);
		l0l1l2l3 = sse_128_hadd_i32(l0l1, l2l3);
		l4 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[24], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[25], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		l5 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[26], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[27], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		l4l5 = sse_128_hadd_i32(l4, l5);					
		l6 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[28], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[29], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		l7 = sse_128_hadd_i32(sse_128_madd_i16_i32(aux[30], sse_128_load_vector_a(dct16x16_coeffs[i])), sse_128_madd_i16_i32(aux[31], sse_128_load_vector_a(dct16x16_coeffs[i]+8)));	
		l6l7 = sse_128_hadd_i32(l6, l7);
		l4l5l6l7 = sse_128_hadd_i32(l4l5, l6l7);
		aux1 = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(l0l1l2l3, round[1]), 10), sse_128_shift_r_i32(sse_128_add_i32(l4l5l6l7, round[1]), 10));	
		sse_128_store_vector_a(dst+8, aux1);
		dst+=16;
	}
}




ALIGN(16) const int16_t inv_dct16x16_odd[8][8] =
{
  { 90, 87, 80, 70, 57, 43, 25,  9},
  { 87, 57,  9,-43,-80,-90,-70,-25},
  { 80,  9,-70,-87,-25, 57, 90, 43},
  { 70,-43,-87,  9, 90, 25,-80,-57},
  { 57,-80,-25, 90, -9,-87, 43, 70},
  { 43,-90, 57, 25,-87, 70,  9,-80},
  { 25,-70, 90,-80, 43,  9,-57, 87},
  {  9,-25, 43,-57, 70,-80, 87,-90}
};

ALIGN(16) const int16_t inv_dct16x16_eo[][8] =
{
  { 89, 75, 50, 18, 75,-18,-89,-50},
  { 50,-89, 18, 75, 18,-50, 75,-89},
};

ALIGN(16) const int16_t inv_dct16x16_ee[][8] =
{
  { 64, 83, 64, 36, 64, 36,-64,-83},
  { 64,-36,-64, 83, 64,-83, 64,-36},
};


void sse_aligned_inv_dct_16x16(int16_t *src, int16_t *dst, int dst_stride, int16_t *aux) 
{
	int shift;
	int iround[2] = {64, 2048};
	int ishift[2] = {7, 12};
	__m128i round;
	__m128i O[4], EO[2], EE[2], E[4]; 

	__m128i odd_c[16];
	__m128i ee_c[8], eo_c[8];
	__m128i shuff_mask_0 = sse_128_load_vector_a(shuffle_mask_dct_16_0);//0,1,2,3,4,5,6,7 -> 7,6,5,4,3,2,1,0
	int16_t *pdst, *psrc;
	int16_t *adst[2] = {aux, dst};
	int16_t *asrc[2] = {src, aux};
	int astride[2] = {16,dst_stride}; 
	int stride;
	int i, k;

	for (k = 0; k < 2; k++)
	{
/*		psrc = asrc[k];
		pdst = adst[k];
		stride = astride[k];
*/
		__m128i l1lowl3low_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+2*8),sse_128_load_vector_a(asrc[k]+6*8));			//ln10,ln30,ln11,ln31,ln12,ln32,ln13,ln33 
		__m128i l1lowl3low_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+2*8),sse_128_load_vector_a(asrc[k]+6*8));			//ln14,ln34,ln15,ln35,ln16,ln36,ln17,ln37 
		__m128i l1highl3high_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+3*8),sse_128_load_vector_a(asrc[k]+7*8));		//ln18,ln38,ln19,ln39,ln110,ln310,ln111,ln311
		__m128i l1highl3high_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+3*8),sse_128_load_vector_a(asrc[k]+7*8));		//ln112,ln312,ln113,ln313,ln114,ln314,ln115,ln315
		__m128i l5lowl7low_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+10*8),sse_128_load_vector_a(asrc[k]+14*8));
		__m128i l5lowl7low_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+10*8),sse_128_load_vector_a(asrc[k]+14*8));
		__m128i l5highl7high_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+11*8),sse_128_load_vector_a(asrc[k]+15*8));
		__m128i l5highl7high_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+11*8),sse_128_load_vector_a(asrc[k]+15*8));

		__m128i l1357_low_low_low = sse128_unpacklo_u32(l1lowl3low_low,l5lowl7low_low);			//ln10,ln30,ln50,ln70,ln11,ln31,ln51,ln71
		__m128i l1357_low_low_high = sse128_unpackhi_u32(l1lowl3low_low,l5lowl7low_low);		//ln12,ln32,ln52,ln72,ln13,ln33,ln53,ln73
		__m128i l1357_low_high_low = sse128_unpacklo_u32(l1lowl3low_high,l5lowl7low_high);		//ln14,ln34,ln54,ln74,ln15,ln35,ln55,ln75
		__m128i l1357_low_high_high = sse128_unpackhi_u32(l1lowl3low_high,l5lowl7low_high);		//ln16,ln36,ln56,ln76,ln17,ln37,ln57,ln77
		__m128i l1357_high_low_low = sse128_unpacklo_u32(l1highl3high_low,l5highl7high_low);	//ln18,ln38,ln58,ln78,ln19,ln39,ln59,ln79
		__m128i l1357_high_low_high = sse128_unpackhi_u32(l1highl3high_low,l5highl7high_low);	//ln110,ln310,ln510,ln710,ln111,ln311,ln511,ln711
		__m128i l1357_high_high_low = sse128_unpacklo_u32(l1highl3high_high,l5highl7high_high);	//ln112,ln312,ln512,ln712,ln113,ln313,ln513,ln713
		__m128i l1357_high_high_high = sse128_unpackhi_u32(l1highl3high_high,l5highl7high_high);//ln114,ln314,ln514,ln714,ln1115,ln3115,ln5115,ln7115


		__m128i l9lowl11low_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+18*8),sse_128_load_vector_a(asrc[k]+22*8));
		__m128i l9lowl11low_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+18*8),sse_128_load_vector_a(asrc[k]+22*8));
		__m128i l9highl11high_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+19*8),sse_128_load_vector_a(asrc[k]+23*8));
		__m128i l9highl11high_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+19*8),sse_128_load_vector_a(asrc[k]+23*8));
		__m128i l13lowl15low_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+26*8),sse_128_load_vector_a(asrc[k]+30*8));
		__m128i l13lowl15low_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+26*8),sse_128_load_vector_a(asrc[k]+30*8));
		__m128i l13highl15high_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+27*8),sse_128_load_vector_a(asrc[k]+31*8));
		__m128i l13highl15high_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+27*8),sse_128_load_vector_a(asrc[k]+31*8));

		__m128i l9111315_low_low_low = sse128_unpacklo_u32(l9lowl11low_low,l13lowl15low_low);			//ln90,ln110,ln130,ln150,ln91,ln111,ln131,ln151
		__m128i l9111315_low_low_high = sse128_unpackhi_u32(l9lowl11low_low,l13lowl15low_low);			//ln92,ln112,ln132,ln152,ln93,ln113,ln133,ln153
		__m128i l9111315_low_high_low = sse128_unpacklo_u32(l9lowl11low_high,l13lowl15low_high);		//ln94,ln114,ln134,ln154,ln95,ln115,ln135,ln155
		__m128i l9111315_low_high_high = sse128_unpackhi_u32(l9lowl11low_high,l13lowl15low_high);		//ln96,ln116,ln136,ln156,ln97,ln117,ln137,ln157
		__m128i l9111315_high_low_low = sse128_unpacklo_u32(l9highl11high_low,l13highl15high_low);		//ln98,ln118,ln138,ln158,ln99,ln119,ln139,ln159
		__m128i l9111315_high_low_high = sse128_unpackhi_u32(l9highl11high_low,l13highl15high_low);		//ln910,ln1110,ln1310,ln1510,ln911,ln1111,ln1311,ln1511
		__m128i l9111315_high_high_low = sse128_unpacklo_u32(l9highl11high_high,l13highl15high_high);	//ln912,ln1112,ln1312,ln1512,ln913,ln1113,ln1313,ln1513
		__m128i l9111315_high_high_high = sse128_unpackhi_u32(l9highl11high_high,l13highl15high_high);	//ln914,ln1114,ln1314,ln1514,ln9115,ln11115,ln1315,ln1515


		//lines 0,4,8,12 for EE
		__m128i l0lowl4low_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+0*8),sse_128_load_vector_a(asrc[k]+8*8));			//ln00,ln40,ln01,ln41,ln02,ln42,ln03,ln43 
		__m128i l0lowl4low_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+0*8),sse_128_load_vector_a(asrc[k]+8*8));			//ln04,ln44,ln05,ln45,ln06,ln46,ln07,ln47 
		__m128i l0highl4high_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+1*8),sse_128_load_vector_a(asrc[k]+9*8));		//ln08,ln48,ln09,ln49,ln010,ln410,ln011,ln411
		__m128i l0highl4high_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+1*8),sse_128_load_vector_a(asrc[k]+9*8));		//ln012,ln412,ln013,ln413,ln014,ln414,ln015,ln415
		__m128i l8lowl12low_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+16*8),sse_128_load_vector_a(asrc[k]+24*8));
		__m128i l8lowl12low_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+16*8),sse_128_load_vector_a(asrc[k]+24*8));
		__m128i l8highl12high_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+17*8),sse_128_load_vector_a(asrc[k]+25*8));
		__m128i l8highl12high_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+17*8),sse_128_load_vector_a(asrc[k]+25*8));


		//lines 2,6,10,14 for EO
		__m128i l2lowl6low_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+4*8),sse_128_load_vector_a(asrc[k]+12*8));			//ln20,ln60,ln21,ln61,ln22,ln62,ln23,ln63 
		__m128i l2lowl6low_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+4*8),sse_128_load_vector_a(asrc[k]+12*8));		//ln24,ln64,ln25,ln65,ln26,ln66,ln27,ln67 
		__m128i l2highl6high_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+5*8),sse_128_load_vector_a(asrc[k]+13*8));		//ln28,ln68,ln29,ln69,ln210,ln610,ln211,ln611
		__m128i l2highl6high_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+5*8),sse_128_load_vector_a(asrc[k]+13*8));		//ln212,ln612,ln213,ln613,ln214,ln614,ln215,ln615
		__m128i l10lowl14low_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+20*8),sse_128_load_vector_a(asrc[k]+28*8));
		__m128i l10lowl14low_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+20*8),sse_128_load_vector_a(asrc[k]+28*8));
		__m128i l10highl14high_low = sse128_unpacklo_u16(sse_128_load_vector_a(asrc[k]+21*8),sse_128_load_vector_a(asrc[k]+29*8));
		__m128i l10highl14high_high = sse128_unpackhi_u16(sse_128_load_vector_a(asrc[k]+21*8),sse_128_load_vector_a(asrc[k]+29*8));


		odd_c[0] = sse128_unpacklo_u64(l1357_low_low_low,l9111315_low_low_low);		//odd members of colum 0 in a row (10,30,50,70,90,110,130,150)
		odd_c[1] = sse128_unpackhi_u64(l1357_low_low_low,l9111315_low_low_low);		//odd members of colum 1 in a row (11,31,51,71,91,111,131,151)
		odd_c[2] = sse128_unpacklo_u64(l1357_low_low_high,l9111315_low_low_high);		//odd members of colum 2 in a row (12,32,52,72,92,112,132,152)
		odd_c[3] = sse128_unpackhi_u64(l1357_low_low_high,l9111315_low_low_high);		//odd members of colum 3 in a row (13,33,53,73,93,113,133,153)
		odd_c[4] = sse128_unpacklo_u64(l1357_low_high_low,l9111315_low_high_low);		//odd members of colum 4 in a row 
		odd_c[5] = sse128_unpackhi_u64(l1357_low_high_low,l9111315_low_high_low);		//odd members of colum 5 in a row 
		odd_c[6] = sse128_unpacklo_u64(l1357_low_high_high,l9111315_low_high_high);	//odd members of colum 6 in a row 
		odd_c[7] = sse128_unpackhi_u64(l1357_low_high_high,l9111315_low_high_high);	//odd members of colum 7 in a row 
		odd_c[8] = sse128_unpacklo_u64(l1357_high_low_low,l9111315_high_low_low);		//odd members of colum 8 in a row 
		odd_c[9] = sse128_unpackhi_u64(l1357_high_low_low,l9111315_high_low_low);		//odd members of colum 9 in a row 
		odd_c[10] = sse128_unpacklo_u64(l1357_high_low_high,l9111315_high_low_high);	//odd members of colum 10 in a row 
		odd_c[11] = sse128_unpackhi_u64(l1357_high_low_high,l9111315_high_low_high);	//odd members of colum 11 in a row 
		odd_c[12] = sse128_unpacklo_u64(l1357_high_high_low,l9111315_high_high_low);	//odd members of colum 12 in a row 
		odd_c[13] = sse128_unpackhi_u64(l1357_high_high_low,l9111315_high_high_low);	//odd members of colum 13 in a row 
		odd_c[14] = sse128_unpacklo_u64(l1357_high_high_high,l9111315_high_high_high);//odd members of colum 14 in a row 
		odd_c[15] = sse128_unpackhi_u64(l1357_high_high_high,l9111315_high_high_high);//odd members of colum 15 in a row 


		ee_c[0] = sse128_unpacklo_u32(l0lowl4low_low,l8lowl12low_low);					//ln00,ln40,ln80,ln120,ln01,ln41,ln81,ln121
		ee_c[1] = sse128_unpackhi_u32(l0lowl4low_low,l8lowl12low_low);					//ln02,ln42,ln82,ln122,ln03,ln43,ln83,ln123
		ee_c[2] = sse128_unpacklo_u32(l0lowl4low_high,l8lowl12low_high);				//ln04,ln44,ln84,ln124,ln05,ln45,ln85,ln125
		ee_c[3] = sse128_unpackhi_u32(l0lowl4low_high,l8lowl12low_high);				//ln06,ln46,ln86,ln126,ln07,ln47,ln87,ln127
		ee_c[4] = sse128_unpacklo_u32(l0highl4high_low,l8highl12high_low);				//ln08,ln48,ln88,ln128,ln09,ln49,ln89,ln129
		ee_c[5] = sse128_unpackhi_u32(l0highl4high_low,l8highl12high_low);				//ln010,ln410,ln810,ln1210,ln011,ln411,ln811,ln1211
		ee_c[6] = sse128_unpacklo_u32(l0highl4high_high,l8highl12high_high);			//ln012,ln412,ln812,ln1212,ln013,ln413,ln813,ln1213
		ee_c[7] = sse128_unpackhi_u32(l0highl4high_high,l8highl12high_high);			//ln014,ln414,ln814,ln1214,ln015,ln415,ln815,ln1215

		eo_c[0] = sse128_unpacklo_u32(l2lowl6low_low,l10lowl14low_low);					//ln20,ln60,ln100,ln140,ln21,ln61,ln101,ln141
		eo_c[1] = sse128_unpackhi_u32(l2lowl6low_low,l10lowl14low_low);					//ln22,ln62,ln102,ln142,ln23,ln63,ln103,ln143
		eo_c[2] = sse128_unpacklo_u32(l2lowl6low_high,l10lowl14low_high);				//ln24,ln64,ln104,ln144,ln25,ln65,ln105,ln145
		eo_c[3] = sse128_unpackhi_u32(l2lowl6low_high,l10lowl14low_high);				//ln26,ln66,ln106,ln146,ln27,ln67,ln107,ln147
		eo_c[4] = sse128_unpacklo_u32(l2highl6high_low,l10highl14high_low);				//ln28,ln68,ln108,ln148,ln29,ln69,ln109,ln149
		eo_c[5] = sse128_unpackhi_u32(l2highl6high_low,l10highl14high_low);				//ln210,ln610,ln1010,ln1410,ln211,ln611,ln1011,ln1411
		eo_c[6] = sse128_unpacklo_u32(l2highl6high_high,l10highl14high_high);			//ln212,ln612,ln1012,ln1412,ln213,ln613,ln1013,ln1413
		eo_c[7] = sse128_unpackhi_u32(l2highl6high_high,l10highl14high_high);			//ln214,ln614,ln1014,ln1414,ln215,ln615,ln1015,ln1415


		round = sse_128_vector_i32(iround[k]);
		shift = ishift[k];
		pdst = adst[k];
		stride = astride[k];

		for (i = 0; i < 8; i++)
		{
			__m128i aux0 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[2*i],sse_128_load_vector_a(inv_dct16x16_odd[0])), sse_128_madd_i16_i32(odd_c[2*i],sse_128_load_vector_a(inv_dct16x16_odd[1])));
			__m128i aux1 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[2*i],sse_128_load_vector_a(inv_dct16x16_odd[2])), sse_128_madd_i16_i32(odd_c[2*i],sse_128_load_vector_a(inv_dct16x16_odd[3])));
			__m128i aux2 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[2*i],sse_128_load_vector_a(inv_dct16x16_odd[4])), sse_128_madd_i16_i32(odd_c[2*i],sse_128_load_vector_a(inv_dct16x16_odd[5])));
			__m128i aux3 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[2*i],sse_128_load_vector_a(inv_dct16x16_odd[6])), sse_128_madd_i16_i32(odd_c[2*i],sse_128_load_vector_a(inv_dct16x16_odd[7])));
			O[0] = sse_128_hadd_i32(aux0, aux1);
			O[1] = sse_128_hadd_i32(aux2, aux3);
			aux0 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[2*i+1],sse_128_load_vector_a(inv_dct16x16_odd[0])), sse_128_madd_i16_i32(odd_c[2*i+1],sse_128_load_vector_a(inv_dct16x16_odd[1])));
			aux1 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[2*i+1],sse_128_load_vector_a(inv_dct16x16_odd[2])), sse_128_madd_i16_i32(odd_c[2*i+1],sse_128_load_vector_a(inv_dct16x16_odd[3])));
			aux2 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[2*i+1],sse_128_load_vector_a(inv_dct16x16_odd[4])), sse_128_madd_i16_i32(odd_c[2*i+1],sse_128_load_vector_a(inv_dct16x16_odd[5])));
			aux3 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[2*i+1],sse_128_load_vector_a(inv_dct16x16_odd[6])), sse_128_madd_i16_i32(odd_c[2*i+1],sse_128_load_vector_a(inv_dct16x16_odd[7])));
			O[2] = sse_128_hadd_i32(aux0, aux1);
			O[3] = sse_128_hadd_i32(aux2, aux3);

			aux0 = sse128_unpacklo_u64(eo_c[i], eo_c[i]); //[c2i,c2i]
			EO[0] = sse_128_hadd_i32(sse_128_madd_i16_i32(aux0,sse_128_load_vector_a(inv_dct16x16_eo[0])), sse_128_madd_i16_i32(aux0,sse_128_load_vector_a(inv_dct16x16_eo[1])));
			aux0 = sse128_unpackhi_u64(eo_c[i], eo_c[i]); //[c2i+1,c2i+1]
			EO[1] = sse_128_hadd_i32(sse_128_madd_i16_i32(aux0,sse_128_load_vector_a(inv_dct16x16_eo[0])), sse_128_madd_i16_i32(aux0,sse_128_load_vector_a(inv_dct16x16_eo[1])));

			aux0 = sse128_unpacklo_u64(ee_c[i], ee_c[i]); //[c2i,c2i]
			EE[0] = sse_128_hadd_i32(sse_128_madd_i16_i32(aux0,sse_128_load_vector_a(inv_dct16x16_ee[0])), sse_128_madd_i16_i32(aux0,sse_128_load_vector_a(inv_dct16x16_ee[1])));
			aux0 = sse128_unpackhi_u64(ee_c[i], ee_c[i]); //[c2i+1,c2i+1]
			EE[1] = sse_128_hadd_i32(sse_128_madd_i16_i32(aux0,sse_128_load_vector_a(inv_dct16x16_ee[0])), sse_128_madd_i16_i32(aux0,sse_128_load_vector_a(inv_dct16x16_ee[1])));

			E[0] = sse_128_add_i32(EE[0], EO[0]);//c2*i l
			E[1] = sse_128_sub_i32(sse_128_shuffle_32(EE[0], 0x1b), sse_128_shuffle_32(EO[0], 0x1b));	//c2*i h
			E[2] = sse_128_add_i32(EE[1], EO[1]);//c2*i+1 l
			E[3] = sse_128_sub_i32(sse_128_shuffle_32(EE[1], 0x1b), sse_128_shuffle_32(EO[1], 0x1b));	//c2*i+1 h


			aux0 = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E[0], O[0]), round), shift), sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E[1], O[1]), round), shift));	
			aux1 = sse_128_shuffle_8(sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E[0], O[0]), round), shift), sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E[1], O[1]), round), shift)), shuff_mask_0);	
			aux2 = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E[2], O[2]), round), shift), sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E[3], O[3]), round), shift));	
			aux3 = sse_128_shuffle_8(sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E[2], O[2]), round), shift), sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E[3], O[3]), round), shift)), shuff_mask_0);		

			sse_128_store_vector_a(pdst, aux0);
			sse_128_store_vector_a(pdst+8, aux1);
			pdst+=stride;
			sse_128_store_vector_a(pdst, aux2);
			sse_128_store_vector_a(pdst+8, aux3);
			pdst+=stride;
		}
	}
}




//---------------------------------------DCT-32x32 ----------------------------------------------


//EEEE and EEEO (2 coeffs per line). Lines 0,8,16,24
ALIGN(16) const int16_t dct32x32_butterfly_coeffs_0[4][8] =
{
	{ 64, 64, 64, 64, 64, 64, 64, 64},//0
	{ 83, 36, 83, 36, 83, 36, 83, 36},//8
	{ 64,-64, 64,-64, 64,-64, 64,-64},//16
	{ 36,-83, 36,-83, 36,-83, 36,-83}//24
};

//EEO (4 coeffs per line). Lines 4,12,20,28
ALIGN(16) const int16_t dct32x32_butterfly_coeffs_1[4][8] =
{
	{ 89, 75, 50, 18, 89, 75, 50, 18},//4
	{ 75,-18,-89,-50, 75,-18,-89,-50},//12
	{ 50,-89, 18, 75, 50,-89, 18, 75},//20
	{ 18,-50, 75,-89, 18,-50, 75,-89}//28
};

//EO (8 coeffs per line). Lines 2,6,10,14,18,22,26,30
ALIGN(16) const int16_t dct32x32_butterfly_coeffs_2[8][8] =
{
	{ 90, 87, 80, 70, 57, 43, 25,  9},//2
	{ 87, 57,  9,-43,-80,-90,-70,-25},//6
	{ 80,  9,-70,-87,-25, 57, 90, 43},//10
	{ 70,-43,-87,  9, 90, 25,-80,-57},//14
	{ 57,-80,-25, 90, -9,-87, 43, 70},//18
	{ 43,-90, 57, 25,-87, 70,  9,-80},//22
	{ 25,-70, 90,-80, 43,  9,-57, 87},//26
	{  9,-25, 43,-57, 70,-80, 87,-90},//30
};


//O (16 coeffs per line). Lines 1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31
ALIGN(16) const int16_t dct32x32_butterfly_coeffs_3[16][16] =
{
	{ 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13,  4},//1
	{ 90, 82, 67, 46, 22, -4,-31,-54,-73,-85,-90,-88,-78,-61,-38,-13},//3
	{ 88, 67, 31,-13,-54,-82,-90,-78,-46, -4, 38, 73, 90, 85, 61, 22},//5
	{ 85, 46,-13,-67,-90,-73,-22, 38, 82, 88, 54, -4,-61,-90,-78,-31},//7
	{ 82, 22,-54,-90,-61, 13, 78, 85, 31,-46,-90,-67,  4, 73, 88, 38},//9
	{ 78, -4,-82,-73, 13, 85, 67,-22,-88,-61, 31, 90, 54,-38,-90,-46},//11
	{ 73,-31,-90,-22, 78, 67,-38,-90,-13, 82, 61,-46,-88, -4, 85, 54},//13
	{ 67,-54,-78, 38, 85,-22,-90,  4, 90, 13,-88,-31, 82, 46,-73,-61},//15
	{ 61,-73,-46, 82, 31,-88,-13, 90, -4,-90, 22, 85,-38,-78, 54, 67},//17
	{ 54,-85, -4, 88,-46,-61, 82, 13,-90, 38, 67,-78,-22, 90,-31,-73},//19
	{ 46,-90, 38, 54,-90, 31, 61,-88, 22, 67,-85, 13, 73,-82,  4, 78},//21
	{ 38,-88, 73, -4,-67, 90,-46,-31, 85,-78, 13, 61,-90, 54, 22,-82},//23
	{ 31,-78, 90,-61,  4, 54,-88, 82,-38,-22, 73,-90, 67,-13,-46, 85},//25
	{ 22,-61, 85,-90, 73,-38, -4, 46,-78, 90,-82, 54,-13,-31, 67,-88},//27
	{ 13,-38, 61,-78, 88,-90, 85,-73, 54,-31,  4, 22,-46, 67,-82, 90},//29
	{  4,-13, 22,-31, 38,-46, 54,-61, 67,-73, 78,-82, 85,-88, 90,-90},//31
};


//EEEE and EEEO (2 coeffs per line). Lines 0,8,16,24
ALIGN(16) const int32_t dct32x32_butterfly_coeffs_4[][4] =
{
	{ 64, 64, 64, 64},//0
	{ 83, 36, 83, 36},//8
	{ 64,-64, 64,-64},//16
	{ 36,-83, 36,-83}//24
};

//EEO (4 coeffs per line). Lines 4,12,20,28
ALIGN(16) const int32_t dct32x32_butterfly_coeffs_5[][4] =
{
	{ 89, 75, 50, 18},//4
	{ 75,-18,-89,-50},//12
	{ 50,-89, 18, 75},//20
	{ 18,-50, 75,-89}//28
};

//EO (8 coeffs per line). Lines 2,6,10,14,18,22,26,30
ALIGN(16) const int32_t dct32x32_butterfly_coeffs_6[][8] =
{
	{ 90, 87, 80, 70, 57, 43, 25,  9},//2
	{ 87, 57,  9,-43,-80,-90,-70,-25},//6
	{ 80,  9,-70,-87,-25, 57, 90, 43},//10
	{ 70,-43,-87,  9, 90, 25,-80,-57},//14
	{ 57,-80,-25, 90, -9,-87, 43, 70},//18
	{ 43,-90, 57, 25,-87, 70,  9,-80},//22
	{ 25,-70, 90,-80, 43,  9,-57, 87},//26
	{  9,-25, 43,-57, 70,-80, 87,-90},//30
};


//O (16 coeffs per line). Lines 1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31
ALIGN(16) const int32_t dct32x32_butterfly_coeffs_7[16][16] =
{
	{ 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13,  4},//1
	{ 90, 82, 67, 46, 22, -4,-31,-54,-73,-85,-90,-88,-78,-61,-38,-13},//3
	{ 88, 67, 31,-13,-54,-82,-90,-78,-46, -4, 38, 73, 90, 85, 61, 22},//5
	{ 85, 46,-13,-67,-90,-73,-22, 38, 82, 88, 54, -4,-61,-90,-78,-31},//7
	{ 82, 22,-54,-90,-61, 13, 78, 85, 31,-46,-90,-67,  4, 73, 88, 38},//9
	{ 78, -4,-82,-73, 13, 85, 67,-22,-88,-61, 31, 90, 54,-38,-90,-46},//11
	{ 73,-31,-90,-22, 78, 67,-38,-90,-13, 82, 61,-46,-88, -4, 85, 54},//13
	{ 67,-54,-78, 38, 85,-22,-90,  4, 90, 13,-88,-31, 82, 46,-73,-61},//15
	{ 61,-73,-46, 82, 31,-88,-13, 90, -4,-90, 22, 85,-38,-78, 54, 67},//17
	{ 54,-85, -4, 88,-46,-61, 82, 13,-90, 38, 67,-78,-22, 90,-31,-73},//19
	{ 46,-90, 38, 54,-90, 31, 61,-88, 22, 67,-85, 13, 73,-82,  4, 78},//21
	{ 38,-88, 73, -4,-67, 90,-46,-31, 85,-78, 13, 61,-90, 54, 22,-82},//23
	{ 31,-78, 90,-61,  4, 54,-88, 82,-38,-22, 73,-90, 67,-13,-46, 85},//25
	{ 22,-61, 85,-90, 73,-38, -4, 46,-78, 90,-82, 54,-13,-31, 67,-88},//27
	{ 13,-38, 61,-78, 88,-90, 85,-73, 54,-31,  4, 22,-46, 67,-82, 90},//29
	{  4,-13, 22,-31, 38,-46, 54,-61, 67,-73, 78,-82, 85,-88, 90,-90},//31
};


#define DCT32x32_EXTRACT_1_LINE_16BITS(src, E, O)									\
			aux0 = sse_128_load_vector_a(psrc);										\
			aux1 = sse_128_shuffle_8(sse_128_load_vector_a(psrc+24), shuff_mask_0);		\
			E[0] = sse_128_add_i16(aux0, aux1);	/*E[0..7]*/							\
			O[0] = sse_128_sub_i16(aux0, aux1);	/*O[0..7]*/							\
			aux0 = sse_128_load_vector_a(psrc+8);									\
			aux1 = sse_128_shuffle_8(sse_128_load_vector_a(psrc+16), shuff_mask_0);		\
			E[1] = sse_128_add_i16(aux0, aux1);/*E[8..15]*/							\
			O[1] = sse_128_sub_i16(aux0, aux1);/*O[8..15]*/							


#define DCT32x32_EXTRACT_1_LINE_32BITS(src, E, O)																						\
			aux0 = sse_128_load_vector_a(psrc);																							\
			aux1 = sse_128_shuffle_8(sse_128_load_vector_a(psrc+24), shuff_mask_0);/*0,1,2,3,4,5,6,7 -> 7,6,5,4,3,2,1,0 - 16 bit data*/		\
			aux2 = sse_128_convert_i16_i32(aux0);																							\
			aux3 = sse_128_convert_i16_i32(aux1);																							\
			E[0] = sse_128_add_i32(aux2, aux3);/*E[0..3]*/																				\
			O[0] = sse_128_sub_i32(aux2, aux3);/*O[0..3]*/																				\
			aux2 = sse_128_convert_i16_i32(sse128_unpackhi_u64(aux0, aux0));																	\
			aux3 = sse_128_convert_i16_i32(sse128_unpackhi_u64(aux1, aux1));																	\
			E[1] = sse_128_add_i32(aux2, aux3);/*E[4..7]*/																				\
			O[1] = sse_128_sub_i32(aux2, aux3);/*O[4..7]*/																				\
			aux0 = sse_128_load_vector_a(psrc+8);																						\
			aux1 = sse_128_shuffle_8(sse_128_load_vector_a(psrc+16), shuff_mask_0);/*0,1,2,3,4,5,6,7 -> 7,6,5,4,3,2,1,0 - 16 bit data*/		\
			aux2 = sse_128_convert_i16_i32(aux0);																							\
			aux3 = sse_128_convert_i16_i32(aux1);																							\
			E[2] = sse_128_add_i32(aux2, aux3);/*E[8..11]*/																				\
			O[2] = sse_128_sub_i32(aux2, aux3);/*O[8..11]*/																				\
			aux2 = sse_128_convert_i16_i32(sse128_unpackhi_u64(aux0, aux0));																	\
			aux3 = sse_128_convert_i16_i32(sse128_unpackhi_u64(aux1, aux1));																	\
			E[3] = sse_128_add_i32(aux2, aux3);/*E[12..15]*/																			\
			O[3] = sse_128_sub_i32(aux2, aux3);/*O[12..15]*/																		


void sse_aligned_dct_32x32(int16_t *src, int16_t *dst, int stride, int16_t *aux) 
{
	int i, j;
	__m128_	zero = sse_128_zero_vector();
	__m128i round[2];

	__m128i O[8][4], E[4];
	__m128i EO[8][2], EE[2];
	__m128i EEO[4][2], EEE[4];
	__m128i EEEO[4], EEEE[4];
	__m128i  aux1, aux0, aux2, aux3;

	__m128i shuff_mask_0 = sse_128_load_vector_a(shuffle_mask_dct_16_0);//0,1,2,3,4,5,6,7 -> 7,6,5,4,3,2,1,0
	__m128i shuff_mask_1 = sse_128_load_vector_a(shuffle_mask_dct_16_1);//0,1,2,3,4,5,6,7 -> 3,2,1,0,7,6,5,4
	__m128i shuff_mask_2 = sse_128_load_vector_a(shuffle_mask_dct_16_2);//0,1,2,3,4,5,6,7 -> 1,0,3,2,5,4,7,6
	int16_t *psrc, *pdst;

	round[0] = sse_128_vector_i32(8);
	round[1] = sse_128_vector_i32(1024);

	//1st pass with 16bit precission
	for (j=0; j<4; j++)//8 lines each 
	{
		psrc = src+j*8*stride;
		pdst = aux+j*8;

		for (i=0; i<2; i++)
		{
			//line 1 of 4
			DCT32x32_EXTRACT_1_LINE_16BITS(psrc, E, O[i*4])
			aux1 = sse_128_shuffle_8(E[1], shuff_mask_0);
			EE[0] = sse_128_add_i16(E[0], aux1);
			EO[2*i][0] = sse_128_sub_i16(E[0], aux1);

			//line 2 of 4
			psrc += stride;
			DCT32x32_EXTRACT_1_LINE_16BITS(psrc, E, O[i*4+1])

			aux1 = sse_128_shuffle_8(E[1], shuff_mask_0);
			EE[1] = sse_128_add_i16(E[0], aux1);
			EO[2*i][1] = sse_128_sub_i16(E[0], aux1);

			aux0 = sse128_unpacklo_u64(EE[0],EE[1]);
			aux1 = sse_128_shuffle_8(sse128_unpackhi_u64(EE[0],EE[1]), shuff_mask_1);//0,1,2,3,4,5,6,7 -> 3,2,1,0,7,6,5,4

			EEE[0] = sse_128_add_i16(aux0, aux1);//2 lines, 4 coeffs each. 
			EEO[i][0] = sse_128_sub_i16(aux0, aux1);//2 lines, 4 coeffs each

			//line 3 of 4
			psrc += stride;
			DCT32x32_EXTRACT_1_LINE_16BITS(psrc, E, O[i*4+2])

			aux1 = sse_128_shuffle_8(E[1], shuff_mask_0);
			EE[0] = sse_128_add_i16(E[0], aux1);
			EO[2*i+1][0] = sse_128_sub_i16(E[0], aux1);

			//line 4 of 4
			psrc += stride;
			DCT32x32_EXTRACT_1_LINE_16BITS(psrc, E, O[i*4+3])

			aux1 = sse_128_shuffle_8(E[1], shuff_mask_0);
			EE[1] = sse_128_add_i16(E[0], aux1);
			EO[2*i+1][1] = sse_128_sub_i16(E[0], aux1);

			aux0 = sse128_unpacklo_u64(EE[0],EE[1]);
			aux1 = sse_128_shuffle_8(sse128_unpackhi_u64(EE[0],EE[1]), shuff_mask_1);//0,1,2,3,4,5,6,7 -> 3,2,1,0,7,6,5,4

			EEE[1] = sse_128_add_i16(aux0, aux1);//2 lines, 4 coeffs each
			EEO[i][1] = sse_128_sub_i16(aux0, aux1);//2 lines, 4 coeffs each

			aux0 = sse128_unpacklo_u32(EEE[0],EEE[1]);
			aux1 = sse128_unpackhi_u32(EEE[0],EEE[1]);
			aux2 = sse128_unpacklo_u32(aux0,aux1);
			aux3 = sse_128_shuffle_8(sse128_unpackhi_u32(aux0,aux1), shuff_mask_2);//0,1,2,3,4,5,6,7 -> 1,0,3,2,5,4,7,6

			EEEE[i] = sse_128_add_i16(aux2, aux3);
			EEEO[i] = sse_128_sub_i16(aux2, aux3);
			psrc += stride;
		}

		//lines 0,8,16,32
		for (i = 0; i < 2; i++)
		{
			__m128i  dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_madd_i16_i32(EEEE[0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_0[2*i])), round[0]), 4);
			__m128i  dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_madd_i16_i32(EEEE[1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_0[2*i])), round[0]), 4);
			sse_128_store_vector_a(pdst+(16*i)*32, sse128_packs_u32_u16(dstl, dsth));
			dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_madd_i16_i32(EEEO[0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_0[2*i+1])), round[0]), 4);
			dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_madd_i16_i32(EEEO[1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_0[2*i+1])), round[0]), 4);
			sse_128_store_vector_a(pdst+(16*i+8)*32, sse128_packs_u32_u16(dstl, dsth));
		}

		//lines 4,12,20,28
		for (i = 0; i < 4; i++)
		{
			__m128i  dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_madd_i16_i32(EEO[0][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_1[i])), sse_128_madd_i16_i32(EEO[0][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_1[i]))), round[0]), 4);
			__m128i  dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_madd_i16_i32(EEO[1][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_1[i])), sse_128_madd_i16_i32(EEO[1][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_1[i]))), round[0]), 4);
			sse_128_store_vector_a(pdst+(8*i+4)*32, sse128_packs_u32_u16(dstl, dsth));
		}

		//lines 2,6,10,14,18,22,26,30
		for (i = 0; i < 8; i++)
		{
			__m128i  dst0 = sse_128_hadd_i32(sse_128_madd_i16_i32(EO[0][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_2[i])), sse_128_madd_i16_i32(EO[0][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_2[i])));
			__m128i  dst1 = sse_128_hadd_i32(sse_128_madd_i16_i32(EO[1][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_2[i])), sse_128_madd_i16_i32(EO[1][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_2[i])));
			__m128i  dsth, dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(dst0, dst1), round[0]), 4);
			dst0 = sse_128_hadd_i32(sse_128_madd_i16_i32(EO[2][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_2[i])), sse_128_madd_i16_i32(EO[2][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_2[i])));
			dst1 = sse_128_hadd_i32(sse_128_madd_i16_i32(EO[3][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_2[i])), sse_128_madd_i16_i32(EO[3][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_2[i])));
			dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(dst0, dst1), round[0]), 4);
			sse_128_store_vector_a(pdst+(4*i+2)*32, sse128_packs_u32_u16(dstl, dsth));
		}

		for (i = 0; i < 16; i++)
		{
			__m128i  dst0 = sse_128_hadd_i32(sse_128_madd_i16_i32(O[0][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i])), sse_128_madd_i16_i32(O[0][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i]+8)));
			__m128i  dst1 = sse_128_hadd_i32(sse_128_madd_i16_i32(O[1][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i])), sse_128_madd_i16_i32(O[1][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i]+8)));
			__m128i  dst2 = sse_128_hadd_i32(sse_128_madd_i16_i32(O[2][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i])), sse_128_madd_i16_i32(O[2][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i]+8)));
			__m128i  dst3 = sse_128_hadd_i32(sse_128_madd_i16_i32(O[3][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i])), sse_128_madd_i16_i32(O[3][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i]+8)));
			__m128i  dsth, dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_hadd_i32(dst0, dst1), sse_128_hadd_i32(dst2, dst3)), round[0]), 4);
			dst0 = sse_128_hadd_i32(sse_128_madd_i16_i32(O[4][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i])), sse_128_madd_i16_i32(O[4][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i]+8)));
			dst1 = sse_128_hadd_i32(sse_128_madd_i16_i32(O[5][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i])), sse_128_madd_i16_i32(O[5][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i]+8)));
			dst2 = sse_128_hadd_i32(sse_128_madd_i16_i32(O[6][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i])), sse_128_madd_i16_i32(O[6][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i]+8)));
			dst3 = sse_128_hadd_i32(sse_128_madd_i16_i32(O[7][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i])), sse_128_madd_i16_i32(O[7][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_3[i]+8)));
			dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_hadd_i32(dst0, dst1), sse_128_hadd_i32(dst2, dst3)), round[0]), 4);
			sse_128_store_vector_a(pdst+(2*i+1)*32, sse128_packs_u32_u16(dstl, dsth));
		}
	}


	shuff_mask_1 = sse_128_load_vector_a(shuffle_mask_dct_32_0);//0,1,2,3 -> 3,2,1,0
	shuff_mask_2 = sse_128_load_vector_a(shuffle_mask_dct_32_1);//0,1,2,3 -> 1,0,3,2

	//2nd pass with 32bit precission
	for (j=0; j<4; j++)//8 lines each 
	{
		psrc = aux+j*8*32;
		pdst = dst+j*8;
		
		for (i=0; i<2; i++)//2 iterations = 8 lines
		{
			//line 1 of 4
			DCT32x32_EXTRACT_1_LINE_32BITS(psrc, E, O[4*i])
			aux1 = sse_128_shuffle_8(E[3], shuff_mask_1);
			EE[0] = sse_128_add_i32(E[0], aux1);		//[0..3]
			EO[4*i][0] = sse_128_sub_i32(E[0], aux1);	//[0..3]
			aux1 = sse_128_shuffle_8(E[2], shuff_mask_1);
			EE[1] = sse_128_add_i32(E[1], aux1);		//[4..7]
			EO[4*i][1] = sse_128_sub_i32(E[1], aux1);	//[4..7]

			aux0 = sse_128_shuffle_8(EE[1],shuff_mask_1);//0,1,2,3 -> 3,2,1,0
			EEE[0] = sse_128_add_i32(EE[0], aux0);//1 line, 4 coeffs each. 
			EEO[2*i][0] = sse_128_sub_i32(EE[0], aux0);//1 line, 4 coeffs each

			//line 2 of 4
			psrc += 32;
			DCT32x32_EXTRACT_1_LINE_32BITS(psrc, E, O[4*i+1])
			aux1 = sse_128_shuffle_8(E[3], shuff_mask_1);
			EE[0] = sse_128_add_i32(E[0], aux1);
			EO[4*i+1][0] = sse_128_sub_i32(E[0], aux1);//1 line, 8 coeffs in E=[0,1]
			aux1 = sse_128_shuffle_8(E[2], shuff_mask_1);
			EE[1] = sse_128_add_i32(E[1], aux1);
			EO[4*i+1][1] = sse_128_sub_i32(E[1], aux1);

			aux0 = sse_128_shuffle_8(EE[1],shuff_mask_1);//0,1,2,3 -> 3,2,1,0
			EEE[1] = sse_128_add_i32(EE[0], aux0);//1 line, 4 coeffs each. 
			EEO[2*i][1] = sse_128_sub_i32(EE[0], aux0);//1 line, 4 coeffs each

			aux0 = sse128_unpacklo_u64(EEE[0],EEE[1]);
			aux1 = sse_128_shuffle_8(sse128_unpackhi_u64(EEE[0],EEE[1]), shuff_mask_2);//0,1,2,3 -> 1,0,3,2

			EEEE[2*i] = sse_128_add_i32(aux0, aux1);//2 lines, 2 coeffs each
			EEEO[2*i] = sse_128_sub_i32(aux0, aux1);//2 lines, 2 coeffs each

			//line 3 of 4
			psrc += 32;
			DCT32x32_EXTRACT_1_LINE_32BITS(psrc, E, O[4*i+2])
			aux1 = sse_128_shuffle_8(E[3], shuff_mask_1);
			EE[0] = sse_128_add_i32(E[0], aux1);		//[0..3]
			EO[4*i+2][0] = sse_128_sub_i32(E[0], aux1);	//[0..3]
			aux1 = sse_128_shuffle_8(E[2], shuff_mask_1);
			EE[1] = sse_128_add_i32(E[1], aux1);		//[4..7]
			EO[4*i+2][1] = sse_128_sub_i32(E[1], aux1);	//[4..7]

			aux0 = sse_128_shuffle_8(EE[1],shuff_mask_1);//0,1,2,3 -> 3,2,1,0
			EEE[0] = sse_128_add_i32(EE[0], aux0);//1 line, 4 coeffs each. 
			EEO[2*i+1][0] = sse_128_sub_i32(EE[0], aux0);//1 line, 4 coeffs each

			//line 4 of 4
			psrc += 32;
			DCT32x32_EXTRACT_1_LINE_32BITS(psrc, E, O[4*i+3])
			aux1 = sse_128_shuffle_8(E[3], shuff_mask_1);
			EE[0] = sse_128_add_i32(E[0], aux1);
			EO[4*i+3][0] = sse_128_sub_i32(E[0], aux1);//1 line, 8 coeffs in E=[0,1]
			aux1 = sse_128_shuffle_8(E[2], shuff_mask_1);
			EE[1] = sse_128_add_i32(E[1], aux1);
			EO[4*i+3][1] = sse_128_sub_i32(E[1], aux1);

			aux0 = sse_128_shuffle_8(EE[1],shuff_mask_1);//0,1,2,3 -> 3,2,1,0
			EEE[1] = sse_128_add_i32(EE[0], aux0);//1 line, 4 coeffs each. 
			EEO[2*i+1][1] = sse_128_sub_i32(EE[0], aux0);//1 line, 4 coeffs each

			aux0 = sse128_unpacklo_u64(EEE[0],EEE[1]);
			aux1 = sse_128_shuffle_8(sse128_unpackhi_u64(EEE[0],EEE[1]), shuff_mask_2);//0,1,2,3 -> 1,0,3,2

			EEEE[2*i+1] = sse_128_add_i32(aux0, aux1);//2 lines, 2 coeffs each
			EEEO[2*i+1] = sse_128_sub_i32(aux0, aux1);//2 lines, 2 coeffs each

			psrc += 32;
		}


		//lines 0,8,16,32
		for (i = 0; i < 2; i++)
		{
			__m128i  dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_mul_i32(EEEE[0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_4[2*i])), sse_128_mul_i32(EEEE[1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_4[2*i]))), round[1]), 11);
			__m128i  dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_mul_i32(EEEE[2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_4[2*i])), sse_128_mul_i32(EEEE[3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_4[2*i]))), round[1]), 11);
			sse_128_store_vector_a(pdst+(16*i)*32, sse128_packs_u32_u16(dstl, dsth));

			dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_mul_i32(EEEO[0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_4[2*i+1])), sse_128_mul_i32(EEEO[1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_4[2*i+1]))), round[1]), 11);
			dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_mul_i32(EEEO[2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_4[2*i+1])), sse_128_mul_i32(EEEO[3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_4[2*i+1]))), round[1]), 11);
			sse_128_store_vector_a(pdst+(16*i+8)*32, sse128_packs_u32_u16(dstl, dsth));
		}


		//voy por aqui - este bucle esta cambiado pero no validado
		//lines 4,12,20,28
		for (i = 0; i < 4; i++)
		{
			__m128i  dst0 = sse_128_hadd_i32(sse_128_mul_i32(EEO[0][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_5[i])), sse_128_mul_i32(EEO[0][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_5[i])));
			__m128i  dst1 = sse_128_hadd_i32(sse_128_mul_i32(EEO[1][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_5[i])), sse_128_mul_i32(EEO[1][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_5[i])));
			__m128i  dst2 = sse_128_hadd_i32(sse_128_mul_i32(EEO[2][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_5[i])), sse_128_mul_i32(EEO[2][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_5[i])));
			__m128i  dst3 = sse_128_hadd_i32(sse_128_mul_i32(EEO[3][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_5[i])), sse_128_mul_i32(EEO[3][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_5[i])));
			__m128i  dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(dst0, dst1), round[1]), 11);
			__m128i  dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(dst2, dst3), round[1]), 11);

			sse_128_store_vector_a(pdst+(8*i+4)*32, sse128_packs_u32_u16(dstl, dsth));
		}

		//lines 2,6,10,14,18,22,26,30
		for (i = 0; i < 8; i++)
		{
			__m128i  dst0 = sse_128_hadd_i32(sse_128_mul_i32(EO[0][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i])), sse_128_mul_i32(EO[0][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i]+4)));
			 __m128i  dst1 = sse_128_hadd_i32(sse_128_mul_i32(EO[1][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i])), sse_128_mul_i32(EO[1][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i]+4)));
			__m128i  dst2 = sse_128_hadd_i32(sse_128_mul_i32(EO[2][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i])), sse_128_mul_i32(EO[2][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i]+4)));
			__m128i  dst3 = sse_128_hadd_i32(sse_128_mul_i32(EO[3][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i])), sse_128_mul_i32(EO[3][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i]+4)));
			__m128i  dst4 = sse_128_hadd_i32(sse_128_mul_i32(EO[4][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i])), sse_128_mul_i32(EO[4][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i]+4)));
			__m128i  dst5 = sse_128_hadd_i32(sse_128_mul_i32(EO[5][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i])), sse_128_mul_i32(EO[5][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i]+4)));
			__m128i  dst6 = sse_128_hadd_i32(sse_128_mul_i32(EO[6][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i])), sse_128_mul_i32(EO[6][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i]+4)));
			__m128i  dst7 = sse_128_hadd_i32(sse_128_mul_i32(EO[7][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i])), sse_128_mul_i32(EO[7][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_6[i]+4)));

			__m128i  dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_hadd_i32(dst0, dst1), sse_128_hadd_i32(dst2, dst3)), round[1]), 11);
			__m128i  dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_hadd_i32(dst4, dst5), sse_128_hadd_i32(dst6, dst7)), round[1]), 11);
			sse_128_store_vector_a(pdst+(4*i+2)*32, sse128_packs_u32_u16(dstl, dsth));
		}

		for (i = 0; i < 16; i++)
		{
			__m128i  dst0 = sse_128_hadd_i32(sse_128_hadd_i32(sse_128_mul_i32(O[0][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i])), sse_128_mul_i32(O[0][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+4))), sse_128_hadd_i32(sse_128_mul_i32(O[0][2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+8)), sse_128_mul_i32(O[0][3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+12))));
			__m128i  dst1 = sse_128_hadd_i32(sse_128_hadd_i32(sse_128_mul_i32(O[1][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i])), sse_128_mul_i32(O[1][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+4))), sse_128_hadd_i32(sse_128_mul_i32(O[1][2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+8)), sse_128_mul_i32(O[1][3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+12))));
			__m128i  dst2 = sse_128_hadd_i32(sse_128_hadd_i32(sse_128_mul_i32(O[2][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i])), sse_128_mul_i32(O[2][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+4))), sse_128_hadd_i32(sse_128_mul_i32(O[2][2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+8)), sse_128_mul_i32(O[2][3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+12))));
			__m128i  dst3 = sse_128_hadd_i32(sse_128_hadd_i32(sse_128_mul_i32(O[3][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i])), sse_128_mul_i32(O[3][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+4))), sse_128_hadd_i32(sse_128_mul_i32(O[3][2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+8)), sse_128_mul_i32(O[3][3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+12))));
			__m128i  dst4 = sse_128_hadd_i32(sse_128_hadd_i32(sse_128_mul_i32(O[4][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i])), sse_128_mul_i32(O[4][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+4))), sse_128_hadd_i32(sse_128_mul_i32(O[4][2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+8)), sse_128_mul_i32(O[4][3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+12))));
			__m128i  dst5 = sse_128_hadd_i32(sse_128_hadd_i32(sse_128_mul_i32(O[5][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i])), sse_128_mul_i32(O[5][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+4))), sse_128_hadd_i32(sse_128_mul_i32(O[5][2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+8)), sse_128_mul_i32(O[5][3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+12))));
			__m128i  dst6 = sse_128_hadd_i32(sse_128_hadd_i32(sse_128_mul_i32(O[6][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i])), sse_128_mul_i32(O[6][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+4))), sse_128_hadd_i32(sse_128_mul_i32(O[6][2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+8)), sse_128_mul_i32(O[6][3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+12))));
			__m128i  dst7 = sse_128_hadd_i32(sse_128_hadd_i32(sse_128_mul_i32(O[7][0], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i])), sse_128_mul_i32(O[7][1], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+4))), sse_128_hadd_i32(sse_128_mul_i32(O[7][2], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+8)), sse_128_mul_i32(O[7][3], sse_128_load_vector_a(dct32x32_butterfly_coeffs_7[i]+12))));

			__m128i  dstl = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_hadd_i32(dst0, dst1), sse_128_hadd_i32(dst2, dst3)), round[1]), 11);
			__m128i  dsth = sse_128_shift_r_i32(sse_128_add_i32(sse_128_hadd_i32(sse_128_hadd_i32(dst4, dst5), sse_128_hadd_i32(dst6, dst7)), round[1]), 11);
			sse_128_store_vector_a(pdst+(2*i+1)*32, sse128_packs_u32_u16(dstl, dsth));
		}
	}
}



ALIGN(16) const int16_t inv_dct32x32_odd[16][16] =
{
  { 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13,  4},//0
  { 90, 82, 67, 46, 22, -4,-31,-54,-73,-85,-90,-88,-78,-61,-38,-13},//1
  { 88, 67, 31,-13,-54,-82,-90,-78,-46, -4, 38, 73, 90, 85, 61, 22},//2
  { 85, 46,-13,-67,-90,-73,-22, 38, 82, 88, 54, -4,-61,-90,-78,-31},//3
  { 82, 22,-54,-90,-61, 13, 78, 85, 31,-46,-90,-67,  4, 73, 88, 38},//4
  { 78, -4,-82,-73, 13, 85, 67,-22,-88,-61, 31, 90, 54,-38,-90,-46},//5
  { 73,-31,-90,-22, 78, 67,-38,-90,-13, 82, 61,-46,-88, -4, 85, 54},//6
  { 67,-54,-78, 38, 85,-22,-90,  4, 90, 13,-88,-31, 82, 46,-73,-61},//7
  { 61,-73,-46, 82, 31,-88,-13, 90, -4,-90, 22, 85,-38,-78, 54, 67},//8
  { 54,-85, -4, 88,-46,-61, 82, 13,-90, 38, 67,-78,-22, 90,-31,-73},//9
  { 46,-90, 38, 54,-90, 31, 61,-88, 22, 67,-85, 13, 73,-82,  4, 78},//10
  { 38,-88, 73, -4,-67, 90,-46,-31, 85,-78, 13, 61,-90, 54, 22,-82},//11
  { 31,-78, 90,-61,  4, 54,-88, 82,-38,-22, 73,-90, 67,-13,-46, 85},//12
  { 22,-61, 85,-90, 73,-38, -4, 46,-78, 90,-82, 54,-13,-31, 67,-88},//13
  { 13,-38, 61,-78, 88,-90, 85,-73, 54,-31,  4, 22,-46, 67,-82, 90},//14
  {  4,-13, 22,-31, 38,-46, 54,-61, 67,-73, 78,-82, 85,-88, 90,-90},//15
};

ALIGN(16) const int16_t inv_dct32x32_eo[8][8] =
{
  { 90, 87, 80, 70, 57, 43, 25,  9},
  { 87, 57,  9,-43,-80,-90,-70,-25},
  { 80,  9,-70,-87,-25, 57, 90, 43},
  { 70,-43,-87,  9, 90, 25,-80,-57},
  { 57,-80,-25, 90, -9,-87, 43, 70},
  { 43,-90, 57, 25,-87, 70,  9,-80},
  { 25,-70, 90,-80, 43,  9,-57, 87},
  {  9,-25, 43,-57, 70,-80, 87,-90}
};

//2 columns per line
ALIGN(16) const int16_t inv_dct32x32_eeo[2][8] =
{
  { 89, 75, 50, 18, 75,-18,-89,-50},
  { 50,-89, 18, 75, 18,-50, 75,-89}
//  {-18, 50,-75, 89,-50, 89,-18,-75},
//  {-75, 18, 89, 50,-89,-75,-50,-18}
};

//2 columns per line
ALIGN(16) const int16_t inv_dct32x32_eee[2][8] =
{
  { 64, 83, 64, 36, 64, 36,-64,-83},
  { 64,-36,-64, 83, 64,-83, 64,-36}
};


void sse_aligned_inv_dct_32x32(int16_t *src, int16_t *dst, int dst_stride, int16_t *aux) 
{
	int i, j, k, shift;
	int iround[2] = {64, 2048};
	int ishift[2] = {7, 12};
	__m128i round;
	__m128i O[4], E[4], EO[2], EE[2], EEO, EEE; 

	__m128i odd_c[8][2];			//16 samples per column. one column = 2 registers
	__m128i eo_c[8];		//8 samples per column. one column per register
	__m128i eeo_c[8];				//4 samples per column. duplicated column per register
	__m128i eee_c[8];				//4 samples per column. duplicated column per register
//	__m128i eeeo_c[2], eeee_c[2];	//2 samples per column. 4 columns per register
	__m128i shuff_mask_0 = sse_128_load_vector_a(shuffle_mask_dct_16_0);//0,1,2,3,4,5,6,7 -> 7,6,5,4,3,2,1,0

	int stride; 
	int astride[2] = {32, dst_stride};
	int16_t *psrc, *pdst;
	int16_t *asrc[2] = {src, aux};
	int16_t *adst[2] = {aux, dst};

	for (k = 0; k < 2; k++)
	{
		stride  = astride[k];
		for (i = 0; i < 4; i++)//eight samples of each line per iteration
		{
			__m128i l2l6low, l2l6high, l10l14low, l10l14high, l18l22low, l18l22high, l26l30low, l26l30high;
			__m128i l261014low_low, l261014low_high, l261014high_low, l261014high_high, l18222630low_low, l18222630low_high, l18222630high_low, l18222630high_high;
			__m128i l4l12low, l4l12high, l20l28low, l20l28high, _aux;
			__m128i l0l8low, l0l8high, l16l24low, l16l24high;

			psrc = asrc[k]+i*8;
			pdst = adst[k]+i*8*stride;
			for (j = 0; j < 2; j++)//eight samples of each line per iteration
			{
				__m128i l1l3low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+1*32),sse_128_load_vector_a(psrc+3*32));				//ln10,ln30,ln11,ln31,ln12,ln32,ln13,ln33 
				__m128i l1l3high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+1*32),sse_128_load_vector_a(psrc+3*32));				//ln14,ln34,ln15,ln35,ln16,ln36,ln17,ln37 
				__m128i l5l7low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+5*32),sse_128_load_vector_a(psrc+7*32));				//ln50,ln70,ln51,ln71,ln52,ln72,ln53,ln73 
				__m128i l5l7high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+5*32),sse_128_load_vector_a(psrc+7*32));				//ln54,ln74,ln55,ln75,ln56,ln76,ln57,ln77 
				__m128i l9l11low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+9*32),sse_128_load_vector_a(psrc+11*32));				//ln90,ln110,ln91,ln111,ln92,ln112,ln93,ln113 
				__m128i l9l11high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+9*32),sse_128_load_vector_a(psrc+11*32));			//ln94,ln114,ln95,ln115,ln96,ln116,ln97,ln117 
				__m128i l13l15low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+13*32),sse_128_load_vector_a(psrc+15*32));			//ln130,ln150,ln131,ln151,ln132,ln152,ln133,ln153 
				__m128i l13l15high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+13*32),sse_128_load_vector_a(psrc+15*32));			//ln134,ln154,ln135,ln155,ln136,ln156,ln137,ln157 

				__m128i l1357low_low = sse128_unpacklo_u32(l1l3low,l5l7low);													//ln10,ln30,ln50,ln70,ln11,ln31,ln51,ln71
				__m128i l1357low_high = sse128_unpackhi_u32(l1l3low,l5l7low);													//ln12,ln32,ln52,ln72,ln13,ln33,ln53,ln73
				__m128i l1357high_low = sse128_unpacklo_u32(l1l3high,l5l7high);													//ln14,ln34,ln54,ln74,ln15,ln35,ln55,ln75
				__m128i l1357high_high = sse128_unpackhi_u32(l1l3high,l5l7high);												//ln16,ln36,ln56,ln76,ln17,ln37,ln57,ln77
				__m128i l9111315low_low = sse128_unpacklo_u32(l9l11low,l13l15low);												//ln90,ln110,ln130,ln150,ln91,ln111,ln131,ln151
				__m128i l9111315low_high = sse128_unpackhi_u32(l9l11low,l13l15low);												//ln92,ln112,ln132,ln152,ln93,ln113,ln133,ln153
				__m128i l9111315high_low = sse128_unpacklo_u32(l9l11high,l13l15high);											//ln94,ln114,ln134,ln154,ln95,ln115,ln135,ln155
				__m128i l9111315high_high = sse128_unpackhi_u32(l9l11high,l13l15high);											//ln96,ln116,ln136,ln156,ln97,ln117,ln137,ln157

				//16 samples per column:(1,3,5,7,...,31). eigth per iteration.  first iteration = (1,3,5,7,9,11,13,15), second iteration = (17,19,21,23,25,27,29,31)  
				odd_c[0][j] = sse128_unpacklo_u64(l1357low_low,l9111315low_low);												//ln10,ln30,ln50,ln70,ln90,ln110,ln130,ln150
				odd_c[1][j] = sse128_unpackhi_u64(l1357low_low,l9111315low_low);												//ln11,ln31,ln51,ln71,ln91,ln111,ln131,ln151
				odd_c[2][j] = sse128_unpacklo_u64(l1357low_high,l9111315low_high);												//ln12,ln32,ln52,ln72,ln92,ln112,ln132,ln152
				odd_c[3][j] = sse128_unpackhi_u64(l1357low_high,l9111315low_high);												//ln13,ln33,ln53,ln73,ln93,ln113,ln133,ln153
				odd_c[4][j] = sse128_unpacklo_u64(l1357high_low,l9111315high_low);												//ln14,ln34,ln54,ln74,ln94,ln114,ln134,ln154
				odd_c[5][j] = sse128_unpackhi_u64(l1357high_low,l9111315high_low);												//ln15,ln35,ln55,ln75,ln95,ln115,ln135,ln155
				odd_c[6][j] = sse128_unpacklo_u64(l1357high_high,l9111315high_high);											//ln14,ln34,ln54,ln74,ln94,ln114,ln134,ln154
				odd_c[7][j] = sse128_unpackhi_u64(l1357high_high,l9111315high_high);											//ln15,ln35,ln55,ln75,ln95,ln115,ln135,ln155
				
				psrc+=16*32;
			}

			//reinit psrc
			psrc = asrc[k]+i*8;
			
			l2l6low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+2*32),sse_128_load_vector_a(psrc+6*32));				//ln20,ln60,ln21,ln61,ln22,ln62,ln23,ln63 
			l2l6high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+2*32),sse_128_load_vector_a(psrc+6*32));				//ln24,ln64,ln25,ln65,ln26,ln66,ln27,ln67 
			l10l14low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+10*32),sse_128_load_vector_a(psrc+14*32));			//ln100,ln140,ln101,ln141,ln102,ln142,ln103,ln143 
			l10l14high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+10*32),sse_128_load_vector_a(psrc+14*32));			//ln104,ln144,ln105,ln145,ln106,ln146,ln107,ln147 
			l18l22low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+18*32),sse_128_load_vector_a(psrc+22*32));			//ln180,ln220,ln181,ln221,ln182,ln222,ln183,ln223 
			l18l22high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+18*32),sse_128_load_vector_a(psrc+22*32));			//ln184,ln224,ln185,ln225,ln186,ln226,ln187,ln227 
			l26l30low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+26*32),sse_128_load_vector_a(psrc+30*32));			//ln260,ln300,ln261,ln301,ln262,ln302,ln263,ln303 
			l26l30high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+26*32),sse_128_load_vector_a(psrc+30*32));			//ln264,ln304,ln265,ln305,ln266,ln306,ln267,ln307 

			l261014low_low = sse128_unpacklo_u32(l2l6low,l10l14low);													//ln20,ln60,ln100,ln140,ln21,ln61,ln101,ln141
			l261014low_high = sse128_unpackhi_u32(l2l6low,l10l14low);
			l261014high_low = sse128_unpacklo_u32(l2l6high,l10l14high);
			l261014high_high = sse128_unpackhi_u32(l2l6high,l10l14high);
			l18222630low_low = sse128_unpacklo_u32(l18l22low,l26l30low);
			l18222630low_high = sse128_unpackhi_u32(l18l22low,l26l30low);
			l18222630high_low = sse128_unpacklo_u32(l18l22high,l26l30high);
			l18222630high_high = sse128_unpackhi_u32(l18l22high,l26l30high);

			//8 samples per column:(2,6,10,14,...,30). 
			eo_c[0] = sse128_unpacklo_u64(l261014low_low,l18222630low_low);												//ln20,ln60,ln100,ln140,ln180,ln220,ln260,ln300
			eo_c[1] = sse128_unpackhi_u64(l261014low_low,l18222630low_low);
			eo_c[2] = sse128_unpacklo_u64(l261014low_high,l18222630low_high);
			eo_c[3] = sse128_unpackhi_u64(l261014low_high,l18222630low_high);
			eo_c[4] = sse128_unpacklo_u64(l261014high_low,l18222630high_low);
			eo_c[5] = sse128_unpackhi_u64(l261014high_low,l18222630high_low);	
			eo_c[6] = sse128_unpacklo_u64(l261014high_high,l18222630high_high);	
			eo_c[7] = sse128_unpackhi_u64(l261014high_high,l18222630high_high);			

			l4l12low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+4*32),sse_128_load_vector_a(psrc+12*32));				//ln40,ln120,ln41,ln121,ln42,ln122,ln43,ln123 
			l4l12high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+4*32),sse_128_load_vector_a(psrc+12*32));			//ln44,ln124,ln45,ln125,ln46,ln126,ln47,ln127 
			l20l28low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+20*32),sse_128_load_vector_a(psrc+28*32));			//ln200,ln280,ln201,ln281,ln202,ln282,ln203,ln283 
			l20l28high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+20*32),sse_128_load_vector_a(psrc+28*32));			//ln204,ln284,ln205,ln285,ln206,ln286,ln207,ln287 

			_aux = sse128_unpacklo_u32(l4l12low,l20l28low);						//ln40,ln120,ln200,ln280,ln41,ln121,ln201,ln281
			eeo_c[0] = sse128_unpacklo_u64(_aux,_aux);							//ln40,ln120,ln200,ln280,ln40,ln120,ln200,ln280
			eeo_c[1] = sse128_unpackhi_u64(_aux,_aux);							//ln41,ln121,ln201,ln281,ln41,ln121,ln201,ln281
			_aux = sse128_unpackhi_u32(l4l12low,l20l28low);						//ln42,ln122,ln202,ln282,ln43,ln123,ln203,ln283																		
			eeo_c[2] = sse128_unpacklo_u64(_aux,_aux);							//ln42,ln122,ln202,ln282,ln42,ln122,ln202,ln282
			eeo_c[3] = sse128_unpackhi_u64(_aux,_aux);							//ln43,ln123,ln203,ln283,ln43,ln123,ln203,ln283
			_aux = sse128_unpacklo_u32(l4l12high,l20l28high);
			eeo_c[4] = sse128_unpacklo_u64(_aux,_aux);
			eeo_c[5] = sse128_unpackhi_u64(_aux,_aux);
			_aux = sse128_unpackhi_u32(l4l12high,l20l28high);
			eeo_c[6] = sse128_unpacklo_u64(_aux,_aux);
			eeo_c[7] = sse128_unpackhi_u64(_aux,_aux);

			l0l8low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc),sse_128_load_vector_a(psrc+8*32));					//ln00,ln80,ln01,ln81,ln02,ln82,ln03,ln83 
			l0l8high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc),sse_128_load_vector_a(psrc+8*32));					//ln04,ln84,ln05,ln85,ln06,ln86,ln07,ln87 
			l16l24low = sse128_unpacklo_u16(sse_128_load_vector_a(psrc+16*32),sse_128_load_vector_a(psrc+24*32));			//ln160,ln240,ln161,ln241,ln162,ln242,ln163,ln243 
			l16l24high = sse128_unpackhi_u16(sse_128_load_vector_a(psrc+16*32),sse_128_load_vector_a(psrc+24*32));			//ln164,ln244,ln165,ln245,ln166,ln246,ln167,ln247 

			_aux = sse128_unpacklo_u32(l0l8low,l16l24low);						//ln00,ln80,ln160,ln240,ln01,ln81,ln161,ln241
			eee_c[0] = sse128_unpacklo_u64(_aux,_aux);							//ln00,ln80,ln160,ln240,ln00,ln80,ln160,ln240
			eee_c[1] = sse128_unpackhi_u64(_aux,_aux);							//ln01,ln81,ln161,ln241,ln01,ln81,ln161,ln241
			_aux = sse128_unpackhi_u32(l0l8low,l16l24low);						//ln02,ln82,ln162,ln242,ln03,ln83,ln163,ln243																		
			eee_c[2] = sse128_unpacklo_u64(_aux,_aux);							//ln02,ln82,ln162,ln242,ln02,ln82,ln162,ln242
			eee_c[3] = sse128_unpackhi_u64(_aux,_aux);							//ln03,ln83,ln163,ln243,ln03,ln83,ln163,ln243
			_aux = sse128_unpacklo_u32(l0l8high,l16l24high);				
			eee_c[4] = sse128_unpacklo_u64(_aux,_aux);
			eee_c[5] = sse128_unpackhi_u64(_aux,_aux);
			_aux = sse128_unpackhi_u32(l0l8high,l16l24high);
			eee_c[6] = sse128_unpacklo_u64(_aux,_aux);
			eee_c[7] = sse128_unpackhi_u64(_aux,_aux);

			round = sse_128_vector_i32(iround[k]);
			shift = ishift[k];

			for (j = 0; j < 8; j++)
			{
				//una linea
				__m128i aux0 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[0])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[0]+8)));
				__m128i aux1 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[1])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[1]+8)));
				__m128i aux2 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[2])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[2]+8)));
				__m128i aux3 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[3])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[3]+8)));
				O[0] = sse_128_hadd_i32(sse_128_hadd_i32(aux0, aux1), sse_128_hadd_i32(aux2, aux3));
				aux0 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[4])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[4]+8)));
				aux1 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[5])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[5]+8)));
				aux2 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[6])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[6]+8)));
				aux3 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[7])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[7]+8)));
				O[1] = sse_128_hadd_i32(sse_128_hadd_i32(aux0, aux1), sse_128_hadd_i32(aux2, aux3));
				aux0 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[8])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[8]+8)));
				aux1 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[9])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[9]+8)));
				aux2 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[10])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[10]+8)));
				aux3 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[11])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[11]+8)));
				O[2] = sse_128_hadd_i32(sse_128_hadd_i32(aux0, aux1), sse_128_hadd_i32(aux2, aux3));							
				aux0 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[12])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[12]+8)));
				aux1 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[13])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[13]+8)));
				aux2 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[14])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[14]+8)));
				aux3 = sse_128_hadd_i32(sse_128_madd_i16_i32(odd_c[j][0],sse_128_load_vector_a(inv_dct32x32_odd[15])), sse_128_madd_i16_i32(odd_c[j][1],sse_128_load_vector_a(inv_dct32x32_odd[15]+8)));
				O[3] = sse_128_hadd_i32(sse_128_hadd_i32(aux0, aux1), sse_128_hadd_i32(aux2, aux3));							

				//una linea
				aux0 = sse_128_hadd_i32(sse_128_madd_i16_i32(eo_c[j],sse_128_load_vector_a(inv_dct32x32_eo[0])), sse_128_madd_i16_i32(eo_c[j],sse_128_load_vector_a(inv_dct32x32_eo[1])));
				aux1 = sse_128_hadd_i32(sse_128_madd_i16_i32(eo_c[j],sse_128_load_vector_a(inv_dct32x32_eo[2])), sse_128_madd_i16_i32(eo_c[j],sse_128_load_vector_a(inv_dct32x32_eo[3])));
				aux2 = sse_128_hadd_i32(sse_128_madd_i16_i32(eo_c[j],sse_128_load_vector_a(inv_dct32x32_eo[4])), sse_128_madd_i16_i32(eo_c[j],sse_128_load_vector_a(inv_dct32x32_eo[5])));
				aux3 = sse_128_hadd_i32(sse_128_madd_i16_i32(eo_c[j],sse_128_load_vector_a(inv_dct32x32_eo[6])), sse_128_madd_i16_i32(eo_c[j],sse_128_load_vector_a(inv_dct32x32_eo[7])));				
				EO[0] = sse_128_hadd_i32(aux0, aux1);
				EO[1] = sse_128_hadd_i32(aux2, aux3);

				//una linea
				EEO = sse_128_hadd_i32(sse_128_madd_i16_i32(eeo_c[j],sse_128_load_vector_a(inv_dct32x32_eeo[0])), sse_128_madd_i16_i32(eeo_c[j],sse_128_load_vector_a(inv_dct32x32_eeo[1])));
				EEE = sse_128_hadd_i32(sse_128_madd_i16_i32(eee_c[j],sse_128_load_vector_a(inv_dct32x32_eee[0])), sse_128_madd_i16_i32(eee_c[j],sse_128_load_vector_a(inv_dct32x32_eee[1])));
				
				EE[0] = sse_128_add_i32(EEE,EEO);
				EE[1] = sse_128_shuffle_32(sse_128_sub_i32(EEE,EEO), 0x1b);

				E[0] = sse_128_add_i32(EE[0],EO[0]);
				E[1] = sse_128_add_i32(EE[1],EO[1]);
				E[2] = sse_128_shuffle_32(sse_128_sub_i32(EE[1],EO[1]), 0x1b);
				E[3] = sse_128_shuffle_32(sse_128_sub_i32(EE[0],EO[0]), 0x1b);

				aux0 = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E[0], O[0]), round), shift), sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E[1], O[1]), round), shift));	
				aux1 = sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E[2], O[2]), round), shift), sse_128_shift_r_i32(sse_128_add_i32(sse_128_add_i32(E[3], O[3]), round), shift));	
				aux2 = sse_128_shuffle_8(sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E[2], O[2]), round), shift), sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E[3], O[3]), round), shift)), shuff_mask_0);
				aux3 = sse_128_shuffle_8(sse128_packs_u32_u16(sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E[0], O[0]), round), shift), sse_128_shift_r_i32(sse_128_add_i32(sse_128_sub_i32(E[1], O[1]), round), shift)), shuff_mask_0);	

				sse_128_store_vector_a(pdst, aux0);
				sse_128_store_vector_a(pdst+8, aux1);
				sse_128_store_vector_a(pdst+16, aux2);
				sse_128_store_vector_a(pdst+24, aux3);
				pdst+=stride;
			}
		}
	}
}


//-----------------------interface functions-----------------------------
#ifndef REG_DCT
	#define REG_DCT					65535
#endif


void sse_transform(int bitDepth, int16_t *block,int16_t *coeff, int block_size, int iWidth, int iHeight, int width_shift, int height_shift, unsigned short uiMode, int16_t *aux)
{
	if(iWidth == 4 && iHeight == 4)
	{
		if (uiMode != REG_DCT)
		{
			sse_aligned_dst_4x4(block, coeff, block_size);
		}
		else
		{
			sse_aligned_dct_4x4(block, coeff, block_size);
		}
	}
	else if(iWidth == 8 && iHeight == 8)
	{
		sse_aligned_dct_8x8(block, coeff, block_size);
	}
	else if(iWidth == 16 && iHeight == 16)
	{
		sse_aligned_dct_16x16(block, coeff, block_size);
	}
	else if (iWidth == 32 && iHeight == 32)
	{
		sse_aligned_dct_32x32(block, coeff, block_size, aux);
	}	

}
#define SHIFT_INV_1ST          7 // Shift after first inverse transform stage
#define SHIFT_INV_2ND         12 // Shift after second inverse transform stage

void sse_itransform(int bitDepth, short *block,short *coeff, int block_size, int iWidth, int iHeight, unsigned int uiMode, short *aux)
{
	int shift_1st = SHIFT_INV_1ST;//g_aucConvertToBit[iWidth]  + 1 + bitDepth-8; // log2(iWidth) - 1 + g_bitDepth - 8
	int shift_2nd = SHIFT_INV_2ND - (bitDepth-8);//g_aucConvertToBit[iHeight]  + 8;                   // log2(iHeight) + 6

	if(iWidth == 4 && iHeight == 4)
	{
		if (uiMode != REG_DCT)
		{
			sse_aligned_inv_dst_4x4(coeff, block, block_size);
		}
		else
		{
			sse_aligned_inv_dct_4x4(coeff, block, block_size); 
		}
	}
	else if(iWidth == 8 && iHeight == 8)
	{
		sse_aligned_inv_dct_8x8(coeff, block, block_size);
	}
	else if(iWidth == 16 && iHeight == 16)
	{
		sse_aligned_inv_dct_16x16(coeff, block, block_size, aux);
	}
	else if (iWidth == 32 && iHeight == 32)
	{
		sse_aligned_inv_dct_32x32(coeff, block, block_size, aux);
	}
}
