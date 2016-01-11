/*****************************************************************************
 * hmr_sse42_functions_pixel.c : homerHEVC encoding library
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

#include "hmr_os_primitives.h"
#include "hmr_sse42_primitives.h"
#include "hmr_sse42_macros.h"
#include "hmr_sse42_functions.h"

//#define EXTRA_OPTIMIZATION	1

void copy_16_16(void* vsrc, uint32_t src_stride, void* vdst, uint32_t dst_stride, int height, int width)
{
	int16_t *src = (int16_t*)vsrc;
	int16_t *dst = (int16_t*)vdst;
	int i,j;
	for(j=0;j<height;j++)
	{
		//for(i=0;i<size;i++)
			//dst[i] = src[i];
		memcpy(dst,src,width*sizeof(src[0]));

		src += src_stride;
		dst += dst_stride;
	}
}


void sse_copy_16_16_4xn(int16_t* src, uint32_t src_stride, int16_t* dst, uint32_t dst_stride, int height)
{
	int i,j;
	for(j=0;j<height;j++)
	{
		__m128_i16	_128_i16_src = sse_128_loadlo_vector64(src);
		sse_64_storel_vector_u(dst, _128_i16_src);

		src += src_stride;
		dst += dst_stride;
	}
}

void sse_copy_16_16_8xn(int16_t* src, uint32_t src_stride, int16_t * dst, uint32_t dst_stride, int height)
{
	int i,j;
	for(j=0;j<height;j++)
	{
		__m128_i16	_128_i16_src = sse_128_load_vector_u(src);
		sse_128_store_vector_u(dst, _128_i16_src);

		src += src_stride;
		dst += dst_stride;
	}
}

void sse_copy_16_16_16xn(int16_t* src, uint32_t src_stride, int16_t * dst, uint32_t dst_stride, int height)
{
	int i,j;
	for(j=0;j<height;j++)
	{
		sse_128_store_vector_u(dst, sse_128_load_vector_u(src));
		sse_128_store_vector_u(dst+8, sse_128_load_vector_u(src+8));

		src += src_stride;
		dst += dst_stride;
	}
}


void sse_copy_16_16_16nxn(int16_t* src, uint32_t src_stride, int16_t * dst, uint32_t dst_stride, int height, int width)
{
	int i,j;
	for(j=0;j<height;j++)
	{
		for(i=0;i<width;i+=16)
		{
			__m128_i16	_128_i16_src = sse_128_load_vector_u(&src[i]);
			__m128_i16	_128_i16_src1 = sse_128_load_vector_u(&src[i+8]);
			sse_128_store_vector_u(&dst[i], _128_i16_src);
			sse_128_store_vector_u(&dst[i+8], _128_i16_src1);
		}

		src += src_stride;
		dst += dst_stride;
	}
}

void sse_copy_16_16(void* vsrc, uint32_t src_stride, void* vdst, uint32_t dst_stride, int height, int width)
{
	int i,j;
	int16_t* src  = (int16_t*)vsrc;
	int16_t* dst  = (int16_t*)vdst;

	if(width==4)
		return sse_copy_16_16_4xn(src, src_stride, dst, dst_stride, height);
	else if(width==8)
		return sse_copy_16_16_8xn(src, src_stride, dst, dst_stride, height);
//	else if(size==16)
//		return sse_copy_16_16_16xn(src, src_stride, dst, dst_stride, size);
	else
		return sse_copy_16_16_16nxn(src, src_stride, dst, dst_stride, height, width);
}


//void copy_8_16(uint8_t* src, uint32_t src_stride, int16_t * dst, uint32_t dst_stride, int size)
void copy_8_16(void* vsrc, uint32_t src_stride, void* vdst, uint32_t dst_stride, int height, int width)
{
	uint8_t *src = (uint8_t*)vsrc;
	int16_t *dst = (int16_t*)vdst;
	int i,j;

	for(j=0;j<height;j++)
	{
		for(i=0;i<width;i++)
			dst[i] = (int16_t)src[i];

		src += src_stride;
		dst += dst_stride;
	}
}


void sse_copy_8_16_4xn(uint8_t* src, uint32_t src_stride, int16_t * dst, uint32_t dst_stride, int height)
{
	int i,j;
	for(j=0;j<height;j++)
	{
		__m128_u8	_128_u8_src0 = sse_128_loadlo_vector64(src);
		__m128_i16	_128_i16_src = sse_128_convert_u8_i16(_128_u8_src0);//
		sse_64_storel_vector_u(dst, _128_i16_src);

		src += src_stride;
		dst += dst_stride;
	}
}

void sse_copy_8_16_8xn(uint8_t* src, uint32_t src_stride, int16_t * dst, uint32_t dst_stride, int height)
{
	int i,j;
	for(j=0;j<height;j++)
	{
		__m128_u8	_128_u8_src0 = sse_128_loadlo_vector64(src);
		__m128_i16	_128_i16_src = sse_128_convert_u8_i16(_128_u8_src0);//
		sse_128_store_vector_u(dst, _128_i16_src);

		src += src_stride;
		dst += dst_stride;
	}
}

void sse_copy_8_16_16nx16n(uint8_t* src, uint32_t src_stride, int16_t * dst, uint32_t dst_stride, int height, int width)
{
	int i,j;
	for(j=0;j<height;j++)
	{
		for(i=0;i<width;i+=16)
		{
			__m128_u8	_128_u8_src = sse_128_load_vector_u(&src[i]);
			__m128_i16	_128_i16_src0 = sse_128_convert_u8_i16(_128_u8_src);//
			__m128_i16	_128_i16_src1 = sse_128_convert_u8_i16(sse128_unpackhi_u64(_128_u8_src, _128_u8_src));//
			sse_128_store_vector_u(&dst[i], _128_i16_src0);
			sse_128_store_vector_u(&dst[i+8], _128_i16_src1);
		}

		src += src_stride;
		dst += dst_stride;
	}
}



void sse_copy_8_16(void* vsrc, uint32_t src_stride, void* vdst, uint32_t dst_stride, int height, int width)
{
	uint8_t *src = (uint8_t*)vsrc;
	int16_t *dst = (int16_t*)vdst;
	if(width==4)
		return sse_copy_8_16_4xn(src, src_stride, dst, dst_stride, height);
	else if(width==8)
		return sse_copy_8_16_8xn(src, src_stride, dst, dst_stride, height);
	else //if(size==16)
		return sse_copy_8_16_16nx16n(src, src_stride, dst, dst_stride, height, width);
}



void copy_16_8(void* vsrc, uint32_t src_stride, void* vdst, uint32_t dst_stride, int height, int width)
{
	int16_t *src = (int16_t*)vsrc;
	uint8_t *dst = (uint8_t*)vdst;

	int i,j;
	for(j=0;j<width;j++)
	{
		for(i=0;i<height;i++)
			dst[i] = (uint8_t)src[i];

		src += src_stride;
		dst += dst_stride;
	}
}

void sse_copy_16_8_4xn(int16_t* src, uint32_t src_stride, uint8_t* dst, uint32_t dst_stride, int height)
{
	int i,j;
	for(j=0;j<height;j++)
	{
//		for(i=0;i<size;i+=4)
		{
			__m128_i16	_128_i16_src0 = sse_128_loadlo_vector64(src);
			__m128_u8	_128_u8_src = sse128_packs_i16_u8(_128_i16_src0, _128_i16_src0);//
			sse_32_store_vector0_u(dst, _128_u8_src);
		}

		src += src_stride;
		dst += dst_stride;
	}
}

void sse_copy_16_8_8xn(int16_t* src, uint32_t src_stride, uint8_t* dst, uint32_t dst_stride, int height)
{
	int i,j;
	for(j=0;j<height;j++)
	{
//		for(i=0;i<size;i+=8)
		{
			__m128_i16	_128_i16_src0 = sse_128_load_vector_u(src);
			__m128_u8	_128_u8_src = sse128_packs_i16_u8(_128_i16_src0, _128_i16_src0);//
			sse_64_storel_vector_u(dst, _128_u8_src);
		}

		src += src_stride;
		dst += dst_stride;
	}
}

void sse_copy_16_8_16nx16n(int16_t* src, uint32_t src_stride, uint8_t* dst, uint32_t dst_stride, int height, int width)
{
	int i,j;
	for(j=0;j<height;j++)
	{
		for(i=0;i<width;i+=16)
		{
			__m128_i16	_128_i16_src0 = sse_128_load_vector_u(&src[i]);
			__m128_i16	_128_i16_src1 = sse_128_load_vector_u(&src[i+8]);
			__m128_u8	_128_u8_src = sse128_packs_i16_u8(_128_i16_src0, _128_i16_src1);//
			sse_128_store_vector_u(&dst[i], _128_u8_src);
		}

		src += src_stride;
		dst += dst_stride;
	}
}


void sse_copy_16_8(void* vsrc, uint32_t src_stride, void* vdst, uint32_t dst_stride, int height, int width)
{
	int16_t *src = (int16_t*)vsrc;
	uint8_t *dst = (uint8_t*)vdst;
	if(width==4)
		return sse_copy_16_8_4xn(src, src_stride, dst, dst_stride, height);
	else if(width==8)
		return sse_copy_16_8_8xn(src, src_stride, dst, dst_stride, height);
	else //if(size==16)
		return sse_copy_16_8_16nx16n(src, src_stride, dst, dst_stride, height, width);
}








uint32_t sse_aligned_sad_4x4(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	uint32_t sad = 0;
//	__m128_u32	_128u32_zero = sse_128_zero_vector();

	//lines 0 and 1
	__m128_i16	_128_i16_src0 = sse128_unpacklo_u64(sse_128_loadlo_vector64(src), sse_128_loadlo_vector64(src+src_stride));//
	__m128_i16	_128_i16_pred0 = sse128_unpacklo_u64(sse_128_loadlo_vector64(pred),sse_128_loadlo_vector64(pred+pred_stride)); 	
	__m128_i16	_128_u16_abs0 = sse_128_sad_i16(_128_i16_src0, _128_i16_pred0);//sse_128_abs_i16(sse_128_sub_i16(_128_i16_src0, _128_i16_pred0));

	//lines 2 and 3
	__m128_i16	_128_i16_src1 = sse128_unpacklo_u64(sse_128_loadlo_vector64(src+2*src_stride), sse_128_loadlo_vector64(src+3*src_stride));//
	__m128_i16	_128_i16_pred1 = sse128_unpacklo_u64(sse_128_loadlo_vector64(pred+2*pred_stride),sse_128_loadlo_vector64(pred+3*pred_stride));
	__m128_i16	_128_u16_abs1 = sse_128_sad_i16(_128_i16_src1, _128_i16_pred1);

	__m128_i16	_128_u16_acc0 = sse_128_add_i16(_128_u16_abs0, _128_u16_abs1);

	return sad = sse_128_hacc_i16(_128_u16_acc0);
}

uint32_t sse_aligned_sad_8x8(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	int i;
	uint32_t sad = 0;
	__m128_u32	_128u32_result = sse_128_zero_vector();

	for(i=0;i<8;i++)
	{
		__m128_i16	_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride, pred+i*pred_stride);
		_128u32_result = sse_128_add_i16(_128u32_result, _128_u16_sad0);
	}

	return sad = sse_128_hacc_i16(_128u32_result);
}


uint32_t sse_aligned_sad_16x16(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	int i;
	uint32_t sad = 0;
	__m128_u32	_128u32_result = sse_128_zero_vector();

	for(i=0;i<16;i++)
	{
		__m128_i16	_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride, pred+i*pred_stride);
		_128u32_result = sse_128_add_i16(_128u32_result, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+8, pred+i*pred_stride+8);
		_128u32_result = sse_128_add_i16(_128u32_result, _128_u16_sad0);
	}
	return sad = sse_128_hacc_i16(_128u32_result);
}


uint32_t sse_aligned_sad_32x32(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	int i,j;
	uint32_t sad = 0;
	__m128_u32	_128u32_acc = sse_128_zero_vector();
	__m128_u16	_128u32_result0 = sse_128_zero_vector();
	__m128_u16	_128u32_result1 = sse_128_zero_vector();
	for(i=0;i<32;i++)
	{
		__m128_u16	_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride, pred+i*pred_stride);
		_128u32_result0 = sse_128_adds_u16(_128u32_result0, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+8, pred+i*pred_stride+8);
		_128u32_result0 = sse_128_adds_u16(_128u32_result0, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+16, pred+i*pred_stride+16);
		_128u32_result1 = sse_128_adds_u16(_128u32_result1, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+24, pred+i*pred_stride+24);
		_128u32_result1 = sse_128_adds_u16(_128u32_result1, _128_u16_sad0);
	}

	_128u32_acc = sse_128_convert_u16_i32(_128u32_result0);
	_128u32_acc = sse_128_add_i32(_128u32_acc, sse_128_convert_u16_i32(sse128_unpackhi_u64(_128u32_result0,_128u32_result0)));
	_128u32_acc = sse_128_add_i32(_128u32_acc, sse_128_convert_u16_i32(_128u32_result1));
	_128u32_acc = sse_128_add_i32(_128u32_acc, sse_128_convert_u16_i32(sse128_unpackhi_u64(_128u32_result1,_128u32_result1)));

	return sad = sse_128_hacc_i32(_128u32_acc);
}


uint32_t sse_aligned_sad_64x64(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	int i,j;
	uint32_t sad = 0;
	__m128_u32	_128u32_acc = sse_128_zero_vector();

	for(i=0;i<64;i++)
	{
		__m128_i16	_128u32_result = SSE42_128_SAD_i16(src+i*src_stride, pred+i*pred_stride);

		__m128_i16	_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+8, pred+i*pred_stride+8);

		_128u32_result = sse_128_adds_u16(_128u32_result, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+16, pred+i*pred_stride+16);
		_128u32_result = sse_128_adds_u16(_128u32_result, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+24, pred+i*pred_stride+24);
		_128u32_result = sse_128_adds_u16(_128u32_result, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+32, pred+i*pred_stride+32);
		_128u32_result = sse_128_adds_u16(_128u32_result, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+40, pred+i*pred_stride+40);
		_128u32_result = sse_128_adds_u16(_128u32_result, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+48, pred+i*pred_stride+48);
		_128u32_result = sse_128_adds_u16(_128u32_result, _128_u16_sad0);

		_128_u16_sad0 = SSE42_128_SAD_i16(src+i*src_stride+56, pred+i*pred_stride+56);
		_128u32_result = sse_128_adds_u16(_128u32_result, _128_u16_sad0);

		_128u32_acc = sse_128_add_i32(_128u32_acc, sse_128_convert_u16_i32(_128u32_result));
		_128u32_acc = sse_128_add_i32(_128u32_acc, sse_128_convert_u16_i32(sse128_unpackhi_u64(_128u32_result,_128u32_result)));
	}
	return sad = sse_128_hacc_i32(_128u32_acc);
}


uint32_t sse_aligned_sad(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride, int size)
{
	if(size==4)
		return sse_aligned_sad_4x4(src, src_stride, pred, pred_stride);
	else if(size==8)
		return sse_aligned_sad_8x8(src, src_stride, pred, pred_stride);
	else if(size==16)
		return sse_aligned_sad_16x16(src, src_stride, pred, pred_stride);
	else if(size==32)
		return sse_aligned_sad_32x32(src, src_stride, pred, pred_stride);
	else// if(size==64)
		return sse_aligned_sad_64x64(src, src_stride, pred, pred_stride);
}



//---------------------------------------------ssd ------------------------------------------------------------------
/*uint32_t sse_ssd_nxn_16x16(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride, uint32_t size)
{
	uint32_t i,j,n;
	uint32_t ssd = 0;
	__m128_u8	_128_aux;
	__m128_u16	_128u32_result = sse_128_zero_vector();


	int16_t *psrc = src;
	int16_t *ppred = pred;

	src_stride-=size;
	pred_stride-=size;

	n = size>>4;
	for(j=0;j<size;j++)
	{
		for(i=0;i<n;i++)
		{
			CALC_ALIGNED_SSD_16(_128u32_result, psrc, ppred, _128_aux)

			psrc+=size;
			ppred+=size;
		}
		psrc+=src_stride;
		ppred+=pred_stride;

	}
	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);
}

uint32_t sse_aligned_ssd_4x4(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	uint32_t ssd = 0;
	__m128_u32	_128_aux;
	__m128_u32	_128u32_result = sse_128_zero_vector();


	CALC_ALIGNED_SSD_2x4(_128u32_result, src, src+src_stride, pred, pred+pred_stride, _128_aux)	
	CALC_ALIGNED_SSD_2x4(_128u32_result, src+2*src_stride, src+3*src_stride, pred+2*pred_stride, pred+3*pred_stride, _128_aux)	

	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);
}

uint32_t sse_aligned_ssd_8x8(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	uint32_t ssd = 0;
	__m128_u32	_128_aux;
	__m128_u32	_128u32_result = sse_128_zero_vector();

	CALC_ALIGNED_SSD_2x8(_128u32_result, src, src+src_stride, pred, pred+pred_stride, _128_aux)
	CALC_ALIGNED_SSD_2x8(_128u32_result, src+2*src_stride, src+3*src_stride, pred+2*pred_stride, pred+3*pred_stride, _128_aux)
	CALC_ALIGNED_SSD_2x8(_128u32_result, src+4*src_stride, src+5*src_stride, pred+4*pred_stride, pred+5*pred_stride, _128_aux)
	CALC_ALIGNED_SSD_2x8(_128u32_result, src+6*src_stride, src+7*src_stride, pred+6*pred_stride, pred+7*pred_stride, _128_aux)

	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);
}


uint32_t sse_aligned_ssd_16x16(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	uint32_t ssd = 0;
	__m128_u32	_128_aux;
	__m128_u32	_128u32_result = sse_128_zero_vector();
	int16_t *psrc = src;
	int16_t *ppred = pred;
	int i;

	for(i=0;i<16;i++)
	{
		CALC_ALIGNED_SSD_16(_128u32_result, psrc, ppred, _128_aux);

		psrc+=src_stride;
		ppred+=pred_stride;
	}

	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);
}



uint32_t sse_aligned_ssd_32x32(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	uint32_t ssd = 0;
	__m128_u8	_128_aux;
	__m128_u32	_128u32_result = sse_128_zero_vector();
	int16_t *psrc = src;
	int16_t *ppred = pred;
	int i;

	for(i=0;i<32;i++)
	{
		CALC_ALIGNED_SSD_32(_128u32_result, psrc, ppred, _128_aux);

		psrc+=src_stride;
		ppred+=pred_stride;
	}

	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);
}


uint32_t sse_aligned_ssd_64x64(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	uint32_t ssd = 0;
	__m128_u8	_128_aux;
	__m128_u32	_128u32_result = sse_128_zero_vector();

	int16_t *psrc = src;
	int16_t *ppred = pred;
	int i;

	for(i=0;i<64;i++)
	{
		CALC_ALIGNED_SSD_64(_128u32_result, psrc, ppred, _128_aux);

		psrc+=src_stride;
		ppred+=pred_stride;
	}

	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);

}


uint32_t sse_aligned_ssd(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride, int size)
{
	if(size==4)
		return sse_aligned_ssd_4x4(src, src_stride, pred, pred_stride);
	else if(size==8)
		return sse_aligned_ssd_8x8(src, src_stride, pred, pred_stride);
	else if(size==16)
		return sse_aligned_ssd_16x16(src, src_stride, pred, pred_stride);
	else if(size==32)
		return sse_aligned_ssd_32x32(src, src_stride, pred, pred_stride);
	else// if(size==64)
		return sse_aligned_ssd_64x64(src, src_stride, pred, pred_stride);

}
*/
//---------------------------------------------ssd16b ------------------------------------------------------------------

uint32_t sse_aligned_ssd16b_4x4(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	uint32_t ssd = 0;
//	__m128_u32	_128u32_zero = sse_128_zero_vector();

	//lines 0 and 1
	__m128_i16	_128_i16_src0 = sse128_unpacklo_u64(sse_128_loadlo_vector64(src), sse_128_loadlo_vector64(src+src_stride));//
	__m128_i16	_128_i16_pred0 = sse128_unpacklo_u64(sse_128_loadlo_vector64(pred),sse_128_loadlo_vector64(pred+pred_stride)); 	

	__m128_i16 _128_i16_sub0 = sse_128_sub_i16(_128_i16_src0, _128_i16_pred0);
	__m128_i32 _128_i32_acc0 = sse_128_madd_i16_i32(_128_i16_sub0, _128_i16_sub0);

	//lines 2 and 3
	__m128_i16	_128_i16_src1 = sse128_unpacklo_u64(sse_128_loadlo_vector64(src+2*src_stride), sse_128_loadlo_vector64(src+3*src_stride));//
	__m128_i16	_128_i16_pred1 = sse128_unpacklo_u64(sse_128_loadlo_vector64(pred+2*pred_stride),sse_128_loadlo_vector64(pred+3*pred_stride));
	__m128_i16 _128_i16_sub1 = sse_128_sub_i16(_128_i16_src1, _128_i16_pred1);
	__m128_i32 _128_i32_acc1 = sse_128_madd_i16_i32(_128_i16_sub1, _128_i16_sub1);

	_128_i32_acc0 = sse_128_add_i32(_128_i32_acc0, _128_i32_acc1);
	return ssd = sse_128_hacc_i32(_128_i32_acc0);
}



uint32_t sse_aligned_ssd16b_8x8(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	int i;
	uint32_t ssd = 0;
	__m128_u32	_128u32_result = sse_128_zero_vector();

	for(i=0;i<8;i++)
	{
		__m128_i32	_128_i32_ssd0 = SSE42_128_SSD_i16_i32(src+i*src_stride, pred+i*pred_stride);
		_128u32_result = sse_128_add_i32(_128u32_result, _128_i32_ssd0);
	}

//	return ssd = sse_128_hacc_i32(_128u32_result);
	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);
}


uint32_t sse_aligned_ssd16b_16x16(int16_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride)
{
	int i;
	uint32_t ssd = 0;
	__m128_u32	_128u32_result = sse_128_zero_vector();

	for(i=0;i<16;i++)
	{
		__m128_i32	_128_i32_ssd0 = SSE42_128_SSD_i16_i32(src+i*src_stride, pred+i*pred_stride);
		_128u32_result = sse_128_add_i32(_128u32_result, _128_i32_ssd0);

		_128_i32_ssd0 = SSE42_128_SSD_i16_i32(src+i*src_stride+8, pred+i*pred_stride+8);
		_128u32_result = sse_128_add_i32(_128u32_result, _128_i32_ssd0);
	}

//	return ssd = sse_128_hacc_i32(_128u32_result);
	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);
}




uint32_t sse_aligned_ssd16b_32x32(int16_t *src, uint32_t src_stride, int16_t *pred, uint32_t pred_stride)
{
	uint32_t ssd = 0;
	__m128_u8	_128_aux;
	__m128_u32	_128u32_result = sse_128_zero_vector();
	int16_t *psrc = src;
	int16_t *ppred = pred;
	int i;

	for(i=0;i<32;i++)
	{
		CALC_ALIGNED_SSD16b_32(_128u32_result, psrc, ppred, _128_aux);

		psrc+=src_stride;
		ppred+=pred_stride;
	}

	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);
}


uint32_t sse_aligned_ssd16b_64x64(int16_t *src, uint32_t src_stride, int16_t *pred, uint32_t pred_stride)
{
	uint32_t ssd = 0;
	__m128_u8	_128_aux;
	__m128_u32	_128u32_result = sse_128_zero_vector();

	int16_t *psrc = src;
	int16_t *ppred = pred;
	int i;

	for(i=0;i<64;i++)
	{
		CALC_ALIGNED_SSD16b_64(_128u32_result, psrc, ppred, _128_aux);

		psrc+=src_stride;
		ppred+=pred_stride;
	}

	return ssd = sse_128_get_data_u32(_128u32_result,0)+sse_128_get_data_u32(_128u32_result,1)+sse_128_get_data_u32(_128u32_result,2)+sse_128_get_data_u32(_128u32_result,3);

}


uint32_t sse_aligned_ssd16b(int16_t *src, uint32_t src_stride, int16_t *pred, uint32_t pred_stride, int size)
{
	if(size==4)
		return sse_aligned_ssd16b_4x4(src, src_stride, pred, pred_stride);
	else if(size==8)
		return sse_aligned_ssd16b_8x8(src, src_stride, pred, pred_stride);
	else if(size==16)
		return sse_aligned_ssd16b_16x16(src, src_stride, pred, pred_stride);
	else if(size==32)
		return sse_aligned_ssd16b_32x32(src, src_stride, pred, pred_stride);
	else// if(size==64)
		return sse_aligned_ssd16b_64x64(src, src_stride, pred, pred_stride);
}


//---------------------------------------------predict ------------------------------------------------------------------

void sse_aligned_predict_4x4(int16_t *orig, int orig_stride, int16_t *pred, int pred_stride, int16_t *residual, int residual_stride)
{
	__m128_i16	_128_i16_src0 = sse128_unpacklo_u64(sse_128_loadlo_vector64(orig), sse_128_loadlo_vector64(orig+orig_stride));//
	__m128_i16	_128_i16_pred0 = sse128_unpacklo_u64(sse_128_loadlo_vector64(pred),sse_128_loadlo_vector64(pred+pred_stride)); 	
	__m128_i16	_128_residual0 = sse_128_sub_i16(_128_i16_src0, _128_i16_pred0);//_128_u16_abs0 = sse_128_sad_i16(_128_i16_src0, _128_i16_pred0);//sse_128_abs_i16(sse_128_sub_i16(_128_i16_src0, _128_i16_pred0));
	sse_64_storel_vector_u(residual, _128_residual0);
	sse_64_storeh_vector_u(residual+residual_stride, _128_residual0);

	_128_i16_src0 = sse128_unpacklo_u64(sse_128_loadlo_vector64(orig+2*orig_stride), sse_128_loadlo_vector64(orig+3*orig_stride));//
	_128_i16_pred0 = sse128_unpacklo_u64(sse_128_loadlo_vector64(pred+2*pred_stride),sse_128_loadlo_vector64(pred+3*pred_stride)); 	
	_128_residual0 = sse_128_sub_i16(_128_i16_src0, _128_i16_pred0);//_128_u16_abs0 = sse_128_sad_i16(_128_i16_src0, _128_i16_pred0);//sse_128_abs_i16(sse_128_sub_i16(_128_i16_src0, _128_i16_pred0));
	sse_64_storel_vector_u(residual+2*residual_stride, _128_residual0);
	sse_64_storeh_vector_u(residual+3*residual_stride, _128_residual0);

}



void sse_aligned_predict_8x8(int16_t *orig, int orig_stride, int16_t *pred, int pred_stride, int16_t *residual, int residual_stride)
{
	__m128_		_128_zero = sse_128_zero_vector();
	int j;
	for(j=0;j<8;j++)
	{
		CALC_ALIGNED_PREDICT_8(orig, pred, residual)
		orig+=orig_stride;
		pred+=pred_stride;
		residual+=residual_stride;
	}
}


void sse_aligned_predict_16x16(int16_t *orig, int orig_stride, int16_t *pred, int pred_stride, int16_t *residual, int residual_stride)
{
	int j;
	for(j=0;j<16;j++)
	{
		CALC_ALIGNED_PREDICT_16(orig, pred, residual)
		orig+=orig_stride;
		pred+=pred_stride;
		residual+=residual_stride;
	}
}



void sse_aligned_predict_32x32(int16_t *orig, int orig_stride, int16_t *pred, int pred_stride, int16_t *residual, int residual_stride)
{
	int j;
	for(j=0;j<32;j++)
	{
		CALC_ALIGNED_PREDICT_32(orig, pred, residual)
		orig+=orig_stride;
		pred+=pred_stride;
		residual+=residual_stride;
	}
}

void sse_aligned_predict_64x64(int16_t *orig, int orig_stride, int16_t *pred, int pred_stride, int16_t *residual, int residual_stride)
{
	int j;
	for(j=0;j<64;j++)
	{
		CALC_ALIGNED_PREDICT_64(orig, pred, residual)
		orig+=orig_stride;
		pred+=pred_stride;
		residual+=residual_stride;
	}
}



void sse_aligned_predict(int16_t *orig, int orig_stride, int16_t *pred, int pred_stride, int16_t *residual, int residual_stride, int size)
{
	if(size==4)
		sse_aligned_predict_4x4(orig, orig_stride, pred, pred_stride, residual, residual_stride);
	else if(size==8)
		sse_aligned_predict_8x8(orig, orig_stride, pred, pred_stride, residual, residual_stride);
	else if(size==16)
		sse_aligned_predict_16x16(orig, orig_stride, pred, pred_stride, residual, residual_stride);
	else if(size==32)
		sse_aligned_predict_32x32(orig, orig_stride, pred, pred_stride, residual, residual_stride);
	else if(size==64)
		sse_aligned_predict_64x64(orig, orig_stride, pred, pred_stride, residual, residual_stride);
}


//---------------------------------------------reconst ------------------------------------------------------------------


void sse_aligned_reconst_4x4(int16_t *pred, int pred_stride, int16_t *residual, int residual_stride, int16_t *decoded, int decoded_stride)
{
	int bit_depth = 8;
	__m128_i16	min_limit = sse_128_zero_vector();
	__m128_i16	max_limit = sse_128_vector_i16((1<<bit_depth)-1);

	__m128_i16 _128_aux1 = sse_128_clip_16(sse_128_adds_i16(sse128_unpacklo_u64(sse_128_loadlo_vector64(pred),sse_128_loadlo_vector64(pred+pred_stride)), sse128_unpacklo_u64(sse_128_loadlo_vector64(residual),sse_128_loadlo_vector64(residual+residual_stride))), min_limit, max_limit);
	__m128_i16 _128_aux2 = sse_128_clip_16(sse_128_adds_i16(sse128_unpacklo_u64(sse_128_loadlo_vector64(pred+2*pred_stride),sse_128_loadlo_vector64(pred+3*pred_stride)), sse128_unpacklo_u64(sse_128_loadlo_vector64(residual+2*residual_stride),sse_128_loadlo_vector64(residual+3*residual_stride))), min_limit, max_limit);

	sse_64_storel_vector_u(decoded, _128_aux1);
	sse_64_storeh_vector_u(decoded+decoded_stride, _128_aux1);
	sse_64_storel_vector_u(decoded+2*decoded_stride, _128_aux2);
	sse_64_storeh_vector_u(decoded+3*decoded_stride, _128_aux2);
}


void sse_aligned_reconst_8x8(int16_t *pred, int pred_stride, int16_t *residual, int residual_stride, int16_t *decoded, int decoded_stride)
{
	int bit_depth = 8;
	__m128_i16	min_limit = sse_128_zero_vector();
	__m128_i16	max_limit = sse_128_vector_i16((1<<bit_depth)-1);
	int j;
	for (j=0;j<8;j++)
	{
		CALC_ALIGNED_RECONST_8(pred+j*pred_stride, residual+j*residual_stride, decoded+j*decoded_stride, min_limit, max_limit)	

//		CALC_ALIGNED_RECONST_8(pred, residual, decoded, min_limit, max_limit)	
//		decoded += decoded_stride;
//		residual += residual_stride;//este es 2D.Podria ser lineal
//		pred += pred_stride;
	}
}

void sse_aligned_reconst_16x16(int16_t *pred, int pred_stride, int16_t *residual, int residual_stride, int16_t *decoded, int decoded_stride)
{
	int j;
	int bit_depth = 8;
	__m128_i16	min_limit = sse_128_zero_vector();
	__m128_i16	max_limit = sse_128_vector_i16((1<<bit_depth)-1);

	for (j=0;j<16;j++)
	{
		CALC_ALIGNED_RECONST_16(pred, residual, decoded, min_limit, max_limit)

		decoded += decoded_stride;
		residual += residual_stride;//este es 2D.Podria ser lineal
		pred += pred_stride;
	}
}

void sse_aligned_reconst_32x32(int16_t *pred, int pred_stride, int16_t *residual, int residual_stride, int16_t *decoded, int decoded_stride)
{
	int j;
	int bit_depth = 8;
	__m128_i16	min_limit = sse_128_zero_vector();
	__m128_i16	max_limit = sse_128_vector_i16((1<<bit_depth)-1);

	for (j=0;j<32;j++)
	{
		CALC_ALIGNED_RECONST_32(pred, residual, decoded, min_limit, max_limit)													


		decoded += decoded_stride;
		residual += residual_stride;//este es 2D.Podria ser lineal
		pred += pred_stride;
	}
}

void sse_aligned_reconst_64x64(int16_t *pred, int pred_stride, int16_t *residual, int residual_stride, int16_t *decoded, int decoded_stride)
{
	int j;
	int bit_depth = 8;
	__m128_i16	min_limit = sse_128_zero_vector();
	__m128_i16	max_limit = sse_128_vector_i16((1<<bit_depth)-1);

	for (j=0;j<64;j++)
	{
		CALC_ALIGNED_RECONST_64(pred, residual, decoded, min_limit, max_limit)													

		decoded += decoded_stride;
		residual += residual_stride;//este es 2D.Podria ser lineal
		pred += pred_stride;
	}
}
void sse_aligned_reconst(int16_t *pred, int pred_stride, int16_t *residual, int residual_stride, int16_t *decoded, int decoded_stride, int size)
{
	if(size==4)
		sse_aligned_reconst_4x4(pred, pred_stride, residual, residual_stride, decoded, decoded_stride);
	else if(size==8)
		sse_aligned_reconst_8x8(pred, pred_stride, residual, residual_stride, decoded, decoded_stride);
	else if(size==16)
		sse_aligned_reconst_16x16(pred, pred_stride, residual, residual_stride, decoded, decoded_stride);
	else if(size==32)
		sse_aligned_reconst_32x32(pred, pred_stride, residual, residual_stride, decoded, decoded_stride);
	else if(size==64)
		sse_aligned_reconst_64x64(pred, pred_stride, residual, residual_stride, decoded, decoded_stride);
}


//--------------------------------------- variance -----------------------------------------

uint32_t sse_variance_16nx16n(int16_t *__restrict p, int size, int stride, int modif)
{
	int i,j;
	unsigned int s=0;
	int16_t *__restrict paux = p;
	__m128_i16	_128_one = sse_128_vector_i16(1);
	__m128_i16	_128_modif = sse_128_vector_i16(modif);
	__m128_i32	acc = sse_128_zero_vector();
	__m128_i16	avg_128;
	__m128_i32	var_128;

	//media 
	for (j=0; j<size; j++)
	{
		__m128_i16	s_128 = sse_128_zero_vector();
		for (i=0; i<size; i+=16)
		{
			__m128_i16	v0 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+i));
			__m128_i16	v1 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+i+8));
			s_128 = sse_128_add_i16(s_128, sse_128_add_i16(v0,v1));
		}
		acc = sse_128_add_i32(acc, sse_128_add_i32(sse_128_convert_i16_i32(s_128), sse_128_convert_i16_i32(sse128_unpackhi_u64(s_128,s_128))));

		paux+= stride;
	}
	acc = sse_128_hadd_i32(acc,sse_128_zero_vector());
	acc = sse_128_hadd_i32(acc,sse_128_zero_vector());

	s = sse_128_get_data_u32(acc,0);
	s/=(size*size);
	
	paux = p;

	avg_128 = sse_128_vector_i16(s);
	var_128 = sse_128_zero_vector();
	paux = p;
	for (j=0; j<size; j++)
	{
		__m128_i16	v_128 = sse_128_zero_vector();
		for (i=0; i<size; i+=16)
		{
			__m128_i16	v0 = sse_128_sub_i16(sse_128_convert_u8_i16(sse_128_load_vector_u(paux+i)),avg_128) ;
			__m128_i16	v1 = sse_128_sub_i16(sse_128_convert_u8_i16(sse_128_load_vector_u(paux+i+8)),avg_128);
			v0 = sse_128_add_i16(_128_one,sse_128_mul_i16(v0,_128_modif));
			v1 = sse_128_add_i16(_128_one,sse_128_mul_i16(v1,_128_modif));

			v_128 = sse_128_add_i32(v_128, sse_128_add_i32(sse_128_madd_i16_i32(v0,v0),sse_128_madd_i16_i32(v1,v1)));
		}
		var_128 = sse_128_add_i32(var_128, v_128);
		paux+= stride;
	}
	var_128 = sse_128_hadd_i32(var_128,sse_128_zero_vector());
	var_128 = sse_128_hadd_i32(var_128,sse_128_zero_vector());

	return 	sse_128_get_data_u32(var_128,0);
}

ALIGN(16) static const int8_t shuffle_mask_variance_16_0[16] ={ 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};//0,1,2,3,4,5,6,7 -> 7,6,5,4,3,2,1,0

uint32_t sse_variance_8x8(int16_t *__restrict p, int stride, int modif)
{
	int16_t *__restrict paux = p;
	__m128_i16	_128_modif = sse_128_vector_i16(modif);
	__m128_i16	_128_one = sse_128_vector_i16(1);

	__m128_i16	zero = sse_128_zero_vector();
	__m128_i32	acc = sse_128_zero_vector();
	__m128_i32	var_128;
	__m128_i16	avg_128;

	__m128_i16	v0, v1, v2, v3, v4, v5, v6, v7;
	v0 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux));
	v1 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+stride));
	v2 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+2*stride));
	v3 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+3*stride));
	v4 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+4*stride));
	v5 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+5*stride));
	v6 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+6*stride));
	v7 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+7*stride));

	acc = sse_128_add_i16(sse_128_add_i16(v2, v3), sse_128_add_i16(v0, v1));
	acc = sse_128_add_i16(acc, sse_128_add_i16(sse_128_add_i16(v6, v7), sse_128_add_i16(v4, v5)));
	acc = sse_128_hadd_i16(acc, zero);
	acc = sse_128_hadd_i16(acc, zero);
	acc = sse_128_hadd_i16(acc, zero);

	avg_128 = sse_128_shift_r_i16(sse_128_shuffle_8(acc, sse_128_load_vector_u(shuffle_mask_variance_16_0)),6);

	v0 = sse_128_sub_i16(v0, avg_128);
	v0 = sse_128_add_i16(_128_one,sse_128_mul_i16(v0,_128_modif));
	v1 = sse_128_sub_i16(v1, avg_128);
	v1 = sse_128_add_i16(_128_one,sse_128_mul_i16(v1,_128_modif));
	v2 = sse_128_sub_i16(v2, avg_128);
	v2 = sse_128_add_i16(_128_one,sse_128_mul_i16(v2,_128_modif));
	v3 = sse_128_sub_i16(v3, avg_128);
	v3 = sse_128_add_i16(_128_one,sse_128_mul_i16(v3,_128_modif));
	v4 = sse_128_sub_i16(v4, avg_128);
	v4 = sse_128_add_i16(_128_one,sse_128_mul_i16(v4,_128_modif));
	v5 = sse_128_sub_i16(v5, avg_128);
	v5 = sse_128_add_i16(_128_one,sse_128_mul_i16(v5,_128_modif));
	v6 = sse_128_sub_i16(v6, avg_128);
	v6 = sse_128_add_i16(_128_one,sse_128_mul_i16(v6,_128_modif));
	v7 = sse_128_sub_i16(v7, avg_128);
	v7 = sse_128_add_i16(_128_one,sse_128_mul_i16(v7,_128_modif));


	var_128 = sse_128_madd_i16_i32(v0,v0);
	var_128 = sse_128_add_i32(var_128, sse_128_madd_i16_i32(v1,v1));
	var_128 = sse_128_add_i32(var_128, sse_128_madd_i16_i32(v2,v2));
	var_128 = sse_128_add_i32(var_128, sse_128_madd_i16_i32(v3,v3));
	var_128 = sse_128_add_i32(var_128, sse_128_madd_i16_i32(v4,v4));
	var_128 = sse_128_add_i32(var_128, sse_128_madd_i16_i32(v5,v5));
	var_128 = sse_128_add_i32(var_128, sse_128_madd_i16_i32(v6,v6));
	var_128 = sse_128_add_i32(var_128, sse_128_madd_i16_i32(v7,v7));

	var_128 = sse_128_hadd_i32(var_128,zero);
	var_128 = sse_128_hadd_i32(var_128,zero);

	return 	sse_128_get_data_u32(var_128,0);
}


uint32_t sse_variance_4x4(int16_t *__restrict p, int stride, int modif)
{
	int16_t *__restrict paux = p;
	__m128_i16	_128_modif = sse_128_vector_i16(modif);
	__m128_i16	_128_one = sse_128_vector_i16(1);
	__m128_i32	acc = sse_128_zero_vector();
	__m128_i16	v0 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux));
	__m128_i16	v1 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+stride));
	__m128_i16	v2 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+2*stride));
	__m128_i16	v3 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+3*stride));
	__m128_i16	v4, v5;
	__m128_i16	avg_128, var_128;


	acc = sse_128_add_i16(sse_128_add_i16(v2, v3), sse_128_add_i16(v0, v1));
	acc = sse_128_hadd_i16(acc, acc);
	acc = sse_128_hadd_i16(acc, acc);

	avg_128 = sse_128_shift_r_i16(sse_128_shuffle_8(acc, sse_128_load_vector_u(shuffle_mask_variance_16_0)),4);

	v4 = sse_128_sub_i16(sse128_unpacklo_u16(v0,v1), avg_128);
	v4 = sse_128_add_i16(_128_one,sse_128_mul_i16(v4,_128_modif));
	var_128 = sse_128_madd_i16_i32(v4,v4);

	v5 = sse_128_sub_i16(sse128_unpacklo_u16(v2,v3), avg_128);
	v5 = sse_128_add_i16(_128_one,sse_128_mul_i16(v5,_128_modif));

	var_128 = sse_128_add_i32(var_128, sse_128_madd_i16_i32(v5,v5));

	var_128 = sse_128_hadd_i32(var_128,sse_128_zero_vector());
	var_128 = sse_128_hadd_i32(var_128,sse_128_zero_vector());

	return 	sse_128_get_data_u32(var_128,0);
}

uint32_t sse_variance_2x2(int16_t *__restrict p, int stride, int modif)
{
	int16_t *__restrict paux = p;

	__m128_i16	_128_modif = sse_128_vector_i16(modif);
	__m128_i16	_128_one = sse_128_vector_i16(1);

	__m128_i32	acc = sse_128_zero_vector();

	__m128_i16	v0 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux));
	__m128_i16	v1 = sse_128_convert_u8_i16(sse_128_load_vector_u(paux+stride));
	__m128_i16	avg_128;
	__m128_i16	v3;
	__m128_i16	var_128;

	acc = sse_128_add_i16(v0, v1);
	acc = sse_128_hadd_i16(acc, acc);

	avg_128 = sse_128_shift_r_i16(sse_128_shuffle_8(acc, sse_128_load_vector_u(shuffle_mask_variance_16_0)),2);

	v3 = sse_128_sub_i16(sse128_unpacklo_u16(v0,v1), avg_128);
	v3 = sse_128_add_i16(_128_one,sse_128_mul_i16(v3,_128_modif));
	var_128 = sse_128_madd_i16_i32(v3,v3);

	var_128 = sse_128_hadd_i32(var_128, var_128);
	return 	sse_128_get_data_u32(var_128,0);
}



uint32_t sse_modified_variance(int16_t *p, int size, int stride, int modif)
{
	if(size==2)
		return sse_variance_2x2(p, stride, modif);
	if(size==4)
		return sse_variance_4x4(p, stride, modif);
	else if(size==8)
		return sse_variance_8x8(p, stride, modif);
	else //if(size>8)
		return sse_variance_16nx16n(p, size, stride, modif);
}



