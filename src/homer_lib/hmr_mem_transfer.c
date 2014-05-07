/*****************************************************************************
 * hmr_mem_transfer.c : homerHEVC encoding library
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
#include <memory.h>
#include <malloc.h>

#include "hmr_common.h"
#include "hmr_private.h"

#define ALIGNMENT 64

void *aligned_alloc(int num, int size)
{
    unsigned char *mem = (unsigned char *)malloc(num*size+ALIGNMENT+sizeof(void*));
    void **ptr = (void**)((uint64_t)(mem+ALIGNMENT+sizeof(void*)) & ~(ALIGNMENT-1));
    ptr[-1] = mem;
    return ptr;
}

void aligned_free(void *ptr) {
    free(((void**)ptr)[-1]);
}



//Alignes by increasing offset_x and size_x to match the boundary. 
void wnd_alloc(wnd_t *wnd, int size_x, int size_y, int offset_x, int offset_y, int pix_size)
{
	int i=0;
	int offset_aligment = ((offset_x*pix_size)%16)?16:0;//((16-((offset_x*pix_size)%16))&0xf);//left boundary

	wnd->pix_size = pix_size;

	//aligment - we want the boundaries of ctus alingned to 16 bytes for SSE2-4
	offset_aligment = offset_aligment/wnd->pix_size;//Normalize


	wnd->data_offset_x = offset_aligment;
	wnd->data_offset_y = offset_y;

	for(i=0;i<3;i++)
	{
		int width = (i==Y_COMP)?size_x:(size_x>>1);//420
		int height = (i==Y_COMP)?size_y:(size_y>>1);//420
		int data_aligment = ((width*pix_size)%16)?16:0;//((16-((size_x*pix_size)%16))&0xf);//right boundary - depends on chroma subsampling type
		data_aligment = data_aligment/wnd->pix_size;

		wnd->window_size_x[i] = wnd->data_offset_x+width+data_aligment;
		wnd->window_size_y[i] = height+wnd->data_offset_y;

		if((wnd->pwnd[i] = (void*)aligned_alloc(wnd->window_size_x[i]*wnd->window_size_y[i]*pix_size, sizeof(byte)))==NULL)
			printf("wnd_alloc - unable to allocate memory for wnd->pwnd[%d]\r\n", i);
	}
}

void wnd_delete(wnd_t *wnd)
{
	int i;
	for(i=0;i<3;i++)
	{
		if(wnd->pwnd[i] != NULL)
			aligned_free(wnd->pwnd[i]);
	}
}


void wnd_realloc(wnd_t *wnd, int size_x, int size_y, int offset_x, int offset_y, int pix_size)
{
	wnd_delete(wnd);
	wnd_alloc(wnd, size_x, size_y, offset_x, offset_y, pix_size);
}

#ifdef WRITE_REF_FRAMES
void wnd_write2file(wnd_t *wnd)
{
	byte * __restrict src;
	int component;

	for(component=Y_COMP;component<=V_COMP;component++)
	{
		src = WND_DATA_PTR(byte*, *wnd, component);
		fwrite(src, sizeof(byte), (wnd->window_size_x[component]*wnd->window_size_y[component]), wnd->out_file); 
	}
}
#endif


void mem_transfer_move_curr_ctu_group(henc_thread_t* et, int i, int j)//i,j are cu indexes
{
	int width, height, component;
	wnd_t* dst_wnd = &et->curr_mbs_wnd;
	byte * src;
	byte * dst;
	int src_stride, dst_stride;

	for(component=Y_COMP;component<=V_COMP;component++)
	{
		src = WND_POSITION_2D(byte *, et->ed->current_pict.img2encode->img, component, (i*et->ctu_width[component]), (j*et->ctu_height[component]), 0, et->ctu_width);
		dst = WND_POSITION_2D(byte *, *dst_wnd, component, 0, 0, 0, et->ctu_width);

		src_stride =  WND_STRIDE_2D(et->ed->current_pict.img2encode->img, component);
		dst_stride =  WND_STRIDE_2D(*dst_wnd, component);

		width = ((i+1)*et->ctu_width[component]<et->pict_width[component])?et->ctu_width[component]:(et->pict_width[component]-(i*et->ctu_width[component]));
		height = ((j+1)*et->ctu_height[component]<et->pict_height[component])?et->ctu_height[component]:(et->pict_height[component]-(j*et->ctu_height[component]));

		mem_transfer_2d2d(src, dst, width, height, src_stride, dst_stride);
	}
}


void mem_transfer_decoded_blocks(henc_thread_t* et, ctu_info_t* ctu)
{
//	int l;//, i;
	wnd_t *decoded_src_wnd = &et->decoded_mbs_wnd[0];
	wnd_t *decoded_dst_wnd = &et->ed->curr_ref_wnd->img;
	int component = Y_COMP;
	int src_stride;
	int dst_stride;
	byte * decoded_buff_src;
	byte * decoded_buff_dst;
	int copy_width, copy_height, decoded_frame_width, decoded_frame_height;

	for(component=Y_COMP;component<=V_COMP;component++)
	{
		decoded_buff_src = WND_POSITION_2D(byte *__restrict, *decoded_src_wnd, component, 0, 0, 0, et->ctu_width);
		decoded_buff_dst = WND_POSITION_2D(byte *__restrict, *decoded_dst_wnd, component, ctu->x[component], ctu->y[component], 0, et->ctu_width);
		src_stride =  WND_STRIDE_2D(*decoded_src_wnd, component);
		dst_stride =  WND_STRIDE_2D(*decoded_dst_wnd, component);

		decoded_frame_width = decoded_dst_wnd->window_size_x[component];
		decoded_frame_height = decoded_dst_wnd->window_size_y[component];
		copy_width = ((ctu->x[component]+et->ctu_width[component]*et->ctu_group_size)<decoded_frame_width)?(et->ctu_width[component]*et->ctu_group_size):(decoded_frame_width-ctu->x[component]);
		copy_height = ((ctu->y[component]+et->ctu_height[component])<decoded_frame_height)?(et->ctu_height[component]):(decoded_frame_height-(ctu->y[component]));
		mem_transfer_2d2d(decoded_buff_src, decoded_buff_dst, copy_width, copy_height, src_stride, dst_stride);
	}
}


void mem_transfer_intra_refs(henc_thread_t* et, ctu_info_t* ctu)
{
	int l, j;
	wnd_t *decoded_src_wnd = &et->ed->curr_ref_wnd->img;
	wnd_t *decoded_dst_wnd = &et->decoded_mbs_wnd[0];
	int component = Y_COMP;
	int src_stride = et->pict_width[component];
	int dst_stride;
	byte * __restrict decoded_buff_src;
	byte * __restrict decoded_buff_dst;
	int cu_size = ctu->size;
	int left_copy = 0, top_copy = 0;


	if(!ctu->ctu_left && !ctu->ctu_top)
		return;

	for(l=0;l<NUM_DECODED_WNDS;l++)
	{
		for(component=Y_COMP;component<=V_COMP;component++)
		{
			decoded_dst_wnd = &et->decoded_mbs_wnd[l];
			src_stride = WND_STRIDE_2D(*decoded_src_wnd, component);
			dst_stride = WND_STRIDE_2D(*decoded_dst_wnd, component);
			decoded_buff_src = WND_POSITION_2D(byte *__restrict, *decoded_src_wnd, component, ctu->x[component], ctu->y[component], 0, et->ctu_width);
			decoded_buff_dst = WND_POSITION_2D(byte *__restrict, *decoded_dst_wnd, component, 0, 0, 0, et->ctu_width);
			cu_size = et->ctu_width[component];
			left_copy = top_copy = 0;
			if(ctu->ctu_left)
				left_copy += cu_size;
			if(ctu->ctu_left_bottom)
				left_copy += cu_size;
			if(ctu->ctu_top)
				top_copy += cu_size;
			if(ctu->ctu_top_right)
				top_copy += cu_size;

			if(left_copy>(et->pict_height[component]-ctu->y[component]))
				left_copy=(et->pict_height[component]-ctu->y[component]);
			if(top_copy>(et->pict_width[component]-ctu->x[component]))
				top_copy=(et->pict_width[component]-ctu->x[component]);

			decoded_buff_src-=(src_stride+1);
			decoded_buff_dst-=(dst_stride+1);
	
			//top-left square
			if(ctu->ctu_left && ctu->ctu_top)
			{
				*decoded_buff_dst = *decoded_buff_src;
			}
			decoded_buff_src++;
			decoded_buff_dst++;
			//bottom line
			memcpy(decoded_buff_dst, decoded_buff_src, top_copy*sizeof(decoded_buff_src[0]));

			//right column
			decoded_buff_src+=src_stride-1;
			decoded_buff_dst+=dst_stride-1;
			for(j=0;j<left_copy;j++)
			{
				decoded_buff_dst[j*dst_stride] = decoded_buff_src[j*src_stride];
			}
		}
	}
}


void mem_transfer_1d1d(unsigned char *src, unsigned char *dst, unsigned int width, unsigned int height)
{
	uint l;
	for(l=0;l<height*width;l+=16)
	{
#ifdef SSE_FUNCS_
		sse_128_store_vector_a(dst+l, sse_128_load_vector_u(src+l));
#else
		memcpy(dst,src,height*width);
#endif
	}		
}


void mem_transfer_1d2d(unsigned char *src, unsigned char *dst, unsigned int width, unsigned int height, unsigned int dst_stride)
{
	uint l;
	for(l=0;l<height;l++)
	{
#ifdef SSE_FUNCS_
		for(uint ll=0;ll<width;ll+=16)
		{
			sse_128_store_vector_a(dst+ll, sse_128_load_vector_u(src));
			src +=16;
		}
		dst  += dst_stride;
#else
		memcpy(dst,src,width);
		dst += dst_stride;
		src += width;
#endif
	}		
}



void mem_transfer_2d1d(unsigned char *src, unsigned char *dst, unsigned int width, unsigned int height, unsigned int src_stride)
{
	uint l;
	for(l=0;l<height;l++)
	{
#ifdef SSE_FUNCS_
		for(uint ll=0;ll<width;ll+=16)
		{
			sse_128_store_vector_u(dst+ll, sse_128_load_vector_u(src+ll));
			dst+=16;
		}
		src  += src_stride;
#else
		memcpy(dst,src,width);
		dst+=width;
		src  += src_stride;
#endif
	}
}

void mem_transfer_2d2d(unsigned char *src, unsigned char *dst, unsigned int width, unsigned int height, unsigned int src_stride, unsigned int dst_stride)
{
	uint l;
	for(l=0;l<height;l++)
	{
#ifdef SSE_FUNCS_
		for(uint ll=0;ll<width;ll+=16)
		{
			sse_128_store_vector_a(dst+ll, sse_128_load_vector_u(src+ll));
		}
		dst += dst_stride;
		src  += src_stride;
#else
		memcpy(dst,src,width);
		dst += dst_stride;
		src  += src_stride;
#endif
	}
}

