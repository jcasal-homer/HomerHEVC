/*****************************************************************************
 * hmr_bitstream.c : homerHEVC encoding library
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
#include <malloc.h>
#include "memory.h"
#include "hmr_private.h"
#include "hmr_common.h"

#define ENDIANESS_CHANGE(n)			((((n) & 0x000000FF) << 24)     |   \
									(((n) & 0x0000FF00) << 8)		|   \
									(((n) & 0x00FF0000) >> 8)		|   \
									(((n) & 0xFF000000) >> 24))



void hmr_bitstream_alloc(bitstream_t* bs, int size)
{
	bs->bitstream = (unsigned char*)calloc(size, 1);
	bs->streamsize = size;
	bs->streambitcnt  = 0;
	bs->streambytecnt = 0;
}


void hmr_bitstream_free(bitstream_t* bs)
{
	if(bs->bitstream != NULL)
		free(bs->bitstream);

	bs->bitstream = NULL;
	bs->streamsize = 0;
	bs->streambitcnt  = 0;
	bs->streambytecnt = 0;
}


void hmr_bitstream_init(bitstream_t* bs)
{
  memset(bs->bitstream, 0, bs->streambytecnt+(bs->streambitcnt!=0));
 
  bs->streambitcnt  = 0;
  bs->streambytecnt = 0;
}

void hmr_bitstream_write_bits(bitstream_t* bs, unsigned int val,int n)
{
	unsigned int *stream = (unsigned int *)&bs->bitstream[bs->streambytecnt];
	*stream = *stream | ENDIANESS_CHANGE(((val<<(32-n))>>bs->streambitcnt));

	if((n+bs->streambitcnt) > 32)
	{
		stream = (unsigned int *)&bs->bitstream[bs->streambytecnt+4];
		*stream = ENDIANESS_CHANGE(val<<(n-bs->streambitcnt));	
	}
	
	bs->streambitcnt  += n;
	bs->streambytecnt += bs->streambitcnt>>3;
	bs->streambitcnt  &= 0x7;
}

void hmr_bitstream_write_bits_uvlc(bitstream_t* bs, unsigned int val)
{
	unsigned int length = 1;
	unsigned int temp = ++val;

	while( 1 != temp )
	{
		temp >>= 1;
		length += 2;
	}
	hmr_bitstream_write_bits(bs, 0,length>>1);
	hmr_bitstream_write_bits(bs, val,(length+1)>>1);
}

void hmr_bitstream_write_bits_svlc(bitstream_t* bs, int val)
{
	unsigned int uval = ( val <= 0) ? -val<<1 : (val<<1)-1;

	hmr_bitstream_write_bits_uvlc(bs, uval);
}



void hmr_bitstream_align_bits_1(bitstream_t* bs)
{
  if (bs->streambitcnt != 0)
    hmr_bitstream_write_bits(bs, 0xff,(8-bs->streambitcnt));
}


void hmr_bitstream_align_bits_0(bitstream_t* bs)
{
  if (bs->streambitcnt != 0)
    hmr_bitstream_write_bits(bs, 0,(8-bs->streambitcnt));
}

void hmr_bitstream_rbsp_trailing_bits(bitstream_t* bs)
{
	hmr_bitstream_write_bits(bs, 1,1);
	hmr_bitstream_align_bits_0(bs);
}

byte emulation_prevention_three_byte = 3;


void hmr_bitstream_nalu_ebsp(bitstream_t* in_bs, bitstream_t* out_bs)
{
	unsigned char* ptr = in_bs->bitstream;
	unsigned char* out_ptr = &out_bs->bitstream[out_bs->streambytecnt];
	int size = in_bs->streambytecnt;
	int i = 0;

	while(i<size)
	{
		while(ptr[i]!=0 || ptr[i+1]!=0) //search 00,00
		{
			*out_ptr++ = ptr[i];
			if(i++==size)
				break;
		}

		if(i++>=size)
			break;
		*out_ptr++ = 0;
		*out_ptr++ = 0;
		if(ptr[++i] <=3)//if the third element after two 0s is less than three
			*out_ptr++ = emulation_prevention_three_byte;
	}

	//in case of cabac_zero_word
	if(ptr[size-1] == 0)
		*out_ptr = emulation_prevention_three_byte;

	out_bs->streambytecnt =  out_ptr - out_bs->bitstream;//update byte count value
}


void hmr_bitstream_put_nal_unit_header(bitstream_t* bs, unsigned int nalu_type, ushort temporal_id, ushort rsvd_zero6bits)
{
	hmr_bitstream_write_bits(bs,0,1);                    // forbidden_zero_bit
	hmr_bitstream_write_bits(bs,nalu_type, 6);  // nal_unit_type
	hmr_bitstream_write_bits(bs,rsvd_zero6bits, 6);                   // nuh_reserved_zero_6bits
	hmr_bitstream_write_bits(bs,temporal_id+1, 3); // nuh_temporal_id_plus1
}


int hmr_bitstream_bitcount(bitstream_t* bs)
{
  return (bs->streambytecnt<<3) + bs->streambitcnt;
}


void hmr_bitstream_write2file(bitstream_t* bs)
{
	FILE *f_annexb;
	if (((f_annexb) = fopen ("prueba.bin", "wb")) == NULL)
	{
		printf("Fatal: cannot open Annex B bytestream file 'prueba.bin', exit (-1)\n");
	}

	fwrite (bs->bitstream, 1, bs->streambytecnt, f_annexb);
	fflush (f_annexb);
	fclose(f_annexb);

}
