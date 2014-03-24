/*****************************************************************************
 * hmr_metrics.c : homerHEVC encoding library
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
/*
 * some of the work below is derived from HM HEVC reference software where 
 * the following license applies
 */
/****************************************************************************
/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.  
 *
 * Copyright (c) 2010-2014, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *****************************************************************************/

#include "hmr_private.h"
#include "hmr_common.h"
#include <math.h>


void homer_psnr(picture_t *orig, wnd_t* decoded, int pic_width[3], int pic_height[3], double psnr[3])
{
	int     x, y;

	uint64_t ssd[NUM_PICT_COMPONENTS]={0,0,0};

	double psnry  = 0;
	double  psnru  = 0;
	double  psnrv  = 0;

	int width;
	int height;

	int size_luma = pic_width[Y_COMP]*pic_height[Y_COMP];
	int	size_chroma   = pic_width[U_COMP]*pic_height[U_COMP];

	double	luma_ref = (double) 255 * 255 * size_luma;
	double	chroma_ref = (double) 255 * 255 * size_chroma;
	int component;

	for (component = 0; component < NUM_PICT_COMPONENTS; component++)
	{

		//===== calculate PSNR =====
		uint8_t*  src = WND_DATA_PTR(byte*, orig->img2encode->img, component);
		uint8_t*  dec = WND_DATA_PTR(byte*, *decoded, component);
		int   dec_stride = WND_STRIDE_2D(*decoded, component);
		int   src_stride;// = width;

		width  = pic_width[component];
		height = pic_height[component];

		src_stride = width;

		for( y = 0; y < height; y++ )
		{
			for( x = 0; x < width; x++ )
			{
				int iDiff = (int)( src[x] - dec[x] );
				ssd[component]   += iDiff * iDiff;
			}
			src += src_stride;
			dec += dec_stride;
		}
	}

	psnry            = ( ssd[Y_COMP] ? 10.0 * log10( luma_ref / (double)ssd[Y_COMP] ) : 99.99 );
	psnru            = ( ssd[U_COMP] ? 10.0 * log10( chroma_ref / (double)ssd[U_COMP] ) : 99.99 );
	psnrv            = ( ssd[V_COMP] ? 10.0 * log10( chroma_ref / (double)ssd[V_COMP] ) : 99.99 );

	psnr[Y_COMP]= psnry;
	psnr[U_COMP]= psnru;
	psnr[V_COMP]= psnrv;

}
