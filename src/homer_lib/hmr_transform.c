/*****************************************************************************
 * hmr_transform.c : homerHEVC encoding library
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
 * some of the work below is derived from HM HEVC reference code where 
 * the following license applies
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

#include <memory.h>
#include <limits.h>

#include "hmr_common.h"
#include "hmr_private.h"


const short g_aiT4[4][4] =
{
  { 64, 64, 64, 64},
  { 83, 36,-36,-83},
  { 64,-64,-64, 64},
  { 36,-83, 83,-36}
};

const short g_aiT8[8][8] =
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

const short g_aiT16[16][16] =
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

const short g_aiT32[32][32] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},//0
  { 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13,  4, -4,-13,-22,-31,-38,-46,-54,-61,-67,-73,-78,-82,-85,-88,-90,-90},//1
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90,-90,-87,-80,-70,-57,-43,-25, -9,  9, 25, 43, 57, 70, 80, 87, 90},//2
  { 90, 82, 67, 46, 22, -4,-31,-54,-73,-85,-90,-88,-78,-61,-38,-13, 13, 38, 61, 78, 88, 90, 85, 73, 54, 31,  4,-22,-46,-67,-82,-90},//3
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89, 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},//4
  { 88, 67, 31,-13,-54,-82,-90,-78,-46, -4, 38, 73, 90, 85, 61, 22,-22,-61,-85,-90,-73,-38,  4, 46, 78, 90, 82, 54, 13,-31,-67,-88},//5
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87,-87,-57, -9, 43, 80, 90, 70, 25,-25,-70,-90,-80,-43,  9, 57, 87},//6
  { 85, 46,-13,-67,-90,-73,-22, 38, 82, 88, 54, -4,-61,-90,-78,-31, 31, 78, 90, 61,  4,-54,-88,-82,-38, 22, 73, 90, 67, 13,-46,-85},//7
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},//8
  { 82, 22,-54,-90,-61, 13, 78, 85, 31,-46,-90,-67,  4, 73, 88, 38,-38,-88,-73, -4, 67, 90, 46,-31,-85,-78,-13, 61, 90, 54,-22,-82},//9
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80,-80, -9, 70, 87, 25,-57,-90,-43, 43, 90, 57,-25,-87,-70,  9, 80},//10
  { 78, -4,-82,-73, 13, 85, 67,-22,-88,-61, 31, 90, 54,-38,-90,-46, 46, 90, 38,-54,-90,-31, 61, 88, 22,-67,-85,-13, 73, 82,  4,-78},//11
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75, 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},//12
  { 73,-31,-90,-22, 78, 67,-38,-90,-13, 82, 61,-46,-88, -4, 85, 54,-54,-85,  4, 88, 46,-61,-82, 13, 90, 38,-67,-78, 22, 90, 31,-73},//13
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70,-70, 43, 87, -9,-90,-25, 80, 57,-57,-80, 25, 90,  9,-87,-43, 70},//14
  { 67,-54,-78, 38, 85,-22,-90,  4, 90, 13,-88,-31, 82, 46,-73,-61, 61, 73,-46,-82, 31, 88,-13,-90, -4, 90, 22,-85,-38, 78, 54,-67},//15
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},//16
  { 61,-73,-46, 82, 31,-88,-13, 90, -4,-90, 22, 85,-38,-78, 54, 67,-67,-54, 78, 38,-85,-22, 90,  4,-90, 13, 88,-31,-82, 46, 73,-61},//17
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57,-57, 80, 25,-90,  9, 87,-43,-70, 70, 43,-87, -9, 90,-25,-80, 57},//18
  { 54,-85, -4, 88,-46,-61, 82, 13,-90, 38, 67,-78,-22, 90,-31,-73, 73, 31,-90, 22, 78,-67,-38, 90,-13,-82, 61, 46,-88,  4, 85,-54},//19
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50, 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},//20
  { 46,-90, 38, 54,-90, 31, 61,-88, 22, 67,-85, 13, 73,-82,  4, 78,-78, -4, 82,-73,-13, 85,-67,-22, 88,-61,-31, 90,-54,-38, 90,-46},//21
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43,-43, 90,-57,-25, 87,-70, -9, 80,-80,  9, 70,-87, 25, 57,-90, 43},//22
  { 38,-88, 73, -4,-67, 90,-46,-31, 85,-78, 13, 61,-90, 54, 22,-82, 82,-22,-54, 90,-61,-13, 78,-85, 31, 46,-90, 67,  4,-73, 88,-38},//23
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},//24
  { 31,-78, 90,-61,  4, 54,-88, 82,-38,-22, 73,-90, 67,-13,-46, 85,-85, 46, 13,-67, 90,-73, 22, 38,-82, 88,-54, -4, 61,-90, 78,-31},//25
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25,-25, 70,-90, 80,-43, -9, 57,-87, 87,-57,  9, 43,-80, 90,-70, 25},//26
  { 22,-61, 85,-90, 73,-38, -4, 46,-78, 90,-82, 54,-13,-31, 67,-88, 88,-67, 31, 13,-54, 82,-90, 78,-46,  4, 38,-73, 90,-85, 61,-22},//27
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18, 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},//28
  { 13,-38, 61,-78, 88,-90, 85,-73, 54,-31,  4, 22,-46, 67,-82, 90,-90, 82,-67, 46,-22, -4, 31,-54, 73,-85, 90,-88, 78,-61, 38,-13},//29
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9, -9, 25,-43, 57,-70, 80,-87, 90,-90, 87,-80, 70,-57, 43,-25,  9},//30
  {  4,-13, 22,-31, 38,-46, 54,-61, 67,-73, 78,-82, 85,-88, 90,-90, 90,-90, 88,-85, 82,-78, 73,-67, 61,-54, 46,-38, 31,-22, 13, -4} //31
};


// Fast DST Algorithm. Full matrix multiplication for DST and Fast DST algorithm 
// give identical results
void fastForwardDst(short *block, short *coeff, int shift, int src_stride)  // input block, output coeff
{
  int i, c[4];
  int rnd_factor = 1<<(shift-1);
  for (i=0; i<4; i++)
  {
    // Intermediate Variables
    c[0] = block[src_stride*i+0] + block[src_stride*i+3];
    c[1] = block[src_stride*i+1] + block[src_stride*i+3];
    c[2] = block[src_stride*i+0] - block[src_stride*i+1];
    c[3] = 74* block[src_stride*i+2];

    coeff[   i] =  ( 29 * c[0] + 55 * c[1]         + c[3]               + rnd_factor ) >> shift;
    coeff[ 4+i] =  ( 74 * (block[src_stride*i+0]+ block[src_stride*i+1] - block[src_stride*i+3])   + rnd_factor ) >> shift;
    coeff[ 8+i] =  ( 29 * c[2] + 55 * c[0]         - c[3]               + rnd_factor ) >> shift;
    coeff[12+i] =  ( 55 * c[2] - 29 * c[1]         + c[3]               + rnd_factor ) >> shift;
  }
}

void fastInverseDst(short *tmp, short *block, int shift, int dst_stride)  // input tmp, output block
{
  int i, c[4];
  int rnd_factor = 1<<(shift-1);
  for (i=0; i<4; i++)
  {  
    // Intermediate Variables
    c[0] = tmp[  i] + tmp[ 8+i];
    c[1] = tmp[8+i] + tmp[12+i];
    c[2] = tmp[  i] - tmp[12+i];
    c[3] = 74* tmp[4+i];

    block[dst_stride*i+0] = clip(( 29 * c[0] + 55 * c[1]     + c[3]		+ rnd_factor ) >> shift, -32768, 32767);
    block[dst_stride*i+1] = clip(( 55 * c[2] - 29 * c[1]     + c[3]		+ rnd_factor ) >> shift, -32768, 32767);
    block[dst_stride*i+2] = clip(( 74 * (tmp[i] - tmp[8+i]  + tmp[12+i]) + rnd_factor ) >> shift, -32768, 32767);
    block[dst_stride*i+3] = clip(( 55 * c[0] + 29 * c[2]     - c[3]      + rnd_factor ) >> shift, -32768, 32767);
  }
}


void partialButterfly4(short *src, short *dst,int shift, int src_stride)//(Short *src,Short *dst,Int shift, Int line)
{
  int j;
  int E[2],O[2];
  int add = 1<<(shift-1);

  for (j=0; j<4; j++)
  {    
    /* E and O */
    E[0] = src[0] + src[3];
    O[0] = src[0] - src[3];
    E[1] = src[1] + src[2];
    O[1] = src[1] - src[2];

    dst[0] = (g_aiT4[0][0]*E[0] + g_aiT4[0][1]*E[1] + add)>>shift;
    dst[2*4] = (g_aiT4[2][0]*E[0] + g_aiT4[2][1]*E[1] + add)>>shift;
    dst[4] = (g_aiT4[1][0]*O[0] + g_aiT4[1][1]*O[1] + add)>>shift;
    dst[3*4] = (g_aiT4[3][0]*O[0] + g_aiT4[3][1]*O[1] + add)>>shift;

    src += src_stride;
    dst++;
  }
}

void partialButterflyInverse4(short *src, short *dst,int shift, int dst_stride)//(Short *src,Short *dst,Int shift, Int line)
{
  int j;
  int E[2],O[2];
  int add = 1<<(shift-1);

  for (j=0; j<4; j++)
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */    
    O[0] = g_aiT4[1][0]*src[4] + g_aiT4[3][0]*src[3*4];
    O[1] = g_aiT4[1][1]*src[4] + g_aiT4[3][1]*src[3*4];
    E[0] = g_aiT4[0][0]*src[0] + g_aiT4[2][0]*src[2*4];
    E[1] = g_aiT4[0][1]*src[0] + g_aiT4[2][1]*src[2*4];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    dst[0] = clip((E[0] + O[0] + add)>>shift, -32768, 32767);
    dst[1] = clip((E[1] + O[1] + add)>>shift, -32768, 32767);
    dst[2] = clip((E[1] - O[1] + add)>>shift, -32768, 32767);
    dst[3] = clip((E[0] - O[0] + add)>>shift, -32768, 32767);
            
    src++;
    dst += dst_stride;
  }
}

void partialButterfly8(short *src, short *dst,int shift, int src_stride)//(Short *src,Short *dst,int shift, int 8)
{
  int j,k;
  int E[4],O[4];
  int EE[2],EO[2];
  int add = 1<<(shift-1);

  for (j=0; j<8; j++)
  {  
    /* E and O*/
    for (k=0;k<4;k++)
    {
      E[k] = src[k] + src[7-k];
      O[k] = src[k] - src[7-k];
    }    
    /* EE and EO */
    EE[0] = E[0] + E[3];    
    EO[0] = E[0] - E[3];
    EE[1] = E[1] + E[2];
    EO[1] = E[1] - E[2];

    dst[0] = (g_aiT8[0][0]*EE[0] + g_aiT8[0][1]*EE[1] + add)>>shift;
    dst[4*8] = (g_aiT8[4][0]*EE[0] + g_aiT8[4][1]*EE[1] + add)>>shift; 
    dst[2*8] = (g_aiT8[2][0]*EO[0] + g_aiT8[2][1]*EO[1] + add)>>shift;
    dst[6*8] = (g_aiT8[6][0]*EO[0] + g_aiT8[6][1]*EO[1] + add)>>shift; 

    dst[8] = (g_aiT8[1][0]*O[0] + g_aiT8[1][1]*O[1] + g_aiT8[1][2]*O[2] + g_aiT8[1][3]*O[3] + add)>>shift;
    dst[3*8] = (g_aiT8[3][0]*O[0] + g_aiT8[3][1]*O[1] + g_aiT8[3][2]*O[2] + g_aiT8[3][3]*O[3] + add)>>shift;
    dst[5*8] = (g_aiT8[5][0]*O[0] + g_aiT8[5][1]*O[1] + g_aiT8[5][2]*O[2] + g_aiT8[5][3]*O[3] + add)>>shift;
    dst[7*8] = (g_aiT8[7][0]*O[0] + g_aiT8[7][1]*O[1] + g_aiT8[7][2]*O[2] + g_aiT8[7][3]*O[3] + add)>>shift;

    src += src_stride;//8;
    dst++;
  }
}


void partialButterflyInverse8(short *src, short *dst,int shift, int dst_stride)//(Short *src,Short *dst,int shift, int line)
{
  int j,k;
  int E[4],O[4];
  int EE[2],EO[2];
  int add = 1<<(shift-1);

  for (j=0; j<8; j++) 
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<4;k++)
    {
      O[k] = g_aiT8[ 1][k]*src[8] + g_aiT8[ 3][k]*src[3*8] + g_aiT8[ 5][k]*src[5*8] + g_aiT8[ 7][k]*src[7*8];
    }

    EO[0] = g_aiT8[2][0]*src[ 2*8 ] + g_aiT8[6][0]*src[ 6*8 ];
    EO[1] = g_aiT8[2][1]*src[ 2*8 ] + g_aiT8[6][1]*src[ 6*8 ];
    EE[0] = g_aiT8[0][0]*src[ 0      ] + g_aiT8[4][0]*src[ 4*8 ];
    EE[1] = g_aiT8[0][1]*src[ 0      ] + g_aiT8[4][1]*src[ 4*8 ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */ 
    E[0] = EE[0] + EO[0];
    E[3] = EE[0] - EO[0];
    E[1] = EE[1] + EO[1];
    E[2] = EE[1] - EO[1];
    for (k=0;k<4;k++)
    {
      dst[ k   ] = clip((E[k] + O[k] + add)>>shift, -32768, 32767);
      dst[ k+4 ] = clip((E[3-k] - O[3-k] + add)>>shift, -32768, 32767);
    }   
    src++;
    dst += dst_stride;//8;
  }
}

void partialButterfly16(short *src, short *dst,int shift, int src_stride)//(Short *src,Short *dst,int shift, int line)
{
  int j,k;
  int E[8],O[8];
  int EE[4],EO[4];
  int EEE[2],EEO[2];
  int add = 1<<(shift-1);

  for (j=0; j<16; j++) 
  {    
    /* E and O*/
    for (k=0;k<8;k++)
    {
      E[k] = src[k] + src[15-k];
      O[k] = src[k] - src[15-k];
    } 
    /* EE and EO */
    for (k=0;k<4;k++)
    {
      EE[k] = E[k] + E[7-k];
      EO[k] = E[k] - E[7-k];
    }
    /* EEE and EEO */
    EEE[0] = EE[0] + EE[3];    
    EEO[0] = EE[0] - EE[3];
    EEE[1] = EE[1] + EE[2];
    EEO[1] = EE[1] - EE[2];

    dst[ 0      ] = (g_aiT16[ 0][0]*EEE[0] + g_aiT16[ 0][1]*EEE[1] + add)>>shift;        
    dst[ 8*16 ] = (g_aiT16[ 8][0]*EEE[0] + g_aiT16[ 8][1]*EEE[1] + add)>>shift;    
    dst[ 4*16 ] = (g_aiT16[ 4][0]*EEO[0] + g_aiT16[ 4][1]*EEO[1] + add)>>shift;        
    dst[ 12*16] = (g_aiT16[12][0]*EEO[0] + g_aiT16[12][1]*EEO[1] + add)>>shift;

    for (k=2;k<16;k+=4)
    {
      dst[ k*16 ] = (g_aiT16[k][0]*EO[0] + g_aiT16[k][1]*EO[1] + g_aiT16[k][2]*EO[2] + g_aiT16[k][3]*EO[3] + add)>>shift;      
    }

    for (k=1;k<16;k+=2)
    {
      dst[ k*16 ] = (g_aiT16[k][0]*O[0] + g_aiT16[k][1]*O[1] + g_aiT16[k][2]*O[2] + g_aiT16[k][3]*O[3] + 
        g_aiT16[k][4]*O[4] + g_aiT16[k][5]*O[5] + g_aiT16[k][6]*O[6] + g_aiT16[k][7]*O[7] + add)>>shift;
    }

    src += src_stride;
    dst ++; 

  }
}

void partialButterflyInverse16(short *src, short *dst,int shift, int dst_stride)
{
  int j,k;
  int E[8],O[8];
  int EE[4],EO[4];
  int EEE[2],EEO[2];
  int add = 1<<(shift-1);

  for (j=0; j<16; j++)
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<8;k++)
    {
      O[k] = g_aiT16[ 1][k]*src[ 16] + g_aiT16[ 3][k]*src[ 3*16] + g_aiT16[ 5][k]*src[ 5*16] + g_aiT16[ 7][k]*src[ 7*16] + 
        g_aiT16[ 9][k]*src[ 9*16] + g_aiT16[11][k]*src[11*16] + g_aiT16[13][k]*src[13*16] + g_aiT16[15][k]*src[15*16];
    }
    for (k=0;k<4;k++)
    {
      EO[k] = g_aiT16[ 2][k]*src[ 2*16] + g_aiT16[ 6][k]*src[ 6*16] + g_aiT16[10][k]*src[10*16] + g_aiT16[14][k]*src[14*16];
    }
    EEO[0] = g_aiT16[4][0]*src[ 4*16 ] + g_aiT16[12][0]*src[ 12*16 ];
    EEE[0] = g_aiT16[0][0]*src[ 0      ] + g_aiT16[ 8][0]*src[ 8*16  ];
    EEO[1] = g_aiT16[4][1]*src[ 4*16 ] + g_aiT16[12][1]*src[ 12*16 ];
    EEE[1] = g_aiT16[0][1]*src[ 0      ] + g_aiT16[ 8][1]*src[ 8*16  ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */ 
    for (k=0;k<2;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+2] = EEE[1-k] - EEO[1-k];
    }    
    for (k=0;k<4;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+4] = EE[3-k] - EO[3-k];
    }    
    for (k=0;k<8;k++)
    {
	clip((E[k] + O[k] + add)>>shift , -32768, 32767);
      dst[k]   = clip((E[k] + O[k] + add)>>shift, -32768, 32767);
      dst[k+8] = clip((E[7-k] - O[7-k] + add)>>shift, -32768, 32767);
    }   
    src ++; 
    dst += dst_stride;
  }
}



void partialButterfly32(short *src, short *dst,int shift, int src_stride)
{
  int j,k;
  int E[16],O[16];
  int EE[8],EO[8];
  int EEE[4],EEO[4];
  int EEEE[2],EEEO[2];
  int add = 1<<(shift-1);

  for (j=0; j<32; j++)
  {    
    /* E and O*/
    for (k=0;k<16;k++)
    {
      E[k] = src[k] + src[31-k];
      O[k] = src[k] - src[31-k];
    } 
    /* EE and EO */
    for (k=0;k<8;k++)
    {
      EE[k] = E[k] + E[15-k];
      EO[k] = E[k] - E[15-k];
    }
    /* EEE and EEO */
    for (k=0;k<4;k++)
    {
      EEE[k] = EE[k] + EE[7-k];
      EEO[k] = EE[k] - EE[7-k];
    }
    /* EEEE and EEEO */
    EEEE[0] = EEE[0] + EEE[3];    
    EEEO[0] = EEE[0] - EEE[3];
    EEEE[1] = EEE[1] + EEE[2];
    EEEO[1] = EEE[1] - EEE[2];

    dst[ 0       ] = (g_aiT32[ 0][0]*EEEE[0] + g_aiT32[ 0][1]*EEEE[1] + add)>>shift;
    dst[ 16*32 ] = (g_aiT32[16][0]*EEEE[0] + g_aiT32[16][1]*EEEE[1] + add)>>shift;
    dst[ 8*32  ] = (g_aiT32[ 8][0]*EEEO[0] + g_aiT32[ 8][1]*EEEO[1] + add)>>shift; 
    dst[ 24*32 ] = (g_aiT32[24][0]*EEEO[0] + g_aiT32[24][1]*EEEO[1] + add)>>shift;
    for (k=4;k<32;k+=8)
    {
      dst[ k*32 ] = (g_aiT32[k][0]*EEO[0] + g_aiT32[k][1]*EEO[1] + g_aiT32[k][2]*EEO[2] + g_aiT32[k][3]*EEO[3] + add)>>shift;
    }       
    for (k=2;k<32;k+=4)
    {
      dst[ k*32 ] = (g_aiT32[k][0]*EO[0] + g_aiT32[k][1]*EO[1] + g_aiT32[k][2]*EO[2] + g_aiT32[k][3]*EO[3] + 
        g_aiT32[k][4]*EO[4] + g_aiT32[k][5]*EO[5] + g_aiT32[k][6]*EO[6] + g_aiT32[k][7]*EO[7] + add)>>shift;
    }       
    for (k=1;k<32;k+=2)
    {
      dst[ k*32 ] = (g_aiT32[k][ 0]*O[ 0] + g_aiT32[k][ 1]*O[ 1] + g_aiT32[k][ 2]*O[ 2] + g_aiT32[k][ 3]*O[ 3] + 
        g_aiT32[k][ 4]*O[ 4] + g_aiT32[k][ 5]*O[ 5] + g_aiT32[k][ 6]*O[ 6] + g_aiT32[k][ 7]*O[ 7] +
        g_aiT32[k][ 8]*O[ 8] + g_aiT32[k][ 9]*O[ 9] + g_aiT32[k][10]*O[10] + g_aiT32[k][11]*O[11] + 
        g_aiT32[k][12]*O[12] + g_aiT32[k][13]*O[13] + g_aiT32[k][14]*O[14] + g_aiT32[k][15]*O[15] + add)>>shift;
    }
    src += src_stride;
    dst ++;
  }
}

void partialButterflyInverse32(short *src, short *dst,int shift, int dst_stride)
{
  int j,k;
  int E[16],O[16];
  int EE[8],EO[8];
  int EEE[4],EEO[4];
  int EEEE[2],EEEO[2];
  int add = 1<<(shift-1);

  for (j=0; j<32; j++)
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<16;k++)
    {
      O[k] = g_aiT32[ 1][k]*src[ 32  ] + g_aiT32[ 3][k]*src[ 3*32  ] + g_aiT32[ 5][k]*src[ 5*32  ] + g_aiT32[ 7][k]*src[ 7*32  ] + 
        g_aiT32[ 9][k]*src[ 9*32  ] + g_aiT32[11][k]*src[ 11*32 ] + g_aiT32[13][k]*src[ 13*32 ] + g_aiT32[15][k]*src[ 15*32 ] + 
        g_aiT32[17][k]*src[ 17*32 ] + g_aiT32[19][k]*src[ 19*32 ] + g_aiT32[21][k]*src[ 21*32 ] + g_aiT32[23][k]*src[ 23*32 ] + 
        g_aiT32[25][k]*src[ 25*32 ] + g_aiT32[27][k]*src[ 27*32 ] + g_aiT32[29][k]*src[ 29*32 ] + g_aiT32[31][k]*src[ 31*32 ];
    }
    for (k=0;k<8;k++)
    {
      EO[k] = g_aiT32[ 2][k]*src[ 2*32  ] + g_aiT32[ 6][k]*src[ 6*32  ] + g_aiT32[10][k]*src[ 10*32 ] + g_aiT32[14][k]*src[ 14*32 ] + 
        g_aiT32[18][k]*src[ 18*32 ] + g_aiT32[22][k]*src[ 22*32 ] + g_aiT32[26][k]*src[ 26*32 ] + g_aiT32[30][k]*src[ 30*32 ];
    }
    for (k=0;k<4;k++)
    {
      EEO[k] = g_aiT32[4][k]*src[ 4*32 ] + g_aiT32[12][k]*src[ 12*32 ] + g_aiT32[20][k]*src[ 20*32 ] + g_aiT32[28][k]*src[ 28*32 ];
    }
    EEEO[0] = g_aiT32[8][0]*src[ 8*32 ] + g_aiT32[24][0]*src[ 24*32 ];
    EEEO[1] = g_aiT32[8][1]*src[ 8*32 ] + g_aiT32[24][1]*src[ 24*32 ];
    EEEE[0] = g_aiT32[0][0]*src[ 0      ] + g_aiT32[16][0]*src[ 16*32 ];    
    EEEE[1] = g_aiT32[0][1]*src[ 0      ] + g_aiT32[16][1]*src[ 16*32 ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    EEE[0] = EEEE[0] + EEEO[0];
    EEE[3] = EEEE[0] - EEEO[0];
    EEE[1] = EEEE[1] + EEEO[1];
    EEE[2] = EEEE[1] - EEEO[1];    
    for (k=0;k<4;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+4] = EEE[3-k] - EEO[3-k];
    }    
    for (k=0;k<8;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+8] = EE[7-k] - EO[7-k];
    }    
    for (k=0;k<16;k++)
    {
      dst[k]    = clip((E[k] + O[k] + add)>>shift , -32768, 32767);
      dst[k+16] = clip((E[15-k] - O[15-k] + add)>>shift , -32768, 32767);
    }
    src ++;
    dst += dst_stride;
  }
}





				
void transform(int bitDepth, int16_t *block,int16_t *coeff, int block_size, int iWidth, int iHeight, int width_shift, int height_shift, unsigned short uiMode, int16_t *aux)
{
	int shift_1st = width_shift - 2  + 1 + bitDepth-8;//g_aucConvertToBit[iWidth]  + 1 + bitDepth-8; // log2(iWidth) - 1 + g_bitDepth - 8
	int shift_2nd = height_shift -2  + 8;//g_aucConvertToBit[iHeight]  + 8;                   // log2(iHeight) + 6

	if(iWidth == 4 && iHeight == 4)
	{
		if (uiMode != REG_DCT)
		{
		  fastForwardDst(block, aux, shift_1st, block_size); // Forward DST BY FAST ALGORITHM, block input, tmp output
		  fastForwardDst(aux, coeff, shift_2nd, 4);// Forward DST BY FAST ALGORITHM, tmp input, coeff output
		}
		else
		{
			partialButterfly4( block, aux, shift_1st, block_size);
			partialButterfly4( aux, coeff, shift_2nd, 4);
		}
	}
	else if(iWidth == 8 && iHeight == 8)
	{
		partialButterfly8( block, aux, shift_1st, block_size);
		partialButterfly8( aux, coeff, shift_2nd, 8);
	}
	else if(iWidth == 16 && iHeight == 16)
	{
		partialButterfly16( block, aux, shift_1st, block_size);
		partialButterfly16( aux, coeff, shift_2nd, 16);
	}
	else if (iWidth == 32 && iHeight == 32)
	{
		//utilizamos estas operaciones para pasar de un buffer 2D del CU (64x64) a uno 1D de la particion 32x32
		partialButterfly32( block, aux, shift_1st, block_size);
		partialButterfly32( aux, coeff, shift_2nd, 32);	
	}	

}
#define SHIFT_INV_1ST          7 // Shift after first inverse transform stage
#define SHIFT_INV_2ND         12 // Shift after second inverse transform stage

void itransform(int bitDepth, int16_t *block,int16_t *coeff, int block_size, int iWidth, int iHeight, unsigned int uiMode, int16_t *aux)
{
	int shift_1st = SHIFT_INV_1ST;//g_aucConvertToBit[iWidth]  + 1 + bitDepth-8; // log2(iWidth) - 1 + g_bitDepth - 8
	int shift_2nd = SHIFT_INV_2ND - (bitDepth-8);//g_aucConvertToBit[iHeight]  + 8;                   // log2(iHeight) + 6

	if(iWidth == 4 && iHeight == 4)
	{
		if (uiMode != REG_DCT)
		{
		  fastInverseDst(coeff, aux, shift_1st, 4); // Forward DST BY FAST ALGORITHM, block input, tmp output
		  fastInverseDst(aux, block, shift_2nd, block_size); // Forward DST BY FAST ALGORITHM, tmp input, coeff output
		}
		else
		{
			partialButterflyInverse4(coeff, aux, shift_1st, 4);
			partialButterflyInverse4(aux, block, shift_2nd, block_size);//pasamos a 2D
		}
	}
	else if(iWidth == 8 && iHeight == 8)
	{
	    partialButterflyInverse8(coeff, aux, shift_1st, 8);
		partialButterflyInverse8(aux, block, shift_2nd, block_size);//pasamos a 2D
	}
	else if(iWidth == 16 && iHeight == 16)
	{
	    partialButterflyInverse16(coeff, aux, shift_1st, 16);
		partialButterflyInverse16(aux, block, shift_2nd, block_size);//pasamos a 2D
	}
	else if (iWidth == 32 && iHeight == 32)
	{
	    partialButterflyInverse32(coeff, aux, shift_1st, 32);
		partialButterflyInverse32(aux, block, shift_2nd, block_size);//pasamos a 2D
//		partialButterflyInverse32(aux, block, shift_2nd, 32, iHeight);//lo dejamos lineal
	}
}
