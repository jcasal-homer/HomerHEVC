/*****************************************************************************
 * hmr_ctx_tables.h : homerHEVC encoding library
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
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
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

#ifndef __HOMER_HEVC_CTX_TABLES_H__
#define __HOMER_HEVC_CTX_TABLES_H__


// initial probability for cu_transquant_bypass flag
static const byte
INIT_CU_TRANSQUANT_BYPASS_FLAG[3][NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX] =
{
  { 154 }, 
  { 154 }, 
  { 154 }, 
};

// initial probability for split flag
static const byte 
INIT_SPLIT_FLAG[3][NUM_SPLIT_FLAG_CTX] =  
{
  { 107,  139,  126, },
  { 107,  139,  126, }, 
  { 139,  141,  157, }, 
};

static const byte 
INIT_SKIP_FLAG[3][NUM_SKIP_FLAG_CTX] =  
{
  { 197,  185,  201, }, 
  { 197,  185,  201, }, 
  { CNU,  CNU,  CNU, }, 
};

static const byte
INIT_MERGE_FLAG_EXT[3][NUM_MERGE_FLAG_EXT_CTX] = 
{
  { 154, }, 
  { 110, }, 
  { CNU, }, 
};

static const byte 
INIT_MERGE_IDX_EXT[3][NUM_MERGE_IDX_EXT_CTX] =  
{
  { 137, }, 
  { 122, }, 
  { CNU, }, 
};

static const byte 
INIT_PART_SIZE[3][NUM_PART_SIZE_CTX] =  
{
  { 154,  139,  154,  154 },
  { 154,  139,  154,  154 },
  { 184,  CNU,  CNU,  CNU },
};

/*static const byte 
INIT_CU_AMP_POS[3][NUM_CU_AMP_CTX] =  
{
  { 154, }, 
  { 154, }, 
  { CNU, }, 
};
*/
static const byte 
INIT_PRED_MODE[3][NUM_PRED_MODE_CTX] = 
{
  { 134, }, 
  { 149, }, 
  { CNU, }, 
};

static const byte 
INIT_INTRA_PRED_MODE[3][NUM_ADI_CTX] = 
{
  { 183, }, 
  { 154, }, 
  { 184, }, 
};

static const byte 
INIT_CHROMA_PRED_MODE[3][NUM_CHROMA_PRED_CTX] = 
{
  { 152,  139, }, 
  { 152,  139, }, 
  {  63,  139, }, 
};

static const byte 
INIT_INTER_DIR[3][NUM_INTER_DIR_CTX] = 
{
  {  95,   79,   63,   31,  31, }, 
  {  95,   79,   63,   31,  31, }, 
  { CNU,  CNU,  CNU,  CNU, CNU, }, 
};

static const byte 
INIT_MVD[3][NUM_MV_RES_CTX] =  
{
  { 169,  198, }, 
  { 140,  198, }, 
  { CNU,  CNU, }, 
};

static const byte 
INIT_REF_PIC[3][NUM_REF_NO_CTX] =  
{
  { 153,  153 }, 
  { 153,  153 }, 
  { CNU,  CNU }, 
};

static const byte 
INIT_DQP[3][NUM_DELTA_QP_CTX] = 
{
  { 154,  154,  154, }, 
  { 154,  154,  154, }, 
  { 154,  154,  154, }, 
};

static const byte 
INIT_QT_CBF[3][2*NUM_QT_CBF_CTX] =  
{
  { 153,  111,  CNU,  CNU,   149,   92,  167,  154 },
  { 153,  111,  CNU,  CNU,   149,  107,  167,  154 },
  { 111,  141,  CNU,  CNU,    94,  138,  182,  154 },
};

static const byte 
INIT_QT_ROOT_CBF[3][NUM_QT_ROOT_CBF_CTX] = 
{
  {  79, }, 
  {  79, }, 
  { CNU, }, 
};

static const byte 
INIT_LAST[3][2*NUM_CTX_LAST_FLAG_XY] =  
{
  { 125,  110,  124,  110,   95,   94,  125,  111,  111,   79,  125,  126,  111,  111,   79,
    108,  123,   93,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU, 
  }, 
  { 125,  110,   94,  110,   95,   79,  125,  111,  110,   78,  110,  111,  111,   95,   94,
    108,  123,  108,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,
  }, 
  { 110,  110,  124,  125,  140,  153,  125,  127,  140,  109,  111,  143,  127,  111,   79, 
    108,  123,   63,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU, 
  }, 
};

static const byte 
INIT_SIG_CG_FLAG[3][2 * NUM_SIG_CG_FLAG_CTX] =  
{
  { 121,  140,  
    61,  154, 
  }, 
  { 121,  140, 
    61,  154, 
  }, 
  {  91,  171,  
    134,  141, 
  }, 
};

static const byte 
INIT_SIG_FLAG[3][NUM_SIG_FLAG_CTX] = 
{
  { 170,  154,  139,  153,  139,  123,  123,   63,  124,  166,  183,  140,  136,  153,  154,  166,  183,  140,  136,  153,  154,  166,  183,  140,  136,  153,  154,  170,  153,  138,  138,  122,  121,  122,  121,  167,  151,  183,  140,  151,  183,  140,  }, 
  { 155,  154,  139,  153,  139,  123,  123,   63,  153,  166,  183,  140,  136,  153,  154,  166,  183,  140,  136,  153,  154,  166,  183,  140,  136,  153,  154,  170,  153,  123,  123,  107,  121,  107,  121,  167,  151,  183,  140,  151,  183,  140,  }, 
  { 111,  111,  125,  110,  110,   94,  124,  108,  124,  107,  125,  141,  179,  153,  125,  107,  125,  141,  179,  153,  125,  107,  125,  141,  179,  153,  125,  140,  139,  182,  182,  152,  136,  152,  136,  153,  136,  139,  111,  136,  139,  111,  }, 
};

static const byte 
INIT_ONE_FLAG[3][NUM_ONE_FLAG_CTX] = 
{
  { 154,  196,  167,  167,  154,  152,  167,  182,  182,  134,  149,  136,  153,  121,  136,  122,  169,  208,  166,  167,  154,  152,  167,  182, }, 
  { 154,  196,  196,  167,  154,  152,  167,  182,  182,  134,  149,  136,  153,  121,  136,  137,  169,  194,  166,  167,  154,  167,  137,  182, }, 
  { 140,   92,  137,  138,  140,  152,  138,  139,  153,   74,  149,   92,  139,  107,  122,  152,  140,  179,  166,  182,  140,  227,  122,  197, }, 
};

static const byte 
INIT_ABS_FLAG[3][NUM_ABS_FLAG_CTX] =  
{
  { 107,  167,   91,  107,  107,  167, }, 
  { 107,  167,   91,  122,  107,  167, }, 
  { 138,  153,  136,  167,  152,  152, }, 
};

static const byte 
INIT_MVP_IDX[3][NUM_MVP_IDX_CTX] =  
{
  { 168 },
  { 168 },
  { CNU }, 
};

static const byte 
INIT_SAO_MERGE_FLAG[3][NUM_SAO_MERGE_FLAG_CTX] = 
{
  { 153,  }, 
  { 153,  }, 
  { 153,  }, 
};

static const byte 
INIT_SAO_TYPE_IDX[3][NUM_SAO_TYPE_IDX_CTX] = 
{
  { 160, },
  { 185, },
  { 200, },
};

static const byte
INIT_TRANS_SUBDIV_FLAG[3][NUM_TRANS_SUBDIV_FLAG_CTX] =
{
  { 224,  167,  122, },
  { 124,  138,   94, },
  { 153,  138,  138, },
};

static const byte
INIT_TRANSFORMSKIP_FLAG[3][2*NUM_TRANSFORMSKIP_FLAG_CTX] = 
{
  { 139,  139}, 
  { 139,  139}, 
  { 139,  139}, 
};
//! \}

#endif //__HOMER_HEVC_CTX_TABLES_H__

