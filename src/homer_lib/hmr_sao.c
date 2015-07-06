/*****************************************************************************
* hmr_sao.c : homerHEVC encoding library
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
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
*****************************************************************************/
/*
* some of the work below is derived from HM HEVC reference code where 
* the following license applies
*****************************************************************************
* The copyright in this software is being made available under the BSD
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

#include <math.h>
#include <memory.h>
#include "hmr_common.h"
#include "hmr_private.h"





uint g_saoMaxOffsetQVal[NUM_PICT_COMPONENTS];
uint m_offsetStepLog2[NUM_PICT_COMPONENTS];


int skiped_lines_r[NUM_PICT_COMPONENTS] = {5,3,3};
int skiped_lines_b[NUM_PICT_COMPONENTS] = {4,2,2};


