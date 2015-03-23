/*****************************************************************************
 * hmr_profiler.h : homerHEVC encoding library
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

#ifndef __HOMER_HEVC_PROFILER_H__
#define __HOMER_HEVC_PROFILER_H__

#include "hmr_os_primitives.h"

//enable - disable profilers on compilation
#ifdef _TIME_PROFILING_
#define PROFILER_INIT(name)		{name,0,0,0}
#define	PROFILER_START(a)		profiler_start(&a);
#define	PROFILER_RESET(a)		profiler_reset_counter(&a);
#define	PROFILER_ACCUMULATE(a)	profiler_accumulate(&a);
#define	PROFILER_PRINT(a)		profiler_print_result(&a);
#else
#define PROFILER_INIT(name)		{name,0,0,0}//#define PROFILER_INIT(name)
#define	PROFILER_START(a)
#define	PROFILER_RESET(a)		
#define	PROFILER_ACCUMULATE(a)	
#define	PROFILER_PRINT(a)
#endif // _TIME_PROFILING_		


typedef struct profiler_t profiler_t;

struct profiler_t
{
	char		name[256];
	int64_t		pc_freq;
	int64_t		init_count;
	int64_t		count;
};


void profiler_print_result(profiler_t* p);
double profiler_get_result(profiler_t* p);
void profiler_accumulate(profiler_t* p);
void profiler_start(profiler_t* p);
void profiler_reset_counter(profiler_t* p);
void init_profiler(profiler_t* p, char* str);

#endif /*__HOMER_HEVC_PROFILER_H__*/