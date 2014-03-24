/*****************************************************************************
 * hmr_profiler.c : homerHEVC encoding library
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
#include <stdio.h>
#include "hmr_profiler.h"

#if _MSC_VER
#include <windows.h>
void profiler_reset_counter(profiler_t* p)
{
    LARGE_INTEGER li;

    QueryPerformanceCounter(&li);
    p->init_count = li.QuadPart;
}


void profiler_start(profiler_t* p)
{
	LARGE_INTEGER li;
	profiler_reset_counter(p);
	p->count = 0;
	QueryPerformanceFrequency(&li);
	p->pc_freq = li.QuadPart;
}


void profiler_accumulate(profiler_t* p)
{
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
	p->count += li.QuadPart-p->init_count;
}


double profiler_get_result(profiler_t* p)
{
	return ((double)p->count)/p->pc_freq;	
}


void profiler_print_result(profiler_t* p)
{
	printf("Counter %s, time = %f ms\r\n", p->name, (double)p->count/p->pc_freq);	
}

#else

#include <sys/time.h>

static const unsigned cycles_per_sec = 1000000;


void QueryPerformanceFrequency(int64_t *freq)
{
	*freq = cycles_per_sec;
}

void QueryPerformanceCounter(int64_t *clock)
{
	struct timeval time;

	gettimeofday(&time, NULL);
	*clock = time.tv_usec + time.tv_sec * cycles_per_sec;
}


void profiler_reset_counter(profiler_t* p)
{
	int64_t li;

	QueryPerformanceCounter(&li);
	p->init_count = li;
}


void profiler_start(profiler_t* p)
{
	profiler_reset_counter(p);
	p->count = 0;
}


void profiler_accumulate(profiler_t* p)
{
	int64_t li;

	QueryPerformanceCounter(&li);
	p->count += li-p->init_count;
}


double profiler_get_result(profiler_t* p)
{
	return ((double)p->count)/p->pc_freq;
}


void profiler_print_result(profiler_t* p)
{
	printf("Counter %s, time = %f ms\r\n", p->name, (double)p->count/p->pc_freq);
}

#endif
