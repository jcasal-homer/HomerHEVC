/*****************************************************************************
 * hmr_os_primitives.h : homerHEVC encoding library
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

#ifndef __HOMER_HEVC_OS_PRIMITIVES_H__
#define __HOMER_HEVC_OS_PRIMITIVES_H__


#ifdef _MSC_VER
#include	<windows.h>
#include	<intrin.h>//for cpuid in windows
#include	"stdint.h"
//-----------------------------------------------------------threads---------------------------------------------------

//data alignment
#define ALIGN(a)	__declspec(align(a))

//get cpu info
#define GET_CPU_ID(cpu_info)	__cpuid(cpu_info, 1);

//-----------------------------------------------------------mutex---------------------------------------------------
typedef CRITICAL_SECTION hmr_mutex;

#define MUTEX_INIT(mutex)									InitializeCriticalSectionAndSpinCount(&mutex, 0)
#define MUTEX_CLOSE(mutex)									DeleteCriticalSection(&mutex)
#define MUTEX_LOCK(mutex)									EnterCriticalSection(&mutex)
#define MUTEX_UNLOCK(mutex)									LeaveCriticalSection(&mutex)
//-----------------------------------------------------------semaphores---------------------------------------------------

typedef void* hmr_sem_t;
typedef void* hmr_sem_ptr;

#define SEM_COPY(a,b) a=b

#define SEM_INIT(sem, count, max_count)									sem = CreateSemaphore(NULL,count,max_count,NULL);
#define SEM_POST(sem)													ReleaseSemaphore(sem,1,NULL);
#define SEM_WAIT(sem)													WaitForSingleObject(sem, INFINITE)
#define SEM_WAIT_MULTIPLE(semaphores, n)								WaitForMultipleObjects(n, semaphores, TRUE, INFINITE)
#define SEM_RESET(sem)																		\
{																							\
	while(WaitForSingleObject(sem, 0)!=WAIT_TIMEOUT);										\
}

#define SEM_DESTROY(sem)												CloseHandle(sem);

//-----------------------------------------------------------thread---------------------------------------------------

//thread function return format
#define THREAD_RETURN_TYPE	void
#define THREAD_RETURN	

//thread handle format
typedef void* hmr_thread_t;

#define CREATE_THREAD(thread, func, param)					thread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func,	(LPVOID)param, 0, NULL);					\

#define CREATE_THREADS(threads, func, param, nthreads)																			\
{																																\
	int nFork;																													\
	for(nFork=0;nFork<nthreads-1;nFork++)																						\
		CREATE_THREAD(threads[nFork], func, param[nFork]);																		\
	func(param[nFork]);																											\
}

#define JOIN_THREAD(thread)	SEM_WAIT(thread);CloseHandle(thread);						

#define JOIN_THREADS(threads, nthreads)									\
{																		\
	int nFork;															\
	SEM_WAIT_MULTIPLE(threads, nthreads);								\
	for(nFork=0;nFork<nthreads;nFork++)									\
		CloseHandle(threads[nFork]);									\
}


#else	//#elif defined(__GNUC__) || defined(__clang__)


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <pthread.h>
#include <semaphore.h>

//-----------------------------------------------------------threads---------------------------------------------------

//data alignment
#define ALIGN(a)	__attribute__((aligned(a)))

//get cpu info
#define GET_CPU_ID(cpu_info)																			\
{																									\
	int a,b,c,d;																					\
	__asm__ __volatile__ ("cpuid":	"=a" (a), "=b" (b), "=c" (c), "=d" (d) : "a" (1));				\
	cpu_info[0] = a;																				\
	cpu_info[1] = b;																				\
	cpu_info[2] = c;																				\
	cpu_info[3] = d;																				\
}

//-----------------------------------------------------------mutex---------------------------------------------------
typedef pthread_mutex_t hmr_mutex;

#define MUTEX_INIT(mutex)									pthread_mutex_init(&mutex, NULL)
#define MUTEX_CLOSE(mutex)									pthread_mutex_destroy(&mutex)
#define MUTEX_LOCK(mutex)									pthread_mutex_lock(&mutex)
#define MUTEX_UNLOCK(mutex)									pthread_mutex_unlock(&mutex)

//-----------------------------------------------------------semaphores---------------------------------------------------

typedef sem_t 	hmr_sem_t;
typedef sem_t* 	hmr_sem_ptr;

#define SEM_COPY(a,b) a=&b
#define SEM_INIT(sem, count, max_count)									sem_init(&sem, 0, count)
#define SEM_POST(sem)													sem_post(sem);
#define SEM_WAIT(sem)													sem_wait(sem)

#define SEM_WAIT_MULTIPLE(semaphores, n)								\
{																		\
	int nSems;															\
	for(nSems=0;nSems<n;nSems++)										\
		SEM_WAIT(semaphores[nSems]);									\
}

#define SEM_RESET(sem)																											\
{																																\
	int rc=0;																													\
	while(rc==0 && errno != EAGAIN)																								\
		rc = sem_trywait(sem);																									\
}

#define SEM_DESTROY(sem)												sem_destroy(sem);


//-----------------------------------------------------------threads---------------------------------------------------
//thread function return format
#define THREAD_RETURN_TYPE	void *
#define THREAD_RETURN	0

//thread handle format
typedef pthread_t hmr_thread_t;

#define CREATE_THREAD(thread, func, param)					pthread_create( &thread, NULL, func, (void*) param);


#define CREATE_THREADS(handle, func, param, nthreads)																			\
{																																\
	int nFork;																													\
	for(nFork=0;nFork<nthreads-1;nFork++)																						\
	{																															\
		pthread_create( &handle[nFork], NULL, func, (void*) param[nFork]);														\
	}																															\
	func(param[nFork]);																											\
}

#define JOIN_THREAD(thread)				pthread_join(thread, NULL);

#define JOIN_THREADS(thread, nthreads)									\
{																		\
	int nFork;															\
	for(nFork=0;nFork<nthreads;nFork++)									\
		JOIN_THREAD(thread[nFork]);										\
}

#endif	//_MSC_VER


#endif /* __HOMER_HEVC_OS_PRIMITIVES_H__ */
