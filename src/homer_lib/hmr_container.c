/*****************************************************************************
 * hmr_container.c : homerHEVC encoding library
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

#include <memory.h>

#include "hmr_private.h"
#include "hmr_os_primitives.h"
#include "hmr_container.h"

//-----------------------------------------------SIMPLE hmr_container --------------------------------------------------------------
typedef struct hmr_container_t hmr_container_t;

struct hmr_container_t
{
	void* list[MAX_CONT_ELEMENTS];
	int read;
	int write;
//	int mask;	
	int num_buffs_in_list;
};



void cont_init(void** h)
{		
	*h = (hmr_container_t *)calloc(1,sizeof(hmr_container_t));
}

void cont_delete(void* h)
{
	free(h);
}

//empty hmr_container
void cont_reset(void *cont)
{
	memset(cont,0,sizeof(hmr_container_t));
}

void cont_put(void* h,void *ppbuff)
{
	hmr_container_t *cont = (hmr_container_t *)h;

	cont->list[cont->write++]=ppbuff;
	cont->write &= MAX_CONT_ELEMENTS_MASK;
	cont->num_buffs_in_list++;
//	*ppbuff = NULL;
}

void cont_get(void* h,void **ppbuff)
{
	hmr_container_t *cont = (hmr_container_t *)h;

	*ppbuff=cont->list[cont->read];
	cont->list[cont->read++] = NULL;
	cont->read &= MAX_CONT_ELEMENTS_MASK;
	cont->num_buffs_in_list--;
}

int get_num_elements(void* h)
{
	hmr_container_t *cont = (hmr_container_t *)h;
	return (cont->num_buffs_in_list);
}



//------------------------------------------------------ SYNC hmr_container ------------------------------------------------
//in this hmr_container it is supposed to be called only by one writer and one reader
typedef struct sync_list_t sync_list_t;
struct sync_list_t
{
	void* list[MAX_CONT_ELEMENTS];
	int read;
	int write;
	hmr_sem_t semaphore;
	hmr_sem_ptr sem;
	int num_buffs_in_list;
};

typedef struct sync_hmr_container_t sync_hmr_container_t;
struct sync_hmr_container_t
{
	sync_list_t buffs_filled;
	sync_list_t buffs_clean;
};

void sync_cont_init(void** h)
{
	sync_hmr_container_t* cont = (sync_hmr_container_t *)calloc(1,sizeof(sync_hmr_container_t));
	SEM_INIT(cont->buffs_clean.semaphore, 0,MAX_CONT_ELEMENTS);
	SEM_COPY(cont->buffs_clean.semaphore, cont->buffs_clean.sem);
	SEM_INIT(cont->buffs_filled.semaphore, 0,MAX_CONT_ELEMENTS);
	SEM_COPY(cont->buffs_filled.semaphore, cont->buffs_filled.sem);

	*h = cont; 
}

void sync_cont_delete(void* h)
{
	sync_hmr_container_t * cont = (sync_hmr_container_t *)h;
	SEM_DESTROY(cont->buffs_clean.sem);
	SEM_DESTROY(cont->buffs_filled.sem);

	free(h);
}

void sync_cont_reset(void* h)
{
	sync_hmr_container_t * cont = (sync_hmr_container_t *)h;
	memset(cont->buffs_clean.list,0,MAX_CONT_ELEMENTS);
	cont->buffs_clean.num_buffs_in_list = 0;
	cont->buffs_clean.read = 0;
	cont->buffs_clean.write = 0;
	cont->buffs_filled.num_buffs_in_list = 0;
	cont->buffs_filled.read = 0;
	cont->buffs_filled.write = 0;
	SEM_RESET(cont->buffs_clean.sem);
	SEM_RESET(cont->buffs_filled.sem);
}

void sync_cont_put_empty(void* h,void *p)
{
	sync_hmr_container_t *cont = (sync_hmr_container_t *)h;
	cont->buffs_clean.list[cont->buffs_clean.write++] = p;	
	cont->buffs_clean.write &= MAX_CONT_ELEMENTS_MASK;
	SEM_POST(cont->buffs_clean.sem);
}

void sync_cont_get_empty(void* h,void **p)
{
	sync_hmr_container_t *cont = (sync_hmr_container_t *)h;
	SEM_WAIT(cont->buffs_clean.sem);
	*p = cont->buffs_clean.list[cont->buffs_clean.read++];
	cont->buffs_clean.read &= MAX_CONT_ELEMENTS_MASK;
}

void sync_cont_put_filled(void* h,void *p)
{
	sync_hmr_container_t *cont = (sync_hmr_container_t *)h;
	cont->buffs_filled.list[cont->buffs_filled.write++] = p;	
	cont->buffs_filled.write &= MAX_CONT_ELEMENTS_MASK;
	SEM_POST(cont->buffs_filled.sem);
}

void sync_cont_get_filled(void* h,void **p)
{
	sync_hmr_container_t *cont = (sync_hmr_container_t *)h;
	SEM_WAIT(cont->buffs_filled.sem);
	*p = cont->buffs_filled.list[cont->buffs_filled.read++];
	cont->buffs_filled.read &= MAX_CONT_ELEMENTS_MASK;
}

int sync_cont_is_empty(void* h)
{
	sync_hmr_container_t *cont = (sync_hmr_container_t *)h;
	return cont->buffs_filled.read == cont->buffs_filled.write;
}
