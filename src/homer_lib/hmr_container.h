/*****************************************************************************
 * hmr_container.h : homerHEVC encoding library
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


#ifndef __HOMER_HEVC_CONTAINER_H__
#define __HOMER_HEVC_CONTAINER_H__



//handles up to 32 elements; has to be 2^n
#define MAX_CONT_ELEMENTS		32
#define MAX_CONT_ELEMENTS_MASK	(MAX_CONT_ELEMENTS-1)

//hmr_container without sync - not atomic. All functions must be called in the same threads
void	cont_init(void** h);
void	cont_delete(void* h);
void	cont_reset(void* cont);
void	cont_put(void* cont,void *ppelement);
void	cont_get(void* cont,void **ppelement);
int		get_num_elements(void* h);

//sync hmr_container - not atomic. Just one producer and one consumer
void sync_cont_init(void** h);
void sync_cont_delete(void* h);
void sync_cont_reset(void* h);
void sync_cont_put_empty(void* h,void *p);
void sync_cont_get_empty(void* h,void **p);
void sync_cont_put_filled(void* h,void *p);
void sync_cont_get_filled(void* h,void **p);

#endif //__HOMER_HEVC_CONTAINER_H__
