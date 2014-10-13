/*****************************************************************************
 * hmr_rate_control.c : homerHEVC encoding library
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

#include "hmr_private.h"
#include "hmr_common.h"
#include <math.h>




void hmr_rc_init(hvenc_t* ed)
{
	ed->rc.vbv_size = ed->vbv_size*1000;
	ed->rc.vbv_fullness = ed->vbv_init*1000;
	ed->rc.average_pict_size = ed->bitrate/ed->frame_rate;
	ed->rc.average_bits_per_ctu = ed->rc.average_pict_size/ed->pict_total_ctu;
}

void hmr_rc_init_seq(hvenc_t* ed)//seq goes from intra to intra
{
/*	ed->r = (int)floor(2.0*ed->bit_rate/ed->frame_rate + 0.5);//floor function returns a floating-point value representing the largest integer that is less than or equal to x

	// average activity 
	ed->avg_act = 1000;//pq el bitrate se mide en unidades de 400 bps
	ed->lastIntraFrame = -20;

	ed->d=(ed->vbv_buffer_size*512);//*1024 = vbv_buff
	ed->dAcc = 0;//(ed->vbv_buffer_size*256);
	ed->dAccRestar = 0;//ed->dAcc>>4;
*/
}


void hmr_rc_gop(hvenc_t* ed)//, int np, int nb)
{
//	ed->Np = ed->fieldpic ? 2*np+1 : np;
//	ed->Nb = ed->fieldpic ? 2*nb : nb;
}


void hmr_rc_init_pic(hvenc_t* ed)
{
	int ithreads;

	ed->rc.consumed_bitrate = 0;
	ed->rc.consumed_ctus = 0;

	switch(ed->slice_type)
	{
	case  I_SLICE:
		ed->rc.target_pict_size = ed->rc.average_bits_per_ctu*sqrt((double)ed->intra_period);
		break;
	case  P_SLICE:
		ed->rc.target_pict_size = ed->rc.average_bits_per_ctu;
		break;	
	case  B_SLICE:
		ed->rc.target_pict_size = ed->rc.average_bits_per_ctu/2;
		break;	
	}

	for(ithreads=0;ithreads<ed->wfpp_num_threads;ithreads++)
	{
		henc_thread_t* henc_th = ed->thread[ithreads];
		
		henc_th->target_pict_size = ed->rc.target_pict_size;
		henc_th->num_encoded_ctus = 0;
		henc_th->num_bits = 0;
	}
}



void hmr_rc_end_pic(hvenc_t* ed)
{
	int ithreads;

	for(ithreads=0;ithreads<ed->wfpp_num_threads;ithreads++)
	{
		henc_thread_t* henc_th = ed->thread[ithreads];
		
		ed->rc.consumed_bitrate += henc_th->num_bits;
		ed->rc.consumed_ctus += henc_th->num_encoded_ctus;
	}

	ed->rc.vbv_fullness += ed->rc.average_pict_size;
	ed->rc.vbv_fullness -= ed->rc.consumed_bitrate;

	if(ed->rc.vbv_fullness>ed->rc.vbv_size)
	{
		printf("HomerHEVC - vbv_overflow: bitrate consumed is lower than expected\r\n");
		ed->rc.vbv_fullness>ed->rc.vbv_size;
	}

	if(ed->rc.vbv_fullness<0)
	{
		printf("HomerHEVC - vbv_underflow: bitrate consumed is higher than expected\r\n");
		ed->rc.vbv_fullness>ed->rc.vbv_size;
	}
}


int hmr_rc_calc_cu_qp(henc_thread_t* curr_thread)
{
	hvenc_t* ed = curr_thread->ed;
	int ithreads;
	double buff_corrector, entropy_corrector;

	for(ithreads=0;ithreads<ed->wfpp_num_threads;ithreads++)
	{
		henc_thread_t* henc_th = ed->thread[ithreads];
		
		ed->rc.consumed_bitrate += henc_th->num_bits;
		ed->rc.consumed_ctus += henc_th->num_encoded_ctus;
	}
	
	buff_corrector = clip((ed->rc.vbv_fullness+(ed->rc.average_bits_per_ctu*ed->rc.consumed_ctus)-ed->rc.consumed_bitrate)/ed->rc.vbv_size, 0.0, 1.0);

	return (int)(buff_corrector*MAX_QP);
}


int hmr_rc_get_cu_qp(henc_thread_t* et, cu_partition_info_t	*curr_cu_info)
{
	uint qp;
	if(et->ed->bitrate_mode == BR_FIXED_QP)
	{
		qp = et->ed->pict_qp;
	}
	else//cbr, vbr
	{
		if(et->ed->qp_depth==0 )
			qp = hmr_rc_calc_cu_qp(et);
		else if(et->ed->qp_depth <= curr_cu_info->depth)
			qp = hmr_rc_calc_cu_qp(et);
		else
			qp = curr_cu_info->parent->qp;
	}
	return qp;
}
