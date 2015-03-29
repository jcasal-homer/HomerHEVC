/*****************************************************************************
 * hmr_encoder_lib.c : homerHEVC encoding library
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
#include <stdio.h>
#include	<malloc.h>
#include	<memory.h>
#include	<math.h>
#include	<limits.h>

#include	"hmr_os_primitives.h"
#include 	"homer_hevc_enc_api.h"
#include 	"hmr_common.h"
//#include	"hmr_tables.h"
#include	"hmr_profiler.h"
#include	"hmr_sse42_functions.h"


static const uint16_t ang_table[] = {0,    2,    5,   9,  13,  17,  21,  26,  32}; 
static const uint16_t inv_ang_table[] = {0, 4096, 1638, 910, 630, 482, 390, 315, 256}; // (256 * 32) / Angle

#ifdef COMPUTE_METRICS
profiler_t frame_metrics = PROFILER_INIT("frame_metrics");
#endif

#ifdef _TIME_PROFILING_
profiler_t intra_luma = PROFILER_INIT("intra_luma");
profiler_t intra_chroma = PROFILER_INIT("intra_chroma");
profiler_t intra = PROFILER_INIT("intra");
profiler_t cabac = PROFILER_INIT("cabac");
//profiling intra_luma
profiler_t intra_luma_bucle1 = PROFILER_INIT("intra_luma.bucle1");
profiler_t intra_luma_bucle1_sad = PROFILER_INIT("intra_luma.bucle1.sad");
profiler_t intra_luma_bucle2 = PROFILER_INIT("intra_luma.bucle2");
profiler_t intra_luma_bucle3 = PROFILER_INIT("intra_luma.bucle3");
profiler_t intra_luma_generate_prediction = PROFILER_INIT("intra_luma.generate_prediction");
profiler_t intra_luma_predict = PROFILER_INIT("intra_luma_predict");
profiler_t intra_luma_tr = PROFILER_INIT("intra_luma.tr");
profiler_t intra_luma_q = PROFILER_INIT("intra_luma.q");
profiler_t intra_luma_iq = PROFILER_INIT("intra_luma.iq");
profiler_t intra_luma_itr = PROFILER_INIT("intra_luma.itr");
profiler_t intra_luma_recon_ssd = PROFILER_INIT("intra_luma.recon+ssd");
#endif




static int num_scaling_list[NUM_SCALING_MODES]={6,6,6,2};

void *HOMER_enc_init()
{
	int i, size, size_index, list_index, qp;
//	hvenc_t* phvenc = (hvenc_t*)calloc(1,sizeof(hvenc_t));
	hvenc_enc_t* hvenc = (hvenc_enc_t*)calloc(1,sizeof(hvenc_enc_t));
	unsigned short* aux_ptr;
	int cpu_info[4];

	//int max_width = 1920, max_height = 1080;

	hvenc->num_encoded_frames = 0;
	hvenc->num_encoder_modules = 0;

//	hvenc->ctu_group_size = MAX_MB_GROUP_SIZE;
	hvenc->ctu_width[0] = hvenc->ctu_height[0] = 64;
	hvenc->ctu_width[1] = hvenc->ctu_width[2] = 32;
	hvenc->ctu_height[1] = hvenc->ctu_height[2] = 32;
	hvenc->bit_depth = 8;

	//---------------------------------- general tables ------------------------------------------------
	//for partition order inside a ctu
	hvenc->abs2raster_table = (unsigned short*)calloc (hvenc->ctu_width[0] * hvenc->ctu_height[0], sizeof(unsigned short));//number of elements in CU
	hvenc->raster2abs_table = (unsigned short*)calloc (hvenc->ctu_width[0] * hvenc->ctu_height[0], sizeof(unsigned short));//number of elements in CU
	aux_ptr = hvenc->abs2raster_table;
	create_abs2raster_tables(&hvenc->abs2raster_table, 5, 1, 0);
	hvenc->abs2raster_table = aux_ptr;
	create_raster2abs_tables( hvenc->abs2raster_table, hvenc->raster2abs_table, hvenc->ctu_width[0], hvenc->ctu_height[0], 5);


	size=2;
	for ( i=0; i<MAX_CU_DEPTHS; i++ ) //scan block size (2x2, ....., 128x128)
	{
		hvenc->scan_pyramid[0][i] = (uint*) hmr_aligned_alloc (size*size, sizeof(uint32_t));
		hvenc->scan_pyramid[1][i] = (uint*) hmr_aligned_alloc (size*size, sizeof(uint32_t));
		hvenc->scan_pyramid[2][i] = (uint*) hmr_aligned_alloc (size*size, sizeof(uint32_t));
		hvenc->scan_pyramid[3][i] = (uint*) hmr_aligned_alloc (size*size, sizeof(uint32_t));
		init_scan_pyramid( hvenc, hvenc->scan_pyramid[0][i], hvenc->scan_pyramid[1][i], hvenc->scan_pyramid[2][i], hvenc->scan_pyramid[3][i], size, size, i);

		size <<= 1;
	}  

	size=4;
	for ( size_index=0; size_index<NUM_SCALING_MODES; size_index++ )//size_index (4x4,8x8,16x16,32x32)
	{
		int ratio = size/min(NUM_MAX_MATRIX_SIZE,size);
		for ( list_index=0; list_index<num_scaling_list[size_index]; list_index++ )//list_index
		{
			short *quant_def_table = get_default_qtable(size_index, list_index);
			for ( qp=0; qp<NUM_SCALING_REM_LISTS; qp++ )//qp
			{
				hvenc->quant_pyramid[size_index][list_index][qp] = (int*) hmr_aligned_alloc (size*size, sizeof(uint32_t));
				hvenc->dequant_pyramid[size_index][list_index][qp] = (int*) hmr_aligned_alloc (size*size, sizeof(uint32_t));
				hvenc->scaling_error_pyramid[size_index][list_index][qp] = (double*) hmr_aligned_alloc (size*size, sizeof(double));
				init_quant_pyramids( hvenc, hvenc->quant_pyramid[size_index][list_index][qp], hvenc->dequant_pyramid[size_index][list_index][qp], hvenc->scaling_error_pyramid[size_index][list_index][qp],
									quant_def_table, size, size, ratio, min(NUM_MAX_MATRIX_SIZE, size), QUANT_DEFAULT_DC, size_index+2, qp);

//				init_flat_quant_pyramids( hvenc, hvenc->quant_pyramid[size_index][list_index][qp], hvenc->dequant_pyramid[size_index][list_index][qp], hvenc->scaling_error_pyramid[size_index][list_index][qp], size*size, size_index+2, qp);
			}  
		}
		size <<= 1;
	}

	for ( qp=0; qp<NUM_SCALING_REM_LISTS; qp++ )//qp
	{
		hvenc->quant_pyramid[SCALING_MODE_32x32][3][qp] = hvenc->quant_pyramid[SCALING_MODE_32x32][1][qp];
		hvenc->dequant_pyramid[SCALING_MODE_32x32][3][qp] = hvenc->dequant_pyramid[SCALING_MODE_32x32][1][qp];
		hvenc->scaling_error_pyramid[SCALING_MODE_32x32][3][qp] = hvenc->scaling_error_pyramid[SCALING_MODE_32x32][1][qp];
	}  

	//angular intra table
	hvenc->ang_table = ang_table;//(ushort*)calloc (9, sizeof(ushort));//number of elements in CU
	hvenc->inv_ang_table = inv_ang_table;//(ushort*)calloc (9, sizeof(ushort));//number of elements in CU

//	memcpy(hvenc->ang_table, ang_table, sizeof(ang_table));
//	memcpy(hvenc->inv_ang_table, inv_ang_table, sizeof(inv_ang_table));

	sync_cont_init(&hvenc->input_hmr_container);	
	cont_init(&hvenc->output_hmr_container);
	cont_init(&hvenc->cont_empty_reference_wnds);

	GET_CPU_ID(cpu_info)
#ifndef COMPUTE_SSE_FUNCS
	cpu_info[2] = 0;
#endif
	if(cpu_info[2] & 0x100000)////
	{
		printf("SSE42 avaliable!!");

		hvenc->funcs.sad = sse_aligned_sad;
		hvenc->funcs.ssd = sse_aligned_ssd;
		hvenc->funcs.ssd16b = sse_aligned_ssd16b;
		hvenc->funcs.modified_variance = sse_modified_variance;
		hvenc->funcs.predict = sse_aligned_predict;
		hvenc->funcs.reconst = sse_aligned_reconst;
		hvenc->funcs.create_intra_planar_prediction = sse_create_intra_planar_prediction;
		hvenc->funcs.create_intra_angular_prediction = sse_create_intra_angular_prediction;
				
		hvenc->funcs.interpolate_luma_m_compensation = sse_interpolate_luma;
		hvenc->funcs.interpolate_chroma_m_compensation = sse_interpolate_chroma;
		hvenc->funcs.interpolate_luma_m_estimation = sse_interpolate_luma;//hmr_fake_interpolate_luma;//

		hvenc->funcs.quant = sse_aligned_quant;
		hvenc->funcs.inv_quant = sse_aligned_inv_quant;

		hvenc->funcs.transform = sse_transform;
		hvenc->funcs.itransform = sse_itransform;
	}
	else
	{
		hvenc->funcs.sad = sad;
		hvenc->funcs.ssd = ssd;
		hvenc->funcs.ssd16b = ssd16b;
		hvenc->funcs.modified_variance = modified_variance;
		hvenc->funcs.predict = predict;
		hvenc->funcs.reconst = reconst;
		hvenc->funcs.create_intra_planar_prediction = create_intra_planar_prediction;
		hvenc->funcs.create_intra_angular_prediction = create_intra_angular_prediction;

		hvenc->funcs.interpolate_luma_m_compensation = hmr_interpolate_luma;
		hvenc->funcs.interpolate_chroma_m_compensation = hmr_interpolate_chroma;
		hvenc->funcs.interpolate_luma_m_estimation = hmr_interpolate_luma;

		hvenc->funcs.quant = quant;
		hvenc->funcs.inv_quant = iquant;

		hvenc->funcs.transform = transform;
		hvenc->funcs.itransform = itransform;
	}

	if(!InitializeCriticalSectionAndSpinCount(&hvenc->CriticalSection, 0))
		return NULL;
	if(!InitializeCriticalSectionAndSpinCount(&hvenc->CriticalSection2, 0))
		return NULL;
	return hvenc;
}


void put_frame_to_encode(hvenc_enc_t* ed, encoder_in_out_t* input_frame)
{
	video_frame_t	*p;
	uint8_t *src, *dst;
	int stride_dst, stride_src;
	int comp, j;

	sync_cont_get_empty(ed->input_hmr_container, (void**)&p);

	p->temp_info.pts = input_frame->pts;
	p->img_type = input_frame->image_type;
	for(comp=Y_COMP;comp<=V_COMP;comp++)
	{
		src = input_frame->stream.streams[comp];
		dst = WND_DATA_PTR(uint8_t*, p->img, comp);
		stride_src = input_frame->stream.data_stride[comp];//,  ed->pict_width[comp];
		stride_dst = WND_STRIDE_2D(p->img, comp);

		for(j=0;j<ed->pict_height[comp];j++)
		{
			memcpy(dst, src, ed->pict_width[comp]);
			src += stride_src;
			dst += stride_dst;
		}
	}
	sync_cont_put_filled(ed->input_hmr_container, p);
}

int get_frame_to_encode(hvenc_enc_t* venc, video_frame_t **picture)
{
	sync_cont_get_filled(venc->input_hmr_container, (void**)picture);

	return (*picture)!=NULL;
}

void put_avaliable_frame(hvenc_enc_t * venc, video_frame_t *picture)
{
	sync_cont_put_empty(venc->input_hmr_container, picture);
}


#define HMR_FREE(a) if(a!=NULL)free(a);(a)=NULL;
#define HMR_ALIGNED_FREE(a) if(a!=NULL)hmr_aligned_free(a);(a)=NULL;

void HOMER_enc_close(void* h)
{
	hvenc_enc_t* hvenc = (hvenc_enc_t*)h;
	hvenc_t* phvenc_mod = (hvenc_t*)h;
	int i,j;
	int ithreads;
	int size_index;
	if(hvenc->run==1)
	{
		hvenc->run = 0;

		if(phvenc_mod->encoder_thread!=NULL)
		{
			sync_cont_put_filled(hvenc->input_hmr_container, NULL);//wake encoder_thread if it is waiting
			JOINT_THREAD(phvenc_mod->encoder_thread);
		}
	}

	HMR_FREE(hvenc->ref_pic_set_list)

	//for all modules

	for(ithreads=0;ithreads<phvenc_mod->wfpp_num_threads;ithreads++)
	{
		henc_thread_t* henc_th = phvenc_mod->thread[ithreads];
		if(henc_th==NULL)
			break;

		HMR_FREE(henc_th->ctu_rd->merge_idx)
		HMR_FREE(henc_th->ctu_rd->merge)
		HMR_FREE(henc_th->ctu_rd->skipped)
		HMR_FREE(henc_th->ctu_rd->pred_mode)
		HMR_FREE(henc_th->ctu_rd->part_size_type)
		HMR_FREE(henc_th->ctu_rd)
	
		HMR_ALIGNED_FREE(henc_th->cbf_buffs_chroma[V_COMP])
		HMR_ALIGNED_FREE(henc_th->cbf_buffs_chroma[U_COMP])

		for(i=0;i<NUM_CBF_BUFFS;i++)
		{
			HMR_ALIGNED_FREE(henc_th->tr_idx_buffs[i])
			HMR_ALIGNED_FREE(henc_th->intra_mode_buffs[V_COMP][i])
			HMR_ALIGNED_FREE(henc_th->intra_mode_buffs[U_COMP][i])
			HMR_ALIGNED_FREE(henc_th->intra_mode_buffs[Y_COMP][i])
			HMR_ALIGNED_FREE(henc_th->cbf_buffs[V_COMP][i])
			HMR_ALIGNED_FREE(henc_th->cbf_buffs[U_COMP][i])
			HMR_ALIGNED_FREE(henc_th->cbf_buffs[Y_COMP][i])			
		}

		HMR_ALIGNED_FREE(henc_th->cabac_aux_buff)

		HMR_FREE(henc_th->deblock_partition_info)

		HMR_ALIGNED_FREE(henc_th->deblock_edge_filter[EDGE_HOR])
		HMR_ALIGNED_FREE(henc_th->deblock_edge_filter[EDGE_VER])
		HMR_ALIGNED_FREE(henc_th->deblock_filter_strength_bs[EDGE_HOR])
		HMR_ALIGNED_FREE(henc_th->deblock_filter_strength_bs[EDGE_VER])


		for(j=0;j<4;j++)
		{
			wnd_delete(&henc_th->filtered_block_temp_wnd[j]);
			for(i=0;i<4;i++)
			{
				wnd_delete(&henc_th->filtered_block_wnd[j][i]);
			}
		}		

		for(i=0;i<NUM_DECODED_WNDS;i++)
			wnd_delete(&henc_th->decoded_mbs_wnd_[i]);

		wnd_delete(&henc_th->itransform_iquant_wnd);

		for(i=0;i<NUM_QUANT_WNDS;i++)
			wnd_delete(&henc_th->transform_quant_wnd_[i]);

		HMR_ALIGNED_FREE(henc_th->aux_buff)

		wnd_delete(&henc_th->residual_dec_wnd);
		wnd_delete(&henc_th->residual_wnd);
		wnd_delete(&henc_th->prediction_wnd);

		HMR_ALIGNED_FREE(henc_th->pred_aux_buff)

		wnd_delete(&henc_th->curr_mbs_wnd);

		HMR_ALIGNED_FREE(henc_th->adi_pred_buff)
		HMR_ALIGNED_FREE(henc_th->adi_filtered_pred_buff)
		HMR_ALIGNED_FREE(henc_th->top_pred_buff)
		HMR_ALIGNED_FREE(henc_th->left_pred_buff)
		HMR_ALIGNED_FREE(henc_th->bottom_pred_buff)
		HMR_ALIGNED_FREE(henc_th->right_pred_buff)

		HMR_FREE(henc_th->partition_info)

		SEM_DESTROY(henc_th->synchro_signal);

		free(henc_th);
	}
	SEM_DESTROY(phvenc_mod->deblock_filter_sem);

	for(i=0;i<phvenc_mod->pict_total_ctu;i++)
	{
		HMR_FREE(phvenc_mod->ctu_info[i].mv_diff_ref_idx[REF_PIC_LIST_0])
		HMR_FREE(phvenc_mod->ctu_info[i].mv_diff[REF_PIC_LIST_0])
		HMR_FREE(phvenc_mod->ctu_info[i].mv_ref_idx[REF_PIC_LIST_0])
		HMR_FREE(phvenc_mod->ctu_info[i].mv_ref[REF_PIC_LIST_0])
		HMR_FREE(phvenc_mod->ctu_info[i].qp)
		HMR_FREE(phvenc_mod->ctu_info[i].merge)
		HMR_FREE(phvenc_mod->ctu_info[i].merge_idx)
		HMR_FREE(phvenc_mod->ctu_info[i].skipped)
		HMR_FREE(phvenc_mod->ctu_info[i].pred_mode)
		HMR_FREE(phvenc_mod->ctu_info[i].part_size_type)
		HMR_FREE(phvenc_mod->ctu_info[i].pred_depth)
		HMR_FREE(phvenc_mod->ctu_info[i].tr_idx)
		HMR_FREE(phvenc_mod->ctu_info[i].inter_mode)
		HMR_FREE(phvenc_mod->ctu_info[i].intra_mode[Y_COMP])
		HMR_FREE(phvenc_mod->ctu_info[i].cbf[Y_COMP])
	}

	HMR_FREE(phvenc_mod->ctu_info)

	for(i=0;i<phvenc_mod->num_ec;i++)
	{
		HMR_FREE(phvenc_mod->ec_list[i].b_ctx)
		HMR_FREE(phvenc_mod->ec_list[i].contexts)
		HMR_FREE(phvenc_mod->ec_list[i].e_ctx)
	}
	HMR_FREE(phvenc_mod->ec_list)

	for(i=0;i<phvenc_mod->num_ee;i++)
	{
		HMR_FREE(phvenc_mod->ee_list[i]->b_ctx)
		HMR_FREE(phvenc_mod->ee_list[i]->contexts)
		HMR_FREE(phvenc_mod->ee_list[i]->e_ctx)
		HMR_FREE(phvenc_mod->ee_list[i]);
	}
	HMR_FREE(phvenc_mod->ee_list)

	for(i=0;i<phvenc_mod->num_sub_streams;i++)
		hmr_bitstream_free(&phvenc_mod->aux_bs[i]);
	HMR_FREE(phvenc_mod->aux_bs)

	HMR_FREE(phvenc_mod->sub_streams_entry_point_list)

	for(i=0;i<2*MAX_NUM_REF;i++)
	{
		wnd_delete(&hvenc->ref_wnds[i].img);
	}

	for(i=0;i<NUM_INPUT_FRAMES;i++)
	{
		wnd_delete(&hvenc->input_frames[i].img);
	}

	//--------------------------------------------------------HOMER_enc_init-----------------------------------------------------------------------------

	for(i=0;i<NUM_OUTPUT_NALUS;i++)
		hmr_bitstream_free(&phvenc_mod->slice_nalu_list[i].bs);

	hmr_bitstream_free(&hvenc->pps_nalu.bs);
	hmr_bitstream_free(&hvenc->sps_nalu.bs);
	hmr_bitstream_free(&hvenc->vps_nalu.bs);
	hmr_bitstream_free(&phvenc_mod->slice_bs);

	sync_cont_delete(hvenc->input_hmr_container);
	cont_delete(phvenc_mod->hvenc->output_hmr_container);
	cont_delete(hvenc->cont_empty_reference_wnds);

//	HMR_FREE(phvenc_mod->ang_table)
//	HMR_FREE(phvenc_mod->inv_ang_table)

	HMR_FREE(phvenc_mod)


	// end for all modules


	for ( size_index=0; size_index<NUM_SCALING_MODES; size_index++ )
	{
		int list_index;
		for ( list_index=0; list_index<num_scaling_list[size_index]; list_index++ )//list_index
		{
			int qp;
			short *quant_def_table = get_default_qtable(size_index, list_index);
			for ( qp=0; qp<NUM_SCALING_REM_LISTS; qp++ )//qp
			{
				HMR_ALIGNED_FREE(hvenc->scaling_error_pyramid[size_index][list_index][qp])
				HMR_ALIGNED_FREE(hvenc->dequant_pyramid[size_index][list_index][qp])
				HMR_ALIGNED_FREE(hvenc->quant_pyramid[size_index][list_index][qp])
			}  
		}
	}

	for ( i=0; i<MAX_CU_DEPTHS; i++ ) 
	{
		HMR_ALIGNED_FREE(hvenc->scan_pyramid[3][i])
		HMR_ALIGNED_FREE(hvenc->scan_pyramid[2][i])
		HMR_ALIGNED_FREE(hvenc->scan_pyramid[1][i])
		HMR_ALIGNED_FREE(hvenc->scan_pyramid[0][i])
	}  

	HMR_FREE(hvenc->raster2abs_table)
	HMR_FREE(hvenc->abs2raster_table)

	HMR_FREE(hvenc)
}



static const int pad_unit_x[]={1,2,2,1};
static int pad_unit_y[]={1,2,1,1};

int HOMER_enc_control(void *h, int cmd, void *in)
{
	hvenc_enc_t* hvenc = (hvenc_enc_t*)h;
	int err=0;
	int i, aux;
	unsigned short* aux_ptr;
//	int cpu_info[4];

	switch (cmd)
	{
		case HENC_SETCFG :
		{
			HVENC_Cfg *cfg = (HVENC_Cfg *)in;
			int n_enc_mod;
			hvenc_t*  phvenc_mod;
			int num_merge_candidates = 2;

#ifdef COMPUTE_AS_HM
				cfg->rd_mode = RD_DIST_ONLY;    //0 only distortion 
				cfg->bitrate_mode = BR_FIXED_QP;//0=fixed qp, 1=cbr (constant bit rate)
				cfg->performance_mode = PERF_FULL_COMPUTATION;//0 full computation(HM)
				cfg->reinit_gop_on_scene_change = 0;
				cfg->chroma_qp_offset = 0;
				cfg->wfpp_num_threads = 1;
				cfg->intra_period = 20;
				num_merge_candidates = MERGE_MVP_MAX_NUM_CANDS;
#endif
			hvenc->num_encoder_modules = 2;
			hvenc->max_sublayers = 1;//TLayers en HM
			hvenc->max_layers = 1;

			hvenc->profile = cfg->profile;
			hvenc->intra_period = cfg->intra_period;
			hvenc->gop_size = hvenc->intra_period==1?1:cfg->gop_size;
			hvenc->num_b = hvenc->intra_period==1?0:0;
			hvenc->num_ref_frames = hvenc->gop_size>0?cfg->num_ref_frames:0;	
			hvenc->ctu_height[0] = hvenc->ctu_width[0] = cfg->cu_size;
			hvenc->ctu_height[1] = hvenc->ctu_width[1] = hvenc->ctu_height[2] = hvenc->ctu_width[2] = cfg->cu_size>>1;

			//bitstreams
			hmr_bitstream_alloc(&hvenc->aux_bs, 256);
			hmr_bitstream_alloc(&hvenc->vps_nalu.bs, 256);
			hmr_bitstream_alloc(&hvenc->sps_nalu.bs, 256);
			hmr_bitstream_alloc(&hvenc->pps_nalu.bs, 256);


			hvenc->num_short_term_ref_pic_sets = hvenc->gop_size+1;
			if(hvenc->ref_pic_set_list)
				memset(hvenc->ref_pic_set_list, 0, sizeof(ref_pic_set_t));
			else
				hvenc->ref_pic_set_list = (ref_pic_set_t*)calloc (hvenc->num_short_term_ref_pic_sets, sizeof(ref_pic_set_t));

			for(i=0;i<hvenc->num_short_term_ref_pic_sets-1;i++)
			{
				if(hvenc->intra_period==1)
					hvenc->ref_pic_set_list[i].num_negative_pics = hvenc->ref_pic_set_list[i].num_positive_pics = hvenc->ref_pic_set_list[i].inter_ref_pic_set_prediction_flag = 0;
				else if(hvenc->num_b == 0)
				{
					int j;
					hvenc->ref_pic_set_list[i].num_negative_pics = hvenc->num_ref_frames;
					for(j=0;j<hvenc->ref_pic_set_list[i].num_negative_pics;j++)
					{
						hvenc->ref_pic_set_list[i].delta_poc_s0[j] = -(j+1);//use the last n pictures
						hvenc->ref_pic_set_list[i].used_by_curr_pic_S0_flag[j] = 1;
					}
					hvenc->ref_pic_set_list[i].num_positive_pics = 0;
					hvenc->ref_pic_set_list[i].inter_ref_pic_set_prediction_flag = 0;					
				}
				else
				{
					hvenc->ref_pic_set_list[i].num_positive_pics = hvenc->ref_pic_set_list[i].num_negative_pics = hvenc->num_ref_frames;
					hvenc->ref_pic_set_list[i].inter_ref_pic_set_prediction_flag = 0;
				}
			}


			for(n_enc_mod=0;n_enc_mod<hvenc->num_encoder_modules;n_enc_mod++)
			{
				int ithreads;

				int prev_num_sub_streams, prev_num_ee, prev_num_ec;
				unsigned int min_cu_size, min_cu_size_mask;

				phvenc_mod = (hvenc_t*)calloc(1, sizeof(hvenc_t));
				phvenc_mod->hvenc = hvenc;
				phvenc_mod->index = n_enc_mod;
				hvenc->encoder_module[n_enc_mod] = phvenc_mod;

				if(hvenc->run==1)
				{
					hvenc->run = 0;

					if(phvenc_mod->encoder_thread!=NULL)
					{
						sync_cont_put_filled(phvenc_mod->hvenc->input_hmr_container, NULL);//wake encoder_thread if it is waiting
						JOINT_THREAD(phvenc_mod->encoder_thread);
					}
				}
				phvenc_mod->wfpp_enable = 0;
				phvenc_mod->num_sub_streams = 0;
				phvenc_mod->wfpp_num_threads = 0;

				//bitstreams
				hmr_bitstream_alloc(&phvenc_mod->slice_bs, 0x8000000);
				for(i=0;i<NUM_OUTPUT_NALUS;i++)
					hmr_bitstream_alloc(&phvenc_mod->slice_nalu_list[i].bs, 0x8000000);

				phvenc_mod->avg_dist = 1000;
				phvenc_mod->ctu_width[0] = hvenc->ctu_width[0];
				phvenc_mod->ctu_height[0] = hvenc->ctu_height[0];
				phvenc_mod->ctu_width[1] = phvenc_mod->ctu_width[2] = hvenc->ctu_width[1];
				phvenc_mod->ctu_height[1] = phvenc_mod->ctu_height[2] = hvenc->ctu_height[1];

				phvenc_mod->performance_mode = clip(cfg->performance_mode,0,NUM_PERF_MODES-1);
				phvenc_mod->rd_mode = clip(cfg->rd_mode,0,NUM_RD_MODES-1);
				phvenc_mod->bitrate_mode = clip(cfg->bitrate_mode,0,NUM_BR_MODES-1);
				phvenc_mod->bitrate = cfg->bitrate;
				if(phvenc_mod->bitrate_mode == BR_VBR)
				{
					phvenc_mod->vbv_size = 40*phvenc_mod->bitrate;
					phvenc_mod->vbv_init = .5*phvenc_mod->vbv_size;
					phvenc_mod->qp_min = 15;
				}
				else
				{
					phvenc_mod->vbv_size = cfg->vbv_size;
					phvenc_mod->vbv_init = cfg->vbv_init;
				}
				phvenc_mod->qp_depth = 0;//cfg->qp_depth;//if rc enabled qp_depth == 0

				phvenc_mod->pict_qp = cfg->qp;
				phvenc_mod->chroma_qp_offset = cfg->chroma_qp_offset;

				phvenc_mod->max_cu_size = cfg->cu_size;//MAX_CU_SIZE;
				phvenc_mod->max_cu_size_shift = 0;//MAX_CU_SIZE_SHIFT;
				while(phvenc_mod->max_cu_size>(1<<phvenc_mod->max_cu_size_shift))phvenc_mod->max_cu_size_shift++;

				phvenc_mod->max_pred_partition_depth = (cfg->max_pred_partition_depth>(phvenc_mod->max_cu_size_shift-MIN_TU_SIZE_SHIFT))?(phvenc_mod->max_cu_size_shift-MIN_TU_SIZE_SHIFT):cfg->max_pred_partition_depth;
				if(cfg->motion_estimation_precision==PEL)
					phvenc_mod->motion_estimation_precision=MOTION_PEL_MASK;
				else if(cfg->motion_estimation_precision==HALF_PEL)
					phvenc_mod->motion_estimation_precision=MOTION_HALF_PEL_MASK;
				else
					phvenc_mod->motion_estimation_precision = MOTION_QUARTER_PEL_MASK;

				phvenc_mod->num_merge_mvp_candidates = num_merge_candidates;//MERGE_MVP_MAX_NUM_CANDS;
				if(cfg->width%(phvenc_mod->max_cu_size>>(phvenc_mod->max_pred_partition_depth-1)) || cfg->height%(phvenc_mod->max_cu_size>>(phvenc_mod->max_pred_partition_depth-1)))
				{
					printf("HENC_SETCFG Error- size is not multiple of minimum cu size\r\n");
					goto config_error;
				}
	//			phvenc_mod->max_inter_pred_depth = 0;
	//			phvenc_mod->max_inter_pred_depth = phvenc_mod->max_pred_partition_depth;

				//depth of TU tree 
				phvenc_mod->max_intra_tr_depth = (cfg->max_intra_tr_depth>(phvenc_mod->max_cu_size_shift-MIN_TU_SIZE_SHIFT+1))?(phvenc_mod->max_cu_size_shift-MIN_TU_SIZE_SHIFT+1):cfg->max_intra_tr_depth; 
				phvenc_mod->max_inter_tr_depth = (cfg->max_inter_tr_depth>(phvenc_mod->max_cu_size_shift-MIN_TU_SIZE_SHIFT+1))?(phvenc_mod->max_cu_size_shift-MIN_TU_SIZE_SHIFT+1):cfg->max_inter_tr_depth; 
				phvenc_mod->max_cu_size_shift_chroma = phvenc_mod->max_cu_size_shift-1;
				phvenc_mod->min_tu_size_shift = MIN_TU_SIZE_SHIFT;
				phvenc_mod->max_tu_size_shift = phvenc_mod->max_cu_size_shift<MAX_TU_SIZE_SHIFT?phvenc_mod->max_cu_size_shift:MAX_TU_SIZE_SHIFT;//

				//--------------------------------these are default values for coding cu and tu sizes and herarchy depth-----------------
				phvenc_mod->mincu_mintr_shift_diff = (phvenc_mod->max_cu_size_shift-phvenc_mod->max_pred_partition_depth) - phvenc_mod->min_tu_size_shift;
				phvenc_mod->max_cu_depth = phvenc_mod->max_pred_partition_depth+phvenc_mod->mincu_mintr_shift_diff;
				phvenc_mod->mincu_mintr_shift_diff++;

				min_cu_size = phvenc_mod->max_cu_size  >> ( phvenc_mod->max_cu_depth-phvenc_mod->mincu_mintr_shift_diff);
				min_cu_size_mask = min_cu_size-1;

				phvenc_mod->min_cu_size = min_cu_size;
				aux = phvenc_mod->min_cu_size;
				phvenc_mod->min_cu_size_shift = 0;
				while (aux>1)
				{
					aux>>=1;
					phvenc_mod->min_cu_size_shift++;
				}

				phvenc_mod->num_partitions_in_cu_shift = 4;//each partition is a 4x4 square
				phvenc_mod->num_partitions_in_cu = ((phvenc_mod->max_cu_size*phvenc_mod->max_cu_size)>>phvenc_mod->num_partitions_in_cu_shift);

				phvenc_mod->abs2raster_table = hvenc->abs2raster_table;
				phvenc_mod->raster2abs_table  = hvenc->raster2abs_table;

				phvenc_mod->ang_table = hvenc->ang_table;
				phvenc_mod->inv_ang_table  = hvenc->inv_ang_table;

				phvenc_mod->intra_period = hvenc->intra_period;
				phvenc_mod->last_intra = -cfg->intra_period;
				phvenc_mod->gop_reinit_on_scene_change = cfg->reinit_gop_on_scene_change;
				phvenc_mod->gop_size = hvenc->gop_size;
				phvenc_mod->num_b = hvenc->num_b;
				phvenc_mod->num_ref_frames = hvenc->num_ref_frames;
				//conformance wnd
				phvenc_mod->pad_left = 0;
				if ((cfg->width & min_cu_size_mask) != 0)
					 phvenc_mod->pad_right = (min_cu_size - (cfg->width & min_cu_size_mask));
				else
					phvenc_mod->pad_left = phvenc_mod->pad_right = 0;

				phvenc_mod->pad_top = 0;
				if ((cfg->height & min_cu_size_mask) != 0)
					phvenc_mod->pad_bottom = (min_cu_size - (cfg->height & min_cu_size_mask));
				else
					phvenc_mod->pad_bottom = 0;

				phvenc_mod->conformance_mode = 1;//(0: no conformance, 1:automatic padding
				phvenc_mod->pict_width[0] = cfg->width+phvenc_mod->pad_left + phvenc_mod->pad_right;
				phvenc_mod->pict_height[0] = cfg->height+phvenc_mod->pad_top + phvenc_mod->pad_bottom;
				phvenc_mod->pict_width[1] = phvenc_mod->pict_width[2] = phvenc_mod->pict_width[0]>>1;
				phvenc_mod->pict_height[1] = phvenc_mod->pict_height[2] = phvenc_mod->pict_height[0]>>1;

				//------------------------processing elements------------------------------------
				phvenc_mod->wfpp_enable = cfg->wfpp_enable;
				prev_num_sub_streams = phvenc_mod->num_sub_streams;

				//bitstreams (one per ctu line y wfpp)
				phvenc_mod->num_sub_streams = (phvenc_mod->wfpp_enable==0)?1:(phvenc_mod->pict_height[0] + phvenc_mod->ctu_height[0]-1)/phvenc_mod->ctu_height[0];
				if(prev_num_sub_streams!=phvenc_mod->num_sub_streams)
				{

					//delete previous streams
					for(i=0;i<prev_num_sub_streams;i++)
						hmr_bitstream_free(&phvenc_mod->aux_bs[i]);

					if(phvenc_mod->aux_bs!=NULL)
						free(phvenc_mod->aux_bs);
					if(phvenc_mod->sub_streams_entry_point_list)
						free(phvenc_mod->sub_streams_entry_point_list);

					phvenc_mod->sub_streams_entry_point_list = (uint*)calloc (phvenc_mod->num_sub_streams, sizeof(uint32_t));
					//create new streams
					phvenc_mod->aux_bs = (bitstream_t	*)calloc (phvenc_mod->num_sub_streams, sizeof(bitstream_t));
					for(i=0;i<phvenc_mod->num_sub_streams;i++)
						hmr_bitstream_alloc(&phvenc_mod->aux_bs[i], 0x8000000/phvenc_mod->num_sub_streams);
				}

				//encoding enviroments and rd enviroments (one per thread if wfpp)
				phvenc_mod->wfpp_num_threads = (phvenc_mod->wfpp_enable)?((cfg->wfpp_num_threads<=phvenc_mod->num_sub_streams)?cfg->wfpp_num_threads:phvenc_mod->num_sub_streams):1;
				prev_num_ee = phvenc_mod->num_ee;
				prev_num_ec = phvenc_mod->num_ec;
				phvenc_mod->num_ee = (phvenc_mod->wfpp_enable)?2*phvenc_mod->wfpp_num_threads:1;
				phvenc_mod->num_ec = (phvenc_mod->wfpp_enable)?phvenc_mod->wfpp_num_threads:1;

				if(prev_num_ee != phvenc_mod->num_ee)
				{
					//delete previous enviroments
					for(i=0;i<prev_num_ee;i++)
					{
						free(phvenc_mod->ee_list[i]->e_ctx);phvenc_mod->ee_list[i]->e_ctx=NULL;
						free(phvenc_mod->ee_list[i]->contexts);phvenc_mod->ee_list[i]->contexts=NULL;
						free(phvenc_mod->ee_list[i]->b_ctx);phvenc_mod->ee_list[i]->b_ctx=NULL;
						phvenc_mod->ee_list[i]->type = EE_INVALID;
					}
					if(phvenc_mod->ee_list!=NULL)
						free(phvenc_mod->ee_list);

					//create new enviroments
					phvenc_mod->ee_list = (enc_env_t **)calloc (phvenc_mod->num_ee, sizeof(enc_env_t*));
					for(i=0;i<phvenc_mod->num_ee;i++)
					{
						phvenc_mod->ee_list[i] = (enc_env_t *)calloc (1, sizeof(enc_env_t));
						phvenc_mod->ee_list[i]->e_ctx = (entropy_model_t*)calloc(1, sizeof(entropy_model_t)); 	
						phvenc_mod->ee_list[i]->contexts = (context_model_t*)calloc(NUM_CTXs, sizeof(context_model_t));
						phvenc_mod->ee_list[i]->b_ctx = (binary_model_t*)calloc(1, sizeof(binary_model_t));//calloc(NUM_CTXs, sizeof(binary_model_t));
						phvenc_mod->ee_list[i]->type = EE_ENCODER;
						ee_init_contexts(phvenc_mod->ee_list[i]);
						bm_map_funcs(phvenc_mod->ee_list[i]);
					}
				}

				if(prev_num_ec != phvenc_mod->num_ec)
				{
					//delete previous rd counters
					for(i=0;i<prev_num_ec;i++)
					{
						free(phvenc_mod->ec_list[i].e_ctx);phvenc_mod->ec_list[i].e_ctx=NULL;
						free(phvenc_mod->ec_list[i].contexts);phvenc_mod->ec_list[i].contexts=NULL;
						free(phvenc_mod->ec_list[i].b_ctx);phvenc_mod->ec_list[i].b_ctx=NULL;
						phvenc_mod->ec_list[i].type = EE_INVALID;
					}
					if(phvenc_mod->ec_list!=NULL)
						free(phvenc_mod->ec_list);

					//create new rd counters
					phvenc_mod->ec_list = (enc_env_t *)calloc (phvenc_mod->num_ec, sizeof(enc_env_t));

					for(i=0;i<phvenc_mod->num_ec;i++)
					{
						phvenc_mod->ec_list[i].e_ctx = (entropy_model_t*)calloc(1, sizeof(entropy_model_t)); 	
						phvenc_mod->ec_list[i].contexts = (context_model_t*)calloc(NUM_CTXs, sizeof(context_model_t));
						phvenc_mod->ec_list[i].b_ctx = (binary_model_t*)calloc(1, sizeof(binary_model_t));//calloc(NUM_CTXs, sizeof(binary_model_t));
						phvenc_mod->ec_list[i].type = EE_COUNTER;
						ee_init_contexts(&phvenc_mod->ec_list[i]);
						bm_map_funcs(&phvenc_mod->ec_list[i]);
					}
				}

				phvenc_mod->frame_rate = cfg->frame_rate;
	//			phvenc_mod->pic_interlaced = 0;
	//			phvenc_mod->mb_interlaced = 0;

				phvenc_mod->bit_depth = hvenc->bit_depth;

				phvenc_mod->pict_width_in_ctu = (phvenc_mod->pict_width[0]>>phvenc_mod->max_cu_size_shift) + ((phvenc_mod->pict_width[0]%phvenc_mod->max_cu_size)!=0);
				phvenc_mod->pict_height_in_ctu = (phvenc_mod->pict_height[0]>>phvenc_mod->max_cu_size_shift) + ((phvenc_mod->pict_height[0]%phvenc_mod->max_cu_size)!=0);

				phvenc_mod->pict_total_ctu = phvenc_mod->pict_width_in_ctu*phvenc_mod->pict_height_in_ctu;

				if(phvenc_mod->ctu_info!=NULL)
				{
					for(i=0;i<phvenc_mod->pict_total_ctu;i++)
					{
						HMR_FREE(phvenc_mod->ctu_info[i].cbf[Y_COMP]);
						//intra mode
						HMR_FREE(phvenc_mod->ctu_info[i].intra_mode[Y_COMP]);
						HMR_FREE(phvenc_mod->ctu_info[i].inter_mode);
						//tr_idx, pred_depth, part_size_type, pred_mode
						HMR_FREE(phvenc_mod->ctu_info[i].tr_idx);
						HMR_FREE(phvenc_mod->ctu_info[i].pred_depth);
						HMR_FREE(phvenc_mod->ctu_info[i].part_size_type);
						HMR_FREE(phvenc_mod->ctu_info[i].pred_mode);

						HMR_FREE(phvenc_mod->ctu_info[i].skipped);
						HMR_FREE(phvenc_mod->ctu_info[i].merge);
						HMR_FREE(phvenc_mod->ctu_info[i].merge_idx);
						//inter
						HMR_FREE(phvenc_mod->ctu_info[i].mv_ref[REF_PIC_LIST_0]);
						HMR_FREE(phvenc_mod->ctu_info[i].mv_ref_idx[REF_PIC_LIST_0]);
						HMR_FREE(phvenc_mod->ctu_info[i].mv_diff[REF_PIC_LIST_0]);
						HMR_FREE(phvenc_mod->ctu_info[i].mv_diff_ref_idx[REF_PIC_LIST_0]);
	//					free(phvenc_mod->ctu_info[i].mv_ref1);
					
	//					free(phvenc_mod->ctu_info[i].ref_idx1);
						free(phvenc_mod->ctu_info[i].qp);
					}
					free(phvenc_mod->ctu_info);
				}
				phvenc_mod->ctu_info = (ctu_info_t*)calloc (phvenc_mod->pict_total_ctu, sizeof(ctu_info_t));

				for(i=0;i<phvenc_mod->pict_total_ctu;i++)
				{
					//------- ctu encoding info -------
					//cbf
					phvenc_mod->ctu_info[i].cbf[Y_COMP] = (uint8_t*)calloc (NUM_PICT_COMPONENTS*MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].cbf[U_COMP] = phvenc_mod->ctu_info[i].cbf[Y_COMP]+MAX_NUM_PARTITIONS;
					phvenc_mod->ctu_info[i].cbf[V_COMP] = phvenc_mod->ctu_info[i].cbf[U_COMP]+MAX_NUM_PARTITIONS;
					//intra mode
					phvenc_mod->ctu_info[i].intra_mode[Y_COMP] = (uint8_t*)calloc (2*MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].intra_mode[CHR_COMP] = phvenc_mod->ctu_info[i].intra_mode[Y_COMP]+MAX_NUM_PARTITIONS;
					//inter mode
					phvenc_mod->ctu_info[i].inter_mode = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					//tr_idx, pred_depth, part_size_type, pred_mode, skipped
					phvenc_mod->ctu_info[i].tr_idx = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].pred_depth = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].part_size_type = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].pred_mode = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].skipped = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].merge = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].merge_idx = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].qp = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));

					//inter
					phvenc_mod->ctu_info[i].mv_ref[REF_PIC_LIST_0] = (motion_vector_t*)calloc (2*MAX_NUM_PARTITIONS, sizeof(motion_vector_t));
					phvenc_mod->ctu_info[i].mv_ref[REF_PIC_LIST_1] = phvenc_mod->ctu_info[i].mv_ref[REF_PIC_LIST_0]+MAX_NUM_PARTITIONS;
					phvenc_mod->ctu_info[i].mv_ref_idx[REF_PIC_LIST_0] = (int8_t*)calloc (2*MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].mv_ref_idx[REF_PIC_LIST_1] = phvenc_mod->ctu_info[i].mv_ref_idx[REF_PIC_LIST_0]+MAX_NUM_PARTITIONS;

					phvenc_mod->ctu_info[i].mv_diff[REF_PIC_LIST_0] = (motion_vector_t*)calloc (2*MAX_NUM_PARTITIONS, sizeof(motion_vector_t));
					phvenc_mod->ctu_info[i].mv_diff[REF_PIC_LIST_1] = phvenc_mod->ctu_info[i].mv_diff[REF_PIC_LIST_0]+MAX_NUM_PARTITIONS;
					phvenc_mod->ctu_info[i].mv_diff_ref_idx[REF_PIC_LIST_0] = (uint8_t*)calloc (2*MAX_NUM_PARTITIONS, sizeof(uint8_t));
					phvenc_mod->ctu_info[i].mv_diff_ref_idx[REF_PIC_LIST_1] = phvenc_mod->ctu_info[i].mv_diff_ref_idx[REF_PIC_LIST_0]+MAX_NUM_PARTITIONS;
				}


				phvenc_mod->max_sublayers = hvenc->max_sublayers;//TLayers en HM
				phvenc_mod->max_layers = hvenc->max_layers;

				if(phvenc_mod->deblock_filter_sem!=NULL)
					SEM_DESTROY(phvenc_mod->deblock_filter_sem);
				SEM_INIT(phvenc_mod->deblock_filter_semaphore, 0,1000);
				SEM_COPY(phvenc_mod->deblock_filter_sem, phvenc_mod->deblock_filter_semaphore);
				for(ithreads=0;ithreads<phvenc_mod->wfpp_num_threads;ithreads++)
				{
					int depth_aux;
					int j;
					henc_thread_t* henc_th = (henc_thread_t*)calloc(1, sizeof(henc_thread_t));
					int filter_buff_width, filter_buff_height;

	//				henc_th = &phvenc_mod->_thread;
					phvenc_mod->thread[ithreads] = henc_th;
					henc_th->ed = phvenc_mod;
					henc_th->index = ithreads;
					henc_th->funcs = &hvenc->funcs;
					henc_th->wfpp_enable = phvenc_mod->wfpp_enable;
					henc_th->wfpp_num_threads = phvenc_mod->wfpp_num_threads;

					//alloc processing windows (processing buffers attached to windows will be allocated depending on the resolution)
					henc_th->wfpp_num_threads = phvenc_mod->wfpp_num_threads;

					SEM_COPY(henc_th->deblock_filter_sem, phvenc_mod->deblock_filter_semaphore);

					if(henc_th->synchro_signal!=NULL)
						SEM_DESTROY(henc_th->synchro_signal);

					SEM_INIT(henc_th->synchro_sem, 0,1000);
					SEM_COPY(henc_th->synchro_signal, henc_th->synchro_sem);
				
					henc_th->vps = &phvenc_mod->hvenc->vps;
					henc_th->sps = &phvenc_mod->hvenc->sps;
					henc_th->pps = &phvenc_mod->hvenc->pps;

					memcpy(henc_th->pict_width, phvenc_mod->pict_width, sizeof(henc_th->pict_width));
					memcpy(henc_th->pict_height, phvenc_mod->pict_height, sizeof(henc_th->pict_height));
					henc_th->pict_width_in_ctu = phvenc_mod->pict_width_in_ctu; 
					henc_th->pict_height_in_ctu = phvenc_mod->pict_height_in_ctu;			
					henc_th->pict_total_ctu = phvenc_mod->pict_total_ctu;

					memcpy(henc_th->ctu_width, phvenc_mod->ctu_width, sizeof(henc_th->ctu_width));
					memcpy(henc_th->ctu_height, phvenc_mod->ctu_height, sizeof(henc_th->ctu_height));
//					henc_th->ctu_group_size = phvenc_mod->ctu_group_size;

					henc_th->max_cu_size = phvenc_mod->max_cu_size;
					henc_th->max_cu_size_shift = phvenc_mod->max_cu_size_shift;//log2 del tama�o del CU maximo
					henc_th->max_cu_size_shift_chroma = phvenc_mod->max_cu_size_shift_chroma;//log2 del tama�o del CU maximo
					henc_th->max_intra_tr_depth = phvenc_mod->max_intra_tr_depth;
					henc_th->max_inter_tr_depth = phvenc_mod->max_inter_tr_depth;
					henc_th->max_pred_partition_depth = phvenc_mod->max_pred_partition_depth;//max depth for prediction
					henc_th->motion_estimation_precision = phvenc_mod->motion_estimation_precision;
	//				henc_th->max_inter_pred_depth = phvenc_mod->max_inter_pred_depth;//max depth for prediction

					henc_th->num_partitions_in_cu = phvenc_mod->num_partitions_in_cu;
					henc_th->num_partitions_in_cu_shift = phvenc_mod->num_partitions_in_cu_shift;
					henc_th->mincu_mintr_shift_diff = phvenc_mod->mincu_mintr_shift_diff;
					henc_th->max_cu_depth = phvenc_mod->max_cu_depth;
					henc_th->min_cu_size = phvenc_mod->min_cu_size;
					henc_th->min_cu_size_shift = phvenc_mod->min_cu_size_shift;
					henc_th->min_tu_size_shift = phvenc_mod->min_tu_size_shift;
					henc_th->max_tu_size_shift = phvenc_mod->max_tu_size_shift;

					henc_th->profile = hvenc->profile;
					henc_th->bit_depth = phvenc_mod->bit_depth;
					henc_th->performance_mode = phvenc_mod->performance_mode;
					henc_th->rd_mode = phvenc_mod->rd_mode;

					henc_th->partition_depth_start = phvenc_mod->partition_depth_start;
					aux = 1;//el CTU
					for(depth_aux=1;depth_aux<MAX_PARTITION_DEPTH;depth_aux++)//5 = max partition depth 
					{
						aux += (4<<(2*(depth_aux-1)));
					}

					henc_th->partition_info = (cu_partition_info_t*)calloc (aux, sizeof(cu_partition_info_t));
					//init partition list info
					init_partition_info(henc_th, henc_th->partition_info);
					//----------------------------current thread processing buffers allocation	-----
					//alloc processing windows and buffers
					henc_th->adi_size = 2*henc_th->ctu_height[0]+2*henc_th->ctu_width[0]+1;//vecinos de la columna izq y fila superior + la esquina
					henc_th->adi_pred_buff = (short*)hmr_aligned_alloc (henc_th->adi_size, sizeof(int16_t));
					henc_th->adi_filtered_pred_buff = (short*)hmr_aligned_alloc (henc_th->adi_size, sizeof(int16_t));
					henc_th->top_pred_buff = (short*)hmr_aligned_alloc (henc_th->adi_size, sizeof(int16_t));
					henc_th->left_pred_buff = (short*)hmr_aligned_alloc (henc_th->adi_size, sizeof(int16_t));
					henc_th->bottom_pred_buff = (short*)hmr_aligned_alloc (henc_th->adi_size, sizeof(int16_t));
					henc_th->right_pred_buff = (short*)hmr_aligned_alloc (henc_th->adi_size, sizeof(int16_t));


					wnd_realloc(&henc_th->curr_mbs_wnd, (henc_th->ctu_width[0]), henc_th->ctu_height[0], 0, 0, sizeof(uint8_t));

					henc_th->pred_aux_buff_size = MAX_CU_SIZE*MAX_CU_SIZE;//tama�o del buffer auxiliar
					henc_th->pred_aux_buff = (short*) hmr_aligned_alloc (henc_th->pred_aux_buff_size, sizeof(short));

					wnd_realloc(&henc_th->prediction_wnd, (henc_th->ctu_width[0]), henc_th->ctu_height[0], 0, 0, sizeof(int16_t));
					wnd_realloc(&henc_th->residual_wnd, (henc_th->ctu_width[0]), henc_th->ctu_height[0], 0, 0, sizeof(int16_t));
					wnd_realloc(&henc_th->residual_dec_wnd, (henc_th->ctu_width[0]), henc_th->ctu_height[0], 0, 0, sizeof(int16_t));

					henc_th->aux_buff = (short*) hmr_aligned_alloc (MAX_CU_SIZE*MAX_CU_SIZE, sizeof(int));

					wnd_realloc(&henc_th->itransform_iquant_wnd, (henc_th->ctu_width[0]), henc_th->ctu_height[0], 0, 0, sizeof(int16_t));

					for(i=0;i<NUM_DECODED_WNDS;i++)
					{
						wnd_realloc(&henc_th->decoded_mbs_wnd_[i], 2*(henc_th->ctu_width[0]), henc_th->ctu_height[0]*2, 1, 1, sizeof(int16_t));
						henc_th->decoded_mbs_wnd[i] = &henc_th->decoded_mbs_wnd_[i];//use pointers to exchange windows
					}

					for(i=0;i<NUM_QUANT_WNDS;i++)
					{
						wnd_realloc(&henc_th->transform_quant_wnd_[i], (henc_th->ctu_width[0]), henc_th->ctu_height[0], 0, 0, sizeof(int16_t));		
						henc_th->transform_quant_wnd[i] = &henc_th->transform_quant_wnd_[i];//use pointers to exchange windows
					}


					//deblocking filter
					henc_th->deblock_partition_info = (cu_partition_info_t*)calloc (aux, sizeof(cu_partition_info_t));
					init_partition_info(henc_th, henc_th->deblock_partition_info);
					henc_th->deblock_filter_strength_bs[EDGE_VER] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->deblock_filter_strength_bs[EDGE_HOR] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->deblock_edge_filter[EDGE_VER] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->deblock_edge_filter[EDGE_HOR] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));

					//quarter pel interpolation
					filter_buff_width = MAX_CU_SIZE	+ 16;
					filter_buff_height = MAX_CU_SIZE + 2;//MAX_CU_SIZE + 1; - modified for the chroma
					for(j=0;j<4;j++)
					{
						wnd_realloc(&henc_th->filtered_block_temp_wnd[j], filter_buff_width, filter_buff_height+8, 0, 0, sizeof(int16_t));//filter_buff_height+7 - modified for chroma
						for(i=0;i<4;i++)
						{
							wnd_realloc(&henc_th->filtered_block_wnd[j][i], filter_buff_width, filter_buff_height, 0, 0, sizeof(int16_t));				
						}

					}

					henc_th->cabac_aux_buff_size = MAX_CU_SIZE*MAX_CU_SIZE;//MAX_TU_SIZE_SHIFT*MAX_TU_SIZE_SHIFT;//tama�o del buffer auxiliar
					henc_th->cabac_aux_buff = (unsigned char*) hmr_aligned_alloc (henc_th->cabac_aux_buff_size, sizeof(unsigned char));

					//alloc buffers to gather and consolidate information
					for(i=0;i<NUM_CBF_BUFFS;i++)
					{
						henc_th->cbf_buffs[Y_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->cbf_buffs[U_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->cbf_buffs[V_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->intra_mode_buffs[Y_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->intra_mode_buffs[U_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->intra_mode_buffs[V_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->tr_idx_buffs[i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));

						//inter
	/*					henc_th->mv_ref0[Y_COMP][i] = (motion_vector_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(motion_vector_t));
						henc_th->mv_ref0[U_COMP][i] = (motion_vector_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(motion_vector_t));
						henc_th->mv_ref0[V_COMP][i] = (motion_vector_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(motion_vector_t));

						henc_th->mv_ref1[Y_COMP][i] = (motion_vector_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(motion_vector_t));
						henc_th->mv_ref1[U_COMP][i] = (motion_vector_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(motion_vector_t));
						henc_th->mv_ref1[V_COMP][i] = (motion_vector_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(motion_vector_t));

						henc_th->ref_idx0[Y_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->ref_idx0[V_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->ref_idx0[U_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));

						henc_th->ref_idx1[Y_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->ref_idx1[V_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
						henc_th->ref_idx1[U_COMP][i] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
	*/				}

					henc_th->cbf_buffs_chroma[U_COMP] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->cbf_buffs_chroma[V_COMP] = (uint8_t*) hmr_aligned_alloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));

					henc_th->ctu_rd = (ctu_info_t*)calloc (1, sizeof(ctu_info_t));
					henc_th->ctu_rd->part_size_type = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->ctu_rd->pred_mode = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->ctu_rd->skipped = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->ctu_rd->merge = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->ctu_rd->merge_idx = (uint8_t*)calloc (MAX_NUM_PARTITIONS, sizeof(uint8_t));
					henc_th->ee = phvenc_mod->ee_list[2*henc_th->index];
					henc_th->ec = &phvenc_mod->ec_list[henc_th->index];
				}

				//exterchange wait and signal semaphores between sucessive threads
				if(phvenc_mod->wfpp_num_threads==1)
				{
					SEM_COPY(phvenc_mod->thread[0]->synchro_wait, phvenc_mod->thread[0]->synchro_sem);
				}
				else
				{
					for(i=0;i<phvenc_mod->wfpp_num_threads;i++)
					{
						SEM_COPY(phvenc_mod->thread[i]->synchro_wait, phvenc_mod->thread[(i+phvenc_mod->wfpp_num_threads-1)%phvenc_mod->wfpp_num_threads]->synchro_sem);
					}			
				}

				//reference picture lists
				phvenc_mod->num_ref_lists = 2;
				phvenc_mod->num_refs_idx_active_list[REF_PIC_LIST_0] = phvenc_mod->intra_period==1?4:1;//this will need to be a consistent decission taken depending on the configuration
				phvenc_mod->num_refs_idx_active_list[REF_PIC_LIST_1] = phvenc_mod->intra_period==1?4:1;

				hmr_rc_init(phvenc_mod);
			}// end for(n_enc_mod=0;n_enc_mod<hvenc->num_encoder_modules;n_enc_mod++)


			hvenc->pict_width[0] = phvenc_mod->pict_width[0];
			hvenc->pict_height[0] = phvenc_mod->pict_height[0];
			hvenc->pict_width[1] = hvenc->pict_width[2] = phvenc_mod->pict_width[1];
			hvenc->pict_height[1] = hvenc->pict_height[2] = phvenc_mod->pict_height[1];

			sync_cont_reset(hvenc->input_hmr_container);
			for(i=0;i<NUM_INPUT_FRAMES;i++)
			{
				wnd_realloc(&hvenc->input_frames[i].img, hvenc->pict_width[0], hvenc->pict_height[0], 0, 0, sizeof(uint8_t));
				sync_cont_put_empty(hvenc->input_hmr_container, &hvenc->input_frames[i]);
			}

			cont_reset(hvenc->output_hmr_container);

			cont_reset(hvenc->cont_empty_reference_wnds);
			for(i=0;i<2*MAX_NUM_REF;i++)
			{
				wnd_realloc(&hvenc->ref_wnds[i].img, hvenc->pict_width[0], hvenc->pict_height[0], hvenc->ctu_width[Y_COMP]+16, hvenc->ctu_height[Y_COMP]+16, sizeof(int16_t));
				cont_put(hvenc->cont_empty_reference_wnds, &hvenc->ref_wnds[i]);
			}

			hvenc->ptl.generalPTL.profileIdc = hvenc->profile;

			hvenc->ptl.generalPTL.profileCompatibilityFlag[hvenc->profile] = 1;
			if (hvenc->profile == PROFILE_MAIN)// A PROFILE_MAIN10 can decode PROFILE_MAIN
				hvenc->ptl.generalPTL.profileCompatibilityFlag[PROFILE_MAIN10] = 1;

			if (hvenc->profile == PROFILE_MAIN10 && phvenc_mod->bit_depth == 8)// PROFILE_MAIN10 with 8 bits = PROFILE_MAIN
				hvenc->ptl.generalPTL.profileCompatibilityFlag[PROFILE_MAIN] = 1;

			//profiles and levels
			memset(hvenc->ptl.subLayerProfilePresentFlag, 0, sizeof(hvenc->ptl.subLayerProfilePresentFlag));
			memset(hvenc->ptl.subLayerLevelPresentFlag,   0, sizeof(hvenc->ptl.subLayerLevelPresentFlag  ));

			//----------------- start vps ------------------
			hvenc->vps.video_parameter_set_id = 0;
			hvenc->vps.temporal_id_nesting_flag = (hvenc->max_sublayers == 1);
			hvenc->vps.ptl = &hvenc->ptl;

			hvenc->vps.sub_layer_ordering_info_present_flag = 1;
			for(i = 0; i < MAX_TLAYER; i++)
			{
				hvenc->vps.max_num_reorder_pics[i] = 0;
				hvenc->vps.max_dec_pic_buffering[i] = (hvenc->intra_period==1)?1:hvenc->num_ref_frames+1;//m_maxDecPicBuffering[m_GOPList[i].m_temporalId] = m_GOPList[i].m_numRefPics + 1;
				hvenc->vps.max_latency_increase[i] = 0;
			}
	
			hvenc->vps.timing_info_present_flag = 0;

			//----------------- end vps ------------------

			//----------------- start sps ------------------
			hvenc->sps.video_parameter_set_id = 0;
			
			hvenc->sps.ptl = &hvenc->ptl;

			hvenc->sps.seq_parameter_set_id = 0;
			hvenc->sps.chroma_format_idc = CHROMA420;

			hvenc->sps.pic_width_in_luma_samples = hvenc->pict_width[0];
			hvenc->sps.pic_height_in_luma_samples = hvenc->pict_height[0];
			hvenc->sps.conformance_window_flag = phvenc_mod->conformance_mode!=0;
			hvenc->sps.conf_win_left_offset = phvenc_mod->pad_left/pad_unit_x[hvenc->sps.chroma_format_idc];
			hvenc->sps.conf_win_right_offset = phvenc_mod->pad_right/pad_unit_x[hvenc->sps.chroma_format_idc];
			hvenc->sps.conf_win_top_offset = phvenc_mod->pad_top/pad_unit_y[hvenc->sps.chroma_format_idc];
			hvenc->sps.conf_win_bottom_offset = phvenc_mod->pad_bottom/pad_unit_y[hvenc->sps.chroma_format_idc];
			hvenc->sps.bit_depth_luma_minus8 = hvenc->sps.bit_depth_chroma_minus8 = hvenc->bit_depth-8;
			hvenc->sps.pcm_enabled_flag = 0;
			hvenc->sps.log2_max_pic_order_cnt_lsb_minus4 = 0;//4; //bits for poc=8 - must be bigger than idr period
			for(i = 0; i < MAX_TLAYER; i++)
			{
				hvenc->sps.max_num_reorder_pics[i] = 0;
				hvenc->sps.max_dec_pic_buffering[i] = (hvenc->intra_period==1)?1:hvenc->num_ref_frames+1;//m_maxDecPicBuffering[m_GOPList[i].m_temporalId] = m_GOPList[i].m_numRefPics + 1;
				hvenc->sps.max_latency_increase[i] = 0;
			}

			hvenc->sps.restricted_ref_pic_lists_flag = 1;
			hvenc->sps.lists_modification_present_flag = 0;
			hvenc->sps.log2_min_coding_block_size_minus3 = phvenc_mod->min_cu_size_shift-3;
			hvenc->sps.log2_diff_max_min_coding_block_size = phvenc_mod->max_cu_depth-phvenc_mod->mincu_mintr_shift_diff;//hvenc->num_partitions_in_cu_shift;//-hvenc->min_cu_size_shift;
//			hvenc->sps.log2_diff_max_min_coding_block_size = hvenc->max_cu_size_shift-hvenc->min_cu_size_shift;
			hvenc->sps.log2_min_transform_block_size_minus2 = phvenc_mod->min_tu_size_shift-2;
			hvenc->sps.log2_diff_max_min_transform_block_size = phvenc_mod->max_tu_size_shift - phvenc_mod->min_tu_size_shift;
			hvenc->sps.max_transform_hierarchy_depth_inter = phvenc_mod->max_inter_tr_depth-1;
			hvenc->sps.max_transform_hierarchy_depth_intra = phvenc_mod->max_intra_tr_depth-1;
			hvenc->sps.scaling_list_enabled_flag = 1;
			hvenc->sps.scaling_list_data_present_flag = 0;

			hvenc->sps.amp_enabled_flag = 0;
			hvenc->sps.sample_adaptive_offset_enabled_flag = 0;//cfg->UseSAO;
			hvenc->sps.temporal_id_nesting_flag = (hvenc->max_sublayers == 1);
			hvenc->num_long_term_ref_pic_sets = 0;
			hvenc->sps.temporal_mvp_enable_flag = 0;
			hvenc->sps.strong_intra_smooth_enabled_flag = 1;
			hvenc->sps.vui_parameters_present_flag = 0;
			//----------------- end sps ------------------
			
			//----------------- start pps ------------------
			hvenc->pps.pic_parameter_set_id = 0;
			hvenc->pps.seq_parameter_set_id = hvenc->sps.seq_parameter_set_id;
			hvenc->pps.dependent_slice_enabled_flag = 0;
			hvenc->pps.output_flag_present_flag = 0;
			hvenc->pps.num_extra_slice_header_bits = 0;
			hvenc->pps.sign_data_hiding_flag = cfg->sign_hiding;
			hvenc->pps.cabac_init_present_flag = 0;//1

			hvenc->pps.num_ref_idx_l0_default_active_minus1 = 0;
			hvenc->pps.num_ref_idx_l1_default_active_minus1 = 0;
			
#ifdef COMPUTE_AS_HM
			hvenc->pps.pic_init_qp_minus26 = 0;
#else
			hvenc->pps.pic_init_qp_minus26 = phvenc_mod->pict_qp - 26;
#endif
			hvenc->pps.constrained_intra_pred_flag = 0;
			hvenc->pps.transform_skip_enabled_flag = 0;
			hvenc->pps.cu_qp_delta_enabled_flag = (phvenc_mod->bitrate_mode==BR_FIXED_QP)?0:1;
			hvenc->pps.diff_cu_qp_delta_depth = phvenc_mod->qp_depth;

			hvenc->pps.cb_qp_offset = phvenc_mod->chroma_qp_offset ;
			hvenc->pps.cr_qp_offset = phvenc_mod->chroma_qp_offset ;

			hvenc->pps.slice_chroma_qp_offsets_present_flag = 0;
			hvenc->pps.weighted_pred_flag = 0;
			hvenc->pps.weighted_bipred_flag = 0;
			hvenc->pps.output_flag_present_flag = 0;
			hvenc->pps.transquant_bypass_enable_flag = 0;
			hvenc->pps.tiles_enabled_flag = 0;	
			hvenc->pps.entropy_coded_sync_enabled_flag = phvenc_mod->wfpp_enable;

/*			if(hvenc->pps.tiles_enabled_flag)
				//....................		
*/
			hvenc->pps.loop_filter_across_slices_enabled_flag = 1;
			hvenc->pps.deblocking_filter_control_present_flag = 0;
			hvenc->pps.beta_offset_div2 = 0;
			hvenc->pps.tc_offset_div2 = 0;
/*			if(hvenc->pps.deblocking_filter_control_present_flag)
				//.......................
*/			hvenc->pps.pps_scaling_list_data_present_flag = 0;
//			if(hvenc->pps.pps_scaling_list_data_present_flag)
				//.......................
			hvenc->pps.lists_modification_present_flag = 0;
			hvenc->pps.log2_parallel_merge_level_minus2 = 0;
			hvenc->pps.num_extra_slice_header_bits = 0;
			hvenc->pps.slice_header_extension_present_flag = 0;
			//----------------- end pps ------------------
			hvenc->run = 1;
			CREATE_THREAD(hvenc->encoder_module[0]->encoder_thread, encoder_thread, hvenc->encoder_module[0]);//falta
			CREATE_THREAD(hvenc->encoder_module[1]->encoder_thread, encoder_thread, hvenc->encoder_module[1]);//falta

		}	
		break;

	}   
	if(err)     
    	return (FALSE);
	else
		return (TRUE);
config_error:
	return (FALSE);
}

int get_nal_unit_type(hvenc_t* ed, slice_t *curr_slice, int curr_poc)
{
	if(curr_slice->slice_type == I_SLICE)//(curr_poc == 0)
	{
		return NALU_CODED_SLICE_IDR_W_RADL;
	}

	return NALU_CODED_SLICE_TRAIL_R;
}


void reference_picture_border_padding_ctu(wnd_t *wnd, ctu_info_t* ctu)
{
	int component, j, i;
	int pad_left = FALSE, pad_right = FALSE, pad_top = FALSE, pad_bottom = FALSE;
	int frame_width = WND_WIDTH_2D(*wnd, Y_COMP);
	int frame_height = WND_HEIGHT_2D(*wnd, Y_COMP);

	if(ctu->x[Y_COMP] == 0)
		pad_left = TRUE;

	if(ctu->y[Y_COMP] == 0)
		pad_top = TRUE;

	if(ctu->x[Y_COMP] + ctu->size >= frame_width)
		pad_right = TRUE;

	if(ctu->y[Y_COMP] + ctu->size >= frame_height)
		pad_bottom = TRUE;

	for(component = Y_COMP; component <= V_COMP; component++)
	{
		int stride = WND_STRIDE_2D(*wnd, component);
		int padding_x = wnd->data_padding_x[component];
		int padding_y = wnd->data_padding_y[component];
		int cu_width = (component==Y_COMP)?ctu->size:(ctu->size>>1);
		int cu_height = cu_width;
		int16_t *ptr = WND_DATA_PTR(int16_t *, *wnd, component) + ctu->y[component]*stride + ctu->x[component];

		frame_width = WND_WIDTH_2D(*wnd, component);
		frame_height = WND_HEIGHT_2D(*wnd, component);

		if(pad_left)
		{
			int16_t *ptr_left = ptr-padding_x;
			int16_t *ptr_orig = ptr;
			for(j=0;j<cu_height;j++)
			{
				for(i=0;i<padding_x;i++)
				{
					ptr_left[i] = ptr_orig[0];
				}
	//			memset(ptr_left, ptr[0], padding_x);
				ptr_left+=stride;
				ptr_orig+=stride;
			}
		}

		if(pad_right)
		{
			int16_t *ptr_right = ptr + cu_width;
			int16_t *ptr_orig = ptr + cu_width - 1;
			for(j=0;j<cu_height;j++)
			{
				for(i=0;i<padding_x;i++)
				{
					ptr_right[i] = ptr_orig[0];
				}
	//			memset(ptr_left, ptr[0], padding_x);
				ptr_right+=stride;
				ptr_orig+=stride;
			}
		}


		if(pad_top)
		{
			int16_t *ptr_top = ptr - padding_y*stride;
			int16_t *ptr_orig = ptr;
			int padding_size = cu_width;
			if(pad_left)
			{
				ptr_top -= padding_x; 
				ptr_orig -= padding_x; 
				padding_size += padding_x;// += cu_width;
			}
			if(pad_right)
			{
//				ptr_top -= padding_x; 
//				ptr_orig -= padding_x; 
				padding_size += padding_x;
			}
			for(j=0;j<padding_y;j++)
			{
				memcpy(ptr_top, ptr_orig, padding_size*sizeof(ptr_orig[0]));
				ptr_top+=stride;
			}
		}

		if(pad_bottom)
		{
			int padding_height_init = (frame_height-ctu->y[component])<cu_height?(frame_height-ctu->y[component]):cu_height;
			int16_t *ptr_bottom = ptr + (padding_height_init)*stride;
			int16_t *ptr_orig = ptr + (padding_height_init-1)*stride;
			int padding_size = cu_width;

			if(component == Y_COMP)
			{
				int iiiii=0;
			}

			if(pad_left)
			{
				ptr_bottom -= padding_x; 
				ptr_orig -= padding_x; 
				padding_size += padding_x;// += cu_width;
			}
			if(pad_right)
			{
//				ptr_bottom -= padding_x; 
//				ptr_orig -= padding_x; 
				padding_size += padding_x;
			}
			for(j=0;j<padding_y;j++)
			{
				memcpy(ptr_bottom, ptr_orig, padding_size*sizeof(ptr_orig[0]));
				ptr_bottom += stride;
			}
		}
	}
}


void reference_picture_border_padding(wnd_t *wnd)
{
	int   component, j, i;

	for(component = Y_COMP; component <= V_COMP; component++)
	{
		int stride = WND_STRIDE_2D(*wnd, component);
		int padding_x = wnd->data_padding_x[component];
		int padding_y = wnd->data_padding_y[component];
		int data_width = wnd->data_width[component];
		int data_height = wnd->data_height[component];
		int16_t *ptr = WND_DATA_PTR(int16_t *, *wnd, component);
		int16_t *ptr_left = ptr-padding_x;
		int16_t *ptr_right = ptr+data_width;
		int16_t *ptr_top;// = ptr_left-stride;
		int16_t *ptr_bottom;// = ptr_left+(data_height)*stride;
		for(j=0;j<data_height;j++)
		{
			for(i=0;i<padding_x;i++)
			{
				if(ptr_left[i] != ptr[0])//ptr_left[i] = ptr[0];
					printf("left padding error: j:%d, i:%d", j, i);
				if(ptr_right[i]  != ptr[data_width-1])//ptr_right[i]  = ptr[data_width-1];
					printf("right padding error: j:%d, i:%d", j, i);
			}
//			memset(ptr_left, ptr[0], padding_x);
//			memset(ptr_right, ptr[data_width-1], padding_x);
			ptr_left+=stride;
			ptr_right+=stride;
			ptr+=stride;
		}

		ptr = WND_DATA_PTR(int16_t *, *wnd, component);
		ptr += (data_height-1)*stride-padding_x;
		ptr_bottom = ptr+stride;
		
		for(j=0;j<padding_y;j++)
		{
			for(i=0;i<stride;i++)
			{
				if(ptr_bottom[i] != ptr[i])
					printf("bottom padding error, component:%d, j:%d, i:%d", component, j, i);
			}
//			memcpy(ptr_bottom, ptr, stride*sizeof(ptr[0]));
			ptr_bottom+=stride;
		}

		ptr = WND_DATA_PTR(int16_t *, *wnd, component);
		ptr -= padding_x;
		ptr_top = ptr-padding_y*stride;
		for(j=0;j<padding_y;j++)
		{
			for(i=0;i<stride;i++)
			{
				if(ptr_top[i] != ptr[i])
					printf("top padding error, j:%d, i:%d", j, i);
			}
//			memcpy(ptr_top, ptr, stride*sizeof(ptr[0]));
			ptr_top+=stride;
		}
	}
}

//TEncTop::selectReferencePictureSet
void hmr_select_reference_picture_set(hvenc_enc_t* hvenc, slice_t *currslice)
{
	currslice->ref_pic_set_index = 0;
	
	currslice->ref_pic_set = &hvenc->ref_pic_set_list[currslice->ref_pic_set_index];
	currslice->ref_pic_set->num_pics = currslice->ref_pic_set->num_negative_pics+currslice->ref_pic_set->num_positive_pics;
}

void apply_reference_picture_set(hvenc_enc_t* hvenc, slice_t *currslice)
{
	int i, j;
	video_frame_t *refpic;

	currslice->ref_pic_list_cnt[REF_PIC_LIST_0] = currslice->ref_pic_list_cnt[REF_PIC_LIST_1] = 0;

	for(i=0;i<MAX_NUM_REF;i++)
	{
		refpic = hvenc->reference_picture_buffer[(hvenc->reference_list_index-1-i)&MAX_NUM_REF_MASK];
		if(refpic!=NULL)
		{
			refpic->is_reference = FALSE;
			for(j=0;j<currslice->ref_pic_set->num_pics;j++)
			{
				if(currslice->poc + currslice->ref_pic_set->delta_poc_s0[j] == refpic->temp_info.poc )
				{
					refpic->is_reference = currslice->ref_pic_set->used_by_curr_pic_S0_flag[j];

					if(refpic->is_reference)
					{
						if(currslice->ref_pic_set->delta_poc_s0[j] < 0)
							currslice->ref_pic_list[REF_PIC_LIST_0][currslice->ref_pic_list_cnt[REF_PIC_LIST_0]] = refpic;
						else
							currslice->ref_pic_list[REF_PIC_LIST_1][currslice->ref_pic_list_cnt[REF_PIC_LIST_1]] = refpic;
					}
				}
			}
		}
	}
}

void hmr_slice_init(hvenc_t* ed, picture_t *currpict, slice_t *currslice)
{
	int img_type = currpict->img2encode->img_type;
	currslice->qp =  ed->pict_qp;
	currslice->poc = ed->last_poc;
	currslice->sps = &ed->hvenc->sps;
	currslice->pps = &ed->hvenc->pps;
	currslice->slice_index = 0;
	currslice->curr_cu_address = currslice->first_cu_address = 0;
	currslice->last_cu_address = ed->pict_total_ctu*ed->num_partitions_in_cu;

	currslice->num_ref_idx[REF_PIC_LIST_0] = ed->num_refs_idx_active_list[REF_PIC_LIST_0];
	currslice->num_ref_idx[REF_PIC_LIST_1] = ed->num_refs_idx_active_list[REF_PIC_LIST_1];
	currslice->slice_temporal_mvp_enable_flag = ed->hvenc->sps.temporal_mvp_enable_flag;
	currslice->deblocking_filter_disabled_flag = 0;//enabled
	currslice->slice_loop_filter_across_slices_enabled_flag = 1;//disabled
	currslice->slice_beta_offset_div2 = currslice->pps->beta_offset_div2;
	currslice->slice_beta_offset_div2 = currslice->pps->beta_offset_div2;
	currslice->max_num_merge_candidates = ed->num_merge_mvp_candidates;

//	if((currslice->poc%ed->intra_period)==0)
	if(currslice->poc==0 || (ed->intra_period!=0 && currslice->poc==(ed->last_intra + ed->intra_period) && img_type == IMAGE_AUTO) || (ed->intra_period==0 && currslice->poc==0) || img_type == IMAGE_I)
	{
		ed->last_intra = currslice->poc;
		ed->last_gop_reinit = currslice->poc;
		currpict->img2encode->img_type = IMAGE_I;
		currslice->slice_type = I_SLICE;
		currslice->slice_temporal_layer_non_reference_flag = 0;
		currslice->is_dependent_slice = 0;
		currslice->nalu_type = get_nal_unit_type(ed, currslice, currslice->poc);//NALU_CODED_SLICE_IDR;
		currslice->sublayer = 0;
		currslice->depth = 0;
		currslice->qp = ed->pict_qp;
	}
	else if((ed->num_b==0 && img_type == IMAGE_AUTO) || img_type == IMAGE_P || ed->intra_period==0)
	{
		currslice->slice_type = P_SLICE;
		currpict->img2encode->img_type = IMAGE_P;
		currslice->slice_temporal_layer_non_reference_flag = 0;
		currslice->is_dependent_slice = 0;
		currslice->nalu_type = get_nal_unit_type(ed, currslice, currslice->poc);//NALU_CODED_SLICE_IDR;
		currslice->sublayer = 0;
		currslice->depth = 0;	
		currslice->qp = ed->pict_qp;//+2;
	}

	hmr_select_reference_picture_set(ed->hvenc, currslice);

	if(currslice->nalu_type == NALU_CODED_SLICE_IDR_W_RADL || currslice->nalu_type == NALU_CODED_SLICE_IDR_N_LP)
	{
		ed->hvenc->last_idr = currslice->poc;
	}
	ed->last_idr = ed->hvenc->last_idr;

	if(currslice->poc == 0)//if(ed->last_poc == 0)//first image?
	{
//		sublayer_non_reference = 0;		
	}

	switch(currslice->slice_type)
	{
		case I_SLICE:
		case P_SLICE:
		{
			currslice->referenced = 1;
		}
		break;
	}
}



void CuGetNeighbors(henc_thread_t* et, ctu_info_t* ctu)
{
	if(ctu->x[Y_COMP]==0)
		ctu->ctu_left = NULL;
	else
	{
		ctu->ctu_left = &et->ed->ctu_info[ctu->ctu_number-1];
	}
	ctu->ctu_left_bottom = NULL;//en raster order este no existe. En wavefront si

	if(ctu->y[Y_COMP]==0)
	{
		ctu->ctu_top = NULL;	
		ctu->ctu_top_right = NULL;
		ctu->ctu_top_left = NULL;
	}
	else
	{
		ctu->ctu_top = &et->ed->ctu_info[ctu->ctu_number-et->pict_width_in_ctu];	

		if(ctu->x[Y_COMP]==0)
			ctu->ctu_top_left = NULL;
		else
			ctu->ctu_top_left = &et->ed->ctu_info[ctu->ctu_number-et->pict_width_in_ctu-1];	

		if(et->cu_current_y==0 || ((et->cu_current_x % et->pict_width_in_ctu) == (et->pict_width_in_ctu-1)))
			ctu->ctu_top_right = NULL;
		else
		{
			ctu->ctu_top_right = &et->ed->ctu_info[ctu->ctu_number-et->pict_width_in_ctu+1];
		}
	}
}



int HOMER_enc_write_annex_b_output(nalu_t *nalu_out[], unsigned int num_nalus, encoder_in_out_t *vout)
{
	uint nalu_idx, bytes_written=0;
	uint8_t code[4] = {0x0,0x0,0x0,0x1};

	for (nalu_idx=0;nalu_idx<num_nalus;nalu_idx++)
	{
		nalu_t *nalu = nalu_out[nalu_idx];
		uint size = 0;

		if (nalu_idx==0 || nalu->nal_unit_type == NALU_TYPE_SPS || nalu->nal_unit_type == NALU_TYPE_PPS)
		{
			memcpy(&vout->stream.streams[0][bytes_written], code, 4);
			bytes_written += 4;
			size += 4;
		}
		else
		{
			memcpy(&vout->stream.streams[0][bytes_written], code+1, 3);
			bytes_written += 3;
			size += 3;
		}

		//
		memcpy(&vout->stream.streams[0][bytes_written], nalu->bs.bitstream, nalu->bs.streambytecnt);
		bytes_written += nalu->bs.streambytecnt;
		size += nalu->bs.streambytecnt;

	}
	vout->stream.data_size[0] = bytes_written;
	return bytes_written;
}

void copy_ctu(ctu_info_t* src_ctu, ctu_info_t* dst_ctu)
{
	uint8_t *part_size_type = dst_ctu->part_size_type;
	uint8_t *pred_mode = dst_ctu->pred_mode;
	uint8_t *skipped = dst_ctu->skipped;
	uint8_t *merge = dst_ctu->merge;
	uint8_t *merge_idx = dst_ctu->merge_idx;
	memcpy(dst_ctu, src_ctu, sizeof(ctu_info_t));
	dst_ctu->part_size_type = part_size_type;
	dst_ctu->pred_mode = pred_mode;
	dst_ctu->skipped = skipped;
	dst_ctu->merge = merge;
	dst_ctu->merge_idx = merge_idx;
}


const uint8_t chroma_scale_conversion_table[58]=
{
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
  17,18,19,20,21,22,23,24,25,26,27,28,29,29,30,31,32,
  33,33,34,34,35,35,36,36,37,37,38,39,40,41,42,43,44,
  45,46,47,48,49,50,51
};


ctu_info_t* init_ctu(henc_thread_t* et)
{
	ctu_info_t *ctu;
	int ctu_width, ctu_height;
	
	ctu = &et->ed->ctu_info[et->cu_current];//&et->curr_ctu_group_info[0];
	ctu->ctu_number = et->cu_current;
	ctu->x[Y_COMP] = et->cu_current_x*et->ctu_width[Y_COMP];
	ctu->y[Y_COMP] = et->cu_current_y*et->ctu_height[Y_COMP];
	ctu->x[U_COMP] = ctu->x[V_COMP] = et->cu_current_x*et->ctu_width[U_COMP];
	ctu->y[U_COMP] = ctu->y[V_COMP] = et->cu_current_y*et->ctu_height[U_COMP];
	ctu->size = et->max_cu_size;
	ctu->num_part_in_ctu = et->num_partitions_in_cu;
	ctu->num_part_in_ctu = et->num_partitions_in_cu;
	ctu->partition_list = &et->partition_info[0];

	ctu_width = ((ctu->x[Y_COMP]+ctu->size)<et->pict_width[Y_COMP])?(ctu->size):et->pict_width[Y_COMP]-ctu->x[Y_COMP];
	ctu_height = ((ctu->y[Y_COMP]+ctu->size)<et->pict_height[Y_COMP])?(ctu->size):et->pict_height[Y_COMP]-ctu->y[Y_COMP];

	if(ctu_width!=ctu->size || ctu_height!=ctu->size)
	{
		int width_in_partitions = ctu_width>>2;
		int height_in_partitions = ctu_height>>2;
		int cu_size_in_partitions = ctu->size>>2;
//		if(height_in_partitions == cu_size_in_partitions)
//			height_in_partitions = cu_size_in_partitions-1;
//		else if(width_in_partitions == cu_size_in_partitions)
		height_in_partitions -= 1;
		ctu->last_valid_partition =	et->ed->raster2abs_table[height_in_partitions*cu_size_in_partitions+width_in_partitions-1];
	}
	else
	{
		ctu->last_valid_partition = et->num_partitions_in_cu-1;
	}
	
	CuGetNeighbors(et, ctu);//raster order
	return ctu;
}

THREAD_RETURN_TYPE ctu_encoder_thread(void *h)
{
	henc_thread_t* et = (henc_thread_t*)h;
	int gcnt=0;
	picture_t *currpict = &et->ed->current_pict;
	slice_t *currslice = &currpict->slice;
	ctu_info_t* ctu;

	//printf("		+ctu_encoder_thread %d\r\n", et->index);

	et->acc_dist = 0;
	et->cu_current = 0;
	et->cu_current_x = 0;
	et->cu_current_y = et->index;
	et->num_intra_partitions = 0;

	ctu = &et->ed->ctu_info[et->cu_current_y*et->pict_width_in_ctu];


	if(et->index==0)
	{
		et->ec = &et->ed->ec_list[0];
		et->ee = et->ed->ee_list[0];
		et->ee->bs = &et->ed->aux_bs[0];
		hmr_bitstream_init(et->ee->bs);

		//resetEntropy
		ee_start_entropy_model(et->ee, currslice->slice_type, currslice->qp, et->pps->cabac_init_flag);//Init CABAC contexts
		//cabac - reset binary encoder and entropy
		et->ee->ee_start(et->ee->b_ctx);
		et->ee->ee_reset_bits(et->ee->b_ctx);//ee_reset(&ed->ee);
	}

#define  GRAIN					1
#define  GRAIN_MASK				(GRAIN-1)


	while(et->cu_current < et->pict_total_ctu)//all ctus loop
	{
		int bits_allocated;
		et->cu_current = et->pict_width_in_ctu*(et->cu_current_y)+et->cu_current_x;
		et->cu_next = et->cu_current+min(1,et->pict_width_in_ctu-et->cu_current_x);

		if(et->cu_current_y > 0)// && ((et->cu_current_x & GRAIN_MASK) == 0))
		{
			SEM_WAIT(et->synchro_wait);
		}

		if(et->cu_current_x==0 && et->wfpp_enable)
		{
			if(et->cu_current_y > 0)
			{
				ptrswap(enc_env_t*, et->ed->ee_list[(2*et->index)], et->ed->ee_list[(2*et->index+et->ed->num_ee-1)%et->ed->num_ee]);//get inherited enviroment, leave current as avaliable for the previous thread of the list
			}

			et->ec = &et->ed->ec_list[et->index];
			et->ee = et->ed->ee_list[(2*et->index)];
			et->ee->bs = &et->ed->aux_bs[et->cu_current_y];
			hmr_bitstream_init(et->ee->bs);
			et->ee->ee_start(et->ee->b_ctx);
			et->ee->ee_reset_bits(et->ee->b_ctx);//ee_reset(&ed->ee);
		}

		ctu = init_ctu(et);

/*
		if(et->ed->num_encoded_frames == 14 && ctu->ctu_number>=14)
		{
			for(gcnt=0;gcnt<NUM_QUANT_WNDS;gcnt++)
			{
				int iiiiiii=0;
				wnd_realloc(&et->transform_quant_wnd_[gcnt], (et->ctu_width[0]), et->ctu_height[0], 0, 0, sizeof(int16_t));		
				et->transform_quant_wnd[gcnt] = &et->transform_quant_wnd_[gcnt];//use pointers to exchange windows
			}
		}
*/

		//Prepare Memory
		mem_transfer_move_curr_ctu_group(et, et->cu_current_x, et->cu_current_y);	//move MBs from image to currMbWnd
		mem_transfer_intra_refs(et, ctu);//copy left and top info for intra prediction

		copy_ctu(ctu, et->ctu_rd);

		bits_allocated = hmr_bitstream_bitcount(et->ee->bs);

		PROFILER_RESET(intra)

		//map spatial features and neighbours in recursive partition structure
		create_partition_ctu_neighbours(et, ctu, ctu->partition_list);


		if(currslice->slice_type != I_SLICE && !et->ed->is_scene_change)// && (ctu->ctu_number & 0x1) == 0)
		{
			int ll;

			motion_inter(et, ctu);

			for(ll = 0; ll<et->num_partitions_in_cu;ll++)
			{
				if(ctu->pred_mode[ll]==INTRA_MODE)
					et->num_intra_partitions++;
			}
		}
		else
		{
			//make ctu intra prediction
			motion_intra(et, ctu, gcnt);
			et->num_intra_partitions += et->num_partitions_in_cu;
		}
		PROFILER_ACCUMULATE(intra)

		mem_transfer_decoded_blocks(et, ctu);

		if(et->cu_current_x>=2 && et->cu_current_y+1 != et->pict_height_in_ctu)// && ((et->cu_current_x & GRAIN_MASK) == 0))
		{
			SEM_POST(et->synchro_signal);
		}

#ifndef COMPUTE_AS_HM
		if(et->ed->intra_period>1)// && ctu->ctu_number == et->pict_total_ctu-1)
			hmr_deblock_filter_ctu(et, currslice, ctu);
//		else
//			reference_picture_border_padding_ctu(&et->ed->curr_reference_frame->img, ctu);

//		if(et->ed->intra_period>1 && ctu->ctu_number == et->pict_total_ctu-1)
//			hmr_deblock_filter(et->ed, currslice);

#endif

		//cabac - encode ctu
		PROFILER_RESET(cabac)
		ctu->coeff_wnd = et->transform_quant_wnd[0];

		ee_encode_ctu(et, et->ee, currslice, ctu, gcnt);
		PROFILER_ACCUMULATE(cabac)

		et->acc_dist += ctu->partition_list[0].distortion;
		et->num_encoded_ctus++;
		et->num_bits += hmr_bitstream_bitcount(et->ee->bs)-bits_allocated;
		et->cu_current_x++;

		//notify first synchronization as this line must go two ctus ahead from next line in wfpp
		if(et->cu_current_x==2 && et->cu_current_y+1 != et->pict_height_in_ctu)
		{
			if(et->wfpp_enable)
				ee_copy_entropy_model(et->ee, et->ed->ee_list[(2*et->index+1)%et->ed->num_ee]);
			SEM_POST(et->synchro_signal);
		}


		//notify last synchronization as this line goes two ctus ahead from next line in wfpp
		if(et->cu_current_x==et->pict_width_in_ctu && et->cu_current_y+1 != et->pict_height_in_ctu)// && ((et->cu_current_x & GRAIN_MASK) == 0))
		{
			SEM_POST(et->synchro_signal);
		}
		if(et->cu_current_x==et->pict_width_in_ctu)
		{
/*			if(et->cu_current_y!=0)
			{
//				printf("sem_post-a - line processed = %d\r\n",et->cu_current_y);
				SEM_POST(et->deblock_filter_sem);
			}
			if(et->cu_current+1 == et->pict_total_ctu)
			{
//				printf("sem_post-b - line processed = %d\r\n",et->cu_current_y);
				SEM_POST(et->deblock_filter_sem);			
			}
*/
			if(et->wfpp_enable)
				ee_end_slice(et->ee, currslice, ctu);
			et->cu_current_y+=et->wfpp_num_threads;
			et->cu_current_x=0;
		}

		et->cu_current = et->pict_width_in_ctu*(et->cu_current_y)+et->cu_current_x;
	}
	
	if(!et->wfpp_enable)
	{
		ee_end_slice(et->ee, currslice, ctu);
	}

	return THREAD_RETURN;
}

THREAD_RETURN_TYPE deblocking_filter_thread(void *h)
{
	hvenc_t* ed = (hvenc_t*)h;
	picture_t *currpict = &ed->current_pict;
	slice_t *currslice = &currpict->slice;

	if(ed->intra_period>1)
		hmr_deblock_filter(ed, currslice);

	return THREAD_RETURN;
}

int HOMER_enc_encode(void* handle, encoder_in_out_t* input_frame)
{
	hvenc_enc_t* ed = (hvenc_enc_t*)handle;
	put_frame_to_encode(ed, input_frame);

	return 0;
}

int HOMER_enc_get_coded_frame(void* handle, encoder_in_out_t* output_frame, nalu_t *nalu_out[], unsigned int *nalu_list_size)
{
	hvenc_enc_t* hvenc = (hvenc_enc_t*)handle;
	*nalu_list_size = 0;

	if(get_num_elements(hvenc->output_hmr_container))
	{
		int comp, j, i, stride_src;
		uint16_t *src;
		uint8_t *dst;
		output_set_t* ouput_set;
		cont_get(hvenc->output_hmr_container, (void**)&ouput_set);
		memcpy(nalu_out, ouput_set->nalu_list, ouput_set->num_nalus*sizeof(ouput_set->nalu_list[0]));
//		memcpy(nalu_out, ouput_set->nalu_list, ouput_set->num_nalus*sizeof(ouput_set->nalu_list[0]));
		*nalu_list_size = ouput_set->num_nalus;
		output_frame->pts = ouput_set->pts;
		output_frame->image_type = ouput_set->image_type;
		if(output_frame->stream.streams[0]!=NULL && output_frame->stream.streams[1]!=NULL && output_frame->stream.streams[2]!=NULL)
		{
			for(comp=Y_COMP;comp<=V_COMP;comp++)
			{
				src = WND_DATA_PTR(uint16_t*, ouput_set->frame->img, comp);
				dst = output_frame->stream.streams[comp];
				stride_src = WND_STRIDE_2D(ouput_set->frame->img, comp);
				
				for(j=0;j<hvenc->pict_height[comp];j++)
				{
					for(i=0;i<hvenc->pict_width[comp];i++)
					{
						*dst++ = (uint8_t)src[i];						
					}
					src += stride_src;
				}
			}
		}
	}

	return 0;
}

#define SYNC_THREAD_CONTEXT(ed, et)									\
		et->rd = ed->rd;											

THREAD_RETURN_TYPE encoder_thread(void *h)
{
	hvenc_t* ed = (hvenc_t*)h;
	picture_t *currpict = &ed->current_pict;
	slice_t *currslice = &currpict->slice;
	int n, i, num_threads;

	while(ed->hvenc->run)
	{
		output_set_t* ouput_sets;// = &ed->hvenc->output_sets[ed->num_encoded_frames & NUM_OUTPUT_NALUS_MASK];
		int		output_nalu_cnt = 0;
		int		nalu_list_size = NALU_SET_SIZE;
		nalu_t	**output_nalu_list;// = ouput_sets->nalu_list;

		EnterCriticalSection(&ed->hvenc->CriticalSection); 

		//get next image
		if(!get_frame_to_encode(ed->hvenc, &ed->current_pict.img2encode))//get next image to encode and init type
			return THREAD_RETURN;

		//entra seccion critica
//		ed->num_pictures++;
		ed->last_poc = ed->hvenc->poc;
		ed->hvenc->poc++;
		ed->num_encoded_frames = ed->hvenc->num_encoded_frames;
		ed->hvenc->num_encoded_frames++;

		ouput_sets = &ed->hvenc->output_sets[ed->num_encoded_frames & NUM_OUTPUT_NALUS_MASK];
		output_nalu_list = ouput_sets->nalu_list;

		memset(output_nalu_list, 0, (nalu_list_size)*sizeof(output_nalu_list[0]));

		ed->slice_nalu = &ed->slice_nalu_list[ed->num_encoded_frames & NUM_OUTPUT_NALUS_MASK];
		hmr_bitstream_init(&ed->slice_nalu->bs);

		hmr_slice_init(ed, &ed->current_pict, &currpict->slice);
		hmr_rd_init(ed, &currpict->slice);

		if(ed->bitrate_mode != BR_FIXED_QP)
		{
			if(currslice->poc==0)
				hmr_rc_init_seq(ed);

			ed->rc = ed->hvenc->rc;
			hmr_rc_init_pic(ed, &currpict->slice);
		}
		//get free img for decoded blocks
		cont_get(ed->hvenc->cont_empty_reference_wnds,(void**)&ed->curr_reference_frame);
		ed->curr_reference_frame->temp_info.poc = currslice->poc;//assign temporal info to decoding window for future use as reference

		apply_reference_picture_set(ed->hvenc, currslice);		

		//prunning of references must be done in a selective way
		if(ed->hvenc->reference_picture_buffer[ed->hvenc->reference_list_index]!=NULL)
			cont_put(ed->hvenc->cont_empty_reference_wnds,ed->hvenc->reference_picture_buffer[ed->hvenc->reference_list_index]);

		ed->hvenc->reference_picture_buffer[ed->hvenc->reference_list_index] = ed->curr_reference_frame;
		ed->hvenc->reference_list_index = (ed->hvenc->reference_list_index+1)&MAX_NUM_REF_MASK;

		if(currslice->nalu_type == NALU_CODED_SLICE_IDR_W_RADL)
		{
			hmr_bitstream_init(&ed->hvenc->vps_nalu.bs);
			hmr_bitstream_init(&ed->hvenc->sps_nalu.bs);
			hmr_bitstream_init(&ed->hvenc->pps_nalu.bs);

			ed->hvenc->vps_nalu.nal_unit_type = NALU_TYPE_VPS;
			ed->hvenc->vps_nalu.temporal_id = ed->hvenc->vps_nalu.rsvd_zero_bits = 0;
			output_nalu_list[output_nalu_cnt++] = &ed->hvenc->vps_nalu;
			hmr_put_vps_header(ed->hvenc);//vps header

			ed->hvenc->sps_nalu.nal_unit_type = NALU_TYPE_SPS;
			ed->hvenc->sps_nalu.temporal_id = ed->hvenc->sps_nalu.rsvd_zero_bits = 0;
			output_nalu_list[output_nalu_cnt++] = &ed->hvenc->sps_nalu;
			hmr_put_seq_header(ed->hvenc);//seq header

			ed->hvenc->pps_nalu.nal_unit_type = NALU_TYPE_PPS;
			ed->hvenc->pps_nalu.temporal_id = ed->hvenc->pps_nalu.rsvd_zero_bits = 0;
			output_nalu_list[output_nalu_cnt++] = &ed->hvenc->pps_nalu;
			hmr_put_pic_header(ed->hvenc);//pic header
		}

		ed->avg_dist = ed->hvenc->avg_dist;

		for(n = 0; n<ed->wfpp_num_threads;n++)
		{
			SYNC_THREAD_CONTEXT(ed, ed->thread[n]);
			//reset remafore
			SEM_RESET(ed->thread[n]->synchro_wait)
		}

		SEM_RESET(ed->deblock_filter_sem)

#ifdef COMPUTE_AS_HM
		CREATE_THREADS((&ed->hthreads[0]), ctu_encoder_thread, ed->thread, ed->wfpp_num_threads)
		JOIN_THREADS(ed->hthreads, ed->wfpp_num_threads)
#else
		num_threads = ed->wfpp_num_threads;// + 1;//wfpp + deblocking thread

//		CREATE_THREAD(ed->hthreads[0], deblocking_filter_thread, ed);
		CREATE_THREADS((&ed->hthreads[0]), ctu_encoder_thread, ed->thread, ed->wfpp_num_threads)

		JOIN_THREADS(ed->hthreads, num_threads)//ed->wfpp_num_threads)		
#endif
		//calc average distortion
		if(ed->num_encoded_frames == 0 || currslice->slice_type != I_SLICE || ed->intra_period==1)
		{
			ed->avg_dist = 0;
			for(n = 0;n<ed->wfpp_num_threads;n++)
			{
				henc_thread_t* henc_th = ed->thread[n];
		
				 ed->avg_dist+= henc_th->acc_dist;
			}
			ed->avg_dist /= ed->pict_total_ctu*ed->num_partitions_in_cu;
			ed->avg_dist = clip(ed->avg_dist,.1,ed->avg_dist);
			if(currslice->slice_type == I_SLICE)
				ed->avg_dist*=1.5;
			else if(ed->is_scene_change)
				ed->avg_dist*=1.375;
		}
		else if(currslice->slice_type == I_SLICE)
			ed->avg_dist/=1.5;

		if(ed->bitrate_mode != BR_FIXED_QP)
			hmr_rc_end_pic(ed, currslice);

		ed->hvenc->avg_dist = ed->avg_dist;
		ed->is_scene_change = 0;

//#ifdef COMPUTE_AS_HM
//		if(ed->intra_period>1)
//			hmr_deblock_filter(ed, currslice);
//#endif
		LeaveCriticalSection(&ed->hvenc->CriticalSection);

		EnterCriticalSection(&ed->hvenc->CriticalSection2); 
		//slice header
		ed->slice_nalu->nal_unit_type = currslice->nalu_type;
		ed->slice_nalu->temporal_id = ed->slice_nalu->rsvd_zero_bits = 0;
		output_nalu_list[output_nalu_cnt++] = ed->slice_nalu;

		hmr_put_slice_header(ed, currslice);//slice header
		if(ed->wfpp_enable)
			hmr_slice_header_code_wfpp_entry_points(ed);
		hmr_bitstream_rbsp_trailing_bits(&ed->slice_bs);

		for(i=0;i<ed->num_sub_streams;i++)
		{
			memcpy(&ed->slice_bs.bitstream[ed->slice_bs.streambytecnt], ed->aux_bs[i].bitstream, ed->aux_bs[i].streambytecnt);
			ed->slice_bs.streambytecnt += ed->aux_bs[i].streambytecnt;
		}

		//escribimos la nalu
		hmr_bitstream_put_nal_unit_header(&ed->slice_nalu->bs, currslice->nalu_type, 0, 0);
		hmr_bitstream_nalu_ebsp(&ed->slice_bs,&ed->slice_nalu->bs);

#ifdef WRITE_REF_FRAMES
		wnd_write2file(&ed->curr_reference_frame->img);
#endif

//		ed->num_encoded_frames++;
#ifdef DBG_TRACE
		{
			char stringI[] = "I";
			char stringP[] = "P";
			char stringB[] = "B";
			char *frame_type_str;
			frame_type_str=currpict->img2encode->img_type==IMAGE_I?stringI:currpict->img2encode->img_type==IMAGE_P?stringP:stringB;

			printf("\r\nmodule:%d, frame:%d, %s, bits:%d,", ed->index,ed->num_encoded_frames, frame_type_str, ed->slice_bs.streambytecnt*8);
#ifdef COMPUTE_METRICS

			homer_psnr(&ed->current_pict, &ed->curr_reference_frame->img, ed->pict_width, ed->pict_height, ed->current_psnr); 
			ed->accumulated_psnr[0] += ed->current_psnr[Y_COMP];
			ed->accumulated_psnr[1] += ed->current_psnr[U_COMP];
			ed->accumulated_psnr[2] += ed->current_psnr[V_COMP];

			printf("PSNRY: %.2f, PSNRU: %.2f,PSNRV: %.2f, ", ed->current_psnr[Y_COMP], ed->current_psnr[U_COMP], ed->current_psnr[V_COMP]);
			printf("Average PSNRY: %.2f, PSNRU: %.2f,PSNRV: %.2f, ", ed->accumulated_psnr[Y_COMP]/(ed->num_encoded_frames+1), ed->accumulated_psnr[U_COMP]/(ed->num_encoded_frames+1), ed->accumulated_psnr[V_COMP]/(ed->num_encoded_frames+1));
#endif
			printf("vbv: %.2f, avg_dist: %.2f, ", ed->rc.vbv_fullness/ed->rc.vbv_size, ed->avg_dist);
			printf("rc.target_pict_size: %.2f", ed->rc.target_pict_size);
			fflush(stdout);
		}
#endif
		//prunning of references must be done in a selective way
//		if(ed->reference_picture_buffer[ed->reference_list_index]!=NULL)
//			cont_put(ed->cont_empty_reference_wnds,ed->reference_picture_buffer[ed->reference_list_index]);

		//fill padding in reference picture
//		reference_picture_border_padding(&ed->curr_reference_frame->img);
//		ed->reference_picture_buffer[ed->reference_list_index] = ed->curr_reference_frame;
//		ed->reference_list_index = (ed->reference_list_index+1)&MAX_NUM_REF_MASK;
//		ed->last_poc++;


		ouput_sets->pts = ed->current_pict.img2encode->temp_info.pts;
		ouput_sets->image_type = ed->current_pict.img2encode->img_type;
		put_avaliable_frame(ed->hvenc, ed->current_pict.img2encode);

		ouput_sets->num_nalus = output_nalu_cnt;
		ouput_sets->frame = ed->curr_reference_frame;
		cont_put(ed->hvenc->output_hmr_container, ouput_sets);
		LeaveCriticalSection(&ed->hvenc->CriticalSection2);

	}

	return THREAD_RETURN;
}

