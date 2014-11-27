/*****************************************************************************
 * hmr_headers.c : homerHEVC encoding library
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

#include <math.h>
#include "hmr_common.h"
#include "hmr_private.h"



void hmr_put_profile_tier(bitstream_t *bs, profile_tier_t* pt )
{
	int j;

	hmr_bitstream_write_bits(bs, pt->profileSpace, 2);
	hmr_bitstream_write_bits(bs, pt->tierFlag, 1);
	hmr_bitstream_write_bits(bs, pt->profileIdc, 5);
	for(j = 0; j < 32; j++)
	{
		hmr_bitstream_write_bits(bs, pt->profileCompatibilityFlag[j], 1);
	}

	hmr_bitstream_write_bits(bs, pt->general_progressive_source_flag, 1);
	hmr_bitstream_write_bits(bs, pt->general_interlaced_source_flag, 1);
	hmr_bitstream_write_bits(bs, pt->general_non_packed_constraint_flag, 1);
	hmr_bitstream_write_bits(bs, pt->general_frame_only_constraint_flag, 1);
  
	hmr_bitstream_write_bits(bs, 0, 16);//reserved zero bits
	hmr_bitstream_write_bits(bs, 0, 16);//reserved zero bits
	hmr_bitstream_write_bits(bs, 0, 12);//reserved zero bits
}

void hmr_put_pic_bitrate_info(bitstream_t *bs, int temp_level_low, int temp_level_high)
{
	int i;
	for(i = temp_level_low; i <= temp_level_high; i++)
	{
		hmr_bitstream_write_bits(bs, 0, 1);
		hmr_bitstream_write_bits(bs, 0, 1);
	}
}
void hmr_put_profile_and_level(bitstream_t *bs, profile_tier_level_t* ptl, int profilePresentFlag, int maxNumSubLayersMinus1)
{
	int i;
	if(profilePresentFlag)
	{
		hmr_put_profile_tier(bs, &ptl->generalPTL);    // general_...
	}
	hmr_bitstream_write_bits(bs, ptl->generalPTL.levelIdc, 8);
 

	for(i = 0; i < maxNumSubLayersMinus1; i++)
	{
		hmr_bitstream_write_bits(bs, ptl->subLayerProfilePresentFlag[i], 1);
		hmr_bitstream_write_bits(bs, ptl->subLayerLevelPresentFlag[i], 1);
		if( profilePresentFlag && ptl->subLayerProfilePresentFlag[i])
		{
			hmr_put_profile_tier(bs, &ptl->subLayerPTL[i]);  // sub_layer_...
		}
		if( ptl->subLayerLevelPresentFlag[i] )
		{
			hmr_bitstream_write_bits(bs, ptl->subLayerPTL[i].levelIdc, 8);
		}
	}
}

void hmr_put_vps_header(hvenc_t* ed)
{
	int i;
	vps_t	*vps = &ed->vps;
	bitstream_t	*bs = ed->aux_bs;

	hmr_bitstream_init(bs);
	hmr_bitstream_write_bits(bs, vps->video_parameter_set_id, 4);
	hmr_bitstream_write_bits(bs, 3, 2);//rsvd three 2 bits
	hmr_bitstream_write_bits(bs, 0, 6);//rsvd zero 6 bits
	hmr_bitstream_write_bits(bs, ed->max_sublayers-1, 3);
	hmr_bitstream_write_bits(bs, vps->temporal_id_nesting_flag, 1);
	hmr_bitstream_write_bits(bs, 0xffff, 16);//vps_reserved_ffff_16bits

	hmr_put_profile_and_level(bs, vps->ptl, TRUE, ed->max_sublayers-1 );

//	hmr_put_pic_bitrate_info(bs, 0, ed->max_sublayers-1);

	hmr_bitstream_write_bits(bs, vps->sub_layer_ordering_info_present_flag, 1);//

	for(i=0; i <= ed->max_sublayers-1; i++)
	{
		hmr_bitstream_write_bits_uvlc(bs, vps->max_dec_pic_buffering[i]-1);//vps_max_dec_pic_buffering_minus1
		hmr_bitstream_write_bits_uvlc(bs, vps->max_num_reorder_pics[i]);//vps_max_num_reorder_pics
		hmr_bitstream_write_bits_uvlc(bs, vps->max_latency_increase[i]);//vps_max_latency_increase_plus1
		if (!vps->sub_layer_ordering_info_present_flag)
		{
			break;
		}
	}

	hmr_bitstream_write_bits(bs, 0, 6);//vps_max_layer_id u(6) // max_nuh_reserved_zero_layer_id
	hmr_bitstream_write_bits_uvlc(bs, vps->num_layer_sets_minus1);//vps_num_layer_sets_minus1 ue(v) //vps_max_op_sets_minus1

//	for( i = 1; i <= vps.num_layer_sets_minus1; i++ )
//	for( j = 0; j <= vps.max_layer_id; j++ )
//		vps.layer_id_included_flag[ i ][ j ]

	hmr_bitstream_write_bits(bs, vps->timing_info_present_flag, 1);//

	if(vps->timing_info_present_flag)
	{
		//.....

		hmr_bitstream_write_bits_uvlc(bs, 0);//num_hdr_parameters
	}

	hmr_bitstream_write_bits(bs, 0, 1);//extension flag

	hmr_bitstream_rbsp_trailing_bits(bs);

	hmr_bitstream_put_nal_unit_header(&ed->vps_nalu.bs, NALU_TYPE_VPS, 0, 0);
	hmr_bitstream_nalu_ebsp(bs, &ed->vps_nalu.bs);
}

void hmr_short_term_ref_pic_set(bitstream_t	*bs, ref_pic_set_t *rps, int idx)
{
	int j;
	if (idx > 0)
	{
		hmr_bitstream_write_bits(bs, rps->inter_ref_pic_set_prediction_flag, 1);//inter_ref_pic_set_prediction_flag
	}
	if (rps->inter_ref_pic_set_prediction_flag)
	{
	
	}
	else
	{
		int prev = 0;
		hmr_bitstream_write_bits_uvlc(bs, rps->num_negative_pics);
		hmr_bitstream_write_bits_uvlc(bs, rps->num_positive_pics);
		for(j=0;j<rps->num_negative_pics;j++)
		{
			hmr_bitstream_write_bits_uvlc(bs, prev - rps->delta_poc_s0[j]-1);//delta_poc_s0_minus1
			prev = rps->delta_poc_s0[j];
			hmr_bitstream_write_bits(bs, rps->used_by_curr_pic_S0_flag[j], 1);//inter_ref_pic_set_prediction_flag
		}
		for(j=0;j<rps->num_positive_pics;j++)
		{
			//...................
		}
	}
}


void hmr_put_seq_header(hvenc_t* ed)
{
	sps_t	*sps = &ed->sps;
	bitstream_t	*bs = ed->aux_bs;//&ed->sps_nalu.bs;
	int i;
	hmr_bitstream_init(bs);

	hmr_bitstream_write_bits(bs, sps->video_parameter_set_id, 4);
	hmr_bitstream_write_bits(bs, ed->max_sublayers-1, 3);
	hmr_bitstream_write_bits(bs, sps->temporal_id_nesting_flag, 1);
		
	hmr_put_profile_and_level(bs, sps->ptl, TRUE, ed->max_sublayers-1 );

	hmr_bitstream_write_bits_uvlc(bs, sps->seq_parameter_set_id);
	hmr_bitstream_write_bits_uvlc(bs, sps->chroma_format_idc);

	if(sps->chroma_format_idc==3)
		hmr_bitstream_write_bits(bs, 0, 1);//separate_colour_plane_flag

	hmr_bitstream_write_bits_uvlc(bs, sps->pic_width_in_luma_samples);
	hmr_bitstream_write_bits_uvlc(bs, sps->pic_height_in_luma_samples);
	hmr_bitstream_write_bits(bs, sps->conformance_window_flag, 1);
	if(sps->conformance_window_flag)
	{
		hmr_bitstream_write_bits_uvlc(bs, sps->conf_win_left_offset);
		hmr_bitstream_write_bits_uvlc(bs, sps->conf_win_right_offset);
		hmr_bitstream_write_bits_uvlc(bs, sps->conf_win_top_offset);
		hmr_bitstream_write_bits_uvlc(bs, sps->conf_win_bottom_offset);
	}

	hmr_bitstream_write_bits_uvlc(bs, sps->bit_depth_luma_minus8);
	hmr_bitstream_write_bits_uvlc(bs, sps->bit_depth_chroma_minus8);

	hmr_bitstream_write_bits_uvlc(bs, sps->log2_max_pic_order_cnt_lsb_minus4);

	hmr_bitstream_write_bits(bs, 1, 1);//sub_layer_ordering_info_present_flag

	for(i=0; i <= ed->max_sublayers-1; i++)
	{
		hmr_bitstream_write_bits_uvlc(bs, sps->max_dec_pic_buffering[i]-1);//vps_max_dec_pic_buffering_minus1
		hmr_bitstream_write_bits_uvlc(bs, sps->max_num_reorder_pics[i]);
		hmr_bitstream_write_bits_uvlc(bs, sps->max_latency_increase[i]);
	}

	hmr_bitstream_write_bits_uvlc(bs, sps->log2_min_coding_block_size_minus3);
	hmr_bitstream_write_bits_uvlc(bs, sps->log2_diff_max_min_coding_block_size);
	hmr_bitstream_write_bits_uvlc(bs, sps->log2_min_transform_block_size_minus2);
	hmr_bitstream_write_bits_uvlc(bs, sps->log2_diff_max_min_transform_block_size);

	hmr_bitstream_write_bits_uvlc(bs, sps->max_transform_hierarchy_depth_inter);
	hmr_bitstream_write_bits_uvlc(bs, sps->max_transform_hierarchy_depth_intra);

	hmr_bitstream_write_bits(bs, sps->scaling_list_enabled_flag, 1);
	if(sps->scaling_list_enabled_flag)
	{
		hmr_bitstream_write_bits(bs, sps->scaling_list_data_present_flag, 1);
		if(sps->scaling_list_data_present_flag)
		{
			//code_scaling_list
		}
	}

	hmr_bitstream_write_bits(bs, sps->amp_enabled_flag, 1);
	hmr_bitstream_write_bits(bs, sps->sample_adaptive_offset_enabled_flag, 1);

	hmr_bitstream_write_bits(bs, sps->pcm_enabled_flag, 1);

	if(sps->pcm_enabled_flag)
	{
		//
	}

	hmr_bitstream_write_bits_uvlc(bs, ed->num_short_term_ref_pic_sets);
	for(i=0;i<ed->num_short_term_ref_pic_sets;i++)
	{
		hmr_short_term_ref_pic_set(bs, &ed->ref_pic_set_list[i], i);
	}
	
	hmr_bitstream_write_bits(bs, ed->num_long_term_ref_pic_sets?1:0,1);//long_term_ref_pics_present_flag
	if(ed->num_long_term_ref_pic_sets=!0)
	{
/*		hmr_bitstream_write_bits_uvlc(bs, ed->num_long_term_ref_pic_sets);
		for(i=0;i<ed->num_long_term_ref_pic_sets;i++)
		{
			sps->lt_ref_pic_poc_lsb_sps[]
			sps->used_by_curr_pic_lt_sps_flag[]
			//......
		}	
*/	}
	hmr_bitstream_write_bits(bs, sps->temporal_mvp_enable_flag, 1);
	hmr_bitstream_write_bits(bs, sps->strong_intra_smooth_enabled_flag, 1);
	hmr_bitstream_write_bits(bs, sps->vui_parameters_present_flag, 1);

//esto ahora se ha pasado al vui	
//	hmr_bitstream_write_bits(bs, sps->restricted_ref_pic_lists_flag, 1);
//	if(sps->restricted_ref_pic_lists_flag)
//		hmr_bitstream_write_bits(bs, sps->lists_modification_present_flag, 1);


	hmr_bitstream_write_bits(bs, 0, 1);	//sps extension flag

	hmr_bitstream_rbsp_trailing_bits(bs);

	hmr_bitstream_put_nal_unit_header(&ed->sps_nalu.bs, NALU_TYPE_SPS, 0, 0);
	hmr_bitstream_nalu_ebsp(bs, &ed->sps_nalu.bs);
}


void hmr_put_pic_header(hvenc_t* ed)
{
	sps_t	*sps = &ed->sps;
	pps_t	*pps = &ed->pps;
	bitstream_t	*bs = ed->aux_bs;//&ed->pps_nalu.bs;

	hmr_bitstream_init(bs);
	hmr_bitstream_write_bits_uvlc(bs, pps->pic_parameter_set_id);
	hmr_bitstream_write_bits_uvlc(bs, pps->seq_parameter_set_id);
	hmr_bitstream_write_bits(bs, pps->dependent_slice_enabled_flag, 1);
	hmr_bitstream_write_bits(bs, pps->output_flag_present_flag, 1);
	hmr_bitstream_write_bits(bs, pps->num_extra_slice_header_bits, 3);
	hmr_bitstream_write_bits(bs, pps->sign_data_hiding_flag,1);
	hmr_bitstream_write_bits(bs, pps->cabac_init_present_flag,1);

	hmr_bitstream_write_bits_uvlc(bs, ed->num_refs_idx_active_list[REF_PIC_LIST_0]-1);//num_ref_idx_l0_default_active_minus1
	hmr_bitstream_write_bits_uvlc(bs, ed->num_refs_idx_active_list[REF_PIC_LIST_1]-1);//num_ref_idx_l1_default_active_minus1
	
	hmr_bitstream_write_bits_svlc(bs, pps->pic_init_qp_minus26);

	hmr_bitstream_write_bits(bs, pps->constrained_intra_pred_flag,1);
	hmr_bitstream_write_bits(bs, pps->transform_skip_enabled_flag,1);
	hmr_bitstream_write_bits(bs, pps->cu_qp_delta_enabled_flag,1);
	if(pps->cu_qp_delta_enabled_flag)
		hmr_bitstream_write_bits_uvlc(bs, pps->diff_cu_qp_delta_depth);
	
	hmr_bitstream_write_bits_svlc(bs, pps->cb_qp_offset);
	hmr_bitstream_write_bits_svlc(bs, pps->cr_qp_offset);

	hmr_bitstream_write_bits(bs, pps->slice_chroma_qp_offsets_present_flag,1);
	hmr_bitstream_write_bits(bs, pps->weighted_pred_flag,1);
	hmr_bitstream_write_bits(bs, pps->weighted_bipred_flag,1);
	hmr_bitstream_write_bits(bs, pps->transquant_bypass_enable_flag,1);
	hmr_bitstream_write_bits(bs, pps->tiles_enabled_flag,1);
	hmr_bitstream_write_bits(bs, pps->entropy_coded_sync_enabled_flag,1);
	if(pps->tiles_enabled_flag)
	{
		//....................		
	}

	hmr_bitstream_write_bits(bs, pps->loop_filter_across_slices_enabled_flag, 1);
	hmr_bitstream_write_bits(bs, pps->deblocking_filter_control_present_flag, 1);
	if(pps->deblocking_filter_control_present_flag)
	{
		//....................		
	}
	hmr_bitstream_write_bits(bs, pps->pps_scaling_list_data_present_flag, 1);
	if(pps->pps_scaling_list_data_present_flag)
	{
		//....................		
	}
	hmr_bitstream_write_bits(bs, pps->lists_modification_present_flag, 1);
	hmr_bitstream_write_bits_uvlc(bs, pps->log2_parallel_merge_level_minus2);
//	hmr_bitstream_write_bits(bs, pps->num_extra_slice_header_bits, 3);
	hmr_bitstream_write_bits(bs, pps->slice_header_extension_present_flag, 1);
	hmr_bitstream_write_bits(bs, 0, 1);//pps_extension_flag	

	hmr_bitstream_rbsp_trailing_bits(bs);

	hmr_bitstream_put_nal_unit_header(&ed->pps_nalu.bs, NALU_TYPE_PPS, 0, 0);
	hmr_bitstream_nalu_ebsp(bs, &ed->pps_nalu.bs);
}

void hmr_put_slice_header(hvenc_t* ed, slice_t *currslice)
{
	bitstream_t	*bs = &ed->slice_bs;//
	sps_t	*sps = currslice->sps;
	pps_t	*pps = currslice->pps;
	int max_add_outer = ed->pict_total_ctu;
	int bits_outer = 0;
	int max_addr_inner, req_bits_inner = 0;
	int slice_address;
	int inner_address;
	int address;

	hmr_bitstream_init(bs);
	//calculate number of bits required for slice address
	while(max_add_outer>(1<<bits_outer)) 
	{
		bits_outer++;
	}

	max_addr_inner = 1;
	while(max_addr_inner>(1<<req_bits_inner))
	{
		req_bits_inner++;
	}
    //slice address
	slice_address = currslice->first_cu_address/ed->num_partitions_in_cu;
	inner_address = 0;

	address = (slice_address<<req_bits_inner) + inner_address;

	hmr_bitstream_write_bits(bs, address==0, 1);//first_slice_in_pic_flag

	if(rap_pic_flag(currslice->nalu_type))
	{
		hmr_bitstream_write_bits(bs, 0, 1);//no_output_of_prior_pics_flag
	}
	hmr_bitstream_write_bits_uvlc(bs, pps->pic_parameter_set_id);

	if ( pps->dependent_slice_enabled_flag && address!=0 )
		hmr_bitstream_write_bits(bs, currslice->is_dependent_slice, 1);//dependent_slice_flag

	if(address>0)//if(!first_slice_in_pic_flag)
		hmr_bitstream_write_bits(bs, address, bits_outer+req_bits_inner);//slice_address

	if(!currslice->is_dependent_slice)//if(!dependent_slice_flag )
	{
		int i; 
		for (i = 0; i < pps->num_extra_slice_header_bits; i++)
		{
			hmr_bitstream_write_bits(bs, 0, 1);//reserved_undetermined_flag
		}

		hmr_bitstream_write_bits_uvlc(bs, currslice->slice_type);
/*		if(pps->output_flag_present_flag)
		{
			hmr_bitstream_write_bits(bs, pic_output_flag, 1);//pic_output_flag
		}		
		if( separate_colour_plane_flag = = 1 )
			hmr_bitstream_write_bits(bs, colour_plane_id
*/
		if(idr_pic_flag(currslice->nalu_type))
		{
			
		}
		if(sps->sample_adaptive_offset_enabled_flag)
		{
		
		}

		if(!idr_pic_flag(currslice->nalu_type))//K0251
		{	
			int num_bits = 0;
			int bits_for_poc = currslice->sps->log2_max_pic_order_cnt_lsb_minus4+4;
			int poc_lsb = (currslice->poc - ed->last_idr+(1<<bits_for_poc))%(1<<bits_for_poc);

			hmr_bitstream_write_bits(bs, poc_lsb, bits_for_poc);

			if(currslice->ref_pic_set_index<0)
			{
				hmr_bitstream_write_bits(bs, 0, 1);//short_term_ref_pic_set_sps_flag
				//hmr_short_term_ref_pic_set( num_short_term_ref_pic_sets )
			}
			else
			{
				hmr_bitstream_write_bits(bs, 1, 1);//short_term_ref_pic_set_sps_flag
		
				while ((1 << num_bits) < ed->num_short_term_ref_pic_sets)//esto se puede tener en una tabla
				{
				  num_bits++;
				}

				if(num_bits)
					hmr_bitstream_write_bits(bs, currslice->ref_pic_set_index, num_bits);//short_term_ref_pic_set_sps_flag

/*				hmr_bitstream_write_bits_uvlc(bs, ed->num_short_term_ref_pic_sets);
				for(i=0;i<ed->num_short_term_ref_pic_sets;i++)
				{
					hmr_short_term_ref_pic_set(bs, &ed->ref_pic_set_list[i], i);
				}
*/
			}

			//if( long_term_ref_pics_present_flag )
			//........

			if(sps->temporal_mvp_enable_flag)
			{
				hmr_bitstream_write_bits(bs, currslice->slice_temporal_mvp_enable_flag, 1);//slice_temporal_mvp_enable_flag
			}
		}


		//if(use_sao)
		//.....

		if(!isIntra(currslice->slice_type))
		{
			int override_flag = (currslice->num_ref_idx[REF_PIC_LIST_0]!=(pps->num_ref_idx_l0_default_active_minus1+1));// ||(pcSlice->isInterB()&&pcSlice->getNumRefIdx( REF_PIC_LIST_1 )!=pcSlice->getPPS()->getNumRefIdxL1DefaultActive()));
			hmr_bitstream_write_bits(bs, override_flag, 1);//num_ref_idx_active_override_flag
//			WRITE_FLAG( overrideFlag ? 1 : 0,                               "num_ref_idx_active_override_flag");
			if (override_flag) 
			{
				hmr_bitstream_write_bits_uvlc(bs, pps->num_ref_idx_l0_default_active_minus1);
//				WRITE_UVLC( pcSlice->getNumRefIdx( REF_PIC_LIST_0 ) - 1,      "num_ref_idx_l0_active_minus1" );
				if (currslice->slice_type != B_SLICE)
				{
					hmr_bitstream_write_bits_uvlc(bs, pps->num_ref_idx_l0_default_active_minus1);
//					WRITE_UVLC( pcSlice->getNumRefIdx( REF_PIC_LIST_1 ) - 1,    "num_ref_idx_l1_active_minus1" );
				}
				else
				{
					currslice->num_ref_idx[REF_PIC_LIST_1] = 0;
				}
			}			
		}
		else
		{
			currslice->num_ref_idx[REF_PIC_LIST_0] = 0;
			currslice->num_ref_idx[REF_PIC_LIST_1] = 0;
		}

//		if( sps->lists_modification_present_flag )
//			ref_pic_list_modification( )

		if(currslice->slice_type == B_SLICE) 
		{
//			mvd_l1_zero_flag	
		}

		if(!isIntra(currslice->slice_type))
		{
			if(pps->cabac_init_present_flag)
			{
				//	cabac_init_flag
/*				SliceType sliceType   = pcSlice->getSliceType();
				Int  encCABACTableIdx = pcSlice->getPPS()->getEncCABACTableIdx();
				Bool encCabacInitFlag = (sliceType!=encCABACTableIdx && encCABACTableIdx!=I_SLICE) ? true : false;
				pcSlice->setCabacInitFlag( encCabacInitFlag );
				WRITE_FLAG( encCabacInitFlag?1:0, "cabac_init_flag" );
*/			}
		}
		if(currslice->slice_temporal_mvp_enable_flag)
		{
			if(currslice->slice_type == B_SLICE) 
			{
				//collocated_from_l0_flag
			}			
			//......
		}
		if((pps->weighted_pred_flag && currslice->slice_type == P_SLICE) || (pps->weighted_bipred_flag && currslice->slice_type == B_SLICE))
		{
			//pred_weight_table( )
		}

		if(!isIntra(currslice->slice_type))
		{
			hmr_bitstream_write_bits_uvlc(bs, 5 - currslice->max_num_merge_candidates);//five_minus_max_num_merge_cand			
		}
		hmr_bitstream_write_bits_svlc(bs, currslice->qp/*ed->pict_qp*/ - (pps->pic_init_qp_minus26 + 26));//slice_qp_delta
		
		if(pps->slice_chroma_qp_offsets_present_flag)
		{
			//hmr_bitstream_write_bits_svlc(slice_qp_delta_cb
			//hmr_bitstream_write_bits_svlc(slice_qp_delta_cr
		}

		if(pps->deblocking_filter_control_present_flag)
		{
/*			if(pps->deblocking_filter_override_enabled_flag)
			{
				currslice->deblocking_filter_override_flag
			}
			if(currslice->deblocking_filter_override_flag)
			{
//				hmr_bitstream_write_bits(bs, currslice->deblocking_filter_disabled_flag, 1);
			}
			//............
*/		}
		if(pps->loop_filter_across_slices_enabled_flag && ( /*slice_sao_luma_flag || slice_sao_chroma_flag ||*/ !currslice->deblocking_filter_disabled_flag))
		{
			hmr_bitstream_write_bits(bs, currslice->slice_loop_filter_across_slices_enabled_flag, 1);
		}
	}
	if(pps->slice_header_extension_present_flag)
	{
		hmr_bitstream_write_bits_uvlc(bs, 0);
	}
}


uint count_needed_start_codes(bitstream_t	*bs)
{
	uint cnt = 0;
	unsigned char* ptr = bs->bitstream;
	int size = bs->streambytecnt;
	int i = 0;
	int found=0;

	while(i<size)
	{
		while(i<size)
		{
			if(ptr[i]==0 && ptr[i+1]==0)//search 00,00
			{
				i++;
				if(i==size)
					break;
				if(ptr[++i]<=3)
					break;
			}
			else
				i++;
		}
		if(i<size)
			cnt++;

	}
	return cnt;
}


void hmr_slice_header_code_wfpp_entry_points(hvenc_t* ed)
{
	bitstream_t	*bs = &ed->slice_bs;//
	int i;
	uint max_offset = 0, offset_len_minus_1 = 1;
	uint num_entry_point_offsets = ed->num_sub_streams-1;

	for (i=0; i<num_entry_point_offsets; i++)
	{
		ed->sub_streams_entry_point_list[ i ] = ed->aux_bs[i].streambytecnt + count_needed_start_codes(&ed->aux_bs[i]);
		if ( ed->sub_streams_entry_point_list[ i ] > max_offset)
		{
			max_offset = ed->sub_streams_entry_point_list[ i ];
		}
	}

	while (max_offset >= (1u << (offset_len_minus_1 + 1)))
	{
		offset_len_minus_1++;
	}

	hmr_bitstream_write_bits_uvlc(bs, num_entry_point_offsets);//num_entry_point_offsets
	if(num_entry_point_offsets>0)
		hmr_bitstream_write_bits_uvlc(bs, offset_len_minus_1);//offset_len_minus1

	for (i=0; i<num_entry_point_offsets; i++)
	{
		hmr_bitstream_write_bits(bs, ed->sub_streams_entry_point_list[ i ]-1, offset_len_minus_1+1);//entry_point_offset_minus1
	}
}
