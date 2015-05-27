/*****************************************************************************
 * homer_hevc_enc_api.c : homerHEVC encoding library
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

#ifndef __HOMER_HEVC_ENCODER_API_H__
#define __HOMER_HEVC_ENCODER_API_H__

#ifdef _MSC_VER
#include "..\win_includes\stdint.h"
#else
#include <stdint.h>
#endif

#define CHROMA420 1
#define CHROMA422 2
#define CHROMA444 3


#ifdef	__cplusplus
extern "C" {
#endif




enum HOMER_CMD {
	HENC_SETCFG
};


enum HOMER_PROFILE {
	PROFILE_MAIN = 1,
	PROFILE_MAIN10
};

enum HOMER_TIER {
  TIER_MAIN = 0,
  TIER_HIGH
};

enum HOMER_RD_MODES {
	RD_DIST_ONLY = 0,
	RD_FULL,
	RD_FAST,
	NUM_RD_MODES
};

enum HOMER_BR_MODES {
	BR_FIXED_QP = 0,
	BR_CBR,
	BR_VBR,
	NUM_BR_MODES
};

enum HOMER_PERFORMANCE_MODES {
	PERF_FULL_COMPUTATION = 0,
	PERF_FAST_COMPUTATION,
	PERF_UFAST_COMPUTATION,//fastest
	NUM_PERF_MODES
};

enum HOMER_IMG_TYPES {
  IMAGE_AUTO,
  IMAGE_B,
  IMAGE_P,
  IMAGE_I
};

enum HOMER_MOTION_ESTIMATION_PRECISSION {
  PEL,
  HALF_PEL,
  QUARTER_PEL
};

#define MAX_STREAMS	8

typedef struct stream_t stream_t;
struct  stream_t
{
	uint8_t		*streams[MAX_STREAMS];
	int32_t 	stream_size[MAX_STREAMS];
	int32_t 	data_size[MAX_STREAMS];
	int32_t 	data_stride[MAX_STREAMS];
};


typedef struct encoder_in_out_t encoder_in_out_t;
struct encoder_in_out_t
{
	stream_t	stream;
	uint64_t	pts;
	uint32_t	image_type;//HOMER_IMG_TYPES - this field allows to force the type of an image whenever is needed. otherwise use IMAGE_AUTO
	int8_t		user_data;
	uint32_t	user_data_size;
};


typedef struct bitstream_t bitstream_t;
struct bitstream_t 
{
	uint8_t *bitstream;
	int32_t  streamsize;
	int32_t  streambytecnt;
	int32_t  streambitcnt;
//	int32_t  total_bytes_written;
};

typedef struct nalu_t nalu_t;
struct nalu_t
{
//  nalu_types	nal_unit_type;	// nal_unit_type
	uint32_t  	nal_unit_type;	// nal_unit_type
	uint32_t  	temporal_id;	// temporal_id
	uint32_t  	rsvd_zero_bits; // reserved_zero_6bits
	bitstream_t	bs;
};

typedef struct HVENC_Cfg HVENC_Cfg;
struct HVENC_Cfg{
	int32_t  size;
	int32_t  profile;
	int32_t  width, height; // frame size (pixels) 
	float frame_rate;
	int32_t  cu_size;
	int32_t  max_pred_partition_depth;
	int32_t  max_intra_tr_depth;
	int32_t  max_inter_tr_depth;
	int32_t  intra_period;
	int32_t  gop_size; 
	int32_t  num_ref_frames;
	int32_t  motion_estimation_precision;
	int32_t  qp;//for fixed qp mode, or initial qp if cbr or vbr
	int32_t  chroma_qp_offset;
	int32_t  num_enc_engines;
	int32_t  wfpp_enable;
	int32_t  wfpp_num_threads;
	int32_t  sign_hiding;
	int32_t  sample_adaptive_offset;
	int32_t  bitrate_mode;//0=BR_FIXED_QP, 1 = BR_CBR (constant bit rate)
	int32_t  bitrate;//in kbps
	int32_t  vbv_size;//in kbps
	int32_t  vbv_init;//in kbps
	int32_t  reinit_gop_on_scene_change;//
	int32_t  rd_mode;//0=RD_DIST_ONLY,1=RD_FULL,2=RD_FAST
	int32_t  performance_mode;//0=PERF_FULL_COMPUTATION,1=PERF_FAST_COMPUTATION,2=PERF_UFAST_COMPUTATION
};

void *HOMER_enc_init();
void HOMER_enc_close(void* handle);//, nalu_t *nalu_out[], uint32_t   *nalu_list_size)
int32_t  HOMER_enc_encode(void* handle, encoder_in_out_t* input_frame);//uint8_t *picture[]);//, nalu_t *nalu_out[], uint32_t   *nalu_list_size)
int32_t  HOMER_enc_get_coded_frame(void* handle, encoder_in_out_t* output_frame, nalu_t *nalu_out[], uint32_t   *nalu_list_size);
int32_t  HOMER_enc_write_annex_b_output(nalu_t *nalu_out[], uint32_t   num_nalus, encoder_in_out_t *vout);
int32_t  HOMER_enc_control(void *h, int32_t  cmd, void *in);

//void encoder_engine_thread(void* enc_engine);//void *encoder_engine_thread(void *h);//

#ifdef	__cplusplus
}
#endif

#endif  /* __HOMER_HEVC_ENCODER_API_H__ */
