/*****************************************************************************
 * homer_hevc_enc_api.c : homerHEVC encoding library
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

#ifndef __HOMER_HEVC_ENCODER_H__
#define __HOMER_HEVC_ENCODER_H__


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
	NUM_BR_MODES
};

enum HOMER_PERFORMANCE_MODES {
	PERF_FULL_COMPUTATION = 0,
	PERF_FAST_COMPUTATION,
	PERF_UFAST_COMPUTATION,//fastest
	NUM_PERF_MODES
};


#define MAX_STREAMS	8

typedef struct stream_t stream_t;
struct  stream_t
{
	unsigned char	*streams[MAX_STREAMS];
	int				stream_size[MAX_STREAMS];
	int				data_size[MAX_STREAMS];
};


typedef struct encoder_in_out_t encoder_in_out_t;
struct encoder_in_out_t
{
	stream_t stream;
	unsigned int uiPTS;
	unsigned int uiDTS;
	char *pUserdata;
	unsigned int sizeUserdata;
};


typedef struct bitstream_t bitstream_t;
struct bitstream_t 
{
	unsigned char *bitstream;
	int streamsize;
	int streambytecnt;
	int streambitcnt;
//	int total_bytes_written;
};

typedef struct nalu_t nalu_t;
struct nalu_t
{
//  nalu_types	nal_unit_type;	// nal_unit_type
  unsigned int	nal_unit_type;	// nal_unit_type
  unsigned int	temporal_id;	// temporal_id
  unsigned int	rsvd_zero_bits; // reserved_zero_6bits
  bitstream_t	bs;
};

typedef struct HVENC_Cfg HVENC_Cfg;
struct HVENC_Cfg{
	int size;
	int profile;
	int width, height; // frame size (pels) 
	int intra_period;
	int gop_size; 
	int num_b; 
	int frame_rate;
	int num_ref_frames;
	int cu_size;
	int max_pred_partition_depth;
	int max_intra_tr_depth;
	int max_inter_tr_depth;
	int wfpp_enable;
	int wfpp_num_threads;
	int sign_hiding;
	int	rd_mode;//0=RD_DIST_ONLY,1=RD_FULL,2=RD_FAST
	int qp;//for fixed qp mode, or initial qp if cbr or vbr
	int bitrate_mode;//0=BR_FIXED_QP, 1 = BR_CBR (constant bit rate)
	int bitrate;//in kbps
	int vbv_size;//in kbps
	int vbv_init;//in kbps
//	int qp_depth;//
	int chroma_qp_offset;
	int performance_mode;//0=PERF_FULL_COMPUTATION,1=PERF_FAST_COMPUTATION,2=PERF_UFAST_COMPUTATION

};

void *HOMER_enc_init();
void HOMER_enc_close(void* handle);//, nalu_t *nalu_out[], unsigned int *nalu_list_size)
//int HOMER_enc_encode(void* handle, unsigned char *picture[], nalu_t *nalu_out[], unsigned int *nalu_list_size);
int HOMER_enc_encode(void* handle, encoder_in_out_t* input_frame);//unsigned char *picture[]);//, nalu_t *nalu_out[], unsigned int *nalu_list_size)
int HOMER_enc_get_coded_frame(void* handle, encoder_in_out_t* output_frame, nalu_t *nalu_out[], unsigned int *nalu_list_size);
void encoder_thread(void *h);//void encoder_thread(void* ed);
int HOMER_enc_write_annex_b_output(nalu_t *nalu_out[], unsigned int num_nalus, encoder_in_out_t *vout);
int HOMER_enc_control(void *h, int cmd, void *in);


#ifdef	__cplusplus
}
#endif

#endif  /* __HOMER_HEVC_ENCODER_H__ */
