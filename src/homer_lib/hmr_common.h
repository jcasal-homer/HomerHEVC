/*****************************************************************************
 * hmr_common.h : homerHEVC encoding library
/*****************************************************************************
 * Copyright (C) 2014 homerHEVC project
 *
 * Authors: Juan Casal	<jcasal.homer@gmail.com>
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

#ifndef __HOMER_HEVC_COMMON_H__
#define __HOMER_HEVC_COMMON_H__

#include "hmr_private.h"
#include "hmr_container.h"


#define FALSE	0
#define TRUE	1

#ifndef min
#define min(a,b) (a<b?a:b)
#endif // !min

#ifndef max
#define max(a,b) (a>b?a:b)
#endif // !max

#define iswap(a,b) { int aux=a; a=b;b=aux;}
#define ptrswap(type, a,b) { void *_aux_ = (void *)a; a=b;b=(type)_aux_;}

#define clip(value,min,max)	( (value < min) ? min : (value > max) ? max : value )

#define EXCHANGE_PTR(type,a,b) {type* aux;aux=a; a=b; b=aux;}


//nal_unit_type values
//Updated to HM-12.0
typedef enum {
  NALU_CODED_SLICE_TRAIL_N = 0,		// 0
  NALU_CODED_SLICE_TRAIL_R,			// 1
  
  NALU_CODED_SLICE_TSA_N,			// 2
  NALU_CODED_SLICE_TLA_R,				// 3   // Current name in the spec: TSA_R
  
  NALU_CODED_SLICE_STSA_N,			// 4
  NALU_CODED_SLICE_STSA_R,			// 5

  NALU_CODED_SLICE_RADL_N,			// 6
  NALU_CODED_SLICE_RADL_R,				// 7 // Current name in the spec: RADL_R
  
  NALU_CODED_SLICE_RASL_N,			// 8
  NALU_CODED_SLICE_RASL_R,			// 9 // Current name in the spec: RASL_R

  NALU_RESERVED_10,
  NALU_RESERVED_11,
  NALU_RESERVED_12,
  NALU_RESERVED_13,
  NALU_RESERVED_14,
  NALU_RESERVED_15,

  NALU_CODED_SLICE_BLA_W_LP,		// 16   // Current name in the spec: BLA_W_LP
  NALU_CODED_SLICE_BLA_W_RADL,		// 17   // Current name in the spec: BLA_W_DLP
  NALU_CODED_SLICE_BLA_N_LP,		// 18
  NALU_CODED_SLICE_IDR_W_RADL,		// 19  // Current name in the spec: IDR_W_DLP
  NALU_CODED_SLICE_IDR_N_LP,		// 20
  NALU_CODED_SLICE_CRA,				// 21
  NALU_RESERVED_22,
  NALU_RESERVED_23,

  NALU_RESERVED_24,
  NALU_RESERVED_25,
  NALU_RESERVED_26,
  NALU_RESERVED_27,
  NALU_RESERVED_28,
  NALU_RESERVED_29,
  NALU_RESERVED_30,
  NALU_RESERVED_31,

  NALU_TYPE_VPS,							// 32
  NALU_TYPE_SPS,							// 33
  NALU_TYPE_PPS,							// 34
  NALU_TYPE_AUD,							// 35
  NALU_TYPE_EOS,							// 36
  NALU_TYPE_EOB,							// 37
  NALU_TYPE_FILLER_DATA,					// 38
  NALU_TYPE_PREFIX_SEI,						// 39 Prefix SEI
  NALU_TYPE_SUFFIX_SEI,						// 40 Suffix SEI
} nalu_types;

#define rap_pic_flag(nalu_type) (nalu_type >= NALU_CODED_SLICE_BLA_W_LP && nalu_type <= NALU_CODED_SLICE_CRA)//16-21
#define idr_pic_flag(nalu_type) (nalu_type == NALU_CODED_SLICE_IDR_W_RADL || nalu_type == NALU_CODED_SLICE_IDR_N_LP)


//nal_ref_idc values
typedef enum {
 NALU_PRIORITY_LOWEST	   = 0,
 NALU_PRIORITY_LOW         = 1,
 NALU_PRIORITY_HIGH        = 2,
 NALU_PRIORITY_HIGHEST     = 3,
} NalRefIdc;


//hmr_bitstream.c
void hmr_bitstream_alloc(bitstream_t* bs, int size);
void hmr_bitstream_free(bitstream_t* bs);
void hmr_bitstream_init(bitstream_t* bs);
void hmr_bitstream_write_bits(bitstream_t* bs, unsigned int val,int n);
int hmr_bitstream_bitcount(bitstream_t* bs);
void hmr_bitstream_write_bits_uvlc(bitstream_t* bs, unsigned int val);
void hmr_bitstream_write_bits_svlc(bitstream_t* bs, int val);
void hmr_bitstream_rbsp_trailing_bits(bitstream_t* bs);
void hmr_bitstream_put_nal_unit_header(bitstream_t* bs, unsigned int nalu_type, ushort temporal_id, ushort rsvd_zero6bits);
void hmr_bitstream_nalu_ebsp(bitstream_t* in_bs, bitstream_t* out_bs);
void hmr_bitstream_align_bits_1(bitstream_t* bs);
void hmr_bitstream_align_bits_0(bitstream_t* bs);
void hmr_bitstream_write2file(bitstream_t* bs);


//hmr_headers.c
void hmr_put_vps_header(hvenc_t* ed);
void hmr_put_seq_header(hvenc_t* ed);
void hmr_put_pic_header(hvenc_t* ed);
void hmr_put_slice_header(hvenc_t* ed, slice_t *currslice);
void hmr_slice_header_code_wfpp_entry_points(hvenc_t* ed);

//hmr_mem_transfer.c
void *aligned_alloc(int num, int size);
void aligned_free(void *p);
void wnd_alloc(wnd_t *wnd_t, int size_x, int size_y, int offset_x, int offset_y, int pix_size);
void wnd_delete(wnd_t *wnd_t);
void wnd_realloc(wnd_t *wnd_t, int size_x, int size_y, int offset_x, int offset_y, int pix_size);
void wnd_write2file(wnd_t *wnd_t, FILE* file);
void mem_transfer_move_curr_ctu_group(henc_thread_t* et, int i, int j);
void mem_transfer_intra_refs(henc_thread_t* et, ctu_info_t* ctu);
void mem_transfer_decoded_blocks(henc_thread_t* et, ctu_info_t* ctu);
void mem_transfer_1d1d(unsigned char *src, unsigned char *dst, unsigned int width, unsigned int height);
void mem_transfer_1d2d(unsigned char *src, unsigned char *dst, unsigned int width, unsigned int height, unsigned int dst_stride);
void mem_transfer_2d1d(unsigned char *src, unsigned char *dst, unsigned int width, unsigned int height, unsigned int src_stride);
void mem_transfer_2d2d(unsigned char *src, unsigned char *dst, unsigned int width, unsigned int height, unsigned int src_stride, unsigned int dst_stride);

//to access partition list
#define WND_POSITION_1D(type, w, comp, gcnt, ctu_width, offset) ((type)(w).pwnd[comp]+gcnt*ctu_width[comp]+(offset))//((type)(w).pwnd[comp]+((w).data_padding_y[comp])*(w).window_size_x[comp]+(w).data_padding_x[comp]+gcnt*ctu_width[comp]+(offset))
#define WND_POSITION_2D(type, w, comp, part_x, part_y, gcnt, ctu_width) ((type)(w).pwnd[comp]+part_y*(w).window_size_x[comp]+gcnt*ctu_width[comp]+part_x)//((type)(w).pwnd[comp]+((w).data_padding_y[comp]+part_y)*(w).window_size_x[comp]+(w).data_padding_x[comp]+gcnt*ctu_width[comp]+part_x)

//random access
#define WND_DATA_PTR(type, w, comp) ((type)(w).pwnd[comp])//+((w).data_padding_y[comp])*(w).window_size_x[comp]+(w).data_padding_x[comp])
#define WND_STRIDE_2D(w, comp) ((w).window_size_x[comp])

#define CBF(ctu, abs_index, comp, tr_depth) ((ctu->cbf[comp][abs_index]>>(tr_depth))&1)


//init_tables.c
void init_scan_pyramid(hvenc_t* ed, uint* pBuffZ, uint* pBuffH, uint* pBuffV, uint* pBuffD, int iWidth, int iHeight, int iDepth);
void init_flat_quant_pyramids( hvenc_t* ed, uint* quant_pyramid, uint* dequant_pyramid, double* scaling_error_pyramid, uint size, int inv_depth, int qp);
void init_quant_pyramids( hvenc_t* ed, int* quant_pyramid, int* dequant_pyramid, double* scaling_error_pyramid, short* quant_def_table, int width, int height, int ratio, uint sizuNum, uint dc, int inv_depth, int qp);
short* get_default_qtable(int size_mode, int list_index);
void create_abs2raster_tables( unsigned short **zigzag, int total_depth, int depth, int start_value);
void create_raster2abs_tables( unsigned short *zigzag, unsigned short *inv_zigzag, int max_cu_width, int max_cu_height, int total_depth);
void hmr_rd_init(hvenc_t* ed, slice_t *currslice);
int find_scan_mode(int is_intra, int is_luma, int width, int dir_mode, int up_left_luma_dir_mode);


//encodelib.cpp
void copy_ctu(ctu_info_t* src_ctu, ctu_info_t* dst_ctu);
THREAD_RETURN_TYPE encoder_thread(void *h);//void encoder_thread(void* ed);

//hmr_motion_intra.c
#define ADI_POINTER_MIDDLE(ptr_adi_orig, size)  (ptr_adi_orig+(size>>1))	//points to the top left square
void init_partition_info(henc_thread_t* et, cu_partition_info_t *partition_list);
void create_partition_ctu_neighbours(henc_thread_t* et, ctu_info_t *ctu, cu_partition_info_t* curr_partition_info);
uint motion_intra(henc_thread_t* et, ctu_info_t* ctu, int gcnt);
//void cu_partition_get_neighbours(cu_partition_info_t *curr_part, int cu_size);
void cu_partition_get_neighbours(cu_partition_info_t *curr_part, int cu_size, int ctu_valid_colums, int ctu_valid_lines);
uint32_t sad(uint8_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride, int size);
uint32_t ssd(uint8_t * src, uint32_t src_stride, int16_t * pred, uint32_t pred_stride, int size);
uint32_t sad_16(int16_t *src, uint32_t src_stride, int16_t *pred, uint32_t pred_stride, int size);
uint32_t ssd_16(int16_t *src, uint32_t src_stride, int16_t *pred, uint32_t pred_stride, int size);
uint32_t modified_variance(uint8_t *p, int size, int stride, int modif);
void fill_reference_samples(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* partition_info, int adi_size, int16_t* decoded_buff, int decoded_buff_stride, int partition_size, int component, int is_filtered);
void create_intra_angular_prediction(henc_thread_t* et, ctu_info_t* ctu, int16_t* pred, int pred_stride, short  *adi_pred_buff, int adi_size, int cu_size, int cu_mode, int is_luma);
void create_intra_planar_prediction(henc_thread_t* et, int16_t* pred, int pred_stride, int16_t *adi_pred_buff, int adi_size, int cu_size, int cu_size_shift);
void synchronize_reference_buffs_luma(henc_thread_t* et, cu_partition_info_t* curr_part, wnd_t *decoded_src, wnd_t * decoded_dst, int gcnt);
void create_intra_recursive_stop_info(henc_thread_t* et, int gcnt);
void predict(uint8_t * orig_auxptr, int orig_buff_stride, int16_t *pred_auxptr, int pred_buff_stride, int16_t *residual_auxptr, int residual_buff_stride, int curr_part_size);
void reconst(int16_t *pred_auxptr, int pred_buff_stride, int16_t *residual_auxptr, int residual_buff_stride, int16_t *decoded_auxptr, int decoded_buff_stride, int curr_part_size);
void synchronize_motion_buffers_luma(henc_thread_t* et, cu_partition_info_t* curr_part, wnd_t * quant_src, wnd_t * quant_dst, wnd_t *decoded_src, wnd_t * decoded_dst, int gcnt);
//void consolidate_intra_prediction_info(henc_thread_t *et, ctu_info_t *ctu, ctu_info_t *ctu_rd, cu_partition_info_t *parent_part_info, int parent_cost, int children_cost, int is_max_depth, uint *cost_sum);
uint encode_intra_luma(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type);
//uint encode_intra_chroma(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position,  int part_size_type);
void synchronize_reference_buffs(henc_thread_t* et, cu_partition_info_t* curr_part, wnd_t *decoded_src, wnd_t * decoded_dst, int gcnt);
uint encode_intra(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int curr_depth, int position, PartSize part_size_type);
void analyse_recursive_info(henc_thread_t* et, ctu_info_t* ctu, int gcnt);
void consolidate_info_buffers_for_rd(henc_thread_t* et, ctu_info_t* ctu, int dest_depth, int abs_index, int num_part_in_cu);
void synchronize_cu_wnd(henc_thread_t* et, cu_partition_info_t* curr_part, wnd_t * wnd_src, wnd_t * wnd_dst);


//hmr_motion_inter.c
uint motion_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt);
void hmr_motion_inter_uni(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint8_t *orig_buff, int orig_buff_stride, uint8_t *reference_buff, int reference_buff_stride, uint8_t *pred_buff, int pred_buff_stride,  int curr_part_global_x, int curr_part_global_y, int curr_part_size, int curr_part_size_shift, motion_vector_t *mv);
void hmr_interpolate_chroma(int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int fraction, int width, int height, int is_vertical, int is_first, int is_last);
//void consolidate_prediction_info(henc_thread_t *et, ctu_info_t *ctu, ctu_info_t *ctu_rd, cu_partition_info_t *parent_part_info, int parent_cost, int children_cost, int is_max_depth, uint *cost_sum);
void consolidate_prediction_info(henc_thread_t *et, ctu_info_t *ctu, ctu_info_t *ctu_rd, cu_partition_info_t *parent_part_info, uint parent_cost, uint children_cost, int is_max_depth, uint *cost_sum);

//hmr_motion_intra_chroma.c
void create_chroma_dir_list(int* list, int luma_mode);
void synchronize_motion_buffers_chroma(henc_thread_t* et, cu_partition_info_t* curr_part, wnd_t * quant_src, wnd_t * quant_dst, wnd_t *decoded_src, wnd_t * decoded_dst, int gcnt);
void synchronize_reference_buffs_chroma(henc_thread_t* et, cu_partition_info_t* curr_part, wnd_t *decoded_src, wnd_t * decoded_dst, int gcnt);
uint encode_intra_chroma(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position,  int part_size_type);
void homer_update_cand_list( uint uiMode, double Cost, uint uiFastCandNum, int CandModeList[3], double CandCostList[3] );

//hmr_motion_inter_chroma.c
int encode_inter_chroma(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type);

//hmr_transform.c
void transform(int bitDepth, int16_t *block,int16_t *coeff, int block_size, int iWidth, int iHeight, int width_shift, int height_shift, unsigned short uiMode, int16_t *aux);
void itransform(int bitDepth, short *block,short *coeff, int block_size, int iWidth, int iHeight, uint uiMode, short *aux);

//hmr_quant.c
void sign_bit_hidding( short * dst, short * src, uint const *scan, short* deltaU, int width, int height );
void quant(henc_thread_t* et, short * src, short * dst, int scan_mode, int depth, int comp, int cu_mode, int is_intra, int *ac_sum, int cu_size, int per, int rem);
void iquant(henc_thread_t* et, short * src, short * dst, int depth, int comp, int is_intra, int cu_size, int per, int rem);


//hmr_deblocking_filter.c
void hmr_deblock_filter(hvenc_t* ed, slice_t *currslice);

//hmr_arithmetic_encoding.c
void ee_init_contexts(enc_env_t *ee);
void ee_start_entropy_model(enc_env_t *ee, int slice_type, int qp, int cabac_init_flag);
void ee_copy_entropy_model(enc_env_t *ee_src, enc_env_t *ee_dst);
void ee_encode_ctu(henc_thread_t* et, enc_env_t* ee, slice_t *currslice, ctu_info_t* cu, int gcnt);
void ee_encode_coding_unit(henc_thread_t* et, enc_env_t* ee, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, int gcnt);
void ee_end_slice(enc_env_t* ee, slice_t *currslice, ctu_info_t* ctu);
int get_intra_dir_luma_predictor(ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, int* uiIntraDirPred, int* piMode);
uint fast_rd_estimate_bits_intra_luma_mode( henc_thread_t* et, cu_partition_info_t* partition_info, uint pred_depth, int dir, int *preds, int num_preds);
uint rd_estimate_bits_intra_mode( henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* partition_info, uint pred_depth, int is_luma);
uint rd_get_intra_bits_qt( henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* partition_info, uint pred_depth, int is_luma, int gcnt);
ctu_info_t *get_pu_left(ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint *aux_part_idx);
ctu_info_t *get_pu_left_bottom(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint *aux_part_idx);
ctu_info_t *get_pu_top(ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint *aux_part_idx, int planarAtLCUBoundary);
ctu_info_t *get_pu_top_right(ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint *aux_part_idx);
ctu_info_t *get_pu_top_left(ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint *aux_part_idx);

//hmr_binary_encoding.c //bm = binary model, be = bienary encoder, bc = binary counter
void bm_copy_binary_model(binary_model_t *bm_src, binary_model_t *bm_dst);
void bm_map_funcs(enc_env_t* ee);
void bc_init_next_state_table();//() for fast-rd


//hmr_metrics.c 
void homer_psnr(picture_t *orig, wnd_t* decoded, int pic_width[3], int pic_height[3], double psnr[3]);



//hmr_rate_control.c
void hmr_rc_init(hvenc_t* ed);
void hmr_rc_init_seq(hvenc_t* ed);
void hmr_rc_gop(hvenc_t* ed);//, int np, int nb)
void hmr_rc_init_pic(hvenc_t* ed, slice_t *currslice);
void hmr_rc_end_pic(hvenc_t* ed, slice_t *currslice);
int hmr_rc_get_cu_qp(henc_thread_t* et, ctu_info_t *ctu, cu_partition_info_t *curr_cu_info, slice_t *currslice);
void hmr_rc_change_pic_mode(henc_thread_t* ed, slice_t *currslice);
#endif //__HOMER_HEVC_COMMON_H__
