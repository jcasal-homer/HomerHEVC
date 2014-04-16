/*****************************************************************************
 * hmr_private.h : homerHEVC encoding library
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


#ifndef __HOMER_HEVC_PRIVATE_H__
#define __HOMER_HEVC_PRIVATE_H__

#include <stdio.h>
#include "hmr_os_primitives.h"
#include "homer_hevc_enc_api.h"


#define COMPUTE_SSE_FUNCS		1
#define COMPUTE_AS_HM			1
#define COMPUTE_METRICS			1

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		This is the herarchy to follow

		Headers									coding layer
	---------------								------------
	video parameter set								Gop 
		|											 |
	Sequence parameter set							Picture
		|											 |
	picture_t parameter set							Slice
		|
	 slice header
********************************************************************/

#define DOUBLE_MAX              (1.7e+308)    

#define MAX_NUM_CTUs				8160	//1920*1088 resolution
#define MAX_MB_GROUP_SIZE		1

#define     MAX_CU_DEPTHS			7		//from 128x128 to 2x2. For scan piramid. 
//#define     MAX_CU_SIZE					(1<<(MAX_CU_DEPTHS))        // maximum allowable size of CU
#define     MAX_CU_SIZE_SHIFT		6								//MAX_CU_DEPTH	
#define     MAX_CU_SIZE			(1<<(MAX_CU_SIZE_SHIFT))	// current maximum allowable size of CU

#define		MAX_PARTITION_DEPTH		5 //de 64x64 a 4x4 ambos inclusive

#define		MIN_TU_SIZE_SHIFT		2		//min size of Transform is 4x4
#define		MAX_TU_SIZE_SHIFT		5		//min size of Transform is 32x32

#define     MAX_NUM_PARTITIONS		256							//(1<<MAX_CU_PARTITIONS_SHIFT)*(1<<MAX_CU_PARTITIONS_SHIFT) - 16 particiones por eje - se corresponde con el peor caso

#define QUANT_DEFAULT_DC	16

//for intra prediction
#define NUM_INTRA_MODES					35
#define NUM_INTRA_CANDIDATES			3


//scannig modes
#define ZIGZAG_SCAN							0     
#define HOR_SCAN							1
#define VER_SCAN							2
#define DIAG_SCAN							3
#define NUM_SCAN_MODES						4

//scaling list modes
#define SCALING_MODE_4x4					0
#define SCALING_MODE_8x8					1
#define SCALING_MODE_16x16					2
#define SCALING_MODE_32x32					3			
#define NUM_SCALING_MODES					4


#define NUM_SCALING_LISTS					6         //list number for quantization matrix
#define NUM_SCALING_REM_LISTS				6         //remainder of QP/6
#define NUM_MAX_MATRIX_SIZE					8		  //max size number for quantization matrix	//MAX_MATRIX_SIZE_NUM en HM

//from HM
#define QUANT_IQUANT_SHIFT					20 // Q(QP%6) * IQ(QP%6) = 2^20
#define QUANT_SHIFT							14 // Q(4) = 2^14
#define SCALE_BITS							15 // Inherited from TMuC, pressumably for fractional bit estimates in RDOQ
#define MAX_TR_DYNAMIC_RANGE				15 // Maximum transform dynamic range (excluding sign bit)



//-------------------------init binary encode-------------------------------------
#define NUM_SPLIT_FLAG_CTX            3       ///< number of context models for split flag
#define NUM_SKIP_FLAG_CTX             3       ///< number of context models for skip flag

#define NUM_MERGE_FLAG_EXT_CTX        1       ///< number of context models for merge flag of merge extended
#define NUM_MERGE_IDX_EXT_CTX         1       ///< number of context models for merge index of merge extended

#define NUM_PART_SIZE_CTX             4       ///< number of context models for partition size
#define NUM_CU_AMP_CTX                1       ///< number of context models for partition size (AMP)
#define NUM_PRED_MODE_CTX             1       ///< number of context models for prediction mode

#define NUM_ADI_CTX                   1       ///< number of context models for intra prediction

#define NUM_CHROMA_PRED_CTX           2       ///< number of context models for intra prediction (chroma)
#define NUM_INTER_DIR_CTX             5       ///< number of context models for inter prediction direction
#define NUM_MV_RES_CTX                2       ///< number of context models for motion vector difference

#define NUM_REF_NO_CTX                2       ///< number of context models for reference index
#define NUM_TRANS_SUBDIV_FLAG_CTX     3       ///< number of context models for transform subdivision flags
#define NUM_QT_CBF_CTX                5       ///< number of context models for QT CBF
#define NUM_QT_ROOT_CBF_CTX           1       ///< number of context models for QT ROOT CBF
#define NUM_DELTA_QP_CTX              3       ///< number of context models for dQP

#define NUM_SIG_CG_FLAG_CTX           2       ///< number of context models for MULTI_LEVEL_SIGNIFICANCE

#define NUM_SIG_FLAG_CTX              42      ///< number of context models for sig flag
#define NUM_SIG_FLAG_CTX_LUMA         27      ///< number of context models for luma sig flag
#define NUM_SIG_FLAG_CTX_CHROMA       15      ///< number of context models for chroma sig flag

#define NUM_CTX_LAST_FLAG_XY          15      ///< number of context models for last coefficient position

#define NUM_ONE_FLAG_CTX              24      ///< number of context models for greater than 1 flag
#define NUM_ONE_FLAG_CTX_LUMA         16      ///< number of context models for greater than 1 flag of luma
#define NUM_ONE_FLAG_CTX_CHROMA        8      ///< number of context models for greater than 1 flag of chroma
#define NUM_ABS_FLAG_CTX               6      ///< number of context models for greater than 2 flag
#define NUM_ABS_FLAG_CTX_LUMA          4      ///< number of context models for greater than 2 flag of luma
#define NUM_ABS_FLAG_CTX_CHROMA        2      ///< number of context models for greater than 2 flag of chroma

#define NUM_MVP_IDX_CTX               2       ///< number of context models for MVP index

#define NUM_SAO_MERGE_FLAG_CTX        1       ///< number of context models for SAO merge flags
#define NUM_SAO_TYPE_IDX_CTX          1       ///< number of context models for SAO type index

#define NUM_TRANSFORMSKIP_FLAG_CTX    1       ///< number of context models for transform skipping 
#define NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX  1 
#define CNU                          154      ///< dummy initialization value for unused context models 'Context model Not Used'


#define NUM_CTXs	(NUM_SPLIT_FLAG_CTX+NUM_SKIP_FLAG_CTX+NUM_MERGE_FLAG_EXT_CTX+NUM_MERGE_IDX_EXT_CTX	\
					+NUM_PART_SIZE_CTX+NUM_PRED_MODE_CTX+NUM_ADI_CTX+NUM_CHROMA_PRED_CTX+NUM_DELTA_QP_CTX	\
					+NUM_INTER_DIR_CTX+NUM_REF_NO_CTX+NUM_MV_RES_CTX+2*NUM_QT_CBF_CTX+NUM_TRANS_SUBDIV_FLAG_CTX	\
					+NUM_QT_ROOT_CBF_CTX+2*NUM_SIG_CG_FLAG_CTX+NUM_SIG_FLAG_CTX+2*NUM_CTX_LAST_FLAG_XY	\
					+2*NUM_CTX_LAST_FLAG_XY+NUM_ONE_FLAG_CTX+NUM_ABS_FLAG_CTX+NUM_MVP_IDX_CTX	\
					+NUM_CU_AMP_CTX+NUM_SAO_MERGE_FLAG_CTX+NUM_SAO_TYPE_IDX_CTX	\
					+2*NUM_TRANSFORMSKIP_FLAG_CTX+NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX)
//-------------------------end binary encode-------------------------------------


#define MAXnum_ref_frames_in_pic_order_cnt_cycle  256	
typedef int Boolean;
typedef unsigned char	byte;
typedef unsigned short	ushort;
typedef unsigned int	uint;
typedef long			int64;
typedef unsigned long	uint64;

#define SIGN(x)					(x<0?-1:1)
#define QP_BITS                 15

//quantified residual access macros
#define FORMAT_PUT_SIGN(s)	(s<<8)
#define FORMAT_GET_SIGN(s)	((s&0x0000ff00)>>8)
#define FORMAT_PUT_RUN	
#define FORMAT_GET_RUN(s)	(s&0x000000ff)
#define FORMAT_PUT_VAL(s)	(s<<16)	
#define FORMAT_GET_VAL(s)	(s>>16)


//Intra modes 
#define PLANAR_IDX				0
#define VER_IDX					26                    // index for intra VERTICAL   mode
#define HOR_IDX					10                    // index for intra HORIZONTAL mode
#define DC_IDX					1                     // index for intra DC mode
#define NUM_CHROMA_MODE			5                     // total number of chroma modes
#define NUM_LUMA_MODES			35                     // total number of chroma modes
#define DM_CHROMA_IDX			36                    // chroma mode index for derived from luma intra mode

#define REG_DCT					65535


#define INTER_MODE				0           
#define INTRA_MODE				1	


// partition types
typedef enum 
{
  SIZE_2Nx2N,           // symmetric motion partition,  2Nx2N
  SIZE_2NxN,            // symmetric motion partition,  2Nx N
  SIZE_Nx2N,            // symmetric motion partition,   Nx2N
  SIZE_NxN,             // symmetric motion partition,   Nx N
  SIZE_2NxnU,           // asymmetric motion partition, 2Nx( N/2) + 2Nx(3N/2)
  SIZE_2NxnD,           // asymmetric motion partition, 2Nx(3N/2) + 2Nx( N/2)
  SIZE_nLx2N,           // asymmetric motion partition, ( N/2)x2N + (3N/2)x2N
  SIZE_nRx2N,           // asymmetric motion partition, (3N/2)x2N + ( N/2)x2N
  SIZE_NONE = 15
}PartSize;


typedef enum 
{
  REF_PIC_LIST_0 = 0,   // reference list 0
  REF_PIC_LIST_1 = 1,   // reference list 1
  REF_PIC_LIST_C = 2,   // combined reference list for uni-prediction in B-Slices
  REF_PIC_LIST_X = 100  // special mark
}RefPicList;


typedef enum {
  Y_COMP=0,
  U_COMP,
  V_COMP,
  CHR_COMP = 1,//chroma
  NUM_PICT_COMPONENTS = 3
} PictComponents;
//#define NUM_PICT_COMPONENTS			3				//Y,U,V,PTS,private data, size of private data
//#define NUM_PICT_PARAMS		(NUM_PICT_COMPONENTS+3)	//Y,U,V,PTS,private data, size of private data


#define B_BITS         10    // Number of bits to represent the whole coding interval
#define BITS_TO_LOAD   16
#define HALF           0x01FE      //(1 << (B_BITS-1)) - 2


typedef enum {
  P_SLICE,
  B_SLICE,
  I_SLICE,
  NUM_SLICE_TYPES
} SliceTypes;

#define isIntra(slice_type) (slice_type == I_SLICE)



//AVC Profile IDC definitions
typedef enum {
  FREXT_CAVLC444 = 44,       //!< YUV 4:4:4/14 "CAVLC 4:4:4"
  BASELINE       = 66,       //!< YUV 4:2:0/8  "Baseline"
  MAIN           = 77,       //!< YUV 4:2:0/8  "Main"
  EXTENDED       = 88,       //!< YUV 4:2:0/8  "Extended"
  FREXT_HP       = 100,      //!< YUV 4:2:0/8  "High"
  FREXT_Hi10P    = 110,      //!< YUV 4:2:0/10 "High 10"
  FREXT_Hi422    = 122,      //!< YUV 4:2:2/10 "High 4:2:2"
  FREXT_Hi444    = 244,      //!< YUV 4:4:4/14 "High 4:4:4"
  MVC_HIGH       = 118,      //!< YUV 4:2:0/8  "Multiview High"
  STEREO_HIGH    = 128       //!< YUV 4:2:0/8  "Stereo High"
} ProfileIDC;


//  Available MB modes
typedef enum {
  PSKIP        =  0,
  BSKIP_DIRECT =  0,
  P16x16       =  1,
  P16x8        =  2,
  P8x16        =  3,
  SMB8x8       =  4,
  SMB8x4       =  5,
  SMB4x8       =  6,
  SMB4x4       =  7,
  P8x8         =  8,
  I4MB         =  9,
  I16MB        = 10,
  IBLOCK       = 11,
  SI4MB        = 12,
  I8MB         = 13,
  IPCM         = 14,
  MAXMODE      = 15
} MBModes;

typedef enum {
  VERT_PRED            = 0,
  HOR_PRED             = 1,
  DC_PRED              = 2,
  DIAG_DOWN_LEFT_PRED  = 3,
  DIAG_DOWN_RIGHT_PRED = 4,
  VERT_RIGHT_PRED      = 5,
  HOR_DOWN_PRED        = 6,
  VERT_LEFT_PRED       = 7,
  HOR_UP_PRED          = 8
} I4x4PredModes;

// 16x16 intra prediction modes
typedef enum {
  VERT_PRED_16   = 0,
  HOR_PRED_16    = 1,
  DC_PRED_16     = 2,
  PLANE_16       = 3
} I16x16PredModes;

// 8x8 chroma intra prediction modes
typedef enum {
  DC_PRED_8     =  0,
  HOR_PRED_8    =  1,
  VERT_PRED_8   =  2,
  PLANE_8       =  3
} I8x8PredModes;

//! definition of H.264 syntax elements
typedef enum
{
  SE_HEADER,
  SE_PTYPE,
  SE_MBTYPE,
  SE_REFFRAME,
  SE_INTRAPREDMODE,
  SE_MVD,
  SE_CBP,
  SE_LUM_DC_INTRA,
  SE_CHR_DC_INTRA,
  SE_LUM_AC_INTRA,
  SE_CHR_AC_INTRA,
  SE_LUM_DC_INTER,
  SE_CHR_DC_INTER,
  SE_LUM_AC_INTER,
  SE_CHR_AC_INTER,
  SE_DELTA_QUANT,
  SE_BFRAME,
  SE_EOS,
  SE_MAX_ELEMENTS = 20 //!< number of maximum syntax elements
} SE_type;             // substituting the definitions in elements.h


#define IS_INTRA(mode)    (mode==SI4MB || mode==I4MB || mode==I16MB || mode==I8MB || mode==IPCM)

//progressive scan pattern 
static const byte scan_progressive[16][2] =
{
  {0,0},{1,0},{0,1},{0,2},
  {1,1},{2,0},{3,0},{2,1},
  {1,2},{0,3},{1,3},{2,2},
  {3,1},{3,2},{2,3},{3,3}
};

//! interlaced scan pattern
static const byte scan_interlaced[16][2] =
{
  {0,0},{0,1},{1,0},{0,2},
  {0,3},{1,1},{1,2},{1,3},
  {2,0},{2,1},{2,2},{2,3},
  {3,0},{3,1},{3,2},{3,3}
};


/*static const byte ZZ_SCAN[16]  =
{  0,  1,  4,  8,  5,  2,  3,  6,  9, 12, 13, 10,  7, 11, 14, 15
};

static const byte ZZ_SCAN8[64] =
{  0,  1,  8, 16,  9,  2,  3, 10, 17, 24, 32, 25, 18, 11,  4,  5,
   12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13,  6,  7, 14, 21, 28,
   35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
   58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63
};
*/
static __inline short I16Offset (int cbp, short i16mode)
{
  return (short) ((cbp & 15 ? 13 : 1) + i16mode + ((cbp & 0x30) >> 2));
}

static __inline short clip_byte(short x)
{
//	return x&(~255) ? (-x)>>31 : x;
	return ( (x < 0) ? 0 : (x > 255) ? 255 : x );
}

static __inline int is_FREXT_profile(unsigned int profile_idc) 
{
  return ( profile_idc>=FREXT_HP || profile_idc == FREXT_CAVLC444 );
}



#define MAX_TLAYER                  8           // max number of sublayers (temporal layers en el HM)
typedef struct profile_tier_t   profile_tier_t;

struct profile_tier_t
{
  unsigned int		profileSpace;
  unsigned int		tierFlag;
  unsigned int		profileIdc;

  unsigned int		profileCompatibilityFlag[32];

  unsigned int		general_progressive_source_flag;
  unsigned int		general_interlaced_source_flag;
  unsigned int		general_non_packed_constraint_flag;
  unsigned int		general_frame_only_constraint_flag;
  
  unsigned int		levelIdc;
};

typedef struct profile_tier_level_t profile_tier_level_t;

struct profile_tier_level_t
{
  profile_tier_t	generalPTL;
  profile_tier_t	subLayerPTL[6];      // max. value of max_sub_layers_minus1 is 6
  unsigned int		subLayerProfilePresentFlag[6];
  unsigned int		subLayerLevelPresentFlag[6];	
};

typedef struct vps_t vps_t;

struct vps_t
{
	unsigned int video_parameter_set_id;			//u(4)
	unsigned int temporal_id_nesting_flag;			//u(1)
	//reserved_zero									//u(2) 
	//reserved_zero									//u(6)
//	unsigned int	max_sub_layers_minus1;			//u(3)

	profile_tier_level_t *ptl;

	unsigned int sub_layer_ordering_info_present_flag;//u(1)
	//profile and level
	//reserved_zero									//u(12)
	unsigned int max_dec_pic_buffering[MAX_TLAYER];	//ue(v)
	unsigned int max_num_reorder_pics[MAX_TLAYER];	//ue(v)
	unsigned int max_latency_increase[MAX_TLAYER];	//ue(v)

	//unsigned int vps_max_layer_id;				//u(6)
	unsigned int num_layer_sets_minus1;			//ue(v) //vps_max_op_sets_minus1

	unsigned int timing_info_present_flag;		//u(1)

	unsigned int num_hdr_parameters;				//u(3)

};

typedef struct sps_t sps_t;
struct sps_t
{
	Boolean   Valid;                  // indicates the parameter set is valid

	unsigned int video_parameter_set_id;			//u(4)
//	unsigned int sps_max_sub_layers_minus1;         //u(3)
	//reserved_zero									//u(1) 

	profile_tier_level_t *ptl;						

	unsigned int seq_parameter_set_id;				//ue(v)
	unsigned int chroma_format_idc;					//ue(v)
	//separate_colour_plane_flag					//u(1) 
	unsigned int pic_width_in_luma_samples;			//ue(v)
	unsigned int pic_height_in_luma_samples;		//ue(v)
	unsigned int conformance_window_flag;					//u(1)
	unsigned int conf_win_left_offset;				//ue(v)
	unsigned int conf_win_right_offset;				//ue(v)
	unsigned int conf_win_top_offset;				//ue(v)
	unsigned int conf_win_bottom_offset;			//ue(v)
	unsigned int bit_depth_luma_minus8;				//ue(v)
	unsigned int bit_depth_chroma_minus8;			//ue(v)
	unsigned int pcm_enabled_flag;					//u(1)
	unsigned int pcm_bit_depth_luma_minus1;			//u(4)
	unsigned int pcm_bit_depth_chroma_minus1;		//u(4)
	unsigned int log2_max_pic_order_cnt_lsb_minus4;	//ue(v)
	unsigned int max_dec_pic_buffering[MAX_TLAYER];	//ue(v)
	unsigned int max_num_reorder_pics[MAX_TLAYER];	//ue(v)
	unsigned int max_latency_increase[MAX_TLAYER];	//ue(v)
	unsigned int restricted_ref_pic_lists_flag;		//u(1)
	unsigned int lists_modification_present_flag;	//u(1)
	unsigned int log2_min_coding_block_size_minus3; //ue(v)
	unsigned int log2_diff_max_min_coding_block_size; //ue(v)
	unsigned int log2_min_transform_block_size_minus2;//ue(v)
	unsigned int log2_diff_max_min_transform_block_size;//ue(v)
	unsigned int log2_min_pcm_coding_block_size_minus3;//ue(v)
	unsigned int log2_diff_max_min_pcm_coding_block_size;//ue(v)
	unsigned int max_transform_hierarchy_depth_inter;//ue(v)
	unsigned int max_transform_hierarchy_depth_intra;//ue(v)
	unsigned int scaling_list_enabled_flag;			//u(1)
	unsigned int scaling_list_data_present_flag;	//u(1)
	//scaling_list_data()
	unsigned int amp_enabled_flag;					//u(1)	//amp = asymetric motion partitions
	unsigned int sample_adaptive_offset_enabled_flag;//u(1)
	unsigned int pcm_loop_filter_disable_flag;		//u(1)
	unsigned int temporal_id_nesting_flag;			//u(1)
//	unsigned int num_short_term_ref_pic_sets;		//ue(v) - moved to ed
	//hmr_short_term_ref_pic_set(i)								- moved to ed
//	unsigned int long_term_ref_pics_present_flag;	//u(1)	- moved to ed
//	unsigned int num_long_term_ref_pic_sps;			//ue(v)	- moved to ed
//	unsigned int lt_ref_pic_poc_lsb_sps;			//u(v)
	unsigned int used_by_curr_pic_lt_sps_flag;		//u(1)
	unsigned int temporal_mvp_enable_flag;			//u(1)	//mvp = motion vector prediction
	unsigned int strong_intra_smooth_enabled_flag;	//u(1)
	unsigned int vui_parameters_present_flag;		//u(1)
		//vui_parameters()
	//unsigned int sps_extension_flag;				//u(1)
};


typedef struct pps_t pps_t;
struct pps_t
{
//  Boolean   Valid;                  // indicates the parameter set is valid
	unsigned int pic_parameter_set_id;                          // ue(v)
	unsigned int seq_parameter_set_id;                          // ue(v)
	unsigned int dependent_slice_enabled_flag;					// u(1)
	unsigned int output_flag_present_flag;                      // u(3)
	unsigned int num_extra_slice_header_bits;					// u(1)
	unsigned int sign_data_hiding_flag;							// u(1)
	unsigned int cabac_init_present_flag;                       // u(1)
//  unsigned int num_ref_idx_l0_default_active_minus1;          // ue(v)
//  unsigned int num_ref_idx_l1_default_active_minus1;          // ue(v)

	unsigned int pic_init_qp_minus26;							// se(v)

	unsigned int constrained_intra_pred_flag;					// u(1)
	unsigned int transform_skip_enabled_flag;					// u(1)
	unsigned int cu_qp_delta_enabled_flag;						// u(1)
	unsigned int diff_cu_qp_delta_depth;						// ue(v)
	unsigned int pic_cb_qp_offset;								// se(v)
	unsigned int pic_cr_qp_offset;								// se(v)

	unsigned int pic_slice_level_chroma_qp_offsets_present_flag;// u(1)
	unsigned int weighted_pred_flag;							// u(1)
	unsigned int weighted_bipred_flag;							// u(1)
	unsigned int transquant_bypass_enable_flag;					// u(1)
	unsigned int tiles_enabled_flag;							// u(1)
	unsigned int entropy_coded_sync_enabled_flag;				// u(1)
//	unsigned int entropy_slice_enabled_flag;					// u(1) - Removed

//	unsigned int num_tile_columns_minus1;						// ue(v)
//	unsigned int num_tile_rows_minus1;							// ue(v)
//	unsigned int uniform_spacing_flag;							// u(1)
//	...........................
	unsigned int loop_filter_across_slices_enabled_flag;		// u(1)
	unsigned int deblocking_filter_control_present_flag;		// u(1)
	//.......................
	unsigned int pps_scaling_list_data_present_flag;			// u(1)
	unsigned int lists_modification_present_flag;				// u(1)
	unsigned int log2_parallel_merge_level_minus2;				// u(1)
	unsigned int slice_header_extension_present_flag;			// u(1)
};


typedef struct wnd_t wnd_t;
struct wnd_t
{
	void	*pwnd[3];
	int		window_size_x[3];
	int		window_size_y[3];
	int		data_offset_x;//indica el desplazamiento de la ventana respecto a los macrobloques a predecir
	int		data_offset_y;//indica el desplazamiento de la ventana respecto a los macrobloques a predecir
	int		window_orig_x;//indica la coordenada x del inicio de la ventana de busqueda en la imagen origen
	int		window_orig_y;//indica la coordenada x del inicio de la ventana de busqueda en la imagen origen
	int		pix_size;
#ifdef WRITE_REF_FRAMES
	FILE	*out_file;//for debug porposes
#endif
};

typedef struct video_frame_t video_frame_t;
struct video_frame_t
{
	wnd_t			img;
	unsigned int	pts;
};

//picture_t pool for image storage
/*#define PICT_POOL_SIZE 5
typedef struct img_pool
{
	unsigned char	*imgs[PICT_POOL_SIZE][NUM_PICT_PARAMS];
	unsigned int	write_index;
	unsigned int	read_index;
}img_pool_t;

typedef struct symbol
{
	int type;
	int value;
	int value2;
	int context;
	int len;
}symbol_t;
*/
/*
typedef struct bit_counter
{
  int mb_total;
  unsigned short mb_mode;
  unsigned short mb_inter;
  unsigned short mb_cbp;
  unsigned short mb_delta_quant;
  int mb_y_coeff;
  int mb_uv_coeff;
  int mb_cb_coeff;
  int mb_cr_coeff;  
  int mb_stuffing;
}bit_counter_t;
*/

typedef struct cu_partition_info_t cu_partition_info_t;
struct cu_partition_info_t
{
	ushort list_index;
	ushort depth;
	ushort abs_index;
	ushort size;
	ushort size_chroma;
	ushort x_position,x_position_chroma;
	ushort y_position,y_position_chroma;
	ushort num_part_in_cu;//tama�o en particiones minimas

	ushort raster_index;
	ushort left_neighbour;
	ushort top_neighbour; 
	ushort left_bottom_neighbour;
	ushort top_right_neighbour; 
	ushort is_tl_inside_frame, is_b_inside_frame, is_r_inside_frame;
	ushort abs_index_left_partition;
	ushort abs_index_top_partition;

	cu_partition_info_t	*parent;//pointer to parent partition
	cu_partition_info_t	*children[4];//pointers to child partitions
	
	ushort mode;
	ushort mode_chroma;

	uint sum;
	uint distortion_chroma, cost_chroma;
	uint distortion, cost;
	uint variance, variance_chroma;
	uint recursive_split;
};

typedef struct ctu_info_t ctu_info_t ;
struct ctu_info_t 
{
	int				ctu_number;
	int				x[3],y[3];//global in picture_t
	int				size;
	int				num_part_in_ctu;

	cu_partition_info_t	*partition_list;

	byte			*cbf[NUM_PICT_COMPONENTS];//[MAX_NUM_PARTITIONS];
	byte			*intra_mode[NUM_PICT_COMPONENTS-1];//[MAX_NUM_PARTITIONS];
	byte			*tr_idx;//[MAX_NUM_PARTITIONS];
	byte			*pred_depth;//[MAX_NUM_PARTITIONS];
	byte			*part_size_type;//[MAX_NUM_PARTITIONS];
	//byte			*pred_mode;//[MAX_NUM_PARTITIONS];//intra or inter
	int				pred_mode;
	wnd_t			*coeff_wnd;

	ctu_info_t		*ctu_left;
	ctu_info_t		*ctu_left_bottom;
	ctu_info_t		*ctu_top;
	ctu_info_t		*ctu_top_right;
	int				top;
	int				left;
	
	//quant
	int				qp/*, prev_qp, prev_dqp*/;
    int				per;
    int				rem;
//    int				qpbits;
//	uint			variance;
};




#define FAST_BIT_EST	1

//TEncBinCABAC
typedef struct binary_model_t binary_model_t;
struct binary_model_t
{
	uint	m_uiLow;
	uint	m_uiRange;
	uint	m_bufferedByte;
	int		m_numBufferedBytes;
	uint	m_bitsLeft;
	uint	m_uiBinsCoded;
	int		m_binCountIncrement;
#if FAST_BIT_EST
	uint64	m_fracBits;
#endif
};



//ContextModel
typedef struct context_model_t context_model_t;
struct context_model_t
{
	uint		num_bins_coded;
	byte		state;
};

//ContextModel3DBuffer
typedef struct context_model_buff_t context_model_buff_t;
struct context_model_buff_t
{
	context_model_t *ctx;
	int size_x;
	int size_y;
	int size_xy;
	const byte *ref_ctx_model;
};

//estos contextos en el HM estan dentro de TEncSbac
typedef struct entropy_model_t entropy_model_t;
struct entropy_model_t
{
	context_model_buff_t	cu_split_flag_model;
	context_model_buff_t	cu_skip_flag_model;
	context_model_buff_t	cu_merge_flag_model;
	context_model_buff_t	cu_merge_idx_model;
	context_model_buff_t	cu_part_size_model;
	context_model_buff_t	cu_pred_mode_flag_model;
	context_model_buff_t	cu_intra_pred_model;
	context_model_buff_t	cu_chroma_pred_model;
	context_model_buff_t	cu_delta_qp_model;
	context_model_buff_t	cu_inter_dir_model;
	context_model_buff_t	cu_ref_pic_model;
	context_model_buff_t	cu_mvd_model;
	context_model_buff_t	cu_qt_cbf_model;
	context_model_buff_t	cu_trans_subdiv_flag_model;
	context_model_buff_t	cu_qt_root_cbf_model;
	context_model_buff_t	cu_sig_coeff_group_model;
	context_model_buff_t	cu_sig_model;
	context_model_buff_t	cu_ctx_last_x_model;
	context_model_buff_t	cu_ctx_last_y_model;
	context_model_buff_t	cu_one_model;
	context_model_buff_t	cu_abs_model;
	context_model_buff_t	mvp_idx_model;
	context_model_buff_t	cu_amp_model;
	context_model_buff_t	sao_merge_model;
	context_model_buff_t	sao_type_model;
	context_model_buff_t	transform_skip_model;
	context_model_buff_t	transquant_bypass_flag_model;
};


typedef enum 
{
	EE_ENCODER,
	EE_COUNTER,
	EE_INVALID
}enc_env_enum;

//equivale a TEncSbac
typedef struct enc_env_t enc_env_t;
struct enc_env_t
{
	context_model_t	*contexts;
	binary_model_t	*b_ctx;
	entropy_model_t	*e_ctx;
	bitstream_t		*bs;
	enc_env_enum		type;
//	function pointers
	void (*ee_reset_bits)(binary_model_t* bm);
	void (*ee_start)(binary_model_t* bm);
	uint (*ee_bitcnt)(bitstream_t *bs, binary_model_t* bm);
	void (*ee_encode_bin)( enc_env_t* ee, context_model_t *cm, uint binValue);
	void (*ee_encode_bin_TRM)( enc_env_t* ee, uint binValue);
	void (*ee_encode_bins_EP)( enc_env_t* ee, uint binValues, int numBins );
	void (*ee_encode_bin_EP)( enc_env_t* ee, uint binValue);
	void (*ee_finish)( enc_env_t* ee);
};

typedef struct rate_distortion_t rate_distortion_t;
struct rate_distortion_t
{
  double                  lambda;
  double                  sqrt_lambda;
  uint                    lambda_SAD;
  uint                    lambda_SSE;
  double                  frame_lambda;
};


#define MAX_NUM_REF_PICS            16          ///< max. number of pictures used for reference
#define MAX_NUM_REF                 16          ///< max. number of entries in picture reference list

typedef struct ref_pic_set_t ref_pic_set_t;
struct ref_pic_set_t
{
	int inter_ref_pic_set_prediction_flag;
	//..................
	int num_negative_pics;
	int num_positive_pics;
	int delta_poc_s0[MAX_NUM_REF_PICS];
	int used_by_curr_pic_S0_flag[MAX_NUM_REF_PICS];
};

typedef struct slice_t slice_t;
struct slice_t
{
	unsigned int slice_index; //indice en la lista de slices
	unsigned int slice_type;
	unsigned int nalu_type;
	unsigned int first_cu_address;//address of first coding tree block in the slice
	unsigned int curr_cu_address;//address of current coding tree block in the slice
	unsigned int last_cu_address;//m_uiDependentSliceCurEndCUAddr - address of current coding tree block in the slice partition units (256 per CU)
	unsigned int is_dependent_slice;//
	unsigned int slice_temporal_layer_non_reference_flag;//
	unsigned int slice_temporal_mvp_enable_flag;
	unsigned int disable_deblocking_filter_flag;
	unsigned int slice_loop_filter_across_slices_enabled_flag;

	sps_t		*sps;
	pps_t		*pps;
	
	unsigned int qp;
	unsigned int poc;
	unsigned int depth;
	unsigned int sublayer;//TLayer
	unsigned int referenced;
	unsigned int num_ref_idx[2];
};

typedef struct picture_t picture_t;
struct picture_t
{
	slice_t	slice;
	video_frame_t	*img2encode;
};



typedef struct henc_thread_t henc_thread_t;
typedef struct hvenc_t hvenc_t;


typedef struct low_level_funcs_t low_level_funcs_t;
struct low_level_funcs_t
{
	uint32_t (*sad)(uint8_t * src, uint32_t src_stride, uint8_t * pred, uint32_t pred_stride, int size);
	uint32_t (*ssd)(uint8_t * src, uint32_t src_stride, uint8_t * pred, uint32_t pred_stride, int size);
	void (*predict)(uint8_t * __restrict orig, int orig_stride, uint8_t* __restrict pred, int pred_stride, int16_t * __restrict residual, int residual_stride, int size);
	void (*reconst)(uint8_t* pred, int pred_stride, int16_t * residual, int residual_stride, uint8_t* decoded, int decoded_stride, int size);
	uint32_t (*modified_variance)(uint8_t *__restrict p, int size, int stride, int modif);

	void (*create_intra_planar_prediction)(henc_thread_t* et, uint8_t *ref_wnd, int ref_wnd_stride_2D, int16_t  *adi_pred_buff, int adi_size, int cu_size, int cu_size_shift);
	void (*create_intra_angular_prediction)(henc_thread_t* et, ctu_info_t* ctu, uint8_t *ref_wnd, int ref_wnd_stride_2D, int16_t  *adi_pred_buff, int adi_size, int cu_size, int cu_mode, int is_luma);

	void (*quant)(henc_thread_t* et, ctu_info_t *ctu, int16_t* src, int16_t* dst, int scan_mode, int depth, int comp, int cu_mode, int is_intra, int *ac_sum, int cu_size);
	void (*inv_quant)(henc_thread_t* et, ctu_info_t *ctu, short * __restrict src, short * __restrict dst, int depth, int comp, int is_intra, int cu_size);

	void (*transform)(int bitDepth, int16_t *block,int16_t *coeff, int block_size, int iWidth, int iHeight, int width_shift, int height_shift, unsigned short uiMode, int16_t *aux);
	void (*itransform)(int bitDepth, int16_t *block,int16_t *coeff, int block_size, int iWidth, int iHeight, unsigned int uiMode, int16_t *aux);
};





#define NUM_QUANT_WNDS			(MAX_PARTITION_DEPTH+1)
#define NUM_DECODED_WNDS		(MAX_PARTITION_DEPTH+1)
#define NUM_CBF_BUFFS			MAX_PARTITION_DEPTH
#define NUM_CBF_BUFFS_CHROMA	5		//0 not used, U_COMP, V_COMP to consolidate best mode, U_COMP+2, V_COMP+2 for curr operation
#define NUM_REFF_WNDS		5

struct henc_thread_t
{
	hvenc_t		*ed;
	uint		index;
	int			wfpp_enable;
	int			wfpp_num_threads;
	hmr_sem_t	synchro_sem;
	hmr_sem_ptr	synchro_signal;
	hmr_sem_ptr	synchro_wait;

	bitstream_t	*bs;

	//header info
	vps_t	*vps;
	sps_t	*sps;
	pps_t	*pps;

	//Encoder Cfg	
	//Encoding layer
	int				pict_width[3], pict_height[3];
	int				pict_width_in_cu, pict_height_in_cu;
	int				pict_total_cu;
	int				ctu_width[3], ctu_height[3];
	int				ctu_group_size;

	//cfg
	int				max_cu_size;
	int				max_cu_size_shift;//log2 del tama�o del CU maximo
	int				max_cu_size_shift_chroma;//log2 del tama�o del CU maximo
	int				max_intra_tr_depth;
	int				max_inter_tr_depth;
	int				max_pred_partition_depth;//max depth for prediction

	int				num_partitions_in_cu;
	int				num_partitions_in_cu_shift;
	int				mincu_mintr_shift_diff;//log2		//g_uiAddCUDepth
	int				max_cu_depth;							//m_uiMaxCUDepth
	int				min_cu_size;
	int				min_cu_size_shift;//log2
	int				min_tu_size_shift;//log2
	int				max_tu_size_shift;//log2

	int				profile;
	int				bit_depth;
	int				rd_mode;
	int				performance_mode;


	//pointers to common buffers for all threads	 
//	picture_t			*current_pict;
//	wnd_t			*curr_decoding_wnd;

	//-------these are for abs_index partitions---------------------
//	unsigned short		*abs2raster_table; //g_auiZscanToRaster en HM
//	unsigned short		*raster2abs_table; //g_auiRasterToZscan en HM
	//--------------------------------------------------------------
//	uint				*scan_pyramid[NUM_SCAN_MODES][MAX_CU_DEPTHS];//[4][7]
//	int					*quant_pyramid[NUM_SCALING_MODES][NUM_SCALING_LISTS][NUM_SCALING_REM_LISTS];//[4][6][6]
//	int					*dequant_pyramid[NUM_SCALING_MODES][NUM_SCALING_LISTS][NUM_SCALING_REM_LISTS];//[4][6][6]
//	double				*scaling_error_pyramid[NUM_SCALING_MODES][NUM_SCALING_LISTS][NUM_SCALING_REM_LISTS];//[4][6][6]//quizas esto lo tendriamos que pasar a int. tiene valores muy bajos

	int					*partition_depth_start;//start of depths in the partition_info list
	cu_partition_info_t	*partition_info;//recursive structure list to store the state of the recursive computing stages

	//current processing state and buffers
	int				cu_current, cu_next;
	int				cu_current_x, cu_current_y;
//	int				cu_num_pixels, cu_num_pixels_shift;
//	int				frame_num, idr_num;

//	wnd_t			fw_wnd;
//	wnd_t			bw_wnd;

	wnd_t			curr_mbs_wnd;									//original MBs to be coded
//	wnd_t			*curr_prediction_wnd;
	wnd_t			prediction_wnd;									//prediction applied to original MBs
	wnd_t			residual_wnd;									//residual after substracting prediction
	wnd_t			transform_quant_wnd[NUM_QUANT_WNDS];			//for transform coefficients and quantification 
	wnd_t			itransform_iquant_wnd;							//for itransform coefficients and iquantification
	wnd_t			decoded_mbs_wnd[NUM_DECODED_WNDS];

	//intra predicition
	short				(*adi_pred_buff);//this buffer holds the left column and top row for intra pred (bottom2top and left2right)
	short				(*adi_filtered_pred_buff);//this buffer holds the left column and top row for intra pred (bottom2top and left2right)
	short				(*top_pred_buff);//intermediate buffer to calculate intra prediction samples
	short				(*left_pred_buff);//intermediate buffer to calculate intra prediction samples
	short				(*bottom_pred_buff);//intermediate buffer to calculate intra prediction samples
	short				(*right_pred_buff);//intermediate buffer to calculate intra prediction samples
	int					adi_size;//tama�o del patron de prediccion (el buffer es el doble)
	short				(*pred_aux_buff);
	int					pred_aux_buff_size;//tama�o del buffer auxiliar
	short				(*aux_buff);
	byte				(*cabac_aux_buff);
	int					cabac_aux_buff_size;
	ctu_info_t			*curr_ctu_group_info;	//esto esta pensado para poder ser una ventana peque�a, podria haber otra ventana con los macrobloques que sirven como referencia (left, top-left, top, top-right)
	byte				*cbf_buffs[NUM_PICT_COMPONENTS][NUM_CBF_BUFFS];
	byte				*cbf_buffs_chroma[NUM_PICT_COMPONENTS];//processing buffers for iteration buff
	byte				*intra_mode_buffs[NUM_PICT_COMPONENTS][NUM_CBF_BUFFS];
	byte				*tr_idx_buffs[NUM_CBF_BUFFS];

//	ctu_info_t			*ctu_info_t;//[MAX_MB_GROUP_SIZE];
	ctu_info_t			*ctu_rd;//[MAX_MB_GROUP_SIZE];

	enc_env_t			*ee;//encoding enviroment of the processing element
	enc_env_t			*ec;//encoding counter  of the processing element

	//rate distortion
	rate_distortion_t	rd;
	low_level_funcs_t	*funcs;
};


#define MAX_NUM_THREADS			32
#define NUM_INPUT_FRAMES		2
#define NUM_OUTPUT_NALUS		(2*NUM_INPUT_FRAMES)
#define NUM_OUTPUT_NALUS_MASK	(NUM_OUTPUT_NALUS-1)
struct hvenc_t
{
	int num_encoded_frames;

//	henc_thread_t	_thread;//*encoders_list;
	henc_thread_t	*thread[MAX_NUM_THREADS];//*encoders_list;
	hmr_thread_t	hthreads[MAX_NUM_THREADS];
	hmr_thread_t	encoder_thread;
	int				run;
	//nalus
	nalu_t		vps_nalu;
	nalu_t		sps_nalu;
	nalu_t		pps_nalu;
	nalu_t		slice_nalu_list[NUM_OUTPUT_NALUS];//slice
	nalu_t		*slice_nalu;//slice
	bitstream_t	slice_bs;//slice information previous to nalu_ebsp conversion
	bitstream_t	*aux_bs;//list of bitstreams for coef wfpp encoding
	int			num_sub_streams;
	uint		*sub_streams_entry_point_list;

	//header info
	vps_t	vps;
	sps_t	sps;
	pps_t	pps;

	//Encoder Cfg	
	//Encoding layer
	int				gop_size;
	int				num_b;
//	img_pool_t		img_list;
	int				pic_interlaced, mb_interlaced;
	unsigned int	conformance_mode;
	unsigned int	pad_left, pad_right;
	unsigned int	pad_top, pad_bottom;
	int				pict_width[3], pict_height[3];
	int				pict_width_in_cu, pict_height_in_cu;
	int				pict_total_cu;
	int				ctu_width[3], ctu_height[3];
	int				ctu_group_size;
	int				blocks_per_macroblock;

	//cfg
	int				frame_rate;
	int				max_cu_size;
	int				max_cu_size_shift;//log2 del tama�o del CU maximo
	int				max_cu_size_shift_chroma;//log2 del tama�o del CU maximo
	int				max_intra_tr_depth;
	int				max_inter_tr_depth;
	int				max_pred_partition_depth;//max depth for prediction
	int				wfpp_enable;
	uint			wfpp_num_threads;

	//
	int				num_partitions_in_cu;
	int				num_partitions_in_cu_shift;
	int				mincu_mintr_shift_diff;//log2		//g_uiAddCUDepth
	int				max_cu_depth;							//m_uiMaxCUDepth
	int				min_cu_size;
	int				min_cu_size_shift;//log2
	int				min_tu_size_shift;//log2
	int				max_tu_size_shift;//log2

	int				partition_depth_start[MAX_PARTITION_DEPTH];//start of depths in the partition_info list
//	cu_partition_info_t			*partition_info;//information for rd

	int				profile;
	int				bit_depth;
	int				max_sublayers;
	int				max_layers;
		
	profile_tier_level_t	ptl;

	//current picture_t Config
	picture_t		current_pict;
	wnd_t			*curr_ref_wnd;
	wnd_t			ref_wnds[NUM_REFF_WNDS];

	void			*cont_empty_reference_wnds;//for decoding and reference frames
	int				last_poc, last_idr, num_pictures;
	int				num_ref_lists;
	int				num_refs_idx_active_list[2];
	int				num_ref_frames;

	int				slice_type;
	int				pict_qp;

//	int				mb_current, mb_next;
//	int				mb_current_x, mb_current_y;
//	int				mb_num_pixels, mb_num_pixels_shift;
//	int				frame_num, idr_num;

	wnd_t			fw_wnd;
	wnd_t			bw_wnd;

//	wnd_t			curr_mbs_wnd;									//original MBs to be coded
//	wnd_t			*curr_prediction_wnd;
//	wnd_t			prediction_wnd;									//prediction applied to original MBs
//	wnd_t			residual_wnd;									//residual after substracting prediction
//	wnd_t			transform_quant_wnd[NUM_QUANT_WNDS];			//for transform coefficients and quantification 
//	wnd_t			itransform_iquant_wnd;							//for itransform coefficients and iquantification
//	wnd_t			decoded_mbs_wnd[NUM_DECODED_WNDS];

	//intra predicition
	ctu_info_t		*ctu_info_t;//[MAX_MB_GROUP_SIZE];
//	ctu_info_t		*ctu_rd;//[MAX_MB_GROUP_SIZE];
//	short				(*__restrict adi_pred_buff);//this buffer holds the left column and top row for intra pred (bottom2top and left2right)
//	short				(*__restrict adi_filtered_pred_buff);//this buffer holds the left column and top row for intra pred (bottom2top and left2right)
//	short				(*__restrict top_pred_buff);//intermediate buffer to calculate intra prediction samples
//	short				(*__restrict left_pred_buff);//intermediate buffer to calculate intra prediction samples
//	short				(*__restrict bottom_pred_buff);//intermediate buffer to calculate intra prediction samples
//	short				(*__restrict right_pred_buff);//intermediate buffer to calculate intra prediction samples
//	int					adi_size;//tama�o del patron de prediccion (el buffer es el doble)
//	short				(*__restrict pred_aux_buff);
//	int					pred_aux_buff_size;//tama�o del buffer auxiliar
//	short				(*__restrict aux_buff);
//	byte				(*__restrict cabac_aux_buff);
//	int					cabac_aux_buff_size;
//	ctu_info_t			*curr_ctu_group_info;	//esto esta pensado para poder ser una ventana peque�a, podria haber otra ventana con los macrobloques que sirven como referencia (left, top-left, top, top-right)
//	byte				*cbf_buffs[NUM_PICT_COMPONENTS][NUM_CBF_BUFFS];
//	byte				*cbf_buffs_chroma[NUM_PICT_COMPONENTS];//processing buffers for iteration buff
//	byte				*intra_mode_buffs[NUM_PICT_COMPONENTS][NUM_CBF_BUFFS];
//	byte				*tr_idx_buffs[NUM_CBF_BUFFS];
//	byte				*pred_depth_buff;
//	byte				*part_size_buff;
	int					performance_mode;
	int					rd_mode;
	//scan tables	 
	//-------these are for abs_index partitions---------------------
	unsigned short		*abs2raster_table; //g_auiZscanToRaster en HM
	unsigned short		*raster2abs_table; //g_auiRasterToZscan en HM
	//--------------------------------------------------------------
	ushort				*ang_table;//for angular intra prediction    
	ushort				*inv_ang_table;//for angular intra prediction
	uint				*scan_pyramid[NUM_SCAN_MODES][MAX_CU_DEPTHS];//[4][7]
	int					*quant_pyramid[NUM_SCALING_MODES][NUM_SCALING_LISTS][NUM_SCALING_REM_LISTS];//[4][6][6]
	int					*dequant_pyramid[NUM_SCALING_MODES][NUM_SCALING_LISTS][NUM_SCALING_REM_LISTS];//[4][6][6]
	double				*scaling_error_pyramid[NUM_SCALING_MODES][NUM_SCALING_LISTS][NUM_SCALING_REM_LISTS];//[4][6][6]//quizas esto lo tendriamos que pasar a int. tiene valores muy bajos

	//reference pictures
	int					ref_pic_set_index;
	int					num_short_term_ref_pic_sets;
	ref_pic_set_t		*ref_pic_set_list;//[MAX_REF_PIC_SETS];
	int					num_long_term_ref_pic_sets;

	//arithmetic coding
	uint				num_ee;
	enc_env_t				**ee_list;//encoding enviroment list hmr_container 
	uint				num_ec;
	enc_env_t				*ec_list;//encoding enviroment list
//	enc_env_t				*ee;//encoding enviroment of the processing element
//	enc_env_t				*ec;//encoding counter  of the processing element

	//rate distortion
	rate_distortion_t	rd;
	low_level_funcs_t	funcs;

	//input and output
	video_frame_t	input_frames[NUM_INPUT_FRAMES];
	void			*input_hmr_container;
	nalu_set_t		output_nalus[NUM_OUTPUT_NALUS];
	void			*output_hmr_container;

#ifdef COMPUTE_METRICS
	double			current_psnr[3];
	double			accumulated_psnr[3];
//	FILE			*f_psnr;
#endif
};

#endif  /* __HOMER_HEVC_PRIVATE_H__*/
