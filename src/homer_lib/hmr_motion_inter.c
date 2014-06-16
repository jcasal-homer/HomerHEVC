/*****************************************************************************
 * hmr_motion_inter.c : homerHEVC encoding library
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
#include <limits.h>
#include <memory.h>


#include "hmr_private.h"
#include "hmr_common.h"
#include "hmr_profiler.h"

#include "hmr_sse42_functions.h"


int encode_inter_cu(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, int depth, int cu_mode, PartSize part_size_type, int *curr_sum, motion_vector_t *mv, int gcnt)//depth = prediction depth
{		
	int ssd_;
	int pred_buff_stride, orig_buff_stride, residual_buff_stride, residual_dec_buff_stride, decoded_buff_stride;
	uint8_t *orig_buff, *decoded_buff;
	int16_t *pred_buff, *residual_buff, *residual_dec_buff, *quant_buff, *iquant_buff;
	uint8_t *cbf_buff = NULL;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL;
	int inv_depth;//, diff;//, is_filtered;
		
	int curr_depth = curr_partition_info->depth;
	int curr_part_x = curr_partition_info->x_position;
	int curr_part_y = curr_partition_info->y_position;
	int curr_part_size = curr_partition_info->size;
	int curr_part_size_shift = et->max_cu_size_shift-curr_depth;
	int curr_scan_mode = find_scan_mode(TRUE, TRUE, curr_part_size, cu_mode, 0);

	ctu->per = ctu->qp/6;
	ctu->rem = ctu->qp%6;

	quant_wnd = &et->transform_quant_wnd[curr_depth];
	decoded_wnd = &et->decoded_mbs_wnd[curr_depth];
	cbf_buff = et->cbf_buffs[Y_COMP][curr_depth];

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, Y_COMP);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
	orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, Y_COMP);
	residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_dec_buff_stride = WND_STRIDE_2D(et->residual_dec_wnd, Y_COMP);
	residual_dec_buff = WND_POSITION_2D(int16_t *, et->residual_dec_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, Y_COMP, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));
	iquant_buff = WND_POSITION_1D(int16_t  *, et->itransform_iquant_wnd, Y_COMP, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));
	decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, Y_COMP);
	decoded_buff = WND_POSITION_2D(uint8_t *, *decoded_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, Y_COMP, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));

	inv_depth = (et->max_cu_size_shift - curr_depth);
//	diff = min(abs(cu_mode - HOR_IDX), abs(cu_mode - VER_IDX));

	//2d -> 1D buffer
	et->funcs->transform(et->bit_depth, residual_buff, et->pred_aux_buff, residual_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, curr_part_size_shift, cu_mode, quant_buff);//usamos quant buff como auxiliar

	et->funcs->quant(et, ctu, et->pred_aux_buff, quant_buff, curr_scan_mode, curr_depth, Y_COMP, cu_mode, 0, curr_sum, curr_part_size);//Si queremos quitar el bit de signo necesitamos hacerlo en dos arrays distintos
	curr_partition_info->sum = *curr_sum;

	curr_partition_info->inter_cbf[Y_COMP] = (( curr_partition_info->sum ? 1 : 0 ) << (curr_depth-depth+(part_size_type==SIZE_NxN)));
	curr_partition_info->inter_tr_idx[Y_COMP]  = (curr_depth-depth + (part_size_type==SIZE_NxN));
//	memset(&et->cbf_buffs[Y_COMP][curr_depth][curr_partition_info->abs_index], (( curr_partition_info->sum ? 1 : 0 ) << (curr_depth-depth+(part_size_type==SIZE_NxN))), curr_partition_info->num_part_in_cu*sizeof(cbf_buff[0]));//(width*width)>>4 num parts of 4x4 in partition
//	memset(&et->tr_idx_buffs[curr_depth][curr_partition_info->abs_index], (curr_depth-depth+(part_size_type==SIZE_NxN)), curr_partition_info->num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));//(width*width)>>4 num parts of 4x4 in partition
//	memset(&et->intra_mode_buffs[Y_COMP][curr_depth][curr_partition_info->abs_index], cu_mode, curr_partition_info->num_part_in_cu*sizeof(et->intra_mode_buffs[Y_COMP][0][0]));//(width*width)>>4 num parts of 4x4 in partition

	et->funcs->inv_quant(et, ctu, quant_buff, iquant_buff, curr_depth, Y_COMP, 0, curr_part_size);

	//1D ->2D buffer
	et->funcs->itransform(et->bit_depth, residual_dec_buff, iquant_buff, residual_buff_stride, curr_part_size, curr_part_size, cu_mode, et->pred_aux_buff);

	ssd_ = ssd_16(residual_buff, residual_buff_stride, residual_dec_buff, residual_buff_stride, curr_part_size);

	et->funcs->reconst(pred_buff, pred_buff_stride, residual_dec_buff, residual_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
//	ssd_ = et->funcs->ssd(orig_buff, orig_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	return ssd_;
}


int encode_inter_cu_chroma(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, int component, int depth, int cu_mode, PartSize part_size_type, int *curr_sum, motion_vector_t *mv, int gcnt)//depth = prediction depth
{		
	int ssd_;
	int pred_buff_stride, orig_buff_stride, residual_buff_stride, residual_dec_buff_stride, decoded_buff_stride;
	uint8_t *orig_buff, *decoded_buff;
	int16_t *pred_buff, *residual_buff, *residual_dec_buff, *quant_buff, *iquant_buff;
	uint8_t *cbf_buff = NULL;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL;
	int inv_depth, diff;//, is_filtered;
		
	int curr_depth = curr_partition_info->depth;
	int curr_part_x = curr_partition_info->x_position;
	int curr_part_y = curr_partition_info->y_position;
	int curr_part_size = curr_partition_info->size_chroma;
	int curr_part_size_shift = et->max_cu_size_shift-curr_depth-1;//420
	int curr_scan_mode = find_scan_mode(TRUE, TRUE, curr_part_size, cu_mode, 0);
	double weight = pow( 2.0, (ctu->qp-ctu->qp_chroma)/3.0 ); 

	ctu->per = ctu->qp_chroma/6;
	ctu->rem = ctu->qp_chroma%6;

	quant_wnd = &et->transform_quant_wnd[curr_depth];
	decoded_wnd = &et->decoded_mbs_wnd[curr_depth];
	cbf_buff = et->cbf_buffs[component][curr_depth];

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, component);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, component);
	orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, component);
	residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_dec_buff_stride = WND_STRIDE_2D(et->residual_dec_wnd, component);
	residual_dec_buff = WND_POSITION_2D(int16_t *, et->residual_dec_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, component, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));
	iquant_buff = WND_POSITION_1D(int16_t  *, et->itransform_iquant_wnd, component, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));
	decoded_buff_stride = WND_STRIDE_2D(*decoded_wnd, component);
	decoded_buff = WND_POSITION_2D(uint8_t *, *decoded_wnd, component, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	quant_buff = WND_POSITION_1D(int16_t  *, *quant_wnd, component, gcnt, et->ctu_width, (curr_partition_info->abs_index<<et->num_partitions_in_cu_shift));

	inv_depth = (et->max_cu_size_shift - curr_depth);
	diff = min(abs(cu_mode - HOR_IDX), abs(cu_mode - VER_IDX));

	//2d -> 1D buffer
	PROFILER_RESET(inter_luma_tr)
	et->funcs->transform(et->bit_depth, residual_buff, et->pred_aux_buff, residual_buff_stride, curr_part_size, curr_part_size, curr_part_size_shift, curr_part_size_shift, cu_mode, quant_buff);//usamos quant buff como auxiliar
	PROFILER_ACCUMULATE(inter_luma_tr)

	PROFILER_RESET(inter_luma_q)
	et->funcs->quant(et, ctu, et->pred_aux_buff, quant_buff, curr_scan_mode, curr_depth, component, cu_mode, 0, curr_sum, curr_part_size);//Si queremos quitar el bit de signo necesitamos hacerlo en dos arrays distintos
	curr_partition_info->sum = *curr_sum;
	PROFILER_ACCUMULATE(inter_luma_q)


	curr_partition_info->inter_cbf[component] = (( curr_partition_info->sum ? 1 : 0 ) << (curr_depth-depth+(part_size_type==SIZE_NxN)));
	curr_partition_info->inter_tr_idx[component]  = (curr_depth-depth+(part_size_type==SIZE_NxN));
//	memset(&et->cbf_buffs[component][curr_depth][curr_partition_info->abs_index], (( curr_partition_info->sum ? 1 : 0 ) << (curr_depth-depth+(part_size_type==SIZE_NxN))), curr_partition_info->num_part_in_cu*sizeof(cbf_buff[0]));//(width*width)>>4 num parts of 4x4 in partition
//	memset(&et->tr_idx_buffs[curr_depth][curr_partition_info->abs_index], (curr_depth-depth+(part_size_type==SIZE_NxN)), curr_partition_info->num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));//(width*width)>>4 num parts of 4x4 in partition
//	memset(&et->intra_mode_buffs[component][curr_depth][curr_partition_info->abs_index], cu_mode, curr_partition_info->num_part_in_cu*sizeof(et->intra_mode_buffs[component][0][0]));//(width*width)>>4 num parts of 4x4 in partition

	et->funcs->inv_quant(et, ctu, quant_buff, iquant_buff, curr_depth, component, 0, curr_part_size);

	//1D ->2D buffer
	et->funcs->itransform(et->bit_depth, residual_dec_buff, iquant_buff, residual_buff_stride, curr_part_size, curr_part_size, cu_mode, et->pred_aux_buff);
	ssd_ = weight*ssd_16(residual_buff, residual_buff_stride, residual_dec_buff, residual_buff_stride, curr_part_size);

	et->funcs->reconst(pred_buff, pred_buff_stride, residual_dec_buff, residual_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
//	ssd_ = et->funcs->ssd(orig_buff, orig_buff_stride, decoded_buff, decoded_buff_stride, curr_part_size);
	return ssd_;
}


void hmr_motion_inter_uni(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint8_t *orig_buff, int orig_buff_stride, uint8_t *reference_buff, int reference_buff_stride, uint8_t *pred_buff, int pred_buff_stride,  int curr_part_global_x, int curr_part_global_y, int curr_part_size, int curr_part_size_shift, motion_vector_t *mv)
{
	int j;
	mv->hor_vector = 0;
	mv->ver_vector = 0;

	for(j=0;j<curr_part_size;j++)
	{
		memcpy(pred_buff, reference_buff, curr_part_size);
		pred_buff += pred_buff_stride;
		reference_buff += reference_buff_stride;
	}	
}



static const int diamond_small[][2] = {{0,-1},{-1,0},{1,0},{0,1}};
static const int diamond_big[][2] = {{0,-2},{-1,-1},{1,-1},{-2,0},{2,0},{-1,1},{1,1},{0,2}};

//static int MPGVENC_SPC_diamondsearch(const unsigned char *__restrict refImg,const unsigned char *__restrict srcImg, int x_pitch, int y_pitch, int i0, int j0, int sx, int sy,int *iminp,int *jminp, int *v0, MvStream *general_MV, int Niteraciones, int TH0, int pict_type, int M)

uint32_t hmr_motion_estimation(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint8_t *orig_buff, int orig_buff_stride, int16_t *reference_buff, int reference_buff_stride, int curr_part_global_x, 
						int curr_part_global_y, int init_x, int init_y, int curr_part_size, int curr_part_size_shift, int search_range_x, int search_range_y, int frame_size_x, int frame_size_y, motion_vector_t *mv)
{
	int i,j; 
	int xlow, xhigh, ylow, yhigh;
	uint32_t curr_best_sad, best_sad;
	int curr_best_x, curr_best_y, best_x, best_y;
	int dist = 1;

	xlow = curr_part_global_x - search_range_x;
	xhigh = curr_part_global_x + search_range_x;
	ylow = curr_part_global_y - search_range_y;//limite vertical superior
	yhigh = curr_part_global_y + search_range_y;//limite vertical inferior

	if (xlow<0)
		xlow = 0;
	if (xhigh> frame_size_x)
		xhigh = frame_size_x;
	if (ylow<0)
		ylow = 0;
	if (yhigh> frame_size_y)
		yhigh = frame_size_y;

	curr_best_x = init_x;
	curr_best_y = init_y;

	if (curr_best_x<xlow)
		curr_best_x=xlow;
	if (curr_best_x>xhigh-curr_part_size)
		curr_best_x=xhigh-curr_part_size;
	if (curr_best_y<ylow)
		curr_best_y=ylow;
	if (curr_best_y>yhigh-curr_part_size)
		curr_best_y=yhigh-curr_part_size;

	curr_best_sad = sad(orig_buff, orig_buff_stride, reference_buff+curr_best_y*reference_buff_stride+curr_best_x, reference_buff_stride, curr_partition_info->size);

	best_sad = curr_best_sad;
	best_x = curr_best_x;
	best_y = curr_best_y;
//	curr_best_x = 0;
//	curr_best_y = 0;

	for(i=0;i<sizeof(diamond_small)/sizeof(diamond_small[0]);i++)
	{
		int curr_x = best_x+diamond_small[i][0];
		int curr_y = best_y+diamond_small[i][1];
		uint32_t curr_sad;

		if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
		{
			curr_sad = sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_partition_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
			}
		}
	}

/*	best_sad = curr_best_sad;
	best_x = curr_best_x;
	best_y = curr_best_y;
	curr_best_x = 0;
	curr_best_y = 0;
*/
	best_x = curr_best_x-1;
	best_y = curr_best_y-1;
	for(dist = 1; 2*dist<search_range_x; )//dist*=2)
	{
		if(best_x != curr_best_x || best_y != curr_best_y)
		{
			best_sad = curr_best_sad;
			best_x = curr_best_x;
			best_y = curr_best_y;
		}
		else
			dist*=2;

		for(i=0;i<sizeof(diamond_big)/sizeof(diamond_big[0]);i++)
		{
			int curr_x = best_x+diamond_big[i][0]*dist;
			int curr_y = best_y+diamond_big[i][1]*dist;
			uint32_t curr_sad;

			if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
			{
				curr_sad = sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_partition_info->size);
				if(curr_sad < curr_best_sad)
				{
					curr_best_sad = curr_sad;
					curr_best_x = curr_x;
					curr_best_y = curr_y;
				}
			}
		}
	}

	best_sad = curr_best_sad;
	best_x = (curr_best_x & 0xfffffffe)+1;
	best_y = (curr_best_y & 0xfffffffe)+1;
//	curr_best_x = 0;
//	curr_best_y = 0;

/*	for(i=0;i<sizeof(diamond_small)/sizeof(diamond_small[0]);i++)
	{
		int curr_x = best_x+diamond_small[i][0];
		int curr_y = best_y+diamond_small[i][1];
		uint32_t curr_sad;

		if (curr_x>=xlow && curr_x<=xhigh && curr_y>=ylow && curr_y<=yhigh)
		{
			curr_sad = sad(orig_buff, orig_buff_stride, reference_buff+curr_y*reference_buff_stride+curr_x, reference_buff_stride, curr_partition_info->size);
			if(curr_sad < curr_best_sad)
			{
				curr_best_sad = curr_sad;
				curr_best_x = curr_x;
				curr_best_y = curr_y;
			}
		}
	}

	best_sad = curr_best_sad;
	best_x = curr_best_x;
	best_y = curr_best_y;
*/
	mv->hor_vector = 4*best_x;
	mv->ver_vector = 4*best_y;

	return best_sad;
}


//uint32_t hmr_motion_estimation(henc_thread_t* et, ctu_info_t* ctu, cu_partition_info_t* curr_partition_info, uint8_t *orig_buff, int orig_buff_stride, uint8_t *reference_buff, int reference_buff_stride, int curr_part_global_x, 
//							int curr_part_global_y, int init_x, int init_y, int curr_part_size, int curr_part_size_shift, int search_range_x, int search_range_y, int frame_size_x, int frame_size_y, motion_vector_t *mv)
void hmr_motion_compensation_luma(henc_thread_t *et, ctu_info_t *ctu, cu_partition_info_t* curr_partition_info, int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int curr_part_size, int curr_part_size_shift, motion_vector_t *mv)
{
	int x_fraction = mv->hor_vector&0x3;
	int y_fraction = mv->ver_vector&0x3;

	if(x_fraction==0 && y_fraction==0)
	{
		int j, i;
		int x_vect = mv->hor_vector>>2;
		int y_vect = mv->ver_vector>>2;
		reference_buff+= y_vect*reference_buff_stride+x_vect;

		for(j=0;j<curr_part_size;j++)
		{
			for(i=0;i<curr_part_size;i++)
			{	
				pred_buff[i] = reference_buff[i];
			}
			reference_buff+=reference_buff_stride;
			pred_buff+=pred_buff_stride;
		}
		//copy samples
	}
}


#define NTAPS_LUMA        8 ///< Number of taps for luma
#define NTAPS_CHROMA      4 ///< Number of taps for chroma
#define IF_INTERNAL_PREC 14 ///< Number of bits for internal precision
#define IF_FILTER_PREC    6 ///< Log2 of sum of filter taps
#define IF_INTERNAL_OFFS (1<<(IF_INTERNAL_PREC-1)) ///< Offset used internally

int16_t chroma_filter_coeffs[8][NTAPS_CHROMA] =
{
  {  0, 64,  0,  0 },
  { -2, 58, 10, -2 },
  { -4, 54, 16, -2 },
  { -6, 46, 28, -4 },
  { -4, 36, 36, -4 },
  { -4, 28, 46, -6 },
  { -2, 16, 54, -4 },
  { -2, 10, 58, -2 }
};

void hmr_motion_filter_chroma(int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int16_t* coeffs, int width, int height, int is_vertical, int is_first, int is_last)
{
	int bit_depth = 8;
	int num_taps = NTAPS_CHROMA;//argument
	int ref_stride = ( is_vertical ) ? reference_buff_stride : 1;

	int offset;
	short maxVal;
	int headRoom = IF_INTERNAL_PREC - bit_depth;
	int shift = IF_FILTER_PREC;
	int row, col;
	short c[4];

	reference_buff -= ( num_taps/2 - 1 ) * ref_stride;

	c[0] = coeffs[0];
	c[1] = coeffs[1];
	c[2] = coeffs[2];
	c[3] = coeffs[3];

	if ( is_last )
	{
		shift += (is_first) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (is_first) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
		maxVal = (1 << bit_depth) - 1;
	}
	else
	{
		shift -= (is_first) ? headRoom : 0;
		offset = (is_first) ? -IF_INTERNAL_OFFS << shift : 0;
		maxVal = 0;
	}

	for (row = 0; row < height; row++)
	{
		for (col = 0; col < width; col++)
		{
			int sum;
			short val;
			sum  = reference_buff[ col + 0 * ref_stride] * c[0];
			sum += reference_buff[ col + 1 * ref_stride] * c[1];
			sum += reference_buff[ col + 2 * ref_stride] * c[2];
			sum += reference_buff[ col + 3 * ref_stride] * c[3];

			val = ( sum + offset ) >> shift;
			if ( is_last)
			{
				val = clip(val,0,maxVal);
			}
			pred_buff[col] = val;
		}

		reference_buff += reference_buff_stride;
		pred_buff += pred_buff_stride;
	}

}






void hmr_motion_compensation_chroma(henc_thread_t* et, int16_t *reference_buff, int reference_buff_stride, int16_t *pred_buff, int pred_buff_stride, int curr_part_size, int curr_part_size_shift, motion_vector_t *mv)
{
	int is_bi_predict = 0;
	int x_fraction = mv->hor_vector&0x7;
	int y_fraction = mv->ver_vector&0x7;

	int x_vect = mv->hor_vector>>3;
	int y_vect = mv->ver_vector>>3;
	reference_buff+= y_vect*reference_buff_stride+x_vect;


	if(x_fraction==0 && y_fraction==0)
	{
		int j, i;
		//copy samples
		for(j=0 ; j<curr_part_size ; j++)
		{

			for(i=0 ; i<curr_part_size ; i++)
			{
				pred_buff[i] = reference_buff[i];
			}
//			memcpy(pred_buff, reference_buff, curr_part_size);
			reference_buff+=reference_buff_stride;
			pred_buff+=pred_buff_stride;
		}
	}
	else if(x_fraction == 0)
	{
		//vertical filter 
		hmr_motion_filter_chroma(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, chroma_filter_coeffs[y_fraction], curr_part_size, curr_part_size, 1, TRUE, !is_bi_predict);
	}
	else if(y_fraction == 0)
	{
		//horizontal filter 
		hmr_motion_filter_chroma(reference_buff, reference_buff_stride, pred_buff, pred_buff_stride, chroma_filter_coeffs[x_fraction], curr_part_size, curr_part_size, 0, TRUE, !is_bi_predict);	
	}
	else //if(x_fraction!=0 && y_fraction!=0)
	{
		int filter_size = NTAPS_CHROMA;
		int half_filter_size = filter_size>>1;
		int16_t *temp_buff = WND_DATA_PTR(int16_t *, et->filtered_blocks_temp_wnd[0], Y_COMP);
		int temp_buff_stride = WND_STRIDE_2D(et->filtered_blocks_temp_wnd[0], Y_COMP);
		//horizontal
		hmr_motion_filter_chroma(reference_buff - (half_filter_size-1)*reference_buff_stride, reference_buff_stride, temp_buff, temp_buff_stride, chroma_filter_coeffs[x_fraction], curr_part_size, curr_part_size+half_filter_size+1, 0, TRUE, FALSE);
		//vertical filter 
		hmr_motion_filter_chroma(temp_buff + (half_filter_size-1)*temp_buff_stride, temp_buff_stride, pred_buff, pred_buff_stride, chroma_filter_coeffs[y_fraction], curr_part_size, curr_part_size, 1, FALSE, !is_bi_predict);
	}

}


void set_mv_and_ref_idx(henc_thread_t* et, cu_partition_info_t *curr_partition_info, motion_vector_t *mv, int ref_idx)
{
	int i;
	motion_vector_t *aux_mv = (motion_vector_t *)&et->mv_ref0[Y_COMP][curr_partition_info->depth][curr_partition_info->abs_index];
	for(i=0;i<curr_partition_info->num_part_in_cu;i++)
	{
		aux_mv->hor_vector = mv->hor_vector;
		aux_mv->ver_vector = mv->ver_vector;
		aux_mv++;
	}
	memset(&et->ref_idx0[Y_COMP][curr_partition_info->depth][curr_partition_info->abs_index], ref_idx, curr_partition_info->num_part_in_cu);
}

int encode_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt, int depth, int part_position, PartSize part_size_type)
{
	int k;
	int cu_mode;
	double distortion = 0.;

	uint32_t sad, cost, best_cost;
	slice_t *currslice = &et->ed->current_pict.slice;
	ctu_info_t *ctu_rd = et->ctu_rd;
	int pred_buff_stride, orig_buff_stride, reference_buff_stride, residual_buff_stride;
	int pred_buff_stride_chroma, orig_buff_stride_chroma, reference_buff_stride_chroma, residual_buff_stride_chroma;
	uint8_t  *orig_buff, *orig_buff_u, *orig_buff_v;
	int16_t  *reference_buff_cu_position, *reference_buff_cu_position_u, *reference_buff_cu_position_v;
	int16_t *pred_buff, *pred_buff_u, *pred_buff_v, *residual_buff, *residual_buff_u, *residual_buff_v;
	wnd_t *reference_wnd=NULL;//, *resi_wnd = NULL;
	uint8_t *cbf_buff = NULL;
	int curr_part_size, curr_part_size_shift;
	int curr_part_size_chroma, curr_part_size_shift_chroma;
	int curr_part_x, curr_part_y, curr_part_global_x, curr_part_global_y;
	int curr_part_x_chroma, curr_part_y_chroma, curr_part_global_x_chroma, curr_part_global_y_chroma;
	int curr_depth = depth;
	cu_partition_info_t*	parent_part_info;
	cu_partition_info_t*	curr_partition_info;
	int curr_sum = 0, best_sum;
	int num_part_in_cu;
	int partition_cost;
	int cu_min_tu_size_shift;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int max_tr_depth, max_tr_processing_depth;
	int initial_state, end_state;
//	int cbf_split[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int acc_cost[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int bitcost_cu_mode;
	int log2cu_size;
	int cu_x_position;
	motion_vector_t mv;
	int ref_idx = 0;

	//rate-control - en motion_intra_chroma se modifica
	ctu->qp = currslice->qp;
	ctu->per = ctu->qp/6;
	ctu->rem = ctu->qp%6;

	curr_partition_info = &ctu->partition_list[et->partition_depth_start[curr_depth]]+part_position;

	parent_part_info = curr_partition_info->parent;

	curr_depth = curr_partition_info->depth;
	curr_part_x = curr_partition_info->x_position;
	curr_part_y = curr_partition_info->y_position;
	curr_part_global_x = ctu->x[Y_COMP]+curr_part_x;
	curr_part_global_y = ctu->y[Y_COMP]+curr_part_y;
	curr_part_size = curr_partition_info->size;
	curr_part_size_shift = et->max_cu_size_shift-curr_depth;
	curr_part_x_chroma = curr_partition_info->x_position_chroma;
	curr_part_y_chroma = curr_partition_info->y_position_chroma;
	curr_part_global_x_chroma = ctu->x[CHR_COMP]+curr_part_x_chroma;
	curr_part_global_y_chroma = ctu->y[CHR_COMP]+curr_part_y_chroma;
	curr_part_size_chroma = curr_partition_info->size_chroma;
	curr_part_size_shift_chroma = et->max_cu_size_shift-curr_depth-1;//420
//	curr_adi_size = 2*2*curr_part_size+1;

	pred_buff_stride = WND_STRIDE_2D(et->prediction_wnd, Y_COMP);
	pred_buff = WND_POSITION_2D(int16_t *, et->prediction_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	pred_buff_stride_chroma = WND_STRIDE_2D(et->prediction_wnd, CHR_COMP);
	pred_buff_u = WND_POSITION_2D(int16_t *, et->prediction_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	pred_buff_v = WND_POSITION_2D(int16_t *, et->prediction_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

	orig_buff_stride = WND_STRIDE_2D(et->curr_mbs_wnd, Y_COMP);
	orig_buff = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	orig_buff_stride_chroma = WND_STRIDE_2D(et->curr_mbs_wnd, CHR_COMP);
	orig_buff_u = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	orig_buff_v = WND_POSITION_2D(uint8_t *, et->curr_mbs_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);

	residual_buff_stride = WND_STRIDE_2D(et->residual_wnd, Y_COMP);
	residual_buff = WND_POSITION_2D(int16_t *, et->residual_wnd, Y_COMP, curr_part_x, curr_part_y, gcnt, et->ctu_width);
	residual_buff_stride_chroma = WND_STRIDE_2D(et->residual_wnd, CHR_COMP);
	residual_buff_u = WND_POSITION_2D(int16_t *, et->residual_wnd, U_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);
	residual_buff_v = WND_POSITION_2D(int16_t *, et->residual_wnd, V_COMP, curr_part_x_chroma, curr_part_y_chroma, gcnt, et->ctu_width);


	reference_wnd = &currslice->ref_pic_list[REF_PIC_LIST_0][ref_idx]->img;//[0] up to now we only use one reference	
	reference_buff_stride = WND_STRIDE_2D(*reference_wnd, Y_COMP);
	reference_buff_cu_position = WND_POSITION_2D(int16_t *, *reference_wnd, Y_COMP, curr_part_global_x, curr_part_global_y, gcnt, et->ctu_width);
	reference_buff_stride_chroma = WND_STRIDE_2D(*reference_wnd, CHR_COMP);
	reference_buff_cu_position_u = WND_POSITION_2D(int16_t *, *reference_wnd, U_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);
	reference_buff_cu_position_v = WND_POSITION_2D(int16_t *, *reference_wnd, V_COMP, curr_part_global_x_chroma, curr_part_global_y_chroma, gcnt, et->ctu_width);


	best_cost = INT_MAX;
	curr_depth = curr_partition_info->depth;
	memset(depth_state, 0, sizeof(depth_state));

//	hmr_motion_inter_uni(et, ctu, curr_partition_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, pred_buff, pred_buff_stride, curr_part_global_x, curr_part_global_y, curr_part_size, curr_part_size_shift, &mv);
	sad = hmr_motion_estimation(et, ctu, curr_partition_info, orig_buff, orig_buff_stride, reference_buff_cu_position, reference_buff_stride, curr_part_global_x, curr_part_global_y, 0, 0, curr_part_size, curr_part_size_shift, 64, 64, et->pict_width[Y_COMP], et->pict_height[Y_COMP], &mv);	

	//set mvs and ref_idx
	curr_partition_info->inter_mv[REF_PIC_LIST_0] = mv;
	curr_partition_info->inter_ref_index[REF_PIC_LIST_0] = ref_idx;
//	set_mv_and_ref_idx(et, curr_partition_info, &mv, ref_idx);

	hmr_motion_compensation_luma(et, ctu, curr_partition_info, reference_buff_cu_position, reference_buff_stride, pred_buff, pred_buff_stride, curr_part_size, curr_part_size_shift, &mv);
	et->funcs->predict(orig_buff, orig_buff_stride, pred_buff, pred_buff_stride, residual_buff, residual_buff_stride, curr_part_size);

	hmr_motion_compensation_chroma(et, reference_buff_cu_position_u, reference_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv);
	hmr_motion_compensation_chroma(et, reference_buff_cu_position_v, reference_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, curr_part_size_chroma, curr_part_size_shift_chroma, &mv);	
	et->funcs->predict(orig_buff_u, orig_buff_stride_chroma, pred_buff_u, pred_buff_stride_chroma, residual_buff_u, residual_buff_stride_chroma, curr_part_size_chroma);
	et->funcs->predict(orig_buff_v, orig_buff_stride_chroma, pred_buff_v, pred_buff_stride_chroma, residual_buff_v, residual_buff_stride_chroma, curr_part_size_chroma);


	if(depth==0 && et->max_cu_size == MAX_CU_SIZE)
	{
		parent_part_info = &ctu->partition_list[et->partition_depth_start[depth]];
		curr_partition_info = parent_part_info->children[0];
		parent_part_info->cost = INT_MAX;//the parent is consolidated if 2Nx2N
		initial_state = part_position & 0x3;//initial_state = 0;
		end_state = initial_state;//end_state = 1;
	}
	else
	{
		curr_partition_info = &ctu->partition_list[et->partition_depth_start[curr_depth]]+part_position;
		parent_part_info = curr_partition_info->parent;
		initial_state = part_position & 0x3;
		end_state = initial_state+1;
	}

	cu_mode = REG_DCT;//if !IsIntra(abs_index)cu_mode = REG_DCT;

	curr_depth = curr_partition_info->depth;

	max_tr_depth = et->max_inter_tr_depth;

	log2cu_size = et->max_cu_size_shift-(depth-(part_size_type==SIZE_NxN));

	if(log2cu_size < et->min_tu_size_shift+et->max_intra_tr_depth-1+(part_size_type==SIZE_NxN))
		cu_min_tu_size_shift = et->min_tu_size_shift;
	else
	{
		cu_min_tu_size_shift = log2cu_size - (max_tr_depth - 1 + (part_size_type==SIZE_NxN));
		if(cu_min_tu_size_shift>MAX_TU_SIZE_SHIFT)
			cu_min_tu_size_shift=MAX_TU_SIZE_SHIFT;
	}
	max_tr_processing_depth = et->max_cu_size_shift-cu_min_tu_size_shift;
/*	max_tr_processing_depth = 0;

	if(curr_part_size_shift > et->max_tu_size_shift)
	{
		while(curr_part_size_shift > (et->max_tu_size_shift+max_tr_processing_depth)) max_tr_processing_depth++;
	}
	max_tr_processing_depth++;

	while((curr_part_size_shift-max_tr_processing_depth) < (et->max_cu_size_shift-et->max_cu_depth)) max_tr_processing_depth--;
*/


	memset(acc_cost, 0, sizeof(acc_cost));
	memset(depth_state, 0, sizeof(depth_state));

	depth_state[curr_depth] = initial_state;

	if(part_size_type==SIZE_NxN)
	{
		int iiii=0;
	}

	PROFILER_RESET(intra_luma_bucle3)
	while(curr_depth!=depth || depth_state[curr_depth]!=end_state)
	{		
		uint bit_cost = 0;
		uint dist_y, dist_u, dist_v;
		curr_partition_info = (parent_part_info==NULL)?curr_partition_info:parent_part_info->children[depth_state[curr_depth]];//if cu_size=64 we process 4 32x32 partitions, else just the curr_partition
		curr_depth = curr_partition_info->depth;

		dist_y = encode_inter_cu(et, ctu, curr_partition_info, depth, cu_mode, part_size_type, &curr_sum, &mv, gcnt);//depth = prediction depth

		dist_u = encode_inter_cu_chroma(et, ctu, curr_partition_info, U_COMP, depth, cu_mode, part_size_type, &curr_sum, &mv, gcnt);//depth = prediction depth
		dist_v = encode_inter_cu_chroma(et, ctu, curr_partition_info, V_COMP, depth, cu_mode, part_size_type, &curr_sum, &mv, gcnt);//depth = prediction depth

		curr_partition_info->distortion = dist_y+dist_u+dist_v;
		curr_partition_info->cost = curr_partition_info->distortion;

		acc_cost[curr_depth] += curr_partition_info->distortion;

		depth_state[curr_depth]++;

		//HM order
		if(curr_depth < max_tr_processing_depth)//go one level down
		{
			curr_depth++;
			parent_part_info = curr_partition_info;
		}
		else if(depth_state[curr_depth]==4)//consolidate 
		{	
			while(depth_state[curr_depth]==4 && (curr_depth > (depth)))
			{
				distortion = parent_part_info->children[0]->distortion+parent_part_info->children[1]->distortion+parent_part_info->children[2]->distortion+parent_part_info->children[3]->distortion;
				cost = distortion;

				depth_state[curr_depth] = 0;

				if(cost < parent_part_info->cost)
				{
//					acc_cost[parent_part_info->depth]+= cost - parent_part_info->cost;
					parent_part_info->cost = cost;
					parent_part_info->distortion = distortion;

					//consolidate in parent
					memcpy(&et->cbf_buffs[Y_COMP][parent_part_info->depth][parent_part_info->abs_index], &et->cbf_buffs[Y_COMP][curr_depth][parent_part_info->abs_index], parent_part_info->num_part_in_cu*sizeof(cbf_buff[0]));
//					memcpy(&et->intra_mode_buffs[Y_COMP][parent_part_info->depth][parent_part_info->abs_index], &et->intra_mode_buffs[Y_COMP][curr_depth][parent_part_info->abs_index], parent_part_info->num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));
					memcpy(&et->tr_idx_buffs[parent_part_info->depth][parent_part_info->abs_index], &et->tr_idx_buffs[curr_depth][parent_part_info->abs_index], parent_part_info->num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));

					//synchronize buffers for next iterations
//					synchronize_motion_buffers_luma(et, parent_part_info, &et->transform_quant_wnd[curr_depth], &et->transform_quant_wnd[curr_depth-1], &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[curr_depth-1], gcnt);
				}
				else
				{
//					synchronize_reference_buffs(et, parent_part_info, &et->decoded_mbs_wnd[curr_depth-1], &et->decoded_mbs_wnd[curr_depth], gcnt);	
				}

				acc_cost[curr_depth] = 0;
				curr_depth--;
				parent_part_info = parent_part_info->parent;
			}

			if(curr_depth+2 <= max_tr_processing_depth)//sync buffers of higher depth 
			{
				int aux_depth, aux_abs_index, aux_num_part_in_cu;
				cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
				aux_abs_index = aux_partition_info->abs_index;
				aux_num_part_in_cu = aux_partition_info->num_part_in_cu;
				for(aux_depth=curr_depth+2;aux_depth<=max_tr_processing_depth;aux_depth++)//+1 pq el nivel superior esta ya sincronizado y tenemos que sincronizar los siguientes
				{
//					synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);	
					if(et->rd_mode==1)
					{
//						memcpy(&et->intra_mode_buffs[Y_COMP][aux_depth][aux_abs_index], &et->intra_mode_buffs[Y_COMP][curr_depth][aux_abs_index], aux_num_part_in_cu*sizeof(et->intra_mode_buffs[Y_COMP][0][0]));
						memcpy(&et->cbf_buffs[Y_COMP][aux_depth][aux_abs_index], &et->cbf_buffs[Y_COMP][curr_depth][aux_abs_index], aux_num_part_in_cu*sizeof(et->cbf_buffs[Y_COMP][0][0]));
						memcpy(&et->tr_idx_buffs[aux_depth][aux_abs_index], &et->tr_idx_buffs[curr_depth][aux_abs_index], aux_num_part_in_cu*sizeof(et->tr_idx_buffs[0][0]));
					}
				}
			}
		}
	}
	PROFILER_ACCUMULATE(intra_luma_bucle3)

	curr_partition_info = &ctu->partition_list[et->partition_depth_start[depth]]+part_position;
//	curr_part_size_shift = et->max_cu_size_shift-curr_partition_info->depth;

//	if(et->rd_mode == 2)
//		return (curr_partition_info->cost+(bitcost_cu_mode/2)*et->rd.lambda+.5);//+curr_partition_info->size_chroma*et->rd.lambda+.5;
//	else
	return curr_partition_info->cost;
}


int motion_inter(henc_thread_t* et, ctu_info_t* ctu, int gcnt)
{
	double cost, cost_luma, cost_chroma, best_cost;
	int position = 0;
	int curr_depth = 0;
	PartSize part_size_type = SIZE_2Nx2N;
	wnd_t *quant_wnd = NULL, *decoded_wnd = NULL, *ref_wnd;
	ctu_info_t *ctu_rd = et->ctu_rd;
	cu_partition_info_t	*parent_part_info = NULL;
	cu_partition_info_t	*curr_partition_info = ctu->partition_list;
	int depth_state[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	uint cost_sum[MAX_PARTITION_DEPTH] = {0,0,0,0,0};
	int abs_index;
	int num_part_in_cu;
	int ll;

	//init rd auxiliar ctu
/*	if(et->rd_mode==1)
	{
		copy_ctu(ctu, ctu_rd);
		memset(ctu_rd->pred_mode, INTER_MODE, curr_partition_info->num_part_in_cu*sizeof(ctu_rd->pred_mode[0]));//indicamos que todas las codificaciones son intra
	}
*/
	while(curr_depth!=0 || depth_state[curr_depth]!=1)
	{
		curr_depth = curr_partition_info->depth;
		num_part_in_cu = curr_partition_info->num_part_in_cu;
		abs_index = curr_partition_info->abs_index;
		part_size_type = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//
		position = curr_partition_info->list_index - et->partition_depth_start[curr_depth];
		quant_wnd = &et->transform_quant_wnd[curr_depth];
		decoded_wnd = &et->decoded_mbs_wnd[curr_depth];
		
		cost_luma = cost_chroma = 0;

		if(curr_partition_info->is_b_inside_frame && curr_partition_info->is_r_inside_frame)//if br (and tl) are inside the frame, process
		{
			//encode
			PROFILER_RESET(inter_luma)
			cost_luma = encode_inter(et, ctu, gcnt, curr_depth, position, part_size_type);
			PROFILER_ACCUMULATE(inter_luma)

			cost_chroma = 0;
/*				if(part_size_type == SIZE_2Nx2N)
			{
				int iiiiii;
				PROFILER_RESET(inter_chroma)
				cost_chroma = encode_inter_chroma(et, ctu, gcnt, curr_depth, position, part_size_type);
				PROFILER_ACCUMULATE(inter_chroma)
				iiiiii=0;
			}
*/	
		}

		curr_partition_info->cost = cost_luma+cost_chroma;//+cost_bits*et->rd.sqrt_lambda;;

		depth_state[curr_depth]++;

		cost_sum[curr_depth]+=curr_partition_info->cost;

//#ifndef COMPUTE_AS_HM
		//if this matches, it is useless to continue the recursion. the case where part_size_type != SIZE_NxN is checked at the end of the consolidation buffer)
//		if(/*et->performance_mode == 0 && */part_size_type == SIZE_NxN && depth_state[curr_depth]!=4 && cost_sum[curr_depth] > parent_part_info->cost && ctu->partition_list[0].is_b_inside_frame && ctu->partition_list[0].is_r_inside_frame)//parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
//		{
//			depth_state[curr_depth]=4;
//		}

		//stop recursion
		if(et->performance_mode>0 && (curr_partition_info->recursive_split==0 || curr_partition_info->cost == 0) && /*curr_depth && */part_size_type != SIZE_NxN && curr_partition_info->is_b_inside_frame && curr_partition_info->is_r_inside_frame)
		{
			int max_processing_depth;// = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

			//if we fill this in here we don't have to consolidate
			memset(&ctu->pred_depth[abs_index], curr_depth-(part_size_type==SIZE_NxN), num_part_in_cu*sizeof(ctu->pred_depth[0]));
			memset(&ctu->part_size_type[abs_index], part_size_type, num_part_in_cu*sizeof(ctu->part_size_type[0]));
			if(et->rd_mode==1)//rd
			{
				memset(&ctu_rd->pred_depth[abs_index], curr_depth-(part_size_type==SIZE_NxN), num_part_in_cu*sizeof(ctu_rd->pred_depth[0]));
				memset(&ctu_rd->part_size_type[abs_index], part_size_type, num_part_in_cu*sizeof(ctu_rd->part_size_type[0]));
			}
			cost_chroma = 0;//only!=0 if part_size_type == SIZE_NxN

			max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);
			if(curr_depth+1 <= max_processing_depth)//el = es para cuando et->max_intra_tr_depth!=4
			{
				int aux_depth;
				cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
				abs_index = aux_partition_info->abs_index;
				num_part_in_cu  = aux_partition_info->num_part_in_cu;

				for(aux_depth=curr_depth+1;aux_depth<=max_processing_depth;aux_depth++)//+1 pq el nivel superior esta ya sincronizado y tenemos que sincronizar los siguientes
				{
//					synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);	
//					synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);
					//for rd
//					consolidate_recursive_info_buffers(et, gcnt, parent_part_info, curr_depth, aux_depth, abs_index, num_part_in_cu);
				}
//				synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
			}
		}
//#endif 

		if(curr_depth<et->max_pred_partition_depth && curr_partition_info->is_tl_inside_frame)//depth_state[curr_depth]!=4 is for fast skip//if tl is not inside the frame don't process the next depths
		{
			curr_depth++;
			parent_part_info = curr_partition_info;
		}
		else if(/*fast_skip || */depth_state[curr_depth]==4)//la depth =1 lo hemos consolidado antes del bucle
		{
			int max_processing_depth;
			cost_chroma = 0;

			if(cost_sum[curr_depth] < parent_part_info->cost &&  part_size_type == SIZE_NxN && (curr_partition_info->is_b_inside_frame && curr_partition_info->is_r_inside_frame))//esto se tiene que hacer aqui pq si sno se eligen cu_modes distintos para cada una de las 4 particiones
			{
//				position = parent_part_info->children[0]->list_index - et->partition_depth_start[curr_depth];
//				PROFILER_RESET(inter_chroma)
//				cost_chroma = encode_inter_chroma(et, ctu, gcnt, curr_depth, position, part_size_type);//encode_intra_chroma(et, gcnt, depth, position, part_size_type);	
//				PROFILER_ACCUMULATE(inter_chroma)
			}
			while(depth_state[curr_depth]==4 && curr_depth>0)//>0 pq consolidamos sobre el padre, 
			{
/*				cost_luma = parent_part_info->children[0]->cost + parent_part_info->children[1]->cost +parent_part_info->children[2]->cost+parent_part_info->children[3]->cost;

				depth_state[curr_depth] = 0;
				cost = cost_luma+cost_chroma;//cost_chroma=0 after the first iteration 
				cost_sum[curr_depth]+=cost_chroma;
				best_cost = parent_part_info->cost;

				abs_index = parent_part_info->abs_index;
				num_part_in_cu = parent_part_info->num_part_in_cu;

				//choose best
				if(cost<best_cost || !(curr_partition_info->is_b_inside_frame && curr_partition_info->is_r_inside_frame))//if we get here, tl should be inside the frame
				{
					//Aqui consolidamos los resultados bottom-up
					//here we consolidate the bottom-up results for being preferred to the top-down computation
					int part_size_type2 = (curr_depth<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//

					cost_sum[parent_part_info->depth] -= parent_part_info->cost;
					cost_sum[parent_part_info->depth] += cost;
					parent_part_info->cost = cost;
					quant_wnd = &et->transform_quant_wnd[curr_depth];//&et->transform_quant_wnd[curr_depth+1];
					decoded_wnd = &et->decoded_mbs_wnd[curr_depth];//&et->decoded_mbs_wnd[curr_depth+1];
			
					//consolidate in parent
					synchronize_motion_buffers_luma(et, parent_part_info, quant_wnd, &et->transform_quant_wnd[curr_depth-1], decoded_wnd, &et->decoded_mbs_wnd[curr_depth-1], gcnt);
					synchronize_motion_buffers_chroma(et, parent_part_info, quant_wnd, &et->transform_quant_wnd[curr_depth-1], decoded_wnd, &et->decoded_mbs_wnd[curr_depth-1], gcnt);

					consolidate_recursive_info_buffers(et, gcnt, parent_part_info, curr_depth, curr_depth-1, abs_index, num_part_in_cu);

					if(part_size_type2==SIZE_NxN)
					{
						int num_part_in_sub_cu = parent_part_info->children[0]->num_part_in_cu;
						cbf_split[Y_COMP] = (et->cbf_buffs[Y_COMP][curr_depth][abs_index]&2) || (et->cbf_buffs[Y_COMP][curr_depth][abs_index+num_part_in_sub_cu]&2) || (et->cbf_buffs[Y_COMP][curr_depth][abs_index+2*num_part_in_sub_cu]&2) || (et->cbf_buffs[Y_COMP][curr_depth][abs_index+3*num_part_in_sub_cu]&2);
						//consolidate cbf flags
						for(ll=abs_index;ll<abs_index+num_part_in_cu;ll++)
						{
							et->cbf_buffs[Y_COMP][curr_depth-1][ll] |= cbf_split[Y_COMP];
						}
					}

					if(curr_depth==et->max_pred_partition_depth)
					{
						//if we fill this in here we don't have to consolidate
						memset(&ctu->pred_depth[abs_index], curr_depth-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu->pred_depth[0]));
						memset(&ctu->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu->part_size_type[0]));
						if(et->rd_mode==1)
						{
							memset(&ctu_rd->pred_depth[abs_index], curr_depth-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu_rd->pred_depth[0]));
							memset(&ctu_rd->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu_rd->part_size_type[0]));
						}
					}
				}
				else
				{
					//Aqui prevalecen los resultados de la computacion top-down
					//top-down computation results are prefered
					int part_size_type2 = (curr_depth-1<et->max_pred_partition_depth)?SIZE_2Nx2N:SIZE_NxN;//

					//if we fill this in here we don't have to consolidate
					memset(&ctu->pred_depth[abs_index], curr_depth-1-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu->pred_depth[0]));
					memset(&ctu->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu->part_size_type[0]));
					if(et->rd_mode==1)//rd
					{
						memset(&ctu_rd->pred_depth[abs_index], curr_depth-1-(part_size_type2==SIZE_NxN), num_part_in_cu*sizeof(ctu_rd->pred_depth[0]));
						memset(&ctu_rd->part_size_type[abs_index], part_size_type2, num_part_in_cu*sizeof(ctu_rd->part_size_type[0]));
					}
				}
				cost_sum[curr_depth] = 0;
				curr_depth--;
				parent_part_info = parent_part_info->parent;

				if(et->performance_mode == 0 && curr_depth>0 && curr_depth<et->max_pred_partition_depth && depth_state[curr_depth]<4 && ctu->partition_list[0].is_b_inside_frame && ctu->partition_list[0].is_r_inside_frame)
				{
					double totalcost = 0;
					int h;
					for(h=0;h<depth_state[curr_depth];h++)
						totalcost += parent_part_info->children[h]->cost;

					//if(curr_depth>0 && curr_depth<et->max_pred_partition_depth && /*depth_state[curr_depth]!=4 && totalcost > parent_part_info->cost && et->partition_info[0].is_b_inside_frame && et->partition_info[0].is_r_inside_frame)//parent_part_info->is_b_inside_frame && parent_part_info->is_r_inside_frame)
					if(totalcost > parent_part_info->cost)//cost of childs is already higher than cost of parent - stop recursion
					{
						depth_state[curr_depth] = 4;
					}
				}
				cost_chroma = 0;//only!=0 if part_size_type == SIZE_NxN
				*/
			}

/*			max_processing_depth = min(et->max_pred_partition_depth+et->max_intra_tr_depth-1, MAX_PARTITION_DEPTH-1);

			if(curr_depth+1 <= max_processing_depth)//el = es para cuando et->max_intra_tr_depth!=4
			{
				int aux_depth;
				cu_partition_info_t*	aux_partition_info = (parent_part_info!=NULL)?parent_part_info->children[(depth_state[curr_depth]+3)&0x3]:&ctu->partition_list[0];
				abs_index = aux_partition_info->abs_index;
				num_part_in_cu  = aux_partition_info->num_part_in_cu;

				for(aux_depth=curr_depth+1;aux_depth<=max_processing_depth;aux_depth++)//+1 pq el nivel superior esta ya sincronizado y tenemos que sincronizar los siguientes
				{
					synchronize_reference_buffs(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);	
					synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[aux_depth], gcnt);
					//for rd
					consolidate_recursive_info_buffers(et, gcnt, parent_part_info, curr_depth, aux_depth, abs_index, num_part_in_cu);
				}
				synchronize_reference_buffs_chroma(et, aux_partition_info, &et->decoded_mbs_wnd[curr_depth], &et->decoded_mbs_wnd[NUM_DECODED_WNDS-1], gcnt);
			}
//			acc_cost = 0;
*/
		}

		if(parent_part_info!=NULL)
			curr_partition_info = parent_part_info->children[depth_state[curr_depth]];
	}
	
/*	curr_partition_info = &ctu->partition_list[0];
	abs_index = curr_partition_info->abs_index;
	curr_depth = curr_partition_info->depth;
	num_part_in_cu = curr_partition_info->num_part_in_cu;

	//esta informacion o parte de ella se necesita para la codificacion de otros ctus (por lo menos el intra_mode y el pred_mode. Creo que los cbfs y los tr_idx solo se necesitan los actuales, no los de ctus anteriores)
	memcpy(&ctu->cbf[Y_COMP][abs_index], &et->cbf_buffs[Y_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->cbf[Y_COMP][0]));
	memcpy(&ctu->cbf[U_COMP][abs_index], &et->cbf_buffs[U_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->cbf[U_COMP][0]));
	memcpy(&ctu->cbf[V_COMP][abs_index], &et->cbf_buffs[V_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->cbf[V_COMP][0]));
	memcpy(&ctu->intra_mode[Y_COMP][abs_index], &et->intra_mode_buffs[Y_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->intra_mode[Y_COMP][0]));
	memcpy(&ctu->intra_mode[CHR_COMP][abs_index], &et->intra_mode_buffs[CHR_COMP][curr_depth][abs_index], num_part_in_cu*sizeof(ctu->intra_mode[CHR_COMP][0]));
	memcpy(&ctu->tr_idx[abs_index], &et->tr_idx_buffs[curr_depth][abs_index], num_part_in_cu*sizeof(ctu->tr_idx[0]));

	memset(&ctu->pred_mode[abs_index], INTER_MODE, num_part_in_cu*sizeof(ctu->pred_mode[0]));//indicamos que todas las codificaciones son intra
	return best_cost;
*/
}

