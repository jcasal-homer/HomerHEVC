/*****************************************************************************
* homer_app.c : homerHEVC example application
******************************************************************************
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
#include <stdlib.h>
//#include <process.h>
#include <malloc.h>
#include <math.h>
#include <homer_hevc_enc_api.h>
#include <string.h>

#ifdef _MSC_VER
#include <Windows.h>
#include <io.h>
#define fseek_64 _fseeki64
#else
#include <sys/time.h>
#define fseek_64 fseek
#endif

//#define FILE_IN  "C:\\PruebasIntegracionCiresXXI\\DroneChateau.yuv"
#define FILE_IN  "C:\\Patrones\\TestEBU720p50.yuv"//TestEBU1080i50_yuv.YUV"//LolaTest420.yuv"//BrazilianDancer.yuv"//DebugPattern_384x256.yuv"//Flags.yuv"//"C:\\PruebasCiresXXI\\Robots.yuv"//TestEBU720p50_synthetic.yuv"//sinthetic_freeze.yuv"//720p5994_parkrun_ter.yuv"
//#define FILE_IN  "C:\\Patrones\\DebugPattern_384x256.yuv"//demo_pattern_192x128.yuv"//table_tennis_420.yuv"//LolaTest420.yuv"//demo_pattern_192x128.yuv"//"C:\\Patrones\\DebugPattern_248x184.yuv"//"C:\\Patrones\\DebugPattern_384x256.yuv"//DebugPattern_208x144.yuv"//Prueba2_deblock_192x128.yuv"//demo_pattern_192x128.yuv"
//#define FILE_IN  "C:\\Patrones\\LolaTest420.yuv"
//#define FILE_IN  "C:\\Patrones\\1080p_pedestrian_area.yuv"
//#define FILE_IN  "C:\\Patrones\\DebugPattern_248x184.yuv"

#define FILE_OUT	"C:\\Patrones\\homer_development.265"//Flags.265"//"C:\\PruebasCiresXXI\\Robots.265"//Flags_zeros_3.265"//output_Homer_synthetic_full_HM_prueba.265"//DebugPattern_248x184.265"//
#define FILE_REF	"C:\\Patrones\\refs_Homer.bin"


#define HOR_SIZE	1280//624//192//(208)//(384+16)//1280//1920//1280//(2*192)//1280//720//(2*192)//(192+16)//720//320//720
#define VER_SIZE	720//352//128//(144)//(256+16)//720//1080//720//(2*128)//720//576//(2*128)//(128+16)//320//576
#define FPS			50//24//25//50


#ifdef _MSC_VER
unsigned int get_ms()
{
	return GetTickCount();
}
#else
unsigned int get_ms()
{
	struct timeval tv;
	if(gettimeofday(&tv, NULL) != 0)
			return 0;

	return (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
}
#endif


char file_in_name[256] = "";
char file_out_name[256] = "";
char file_ref_name[256] = "";

void print_help()
{
	printf("\r\nhomer_app [-option] [value]...\r\n");

	printf("options:\r\n");
	printf("-h: \t\t\t\t\t help\r\n");
	printf("-i: \t\t\t\t\t input yuv file\r\n");
	printf("-o: \t\t\t\t\t output 265 file\r\n");
	printf("-o-raw: \t\t\t\t output raw frames in yuv format\r\n");
	printf("-widthxheight: \t\t\t\t default = 1280x720\r\n");
	printf("-frame_rate: \t\t\t\t default = 50 fps\r\n");	
	printf("-cu_size: \t\t\t\t cu size [16,32 or 64], default = 64 (only 64 supported for inter prediction)\r\n");
	printf("-intra_period: \t\t\t 0=infinite, default = 30 \r\n");
	printf("-gop_size: \t\t\t\t 0:intra profile, 1: IPPP.. profile, default = 1\r\n");
	printf("-num_ref_frame: \t\t\t default = 1 (only 1 reference currently supported) \r\n");	
	printf("-qp: \t\t\t\t\t qp[0-51], default = 32\r\n");
	printf("-motion_estimation_precision: \t\t 0=pel, 1=half_pel, 2=quarter_pel, default = 2\r\n");
	printf("-chroma_qp_offset: \t\t\t chroma_qp_offset[-12,12], default = 2\r\n");	
	printf("-n_wpp_threads: \t\t\t 0:no wpp, >0-number of wpp threads, default = 10\r\n");	
	printf("-max_pred_depth: \t\t\t [0-4], default = 4\r\n");
	printf("-max_intra_tr_depth: \t\t\t [0-4], default = 2\r\n");
	printf("-max_inter_tr_depth: \t\t\t [0-4], default = 1\r\n");
	printf("-sign_hiding: \t\t\t\t 0=off, 1=on, default = 1\r\n");
	printf("-bitrate_mode: \t\t\t\t 0=fixed qp, 1=CBR (Constant bitrate), default = CBR\r\n");
	printf("-bitrate: \t\t\t\t in kbps when bitrate_mode=CBR, default = 5000\r\n");
	printf("-vbv_size: \t\t\t\t in kbps when bitrate_mode=CBR, default = .5*bitrate\r\n");
	printf("-vbv_init: \t\t\t\t in kbps when bitrate_mode=CBR, default = .1*bitrate\r\n");
	printf("-performance_mode: \t\t\t 0=full computation, 1=fast , 2= ultra fast, default = fast\r\n");
	printf("-scene_change: \t\t\t\t 0=do not reinit, 1=reinit gop on scene change, default = 1\r\n");
	printf("-rd: \t\t\t\t\t 0=off, 1=full rd (only in intra) , 2= fast rd, default = fast\r\n");
	printf("-n_frames: \t\t\t\t default = 15\r\n");
	printf("-skipped_frames: \t\t\t default = 0\r\n");

	printf("\r\nexamples:\r\n\r\n");
	printf("intra:\r\n");
	printf("homer_app -i /home/juan/Patrones/720p5994_parkrun_ter.yuv -o output0.265 -widthxheight 1280x720 -frame_rate 50 -intra_period 1 -gop_size 0 -max_pred_depth 4 -max_intra_tr_depth 3 -bitrate 25000 -vbv_size 1000 -vbv_init 1000 -n_wpp_threads 10 -performance_mode 1 -rd_mode 2 -n_frames 400\r\n\r\n");

	printf("inter:\r\n");
	printf("homer_app -i /home/juan/Patrones/720p5994_parkrun_ter.yuv -o output0.265 -widthxheight 1280x720 -frame_rate 50 -intra_period 100 -gop_size 1 -max_pred_depth 4 -max_intra_tr_depth 3 -max_inter_tr_depth 1 -bitrate 5000 -vbv_size 2500 -vbv_init 750 -n_wpp_threads 10 -performance_mode 1 -rd_mode 2 -n_frames 400\r\n\r\n");
}




void parse_args(int argc, char* argv[], HVENC_Cfg *cfg, int *num_frames, int *skipped_frames)
{
	int args_parsed = 1;

/*	if(argc==1)
	{
		printf ("\r\nno args passed!\r\ntype -h for help\r\n");
		exit(0);
	}
*/
	while(args_parsed<argc)
	{
		if(strcmp(argv[args_parsed] , "-h")==0)//input
		{
			print_help();
			exit(0);
		}
		else if(strcmp(argv[args_parsed] , "-i")==0 && args_parsed+1<argc)//input
		{
			args_parsed++;
			strcpy(file_in_name, argv[args_parsed++]);
		}
		else if(strcmp(argv[args_parsed], "-o")==0 && args_parsed+1<argc)//output
		{
			args_parsed++;
			strcpy(file_out_name, argv[args_parsed++]);
		}
		else if(strcmp(argv[args_parsed], "-o-raw")==0 && args_parsed+1<argc)//output
		{
			args_parsed++;
			strcpy(file_ref_name, argv[args_parsed++]);
		}
		else if(strcmp(argv[args_parsed], "-widthxheight")==0 && args_parsed+1<argc)//720x576, 1280x720, 1920x1080.... Multiple of 16
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%dx%d", &cfg->width, &cfg->height);
		}
		else if(strcmp(argv[args_parsed], "-frame_rate")==0 && args_parsed+1<argc)//frames per second of video
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%f", &cfg->frame_rate);
		}
		else if(strcmp(argv[args_parsed], "-cu_size")==0 && args_parsed+1<argc)//cu_size: 64, 32, 16
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->cu_size);
		}
		else if(strcmp(argv[args_parsed], "-max_pred_depth")==0 && args_parsed+1<argc)//depth of prediction, default 4
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->max_pred_partition_depth);
		}
		else if(strcmp(argv[args_parsed], "-max_intra_tr_depth")==0 && args_parsed+1<argc)//transform of intra prediction, default 2
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->max_intra_tr_depth);
		}
		else if(strcmp(argv[args_parsed], "-max_inter_tr_depth")==0 && args_parsed+1<argc)//transform of inter prediction, default 1
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->max_inter_tr_depth);
		}
		else if(strcmp(argv[args_parsed], "-intra_period")==0 && args_parsed+1<argc)//Period between two I images
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->intra_period);
		}
		else if(strcmp(argv[args_parsed], "-gop_size")==0 && args_parsed+1<argc)//Number of P and B frames repeated in a video sequence
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->gop_size);
		}
		else if(strcmp(argv[args_parsed], "-num_ref_frames")==0 && args_parsed+1<argc)//Number of reference frames to be used as reference for inter prediction
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->num_ref_frames);
		}
		else if(strcmp(argv[args_parsed], "-qp")==0 && args_parsed+1<argc)//quant
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->qp );
		}
		else if(strcmp(argv[args_parsed], "-motion_estimation_precision")==0 && args_parsed+1<argc)//Number of reference frames to be used as reference for inter prediction
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->motion_estimation_precision);
		}
		else if(strcmp(argv[args_parsed], "-chroma_qp_offset")==0 && args_parsed+1<argc)//qp offset to apply for chroma
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->chroma_qp_offset);
		}	
		else if(strcmp(argv[args_parsed], "-n_wpp_threads")==0 && args_parsed+1<argc)//number of threads
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->wfpp_num_threads );
			if(cfg->wfpp_num_threads==0)
				cfg->wfpp_enable = 0;
			else
				cfg->wfpp_enable = 1;
		}
		else if(strcmp(argv[args_parsed], "-sign_hiding")==0 && args_parsed+1<argc)//sign_hiding, default 1
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->sign_hiding);
		}
		else if(strcmp(argv[args_parsed], "-bitrate_mode")==0 && args_parsed+1<argc)//bitrate_mode: 0=BR_FIXED_QP, 1 = BR_CBR (vbv based constant bit rate)
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->bitrate_mode);
		}
		else if(strcmp(argv[args_parsed], "-bitrate")==0 && args_parsed+1<argc)//stream bitrate in kbps when bitrate_mode=BR_CBR 
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->bitrate);
		}
		else if(strcmp(argv[args_parsed], "-vbv_size")==0 && args_parsed+1<argc)//size of VBV (video buffering verifier) in kbps when bitrate_mode=BR_CBR 
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->vbv_size);
		}
		else if(strcmp(argv[args_parsed], "-vbv_init")==0 && args_parsed+1<argc)//initial fullnes of VBV in kbps when bitrate_mode=BR_CBR 
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->vbv_init);
		}
		else if(strcmp(argv[args_parsed], "-performance_mode")==0 && args_parsed+1<argc)//performance_mode:	//0=FULL_COMPUTATION, 1=PERF_FAST_COMPUTATION, 2=FAST_PRED_MODE
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->performance_mode);
		}
		else if(strcmp(argv[args_parsed], "-scene_change")==0 && args_parsed+1<argc)//-scene_change: //0:do not reinit gop on scene change, 1:reinit gop on scene change
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->reinit_gop_on_scene_change);
		}
		else if(strcmp(argv[args_parsed], "-rd_mode")==0 && args_parsed+1<argc)//rd_mode: 0=OFF(DISTORTION_ONLY), 1=FULL_RD, 2=FAST_RD
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", &cfg->rd_mode);
		}
		else if(strcmp(argv[args_parsed], "-n_frames")==0 && args_parsed+1<argc)//number of frames to encode
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", num_frames);
		}
		else if(strcmp(argv[args_parsed], "-skipped_frames")==0 && args_parsed+1<argc)//frames to skip before starting encoding
		{
			args_parsed++;
			sscanf( argv[args_parsed++], "%d", skipped_frames);
		}
		else	//arg not recognized. Continue to next
		{
			printf("\r\nunrecognized argument: %s\r\n",argv[args_parsed]);
			args_parsed++;
		}
	}
}



int main (int argc, char **argv)
{
	int totalbits=0;
	unsigned int msInit=0, msTotal=0;
	int bCoding = 1;
	int bytes_read = 0;
	int frames_read = 0, encoded_frames = 0;
	FILE *infile = NULL, *outfile = NULL, *reffile = NULL;
	int skipped_frames = 00;//2075;//400+1575+25;//25;//1050;//800;//200;//0;
	int num_frames = 2000;//1500;//500;//2200;//100;//700;//15;

	unsigned char *frame[3];
	stream_t stream;
	encoder_in_out_t input_frame, output_stream, output_frame;

	void *pEncoder;

	nalu_t  *nalu_out[8];
	unsigned int num_nalus = 8;

	HVENC_Cfg	HmrCfg;

	strcpy(file_in_name, FILE_IN);
	strcpy(file_out_name, FILE_OUT);
//	strcpy(file_ref_name, FILE_REF);

	HmrCfg.size = sizeof(HmrCfg);
	HmrCfg.width = HOR_SIZE;
	HmrCfg.height = VER_SIZE;
	HmrCfg.profile = PROFILE_MAIN;
	HmrCfg.intra_period = 100;//1;
	HmrCfg.gop_size = 1;//0;
	HmrCfg.motion_estimation_precision = QUARTER_PEL;//QUARTER_PEL;//HALF_PEL;//PEL;//
	HmrCfg.qp = 32;//32;
	HmrCfg.frame_rate = FPS;
	HmrCfg.num_ref_frames = 1;
	HmrCfg.cu_size = 64;
	HmrCfg.max_pred_partition_depth = 4;
	HmrCfg.max_intra_tr_depth = 2;
	HmrCfg.max_inter_tr_depth = 1;
	HmrCfg.wfpp_enable = 1;
	HmrCfg.wfpp_num_threads = 10;
	HmrCfg.sign_hiding = 1;
	HmrCfg.rd_mode = RD_FAST;	  //0 no rd, 1 similar to HM, 2 fast
	HmrCfg.bitrate_mode = BR_VBR;//BR_CBR;//BR_FIXED_QP;//BR_FIXED_QP;//BR_FIXED_QP;//0=fixed qp, 1=cbr (constant bit rate)
	HmrCfg.bitrate = 5000;//in kbps
	HmrCfg.vbv_size = HmrCfg.bitrate*1.;//in kbps
	HmrCfg.vbv_init = HmrCfg.bitrate*0.5;//in kbps
	HmrCfg.chroma_qp_offset = 2;
	HmrCfg.reinit_gop_on_scene_change = 0;
	HmrCfg.performance_mode = PERF_UFAST_COMPUTATION;//PERF_FULL_COMPUTATION;//0=PERF_FULL_COMPUTATION (HM), 1=PERF_FAST_COMPUTATION (rd=1 or rd=2), 2=PERF_UFAST_COMPUTATION (rd=2)

	parse_args(argc, argv, &HmrCfg, &num_frames, &skipped_frames);

	if(!(infile = fopen(file_in_name, "rb")))
	{
		printf("Error opening input file: %s\r\n", file_in_name);
		exit(0);
	}

	if((outfile = fopen(file_out_name, "wb")) == NULL)
	{
		printf("Error opening output file: %s\r\n", file_out_name);
		exit(0);
	}

	if((strlen(file_ref_name)>0) && !(reffile = fopen(file_ref_name, "wb")))
	{
		printf("Error opening raw output file: %s\r\n", file_ref_name);
//		exit(0);
	}

	memset(&input_frame, 0, sizeof(input_frame));
	memset(&output_stream, 0, sizeof(output_stream));
	memset(&output_frame, 0, sizeof(output_frame));
	memset(&stream, 0, sizeof(stream));

	stream.streams[0] = (unsigned char *)calloc(HmrCfg.width*HmrCfg.height, 1);
	stream.streams[1] = (unsigned char *)calloc(HmrCfg.width*HmrCfg.height>>1,1);
	stream.streams[2] = (unsigned char *)calloc(HmrCfg.width*HmrCfg.height>>1,1);
	output_stream.stream.streams[0] = (unsigned char *)calloc(0x8000000,1);

	output_frame.stream.streams[0] = NULL;
	output_frame.stream.streams[1] = NULL;
	output_frame.stream.streams[2] = NULL;

	if(strlen(file_ref_name)>0)//if this is not allocated, the internal copy is not done
	{
		output_frame.stream.streams[0] = (unsigned char *)calloc(HmrCfg.width*HmrCfg.height, 1);
		output_frame.stream.streams[1] = (unsigned char *)calloc(HmrCfg.width*HmrCfg.height>>1,1);
		output_frame.stream.streams[2] = (unsigned char *)calloc(HmrCfg.width*HmrCfg.height>>1,1);
	}

	pEncoder = HOMER_enc_init();

	if(!HOMER_enc_control(pEncoder,HENC_SETCFG,&HmrCfg))
		return -1;

	msInit = get_ms();
	fseek_64(infile, 0, SEEK_SET);
	while(bCoding)
	{
		int frame_size = (int) (HmrCfg.width*HmrCfg.height*1.5);
//		fpos_t pos = frame_size*skipped_frames;
//		long long frame_size = ;
		if(frames_read<skipped_frames)
		{
			fseek_64(infile, frame_size, SEEK_CUR); // move to first frame
			frames_read++;
			msInit = get_ms();
			continue;
		}
		frame[0] = (unsigned char*)stream.streams[0];
		frame[1] = (unsigned char*)stream.streams[1];
		frame[2] = (unsigned char*)stream.streams[2];

		bytes_read = fread(frame[0],HmrCfg.width,HmrCfg.height,infile)*HmrCfg.width;
		bytes_read += fread(frame[1],HmrCfg.width>>1,HmrCfg.height>>1,infile)*(HmrCfg.width>>1);
		bytes_read += fread(frame[2],HmrCfg.width>>1,HmrCfg.height>>1,infile)*(HmrCfg.width>>1);

		stream.data_stride[0] = HmrCfg.width;
		stream.data_stride[1] = stream.data_stride[2]  = stream.data_stride[0]/2;

		input_frame.stream = stream;
		input_frame.pts = frames_read-skipped_frames;
		input_frame.image_type = IMAGE_AUTO;

		if(bCoding)
		{
			int bEndOfFile = bytes_read!=frame_size;
			if(!bEndOfFile)//if EOF is reached don´t try to encode
			{
				num_nalus = 8;
				HOMER_enc_encode(pEncoder, &input_frame);//, nalu_out, &num_nalus);
				Sleep(1000/(2*HmrCfg.frame_rate));
//				encoder_thread(pEncoder);
				frames_read++;			
			}
			else
			{
				int iiiii=0;
			}


			HOMER_enc_get_coded_frame(pEncoder, &output_frame, nalu_out, &num_nalus);

			if(num_nalus>0)
			{
				HOMER_enc_write_annex_b_output(nalu_out, num_nalus, &output_stream);
				fwrite(output_stream.stream.streams[0], sizeof(unsigned char), output_stream.stream.data_size[0], outfile);
				fflush(outfile);
				totalbits+=output_stream.stream.data_size[0];

				if(reffile!=NULL)
				{
					fwrite(output_frame.stream.streams[0], HmrCfg.width, HmrCfg.height, reffile); 
					fwrite(output_frame.stream.streams[1], HmrCfg.width>>1, HmrCfg.height>>1, reffile); 
					fwrite(output_frame.stream.streams[2], HmrCfg.width>>1, HmrCfg.height>>1, reffile); 
				}

				encoded_frames++;
			}
			if(encoded_frames == num_frames || (bEndOfFile && encoded_frames==frames_read-skipped_frames))
			{
				bCoding = 0;
			}
		}

		if(!bCoding)
		{
//			HOMER_enc_close(pEncoder);
			msTotal += get_ms()-msInit;
			printf("\r\n%d frames in %d milliseconds: %f fps", encoded_frames, msTotal, 1000.0*(encoded_frames)/(double)msTotal);

			break;
		}
	}

//	printf("\r\npulse una tecla\r\n");
//	getchar();


	free(stream.streams[0]);
	free(stream.streams[1]);
	free(stream.streams[2]);
	free(output_stream.stream.streams[0]);

	if(strlen(file_ref_name)>0)//if this is not allocated, the internal copy is not done
	{
		free(output_frame.stream.streams[0]);
		free(output_frame.stream.streams[1]);
		free(output_frame.stream.streams[2]);
	}

	fclose(infile);
	fclose(outfile);
	if(reffile!=NULL)
		fclose(reffile);
	return 0;
}
