HomerHEVC
=========

See www.homerhevc.com. 

HOMER (Hevc Open Mpeg EncodeR) is an open-source HEVC encoder to encode YUV420 video sequences to the new HEVC(H.265) stream format. 

It is published under the LPGLv2.1 license, and is therefore Free Software according to the Free Software Foundation.

Development is guided by three main aspects:
- Easy portability.
- Multiplatform (currently tested in Windows and Linux). 
- Performance.

Code style of the development is C'99 and recursive functions have been implemented in a sequential manner to avoid the drawback of recursive calls to complex functions and to ease portability. OS dependant code is isolated in a very simple file, while low level optimizations are handled with an interface for low level functions. 

HomerHEVC is still under development and will improve in quality and performance during the development.

Current Features (latest version is HomerHEVC_v1.2)
--------------------------------------------------------
- Multiplatform (Linux,Windows)
- 8 bit-depth.
- Intra and Baseline profile (I and P images with 1 reference image).
- All intra prediction modes.
- 2Nx2N and NxN inter prediction modes.
- Half pixel and quarter pixel precission motion estimation.
- CTU size 64 (Intra mode supports CTU size 64, 32 and 16).
- All transform sizes (32,16,8,4).
- Fixed QP, CBR (Constant Bit Rate), and VBR (Variable Bit Rate).
- Deblocking filter.
- Wpp and Frame based parallelization for massive parallelism
- Sign hiding bit enabled.
- intra RDO.
- intra-inter fast RD.

Optimizations (SSE42):
- Intra prediction generation.
- Motion estimation.
- Inter prediction.
- Intra Prediction.
- Reconstruction. 
- SAD, SSD.
- Transforms.
- Quantization.


Version Releases and Expected Evolution
----------------------------------------
A new version (v1.0, v2.0, ...) is expected to be released every 4 months including new features and further optimizations. 
Reported bugs and small improvements of the already published versions will be published as soon as they are developed, being signaled by increasing the subversion index (eg. v0.1, v0.2, v0.3...).


Roadmap
-------
Current Version improvements - homerHEVC_V1.0 (November 2014)
- Inter prediction.
- Rate control
- Further SSE/AVX optimizations & quality improvements.

homerHEVC_V2.0 (April 2015)
- SAO.
- Frame based parallelism.
- Half pixel and quarter pixel precission motion estimation.
- Further SSE/AVX optimizations & quality improvements.

homerHEVC_V3.0 (July 2015)
- B frames.

Donwnload
---------
It is recomended to download the latest release from the Release Tab at: https://github.com/jcasal-homer/HomerHEVC/releases


Compiling homerHEVC
-------------------

homerHEVC is composed of a simple example application (homer_app) and the encoder library (homer_lib).

Compile in Windows 
- Visual C++ 2008. Open the solution inside the "build\vc9\" folder and rebuild. Default project is for 64 bits.
- Visual C++ 2012. Open the solution inside the "build\vc11\" folder and rebuild. Default project is for 64 bits.

Compile in Linux
- Open a console. Go to the "build/Linux/" folder and type "make clean all" to execute the makefile.


Configuration:
--------------

These are a list of the configuration variables currently supported:

    call: 

    homer_app [-option] [value]...

    options:

	homer_app [-option] [value]...
	options:
	-h:                             help
	-i:                             input yuv file
	-o:                             output 265 file
	-o-raw:                         output raw frames in yuv format
	-widthxheight:                  default = 1280x720
	-frame_rate:                    default = 25 fps
	-cu_size:                       cu size [16,32 or 64], default = 64 (only 64 supported for inter prediction)
	-intra_preriod:                 [0-...], 0=infinite, default = 100
	-gop_size:                      0:intra profile, 1: IPPP.. profile, default = 1
	-num_ref_frame:                 default = 1 (only 1 reference currently supported)
	-qp:                            qp[0-51], default = 32
	-motion_estimation_precision    0=pel, 1=half_pel, 2=quarter_pel, default = 2
	-chroma_qp_offset:              chroma_qp_offset[-12,12], default = 2
	-n_enc_engines:                 number of encoder engines encoding frames in parallel, default = 2
	-n_wpp_threads:                 0:no wpp, >0-number of wpp threads, default = 10
	-max_pred_depth:                [0-4], default = 4
	-max_intra_tr_depth:            [0-4], default = 2
	-max_inter_tr_depth:            [0-4], default = 1
	-sign_hiding:                   0=off, 1=on, default = 1
	-bitrate_mode:                  0=fixed qp, 1=CBR (Constant bitrate), 2=VBR (Variable bitrate), default = VBR
	-bitrate:                       in kbps when bitrate_mode=CBR, default = 1250
	-vbv_size:                      in kbps when bitrate_mode=CBR, default = 1.*bitrate
	-vbv_init:                      in kbps when bitrate_mode=CBR, default = .25*vbv_size
	-performance_mode:              0=full computation, 1=fast , 2= ultra fast, default = ufast
	-scene_change:                  0=do not reinit, 1=reinit gop on scene change, default = 1;
	-rd:                            0=off, 1=full rd (only in intra) , 2= fast rd, default = fast
	-n_frames:                      default = 1000
	-skipped_frames:                default = 0

examples:

	intra:
	homer_app -i /home/juan/Patrones/720p5994_parkrun_ter.yuv -o output0.265 -widthxheight 1280x720 -frame_rate 50 -intra_period 1 -gop_size 0 -max_pred_depth 4 -max_intra_tr_depth 3 -bitrate 25000 -vbv_size 1000 -vbv_init 1000 -n_wpp_threads 10 -performance_mode 1 -rd_mode 2 -n_frames 400

	inter:
	homer_app -i /home/juan/Patrones/720p5994_parkrun_ter.yuv -o output0.265 -widthxheight 1280x720 -frame_rate 50 -intra_period 100 -gop_size 1 -max_pred_depth 4 -max_intra_tr_depth 3 -max_inter_tr_depth 1 -bitrate 5000 -vbv_size 2500 -vbv_init 750 -n_wpp_threads 10 -performance_mode 1 -rd_mode 2 -n_frames 400

Contribute:
--------------
We encourage developers and users to participate in the project as contributors by: developing new features, reporting bugs, and giving feedback.

HomerHEVC Contribution License Agreement (CLA) must be signed before starting contributing.

if you would like to contribute, please write to jcasal@homerhevc.com.


More Info:
----------
HomerHEVC is led by Juan Casal, an IT and Telecomunications engineer with more than 12 years of experience in video encoding and transmission.

Find Juan Casal's personal linkedIn profile in https://www.linkedin.com/pub/juan-casal/47/373/8a6.

Contact: jcasal@homerhevc.com

Web page: www.homerhevc.com

Please report bugs and different issues in https://github.com/jcasal-homer/HomerHEVC/issues
