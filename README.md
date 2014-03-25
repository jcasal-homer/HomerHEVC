HomerHEVC:

HOMER (Hevc Open Mpeg EncodeR) is an open-source HEVC encoder, to encode YUV420 video sequences to the new HEVC(H.265) stream format.
It is published under the LPGLv2.1 license.
Homer is still under development and will improve in quality and performance while being developed.

Development is guided by three main aspects:
- Easy portability.
- Multiplatform (currently tested in Windows and Linux). 
- Performance.

Thus, code style of the development is C'99 and recursive functions have been implemented in a sequential manner to avoid the drawback of recursive calls to complex functions and to ease portability. OS dependant code is isolated in a very simple file, while low level optimizations are handled with an interface for low level functions. 



Current Features: 

- 8 bit-depth.
- All intra prediction modes.
- All CTU sizes (64, 32, 16).
- All transform sizes (32,16,8,4).
- Sign hiding bit enabled.
- RDO.
- Wpp parallelization.
- Fast decision mode algorithm.

Current SSE42 optimizations (intrinsics):
- Intra prediction generation.
- Prediction.
- Reconstruction. 
- SAD, SSD.
- Transforms.
- Quantization.




Version Releases and Expected Evolution:

A new version (v1.0, v2.0, ...) is expected to be released every 4 months including new features and further optimizations. 
Reported bugs and small improvements of the already published versions will be published as soon as they are developed, being signaled by increasing the subversion index (eg. v0.1, v0.2, v0.3...).


Roadmap:

HomerHEVC_V0.1 (April 2014)
- Intra smooth filter (To be added)
- Intra skipped (To be added)

HomerHEVC_V1.0 (July 2014)
- Inter prediction including interpolation filters.
- Rate control (CBR, VBR, Capped VBR)
- Further SSE/AVX optimizations & quality improvements.

HomerHEVC_V2.0 (November 2014)
- SAO.
- Deblocking.
- Further SSE/AVX optimizations & quality improvements.



Compiling Homer:

HomerHEVC is composed of a simple example application (homer_app) and the encoder library (homer_lib).

Compile in Windows 
- Visual C++ 2008. Open the solution inside the "build\vc9\" folder and rebuild. Default project is for 64 bits.
- Visual C++ 2012. Open the solution inside the "build\vc11\" folder and rebuild. Default project is for 64 bits.
- Cygwin. Open the Cygwin console, go to the "build/Cygwin/" folder and type "make clean all" to execute the makefile.

Compile in Linux
- Open a console. Go to the "build/Linux/" folder and type "make clean all" to execute the makefile.



Configuration:

These are a list of the configuration variables currently supported:

call: homer_app [-option] [value]...

options:


-h:					help


-i:					input yuv file


-o:					output 265 file


-widthxheight:           		frame resolution, default = 1280x720
-cu_size:                		cu size[16,32 or 64], default = 64


-qp:                     		fixed qp[0-51], default = 32


-n_wpp_threads:          		0:no wpp, >0:number of wpp threads, default = 1


-max_intra_pred_depth:   	    	[0-4], default = 4


-max_intra_tr_depth:     	    	[0-4], default = 4


-sign_hiding:            		0=off, 1=on, default = 1


-performance_mode:       	    	0=full computation, 1=fast , 2= ultra fast


-rd:                     		0=off, 1=full rd , 2= fast rd

Configuration examples:

homer_app -i /home/juan/Patrones/720p5994_parkrun_ter.yuv -o output0.265  -widthxheight 1280x720 -n_wpp_threads 10 -performance_mode 2 -rd_mode 2

homer_app -i /home/juan/Patrones/720p5994_parkrun_ter.yuv -o output0.265 -widthxheight 1280x720 -n_wpp_threads 10 -performance_mode 1 -rd_mode 2

homer_app -i /home/juan/Patrones/720p5994_parkrun_ter.yuv -o output0.265 -widthxheight 1280x720 -n_wpp_threads 10 -performance_mode 1 -rd_mode 1





Contributors:

Up to now all code has being developed by Juan Casal, an IT and Telecomunications engineer with more than 12 years of experience in video encoding and transmission.
Please fin Juan Casal's personal linkedIn profile in here.

If you want to contribute to the development of HomerHEVC, please contact me on jcasal.homer@gmail.com.

