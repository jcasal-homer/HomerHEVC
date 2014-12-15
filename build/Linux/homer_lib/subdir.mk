################################################################################
# subdir.mk for homer_lib
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../../../../src/homer_lib/hmr_arithmetic_encoding.c \
../../../../src/homer_lib/hmr_binary_encoding.c \
../../../../src/homer_lib/hmr_bitstream.c \
../../../../src/homer_lib/hmr_container.c \
../../../../src/homer_lib/hmr_encoder_lib.c \
../../../../src/homer_lib/hmr_headers.c \
../../../../src/homer_lib/hmr_mem_transfer.c \
../../../../src/homer_lib/hmr_metics.c \
../../../../src/homer_lib/hmr_motion_intra.c \
../../../../src/homer_lib/hmr_motion_intra_chroma.c \
../../../../src/homer_lib/hmr_motion_inter.c \
../../../../src/homer_lib/hmr_rate_control.c \
../../../../src/homer_lib/hmr_deblocking_filter.c \
../../../../src/homer_lib/hmr_profiler.c \
../../../../src/homer_lib/hmr_quant.c \
../../../../src/homer_lib/hmr_sse42_functions_pixel.c \
../../../../src/homer_lib/hmr_sse42_functions_prediction.c \
../../../../src/homer_lib/hmr_sse42_functions_quant.c \
../../../../src/homer_lib/hmr_sse42_functions_transform.c \
../../../../src/homer_lib/hmr_sse42_functions_inter_prediction.c \
../../../../src/homer_lib/hmr_tables.c \
../../../../src/homer_lib/hmr_transform.c \

OBJS += \
./hmr_arithmetic_encoding.o \
./hmr_binary_encoding.o \
./hmr_bitstream.o \
./hmr_container.o \
./hmr_encoder_lib.o \
./hmr_headers.o \
./hmr_mem_transfer.o \
./hmr_metics.o \
./hmr_motion_intra.o \
./hmr_motion_intra_chroma.o \
./hmr_motion_inter.o \
./hmr_rate_control.o \
./hmr_deblocking_filter.o \
./hmr_profiler.o \
./hmr_quant.o \
./hmr_sse42_functions_pixel.o \
./hmr_sse42_functions_prediction.o \
./hmr_sse42_functions_quant.o \
./hmr_sse42_functions_transform.o \
./hmr_sse42_functions_inter_prediction.o \
./hmr_tables.o \
./hmr_transform.o \


C_DEPS += \
./hmr_arithmetic_encoding.d \
./hmr_binary_encoding.d \
./hmr_bitstream.d \
./hmr_container.d \
./hmr_encoder_lib.d \
./hmr_headers.d \
./hmr_mem_transfer.d \
./hmr_metics.d \
./hmr_motion_intra.d \
./hmr_motion_intra_chroma.d \
./hmr_motion_inter.d \
./hmr_rate_control.d \
./hmr_deblocking_filter.d \
./hmr_profiler.d \
./hmr_quant.d \
./hmr_sse42_functions_pixel.d \
./hmr_sse42_functions_prediction.d \
./hmr_sse42_functions_quant.d \
./hmr_sse42_functions_transform.d \
./hmr_sse42_functions_inter_prediction.d \
./hmr_tables.d \
./hmr_transform.d \



