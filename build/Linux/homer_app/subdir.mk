################################################################################
# subdir.mk for homer_app/Release
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 

LIBS := -lhomer_lib -lpthread -lm

C_SRCS += \
../../../../src/homer_app/homer_app.c 

OBJS += \
./homer_app.o 

C_DEPS += \
./homer_app.d 


