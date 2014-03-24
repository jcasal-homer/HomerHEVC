############################################################################################
#	makefile for homer_app/Debug
############################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include ../subdir.mk
-include rules.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: homer_app

# Tool invocations
homer_app: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -L../../homer_lib/Debug -o "homer_app" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) homer_app
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
