################################################################################
# rules.mk for homer_app/Debug
################################################################################

# rules for building homer_app/Debug sources
%.o: ../../../../src/homer_app/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -I../../../../src/homer_lib/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


