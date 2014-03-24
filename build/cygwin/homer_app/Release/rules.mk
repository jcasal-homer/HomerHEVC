################################################################################
# rules.mk for homer_app/Release
################################################################################

# rules for building homer_app/Release sources 
%.o: ../../../../src/homer_app/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -I../../../../src/homer_lib/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


