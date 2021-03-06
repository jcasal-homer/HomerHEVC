################################################################################
# rules.mk for homer_lib/Release
################################################################################

# rules for building sources of homer_lib/Release
./%.o: ../../../../src/homer_lib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -msse -msse2 -mssse3 -msse4 -msse4.1 -msse4.2 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


