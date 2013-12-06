# Makefile for binary lens code
###############################


ifeq (para,$(firstword $(MAKECMDGOALS)))
	F77 = mpif77
	FFLAGS = -xf77-cpp-input -DMPI_ON
else
	F77 = g77
	FFLAGS = -xf77-cpp-input
endif


para:
	$(F77) $(FFLAGS) parallel_dev.f -o para_lens



serial:
	$(F77) parallel_dev.f -o serial_lens
