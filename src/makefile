# BELOW ARE EXAMPLES OF POTENTIAL HARDWARE PROFILES IF YOU PLAN TO DEPLOY
# THE PROGRAM ON MULTIPLE MACHINES

include include.mk

DIRS=init
DIRS+=penetrance
DIRS+=imputation

LIBDIR=../lib

FORTRAN_OBJS=matrix_completion.o fortran_wrapper.o 
COMMON_OBJS=main.o io_manager.o cl_templates.o
SUB_OBJS=init/*.o penetrance/*.o imputation/*.o

OBJS=$(C_OBJS)
#OBJS=$(C_OBJS) $(FORTRAN_OBJS)

C_COMPILER=g++
F_COMPILER=gfortran

PROGRAM=../bin/impute

impute::
	for dir in $(DIRS); do(cd $$dir; make $1 || exit 1) || exit 1; done

impute:: $(COMMON_OBJS)
	$(LINKER) -o $(PROGRAM) $(LINK_FLAGS) $(SUB_OBJS) $(COMMON_OBJS)  -L. -ltools $(LIBDIR)/libbam.a

matrix_completion.o: matrix_completion.f90
	$(F_COMPILER) $(FINC_FLAGS) -c $<

fortran_wrapper.o: fortran_wrapper.f90
	$(F_COMPILER) $(FINC_FLAGS) -c $<

main.o: main.cpp mendel_gpu.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

io_manager.o: io_manager.cpp io_manager.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

bam_dummyinterface.o: bam_dummyinterface.c bam_interface.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

cl_templates.o: cl_templates.cpp cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

clean:
	for dir in $(DIRS); do(cd $$dir; rm -f *.o  || exit 1) || exit 1; done
	rm -f $(PROGRAM) $(OBJS) *.mod *.o
