# BELOW ARE EXAMPLES OF POTENTIAL HARDWARE PROFILES IF YOU PLAN TO DEPLOY
# THE PROGRAM ON MULTIPLE MACHINES
#include locations.baseline
include locations.gpu

FORTRAN_OBJS=fortran_wrapper.o
CL_OBJS=cl_interface.o clsafe.o mendel_gpu_common.o mendel_gpu_denovo.o mendel_gpu_guide.o

OBJS=$(CL_OBJS) $(FORTRAN_OBJS)

C_COMPILER=g++
F_COMPILER=gfortran

PROGRAM=../bin/impute

impute: $(OBJS)
	$(LINKER) -o $(PROGRAM) $(LINK_FLAGS) $(OBJS)

utility: compress_float.cpp prepare_inputfiles.cpp chunk.cpp prepare_denovo.cpp make_haplo_posteriors.cpp split_file.cpp
	g++ -o ../bin/compress_float compress_float.cpp
	g++ -o ../bin/prepare_inputfiles prepare_inputfiles.cpp
	g++ -o ../bin/chunk chunk.cpp
	g++ -o ../bin/prepare_denovo prepare_denovo.cpp
	g++ -o ../bin/make_haplo_posteriors make_haplo_posteriors.cpp
	g++ -o ../bin/split_file split_file.cpp

fortran_wrapper.o: fortran_wrapper.f90
	$(F_COMPILER) -c $<

cl_interface.o: cl_interface.cpp mendel_gpu.hpp clsafe.h cl_constants.h
	$(C_COMPILER) $(INC_FLAGS) $(PREPROC) -c $<

mendel_gpu_common.o: mendel_gpu_common.cpp clsafe.h cl_constants.h
	$(C_COMPILER) $(INC_FLAGS) $(PREPROC) -c $<

mendel_gpu_denovo.o: mendel_gpu_denovo.cpp clsafe.h cl_constants.h
	$(C_COMPILER) $(INC_FLAGS) $(PREPROC) -c $<

mendel_gpu_guide.o: mendel_gpu_guide.cpp clsafe.h cl_constants.h
	$(C_COMPILER) $(INC_FLAGS) $(PREPROC) -c $<

clsafe.o: clsafe.c clsafe.h
	$(C_COMPILER) $(INC_FLAGS) $(PREPROC) -c $<

clean:
	rm -f $(PROGRAM) $(OBJS) *.mod
