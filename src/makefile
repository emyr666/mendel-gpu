# BELOW ARE EXAMPLES OF POTENTIAL HARDWARE PROFILES IF YOU PLAN TO DEPLOY
# THE PROGRAM ON MULTIPLE MACHINES

#TARGET=baseline
TARGET=epigraph
#TARGET=soto
#TARGET=mec
ifeq ($(TARGET),baseline)
  #CL_ROOT=/opt/firestream
  # 1st arg make sure you see directory CL under include
  # 2nd arg make sure you see directory GL under include 
  #INC_FLAGS=-I$(CL_ROOT)/include -I/opt/X11R6/include
  # 1st arg make sure libOpenCL.so is under lib
  #LINK_FLAGS=-L$(CL_ROOT)/lib -lgfortran -lOpenCL
  LINK_FLAGS=-lgfortran -lstdc++
  LINKER=gfortran
else ifeq ($(TARGET),mec)
  PREPROC=-DUSE_GPU
  # this should point to where your OpenCL installation is (below would be things like lib)
  CL_ROOT=/opt/firestream
  # 1st arg make sure you see directory CL under include
  # 2nd arg make sure you see directory GL under include 
  INC_FLAGS=-I$(CL_ROOT)/include -I/opt/X11R6/include
  # 1st arg make sure libOpenCL.so is under lib
  LINK_FLAGS=-L$(CL_ROOT)/lib -lgfortran -lOpenCL
  LINKER=g++
else ifeq ($(TARGET),soto)
  PREPROC=-DUSE_GPU
  # this should point to where your OpenCL installation is (below would be things like lib)
  CL_ROOT=/opt/AMDAPP
  # 1st arg make sure you see directory CL under include
  # 2nd arg make sure you see directory GL under include 
  INC_FLAGS=-I$(CL_ROOT)/include 
  # 1st arg make sure libOpenCL.so is under lib
  LINK_FLAGS=-L$(CL_ROOT)/lib -lgfortran -lOpenCL
  LINKER=g++
else ifeq ($(TARGET),epigraph)
  PREPROC=-DUSE_GPU
  # this should point to where your OpenCL installation is (below would be things like lib)
  CL_ROOT=/home/garykche/ati-stream-sdk-v2.2-lnx64
  # 1st arg make sure you see directory CL under include
  # 2nd arg make sure you see directory GL under include 
  INC_FLAGS=-I$(CL_ROOT)/include -I/usr/local/cuda_sdk/shared/inc
  # 1st arg make sure libOpenCL.so is under lib
  LINK_FLAGS=-L$(CL_ROOT)/lib/x86_64 -lgfortran -lOpenCL
  LINKER=gfortran
else
  LINK_FLAGS=-lgfortran -lstdc++
  LINKER=gfortran
endif

FORTRAN_OBJS=fortran_wrapper.o
CL_OBJS=cl_interface.o clsafe.o mendel_gpu_common.o mendel_gpu_denovo.o mendel_gpu_guide.o

OBJS=$(CL_OBJS) $(FORTRAN_OBJS)

C_COMPILER=g++
F_COMPILER=gfortran

PROGRAM=../bin/impute

impute: $(OBJS)
	$(LINKER) -o $(PROGRAM) $(LINK_FLAGS) $(OBJS)

utility: compress_float.cpp prepare_inputfiles.cpp chunk.cpp prepare_denovo.cpp make_haplo_posteriors.cpp
	g++ -o ../bin/compress_float compress_float.cpp
	g++ -o ../bin/prepare_inputfiles prepare_inputfiles.cpp
	g++ -o ../bin/chunk chunk.cpp
	g++ -o ../bin/prepare_denovo prepare_denovo.cpp
	g++ -o ../bin/make_haplo_posteriors make_haplo_posteriors.cpp

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
