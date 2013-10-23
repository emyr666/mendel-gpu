#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void GuidedReadsMendelGPU::init_opencl(){
  cerr<<"Initializing OpenCL for GuidedReadsMendelGPU\n";
  if(run_gpu){
#ifdef USE_GPU
    GuidedMendelGPU::init_opencl();
    read_penetrance->init_opencl();
  #endif
  }
}

void GuidedReadsMendelGPU::free_opencl(){
  if(run_gpu){
#ifdef USE_GPU
    read_penetrance->free_opencl();
    GuidedMendelGPU::free_opencl();
#endif
  }
}
