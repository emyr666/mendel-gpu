#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"


void GuidedReadsMendelGPU::compute_penetrance(){
  cerr<<"Entering precompute_penetrance at left marker "<<g_left_marker<<"\n";
  read_penetrance->populate_read_matrices();
  read_penetrance->process_read_matrices();
  return;
}
  
  

