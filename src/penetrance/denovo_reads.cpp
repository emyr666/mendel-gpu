#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif


void DenovoReadsMendelGPU::compute_penetrance(){
  double start = clock();
  read_penetrance->populate_read_matrices();
  read_penetrance->process_read_matrices();
  read_penetrance->prefetch_reads(g_left_marker+g_markers,1);
  cerr<<"Compute denovo reads penentrance in "<<(clock()-start)/CLOCKS_PER_SEC<<" seconds.\n";
}


