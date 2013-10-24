using namespace std;
#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"

GuidedReadsMendelGPU::GuidedReadsMendelGPU(IO_manager * io):GuidedMendelGPU(io){
}

void GuidedReadsMendelGPU::allocate_memory(){
  GuidedMendelGPU::allocate_memory();
  read_penetrance = new ReadPenetrance(this);
  read_penetrance->prefetch_reads(0,2*g_flanking_snps);
}

void GuidedReadsMendelGPU::load_datasets(){
  cerr<<"Loading dataset for guided reads haplotyper\n";
  io_manager->read_input(ref_haplotype,this->informative_snp,this->g_people,this->g_snps, this->ref_haplotypes);
}

GuidedReadsMendelGPU::~GuidedReadsMendelGPU(){
  cerr<<"Entering destructor guided reads haplotyper\n";
  delete read_penetrance;
  cerr<<"Exiting destructor guided reads haplotyper\n";
}

void GuidedReadsMendelGPU::init_window(){
  GuidedMendelGPU::init_window();
  // Load the reads for this current window
  if (run_gpu){
    init_window_opencl();
  }
  return;
}


