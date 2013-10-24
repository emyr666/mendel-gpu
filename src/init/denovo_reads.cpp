#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"

DenovoReadsMendelGPU::DenovoReadsMendelGPU(IO_manager * io):DenovoMendelGPU(io){
}

void DenovoReadsMendelGPU::allocate_memory(){
  DenovoMendelGPU::allocate_memory();
  cerr<<"Initializing variables for denovo reads haplotyper\n";
  cerr<<"Analysis on "<<g_snps<<" SNPs and "<<g_people<<" persons\n";
  read_penetrance = new ReadPenetrance(this);
  read_penetrance->prefetch_reads(0,1);
  cerr<<"Initialized variables for denovo reads haplotyper\n";
}

void DenovoReadsMendelGPU::load_datasets(){
  io_manager->read_input(this->g_people,this->g_snps);
  cerr<<"Loaded input in DenovoReadsMendelGPU\n";
}


DenovoReadsMendelGPU::~DenovoReadsMendelGPU(){
  cerr<<"Entering destructor denovo reads haplotyper\n";
  delete read_penetrance;
  cerr<<"Exiting destructor denovo reads haplotyper\n";

}


void DenovoReadsMendelGPU::init_window(){
  DenovoMendelGPU::init_window();
  // Load the reads for this current window
  if (run_gpu){
    init_window_opencl();
  }
  return;
}
