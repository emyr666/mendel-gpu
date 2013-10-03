#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"

void DenovoGlfMendelGPU::allocate_memory(){
  DenovoMendelGPU::allocate_memory();
  cerr<<"Initializing variables for denovo GLF haplotyper\n";
}

DenovoGlfMendelGPU::DenovoGlfMendelGPU(IO_manager * io):DenovoMendelGPU(io){
}

DenovoGlfMendelGPU::~DenovoGlfMendelGPU(){
  cerr<<"Entering destructor denovo GLF haplotyper\n";
  cerr<<"Exiting destructor denovo GLF haplotyper\n";

}
