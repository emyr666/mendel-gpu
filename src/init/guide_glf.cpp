using namespace std;
#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"

GuidedGlfMendelGPU::GuidedGlfMendelGPU(IO_manager * io):GuidedMendelGPU(io){
}

void GuidedGlfMendelGPU::load_datasets(){
  cerr<<"Loading dataset for guided GLF haplotyper\n";
  io_manager->read_input(ref_haplotype,g_snp_penetrance,this->informative_snp,this->g_people,this->g_snps, this->ref_haplotypes);
}

void GuidedGlfMendelGPU::allocate_memory(){
  GuidedMendelGPU::allocate_memory();
}

GuidedGlfMendelGPU::~GuidedGlfMendelGPU(){
  cerr<<"Entering destructor guided GLF haplotyper\n";
  cerr<<"Exiting destructor guided GLF haplotyper\n";
}
