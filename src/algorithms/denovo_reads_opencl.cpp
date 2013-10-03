#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../clsafe.h"
#endif

void DenovoReadsMendelGPU::init_opencl(){
  DenovoMendelGPU::init_opencl();
  if(run_gpu){
#ifdef USE_GPU
    clSafe(err,"creating kernel impute_genotype_denovo");
    kernel_impute_genotype_denovo_haploid = new cl::Kernel(*program,"impute_genotype_denovo_haploid",&err);
    clSafe(err,"creating kernel impute_genotype_denovo_haplod");
    buffer_beyond_left_edge_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    clSafe(err,"creating buffer beyond left edge dosage");
    buffer_center_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    clSafe(err,"creating buffer center dosage");
    buffer_right_edge_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    clSafe(err,"creating buffer right edge dosage");
    buffer_twin_hap_index = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    clSafe(err,"creating buffer twin index");
    buffer_subject_posterior_prob = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*4, NULL, &err);
    clSafe(err,"creating buffer subject_posterior_prob");
    buffer_subject_genotype = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_people, NULL, &err);
    clSafe(err,"creating buffer subject genotype");
    buffer_subject_dosage = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people, NULL, &err);
    clSafe(err,"creating buffer subject dosage");
    buffer_prev_left_marker = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer prev left marker");
    buffer_center_snp = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer center marker");
  #endif
  }
}
