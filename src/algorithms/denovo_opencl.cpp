#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../clsafe.h"
#endif


void DenovoMendelGPU::init_opencl(){
  MendelGPU::init_opencl();
  if(run_gpu){
#ifdef USE_GPU
    //kernel_precompute_penetrance_fast = new cl::Kernel(*program,"precompute_penetrance_fast",&err);
    //clSafe(err,"creating kernel precompute_fast penetrance");
    //kernel_precompute_penetrance_fast_haploid = new cl::Kernel(*program,"precompute_penetrance_fast_haploid",&err);
    //clSafe(err,"creating kernel precompute_fast penetrance haploid");
    //kernel_impute_penetrance = new cl::Kernel(*program,"impute_penetrance",&err);
    //clSafe(err,"creating kernel impute penetrance");
    //kernel_impute_penetrance_haploid = new cl::Kernel(*program,"impute_penetrance_haploid",&err);
    //clSafe(err,"creating kernel impute penetrance_haploid");
    kernel_impute_genotype_denovo = new cl::Kernel(*program,"impute_genotype_denovo",&err);
    clSafe(err,"creating kernel impute_genotype_denovo");
    kernel_impute_genotype_denovo_haploid = new cl::Kernel(*program,"impute_genotype_denovo_haploid",&err);
    clSafe(err,"creating kernel impute_genotype_denovo_haplod");
    //buffer_beyond_left_edge_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    //clSafe(err,"creating buffer beyond left edge dosage");
    //buffer_center_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    //clSafe(err,"creating buffer center dosage");
    //buffer_right_edge_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    //clSafe(err,"creating buffer right edge dosage");
    //buffer_twin_hap_index = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    //clSafe(err,"creating buffer twin index");
    buffer_subject_posterior_prob = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*4, NULL, &err);
    clSafe(err,"creating buffer subject_posterior_prob");
    buffer_subject_genotype = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_people, NULL, &err);
    clSafe(err,"creating buffer subject genotype");
    buffer_subject_dosage = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people, NULL, &err);
    clSafe(err,"creating buffer subject dosage");
    //buffer_prev_left_marker = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating buffer prev left marker");
    //buffer_center_snp = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating buffer center marker");
  #endif
  }
}

void DenovoMendelGPU::init_window_opencl(){
  if (run_gpu){
    #ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_active_haplotype, CL_TRUE, 0,  sizeof(int) * g_max_haplotypes, g_active_haplotype, NULL, NULL );
    clSafe(err, "write active haplotype indices");
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotypes, CL_TRUE, 0,  sizeof(int), &g_haplotypes, NULL, NULL );
    clSafe(err, "write haplotypes");
    err = commandQueue->enqueueWriteBuffer(*buffer_markers, CL_TRUE, 0,  sizeof(int), &g_markers, NULL, NULL );
    clSafe(err, "write markers");
    err = commandQueue->enqueueWriteBuffer(*buffer_left_marker, CL_TRUE, 0,  sizeof(int), &g_left_marker, NULL, NULL );
    clSafe(err, "write left marker");
    err = commandQueue->enqueueWriteBuffer(*buffer_prev_left_marker, CL_TRUE, 0,  sizeof(int), &g_prev_left_marker, NULL, NULL );
    clSafe(err, "write prev left marker");
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotype, CL_TRUE, sizeof(int)*0,  sizeof(int)*g_max_window*g_max_haplotypes, g_haplotype, NULL, NULL );
    clSafe(err, "write haplotype");
    #endif
    cerr<<"Buffers sent to GPU\n";
  }
}
