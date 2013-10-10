#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif



void DenovoMendelGPU::init_opencl(){
  cerr<<"Initializing OpenCL for DenovoMendelGPU\n";
  MendelGPU::init_opencl();
  if(run_gpu){
#ifdef USE_GPU
    createKernel("impute_genotype_denovo",kernel_impute_genotype_denovo);
    cerr<<"De novo kernels created\n";
    createBuffer<int>(CL_MEM_READ_WRITE,g_people,"buffer_subject_genotype",buffer_subject_genotype);
    createBuffer<float>(CL_MEM_READ_WRITE,g_people,"buffer_subject_dosage",buffer_subject_dosage);
    createBuffer<float>(CL_MEM_READ_WRITE,4*g_people,"buffer_subject_posterior_prob",buffer_subject_posterior_prob);
    createBuffer<int>(CL_MEM_READ_WRITE,g_max_haplotypes,"buffer_center_dosage",buffer_center_dosage);
    createBuffer<int>(CL_MEM_READ_ONLY,g_max_haplotypes,"buffer_twin_hap_index",buffer_twin_hap_index);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_prev_left_marker",buffer_prev_left_marker);
    createBuffer<int>(CL_MEM_READ_ONLY,g_max_haplotypes,"buffer_beyond_left_edge_dosage",buffer_beyond_left_edge_dosage);
    createBuffer<int>(CL_MEM_READ_ONLY,g_max_haplotypes,"buffer_right_edge_dosage",buffer_right_edge_dosage);
    cerr<<"De novo buffers created\n";
    int arg;
    arg = 0;
    setArg(kernel_impute_genotype_denovo,arg,g_max_haplotypes,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,g_max_window,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,gf_epsilon,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,penetrance_matrix_size,"kernel_impute_genotype_denovo");
    int dim = g_genotype_imputation?3:4;
    setArg(kernel_impute_genotype_denovo,arg,dim,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_genotype_imputation,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_markers,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_haplotypes,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_left_marker,"kernel_impute_genotype_denovo");
    if(geno_dim==UNPHASED_INPUT){
      setArg(kernel_impute_genotype_denovo,arg,*buffer_haploid_arr,"kernel_impute_genotype_denovo");
    }
    setArg(kernel_impute_genotype_denovo,arg,*buffer_frequency,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_penetrance_cache,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_subject_genotype,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_subject_dosage,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_subject_posterior_prob,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_active_haplotype,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_center_dosage,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,*buffer_subject_haplotype_weight,"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(float)*g_max_haplotypes),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(float)*4),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(float)*g_max_haplotypes),"kernel_impute_genotype_denovo");
    setArg(kernel_impute_genotype_denovo,arg,cl::__local(sizeof(float)*1),"kernel_impute_genotype_denovo");
    cerr<<"De novo kernel arguments set\n";
  #endif
  }
}

void DenovoMendelGPU::init_window_opencl(){
  MendelGPU::init_window_opencl();
  if (run_gpu){
    #ifdef USE_GPU
    writeToBuffer(buffer_prev_left_marker, 1, &g_prev_left_marker, "buffer_prev_left_marker" );
    #endif
    cerr<<"Buffers sent to GPU for current window\n";
  }
}
