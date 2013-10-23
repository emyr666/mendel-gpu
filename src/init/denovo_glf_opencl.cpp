#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void DenovoGlfMendelGPU::init_opencl(){
  DenovoMendelGPU::init_opencl();
  if(run_gpu){
#ifdef USE_GPU
    createKernel("impute_penetrance",kernel_impute_penetrance);
    createKernel("precompute_penetrance_fast",kernel_precompute_penetrance_fast);
    cerr<<"DenovoGlfMendelGPU kernels created\n";


    createBuffer<float>(CL_MEM_READ_ONLY,g_people*geno_dim*g_snps,"buffer_snp_penetrance",buffer_snp_penetrance);
    writeToBuffer(buffer_snp_penetrance,g_people*geno_dim*g_snps,g_snp_penetrance,"buffer_snp_penetrance");
    cerr<<"DenovoGlfMendelGPU buffers created\n";
    int arg;
    arg = 0;
    setArg(kernel_impute_penetrance,arg,g_max_haplotypes,"kernel_impute_penetrance");
    setArg(kernel_impute_penetrance,arg,penetrance_matrix_size,"kernel_impute_penetrance");
    setArg(kernel_impute_penetrance,arg,*buffer_active_haplotype,"kernel_impute_penetrance");
    setArg(kernel_impute_penetrance,arg,*buffer_twin_hap_index,"kernel_impute_penetrance");
    setArg(kernel_impute_penetrance,arg,*buffer_logpenetrance_cache,"kernel_impute_penetrance");
    setArg(kernel_impute_penetrance,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_impute_penetrance");
    if(geno_dim==UNPHASED_INPUT){
      setArg(kernel_impute_penetrance,arg,cl::__local(sizeof(float)*g_max_haplotypes),"kernel_impute_penetrance");
    }
    setArg(kernel_impute_penetrance,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_impute_penetrance");

    arg = 0;
    setArg(kernel_precompute_penetrance_fast,arg,g_max_haplotypes,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,penetrance_matrix_size,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,g_snps,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,gf_logpen_threshold,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_haplotypes,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_markers,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_prev_left_marker,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_left_marker,"kernel_precompute_penetrance_fast");
    if(geno_dim==UNPHASED_INPUT){
      setArg(kernel_precompute_penetrance_fast,arg,*buffer_haploid_arr,"kernel_precompute_penetrance_fast");
    }
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_beyond_left_edge_dosage,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_right_edge_dosage,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_snp_penetrance,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_logpenetrance_cache,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_penetrance_cache,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,*buffer_active_haplotype,"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,cl::__local(sizeof(float)*geno_dim),"kernel_precompute_penetrance_fast");
    setArg(kernel_precompute_penetrance_fast,arg,cl::__local(sizeof(float)*geno_dim),"kernel_precompute_penetrance_fast");
    if(geno_dim==PHASED_INPUT){
      setArg(kernel_precompute_penetrance_fast,arg,cl::__local(sizeof(float)*2*BLOCK_WIDTH),"kernel_precompute_penetrance_fast");
    }
    setArg(kernel_precompute_penetrance_fast,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_precompute_penetrance_fast");

  #endif
  }
}

void DenovoGlfMendelGPU::free_opencl(){
  if(run_gpu){
#ifdef USE_GPU
    delete buffer_snp_penetrance;
    delete kernel_precompute_penetrance_fast;
    delete kernel_impute_penetrance;
    DenovoMendelGPU::free_opencl();
#endif
  }
}
