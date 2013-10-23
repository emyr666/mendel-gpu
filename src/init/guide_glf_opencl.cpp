#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void GuidedGlfMendelGPU::init_opencl(){
  GuidedMendelGPU::init_opencl();
  if(run_gpu){
#ifdef USE_GPU
    createKernel("precompute_penetrance",kernel_precompute_penetrance);
    cerr<<"GuidedGlfMendelGPU kernels created\n";
    createBuffer<float>(CL_MEM_READ_ONLY,g_people*geno_dim*g_snps,"buffer_snp_penetrance",buffer_snp_penetrance);
    writeToBuffer(buffer_snp_penetrance,g_people*geno_dim*g_snps,g_snp_penetrance,"buffer_snp_penetrance");
    cerr<<"GuidedGlfMendelGPU buffers created\n";
    int arg;
    arg = 0;
    setArg(kernel_precompute_penetrance,arg, g_max_haplotypes,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, g_max_window,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, penetrance_matrix_size,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, g_snps,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, gf_logpen_threshold,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, packedhap_len,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, *buffer_markers,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, *buffer_haplotypes,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, *buffer_left_marker,"kernel_precompute_penetrance");
    if(geno_dim==UNPHASED_INPUT){
      setArg(kernel_precompute_penetrance,arg, *buffer_haploid_arr,"kernel_precompute_penetrance");
    }
    setArg(kernel_precompute_penetrance,arg, *buffer_packedhap,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, *buffer_haplotype,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, *buffer_snp_penetrance,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, *buffer_logpenetrance_cache,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, *buffer_penetrance_cache,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, *buffer_extended_snp_mapping,"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, cl::__local(sizeof(int)*g_max_window),"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, cl::__local(sizeof(packedhap_t)*g_max_haplotypes*packedhap_len),"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, cl::__local(sizeof(int)*g_max_window),"kernel_precompute_penetrance");
    if(geno_dim==UNPHASED_INPUT){
      setArg(kernel_precompute_penetrance,arg, cl::__local(sizeof(int)*g_max_window),"kernel_precompute_penetrance");
    }
    if(geno_dim==PHASED_INPUT){
      setArg(kernel_precompute_penetrance,arg, cl::__local(sizeof(float)*2*BLOCK_WIDTH),"kernel_precompute_penetrance");
    }
    setArg(kernel_precompute_penetrance,arg, cl::__local(sizeof(float)*g_max_window*geno_dim),"kernel_precompute_penetrance");
    setArg(kernel_precompute_penetrance,arg, cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_precompute_penetrance");
    cerr<<"Total args for precompute pen: "<<arg<<endl;
    int kernelWorkGroupSize = kernel_precompute_penetrance->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel precompute_penetrance");
    cerr<<"precompute_penetrance kernel work group size is "<<kernelWorkGroupSize<<endl;
#endif
  }
}

void GuidedGlfMendelGPU::free_opencl(){
  if(run_gpu){
#ifdef USE_GPU
    delete buffer_snp_penetrance;
    delete kernel_precompute_penetrance;
    GuidedMendelGPU::free_opencl();
#endif
  }
}
