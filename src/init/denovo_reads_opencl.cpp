#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void DenovoReadsMendelGPU::init_opencl(){
  cerr<<"Initializing OpenCL for DenovoReadsMendelGPU\n";
  if(run_gpu){
#ifdef USE_GPU
    DenovoMendelGPU::init_opencl();
    createKernel("reads_compute_penetrance",kernel_reads_compute_penetrance);
    createKernel("reads_adjust_penetrance",kernel_reads_adjust_penetrance);
    cerr<<"Read compact matrix size people are "<<read_compact_matrix_size<<","<<g_people<<endl;
    createBuffer<int>(CL_MEM_READ_ONLY,g_people*read_compact_matrix_size,"buffer_read_alleles_mat",buffer_read_alleles_mat);
    createBuffer<float>(CL_MEM_READ_ONLY,g_people*read_compact_matrix_size,"buffer_read_match_logmat",buffer_read_match_logmat);
    createBuffer<float>(CL_MEM_READ_ONLY,g_people*read_compact_matrix_size,"buffer_read_mismatch_logmat",buffer_read_mismatch_logmat);
    createBuffer<int>(CL_MEM_READ_ONLY,g_people,"buffer_mat_rows_by_subject",buffer_mat_rows_by_subject);
    int arg;
    arg = 0;
    setArg(kernel_reads_compute_penetrance,arg,g_max_window,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,g_max_haplotypes,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,penetrance_matrix_size,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,read_compact_matrix_size,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,*buffer_markers,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,*buffer_mat_rows_by_subject,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,*buffer_read_alleles_mat,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,*buffer_read_match_logmat,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,*buffer_read_mismatch_logmat,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,*buffer_haplotype,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,*buffer_logpenetrance_cache,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,*buffer_active_haplotype,"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_compute_penetrance");
    setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_compute_penetrance");
    arg = 0;
    setArg(kernel_reads_adjust_penetrance,arg,g_max_haplotypes,"kernel_reads_adjust_penetrance");
    setArg(kernel_reads_adjust_penetrance,arg,penetrance_matrix_size,"kernel_reads_adjust_penetrance");
    setArg(kernel_reads_adjust_penetrance,arg,gf_logpen_threshold,"kernel_reads_adjust_penetrance");
    setArg(kernel_reads_adjust_penetrance,arg,*buffer_logpenetrance_cache,"kernel_reads_adjust_penetrance");
    setArg(kernel_reads_adjust_penetrance,arg,*buffer_penetrance_cache,"kernel_reads_adjust_penetrance");
    setArg(kernel_reads_adjust_penetrance,arg,*buffer_active_haplotype,"kernel_reads_adjust_penetrance");
    setArg(kernel_reads_adjust_penetrance,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_reads_adjust_penetrance");
    setArg(kernel_reads_adjust_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_adjust_penetrance");
  #endif
  }
}

void DenovoReadsMendelGPU::free_opencl(){
  if(run_gpu){
#ifdef USE_GPU
    DenovoMendelGPU::init_opencl();
#endif
  }
}
