#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void ReadPenetrance::init_opencl(){
  cerr<<"Initializing OpenCL for ReadPenetrance\n";
#ifdef USE_GPU
    mendelgpu->createKernel("reads_compute_penetrance",kernel_reads_compute_penetrance);
    mendelgpu->createKernel("reads_adjust_penetrance",kernel_reads_adjust_penetrance);
    cerr<<"Read compact matrix size people are "<<read_compact_matrix_size<<","<<mendelgpu->g_people<<endl;
    mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people*read_compact_matrix_size,"buffer_read_alleles_mat",buffer_read_alleles_mat);
    mendelgpu->createBuffer<float>(CL_MEM_READ_ONLY,mendelgpu->g_people*read_compact_matrix_size,"buffer_read_match_logmat",buffer_read_match_logmat);
    mendelgpu->createBuffer<float>(CL_MEM_READ_ONLY,mendelgpu->g_people*read_compact_matrix_size,"buffer_read_mismatch_logmat",buffer_read_mismatch_logmat);
    mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people,"buffer_mat_rows_by_subject",buffer_mat_rows_by_subject);
    int arg;
    arg = 0;
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,mendelgpu->g_max_window,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,mendelgpu->g_max_haplotypes,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,mendelgpu->penetrance_matrix_size,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,read_compact_matrix_size,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*mendelgpu->buffer_markers,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_mat_rows_by_subject,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_read_alleles_mat,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_read_match_logmat,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_read_mismatch_logmat,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*mendelgpu->buffer_haplotype,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*mendelgpu->buffer_logpenetrance_cache,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*mendelgpu->buffer_active_haplotype,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_compute_penetrance");

    arg = 0;
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,mendelgpu->g_max_haplotypes,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,mendelgpu->penetrance_matrix_size,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,mendelgpu->gf_logpen_threshold,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,*mendelgpu->buffer_logpenetrance_cache,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,*mendelgpu->buffer_penetrance_cache,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,*mendelgpu->buffer_active_haplotype,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,cl::__local(sizeof(int)*mendelgpu->g_max_haplotypes),"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_adjust_penetrance");
  #endif
}

void ReadPenetrance::free_opencl(){
#ifdef USE_GPU
    delete buffer_read_alleles_mat;
    delete buffer_read_match_logmat;
    delete buffer_read_mismatch_logmat;
    delete buffer_mat_rows_by_subject;
#endif
}
