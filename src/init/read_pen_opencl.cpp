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
    //mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people*read_compact_matrix_size,"buffer_read_alleles_mat",buffer_read_alleles_mat);
    //mendelgpu->createBuffer<float>(CL_MEM_READ_ONLY,mendelgpu->g_people*read_compact_matrix_size,"buffer_read_match_logmat",buffer_read_match_logmat);
    //mendelgpu->createBuffer<float>(CL_MEM_READ_ONLY,mendelgpu->g_people*read_compact_matrix_size,"buffer_read_mismatch_logmat",buffer_read_mismatch_logmat);
    mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people*max_vector_length,"buffer_read_alleles_vec",buffer_read_alleles_vec);
    mendelgpu->createBuffer<float>(CL_MEM_READ_ONLY,mendelgpu->g_people*max_vector_length,"buffer_read_match_logvec",buffer_read_match_logvec);
    mendelgpu->createBuffer<float>(CL_MEM_READ_ONLY,mendelgpu->g_people*max_vector_length,"buffer_read_mismatch_logvec",buffer_read_mismatch_logvec);

    mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people*compact_rows,"buffer_read_lengths",buffer_read_lengths);
    mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people*compact_rows,"buffer_vector_offsets",buffer_vector_offsets);
    mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people*compact_rows,"buffer_haplotype_offsets",buffer_haplotype_offsets);

    mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people,"buffer_mat_rows_by_subject",buffer_mat_rows_by_subject);
    mendelgpu->createBuffer<int>(CL_MEM_READ_ONLY,mendelgpu->g_people*mendelgpu->g_max_haplotypes,"buffer_best_haplotypes",buffer_best_haplotypes);


    int arg;
    arg = 0;
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,mendelgpu->g_max_window,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,mendelgpu->g_max_haplotypes,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,mendelgpu->penetrance_matrix_size,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,read_compact_matrix_size,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,max_vector_length,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,compact_rows,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*mendelgpu->buffer_markers,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_mat_rows_by_subject,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_read_alleles_vec,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_read_match_logvec,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_read_mismatch_logvec,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_vector_offsets,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_read_lengths,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_haplotype_offsets,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*mendelgpu->buffer_haplotype,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*mendelgpu->buffer_logpenetrance_cache,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*mendelgpu->buffer_active_haplotype,"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,*buffer_best_haplotypes,"kernel_reads_compute_penetrance");
    // alleles
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(int)*max_vector_length),"kernel_reads_compute_penetrance");
    // matches
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*max_vector_length),"kernel_reads_compute_penetrance");
    // mismatches
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*max_vector_length),"kernel_reads_compute_penetrance");
    // vector_offsets
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(int)*compact_rows),"kernel_reads_compute_penetrance");
    // read lengths
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(int)*compact_rows),"kernel_reads_compute_penetrance");
    // hap offsets
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(int)*compact_rows),"kernel_reads_compute_penetrance");

    // haplotype 1
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(int)*mendelgpu->g_max_window),"kernel_reads_compute_penetrance");
    // haplotype 2 
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(int)*mendelgpu->g_max_window),"kernel_reads_compute_penetrance");
    // read alleles
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(int)*mendelgpu->g_max_window),"kernel_reads_compute_penetrance");
    // match
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*mendelgpu->g_max_window),"kernel_reads_compute_penetrance");
    // mismatch
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*mendelgpu->g_max_window),"kernel_reads_compute_penetrance");

    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_compute_penetrance");
    mendelgpu->setArg(kernel_reads_compute_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_compute_penetrance");

    arg = 0;
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,mendelgpu->g_max_haplotypes,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,mendelgpu->penetrance_matrix_size,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,mendelgpu->gf_logpen_threshold,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,*mendelgpu->buffer_logpenetrance_cache,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,*mendelgpu->buffer_penetrance_cache,"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,*mendelgpu->buffer_active_haplotype,"kernel_reads_adjust_penetrance");
    //mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,*buffer_best_haplotypes,"kernel_reads_adjust_penetrance");
    //mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,cl::__local(sizeof(int)*mendelgpu->g_max_haplotypes),"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,cl::__local(sizeof(int)*mendelgpu->g_max_haplotypes),"kernel_reads_adjust_penetrance");
    mendelgpu->setArg(kernel_reads_adjust_penetrance,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_reads_adjust_penetrance");
  #endif
}

void ReadPenetrance::free_opencl(){
#ifdef USE_GPU
    delete buffer_read_alleles_vec;
    delete buffer_read_match_logvec;
    delete buffer_read_mismatch_logvec;
    delete buffer_vector_offsets;
    delete buffer_read_lengths;
    delete buffer_haplotype_offsets;
    delete buffer_mat_rows_by_subject;
#endif
}
