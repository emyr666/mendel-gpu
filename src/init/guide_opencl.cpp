#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void GuidedMendelGPU::init_opencl(){
  MendelGPU::init_opencl();
  if(run_gpu){
#ifdef USE_GPU
    // CREATE ALL KERNELS HERE
    createKernel("precompute_penetrance",kernel_precompute_penetrance);
    createKernel("impute_genotype_guide",kernel_impute_genotype_guide);
    // CREATE ALL BUFFERS HERE
    createBuffer<packedhap_t>(CL_MEM_READ_ONLY,g_max_haplotypes*packedhap_len,"buffer_packedhap",buffer_packedhap);
    createBuffer<float>(CL_MEM_READ_ONLY,g_people*geno_dim*g_snps,"buffer_snp_penetrance",buffer_snp_penetrance);
    writeToBuffer(buffer_snp_penetrance,g_people*geno_dim*g_snps,g_snp_penetrance,"buffer_snp_penetrance");
    createBuffer<int>(CL_MEM_READ_ONLY,g_max_window,"buffer_extended_snp_mapping",buffer_extended_snp_mapping);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_extended_haplotypes",buffer_extended_haplotypes);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_center_snp_end",buffer_center_snp_end);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_center_snp_end",buffer_center_snp_end);
    createBuffer<int>(CL_MEM_READ_ONLY,ref_haplotypes,"buffer_extended_root_mapping",buffer_extended_root_mapping);
    createBuffer<float>(CL_MEM_READ_ONLY,ref_haplotypes,"buffer_extended_frequency",buffer_extended_frequency);
    createBuffer<int>(CL_MEM_READ_WRITE,g_people*g_flanking_snps,"buffer_subject_genotype_block",buffer_subject_genotype_block);
    createBuffer<float>(CL_MEM_READ_WRITE,g_people*g_flanking_snps,"buffer_subject_dosage_block",buffer_subject_dosage_block);
    createBuffer<float>(CL_MEM_READ_WRITE,g_people*g_flanking_snps*4,"buffer_subject_posterior_prob_block",buffer_subject_posterior_prob_block);
    createBuffer<packedhap_t>(CL_MEM_READ_ONLY,ref_haplotypes*packedextendedhap_len,"buffer_packedextendedhap",buffer_packedextendedhap);
    cerr<<"GPU Buffers created\n";
    // SET KERNEL ARGUMENTS HERE
    int arg;
      // PRECOMPUTE PENETRANCE FOR THE GUIDE HAP APPROACH
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
      arg = 0;
      setArg(kernel_impute_genotype_guide,arg, g_max_haplotypes,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, gf_epsilon,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, penetrance_matrix_size,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, packedextendedhap_len,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, g_flanking_snps,"kernel_impute_genotype_guide");
      int max_geno = g_genotype_imputation?3:4;
      setArg(kernel_impute_genotype_guide,arg, max_geno,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_genotype_imputation,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_haplotypes,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_extended_haplotypes,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_center_snp_end,"kernel_impute_genotype_guide");
      if(geno_dim==UNPHASED_INPUT){
        setArg(kernel_impute_genotype_guide,arg, *buffer_haploid_arr,"kernel_impute_genotype_guide");
      }
      setArg(kernel_impute_genotype_guide,arg, *buffer_extended_root_mapping,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_extended_frequency,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_penetrance_cache,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_subject_genotype_block,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_subject_dosage_block,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_subject_posterior_prob_block,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_packedextendedhap,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, *buffer_subject_haplotype_weight,"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(int)*ref_haplotypes),"kernel_impute_genotype_guide"); // for root mapping
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(float)*ref_haplotypes),"kernel_impute_genotype_guide"); // for loading haplotype frequencies
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len),"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE),"kernel_impute_genotype_guide"); // prob0
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE),"kernel_impute_genotype_guide"); //prob1
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE),"kernel_impute_genotype_guide"); //prob2
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE),"kernel_impute_genotype_guide"); //prob3
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(float)*4),"kernel_impute_genotype_guide");
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(float)*g_max_haplotypes),"kernel_impute_genotype_guide"); // hap marginals
      setArg(kernel_impute_genotype_guide,arg, cl::__local(sizeof(float)*1),"kernel_impute_genotype_guide"); // normalizing constant
      kernelWorkGroupSize = kernel_impute_genotype_guide->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_guide");
      cerr<<"impute_genotype_guide kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;

    cerr<<"GPU kernel arguments assigned.\n";

    // TRANSFER ANY DATA HERE
    err = commandQueue->enqueueWriteBuffer(*buffer_genotype_imputation, CL_TRUE, 0,sizeof(int), &g_genotype_imputation, NULL, NULL );
    clSafe(err, "write imputation mode");
    err = commandQueue->enqueueWriteBuffer(*buffer_haploid_arr, CL_TRUE, 0,sizeof(int)*g_people, haploid_arr, NULL, NULL );
    clSafe(err, "write haploid arr");
#endif
  }
  cerr<<"Buffers initialized\n";
}

void GuidedMendelGPU::init_window_opencl(){
  MendelGPU::init_window_opencl();
  if (run_gpu){
    #ifdef USE_GPU
    #endif
  }
}
