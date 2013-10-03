#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include<CL/cl.hpp>
#include"../clsafe.h"
#endif

void GuidedMendelGPU::init_opencl(){
  MendelGPU::init_opencl();
  if(run_gpu){
#ifdef USE_GPU
    // CREATE ALL KERNELS HERE
    kernel_precompute_penetrance = new cl::Kernel(*program,"precompute_penetrance",&err);
    clSafe(err,"creating kernel precompute");
    kernel_precompute_penetrance_haploid = new cl::Kernel(*program,"precompute_penetrance_haploid",&err);
    clSafe(err,"creating kernel precompute");

    kernel_impute_genotype_guide = new cl::Kernel(*program,"impute_genotype_guide",&err);
    clSafe(err,"creating kernel impute_genotype_guide");
    kernel_impute_genotype_guide_haploid = new cl::Kernel(*program,"impute_genotype_guide_haploid",&err);
    clSafe(err,"creating kernel impute_genotype_guide_haploid");
    // CREATE ALL BUFFERS HERE
    buffer_markers = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer markers");
    buffer_haplotypes = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer haplotypes");
    buffer_extended_haplotypes = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer extended haplotypes");
    //buffer_left_marker = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating buffer left marker");
    //buffer_prev_left_marker = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating buffer prev left marker");
    buffer_haplotype = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_max_window * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype buffer");
    //buffer_region_snp_offset = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating snp offset buffer");
    //buffer_region_snp_penetrance = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people * geno_dim * g_max_region_size , NULL, &err);
    //clSafe(err,"creating snp penetrance buffer");
    buffer_frequency = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype frequency buffer");
    buffer_subject_haplotype_weight = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating subject haplotype weights buffer");
    buffer_haplotype_weight = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype weights buffer");
    //buffer_center_snp = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating buffer center marker");
    //buffer_center_snp_start = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating buffer center start marker");
    buffer_center_snp_end = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer center end marker");
    //buffer_center_snp_offset = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating buffer center offset marker");
    //buffer_subject_genotype = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_people, NULL, &err);
    //clSafe(err,"creating buffer subject genotype");
    buffer_subject_genotype_block = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_people*g_flanking_snps, NULL, &err);
    clSafe(err,"creating buffer subject genotype block");
    //buffer_subject_dosage = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people, NULL, &err);
    //clSafe(err,"creating buffer subject dosage");
    buffer_subject_dosage_block = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people*g_flanking_snps, NULL, &err);
    clSafe(err,"creating buffer subject dosage block");
    //buffer_subject_posterior_prob = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*4, NULL, &err);
    //clSafe(err,"creating buffer subject_posterior_prob");
    buffer_subject_posterior_prob_block = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*g_flanking_snps*4, NULL, &err);
    clSafe(err,"creating buffer subject_posterior_prob block");
    buffer_haploid_arr = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_people, NULL, &err);
    clSafe(err,"creating buffer haploid_arr");
    if (geno_dim==PHASED_INPUT){
      buffer_logpenetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*2*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer logpenetrance_cache");
      buffer_penetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*2*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer penetrance_cache");
      buffer_frequency_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer frequency_cache");
    }else if (geno_dim==UNPHASED_INPUT){
      buffer_logpenetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer logpenetrance_cache");
      buffer_penetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer penetrance_cache");
      buffer_frequency_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer frequency_cache");
    }
    buffer_max_penetrance = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people, NULL, &err);
    clSafe(err,"creating buffer max penetrance");
    buffer_genotype_imputation = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int), NULL, &err);
    clSafe(err,"creating buffer geno impute");
    //buffer_twin_hap_index = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    //clSafe(err,"creating buffer twin index");
    buffer_active_haplotype = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    clSafe(err,"creating buffer active hap");
    buffer_iteration = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*1, NULL, &err);
    clSafe(err,"creating buffer iteration");
    if (g_haplotype_mode==HAPLOTYPE_MODE_GUIDE){
      buffer_packedhap = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(packedhap_t) *g_max_haplotypes * packedhap_len, NULL, &err);
      clSafe(err,"creating packedhaplotype buffer");
      buffer_packedextendedhap = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(packedhap_t) *ref_haplotypes * packedextendedhap_len, NULL, &err);
      clSafe(err,"creating packedextendedhaplotype buffer");
      buffer_extended_root_mapping = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*ref_haplotypes, NULL, &err);
      clSafe(err,"creating buffer extended root mapping");
      buffer_extended_snp_mapping = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_window, NULL, &err);
      clSafe(err,"creating buffer extended snp mapping");
      buffer_extended_frequency = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * ref_haplotypes, NULL, &err);
      clSafe(err,"creating extended frequency buffer");
    }else if (g_haplotype_mode==HAPLOTYPE_MODE_DENOVO){
    }
    cerr<<"GPU Buffers created\n";
    // SET KERNEL ARGUMENTS HERE
    int arg;
    int kernelWorkGroupSize;
    if (g_haplotype_mode==HAPLOTYPE_MODE_DENOVO){
    } else if (g_haplotype_mode == HAPLOTYPE_MODE_GUIDE){
      // PRECOMPUTE PENETRANCE FOR THE GUIDE HAP APPROACH
      arg = 0;
      err = kernel_precompute_penetrance->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, g_max_window);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      //err = kernel_precompute_penetrance->setArg(arg++, g_max_region_size);
      //clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, logpenetrance_threshold);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, packedhap_len);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      //err = kernel_precompute_penetrance->setArg(arg++, *buffer_region_snp_offset);
      //clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_markers);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_left_marker);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_haploid_arr);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_packedhap);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_haplotype);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      //err = kernel_precompute_penetrance->setArg(arg++, *buffer_region_snp_penetrance);
      //clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_logpenetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_extended_snp_mapping);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(packedhap_t)*g_max_haplotypes*packedhap_len));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(float)*g_max_window*geno_dim));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      cerr<<"Total args for precompute pen: "<<arg<<endl;
      kernelWorkGroupSize = kernel_precompute_penetrance->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel precompute_penetrance");
      cerr<<"precompute_penetrance kernel work group size is "<<kernelWorkGroupSize<<endl;
      arg = 0;
      err = kernel_precompute_penetrance_haploid->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, g_max_window);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      //err = kernel_precompute_penetrance_haploid->setArg(arg++, g_max_region_size);
      //clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, logpenetrance_threshold);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, packedhap_len);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      //err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_region_snp_offset);
      //clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_markers);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_left_marker);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_packedhap);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_haplotype);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      //err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_region_snp_penetrance);
      //clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_logpenetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_extended_snp_mapping);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(packedhap_t)*g_max_haplotypes*packedhap_len));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_window*geno_dim));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(float)*2*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      cerr<<"Total args for precompute pen: "<<arg<<endl;
      kernelWorkGroupSize = kernel_precompute_penetrance_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel precompute_penetrance_haploid");
      cerr<<"precompute_penetrance_haploid kernel work group size is "<<kernelWorkGroupSize<<endl;
      arg = 0;
      err = kernel_impute_genotype_guide->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, epsilon);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, packedextendedhap_len);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, g_flanking_snps);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_genotype_imputation);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_extended_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_center_snp_end);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_haploid_arr);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_extended_root_mapping);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_extended_frequency);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_subject_genotype_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_subject_dosage_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_subject_posterior_prob_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_packedextendedhap);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_subject_haplotype_weight);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(int)*ref_haplotypes)); // for root mapping
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*ref_haplotypes)); // for loading haplotype frequencies
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); // prob0
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob1
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob2
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob3
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*4));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // hap marginals
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*1)); // normalizing constant
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      kernelWorkGroupSize = kernel_impute_genotype_guide->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_guide");
      cerr<<"impute_genotype_guide kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;

      arg = 0;
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, epsilon);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, packedextendedhap_len);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, g_flanking_snps);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_genotype_imputation);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_center_snp_end);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_root_mapping);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_frequency);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_genotype_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_dosage_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_posterior_prob_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_packedextendedhap);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_haplotype_weight);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(int)*ref_haplotypes)); // for root mapping
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*ref_haplotypes)); // for loading haplotype frequencies
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); // prob0
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob1
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob2
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob3
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*4));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // hap marginals
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*1)); // normalizing constant
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      kernelWorkGroupSize = kernel_impute_genotype_guide_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_guide_haploid");
      cerr<<"impute_genotype_guide_haploid kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;

      arg = 0;
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, epsilon);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, packedextendedhap_len);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, g_flanking_snps);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_genotype_imputation);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_center_snp_end);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_root_mapping);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_frequency);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_genotype_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_dosage_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_posterior_prob_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_packedextendedhap);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_haplotype_weight);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(int)*ref_haplotypes)); // for root mapping
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*ref_haplotypes)); // for loading haplotype frequencies
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); // prob0
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob1
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob2
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob3
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*4));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // hap marginals
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*1)); // normalizing constant
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      kernelWorkGroupSize = kernel_impute_genotype_guide_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_guide_haploid");
      cerr<<"impute_genotype_guide_haploid kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;
    }


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
