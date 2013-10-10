#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif


void GuidedMendelGPU::prep_impute_genotypes_guide_opencl(){
  if (run_gpu){
    #ifdef USE_GPU
    writeToBuffer(buffer_haplotypes, 1, &g_haplotypes, "buffer_haplotypes" );
    writeToBuffer(buffer_extended_haplotypes, 1, &extended_haplotypes, "buffer_extended_haplotypes" );
    writeToBuffer(buffer_extended_root_mapping, ref_haplotypes, extended_root_mapping, "buffer_extended_root_mapping" );
    writeToBuffer(buffer_extended_frequency, ref_haplotypes, extended_frequency, "buffer_extended_frequency" );
    #endif
  }
}

void GuidedMendelGPU::impute_genotypes_guide_opencl(){
  int c_snp_offset = g_center_snp_start;
  int c_snp_start = g_center_snp_start - g_left_marker;
  int c_snp_end = c_snp_start+g_flanking_snps;
  bool debug_geno = false;
  bool debug_dosage = c_snp_start==-50;
  bool debug_posterior = false;
  bool debug_pen = false;
  cerr<<"Imputing from "<<c_snp_start<<" to "<<(c_snp_end-1)<<" with offset "<<c_snp_offset<<endl;

  int max_geno = g_genotype_imputation?3:4;
  if (run_gpu){
#ifdef USE_GPU
    double start = clock();
    int last_site = c_snp_end-c_snp_start; 
    writeToBuffer(buffer_center_snp_end, 1, &last_site, "buffer_center_snp_end" );
    writeToBuffer(buffer_packedextendedhap, ref_haplotypes*packedextendedhap_len, packedextendedhap, "buffer_packedextendehap");
    cerr<<"Launching impute geno with max site "<<last_site<<endl;
    runKernel("kernel_impute_genotype_guide",kernel_impute_genotype_guide,g_people*BLOCK_WIDTH_IMPUTE_GUIDE,last_site,1,BLOCK_WIDTH_IMPUTE_GUIDE,1,1);
    // READ GENOTYPE POSTERIORS
    float subject_posterior_prob_block[g_people * g_flanking_snps * 4];
    readFromBuffer(buffer_subject_posterior_prob_block, g_people*g_flanking_snps*4,subject_posterior_prob_block,"buffer_subject_posterior_prob_block");
    int subject_geno_block[g_people*g_flanking_snps];
    readFromBuffer(buffer_subject_genotype_block, g_people*g_flanking_snps,subject_geno_block,"buffer_subject_genotype_block");
    float subject_dosage_block[g_people*g_flanking_snps];
    readFromBuffer(buffer_subject_dosage_block, g_people*g_flanking_snps,subject_dosage_block,"buffer_subject_dosage_block");
    for(int j=g_center_snp_start;j<=g_center_snp_end;++j){
      int col = j-g_center_snp_start;
      float posteriors[g_people*max_geno];
      int genotypes[g_people];
      float dosages[g_people];
      for(int i=0;i<g_people;++i){
        if (debug_posterior) cout<<"GPU_POSTERIOR:\t"<<j<<"\t"<<i;
        for(int k=0;k<max_geno;++k){
          posteriors[i*max_geno+k] = 
          subject_posterior_prob_block[i*g_flanking_snps*4+col*4+k];
          if (debug_posterior) cout<<"\t"<<posteriors[i*max_geno+k];
        }
        if (debug_posterior) cout<<endl;
        genotypes[i] = subject_geno_block[i*g_flanking_snps+col];
        dosages[i] = subject_dosage_block[i*g_flanking_snps+col];
      }
      io_manager->writePosterior(max_geno,j,posteriors,g_people);
      io_manager->writeGenotype(j,genotypes,g_people);
      if (g_genotype_imputation){
        io_manager->writeDosage(j,dosages,g_people);
        float rsq = compute_rsq(dosages,1,0);
        io_manager->writeQuality(j,rsq);
      }
    }
#endif
  }
}
