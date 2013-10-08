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
  int max_geno=g_genotype_imputation?3:4;
  int c_snp_offset = g_center_snp_start;
  int c_snp_start = g_center_snp_start - g_left_marker;
  int c_snp_end = c_snp_start+g_flanking_snps;
  bool debug_geno = false;
  bool debug_dosage = c_snp_start==-50;
  bool debug_posterior = false;
  bool debug_pen = false;
  cerr<<"Imputing from "<<c_snp_start<<" to "<<(c_snp_end-1)<<" with offset "<<c_snp_offset<<endl;

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
    int len = c_snp_end-c_snp_start;
    for(int j=0;j<len;++j){ // loop over SNPs j
      if(debug_posterior){
        for(int i=0;i<g_people;++i){
          cout<<"GPU_POSTERIOR:\t"<<j;
          cout<<"\t"<<i;
          for(int k=0;k<max_geno;++k){
            cout<<"\t"<<
            subject_posterior_prob_block[i*g_flanking_snps*4+j*4+k];
          }
          cout<<endl;
        }
      }
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_posterior_file<<j;
        for(int i=0;i<g_people;++i){
          for(int k=0;k<max_geno;++k){
            ofs_posterior_file<<"\t"<<
            subject_posterior_prob_block[i*g_flanking_snps*4+j*4+k];
          }
        }
        ofs_posterior_file<<endl;
      }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
        for(int i=0;i<g_people;++i){
          ofs_posterior_file<<"GPU_POSTERIOR\t"<<i<<"\t"<<j;
          for(int k=0;k<max_geno;++k){
            ofs_posterior_file<<"\t"<<subject_posterior_prob_block[i*g_flanking_snps*4+j*4+k];
          }
          ofs_posterior_file<<endl;
        }
      }
    }
    int subject_geno_block[g_people*g_flanking_snps];
    readFromBuffer(buffer_subject_genotype_block, g_people*g_flanking_snps,subject_geno_block,"buffer_subject_genotype_block");
    cerr<<"Personloop: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
      len = c_snp_end-c_snp_start;
      for(int j=0;j<len;++j){
        int current_snp = c_snp_offset+c_snp_start+j;
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_genotype_file<<current_snp<<"\t";
        }
        for(int i=0;i<g_people;++i){
          if(debug_geno) cout<<"GPUGENO: "<<(c_snp_offset+c_snp_start+j)<<" "<<i<<" "<<subject_geno_block[i*g_flanking_snps+j]<<endl;
          if (outfile_format.compare(FORMAT_MEC)==0){
            ofs_genotype_file<<subject_geno_block[i*g_flanking_snps+j];
          }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
            ofs_genotype_file<<"GPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_geno_block[i*g_flanking_snps+j]<<endl;
          }else{
       //     ofs_genotype_file<<" "<<outfile_format;
          }
        }
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_genotype_file<<endl;
        }
      }
//    for(int i=0;i<g_people;++i){
//      ofs_genotype_file<<"GPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_genotypes[i]<<endl;
//    }
    if (g_genotype_imputation){
      float subject_dosage_block[g_people*g_flanking_snps];
      readFromBuffer(buffer_subject_dosage_block, g_people*g_flanking_snps,subject_dosage_block,"buffer_subject_dosage_block");
      int len = c_snp_end-c_snp_start;
      for(int j=0;j<len;++j){
        int current_snp = c_snp_offset+c_snp_start+j;
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_dosage_file<<current_snp<<"\t";
        }
        for(int i=0;i<g_people;++i){
          if (debug_dosage) cout<<"GPUDOSE: "<<current_snp<<" "<<i<<" "<<subject_dosage_block[i*g_flanking_snps+j]<<endl;
          if (outfile_format.compare(FORMAT_MEC)==0){
            char chardose = 
            (int)(subject_dosage_block[i*g_flanking_snps+j]*10)+'A';
            ofs_dosage_file<<chardose;
            //ofs_dosage_file<<subject_dosage_block[i*g_flanking_snps+j]<<endl;
          }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
            ofs_dosage_file<<"GPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosage_block[i*g_flanking_snps+j]<<endl;
          }
        }
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_dosage_file<<endl;
        }
        float rsq = compute_rsq(subject_dosage_block, g_flanking_snps,j);
        ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
      }
    }
//      for(int i=0;i<g_people;++i){
//        //float dose = 0;
//        //for(int j=0;j<max_geno;++j){
//        //  dose+=j*subject_posterior[i*4+j]/denom;
//        //}
//        ofs_dosage_file<<"GPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosages[i]<<endl;
//        if (debug_dosage) cout<<"GPU DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosages[i]<<endl;
//      }
//      float rsq = compute_rsq(subject_dosages);
//      ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
//    }
//    cerr<<"Elapsed time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
//    #endif
#endif
  }
}
