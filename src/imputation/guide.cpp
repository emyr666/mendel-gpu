#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"


void GuidedMendelGPU::impute_genotypes(){
  prep_impute_genotypes_guide();
  impute_genotypes_guide();
}

void GuidedMendelGPU::prep_impute_genotypes_guide(){
  int left_marker = g_left_marker;
  cerr<<"begin prep impute_genotypes at left SNP "<<left_marker<<endl;
  int max_geno=g_genotype_imputation?3:4;
  int curhapindex = 0;
  if (debug_mode) ofs_debug_haplotype_file<<"Extended haplotypes for left marker "<<left_marker<<endl;
  for(int i=0;i<g_haplotypes;++i){
    string hapstring = hapstr(g_haplotype+i*g_max_window,g_markers);
    if (full_hap_window.find(hapstring)==full_hap_window.end()){
      cerr<<"Orphaned hapstr "<<hapstring<<endl;
      exit(1);
    }
    set<string> extendedhaps = full_hap_window[hapstring];
    float freq = 1.*g_frequency[i]/(extendedhaps.size());
    if (debug_mode) ofs_debug_haplotype_file<<"Compact haplotype: "<<hapstring<<endl;
    float totalcounts = 0;
    int j = 0;
    float extendedfreq[extendedhaps.size()];
    for(set<string>::iterator it  = extendedhaps.begin();
    it!=extendedhaps.end();it++){
       string str = *it;
       if (ref_hap_counts.find(str)==ref_hap_counts.end()){
         cerr<<"Exception: couldn't find refhap counts for "<<str<<endl;
       }
       extendedfreq[j] = ref_hap_counts[str];
       totalcounts+=1.*extendedfreq[j];
       ++j;
    }
    for(int j=0;j<extendedhaps.size();++j) extendedfreq[j] = extendedfreq[j]/totalcounts*g_frequency[i]*1.;
    j = 0;
    for(set<string>::iterator it  = extendedhaps.begin();
    it!=extendedhaps.end();it++){
       string str = *it;
       for(int k=0;k<extended_markers;++k){
         //cerr<<"Making extended hap for "<< curhapindex*g_max_window+k<<"\n";
         extended_haplotype[curhapindex*g_max_window+k] = hapchar2int(str[k]);
       }
       extended_frequency[curhapindex] = extendedfreq[j];
       if (debug_mode) ofs_debug_haplotype_file<<"Extended hap: "<<str<<" base freq "<<g_frequency[i]<<","<<ref_hap_counts[str]<<", extendedfreq: "<<extendedfreq[j]<<endl;
       extended_root_mapping[curhapindex] = i;
       ++j;
       ++curhapindex;
    }
  }
  if (curhapindex!=extended_haplotypes) {
    cerr<<"Assertion failed curhapindex ("<<curhapindex<<")!=extended_haplotypes("<<extended_haplotypes<<")\n ";
    exit(1);
  }
  if (run_gpu){
    prep_impute_genotypes_guide_opencl();
  }
  cerr<<"done prep impute_genotypes at left SNP "<<left_marker<<endl;
}

void GuidedMendelGPU::impute_genotypes_guide(){
  int max_geno=g_genotype_imputation?3:4;
  int c_snp_offset = g_center_snp_start;
  int c_snp_start = g_center_snp_start - g_left_marker;
  int c_snp_end = c_snp_start+g_flanking_snps;
  bool debug_geno = false;
  bool debug_dosage = c_snp_start==-50;
  bool debug_posterior = false;
  bool debug_pen = false;
  cerr<<"Imputing from "<<c_snp_start<<" to "<<(c_snp_end-1)<<" with offset "<<c_snp_offset<<endl;
  for(int hapindex=0;hapindex<extended_haplotypes;++hapindex){
    compress_extendedhap(hapindex,c_snp_start,c_snp_end,extended_haplotype);
    //compress(hapindex,extended_markers,extended_haplotype);
  }
  cerr<<"Compressed haplotypes\n";
  for(int i=0;i<extended_haplotypes;++i){
    //extended_center_dosage[i] = extended_haplotype[i*g_max_window+c_snp];
  }
  if (run_gpu){
    impute_genotypes_guide_opencl();
  }
  if(run_cpu){
    float subject_dosages[g_people * g_flanking_snps];
    float subject_posterior_prob[g_people * g_flanking_snps * 4];
    memset(subject_posterior_prob,0,sizeof(float)*g_people*g_flanking_snps * 4);
    double start = clock();
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
      for(int j=0;j<extended_haplotypes;++j){
        //cerr<<"subject hap weight for hap "<<extended_root_mapping[j]<<" is "<<subject_haplotype_weight[i*g_max_haplotypes+extended_root_mapping[j]] <<" for subject "<<i<<endl;
        if(1==1){
          //extended_root_mapping[j]]>0){ 
        //if(subject_haplotype_weight[i*g_max_haplotypes+
        //  extended_root_mapping[j]]>0){ 
          //decompress(j,extended_markers,extended_haplotype);
          for(int k=j;k<extended_haplotypes;++k){
            if((!haploid ||  extended_root_mapping[j] == 
            extended_root_mapping[k]) && subject_haplotype_weight[i*
            g_max_haplotypes+extended_root_mapping[k]]>0){ 
              //cerr<<"subject hap weight for hap "<<k<<" is pos for subject "<<i<<endl;
              //decompress(k,extended_markers,extended_haplotype);
              float penetrance =
              penetrance_cache[i*penetrance_matrix_size+extended_root_mapping[j]*
              g_max_haplotypes+extended_root_mapping[k]];
              if (penetrance>0){
                for(int c_snp = c_snp_start;c_snp<c_snp_end; ++c_snp){
                  int packedsite = c_snp - c_snp_start;
                  int current_snp = c_snp_offset + c_snp;
//cerr<<"packed site and current snp: "<<packedsite<<","<<current_snp<<endl;
                  int m;
                  if(g_genotype_imputation){
                    m = (((int)packedextendedhap[j*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) + (((int)packedextendedhap[k*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) ;
                  }else{
                    m = 2*(((int)packedextendedhap[j*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) + (((int)packedextendedhap[k*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) ;
                  }
                  float freq = extended_frequency[j]*extended_frequency[k];
                  if (j!=k) freq*=2;
                  float p = freq*penetrance;
                  //p = 1;
                  subject_posterior_prob[i*g_flanking_snps*4+(packedsite)*4+m]+=p;
                  if (debug_pen) cerr<<"i,j,k,m,pen,freq:"<<i<<","<<j<<","<<k<<","<<m<<","<<penetrance<<","<<freq<<endl;
                  //if (debug_pen) cerr<<"packedsite: "<<packedsite<<endl;
                }
              }
            }
          }
        }
      }
      //cerr<<" "<<i;
    }
    cerr<<"Personloop: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    //cerr<<endl;
      //cout<<"CPU POSTERIOR "<<current_snp<<","<<i;
      //for(int j=0;j<max_geno;++j){
       // cout<<" "<<posterior_prob[j];
      //}
      //cout<<endl;
    for(int c_snp = c_snp_start;c_snp<c_snp_end; ++c_snp){
      int current_snp = c_snp_offset + c_snp;
//      if (outfile_format.compare(FORMAT_MEC)==0){
//        ofs_genotype_file<<current_snp<<"\t";
//        ofs_dosage_file<<current_snp<<"\t";
//        ofs_posterior_file<<current_snp;
//      }
      float dosages[g_people];
      float posteriors[g_people*max_geno];
      int genotypes[g_people];
      for(int i=0;i<g_people;++i){
        // compute normalizing constant
        float denom = gf_epsilon;
        float maxval = 0;
        int genotype = 0;
        for(int m=0;m<max_geno;++m){
          float p = subject_posterior_prob[i*g_flanking_snps*4+
          (c_snp-c_snp_start)*4+m];
          denom+=p;
          if (p>maxval) {
            maxval = p;
            genotype = m;
          }
        }
        if (debug_posterior){
           //cout<<"denom: "<<denom<<endl;
           cout<<"CPU_POSTERIOR: "<<current_snp<<" "<<i;
           for(int k=0;k<max_geno;++k){
             cout<<" "<<subject_posterior_prob[i*g_flanking_snps*4+(c_snp-c_snp_start)*4+k]/denom;
           }
           cout<<endl;
        }
        //cerr<<"Person "<<i<<" denom is "<<denom<<endl;
        if (debug_geno) cout<<"CPUGENO:\t"<<i<<"\t"<<current_snp<<"\t"<<genotype<<endl;
        genotypes[i] = genotype;
//        if (outfile_format.compare(FORMAT_MEC)==0){
//          ofs_genotype_file<<genotype;
//        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
//          ofs_genotype_file<<"CPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<genotype<<endl;
//        }
        dosages[i] = 0;
        for(int j=0;j<max_geno;++j){
          int snpindex = c_snp-c_snp_start;
          float p = subject_posterior_prob[i*g_flanking_snps*4+
          (c_snp-c_snp_start)*4+j]/denom;
          posteriors[i*max_geno+j] = p;
          dosages[i]+=j*p;
        }
        subject_dosages[i*g_flanking_snps+(c_snp-c_snp_start)] = dosages[i];
 

//        float dose=0;
//        if (outfile_format.compare(FORMAT_MEC)==0){
//          for(int j=0;j<max_geno;++j){
//            float p = subject_posterior_prob[i*g_flanking_snps*4+
//            (c_snp-c_snp_start)*4+j]/denom;
//            //ofs_posterior_file<<"\t"<<p;
//            dose+=j*p;
//          }
//          char chardose =  (int)(dose*10)+'A';
//          //ofs_dosage_file<<chardose;
//        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
//          //ofs_posterior_file<<"CPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
//          for(int j=0;j<max_geno;++j){
//            float p = subject_posterior_prob[i*g_flanking_snps*4+
//            (c_snp-c_snp_start)*4+j]/denom;
//            if(debug_mode)ofs_posterior_file<<"\t"<< p;
//            dose+=j*p;
//          }
//          //ofs_dosage_file<<"CPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<dose<<endl;
//          //ofs_posterior_file<<endl;
//        }
      }
      //if (outfile_format.compare(FORMAT_MEC)==0){
        //ofs_genotype_file<<endl;
        //ofs_dosage_file<<endl;
      //}
      float rsq = compute_rsq(subject_dosages,
      g_flanking_snps,(c_snp-c_snp_start));
      if (max_geno==3){
        io_manager->writeDosage(current_snp,dosages,g_people);
        io_manager->writeGenotype(current_snp,genotypes,g_people);
      }
      io_manager->writePosterior(max_geno,current_snp,posteriors,g_people);

      //ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
    }
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_posterior|| debug_dosage||debug_geno) exit(1);
  if (debug_pen) {
    if (debug_mode) ofs_debug_haplotype_file.close();
    exit(1);
  }
}

