#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"


void GuidedMendelGPU::impute_genotypes(){
  prep_impute_genotypes_guide();
  impute_genotypes_guide();
}

void GuidedMendelGPU::prep_impute_genotypes_guide(){
  bool debug_mode = false;
  int left_marker = g_left_marker;
  cerr<<"begin prep impute_genotypes at left SNP "<<left_marker<<endl;
  int curhapindex = 0;
  if (debug_mode) cerr<<"Extended haplotypes for left marker "<<left_marker<<endl;
  for(int i=0;i<g_haplotypes;++i){
    string hapstring = hapstr(g_haplotype+i*g_max_window,g_markers);
    if (full_hap_window.find(hapstring)==full_hap_window.end()){
      cerr<<"Orphaned hapstr "<<hapstring<<endl;
      exit(1);
    }
    set<string> extendedhaps = full_hap_window[hapstring];
    float freq = 1.*g_frequency[i]/(extendedhaps.size());
    if (debug_mode) cerr<<"Compact haplotype: "<<hapstring<<endl;
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
       if (debug_mode) cerr<<"Extended hap: "<<str<<" base freq "<<g_frequency[i]<<","<<ref_hap_counts[str]<<", extendedfreq: "<<extendedfreq[j]<<endl;
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
  int geno_dim = g_genotype_imputation?3:4;
  if(run_cpu){
    float subject_posterior_block[g_people*g_flanking_snps*4];
    for(int i=0;i<g_people*g_flanking_snps*4;++i) subject_posterior_block[i] = gf_epsilon;
    double start = clock();
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
      for(int j=0;j<extended_haplotypes;++j){
        for(int k=j;k<extended_haplotypes;++k){
          if(!haploid ||  extended_root_mapping[j] == 
          extended_root_mapping[k]){ 
            float penetrance =
            penetrance_cache[i*penetrance_matrix_size+extended_root_mapping[j]*
            g_max_haplotypes+extended_root_mapping[k]];
            if (penetrance>0){
              for(int c_snp = c_snp_start;c_snp<c_snp_end; ++c_snp){
                int packedsite = c_snp - c_snp_start;
                int current_snp = c_snp_offset + c_snp;
                int m;
                if(g_genotype_imputation){
                  m = (((int)packedextendedhap[j*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) + (((int)packedextendedhap[k*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) ;
                }else{
                  m = 2*(((int)packedextendedhap[j*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) + (((int)packedextendedhap[k*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) ;
                }
                float freq = extended_frequency[j]*extended_frequency[k];
                if (j!=k) freq*=2;
                float p = freq*penetrance;
                subject_posterior_block[i*g_flanking_snps*4+(packedsite)*4+m]+=p;
                if (debug_pen) cerr<<"i,j,k,m,pen,freq:"<<i<<","<<j<<","<<k<<","<<m<<","<<penetrance<<","<<freq<<endl;
              }
            }
          }
        }
      }
    }
    cerr<<"Personloop: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    float posteriors[g_people*geno_dim];
    float dosages[g_people];
    int genotypes[g_people];
    for(int j=g_center_snp_start;j<=g_center_snp_end;++j){
      int col = j-g_center_snp_start;
      for(int i=0;i<g_people;++i){
        // compute normalizing constant
        float denom = 0;
        float maxval = 0;
        for(int m=0;m<geno_dim;++m){
          float p = subject_posterior_block[i*g_flanking_snps*4+
          col*4+m];
          denom+=p;
          if (p>maxval) {
            maxval = p;
            genotypes[i] = m;
          }
        }
        for(int m=0;m<geno_dim;++m){
          posteriors[i*geno_dim+m] = subject_posterior_block[i*g_flanking_snps*4+col*4+m]/denom;
        }
        if (debug_posterior){
         //cout<<"denom: "<<denom<<endl;
          cout<<"CPU_POSTERIOR: "<<j<<" "<<i;
          for(int k=0;k<geno_dim;++k){
            cout<<" "<<posteriors[i*geno_dim+k];
          }
          cout<<endl;
        }
        dosages[i] = 0;
        for(int j=0;j<geno_dim;++j){
          dosages[i] +=j*posteriors[i*geno_dim+j];
        }
      }
      if (g_genotype_imputation){
        io_manager->writeDosage(j,dosages,g_people);
        float rsq = compute_rsq(dosages,1,0);
        io_manager->writeQuality(j,rsq);
      }
      io_manager->writePosterior(geno_dim,j,posteriors,g_people);
      io_manager->writeGenotype(j,genotypes,g_people);
    }
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_posterior|| debug_dosage||debug_geno) exit(1);
  if (debug_pen) {
    exit(1);
  }
}

