#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"

void DenovoMendelGPU::impute_genotypes(){
  //return;
  if (run_gpu){
    #ifdef USE_GPU
    impute_genotypes_opencl();
    #endif
  }
  int c_snp = g_center_snp_start-g_left_marker;
  int current_snp = g_center_snp_start;
  cerr<<"begin diploid impute_genotypes at SNP "<<current_snp<<", center: "<<c_snp<<endl;
  bool debug_dosage = false;
  bool debug_geno = false;
  bool debug_posterior = current_snp==-94;
  bool debug_pen = false; 
  bool debug_subject = 0;
  int geno_dim = g_genotype_imputation?3:4;
  cerr<<"Max geno "<<geno_dim<<endl;
  for(int i=0;i<g_max_haplotypes;++i){
    center_dosage[i] = g_haplotype[i*g_max_window+c_snp];
  }
  float subject_dosages[g_people];
  if(run_cpu){
    //if (g_likelihood_mode==Config::LIKELIHOOD_MODE_READS && 
    //  window_polymorphisms[c_snp]){
    float dosages[g_people];
    int genotypes[g_people*geno_dim];
    float posteriors[g_people*geno_dim];
//    if (outfile_format.compare(FORMAT_MEC)==0){
//      ofs_dosage_file<<current_snp<<"\t";
//      ofs_genotype_file<<current_snp<<"\t";
//      ofs_posterior_file<<current_snp;
//    }
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
      float best_prob = 0;
      float posterior_prob[4];
      int best_pair[2];
      best_pair[0] = best_pair[1] = 0;
      for(int j=0;j<4;++j) posterior_prob[j] = gf_epsilon;
      for(int j=0;j<g_max_haplotypes;++j){
        if(g_active_haplotype[j] ){
          for(int k=j;k<g_max_haplotypes;++k){
            if(g_active_haplotype[k] && (!haploid || j==k)){
              int m;
              float penetrance = penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
              //float penetrance = 0;
              if (penetrance>0 ){
                //cerr<<"Subject "<<i<<" penetrance at haps "<<j<<" and "<<k<<" :" <<penetrance<<endl;
                float freq = g_frequency[j]*g_frequency[k];
                if (j!=k) freq*=2;
                float p = freq*penetrance;
                if (p>best_prob){
                  best_prob = p;
                  best_pair[0] = j;
                  best_pair[1] = k;
                }
                if(g_genotype_imputation){
                  m = center_dosage[j] + center_dosage[k];
                }else{
                  m = 2*center_dosage[j] + center_dosage[k];
                }
                if (debug_pen && debug_subject==i) cerr<<"i,j,k,m,pen,freq:"<<i<<","<<j<<","<<k<<","<<m<<","<<penetrance<<","<<freq<<endl;
                posterior_prob[m]+=p;
              }
            }
          }
        }
      }
      float denom = 0;
      for(int m=0;m<geno_dim;++m){
        //if (posterior_prob[m]<gf_epsilon) posterior_prob[m] = gf_epsilon;
        denom+=posterior_prob[m];
      }
      for(int j=0;j<geno_dim;++j){
        posterior_prob[j]/=denom;
      }
      if (debug_posterior){
        cout<<"CPU_POSTERIOR "<<current_snp<<","<<i;
        for(int j=0;j<geno_dim;++j){
          cout<<" "<<posterior_prob[j];
        }
        cout<<endl;
      }
      if (denom==0) {
        cerr<<"divide by zero denom\n"; 
        for(int m=0;m<geno_dim;++m){
          cerr<<"posterior for geno "<<m<<": "<<posterior_prob[m]<<endl;
        }
        cerr<<"hap1,hap2:pen,logpen,freq,genotype\n";
        for(int j=0;j<g_max_haplotypes;++j){
          int start = g_genotype_imputation?j:0;
          for(int k=start;k<g_max_haplotypes;++k){
            if(g_active_haplotype[j] && g_active_haplotype[k]){
              int m = g_haplotype[j*g_max_window+ c_snp] + 
              g_haplotype[k*g_max_window+ c_snp];
              float logpenetrance = logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
              float penetrance = penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
              float freq = g_frequency[j]*g_frequency[k];
              cerr<<j<<","<<k<<":"<<penetrance<<","<<logpenetrance<<","<<freq<<","<<m<<endl;
            }
          }
        }
        exit(1);
      }
      dosages[i] = 0;
      for(int j=0;j<geno_dim;++j){
        float & p = posterior_prob[j];
        posteriors[i*geno_dim+j] = p;
        dosages[i]+=j*p;
      }
      if (debug_dosage){
          cout<<"CPU_DOSAGE:\t"<<i<<"\t"<<current_snp<<"\t"<<dosages[i]<<endl;
      }
      //subject_dosages[i] = dose;
      float maxval = 0;
      for(int j=0;j<4;++j){
        if (posterior_prob[j]>maxval) maxval = posterior_prob[j];
      }
      if (maxval>0){
        int geno1= -1;
        float max_prob = posterior_prob[0];
        int geno2 = 0;
        for(int h=1;h<4;++h){
          if (posterior_prob[h]-max_prob>gf_epsilon*max_prob){
            max_prob = posterior_prob[h];
            geno2 = h;
          }
        }  
        if (debug_geno){
          cout<<"CPU GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<geno1<<","<<geno2<<endl;
        }
        genotypes[i] = geno2;
      }
    }
    if (geno_dim==3){
      io_manager->writeDosage(current_snp,dosages,g_people);
      io_manager->writeGenotype(current_snp,genotypes,g_people);
    }
    io_manager->writePosterior(geno_dim,current_snp,posteriors,g_people);
    float rsq = compute_rsq(dosages,1,0);
    io_manager->writeQuality(current_snp,rsq);
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_dosage || debug_posterior|| debug_geno) exit(1);
  if (debug_pen) {
    exit(1);
  }
}
