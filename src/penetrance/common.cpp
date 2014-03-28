#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"


bool MendelGPU::check_mm_converged(){
  bool debug_freq = false;
  if(debug_freq)cerr<<"Doing MM updates\n";
  float d = 0;
  float e = 0;
  float t = 0;
  for(int j=0;j<g_max_haplotypes;++j){
    if (g_active_haplotype[j]){
      if(debug_freq) cerr<<"Current weight of hap "<<j<<" is "<<g_weight[j]<<endl;
      if (g_weight[j]<gf_epsilon) g_weight[j] = gf_epsilon;
      if (g_frequency[j] > g_delta) d+=g_weight[j];
      else e+=g_weight[j];
    }
  }
  t = d+e;
  float q;
  if (g_lambda*d*e>0){
    q = 2*e/(d+e+g_lambda+sqrt((d+e+g_lambda)*(d+e+g_lambda)-4*g_lambda*e));
    for(int j=0;j<g_max_haplotypes;++j){
      if (g_active_haplotype[j]){
        g_frequency[j] = g_frequency[j]>g_delta?g_weight[j]/(t-g_lambda*q):g_weight[j]/(t-g_lambda*q+g_lambda);
        if (debug_freq) cerr<<"1 Assigning new freq "<<j<<" of "<<g_frequency[j]<<endl;
        if(g_frequency[j]<gf_epsilon) g_frequency[j] = gf_epsilon;
        
      }
    }
  }else if(t>0){
    for(int j=0;j<g_max_haplotypes;++j){
      g_frequency[j] = g_weight[j]/t;
        if (debug_freq)cerr<<"2 Assigning new freq "<<j<<" of "<<g_frequency[j]<<endl;
    }
  }else{
    for(int j=0;j<g_max_haplotypes;++j){
      g_frequency[j] = 1./g_max_haplotypes;
        if (debug_freq)cerr<<"3 Assigning new freq "<<j<<" of "<<g_frequency[j]<<endl;
    }
  }
  t = 0;
  for(int j=0;j<g_max_haplotypes;++j){
    if (g_active_haplotype[j]){
      t+=fabs(g_frequency[j] - g_old_frequency[j]);
    }
  }
  t/=g_haplotypes;
  if(debug_freq)cerr<<"MM: T is "<<t<<endl;
  if (t<gf_convergence_criterion) return true;
  else{
    for(int j=0;j<g_max_haplotypes;++j){
      g_old_frequency[j] = g_frequency[j];
    }
  }
  return false;
}

void MendelGPU::compute_haplotype_weights(){
  if (run_gpu){
    compute_haplotype_weights_opencl();
  }
  bool b_impute = g_genotype_imputation;
  cerr<<"begin compute_weights at iteration: "<<gi_iteration<<"\n";
  //cerr<<"Geno dim is "<<geno_dim<<endl;
  //return;
  bool debug_personhap = false;
  //bool debug_personhap = gi_iteration==2;
  //bool debug_personhap = g_haplotypes<-4;
  int debug_person = 0;
  bool debug_weights = false;
  bool debug_freq = false;
  if(run_cpu){
    memset(g_weight,0,sizeof(float)*g_max_haplotypes);
    if (gi_iteration==0){
      memset(subject_haplotype_weight,0,sizeof(float)*g_people*g_max_haplotypes);
    }
    //for(int i=0;i<0;++i){
    double start = clock();
    if (debug_freq) cout<<"CPU freq: "<<endl;
    if (geno_dim==UNPHASED_INPUT){
      for(int j=0;j<g_max_haplotypes;++j){
        int start = b_impute?j:0;
        if (g_active_haplotype[j]){
          if (debug_freq) cout<<"Freq "<<j<<":"<<g_frequency[j];
          for(int k=start;k<g_max_haplotypes;++k){
            if (g_active_haplotype[k]){
              frequency_cache[j*g_max_haplotypes+k] = g_frequency[j]*g_frequency[k];
              if (j!=k) frequency_cache[j*g_max_haplotypes+k]*=2;
            }
          }
        }
        if (debug_freq) cout<<endl;
      }
    }
    for(int i=0;i<g_people;++i){
      //cerr<<"Person "<<i<<endl;
      int haploid = haploid_arr[i];
      float likelihood = 0;
      memset(g_current_weight,0,sizeof(float)*g_max_haplotypes);
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] && (gi_iteration==0 || subject_haplotype_weight[i*g_max_haplotypes+j]>0)){
          if(geno_dim==PHASED_INPUT){
            for(int parent=0;parent<2;++parent){
              float penetrance = penetrance_cache[i*2*g_max_haplotypes+
              2*j+parent];
              if (penetrance>0){
                float freq = g_frequency[j];
                float p = freq * penetrance;
                likelihood+=p;
                g_current_weight[j]+=p;
              }
            }
          }else if (geno_dim==UNPHASED_INPUT){
            for(int k=j;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k] &&
              (gi_iteration==0 || 
              subject_haplotype_weight[i*g_max_haplotypes+k]>0)){
                //float penetrance = (haploid_arr[i] || j==k) ?
                //penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]
                //: 0 ;
                float penetrance = penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
                if (penetrance>0){
   //cerr<<"NONZERO LOG PEN FOR hap pair "<<j<<","<<k<<" is "<<logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]<<" haploidarr: "<<haploid_arr[i]<<" penetrance: "<<penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]<<" pen:"<<penetrance<<" freq: "<<frequency_cache[j*g_max_haplotypes+k]<<" likelihood: "<<likelihood<<endl;
                  //float freq = frequency_cache[j*g_max_haplotypes+k];
                  float freq = 1;
                  float p = freq*penetrance;
                  likelihood+=(2*p);
                  g_current_weight[j]+=p;
                  g_current_weight[k]+=p;
                }
              }
            } 
          }
        }
      }
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] && (gi_iteration==0 || 
        subject_haplotype_weight[i*g_max_haplotypes+j]>0)){
            subject_haplotype_weight[i*g_max_haplotypes+j] = 
            g_current_weight[j];
        }
      }
      char zero='0';
      int offset=(int)zero;
      if (likelihood>0){
        for(int j=0;j<g_max_haplotypes;++j){
          if(g_active_haplotype[j]){
            g_current_weight[j]/=likelihood;   
            g_weight[j]+=g_current_weight[j];
            if (debug_personhap ) {
             cerr<<"CPU person "<<i<<" hap "<<j<<" weight "<<g_current_weight[j]<<" marginal "<<subject_haplotype_weight[i*g_max_haplotypes+j]<<endl;
            }
          }
        }
      }else{
        //cerr<<"Likelihood for person "<<i<<" is zero\n";
        //exit(1);
      }
    }
    cerr<<"Elapsed time for CPU compute haplotype weights: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    if (debug_weights){
      for(int j=0;j<g_max_haplotypes;++j){
        if(g_active_haplotype[j]){
          cout<<"CPU hap "<<j<<" ";
          for(int m2=0;m2<g_markers;++m2){
            cout<<g_haplotype[j*g_max_window+m2];
          } 
          cout<<" weight "<<g_weight[j]<<endl;
        }
      }
    }
  }//END CPU VERSION
  //cerr<<"done do_iteration\n";
  //if (debug_personhap||debug_weights) exit(1);
}

void MendelGPU::assign_penetrance_element(int subject,int hap1, int hap2,float logpenval,float * pen_cache, bool swap){
  float penval = (logpenval>=gf_logpen_threshold)?exp(logpenval):0;
  if(penval>0){
    if(!swap || g_genotype_imputation){
      pen_cache[subject*penetrance_matrix_size+hap1*
      g_max_haplotypes+hap2] = penval;
      pen_cache[subject*penetrance_matrix_size+hap2*
      g_max_haplotypes+hap1] = penval;
      //order_hap_t oh(logpenetrance_cache[subject*penetrance_matrix_size+hap1* g_max_haplotypes+hap2],penval,hap1,hap2);
      //ohs_unphased.insert(oh);
    }else{
      int no_switch_dist=0;
      int switch_dist=0;
      int first_right_dist=0;
      int second_right_dist=0;
      if (g_flanking_snps){
        int right_half_arr[g_flanking_snps];
        for(int i=0;i<g_flanking_snps;++i) right_half_arr[i] = 0;
        // determine whether a switch occured from left half
        //cerr<<"No switch\n";
        no_switch_dist = hamming(g_left_hap_imputed+(subject*2)*g_flanking_snps,g_haplotype+hap1*g_max_window,g_left_flanking_snps) + hamming(g_left_hap_imputed+(subject*2+1)*g_flanking_snps,g_haplotype+hap2*g_max_window,g_left_flanking_snps);
        //cerr<<"switch\n";
        switch_dist = hamming(g_left_hap_imputed+(subject*2)*g_flanking_snps,g_haplotype+hap2*g_max_window,g_left_flanking_snps) + hamming(g_left_hap_imputed+(subject*2+1)*g_flanking_snps,g_haplotype+hap1*g_max_window,g_left_flanking_snps);
        //cerr<<"first right\n";
        first_right_dist = hamming(right_half_arr,g_haplotype+hap1*g_max_window,g_flanking_snps);
        //cerr<<"second right\n";
        second_right_dist = hamming(right_half_arr,g_haplotype+hap2*g_max_window,g_flanking_snps);
      }
      //cerr<<"SUBJECT "<<subject<<" HAP1,2 "<<hap1<<","<<hap2<<": DIST NO_SWITCH "<<no_switch_dist<<" SWITCH "<<switch_dist<<endl;
      if (no_switch_dist<switch_dist){
        pen_cache[subject*penetrance_matrix_size+hap1* g_max_haplotypes+hap2] = penval;
        pen_cache[subject*penetrance_matrix_size+hap2* g_max_haplotypes+hap1] = 0;
        //cerr<<"Not switching\n";
        //order_hap_t oh(logpenetrance_cache[subject*penetrance_matrix_size+hap1* g_max_haplotypes+hap2],penval,hap1,hap2);
        //ohs_unphased.insert(oh);
      }else if (no_switch_dist>switch_dist){
        pen_cache[subject*penetrance_matrix_size+hap2* g_max_haplotypes+hap1] = penval;
        pen_cache[subject*penetrance_matrix_size+hap1* g_max_haplotypes+hap2] = 0;
        //cerr<<"Switching\n";
        //order_hap_t oh(logpenetrance_cache[subject*penetrance_matrix_size+hap1* g_max_haplotypes+hap2],penval,hap2,hap1);
        //ohs_unphased.insert(oh);
      }else if (first_right_dist<second_right_dist){
        //cerr<<"Not switch from right half\n";
        pen_cache[subject*penetrance_matrix_size+hap1* g_max_haplotypes+hap2] = penval;
        pen_cache[subject*penetrance_matrix_size+hap2* g_max_haplotypes+hap1] = 0;
      }else{
        //cerr<<"Will switch from right half\n";
//cerr<<"PENETRANCE_MATRIX_SIZE: "<<penetrance_matrix_size<<" MAX_HAP: "<<g_max_haplotypes<<" HAP1,2: "<<hap1<<","<<hap2<<" INDEX: "<<subject*penetrance_matrix_size+hap2* g_max_haplotypes+hap1<<" VAL "<<penval<<endl;
        pen_cache[subject*penetrance_matrix_size+hap2* g_max_haplotypes+hap1] = penval;
        pen_cache[subject*penetrance_matrix_size+hap1* g_max_haplotypes+hap2] = 0;
      }
    } 
  }else{
    pen_cache[subject*penetrance_matrix_size+hap1*
    g_max_haplotypes+hap2] = 
    pen_cache[subject*penetrance_matrix_size+hap2*
    g_max_haplotypes+hap1] = penval;
  }
}
