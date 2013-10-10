#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"


bool MendelGPU::check_mm_converged(){
  bool debug_freq = false;
  cerr<<"Doing MM updates\n";
  float d = 0;
  float e = 0;
  float t = 0;
  for(int j=0;j<g_max_haplotypes;++j){
    if (g_active_haplotype[j]){
      //cerr<<"Current weight of hap "<<j<<" is "<<g_weight[j]<<endl;
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
  cerr<<"MM: T is "<<t<<endl;
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
  bool debug_personhap = g_left_marker<0;
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
                float penetrance = (haploid_arr[i] || j==k) ?
                penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]
                : 0 ;
                if (penetrance>0){
                  float freq = frequency_cache[j*g_max_haplotypes+k];
                  float p = freq*penetrance;
                  likelihood+=p;
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
             cout<<"CPU person "<<i<<" hap "<<j<<" weight "<<g_current_weight[j]<<" marginal "<<subject_haplotype_weight[i*g_max_haplotypes+j]<<endl;
            }
          }
        }
      }else{
        //cerr<<"Likelihood for person "<<i<<" is zero\n";
        //exit(1);
      }
    }
    cerr<<"Compute haplotype weights: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
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
  if (debug_personhap||debug_weights) exit(1);
}
