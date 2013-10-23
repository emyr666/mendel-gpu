#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"

void DenovoGlfMendelGPU::compute_penetrance(){
  impute_penetrance_matrix();
  precompute_penetrance_fast();
}

void DenovoGlfMendelGPU::precompute_penetrance_fast(){
  cerr<<"Entering precompute_penetrance_fast on center_snp "<<g_center_snp_start<<"\n";
  int last_marker = g_markers-1;
  //cerr<<"Taking hap doses at site "<<last_marker<<" and prev and current left marker: "<<g_prev_left_marker<<","<<g_left_marker<<"\n";
  //cerr<<"right edge dosage: ";
  for(int j=0;j<g_max_haplotypes;++j){
    if(g_active_haplotype[j] ){
      right_edge_dosage[j] = g_haplotype[j*g_max_window+last_marker]; 
      //cerr<<" "<<right_edge_dosage[j];
    } 
  }
  //cerr<<endl;
  if (run_gpu){
#ifdef USE_GPU
    precompute_penetrance_fast_opencl();
#endif
  }
  int debug_penmat_person = 0;
  //bool debug_penmat = g_markers==-6;
  bool debug_penmat = g_left_marker==-91;
  if(run_cpu){
    double start = clock();
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
      float maxlog = -1e10;
      for(int j=0;j<g_max_haplotypes;++j){
        if(g_active_haplotype[j] ){
          if (geno_dim==PHASED_INPUT){
            int dose= right_edge_dosage[j];
            float logpenetrance_new[2];
            float logpenetrance_old[2]={0,0};
            float logpenetrance[2];
            for(int parent=0;parent<2;++parent){
              logpenetrance_new[parent] = g_snp_penetrance[i*(geno_dim*
              g_snps)+ geno_dim*(g_left_marker+last_marker) + 2*parent+dose];
      //        cerr<<"Person "<<i<<" parent "<<parent<<" SNP "<<g_left_marker+last_marker<<" dose "<<dose<<" value "<<
       //       logpenetrance_new[parent]<<endl;
            }
            if (g_prev_left_marker!=g_left_marker){
              dose=beyond_left_edge_dosage[j];
              for(int parent=0;parent<2;++parent){
                logpenetrance_old[parent] = g_snp_penetrance[i*(geno_dim*
                g_snps)+ geno_dim*(g_prev_left_marker) + 2*parent+dose];
              }
            }
            for(int parent=0;parent<2;++parent){
              logpenetrance[parent] = g_haplotypes>2?logpenetrance_cache[i*
              2*g_max_haplotypes+2*j+parent] + logpenetrance_new[parent]-
              logpenetrance_old[parent]: logpenetrance_new[parent];
              if (logpenetrance[parent] > maxlog) 
              maxlog = logpenetrance[parent];
              logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              logpenetrance[parent]; 
            }
          }else if (geno_dim==UNPHASED_INPUT){
            for(int k=j;k<g_max_haplotypes;++k){
              if(g_active_haplotype[k] && (!haploid || k==j)){
                if (g_likelihood_mode==Config::LIKELIHOOD_MODE_GENOTYPES){
                  int dose= right_edge_dosage[j] + right_edge_dosage[k];
                  float logpenetrance_new = g_snp_penetrance[i*(geno_dim*
                  g_snps)+ geno_dim*(g_left_marker+last_marker) + dose];
                  float logpenetrance_old = 0;
                  if (g_prev_left_marker!=g_left_marker){
                    dose=beyond_left_edge_dosage[j] + beyond_left_edge_dosage[k];
                    logpenetrance_old = g_snp_penetrance[i*(geno_dim*
                    g_snps)+ geno_dim*(g_prev_left_marker) + dose];
                  }
                  float logpenetrance = g_haplotypes>2?logpenetrance_cache[i*
                  penetrance_matrix_size+j*g_max_haplotypes+k] + 
                  logpenetrance_new-logpenetrance_old:logpenetrance_new;
                  if (logpenetrance > maxlog) maxlog = logpenetrance;
                  logpenetrance_cache[i*penetrance_matrix_size+j*
                  g_max_haplotypes+k] = logpenetrance; 
                }
              }
            }
          }else{
            cerr<<"Precompute penetrance fast found geno dim of "<<geno_dim<<endl;
          }
        } 
      }
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] ){
          if (geno_dim==PHASED_INPUT){
            for(int parent=0;parent<2;++parent){
              float val = logpenetrance_cache[i*2*g_max_haplotypes+2*
              j+parent]-maxlog;
              logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent] = val;
              penetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              val>=gf_logpen_threshold?exp(val):0;
            }
          }else{
            for(int k=j;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k] && (!haploid || k==j)){
                float val = logpenetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k]-maxlog;
                logpenetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k] = val;
                logpenetrance_cache[i*penetrance_matrix_size+k*
                g_max_haplotypes+j] = val;
                penetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k] = val>=gf_logpen_threshold?exp(val):0;
                penetrance_cache[i*penetrance_matrix_size+k*
                g_max_haplotypes+j] = penetrance_cache[i*
                penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
          }
        }   
      }
    }
    if(debug_penmat){
      if (geno_dim==PHASED_INPUT){
        cout<<"fast CPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<logpenetrance_cache[debug_penmat_person*g_max_haplotypes*2+j*2+k];
            }
            cout<<endl;
          }
        }
      }else if (geno_dim==UNPHASED_INPUT){ 
        cout<<"fast CPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<penetrance_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
                //cout<<" "<<logpenetrance_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
      }
    }
    cerr<<"Precompute penetrance fast time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
  }
  if (debug_penmat) exit(1);
  cerr<<"Exiting precompute_penetrance_fast\n";
}

void DenovoGlfMendelGPU::impute_penetrance_matrix(){
  if (run_gpu){
#ifdef USE_GPU
    impute_penetrance_matrix_opencl();
#endif
  }
  bool debug_penmat = g_left_marker==-3;
  int debug_penmat_person = 10;
  cerr<<"Entering impute penetrance matrix\n";
  if (debug_penmat){
    for(int i=0;i<g_people;++i){
      for(int j=0;j<penetrance_matrix_size;++j){
        //logpenetrance_cache[i*penetrance_matrix_size+j] = j;
      }
    }
    //twin_hap_index[0]=2;
    //twin_hap_index[1]=3;
  }
  if (run_cpu){
    for(int i=0;i<g_people;++i){
      // for diplid case, columns matter
      if (geno_dim==UNPHASED_INPUT){  
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            for(int k=0;k<g_max_haplotypes;++k){ //copy col
              if (g_active_haplotype[k] && twin_hap_index[k]!=k){
                logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+twin_hap_index[k]]; 
              }
            }
          }
        }
      }
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] && twin_hap_index[j]!=j){
          if (geno_dim==PHASED_INPUT){
            // copy penetrance for maternal and paternal haplotypes
            for(int parent=0;parent<2;++parent){
              logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              logpenetrance_cache[i*2*g_max_haplotypes+2*twin_hap_index[j]
              +parent]; 
              //cerr<<"Copied from "<<twin_hap_index[j]<<" to "<<j<<" value "<<logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent]<<endl;
            }
          }else{
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                logpenetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k] = logpenetrance_cache[i*
                penetrance_matrix_size+twin_hap_index[j]*g_max_haplotypes+k]; 
              }
            }
          }
        }
      }
    }
    if(debug_penmat){
      cout<<"CPU for person "<<debug_penmat_person<<":\n";
      if (geno_dim==PHASED_INPUT){
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<logpenetrance_cache[debug_penmat_person*g_max_haplotypes*2+j*2+k];
            }
            cout<<endl;
          }
        }
      }else if (geno_dim==UNPHASED_INPUT){
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<logpenetrance_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
      }
    }
  }
  cerr<<"Leaving impute penetrance matrix\n";
  if (debug_penmat) exit(1);
}

