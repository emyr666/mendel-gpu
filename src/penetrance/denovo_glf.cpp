#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../clsafe.h"
#endif

void DenovoGlfMendelGPU::compute_penetrance(){
  impute_penetrance_matrix();
  precompute_penetrance_fast(g_center_snp_start);
}

void DenovoGlfMendelGPU::precompute_penetrance_fast(int csnp){
  // if we are in sequencing mode, assume non polymorphic sites are ref allele 
  if (g_likelihood_mode == Config::LIKELIHOOD_MODE_READS && !polymorphic_window){
    cerr<<"Skipping non polymorphic site\n";
    return;
  }
  cerr<<"Entering precompute_penetrance_fast on center_snp "<<csnp<<"\n";
  int last_marker = g_markers-1;
  int prev_left_marker = g_prev_left_marker;
  int left_marker = g_left_marker;
  cerr<<"Taking hap doses at site "<<last_marker<<" and prev and current left marker: "<<prev_left_marker<<","<<left_marker<<"\n";
  int debug_penmat_person = 250;
  //bool debug_penmat = false;
  bool debug_penmat = g_markers<-8;
  //bool debug_penmat = g_haplotypes>2;
  if (debug_penmat) cout<<"Debugging precompute_penetrance_fast_()\n";
  cerr<<"right edge dosage: ";
  for(int j=0;j<g_max_haplotypes;++j){
    if(g_active_haplotype[j] ){
      right_edge_dosage[j] = g_haplotype[j*g_max_window+last_marker]; 
      cerr<<" "<<right_edge_dosage[j];
    }
  }
  cerr<<endl;
  if (run_gpu){
#ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_right_edge_dosage, CL_TRUE, 0, sizeof(int)*g_max_haplotypes,right_edge_dosage,NULL,NULL);
    clSafe(err, "write right edge dosage");
    if (prev_left_marker!=left_marker){
      err = commandQueue->enqueueWriteBuffer(*buffer_beyond_left_edge_dosage, CL_TRUE, 0, sizeof(int)*g_max_haplotypes,beyond_left_edge_dosage,NULL,NULL);
      clSafe(err, "write beyond left edge dosage");
    }
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotypes, CL_TRUE, 0, sizeof(int),&g_haplotypes,NULL,NULL);
    clSafe(err, "write total haplotypes");
    if (geno_dim==PHASED_INPUT){
      cerr<<"Launching haploid fast precompute\n";
      err = commandQueue->enqueueNDRangeKernel(*kernel_precompute_penetrance_fast_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch precompute penetrance_fast");
    }else if(geno_dim==UNPHASED_INPUT){
      err = commandQueue->enqueueNDRangeKernel(*kernel_precompute_penetrance_fast,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch precompute penetrance_fast");
    }
    if(debug_penmat){
      if (geno_dim==PHASED_INPUT){
        float * debug_cache = new float[g_people*2*g_max_haplotypes];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*2*g_max_haplotypes,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"fast GPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<debug_cache[debug_penmat_person*g_max_haplotypes*2+j*2+k];
            }
            cout<<endl;
          }
        }
        delete[] debug_cache;
      }else if (geno_dim==UNPHASED_INPUT){ 
        float * debug_cache = new float[g_people*penetrance_matrix_size];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*penetrance_matrix_size,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"fast GPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<debug_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
        delete[] debug_cache;
      }
    }
#endif
  }
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
              g_snps)+ geno_dim*(left_marker+last_marker) + 2*parent+dose];
      //        cerr<<"Person "<<i<<" parent "<<parent<<" SNP "<<left_marker+last_marker<<" dose "<<dose<<" value "<<
       //       logpenetrance_new[parent]<<endl;
            }
            if (prev_left_marker!=left_marker){
              dose=beyond_left_edge_dosage[j];
              for(int parent=0;parent<2;++parent){
                logpenetrance_old[parent] = g_snp_penetrance[i*(geno_dim*
                g_snps)+ geno_dim*(prev_left_marker) + 2*parent+dose];
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
                  g_snps)+ geno_dim*(left_marker+last_marker) + dose];
                  if (debug_penmat && i==debug_penmat_person){
                    cerr<<"logpenetrance new for person "<<i<<" haps "<<j<<","<<k<<": dose "<<dose<<":"<<logpenetrance_new<<endl;
                  }
                  float logpenetrance_old = 0;
                  if (prev_left_marker!=left_marker){
                    dose=beyond_left_edge_dosage[j] + beyond_left_edge_dosage[k];
                    logpenetrance_old = g_snp_penetrance[i*(geno_dim*
                    g_snps)+ geno_dim*(prev_left_marker) + dose];
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
              val>=logpenetrance_threshold?exp(val):0;
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
                g_max_haplotypes+k] = val>=logpenetrance_threshold?exp(val):0;
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
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<"haps "<<j<<","<<k<<": "<<logpenetrance_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k]<<endl;
              }
            }
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
  int left_marker=g_left_marker;
  bool debug_basic = false;
  bool debug_penmat = false;
  int debug_penmat_person = 600;
  cerr<<"Entering impute penetrance matrix\n";
  //return;
  if (run_gpu){
#ifdef USE_GPU
    if (debug_basic){
      cerr<<"Twin indices:";
      for(int i=0;i<g_max_haplotypes;++i){
        g_active_haplotype[i] = (i%2==0);
        twin_hap_index[i] = i<g_max_haplotypes/2 ? i+g_max_haplotypes/2 : i;
        cerr<<" "<<i<<":"<<twin_hap_index[i];
      }
      cerr<<endl;
      for(int i=0;i<g_people;++i){
        for(int j=0;j<g_max_haplotypes;++j){
          for(int k=0;k<g_max_haplotypes;++k){
            logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = (j>=g_max_haplotypes/2 && k>=g_max_haplotypes/2)? j*g_max_haplotypes+k : 0;      
          }
        }
      }
      err = commandQueue->enqueueWriteBuffer(*buffer_active_haplotype, CL_TRUE, 0, sizeof(int)*g_max_haplotypes,g_active_haplotype,NULL,NULL);
      clSafe(err, "write active haplotype");
      err = commandQueue->enqueueWriteBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people * penetrance_matrix_size ,logpenetrance_cache,NULL,NULL);
      clSafe(err, "write logpenetrance cache");
    }
    err = commandQueue->enqueueWriteBuffer(*buffer_twin_hap_index, CL_TRUE, 0, sizeof(int)*g_max_haplotypes,twin_hap_index,NULL,NULL);
    clSafe(err, "write twin indices");
    if (geno_dim==PHASED_INPUT){
      err = commandQueue->enqueueNDRangeKernel(*kernel_impute_penetrance_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch impute penetrance_haploid");
    }else if(geno_dim==UNPHASED_INPUT){
      err = commandQueue->enqueueNDRangeKernel(*kernel_impute_penetrance,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch impute penetrance");
    }
    if (debug_penmat){
      if (geno_dim==PHASED_INPUT){
        float * debug_cache = new float[g_people*2*g_max_haplotypes];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*2*g_max_haplotypes,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"GPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<debug_cache[debug_penmat_person*g_max_haplotypes*2+j*2+k];
            }
            cout<<endl;
          }
        }
        delete[] debug_cache;
      }else if (geno_dim==UNPHASED_INPUT){
        float * debug_cache = new float[g_people*penetrance_matrix_size];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*penetrance_matrix_size,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"GPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<debug_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
        delete[] debug_cache;
      }
    }
#endif
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

