
using namespace std;
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<queue>
#include<tr1/unordered_set>
#include<tr1/unordered_map>
#include<set>
#include<map>
#include<list>
#include<math.h>
#include"../cl_constants.h"
#include<cstdlib>
#include<string.h>
#ifdef USE_GPU
#include<CL/cl.hpp>
#include"../clsafe.h"
#endif
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"


bool MendelGPU::check_mm_converged(){
  cerr<<"Doing MM updates\n";
  float d = 0;
  float e = 0;
  float t = 0;
  for(int j=0;j<g_max_haplotypes;++j){
    if (g_active_haplotype[j]){
      //cerr<<"Current weight of hap "<<j<<" is "<<g_weight[j]<<endl;
      if (g_weight[j]<epsilon) g_weight[j] = epsilon;
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
        //cerr<<"1 Assigning new freq "<<j<<" of "<<g_frequency[j]<<endl;
        if(g_frequency[j]<epsilon) g_frequency[j] = epsilon;
        
      }
    }
  }else if(t>0){
    for(int j=0;j<g_max_haplotypes;++j){
      g_frequency[j] = g_weight[j]/t;
        //cerr<<"2 Assigning new freq "<<j<<" of "<<g_frequency[j]<<endl;
    }
  }else{
    for(int j=0;j<g_max_haplotypes;++j){
      g_frequency[j] = 1./g_max_haplotypes;
        //cerr<<"3 Assigning new freq "<<j<<" of "<<g_frequency[j]<<endl;
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
  if (t<convergence_criterion) return true;
  else{
    for(int j=0;j<g_max_haplotypes;++j){
      g_old_frequency[j] = g_frequency[j];
    }
  }
  return false;
}

void MendelGPU::compute_haplotype_weights_(int * iteration){
  bool b_impute = g_genotype_imputation;
  cerr<<"begin compute_weights at iteration: "<<iteration[0]<<"\n";
  //cerr<<"Geno dim is "<<geno_dim<<endl;
  //return;
  bool debug_personhap = false;
  //bool debug_personhap = iteration[0]==2;
  //bool debug_personhap = g_haplotypes<-4;
  int debug_person = 10;
  bool debug_weights = g_haplotypes<-4;
  bool debug_freq = g_haplotypes<-4;
  if (run_gpu){
    double start = clock();
    #ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_iteration, CL_TRUE, 0,sizeof(int), iteration, NULL, NULL );
    clSafe(err, "write iteration");
    if (geno_dim==PHASED_INPUT){
      // for phased input
      err = commandQueue->enqueueNDRangeKernel(*kernel_compute_weights_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch compute_weights");
    }else if (geno_dim==UNPHASED_INPUT){
      // for unphased input
      err = commandQueue->enqueueNDRangeKernel(*kernel_compute_weights,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch compute_weights");
    }
    cerr<<"Launched compute_haplotype_weights\n";
    if (debug_personhap){
      float subject_haplotype_weight[g_people*g_max_haplotypes];
      err = commandQueue->enqueueReadBuffer(*buffer_subject_haplotype_weight, CL_TRUE, 0, sizeof(float)*g_people*g_max_haplotypes,subject_haplotype_weight);
      clSafe(err, "read subject_hap weight");
      float haplotype_weight[g_max_haplotypes];
      memset(haplotype_weight,0,sizeof(float)*g_max_haplotypes);
      for(int i = 0;i<g_people;++i){
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<"GPU person "<<i<<" hap "<<j<<" weight "<<subject_haplotype_weight[i*g_max_haplotypes+j]<<endl;
            haplotype_weight[j]+=subject_haplotype_weight[i*g_max_haplotypes+j];
          }
        }
      }
    }
    err = commandQueue->enqueueNDRangeKernel(*kernel_reduce_weights2,cl::NullRange,cl::NDRange(BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
    clSafe(err,"launch reduce_weights");
    cerr<<"Launched reduce_weights\n";
    err = commandQueue->enqueueReadBuffer(*buffer_haplotype_weight, CL_TRUE, 0, sizeof(float)*g_max_haplotypes,g_weight);
    clSafe(err, "read hap weight");
    if (debug_weights){
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j]){
          cout<<"GPU hap "<<j<<" weight "<<g_weight[j]<<endl;
        }
      }
    }
    cerr<<"Elapsed time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    #endif
  }
  if(run_cpu){
    memset(g_weight,0,sizeof(float)*g_max_haplotypes);
    if (iteration[0]==0){
      memset(subject_haplotype_weight,0,sizeof(float)*g_people*g_max_haplotypes);
    }
    //for(int i=0;i<0;++i){
    double start = clock();
    if (debug_freq) cout<<"CPU freq: "<<endl;
    if (geno_dim==UNPHASED_INPUT){
      for(int j=0;j<g_max_haplotypes;++j){
        int start = b_impute?j:0;
        if (g_active_haplotype[j]){
          for(int k=start;k<g_max_haplotypes;++k){
            if (g_active_haplotype[k]){
              frequency_cache[j*g_max_haplotypes+k] = g_frequency[j]*g_frequency[k];
              if (j!=k) frequency_cache[j*g_max_haplotypes+k]*=2;
              if (debug_freq) cout<<" "<<frequency_cache[j*g_max_haplotypes+k];
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
        if (g_active_haplotype[j] && (iteration[0]==0 || subject_haplotype_weight[i*g_max_haplotypes+j]>0)){
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
              (iteration[0]==0 || 
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
                //cerr<<"freq is "<<freq<<" and pen is "<<penetrance<<endl;
                //cerr<<"weight for haps "<<j<<","<<k<<" is now "<<g_current_weight[j]<<endl;
                }
              }
            } 
          }
        }
      }
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] && (iteration[0]==0 || 
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
  if (debug_freq) exit(1);
}
