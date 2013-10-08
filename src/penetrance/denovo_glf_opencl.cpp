#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void DenovoGlfMendelGPU::precompute_penetrance_fast_opencl(){
  bool debug_penmat = g_markers==-6;
  //bool debug_penmat = g_left_marker==-2;
  int debug_penmat_person = 3;
  cerr<<endl;
  if (run_gpu){
#ifdef USE_GPU
    writeToBuffer(buffer_right_edge_dosage,g_max_haplotypes,right_edge_dosage,"buffer_right_edge_dosage");
    if (g_prev_left_marker!=g_left_marker){
      writeToBuffer(buffer_beyond_left_edge_dosage,g_max_haplotypes,beyond_left_edge_dosage,"buffer_beyond_left_edge_dosage");
    }
    runKernel("precompute_penetrance_fast",kernel_precompute_penetrance_fast,g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
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
        readFromBuffer(buffer_penetrance_cache, g_people*penetrance_matrix_size,debug_cache,"buffer_penetrance_cache");
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
}

void DenovoGlfMendelGPU::impute_penetrance_matrix_opencl(){
  bool debug_penmat = g_left_marker==-3;
  int debug_penmat_person = 10;
  cerr<<"Entering impute penetrance matrix\n";
  //return;
  if (run_gpu){
#ifdef USE_GPU
    if (debug_penmat){
      float * debug_cache = new float[g_people*penetrance_matrix_size];
      for(int i=0;i<g_people;++i){
        for(int j=0;j<penetrance_matrix_size;++j){
          debug_cache[i*penetrance_matrix_size+j] = j;
        }
      }
      //twin_hap_index[0] =2;
      //twin_hap_index[1] =3;
      //writeToBuffer(buffer_logpenetrance_cache,g_people*penetrance_matrix_size,debug_cache,"buffer_logpenetrance_cache");
      delete[] debug_cache;
    }
    writeToBuffer(buffer_twin_hap_index,g_max_haplotypes,twin_hap_index,"buffer_twin_hap_index");
    runKernel("impute_penetrance",kernel_impute_penetrance,g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
    if (debug_penmat){
      if (geno_dim==PHASED_INPUT){
        float * debug_cache = new float[g_people*2*g_max_haplotypes];
        readFromBuffer(buffer_logpenetrance_cache,g_people*2*g_max_haplotypes,debug_cache,"buffer_logpenetrance_cache");
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
cerr<<"g_max_haplotypes "<<g_max_haplotypes<<" and people "<<g_people<<endl;
        readFromBuffer(buffer_logpenetrance_cache,g_people*penetrance_matrix_size,debug_cache,"buffer_logpenetrance_cache");
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
}

