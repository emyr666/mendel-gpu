#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif


void GuidedMendelGPU::precompute_penetrance_opencl(){
  bool debug_penmat = g_left_marker<0;
  cerr<<"Entering precompute_penetrance at left marker "<<g_left_marker<<"\n";
  int debug_subject = 1;
  if (run_gpu){
#ifdef USE_GPU
    writeToBuffer(buffer_haplotypes, 1,&g_haplotypes,"buffer_haplotypes");
    writeToBuffer(buffer_packedhap, g_max_haplotypes*packedhap_len, packedhap, "buffer_packedhap");
    writeToBuffer(buffer_extended_snp_mapping, g_max_window,extended_snp_mapping,"buffer_extended_snp_mapping");
    runKernel("kernel_precompute_penetrance",kernel_precompute_penetrance,g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
    if(debug_penmat){
      if (geno_dim==PHASED_INPUT){
        float * debug_cache = new float[g_people*g_max_haplotypes*2];
        readFromBuffer(buffer_logpenetrance_cache, g_people*2*g_max_haplotypes,debug_cache,"buffer_logpenetrance_cache");
        cout<<"GPU:\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<debug_cache[debug_subject*2*g_max_haplotypes+2*j+k];
            }
            cout<<endl;
          }
        }
      }else{
        float * debug_cache = new float[g_people*penetrance_matrix_size];
        readFromBuffer(buffer_logpenetrance_cache, g_people*penetrance_matrix_size,debug_cache,"buffer_logpenetrance_cache");
        cout<<"GPU:\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<debug_cache[debug_subject*penetrance_matrix_size+j*g_max_haplotypes+k];
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

