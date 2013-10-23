using namespace std;
#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif


void MendelGPU::compute_haplotype_weights_opencl(){
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
  if (run_gpu){
    double start = clock();
    #ifdef USE_GPU
    cerr<<"Iteration is "<<gi_iteration<<endl;
    writeToBuffer(buffer_iteration, 1, &gi_iteration, "buffer_iteration" );
    writeToBuffer(buffer_frequency,g_max_haplotypes,g_frequency,"buffer_frequency");
    runKernel("kernel_compute_weights",kernel_compute_weights,g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
    if (debug_freq){
      int active[g_max_haplotypes];
      readFromBuffer(buffer_active_haplotype,g_max_haplotypes,active,"buffer_active_haplotype");
      float freq[g_max_haplotypes];
      readFromBuffer(buffer_frequency,g_max_haplotypes,freq,"buffer_frequency");
      for(int j=0;j<g_max_haplotypes;++j){
        if (active[j]){
          if (debug_freq) cout<<"Freq "<<j<<":"<<freq[j]<<endl;
        }
      }
    }
    if (debug_personhap){
    //runKernel("kernel_compute_weights",kernel_compute_weights,g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
      float subject_haplotype_weight[g_people*g_max_haplotypes];
      readFromBuffer(buffer_subject_haplotype_weight, g_people*g_max_haplotypes,subject_haplotype_weight,"buffer_subject_haplotype_weight");
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
    runKernel("kernel_reduce_weights2",kernel_reduce_weights2,BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
    cerr<<"Launched reduce_weights\n";
    readFromBuffer(buffer_haplotype_weight, g_max_haplotypes,g_weight,"buffer_haplotype_weight");
    if (debug_weights){
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j]){
          cout<<"GPU hap "<<j<<" weight "<<g_weight[j]<<endl;
        }
      }
    }
    cerr<<"Elapsed time for GPU compute haplotype weights: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    #endif
  }
  //cerr<<"done do_iteration\n";
}
