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
  bool debug_personhap = false;
  //bool debug_personhap = gi_iteration==2;
  //bool debug_personhap = g_haplotypes<-4;
  int debug_person = 10;
  bool debug_weights = g_haplotypes<-4;
  bool debug_freq = g_haplotypes<-4;
  if (run_gpu){
    double start = clock();
    #ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_iteration, CL_TRUE, 0,sizeof(int), &gi_iteration, NULL, NULL );
    clSafe(err, "write iteration");
    err = commandQueue->enqueueNDRangeKernel(*kernel_compute_weights,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
    clSafe(err,"launch compute_weights");
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
  //cerr<<"done do_iteration\n";
  if (debug_personhap||debug_weights) exit(1);
  if (debug_freq) exit(1);
}
