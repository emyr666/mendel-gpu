#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void DenovoReadsMendelGPU::process_read_matrices_opencl(){
  //bool debug_penmat = g_markers==-6;
  bool debug_penmat = true;
  int debug_penmat_person = 10;
  if (run_gpu){
#ifdef USE_GPU
    writeToBuffer(buffer_read_alleles_mat,g_people*read_compact_matrix_size,read_alleles_mat,"buffer_read_alleles_mat");
    writeToBuffer(buffer_read_match_logmat,g_people*read_compact_matrix_size,read_match_logmat,"buffer_read_match_logmat");
    writeToBuffer(buffer_read_mismatch_logmat,g_people*read_compact_matrix_size,read_mismatch_logmat,"buffer_read_mismatch_logmat");
    writeToBuffer(buffer_mat_rows_by_subject,g_people,mat_rows_by_subject,"buffer_mat_rows_by_subject");
    runKernel("reads_compute_penetrance",kernel_reads_compute_penetrance,g_max_haplotypes*BLOCK_WIDTH,g_max_haplotypes,g_people,BLOCK_WIDTH,1,1);
    runKernel("reads_adjust_penetrance",kernel_reads_adjust_penetrance,g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
    if(debug_penmat){
      float * debug_cache = new float[g_people*penetrance_matrix_size];
      float * debug_logcache = new float[g_people*penetrance_matrix_size];
      readFromBuffer(buffer_penetrance_cache, g_people*penetrance_matrix_size,debug_cache,"buffer_penetrance_cache");
      readFromBuffer(buffer_logpenetrance_cache, g_people*penetrance_matrix_size,debug_logcache,"buffer_logpenetrance_cache");
      cout<<"fast GPU for person "<<debug_penmat_person<<":\n";
      for (int subject=debug_penmat_person;subject<=debug_penmat_person;++subject){
        cout<<"SUBJECT "<<subject<<endl;
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<debug_cache[subject*penetrance_matrix_size+j*g_max_haplotypes+k];
                //cout<<" "<<debug_logcache[subject*penetrance_matrix_size+j*g_max_haplotypes+k]<<","<<debug_cache[subject*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
      }
      delete[] debug_logcache;
      delete[] debug_cache;
    }
#endif
  }
}
