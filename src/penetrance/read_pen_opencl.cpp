#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void ReadPenetrance::process_read_matrices_opencl(){
  
  int penetrance_matrix_size = mendelgpu->g_max_haplotypes*mendelgpu->g_max_haplotypes;
  //bool debug_penmat = g_markers==-6;
  bool debug_penmat = false;
  int debug_penmat_person = 11;
  if (mendelgpu->run_gpu){
#ifdef USE_GPU
  //  double start = clock();
    mendelgpu->writeToBuffer(buffer_read_alleles_vec,mendelgpu->g_people*max_vector_length,read_alleles_vec,"buffer_read_alleles_vec");
    mendelgpu->writeToBuffer(buffer_read_match_logvec,mendelgpu->g_people*max_vector_length,read_match_logvec,"buffer_read_match_logvec");
    mendelgpu->writeToBuffer(buffer_read_mismatch_logvec,mendelgpu->g_people*max_vector_length,read_mismatch_logvec,"buffer_read_mismatch_logvec");
    mendelgpu->writeToBuffer(buffer_vector_offsets,mendelgpu->g_people*compact_rows,vector_offsets,"buffer_vector_offsets");
    mendelgpu->writeToBuffer(buffer_read_lengths,mendelgpu->g_people*compact_rows,read_lengths,"buffer_read_lengths");
    mendelgpu->writeToBuffer(buffer_haplotype_offsets,mendelgpu->g_people*compact_rows,haplotype_offsets,"buffer_haplotype_offsets");
    mendelgpu->writeToBuffer(buffer_mat_rows_by_subject,mendelgpu->g_people,mat_rows_by_subject,"buffer_mat_rows_by_subject");
    mendelgpu->writeToBuffer(buffer_best_haplotypes,mendelgpu->g_people*mendelgpu->g_max_haplotypes,best_haplotypes,"buffer_best_haplotypes");
    int iteration;
    //double start = clock();
    mendelgpu->runKernel("reads_compute_penetrance",kernel_reads_compute_penetrance,mendelgpu->g_max_haplotypes*SMALL_BLOCK_WIDTH,mendelgpu->g_max_haplotypes,mendelgpu->g_people,SMALL_BLOCK_WIDTH,1,1);
    //mendelgpu->readFromBuffer(mendelgpu->buffer_iteration, 1,&iteration,"mendelgpu->buffer_iteration");
    //cerr<<"Elapsed time for GPU read_compute_pen: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    //start = clock();
    mendelgpu->runKernel("reads_adjust_penetrance",kernel_reads_adjust_penetrance,mendelgpu->g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
    //mendelgpu->readFromBuffer(mendelgpu->buffer_iteration, 1,&iteration,"mendelgpu->buffer_iteration");
    //cerr<<"Elapsed time for GPU adjust_pen: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    //int iteration;
    //mendelgpu->readFromBuffer(mendelgpu->buffer_iteration, 1,&iteration,"mendelgpu->buffer_iteration");
    //cerr<<"Elapsed time for GPU process_read_matrices: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    if(debug_penmat){
      float * debug_cache = new float[mendelgpu->g_people*penetrance_matrix_size];
      float * debug_logcache = new float[mendelgpu->g_people*penetrance_matrix_size];
      mendelgpu->readFromBuffer(mendelgpu->buffer_penetrance_cache, mendelgpu->g_people*penetrance_matrix_size,debug_cache,"mendelgpu->buffer_penetrance_cache");
      mendelgpu->readFromBuffer(mendelgpu->buffer_logpenetrance_cache, mendelgpu->g_people*penetrance_matrix_size,debug_logcache,"mendelgpu->buffer_logpenetrance_cache");
      cerr<<"fast GPU for person "<<debug_penmat_person<<":\n";
      for (int subject=debug_penmat_person;subject<=debug_penmat_person;++subject){
        cerr<<"SUBJECT "<<subject<<endl;
        for(int j=0;j<mendelgpu->g_max_haplotypes;++j){
          if (mendelgpu->g_active_haplotype[j]){
            for(int k=0;k<mendelgpu->g_max_haplotypes;++k){
              if (mendelgpu->g_active_haplotype[k]){
                //cerr<<" "<<debug_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k];
                cerr<<" "<<j<<","<<k<<":"<<debug_logcache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k]<<","<<debug_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k];
              }
            }
            cerr<<endl;
          }
        }
      }
      delete[] debug_logcache;
      delete[] debug_cache;
    }
#endif
  }
}
