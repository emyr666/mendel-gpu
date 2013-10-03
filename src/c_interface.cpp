#include<iostream>
#include<fstream>
#include<sstream>
#include<list>
#include<set>
#include<queue>
#include<tr1/unordered_set>
#include<tr1/unordered_map>

using namespace std;

#include<map>
#ifdef USE_GPU
#include<CL/cl.hpp>
#include"clsafe.h"
#endif
#include"cl_constants.h"
#include"config.hpp"
#include"io_manager.hpp"
#include"mendel_gpu.hpp"


MendelGPU * mendelGPU = NULL;
  
extern "C" void load_constants_(int * model, int * people,int * snps, int * total_regions,int * haplotype_mode, int * flanking_snps,int * max_haplotypes,int * max_extended_haplotypes,int * platform_id, int * device_id, float * delta, float * lambda, int * i_geno_dim){
  mendelGPU = new MendelGPU();
  mendelGPU->load_constants_(model, people,snps, total_regions,haplotype_mode, flanking_snps,max_haplotypes,max_extended_haplotypes,platform_id, device_id, delta, lambda, i_geno_dim);
}

extern "C" void init_gpu_(int * platform_id, int * device_id){
  mendelGPU->init_gpu_(platform_id, device_id);
}

extern "C" void run_simple_(int * scaler, int * return_vec){
  mendelGPU->run_simple_(scaler, return_vec);
}

//extern "C" void read_stream_(int * people, int * snps){
 // mendelGPU->read_stream_(people, snps);
//}

extern "C" void init_buffers_(int * active_haplotype,int * max_window,int * max_haplotypes,int * max_extended_haplotypes,int * max_region_size,int * people, int * snps, int * genotype_imputation,int * haplotypes, int * markers, int * window1,int * prev_left_marker, int * window2, int * prev_right_marker,float * frequency, float * weight, int * haplotype,int * flanking_snps){
  mendelGPU->init_buffers_(active_haplotype,max_window,max_haplotypes,max_extended_haplotypes,max_region_size, people, snps, genotype_imputation,haplotypes, markers, window1,prev_left_marker, window2,prev_right_marker,frequency, weight, haplotype, flanking_snps);
}

extern "C" void init_region_buffers_(int * regionstart,int * regionend){
  mendelGPU->init_region_buffers_(regionstart,regionend);
}

extern "C" void double_haplotypes_(){
  mendelGPU->double_haplotypes_();
}

extern "C" void prune_haplotypes_(){
  mendelGPU->prune_haplotypes_();
}

extern "C" void copy_ref_haplotypes_(){
  mendelGPU->copy_ref_haplotypes_();
}

extern "C" void init_window_buffers_(){
  mendelGPU->init_window_buffers_();
}

extern "C" void init_iteration_buffers_(){
  mendelGPU->init_iteration_buffers_();
}

extern "C" void precompute_penetrance_fast_(int * center_snp){
  mendelGPU->precompute_penetrance_fast_( center_snp);
}

extern "C" void impute_penetrance_matrix_(){
  mendelGPU->impute_penetrance_matrix_();
}

extern "C" void precompute_penetrance_(){
  mendelGPU->precompute_penetrance_();
}

extern "C" void compute_haplotype_weights_(int * iteration){
  mendelGPU->compute_haplotype_weights_(iteration);
}

extern "C" void prep_impute_genotypes_guide_(){
  mendelGPU->prep_impute_genotypes_guide_();
}

extern "C" void impute_diploid_genotypes_guide_(int * center_snp_start, int * center_snp_end, int * center_snp_offset){
  mendelGPU->impute_diploid_genotypes_guide_(center_snp_start, center_snp_end, center_snp_offset);
}

extern "C" void impute_haploid_genotypes_denovo_(int * center_snp, int * d_snp){
  mendelGPU->impute_haploid_genotypes_denovo_(center_snp, d_snp);
}

extern "C" void impute_haploid_genotypes_guide_(int * center_snp_start, int * center_snp_end, int * center_snp_offset){
  mendelGPU->impute_haploid_genotypes_guide_(center_snp_start, center_snp_end, center_snp_offset);
}

extern "C" void impute_diploid_genotypes_denovo_(int * center_snp, int * d_snp){
  mendelGPU->impute_diploid_genotypes_denovo_(center_snp, d_snp);
}

extern "C" void cleanup_(){
  delete mendelGPU;
}

