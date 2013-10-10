// this set of kernels compiles successfully at this point
// we need to figure out which data structures to pass in
// from the algorithm

//#include "cl_constants.h"
//

#ifdef phased
#endif

__kernel void simple(
__constant int * scaler,
__global int * outputvec
) {
  int threadindex = get_local_id(0);
  outputvec[threadindex] = *scaler * threadindex;
  return;
}

__kernel void compute_weights(
__const int max_haplotypes,
__const int max_window,
__const int penetrance_matrix_size,
__constant int * iteration,
__constant int * genotype_imputation,
__constant int * markers,
__constant int * haplotypes,
__constant int * left_marker,
#ifdef unphased
__global int * haploid_arr,
#endif
__global int * haplotype,
__global float * frequency,
__global float * penetrance_cache,
__global float * subject_haplotype_weight,
__global int * active_haplotype,
__local int * local_active_haplotype,
__local float * local_cached_marginals,
__local float * local_frequency,
__local float * local_marginals,
__local float * local_penetrance_block,
__local float * local_likelihood
) {
  int subject = get_group_id(0);
#ifdef unphased
  int haploid = haploid_arr[subject];
#endif
  int threadindex = get_local_id(0);

  // INITIALIZE HAPLOTYPE SPECIFIC VARIABLES
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){
      //local_frequency[hapindex] = .1;
      local_frequency[hapindex] = frequency[hapindex];
      local_marginals[hapindex] = 0;
      local_cached_marginals[hapindex] = iteration[0]==0?0:
      subject_haplotype_weight[subject*max_haplotypes+hapindex];
    }
  }
  float likelihood = 0;
  if(threadindex==0) local_likelihood[0] = 0;
  barrier(CLK_LOCAL_MEM_FENCE);
#ifdef unphased
  for(int righthapindex=0;righthapindex<max_haplotypes;++righthapindex){
    if (local_active_haplotype[righthapindex] && 
    (iteration[0]==0 || local_cached_marginals[righthapindex]>0)){
#endif
      local_penetrance_block[threadindex] = 0;
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
#ifdef unphased
        int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
        if (lefthapindex<=righthapindex&&local_active_haplotype[lefthapindex]
        &&(iteration[0]==0||local_cached_marginals[lefthapindex]>0)){
          float penetrance = (haploid || lefthapindex==righthapindex) ?
          penetrance_cache[subject*penetrance_matrix_size+
          righthapindex*max_haplotypes+lefthapindex]: 0;
          if (penetrance>0){
            float prob = (lefthapindex!=righthapindex)?
            local_frequency[lefthapindex]*
            local_frequency[righthapindex]*penetrance*2:
            local_frequency[lefthapindex]*
            local_frequency[righthapindex]*penetrance;
            local_marginals[lefthapindex]+=prob;
            local_penetrance_block[threadindex]+=prob;
          }
        }
#endif
#ifdef phased
    int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
    if (lefthapindex<max_haplotypes&&local_active_haplotype[lefthapindex]
    &&(iteration[0]==0||local_cached_marginals[lefthapindex]>0)){
      for(int gamete=0;gamete<2;++gamete){
        float penetrance = penetrance_cache[subject*max_haplotypes*2+
        lefthapindex*2+gamete];
        if (penetrance>0){
          float prob = local_frequency[lefthapindex]*penetrance;
          local_penetrance_block[threadindex]+=prob;
          local_marginals[lefthapindex] += prob;
        }
      }
    }
#endif
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      // LOG2 REDUCTION
      for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
        if (threadindex<s) {
          local_penetrance_block[threadindex]+=local_penetrance_block[threadindex+s];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if(threadindex==0){
#ifdef unphased
        local_marginals[righthapindex]+=local_penetrance_block[threadindex];
#endif
        local_likelihood[threadindex]+=local_penetrance_block[threadindex];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
#ifdef unphased
    }
  }
#endif
  for(int chunk=0;chunk<max_haplotypes/BLOCK_WIDTH+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      subject_haplotype_weight[subject*max_haplotypes+hapindex] = 
      //frequency[hapindex];
      //local_likelihood[0];
      local_likelihood[0]>0?local_marginals[hapindex] / local_likelihood[0] : 0;
      //local_marginals[hapindex] ;
    }
  }
  return;
}

__kernel void reduce_weights2(
__const int people,
__const int max_haplotypes,
__constant int * haplotypes,
__global float * subject_haplotype_weight,
__global float * haplotype_weight,
__global int * active_haplotype,
__local int * local_active_haplotype,
__local float * local_haplotype_weight
) {
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<max_haplotypes/BLOCK_WIDTH+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_haplotype_weight[hapindex] = 0;
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int subject=0;subject<people;++subject){
    for(int chunk=0;chunk<max_haplotypes/BLOCK_WIDTH+1;++chunk){
      int hapindex = chunk*BLOCK_WIDTH+threadindex;
      if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){
        local_haplotype_weight[hapindex] += 
        subject_haplotype_weight[subject*max_haplotypes+hapindex];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  // store in global memory
  for(int chunk=0;chunk<max_haplotypes/BLOCK_WIDTH+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){
      //haplotype_weight[hapindex] = hapindex;
      haplotype_weight[hapindex] = local_haplotype_weight[hapindex];
    }
  }
  return;
}
