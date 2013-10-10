__kernel void impute_genotype_denovo(
__const int max_haplotypes,
__const int max_window,
__const float epsilon,
__const int penetrance_matrix_size,
__const int geno_dim,
__constant int * genotype_imputation,
__constant int * markers,
__constant int * haplotypes,
__constant int * left_marker,
#ifdef unphased
__global int * haploid_arr,
#endif
__global float * frequency,
__global float * penetrance_cache,
__global int * subject_genotype,
__global float * subject_dosage,
__global float * subject_posterior_prob,
__global int * active_haplotype,
__global int * center_dosage,
__global float * subject_haplotype_weight,
__local int * local_active_haplotype,
__local int * local_center_dosage,
__local float * local_frequency,
__local float * local_posterior_prob0,
__local float * local_posterior_prob1,
__local float * local_posterior_prob2,
__local float * local_posterior_prob3,
__local float * local_posterior_prob,
__local float * local_cached_marginals, 
__local float * local_denom
) {
  int subject = get_group_id(0);
#ifdef unphased
  int haploid = haploid_arr[subject];
#endif
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
  if (threadindex==0) local_denom[threadindex]=0;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){
      local_frequency[hapindex] = frequency[hapindex];
      local_center_dosage[hapindex] = center_dosage[hapindex];
    }
  }
  local_posterior_prob0[threadindex]=0;
  local_posterior_prob1[threadindex]=0;
  local_posterior_prob2[threadindex]=0;
  local_posterior_prob3[threadindex]=0;
  barrier(CLK_LOCAL_MEM_FENCE);
  #ifdef unphased
  for(int righthapindex=0;righthapindex<max_haplotypes;++righthapindex){
    if (local_active_haplotype[righthapindex]){
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
        int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
        if (lefthapindex<=righthapindex&&local_active_haplotype[lefthapindex]
        && (!haploid||lefthapindex==righthapindex)){
          float penetrance = penetrance_cache[subject*penetrance_matrix_size+
          righthapindex*max_haplotypes+lefthapindex];
          if (penetrance>0){
            float prob = (lefthapindex!=righthapindex)?
            local_frequency[lefthapindex]*
            local_frequency[righthapindex]*penetrance*2:
            local_frequency[lefthapindex]*
            local_frequency[righthapindex]*penetrance;


            // COMPUTE POSTERIOR PROBS
            int m = genotype_imputation[0]?local_center_dosage[lefthapindex] +
              local_center_dosage[righthapindex]:local_center_dosage[lefthapindex] +
              2*local_center_dosage[righthapindex];
            local_posterior_prob0[threadindex]+=(m==0)*prob;
            local_posterior_prob1[threadindex]+=(m==1)*prob;
            local_posterior_prob2[threadindex]+=(m==2)*prob;
            local_posterior_prob3[threadindex]+=(m==3)*prob;
          }
        } 
        barrier(CLK_LOCAL_MEM_FENCE);
      }
    }
  }
  #endif
  #ifdef phased
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]) {
      for(int gamete=0;gamete<2;++gamete){
        float penetrance = penetrance_cache[subject*max_haplotypes*2+
        hapindex*2+gamete];
        if (penetrance>0){
          float prob = penetrance * local_frequency[hapindex];
          // COMPUTE POSTERIOR PROBS
          int m = 2 * gamete + local_center_dosage[hapindex] ;
          local_posterior_prob0[threadindex]+=(m==0)*prob;
          local_posterior_prob1[threadindex]+=(m==1)*prob;
          local_posterior_prob2[threadindex]+=(m==2)*prob;
          local_posterior_prob3[threadindex]+=(m==3)*prob;
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  #endif
// LOG2 REDUCTION
//  for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
//    if (threadindex<s && bestpair[threadindex+s].prob>bestpair[threadindex].prob){
//      bestpair[threadindex]= bestpair[threadindex+s];
//    }
//    barrier(CLK_LOCAL_MEM_FENCE);
//  }
  for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
    if (threadindex<s) {
       local_posterior_prob0[threadindex]+=local_posterior_prob0[threadindex+s];
       local_posterior_prob1[threadindex]+=local_posterior_prob1[threadindex+s];
       local_posterior_prob2[threadindex]+=local_posterior_prob2[threadindex+s];
       local_posterior_prob3[threadindex]+=local_posterior_prob3[threadindex+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if (threadindex<4){
    local_posterior_prob[threadindex] = epsilon;
  }
  if (threadindex==0){
//    int geno1 = local_center_dosage[bestpair[0].maternal] +
//              local_center_dosage[bestpair[0].paternal];
//    subject_genotype[subject] = geno1; 
    // MOVE THESE INTO AN ARRAY OF 3 FOR CONVENIENT ACCESS
    local_posterior_prob[0] += local_posterior_prob0[threadindex];
    local_posterior_prob[1] += local_posterior_prob1[threadindex];
    local_posterior_prob[2] += local_posterior_prob2[threadindex];
    local_posterior_prob[3] += local_posterior_prob3[threadindex];
    local_denom[0] = 0;
#ifdef unphased
    float max_prob = 0;
    int geno2 = 0;
    for(int h=0;h<geno_dim;++h){
      if(local_posterior_prob[h]>max_prob){
        max_prob = local_posterior_prob[h];
        geno2 = h;
      }
      local_denom[threadindex]+=local_posterior_prob[h];
    }
    //local_denom[threadindex] = 1;
    subject_genotype[subject] = geno2;
    if (genotype_imputation[0]){
      subject_dosage[subject] = (local_posterior_prob[1]+2*local_posterior_prob[2])/local_denom[threadindex];
    }
#endif
  }
  //barrier(CLK_LOCAL_MEM_FENCE);
  if(threadindex<geno_dim){
    //local_denom[0] = 1;
    subject_posterior_prob[subject*4+threadindex] = 
    //local_posterior_prob[threadindex];
    //local_denom[0];
    //subject_posterior_prob[subject*4+threadindex] = 
    local_posterior_prob[threadindex]/local_denom[0];
    //local_posterior_prob[threadindex];
    //epsilon;
  }
  return;
}
