
__kernel void impute_penetrance(
__const int max_haplotypes,
__const int penetrance_matrix_size,
__global int * active_haplotype,
__global int * twin_hap_index,
__global float * logpenetrance_cache,
__local int * local_active_haplotype,
__local int * local_twin_hap_index,
__local float * local_temp_vec
) {
  int subject = get_group_id(0);
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_temp_vec[hapindex] = 0;
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){
      local_twin_hap_index[hapindex] = twin_hap_index[hapindex];
    }
  }
  for(int rowindex = 0;rowindex<max_haplotypes;++rowindex){
    if(local_active_haplotype[rowindex]){
      // STEP 1 Do coalesced read from penetrance matrix
      for(int colchunk=0;colchunk<(max_haplotypes/BLOCK_WIDTH)+1;++colchunk){
        int colindex = colchunk*BLOCK_WIDTH+threadindex;
        if (colindex<max_haplotypes && local_active_haplotype[colindex]){
          local_temp_vec[colindex] = 
          logpenetrance_cache[subject*penetrance_matrix_size+rowindex*
          max_haplotypes+colindex];
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      // STEP 2 Copy from equivalent columns
      for(int colchunk=0;colchunk<(max_haplotypes/BLOCK_WIDTH)+1;++colchunk){
        int colindex = colchunk*BLOCK_WIDTH+threadindex;
        if (colindex<max_haplotypes && local_active_haplotype[colindex] && 
        local_twin_hap_index[colindex]!=colindex){
          local_temp_vec[colindex] = 
          local_temp_vec[local_twin_hap_index[colindex]];
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      // STEP 3 Do coalesced write to penetrance matrix
      for(int colchunk=0;colchunk<(max_haplotypes/BLOCK_WIDTH)+1;++colchunk){
        int colindex = colchunk*BLOCK_WIDTH+threadindex;
        if (colindex<max_haplotypes && local_active_haplotype[colindex]){
          logpenetrance_cache[subject*penetrance_matrix_size+rowindex*
          max_haplotypes+colindex] = local_temp_vec[colindex];
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  for(int rowindex = 0;rowindex<max_haplotypes;++rowindex){
    if(local_active_haplotype[rowindex]&&
    local_twin_hap_index[rowindex]!=rowindex){
      for(int colchunk=0;colchunk<(max_haplotypes/BLOCK_WIDTH)+1;++colchunk){
        int colindex = colchunk*BLOCK_WIDTH+threadindex;
        if (colindex<max_haplotypes && local_active_haplotype[colindex]){
          logpenetrance_cache[subject*penetrance_matrix_size+rowindex*
          max_haplotypes+colindex] = 
          logpenetrance_cache[subject*penetrance_matrix_size+
          local_twin_hap_index[rowindex]*max_haplotypes+colindex];
        }
      }
    }
  }
  return;
}

__kernel void impute_penetrance_haploid(
__const int max_haplotypes,
__const int penetrance_matrix_size,
__global int * active_haplotype,
__global int * twin_hap_index,
__global float * logpenetrance_cache,
__local int * local_active_haplotype,
__local int * local_twin_hap_index
) {
  int subject = get_group_id(0);
  int threadindex = get_local_id(0);
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
      local_twin_hap_index[hapindex] = twin_hap_index[hapindex];
    }
  }
  for(int rowindex = 0;rowindex<max_haplotypes;++rowindex){
    if(local_active_haplotype[rowindex]&&
    local_twin_hap_index[rowindex]!=rowindex){
      if (threadindex<2){
      //for(int gamete=0;gamete<2;++gamete){
        logpenetrance_cache[subject*2*max_haplotypes+
        rowindex*2+threadindex] =  
        logpenetrance_cache[subject*2*max_haplotypes+
        local_twin_hap_index[rowindex]*2+threadindex];
      }
    }
  }
  return;
}

__kernel void precompute_penetrance_fast(
__const int max_haplotypes,
__const int penetrance_matrix_size,
__const int max_region_size,
__const float logpenetrance_threshold,
__constant int * haplotypes,
__constant int * region_snp_offset,
__constant int * markers,
__constant int * prev_left_marker,
__constant int * left_marker,
__global int * haploid_arr,
__global int * beyond_left_edge_dosage,
__global int * right_edge_dosage,
__global float * region_snp_penetrance,
__global float * logpenetrance_cache,
__global float * penetrance_cache,
__global int * active_haplotype,
__local int * local_active_haplotype,
__local int * local_beyond_left_edge_dosage,
__local int * local_right_edge_dosage,
__local float * local_snp_penetrance_beyond_left,
__local float * local_snp_penetrance_right,
__local float * local_max_penetrance
) {
  int subject = get_group_id(0);
  int haploid = haploid_arr[subject];
  //int haploid = 0;
  int threadindex = get_local_id(0);
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
      local_beyond_left_edge_dosage[hapindex] = 
      beyond_left_edge_dosage[hapindex];
      local_right_edge_dosage[hapindex] = 
      right_edge_dosage[hapindex];
    }
  }
  // COPY SNP PENETRANCES FOR MARKERS IN WINDOW
  if (threadindex<3){
    local_snp_penetrance_beyond_left[threadindex] = region_snp_penetrance[subject*3*max_region_size+ 3*(prev_left_marker[0]-region_snp_offset[0]) + threadindex];
    local_snp_penetrance_right[threadindex] =  region_snp_penetrance[subject*3*max_region_size+ 3*(left_marker[0]+markers[0]-1 -region_snp_offset[0]) + threadindex];
  }
  local_max_penetrance[threadindex] = -1e10;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int righthapindex=0;righthapindex<max_haplotypes;++righthapindex){
    if (local_active_haplotype[righthapindex]){
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
        int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
        if (lefthapindex<=righthapindex&&local_active_haplotype[lefthapindex]
        &&(!haploid || lefthapindex==righthapindex)){
          float logpenetrance_old = prev_left_marker[0]!=left_marker[0] ?
          local_snp_penetrance_beyond_left[local_beyond_left_edge_dosage[
          righthapindex]+ local_beyond_left_edge_dosage[lefthapindex]] : 0;

          float logpenetrance = haplotypes[0]>2 ?
          logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*
          max_haplotypes+lefthapindex] + local_snp_penetrance_right[
          local_right_edge_dosage[righthapindex]+local_right_edge_dosage[lefthapindex]] - 
          logpenetrance_old : 
          local_snp_penetrance_right[local_right_edge_dosage[righthapindex]+
          local_right_edge_dosage[lefthapindex]];

          if (logpenetrance>local_max_penetrance[threadindex])
          local_max_penetrance[threadindex] = logpenetrance;        
          //logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = logpenetrance;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
    }
  }
  for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
    if (threadindex<s && local_max_penetrance[threadindex+s]>local_max_penetrance[threadindex]) {
      local_max_penetrance[threadindex]=local_max_penetrance[threadindex+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  //if (threadindex==0) local_max_penetrance[0] = -local_max_penetrance[0];
  for(int righthapindex=0;righthapindex<max_haplotypes;++righthapindex){
    if(local_active_haplotype[righthapindex]){
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
        int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
        if (lefthapindex<max_haplotypes&&local_active_haplotype[lefthapindex]&&
        (!haploid || lefthapindex==righthapindex)){

          float logpenetrance_old = prev_left_marker[0]!=left_marker[0] ?
          local_snp_penetrance_beyond_left[local_beyond_left_edge_dosage[
          righthapindex]+ local_beyond_left_edge_dosage[lefthapindex]] : 0;

          float logpenetrance = haplotypes[0]>2 ?
          logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*
          max_haplotypes+lefthapindex] + local_snp_penetrance_right[
          local_right_edge_dosage[righthapindex]+local_right_edge_dosage[lefthapindex]] - 
          logpenetrance_old : 
          local_snp_penetrance_right[local_right_edge_dosage[righthapindex]+
          local_right_edge_dosage[lefthapindex]];

          //float rescaled = logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex]+local_max_penetrance[0];
          //logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = rescaled;
          //penetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = rescaled>=logpenetrance_threshold?exp(rescaled):0;
          //float rescaled = logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex]+local_max_penetrance[0];
          logpenetrance-=local_max_penetrance[0];
          logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = logpenetrance;
          //logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = (left_marker[0]==2 && subject==600)?local_snp_penetrance_right[local_right_edge_dosage[righthapindex]+local_right_edge_dosage[lefthapindex]] - logpenetrance_old:logpenetrance;
          penetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = logpenetrance>=logpenetrance_threshold?exp(logpenetrance):0;
        }
      }
    }
  }
  return;
}

__kernel void precompute_penetrance_fast_haploid(
__const int max_haplotypes,
__const int penetrance_matrix_size,
__const int max_region_size,
__const float logpenetrance_threshold,
__constant int * haplotypes,
__constant int * region_snp_offset,
__constant int * markers,
__constant int * prev_left_marker,
__constant int * left_marker,
__global int * beyond_left_edge_dosage,
__global int * right_edge_dosage,
__global float * region_snp_penetrance,
__global float * logpenetrance_cache,
__global float * penetrance_cache,
__global int * active_haplotype,
__local int * local_active_haplotype,
__local int * local_beyond_left_edge_dosage,
__local int * local_right_edge_dosage,
__local float * local_snp_penetrance_beyond_left,
__local float * local_snp_penetrance_right,
__local float * local_max_penetrance,
__local float * local_penetrance_pair
) {
  int subject = get_group_id(0);
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
  for(int gamete=0;gamete<2;++gamete){
    local_penetrance_pair[2*threadindex+gamete] = 0;
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){

      local_beyond_left_edge_dosage[hapindex] = 
      prev_left_marker[0]!=left_marker[0] ?
      beyond_left_edge_dosage[hapindex] : 0 ;

      local_right_edge_dosage[hapindex] = 
      right_edge_dosage[hapindex];
    }
  }
  // COPY SNP PENETRANCES FOR MARKERS IN WINDOW
  if (threadindex<4){
    local_snp_penetrance_beyond_left[threadindex] = 
    prev_left_marker[0]!=left_marker[0] ?
    region_snp_penetrance[subject*4*max_region_size+ 
    4*(prev_left_marker[0]-region_snp_offset[0]) + threadindex] : 0;

    local_snp_penetrance_right[threadindex] =  region_snp_penetrance[subject*4*max_region_size+ 4*(left_marker[0]+markers[0]-1 -region_snp_offset[0]) + threadindex];
  }
  local_max_penetrance[threadindex] = -1e10;
  barrier(CLK_LOCAL_MEM_FENCE);
  //for(int righthapindex=0;righthapindex<max_haplotypes;++righthapindex){
  //  if (local_active_haplotype[righthapindex]){
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes&&local_active_haplotype[hapindex]){
      //int dose = local_beyond_left_edge_dosage[hapindex];
      for(int gamete=0;gamete<2;++gamete){
        if (haplotypes[0]>2){
          local_penetrance_pair[2*threadindex+gamete] = 
          logpenetrance_cache[subject*max_haplotypes*2+hapindex*2
          +gamete] + 
          local_snp_penetrance_right[2*gamete+local_right_edge_dosage[hapindex]] 
          -local_snp_penetrance_beyond_left[
          2*gamete+local_beyond_left_edge_dosage[hapindex]]
          ;
        }else if (haplotypes[0]==16){
          //logpenetrance_cache[subject*max_haplotypes*2+hapindex*2+gamete] = 
          //hapindex;
          //local_right_edge_dosage[hapindex];
        }else{
          local_penetrance_pair[2*threadindex+gamete] = 
          local_snp_penetrance_right[2*gamete+local_right_edge_dosage[hapindex]];
        }
        if (local_penetrance_pair[2*threadindex+gamete]>
        local_max_penetrance[threadindex]) 
        local_max_penetrance[threadindex] = 
        local_penetrance_pair[2*threadindex+gamete];
        //logpenetrance_cache[subject*max_haplotypes*2+hapindex*2+gamete] = 
        //local_right_edge_dosage[hapindex];
        //local_penetrance_pair[2*threadindex+gamete];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  //if (haplotypes[0]==16) return;
  for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
    if (threadindex<s && local_max_penetrance[threadindex+s]>local_max_penetrance[threadindex]) {
      local_max_penetrance[threadindex]=local_max_penetrance[threadindex+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  //return;
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes&&local_active_haplotype[hapindex]){
      for(int gamete=0;gamete<2;++gamete){
        if (haplotypes[0]>2){
          local_penetrance_pair[2*threadindex+gamete] = 
          logpenetrance_cache[subject*max_haplotypes*2+hapindex*2
          +gamete] + local_snp_penetrance_right[
          2*gamete+local_right_edge_dosage[hapindex]] - 
          local_snp_penetrance_beyond_left[
          2*gamete+local_beyond_left_edge_dosage[hapindex]]
          - local_max_penetrance[0];
        }else{
          local_penetrance_pair[2*threadindex+gamete] = 
          local_snp_penetrance_right[2*gamete+local_right_edge_dosage[hapindex]]
          - local_max_penetrance[0];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        logpenetrance_cache[subject*max_haplotypes*2+hapindex*2+gamete] = 
        local_penetrance_pair[2*threadindex+gamete];
        penetrance_cache[subject*max_haplotypes*2+hapindex*2+gamete] =
        local_penetrance_pair[2*threadindex+gamete]>
        logpenetrance_threshold?
        exp(local_penetrance_pair[2*threadindex+gamete]):0;
      }
    }
  }
  return;
}

__kernel void impute_genotype_denovo(
__const int max_haplotypes,
__const int max_window,
__const float epsilon,
__const int penetrance_matrix_size,
__constant int * genotype_imputation,
__constant int * markers,
__constant int * haplotypes,
__constant int * left_marker,
__global int * haploid_arr,
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
__local bestpair_t * bestpair,
__local float * local_posterior_prob0,
__local float * local_posterior_prob1,
__local float * local_posterior_prob2,
__local float * local_posterior_prob3,
__local float * local_posterior_prob,
__local float * local_cached_marginals, 
__local float * local_denom
) {
  int subject = get_group_id(0);
  int haploid = haploid_arr[subject];
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
  if(threadindex==0) local_denom[threadindex]=0;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){
      local_frequency[hapindex] = frequency[hapindex];
      local_center_dosage[hapindex] = center_dosage[hapindex];
      local_cached_marginals[hapindex] = 
      subject_haplotype_weight[subject*max_haplotypes+hapindex];
    }
  }
  local_posterior_prob0[threadindex]=0;
  local_posterior_prob1[threadindex]=0;
  local_posterior_prob2[threadindex]=0;
  local_posterior_prob3[threadindex]=0;
  //bestpair[threadindex].prob = 0;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int righthapindex=0;righthapindex<max_haplotypes;++righthapindex){
    if (local_active_haplotype[righthapindex] && 
    local_cached_marginals[righthapindex]>0){
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
        int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
        if (lefthapindex<=righthapindex&&local_active_haplotype[lefthapindex]
        && local_cached_marginals[lefthapindex]>0 && (!haploid||lefthapindex==righthapindex)){
          float penetrance = penetrance_cache[subject*penetrance_matrix_size+
          righthapindex*max_haplotypes+lefthapindex];
          if (penetrance>0){
//            float prob = penetrance;
            float prob = (lefthapindex!=righthapindex)?
            local_frequency[lefthapindex]*
            local_frequency[righthapindex]*penetrance*2:
            local_frequency[lefthapindex]*
            local_frequency[righthapindex]*penetrance;

            // DETERMINE BEST HAP PAIR 
//            if (prob>bestpair[threadindex].prob){
//              bestpair[threadindex].maternal = righthapindex;
//              bestpair[threadindex].paternal = lefthapindex;
//              bestpair[threadindex].prob = prob;
//            }


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
  if (threadindex==0){
//    int geno1 = local_center_dosage[bestpair[0].maternal] +
//              local_center_dosage[bestpair[0].paternal];
//    subject_genotype[subject] = geno1; 
    // MOVE THESE INTO AN ARRAY OF 3 FOR CONVENIENT ACCESS
    local_posterior_prob[0] = local_posterior_prob0[threadindex];
    local_posterior_prob[1] = local_posterior_prob1[threadindex];
    local_posterior_prob[2] = local_posterior_prob2[threadindex];
    local_posterior_prob[3] = local_posterior_prob3[threadindex];
    float max_prob = 0;
    int geno2 = 0;
    for(int h=0;h<4;++h){
      if(local_posterior_prob[h]>max_prob){
        max_prob = local_posterior_prob[h];
        geno2 = h;
      }
      local_denom[threadindex]+=local_posterior_prob[h];
    }
    subject_genotype[subject] = geno2;
    if (genotype_imputation[0]){
      subject_dosage[subject] = (local_posterior_prob[1]+2*local_posterior_prob[2])/local_denom[threadindex];
    }
  }
  //barrier(CLK_LOCAL_MEM_FENCE);
  if(threadindex<4){
    subject_posterior_prob[subject*4+threadindex] = 
    local_posterior_prob[threadindex]/local_denom[0];
    //local_posterior_prob[threadindex];
  }
  return;
}

__kernel void impute_genotype_denovo_haploid(
__const int max_haplotypes,
__const int max_window,
__const float epsilon,
__const int penetrance_matrix_size,
__constant int * genotype_imputation,
__constant int * markers,
__constant int * haplotypes,
__constant int * left_marker,
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
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
  if(threadindex==0) local_denom[threadindex]=0;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){
      local_frequency[hapindex] = frequency[hapindex];
      local_center_dosage[hapindex] = center_dosage[hapindex];
      local_cached_marginals[hapindex] = 
      subject_haplotype_weight[subject*max_haplotypes+hapindex];
    }
  }
  local_posterior_prob0[threadindex]=0;
  local_posterior_prob1[threadindex]=0;
  local_posterior_prob2[threadindex]=0;
  local_posterior_prob3[threadindex]=0;
  barrier(CLK_LOCAL_MEM_FENCE);
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
// LOG2 REDUCTION
  for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
    if (threadindex<s) {
       local_posterior_prob0[threadindex]+=local_posterior_prob0[threadindex+s];
       local_posterior_prob1[threadindex]+=local_posterior_prob1[threadindex+s];
       local_posterior_prob2[threadindex]+=local_posterior_prob2[threadindex+s];
       local_posterior_prob3[threadindex]+=local_posterior_prob3[threadindex+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if (threadindex==0){
    // MOVE THESE INTO AN ARRAY OF 3 FOR CONVENIENT ACCESS
    local_posterior_prob[0] = local_posterior_prob0[threadindex];
    local_posterior_prob[1] = local_posterior_prob1[threadindex];
    local_posterior_prob[2] = local_posterior_prob2[threadindex];
    local_posterior_prob[3] = local_posterior_prob3[threadindex];
  }
  if(threadindex<4){
    subject_posterior_prob[subject*4+threadindex] = 
    local_posterior_prob[threadindex];
  }
  return;
}
