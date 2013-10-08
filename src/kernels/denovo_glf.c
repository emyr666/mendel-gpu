__kernel void precompute_penetrance_fast(
__const int max_haplotypes,
__const int penetrance_matrix_size,
__const int snps,
__const float logpenetrance_threshold,
__constant int * haplotypes,
__constant int * markers,
__constant int * prev_left_marker,
__constant int * left_marker,
#ifdef unphased
__global int * haploid_arr,
#endif
__global int * beyond_left_edge_dosage,
__global int * right_edge_dosage,
__global float * snp_penetrance,
__global float * logpenetrance_cache,
__global float * penetrance_cache,
__global int * active_haplotype,
__local int * local_active_haplotype,
__local int * local_beyond_left_edge_dosage,
__local int * local_right_edge_dosage,
__local float * local_snp_penetrance_beyond_left,
__local float * local_snp_penetrance_right,
#ifdef phased
__local float * local_penetrance_pair,
#endif
__local float * local_max_penetrance
) {
  int subject = get_group_id(0);
#ifdef unphased
  int haploid = haploid_arr[subject];
#endif
  //int haploid = 0;
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
#ifdef phased
  for(int gamete=0;gamete<2;++gamete){
    local_penetrance_pair[2*threadindex+gamete] = 0;
  }
#endif
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes && local_active_haplotype[hapindex]){
#ifdef unphased
      local_beyond_left_edge_dosage[hapindex] = 
      beyond_left_edge_dosage[hapindex];
#endif
#ifdef phased
      local_beyond_left_edge_dosage[hapindex] = 
      prev_left_marker[0]!=left_marker[0] ?
      beyond_left_edge_dosage[hapindex] : 0 ;
#endif
      local_right_edge_dosage[hapindex] = 
      right_edge_dosage[hapindex];
    }
  }
  // COPY SNP PENETRANCES FOR MARKERS IN WINDOW
#ifdef unphased
  if (threadindex<3){
    local_snp_penetrance_beyond_left[threadindex] = snp_penetrance[subject*3*snps+ 3*(prev_left_marker[0]) + threadindex];
    local_snp_penetrance_right[threadindex] =  snp_penetrance[subject*3*snps+ 3*(left_marker[0]+markers[0]-1 ) + threadindex];
  }
#endif
#ifdef phased
  if (threadindex<4){
    local_snp_penetrance_beyond_left[threadindex] = 
    prev_left_marker[0]!=left_marker[0] ?
    snp_penetrance[subject*4*3+ 
    4*(prev_left_marker[0]) + threadindex] : 0;

    local_snp_penetrance_right[threadindex] =  snp_penetrance[subject*4*snps+ 4*(left_marker[0]+markers[0]-1 ) + threadindex];
  }
#endif
  local_max_penetrance[threadindex] = -1e10;
  barrier(CLK_LOCAL_MEM_FENCE);
#ifdef unphased
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
#endif
#ifdef phased
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
#endif
  for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
    if (threadindex<s && local_max_penetrance[threadindex+s]>local_max_penetrance[threadindex]) {
      local_max_penetrance[threadindex]=local_max_penetrance[threadindex+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  //if (threadindex==0) local_max_penetrance[0] = -local_max_penetrance[0];
#ifdef unphased
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

          logpenetrance-=local_max_penetrance[0];
          logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = logpenetrance;
          penetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = logpenetrance>=logpenetrance_threshold?exp(logpenetrance):0;
        }
      }
    }
  }
#endif
#ifdef phased
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
#endif
  return;
}


__kernel void impute_penetrance(
__const int max_haplotypes,
__const int penetrance_matrix_size,
__global int * active_haplotype,
__global int * twin_hap_index,
__global float * logpenetrance_cache,
__local int * local_active_haplotype,
#ifdef unphased
__local float * local_temp_vec,
#endif
__local int * local_twin_hap_index
) {
  int subject = get_group_id(0);
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
#ifdef unphased
      local_temp_vec[hapindex] = 0;
#endif
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
#ifdef unphased
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
#endif
  for(int rowindex = 0;rowindex<max_haplotypes;++rowindex){
    if(local_active_haplotype[rowindex]&&
    local_twin_hap_index[rowindex]!=rowindex){
#ifdef unphased
      for(int colchunk=0;colchunk<(max_haplotypes/BLOCK_WIDTH)+1;++colchunk){
        int colindex = colchunk*BLOCK_WIDTH+threadindex;
        if (colindex<max_haplotypes && local_active_haplotype[colindex]){
          logpenetrance_cache[subject*penetrance_matrix_size+rowindex*
          max_haplotypes+colindex] = 
          logpenetrance_cache[subject*penetrance_matrix_size+
          local_twin_hap_index[rowindex]*max_haplotypes+colindex];
        }
      }
#endif
#ifdef phased
      if (threadindex<2){
      //for(int gamete=0;gamete<2;++gamete){
        logpenetrance_cache[subject*2*max_haplotypes+
        rowindex*2+threadindex] =  
        logpenetrance_cache[subject*2*max_haplotypes+
        local_twin_hap_index[rowindex]*2+threadindex];
      }
#endif
    }
  }
  return;
}

