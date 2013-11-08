__kernel void reads_compute_penetrance(
__const int max_window,
__const int max_haplotypes,
__const int penetrance_matrix_size,
__const int read_compact_matrix_size,
__const int max_vector_length,
__const int compact_rows,
__constant int * markers,
__global int * mat_rows_by_subject,
__global int * read_alleles_vec,
__global float * read_match_logvec,
__global float * read_mismatch_logvec,
__global int * vector_offsets,
__global int * read_lengths,
__global int * haplotype_offsets,
__global int * haplotype,
__global float * logpenetrance_cache,
__global int * active_haplotype,
__global int * best_haplotypes,
__local int * local_read_alleles_vec,
__local float * local_read_match_logvec,
__local float * local_read_mismatch_logvec,
__local int * local_vector_offsets,
__local int * local_read_lengths,
__local int * local_haplotype_offsets,
__local int * local_x_hap_alleles,
__local int * local_y_hap_alleles,
__local int * local_temp_read_alleles,
__local float * local_temp_match,
__local float * local_temp_mismatch,
__local float * local_x_hap_logprob,
__local float * local_y_hap_logprob
) {
  int x_hap = get_group_id(0);
  int y_hap = get_group_id(1);
  int subject = get_group_id(2);
  int threadindex = get_local_id(0);
  // compute for lower triangle only
  if (x_hap>y_hap) return;
  //if(x_hap!=0 || y_hap!=1) return;
  //if (!active_haplotype[x_hap] || !active_haplotype[y_hap] ){
  if (!active_haplotype[x_hap] || !active_haplotype[y_hap] 
  ||!best_haplotypes[subject*max_haplotypes+x_hap] || !best_haplotypes[subject*max_haplotypes+y_hap]) {
    if (threadindex==0)
      logpenetrance_cache[subject*penetrance_matrix_size+y_hap*max_haplotypes+x_hap] = -999; 
    return ;
  }
  int total_depth = mat_rows_by_subject[subject];
  for(int chunk=0;chunk<(max_window/SMALL_BLOCK_WIDTH)+1;++chunk){
    int col = chunk*SMALL_BLOCK_WIDTH+threadindex;
    if (col<markers[0]){
      local_x_hap_alleles[col] = haplotype[x_hap*max_window+col];
      local_y_hap_alleles[col] = haplotype[y_hap*max_window+col];
    }
  }
  for(int chunk=0;chunk<(compact_rows/SMALL_BLOCK_WIDTH)+1;++chunk){
    int col = chunk*SMALL_BLOCK_WIDTH+threadindex;
    if (col<compact_rows){
      local_vector_offsets[col] = vector_offsets[subject*compact_rows+col];
      local_read_lengths[col] = read_lengths[subject*compact_rows+col];
      local_haplotype_offsets[col] = haplotype_offsets[subject*compact_rows+col];
    }
  }
  for(int chunk=0;chunk<(max_vector_length/SMALL_BLOCK_WIDTH)+1;++chunk){
    int col = chunk*SMALL_BLOCK_WIDTH+threadindex;
    if (col<max_vector_length){
      local_read_alleles_vec[col] = read_alleles_vec[subject*max_vector_length+col];
      local_read_match_logvec[col] = read_match_logvec[subject*max_vector_length+col];
      local_read_mismatch_logvec[col] = read_mismatch_logvec[subject*max_vector_length+col];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  float loglike = 0;
  //if(threadindex==0){
    //logpenetrance_cache[subject*penetrance_matrix_size+y_hap*max_haplotypes+x_hap] =  y_hap*max_haplotypes+x_hap<max_vector_length?local_read_mismatch_logvec[y_hap*max_haplotypes+x_hap]:0;
  //}
  //return;
  for(int row=0;row<total_depth;++row){
    if(local_read_lengths[row]>0 
    && local_vector_offsets[row]<max_vector_length){
      // initialize the local arrays
      local_x_hap_logprob[threadindex] = 0;
      local_y_hap_logprob[threadindex] = 0;
      for(int chunk=0;chunk<(max_window/SMALL_BLOCK_WIDTH)+1;++chunk){
        int col = chunk*SMALL_BLOCK_WIDTH+threadindex;
        if (col<max_window){
          local_temp_read_alleles[col] = 9;
          local_temp_match[col] = 1;
          local_temp_mismatch[col] = 1;
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      // store in local array from vectorized reads, and align with the 
      // reference haplotypes
      for(int chunk=0;chunk<(max_window/SMALL_BLOCK_WIDTH)+1;++chunk){
        int col = chunk*SMALL_BLOCK_WIDTH+threadindex;
        if (col < local_read_lengths[row] 
        && local_haplotype_offsets[row]+col<max_window 
        && local_vector_offsets[row]+col<max_vector_length){
          local_temp_read_alleles[local_haplotype_offsets[row]+col] = local_read_alleles_vec[local_vector_offsets[row]+col];
          local_temp_match[local_haplotype_offsets[row]+col] = local_read_match_logvec[local_vector_offsets[row]+col];
          local_temp_mismatch[local_haplotype_offsets[row]+col] = local_read_mismatch_logvec[local_vector_offsets[row]+col];
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      // accumulate evidence in the local array
      for(int chunk=0;chunk<(max_window/SMALL_BLOCK_WIDTH)+1;++chunk){
        int col = chunk*SMALL_BLOCK_WIDTH+threadindex;
        if (col<markers[0]){
          local_x_hap_logprob[threadindex]+=(local_temp_read_alleles[col]!=9&&local_temp_read_alleles[col]==local_x_hap_alleles[col])*local_temp_match[col];
          local_x_hap_logprob[threadindex]+=(local_temp_read_alleles[col]!=9&&local_temp_read_alleles[col]!=local_x_hap_alleles[col])*local_temp_mismatch[col];
          local_y_hap_logprob[threadindex]+=(local_temp_read_alleles[col]!=9&&local_temp_read_alleles[col]==local_y_hap_alleles[col])*local_temp_match[col];
          local_y_hap_logprob[threadindex]+=(local_temp_read_alleles[col]!=9&&local_temp_read_alleles[col]!=local_y_hap_alleles[col])*local_temp_mismatch[col];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      // reduce the evidence
      for(int s=SMALL_BLOCK_WIDTH/2; s>0; s>>=1) {
        if (threadindex<s) {
          local_x_hap_logprob[threadindex]+=local_x_hap_logprob[threadindex+s];
          local_y_hap_logprob[threadindex]+=local_y_hap_logprob[threadindex+s];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
      }
      if(threadindex==0){
        float mean = (local_x_hap_logprob[0]+local_y_hap_logprob[0])/2;
        loglike+=  LOG_HALF+(log(exp(local_x_hap_logprob[0]-mean)+exp(local_y_hap_logprob[0]-mean))+mean);
      }
    }
  }
  if(threadindex==0){
    logpenetrance_cache[subject*penetrance_matrix_size+y_hap*max_haplotypes+x_hap] =  loglike;
  }
  return;
}

__kernel void reads_adjust_penetrance(
__const int max_haplotypes,
__const int penetrance_matrix_size,
__const float logpenetrance_threshold,
__global float * logpenetrance_cache,
__global float * penetrance_cache,
__global int * active_haplotype,
__local int * local_active_haplotype,
__local float * local_max_penetrance
) {
  int subject = get_group_id(0);
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<max_haplotypes){
      local_active_haplotype[hapindex] = active_haplotype[hapindex];
    }
  }
  local_max_penetrance[threadindex] = -1e10;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int y_hap_index=0;y_hap_index<max_haplotypes;++y_hap_index){
    if (local_active_haplotype[y_hap_index]){
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
        int x_hap_index = chunk*BLOCK_WIDTH+threadindex;
        if (x_hap_index<=y_hap_index&&local_active_haplotype[x_hap_index]){
          float logpenetrance = 
          logpenetrance_cache[subject*penetrance_matrix_size+y_hap_index*
          max_haplotypes+x_hap_index] ;

          if (logpenetrance>local_max_penetrance[threadindex])
          local_max_penetrance[threadindex] = logpenetrance;        
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
  for(int y_hap_index=0;y_hap_index<max_haplotypes;++y_hap_index){
    if(local_active_haplotype[y_hap_index]){
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
        int x_hap_index = chunk*BLOCK_WIDTH+threadindex;
        if (x_hap_index<=y_hap_index&&local_active_haplotype[x_hap_index]){

          float logpenetrance = 
          logpenetrance_cache[subject*penetrance_matrix_size+y_hap_index*
          max_haplotypes+x_hap_index] - local_max_penetrance[0];
          float penetrance = logpenetrance>=logpenetrance_threshold?
          exp(logpenetrance):0;
          logpenetrance_cache[subject*penetrance_matrix_size+y_hap_index*max_haplotypes+x_hap_index] = logpenetrance;
          logpenetrance_cache[subject*penetrance_matrix_size+x_hap_index*max_haplotypes+y_hap_index] = logpenetrance;
          penetrance_cache[subject*penetrance_matrix_size+y_hap_index*max_haplotypes+x_hap_index] = penetrance;
          penetrance_cache[subject*penetrance_matrix_size+x_hap_index*max_haplotypes+y_hap_index] = penetrance;
        }
      }
    }
  }
  return;
}

