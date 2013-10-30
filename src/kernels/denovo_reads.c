__kernel void reads_compute_penetrance(
__const int max_window,
__const int max_haplotypes,
__const int penetrance_matrix_size,
__const int read_compact_matrix_size,
__constant int * markers,
__global int * mat_rows_by_subject,
__global int * read_alleles_mat,
__global float * read_match_logmat,
__global float * read_mismatch_logmat,
__global int * haplotype,
__global float * logpenetrance_cache,
__global int * active_haplotype,
__local int * x_hap_alleles,
__local int * y_hap_alleles,
__local float * x_hap_logprob,
__local float * y_hap_logprob
) {
  int x_hap = get_group_id(0);
  int y_hap = get_group_id(1);
  int subject = get_group_id(2);
  int threadindex = get_local_id(0);
  // compute for lower triangle only
  if (x_hap>y_hap) return;
  if (!active_haplotype[x_hap] || !active_haplotype[y_hap]) return ;
  int total_depth = mat_rows_by_subject[subject];
  float loglike = 0;
  for(int chunk=0;chunk<(max_window/SMALL_BLOCK_WIDTH)+1;++chunk){
    int col = chunk*SMALL_BLOCK_WIDTH+threadindex;
    if (col<markers[0]){
      x_hap_alleles[col] = haplotype[x_hap*max_window+col];
      y_hap_alleles[col] = haplotype[y_hap*max_window+col];
    }
  }
  //if(threadindex==0){
    //logpenetrance_cache[subject*penetrance_matrix_size+y_hap*max_haplotypes+x_hap] =  loglike;
  //}
  //return;
  //for(int row=0;row<1;++row){
  for(int row=0;row<total_depth;++row){
    x_hap_logprob[threadindex] = 0;
    y_hap_logprob[threadindex] = 0;
    barrier(CLK_LOCAL_MEM_FENCE);
    for(int chunk=0;chunk<(max_window/SMALL_BLOCK_WIDTH)+1;++chunk){
      int col = chunk*SMALL_BLOCK_WIDTH+threadindex;
      if (col<markers[0]){
        int read_allele = read_alleles_mat[subject*read_compact_matrix_size+row*max_window+col];
        float match = read_match_logmat[subject*read_compact_matrix_size+row*max_window+col];
        float mismatch = read_mismatch_logmat[subject*read_compact_matrix_size+row*max_window+col];
        x_hap_logprob[threadindex]+=(read_allele!=9&&read_allele==x_hap_alleles[col])*match;
        x_hap_logprob[threadindex]+=(read_allele!=9&&read_allele!=x_hap_alleles[col])*mismatch;
        y_hap_logprob[threadindex]+=(read_allele!=9&&read_allele==y_hap_alleles[col])*match;
        y_hap_logprob[threadindex]+=(read_allele!=9&&read_allele!=y_hap_alleles[col])*mismatch;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
    for(int s=SMALL_BLOCK_WIDTH/2; s>0; s>>=1) {
      if (threadindex<s) {
        x_hap_logprob[threadindex]+=x_hap_logprob[threadindex+s];
        y_hap_logprob[threadindex]+=y_hap_logprob[threadindex+s];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(threadindex==0){
      float mean = (x_hap_logprob[0]+y_hap_logprob[0])/2;
      loglike+=  LOG_HALF+(log(exp(x_hap_logprob[0]-mean)+exp(y_hap_logprob[0]-mean))+mean);
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

