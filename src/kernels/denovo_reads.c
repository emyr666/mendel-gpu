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
__local float * hap1_logprob,
__local float * hap2_logprob
) {
  int hap1 = get_group_id(0);
  int hap2 = get_group_id(1);
  int subject = get_group_id(2);
  int threadindex = get_local_id(0);
  if (!active_haplotype[hap1] || !active_haplotype[hap2]) return ;
  int total_depth = mat_rows_by_subject[subject];
  float loglike = 0;
  if(threadindex==0){
    logpenetrance_cache[subject*penetrance_matrix_size+hap1*max_haplotypes+hap2] =  loglike;
  }
  for(int row=0;row<total_depth;++row){
    hap1_logprob[threadindex] = 0;
    hap2_logprob[threadindex] = 0;
    barrier(CLK_LOCAL_MEM_FENCE);
    for(int chunk=0;chunk<(max_window/BLOCK_WIDTH)+1;++chunk){
      int col = chunk*BLOCK_WIDTH+threadindex;
      if (col<markers[0]){
        int read_allele = read_alleles_mat[subject*read_compact_matrix_size+row*max_window+col];
        if(read_allele!=9){
          float match = read_match_logmat[subject*read_compact_matrix_size+row*max_window+col];
          float mismatch = read_mismatch_logmat[subject*read_compact_matrix_size+row*max_window+col];
          int hap1_allele = haplotype[hap1*max_window+col];
          int hap2_allele = haplotype[hap2*max_window+col];
          hap1_logprob[threadindex]+=read_allele==hap1_allele?match:mismatch;
          hap2_logprob[threadindex]+=read_allele==hap2_allele?match:mismatch;
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
    for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
      if (threadindex<s) {
        hap1_logprob[threadindex]+=hap1_logprob[threadindex+s];
        hap2_logprob[threadindex]+=hap2_logprob[threadindex+s];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(threadindex==0){
      float mean = (hap1_logprob[0]+hap2_logprob[0])/2;
      loglike+=  LOG_HALF+(log(exp(hap1_logprob[0]-mean)+exp(hap2_logprob[0]-mean))+mean);
    }
  }
  if(threadindex==0){
    logpenetrance_cache[subject*penetrance_matrix_size+hap1*max_haplotypes+hap2] =  loglike;
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
  for(int righthapindex=0;righthapindex<max_haplotypes;++righthapindex){
    if (local_active_haplotype[righthapindex]){
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
        int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
        if (lefthapindex<=righthapindex&&local_active_haplotype[lefthapindex]){
          float logpenetrance = 
          logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*
          max_haplotypes+lefthapindex] ;

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
  for(int righthapindex=0;righthapindex<max_haplotypes;++righthapindex){
    if(local_active_haplotype[righthapindex]){
      for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
        int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
        if (lefthapindex<max_haplotypes&&local_active_haplotype[lefthapindex]){

          float logpenetrance = 
          logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*
          max_haplotypes+lefthapindex] - local_max_penetrance[0];
          logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = logpenetrance;
          penetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = logpenetrance>=logpenetrance_threshold?exp(logpenetrance):0;
        }
      }
    }
  }
  return;
}

