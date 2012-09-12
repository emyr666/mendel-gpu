
__kernel void precompute_penetrance(
__const int max_haplotypes,
__const int max_window,
__const int penetrance_matrix_size,
__const int max_region_size,
__const float logpenetrance_threshold,
__const int packedhap_len,
__constant int * region_snp_offset,
__constant int * markers,
__constant int * haplotypes,
__constant int * left_marker,
__global packedhap_t * packedhap,
__global int * haplotype,
__global float * region_snp_penetrance,
__global float * logpenetrance_cache,
__global float * penetrance_cache,
__global int * extended_snp_mapping,
__local int * local_extended_snp_mapping,
__local packedhap_t * local_packedhap,
__local int * local_hap1,
__local int * local_hap2,
__local float * local_snp_penetrance,
__local float * local_max_penetrance
) {
  int subject = get_group_id(0);
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(markers[0]/BLOCK_WIDTH)+1;++chunk){
    int markerindex = chunk*BLOCK_WIDTH+threadindex;
    if(markerindex<markers[0]){
      local_extended_snp_mapping[markerindex] = 
      extended_snp_mapping[markerindex];
    }
  }
  // COPY HAPLOTYPE GROUP VERBATIM
  for(int haploindex = 0;haploindex<haplotypes[0];++haploindex){
    if(threadindex<packedhap_len){
      local_packedhap[haploindex*packedhap_len+threadindex] = packedhap[haploindex*packedhap_len+threadindex];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  // COPY SNP PENETRANCES FOR MARKERS IN WINDOW
  if (threadindex<3){
    for(int markerindex=0;markerindex<markers[0];++markerindex){
      local_snp_penetrance[markerindex*3+threadindex] = region_snp_penetrance[subject*3*max_region_size+ 3*(left_marker[0] + local_extended_snp_mapping[markerindex]-region_snp_offset[0]) + threadindex];
    }
  }
  local_max_penetrance[threadindex] = -999;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int righthapindex=0;righthapindex<haplotypes[0];++righthapindex){
    for(int chunk=0;chunk<(haplotypes[0]/BLOCK_WIDTH)+1;++chunk){
      int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
      if (lefthapindex<haplotypes[0]){
        float logpenetrance = 0;
        for(int windowmarker=0;windowmarker<markers[0];++windowmarker){
          int geno = (((int)local_packedhap[righthapindex*packedhap_len+(windowmarker/32)].octet[(windowmarker%32)/8]) >> ((windowmarker%32)%8) & 1) + (((int)local_packedhap[lefthapindex*packedhap_len+(windowmarker/32)].octet[(windowmarker%32)/8]) >> ((windowmarker%32)%8) & 1) ;
          logpenetrance+=local_snp_penetrance[windowmarker*3+geno];
        }
        if (logpenetrance>local_max_penetrance[threadindex])
          local_max_penetrance[threadindex] = logpenetrance;        
        logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = logpenetrance;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
    if (threadindex<s && local_max_penetrance[threadindex+s]>local_max_penetrance[threadindex]) {
      local_max_penetrance[threadindex]=local_max_penetrance[threadindex+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  for(int righthapindex=0;righthapindex<haplotypes[0];++righthapindex){
    for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
      int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
      if (lefthapindex<haplotypes[0]){
        float rescaled = logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex]-local_max_penetrance[0];
        logpenetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = rescaled;
        penetrance_cache[subject*penetrance_matrix_size+righthapindex*max_haplotypes+lefthapindex] = rescaled>=logpenetrance_threshold?exp(rescaled):0;
      }
    }
  }
  return;
}

__kernel void precompute_penetrance_haploid(
__const int max_haplotypes,
__const int max_window,
__const int penetrance_matrix_size,
__const int max_region_size,
__const float logpenetrance_threshold,
__const int packedhap_len,
__constant int * region_snp_offset,
__constant int * markers,
__constant int * haplotypes,
__constant int * left_marker,
__global packedhap_t * packedhap,
__global int * haplotype,
__global float * region_snp_penetrance,
__global float * logpenetrance_cache,
__global float * penetrance_cache,
__global int * extended_snp_mapping,
__local int * local_extended_snp_mapping,
__local packedhap_t * local_packedhap,
__local int * local_hap1,
__local float * local_snp_penetrance,
__local float * local_max_penetrance,
__local float * local_penetrance_pair
) {
  int subject = get_group_id(0);
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(markers[0]/BLOCK_WIDTH)+1;++chunk){
    int markerindex = chunk*BLOCK_WIDTH+threadindex;
    if(markerindex<markers[0]){
      local_extended_snp_mapping[markerindex] = 
      extended_snp_mapping[markerindex];
    }
  }
  // COPY HAPLOTYPE GROUP VERBATIM
  for(int haploindex = 0;haploindex<haplotypes[0];++haploindex){
    if(threadindex<packedhap_len){
      local_packedhap[haploindex*packedhap_len+threadindex] = packedhap[haploindex*packedhap_len+threadindex];
    }
  }
  for(int gamete=0;gamete<2;++gamete){
    local_penetrance_pair[2*threadindex+gamete] = 0;    
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  // COPY SNP PENETRANCES FOR MARKERS IN WINDOW
  if (threadindex<4){
    for(int markerindex=0;markerindex<markers[0];++markerindex){
      local_snp_penetrance[markerindex*4+threadindex] = region_snp_penetrance[subject*4*max_region_size+ 4*(left_marker[0] + local_extended_snp_mapping[markerindex]-region_snp_offset[0]) + threadindex];
    }
  }
  local_max_penetrance[threadindex] = -999;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int chunk=0;chunk<(haplotypes[0]/BLOCK_WIDTH)+1;++chunk){
    int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
    if (lefthapindex<haplotypes[0]){
      for(int windowmarker=0;windowmarker<markers[0];++windowmarker){
        int dose = (((int)local_packedhap[lefthapindex*packedhap_len+
        (windowmarker/32)].octet[(windowmarker%32)/8]) >> 
        ((windowmarker%32)%8) & 1) ;
        for(int gamete=0;gamete<2;++gamete){
          local_penetrance_pair[2*threadindex+gamete]+=local_snp_penetrance[windowmarker*4+2*gamete+dose];
        }
      }
      for(int gamete=0;gamete<2;++gamete){
          if (local_penetrance_pair[2*threadindex+gamete]>
          local_max_penetrance[threadindex]) 
          local_max_penetrance[threadindex] = 
          local_penetrance_pair[2*threadindex+gamete];
          logpenetrance_cache[subject*2*max_haplotypes+2*lefthapindex+gamete] 
          = local_penetrance_pair[2*threadindex+gamete];
      }
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  //return;
  for(int s=BLOCK_WIDTH/2; s>0; s>>=1) {
    if (threadindex<s && local_max_penetrance[threadindex+s]>local_max_penetrance[threadindex]) {
      local_max_penetrance[threadindex]=local_max_penetrance[threadindex+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  for(int chunk=0;chunk<(max_haplotypes/BLOCK_WIDTH)+1;++chunk){
    int lefthapindex = chunk*BLOCK_WIDTH+threadindex;
    if (lefthapindex<haplotypes[0]){
      for(int gamete=0;gamete<2;++gamete){
        local_penetrance_pair[2*threadindex+gamete] = 
        logpenetrance_cache[subject*2*max_haplotypes+2*lefthapindex + gamete]
        -local_max_penetrance[0];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      for(int gamete=0;gamete<2;++gamete){
        logpenetrance_cache[subject*2*max_haplotypes+2*lefthapindex + gamete] 
        = local_penetrance_pair[2*threadindex+gamete];
        penetrance_cache[subject*2*max_haplotypes+2*lefthapindex + gamete] 
        = local_penetrance_pair[2*threadindex+gamete]>=logpenetrance_threshold?
        exp(local_penetrance_pair[2*threadindex+gamete]):0;
      }
    }
  }
  return;
}

__kernel void impute_genotype_guide(
__const int max_haplotypes,
__const float epsilon,
__const int penetrance_matrix_size,
__const int packedextendedhap_len,
__const int flanking_snps,
__constant int * genotype_imputation,
__constant int * haplotypes,
__constant int * extended_haplotypes,
__constant int * last_site,
__global int * extended_root_mapping,
__global float * extended_frequency,
__global float * penetrance_cache,
__global int * subject_genotype_block,
__global float * subject_dosage_block,
__global float * subject_posterior_prob_block,
__global packedhap_t * packedextendedhap,
__global float * subject_haplotype_weight,
__local int * local_root_mapping,
__local float * local_frequency,
__local packedhap_t * local_packedextendedhap, //ref_haplo*packedextendedhaplen
__local float * local_posterior_prob0,
__local float * local_posterior_prob1,
__local float * local_posterior_prob2,
__local float * local_posterior_prob3,
__local float * local_posterior_prob, // FLANKING_SNPS * 4
__local float * local_cached_marginals,
__local float * local_denom
) {
  int subject = get_group_id(0);
  int packedsite = get_group_id(1);
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(extended_haplotypes[0]/BLOCK_WIDTH_IMPUTE_GUIDE)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH_IMPUTE_GUIDE+threadindex;
    if (hapindex<extended_haplotypes[0]){
      local_root_mapping[hapindex] = extended_root_mapping[hapindex];
      local_frequency[hapindex] = extended_frequency[hapindex];
    }
  }
  local_posterior_prob0[threadindex]=0;
  local_posterior_prob1[threadindex]=0;
  local_posterior_prob2[threadindex]=0;
  local_posterior_prob3[threadindex]=0;
  if (threadindex==0) local_denom[threadindex] = 0;
  // COPY HAPLOTYPE GROUP VERBATIM
  for(int haploindex = 0;haploindex<extended_haplotypes[0];++haploindex){
    if(threadindex<packedextendedhap_len){
      local_packedextendedhap[haploindex*packedextendedhap_len+threadindex] = 
      packedextendedhap[haploindex*packedextendedhap_len+threadindex];
    }
  }
  for(int chunk=0;chunk<(haplotypes[0]/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<haplotypes[0]){
      local_cached_marginals[hapindex] = 
      subject_haplotype_weight[subject*max_haplotypes+hapindex];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int righthapindex=0;righthapindex<extended_haplotypes[0];++righthapindex){
    if (local_cached_marginals[local_root_mapping[righthapindex]]>0){
      for(int chunk=0;chunk<(extended_haplotypes[0]/BLOCK_WIDTH_IMPUTE_GUIDE)+1;++chunk){
        int lefthapindex = chunk*BLOCK_WIDTH_IMPUTE_GUIDE+threadindex;
        //if (lefthapindex<=righthapindex ){
        if (lefthapindex<=righthapindex && local_cached_marginals[local_root_mapping[lefthapindex]]>0){
          // penetrance only needs to be fetched once as it doesn't change 
          // across sites
          //float penetrance = .5;
          float penetrance = penetrance_cache[subject*penetrance_matrix_size+
          local_root_mapping[righthapindex]*max_haplotypes+local_root_mapping[lefthapindex]];
          if (penetrance>0){
            // scale by the haplotype frequencies
            float prob = (lefthapindex!=righthapindex)?
            local_frequency[lefthapindex]*
            local_frequency[righthapindex]*penetrance*2:
            local_frequency[lefthapindex]*
            local_frequency[righthapindex]*penetrance;
            // get the dosage at packedsite
            int m = genotype_imputation[0]?(((int)local_packedextendedhap[righthapindex*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) + (((int)local_packedextendedhap[lefthapindex*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) 
            :2*(((int)local_packedextendedhap[righthapindex*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) + (((int)local_packedextendedhap[lefthapindex*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1); 
            // store the prob in the appropriate array
            local_posterior_prob0[threadindex]+=(m==0)*prob;
            local_posterior_prob1[threadindex]+=(m==1)*prob;
            local_posterior_prob2[threadindex]+=(m==2)*prob;
            local_posterior_prob3[threadindex]+=(m==3)*prob;
          }
        }
      } 
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  //return;
  // reduce
  for(int s=BLOCK_WIDTH_IMPUTE_GUIDE/2; s>0; s>>=1) {
    if (threadindex<s) {
       local_posterior_prob0[threadindex]+=
       local_posterior_prob0[threadindex+s];
       local_posterior_prob1[threadindex]+=
       local_posterior_prob1[threadindex+s];
       local_posterior_prob2[threadindex]+=
       local_posterior_prob2[threadindex+s];
       local_posterior_prob3[threadindex]+=
       local_posterior_prob3[threadindex+s];
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
  if (threadindex==0){
    float max_prob = 0;
    int geno2 = 0;
    for(int h=0;h<4;++h){
      if(local_posterior_prob[h]>max_prob){
        max_prob = local_posterior_prob[h];
        geno2 = h;
      }
      local_denom[threadindex]+=local_posterior_prob[h];
    }
    subject_genotype_block[subject*flanking_snps+packedsite] = geno2;
    if (genotype_imputation[0]){
      subject_dosage_block[subject*flanking_snps+packedsite] = (local_posterior_prob[1]+2*local_posterior_prob[2])/local_denom[threadindex];
    }
  }
  if (threadindex<4){
    subject_posterior_prob_block[subject*flanking_snps*4+packedsite*4+threadindex] = local_posterior_prob[threadindex]/local_denom[0];
  }
  return;
}

__kernel void impute_genotype_guide_haploid(
__const int max_haplotypes,
__const float epsilon,
__const int penetrance_matrix_size,
__const int packedextendedhap_len,
__const int flanking_snps,
__constant int * genotype_imputation,
__constant int * haplotypes,
__constant int * extended_haplotypes,
__constant int * last_site,
__global int * extended_root_mapping,
__global float * extended_frequency,
__global float * penetrance_cache,
__global int * subject_genotype_block,
__global float * subject_dosage_block,
__global float * subject_posterior_prob_block,
__global packedhap_t * packedextendedhap,
__global float * subject_haplotype_weight,
__local int * local_root_mapping,
__local float * local_frequency,
__local packedhap_t * local_packedextendedhap, //ref_haplo*packedextendedhaplen
__local float * local_posterior_prob0,
__local float * local_posterior_prob1,
__local float * local_posterior_prob2,
__local float * local_posterior_prob3,
__local float * local_posterior_prob, // FLANKING_SNPS * 4
__local float * local_cached_marginals,
__local float * local_denom
) {
  int subject = get_group_id(0);
  int packedsite = get_group_id(1);
  int threadindex = get_local_id(0);
  for(int chunk=0;chunk<(extended_haplotypes[0]/BLOCK_WIDTH_IMPUTE_GUIDE)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH_IMPUTE_GUIDE+threadindex;
    if (hapindex<extended_haplotypes[0]){
      local_root_mapping[hapindex] = extended_root_mapping[hapindex];
      local_frequency[hapindex] = extended_frequency[hapindex];
    }
  }
  local_posterior_prob0[threadindex]=0;
  local_posterior_prob1[threadindex]=0;
  local_posterior_prob2[threadindex]=0;
  local_posterior_prob3[threadindex]=0;
  if (threadindex==0) local_denom[threadindex] = 0;
  // COPY HAPLOTYPE GROUP VERBATIM
  for(int haploindex = 0;haploindex<extended_haplotypes[0];++haploindex){
    if(threadindex<packedextendedhap_len){
      local_packedextendedhap[haploindex*packedextendedhap_len+threadindex] = 
      packedextendedhap[haploindex*packedextendedhap_len+threadindex];
    }
  }
  for(int chunk=0;chunk<(haplotypes[0]/BLOCK_WIDTH)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH+threadindex;
    if (hapindex<haplotypes[0]){
      local_cached_marginals[hapindex] = 
      subject_haplotype_weight[subject*max_haplotypes+hapindex];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  //for(int righthapindex=0;righthapindex<extended_haplotypes[0];++righthapindex){
  //if (local_cached_marginals[local_root_mapping[righthapindex]]>0){
  for(int chunk=0;chunk<(extended_haplotypes[0]/BLOCK_WIDTH_IMPUTE_GUIDE)+1;++chunk){
    int hapindex = chunk*BLOCK_WIDTH_IMPUTE_GUIDE+threadindex;
      //if (lefthapindex<=righthapindex ){
    if (hapindex<extended_haplotypes[0] && 
    local_cached_marginals[local_root_mapping[hapindex]]>0){
      for(int gamete=0;gamete<2;++gamete){
        // penetrance only needs to be fetched once as it doesn't change 
        // across sites
        //float penetrance = .5;
        float penetrance = penetrance_cache[subject*max_haplotypes*2+
        local_root_mapping[hapindex]*2+gamete];
        if (penetrance>0){
          // scale by the haplotype frequencies
          float prob = local_frequency[hapindex]*penetrance;
          // get the dosage at packedsite
          int m = 2 * gamete + (((int)local_packedextendedhap[hapindex*
          packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) 
          >> ((packedsite%32)%8) & 1) ;
          // store the prob in the appropriate array
          local_posterior_prob0[threadindex]+=(m==0)*prob;
          local_posterior_prob1[threadindex]+=(m==1)*prob;
          local_posterior_prob2[threadindex]+=(m==2)*prob;
          local_posterior_prob3[threadindex]+=(m==3)*prob;
        }
      } 
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  //return;
  // reduce
  for(int s=BLOCK_WIDTH_IMPUTE_GUIDE/2; s>0; s>>=1) {
    if (threadindex<s) {
       local_posterior_prob0[threadindex]+=
       local_posterior_prob0[threadindex+s];
       local_posterior_prob1[threadindex]+=
       local_posterior_prob1[threadindex+s];
       local_posterior_prob2[threadindex]+=
       local_posterior_prob2[threadindex+s];
       local_posterior_prob3[threadindex]+=
       local_posterior_prob3[threadindex+s];
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
  if (threadindex<4){
    subject_posterior_prob_block[subject*flanking_snps*4+packedsite*4+threadindex] = local_posterior_prob[threadindex];
  }
  return;
}
