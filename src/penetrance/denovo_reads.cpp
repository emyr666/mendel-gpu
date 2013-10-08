#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_penetrance.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif


void DenovoReadsMendelGPU::compute_penetrance(){
  populate_read_matrices();
  process_read_matrices();
}

void DenovoReadsMendelGPU::populate_read_matrices(){
  memset(read_match_logmat,0.,sizeof(float)*g_people*read_full_matrix_size);
  memset(read_mismatch_logmat,0.,sizeof(float)*g_people*read_full_matrix_size);
  for(int subject=0;subject<g_people;++subject){
    bool debug = subject==-50;
    ReadParser & parser = parsers[subject];
    if(debug) cerr<<"Extracting region for subject "<<subject<<" for left marker "<<g_left_marker<<" of markers length "<<g_markers<<"\n";
    parser.extract_region(config->chromosome.data(),g_left_marker,g_markers);
    if(debug) cerr<<"Printing data\n";
    if(debug)parser.print_data();
    // GPU friendly implementation
    // get total depth
    int expanded_row = 0;
    for(uint matrow=0;matrow<parser.matrix_rows;++matrow){
      int offset = parser.offset[matrow];
      if(debug) cerr<<"Working on main read "<<matrow<<" with offset "<<offset<<endl;
      for(uint j=0;j<g_max_window;++j){
        read_alleles_mat[subject*read_compact_matrix_size+matrow*g_max_window+j] = parser.matrix_alleles[matrow*g_max_window+j];
      }
      // first get the depth of main read
      int depth = parser.depth[matrow];
      for(uint readindex=0;readindex<depth;++readindex){
        superread_indices[subject*full_rows+expanded_row] = matrow;
        for(uint baseindex=0;baseindex<parser.read_len[matrow];
        ++baseindex){
          // APPROACH 1: GPU friendly
          int phredscore = parser.base_quality[matrow][readindex][baseindex];
          //if(debug) cerr<<"Phred score is "<<phredscore<<endl;
          read_match_logmat[subject*read_full_matrix_size+expanded_row*g_max_window+offset+baseindex] = parser.logprob_match_lookup[phredscore];
          read_mismatch_logmat[subject*read_full_matrix_size+expanded_row*g_max_window+offset+baseindex] = parser.logprob_mismatch_lookup[phredscore];
          if (debug)cerr<<"POPULATE for row "<<expanded_row<<" col "<<offset+baseindex<<" match,mismatch: "<<read_match_logmat[subject*read_full_matrix_size+expanded_row*g_max_window+offset+baseindex]<<","<<read_mismatch_logmat[subject*read_full_matrix_size+expanded_row*g_max_window+offset+baseindex]<<endl;
        }
        ++expanded_row;
      }   
    } 
    mat_rows_by_subject[subject] = expanded_row;
  }
  return;
}

float DenovoReadsMendelGPU::get_bam_loglikelihood(int len,float *  hap1logprob,float *  hap2logprob){
  float logh1 = 0,logh2 = 0;
  for(int i=0;i<len;++i){
    //cerr<<"debug"<<i<<":"<<hap1logprob[i]<<","<<hap2logprob[i]<<endl;
    logh1+=hap1logprob[i]!=0?(hap1logprob[i]):0;
    logh2+=hap2logprob[i]!=0?(hap2logprob[i]):0;
  }
  float mean = (logh1+logh2)/2.;
  //cerr<<"Mean offset is "<<mean<<endl;
  float loglike = log_half + (log(exp(logh1-mean)+exp(logh2-mean))+mean);
  return loglike;
}

void DenovoReadsMendelGPU::process_read_matrices(){
  for(int subject=0;subject<g_people;++subject){
    bool debug = subject==-50;
    float max_log_pen = -1e10;
    for(int hap1=0;hap1<g_max_haplotypes;++hap1){
      for(int hap2=hap1;hap2<g_max_haplotypes;++hap2){
        if (g_active_haplotype[hap1] && g_active_haplotype[hap2]){
          if(debug) cerr<<"Computing reads penetrance for subject "<<subject<<" haps "<<hap1<<","<<hap2<<endl;
          if(debug) cerr<<hap1<<":";
          for(int i=0;i<g_markers;++i) if(debug) cerr<<g_haplotype[hap1*g_max_window+i];
          if(debug) cerr<<endl;
          if(debug) cerr<<hap2<<":";
          for(int i=0;i<g_markers;++i) if(debug) cerr<<g_haplotype[hap2*g_max_window+i];
          if(debug) cerr<<endl;
          // GPU friendly implementation
          // get total depth
          int total_depth = mat_rows_by_subject[subject];
          int total_width = g_max_window;
          //explore the entire matrix 
          float log_like = 0;
          for(int row=0;row<total_depth;++row){
            float hap1_logprob[g_max_window];
            float hap2_logprob[g_max_window];
            int superread_index = superread_indices[subject*full_rows+row];
            for(int col=0;col<total_width;++col){
              int allele = read_alleles_mat[subject*read_compact_matrix_size+superread_index*g_max_window+col] ;
              if(allele!=9 && col<g_markers){
                hap1_logprob[col] = g_haplotype[hap1*g_max_window+col]==allele?
                read_match_logmat[subject*read_full_matrix_size+row*g_max_window+col]:
                read_mismatch_logmat[subject*read_full_matrix_size+row*g_max_window+col];
                hap2_logprob[col] = g_haplotype[hap2*g_max_window+col]==allele?
                read_match_logmat[subject*read_full_matrix_size+row*g_max_window+col]:
                read_mismatch_logmat[subject*read_full_matrix_size+row*g_max_window+col];
                if (debug)cerr<<"for row,col "<<row<<","<<col<<" GOT HERE with vals "<<hap1_logprob[col]<<","<<hap2_logprob[col]<<"\n";
             
              }else{
                hap1_logprob[col] = 0;
                hap2_logprob[col] = 0;
              }
            }
            float ll = get_bam_loglikelihood(g_max_window,hap1_logprob,hap2_logprob);
            //if (debug) cerr<<"current log like "<<ll<<endl; 
            log_like+=ll;
          }
          if(debug) cerr<<"Log likelihood of all reads for haps "<<hap1<<","<<hap2<<" is "<<log_like<<endl;
          if (log_like>max_log_pen) max_log_pen = log_like;
          logpenetrance_cache[subject*penetrance_matrix_size+hap1*g_max_haplotypes+hap2] = log_like;
        } // active hap condition
      } // second hap loop
    }// first hap loop
    //if(debug) cerr<<"Max log penetrance for subject "<<subject<<" is "<<max_log_pen<<endl;
    for(int hap1=0;hap1<g_max_haplotypes;++hap1){
      for(int hap2=hap1;hap2<g_max_haplotypes;++hap2){
        if (g_active_haplotype[hap1] && g_active_haplotype[hap2]){
          float val = logpenetrance_cache[subject*penetrance_matrix_size+hap1*
          g_max_haplotypes+hap2]-max_log_pen;
          logpenetrance_cache[subject*penetrance_matrix_size+hap1*
          g_max_haplotypes+hap2] = val;
          logpenetrance_cache[subject*penetrance_matrix_size+hap2*
          g_max_haplotypes+hap1] = val;
          penetrance_cache[subject*penetrance_matrix_size+hap1*
          g_max_haplotypes+hap2] = val>=gf_logpen_threshold?exp(val):0;
          penetrance_cache[subject*penetrance_matrix_size+hap2*
          g_max_haplotypes+hap1] = penetrance_cache[subject*
          penetrance_matrix_size+hap1*g_max_haplotypes+hap2];
          //if(debug) cerr<<"Penetrance for subject "<<subject<<" haps "<<hap1<<","<<hap2<<": "<<penetrance_cache[subject*penetrance_matrix_size+hap1*g_max_haplotypes+hap2]<<endl;
        } // active hap condition
      } // second hap loop
    }// first hap loop
  }
  cerr<<"Exiting compute penetrance\n";
  return;
}


