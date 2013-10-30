#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif


void ReadPenetrance::prefetch_reads(int start_snp, int total_snps){
  // delegate call to the reads parser
  if (start_snp<mendelgpu->g_snps){
    double start = clock();
    int snps2fetch = start_snp+total_snps<=mendelgpu->g_snps?total_snps:
    mendelgpu->g_snps-start_snp;
    cerr<<"Prefetching reads for "<<snps2fetch<<" SNPs starting from SNP "<<start_snp<<endl;
    for(int subject_index=0;subject_index<mendelgpu->g_people;++subject_index){
      parser->extract_region(subject_index,start_snp,total_snps,false);
    }
    cerr<<"Prefetch time is "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
  }
}

void ReadPenetrance::populate_read_matrices(){
  Config * config = mendelgpu->config;
  double start = clock();
  cerr<<"Extracting reads at left marker "<<mendelgpu->g_left_marker<<" with region length "<<mendelgpu->g_markers<<"\n";
  for(int subject=0;subject<mendelgpu->g_people;++subject){
    bool debug = (subject==TEST_SUBJECT);
    parser->extract_region(subject,mendelgpu->g_left_marker,mendelgpu->g_markers,true);
    if(debug) cerr<<"Printing subject "<<subject<<"\n";
    if(debug)parser->print_data();
    // GPU friendly implementation
    // get total depth
    for(uint matrow=0;matrow<parser->matrix_rows;++matrow){
      int offset = parser->offset[matrow];
      //if(debug) cerr<<"Working on main read "<<matrow<<" with offset "<<offset<<endl;
      for(uint j=0;j<mendelgpu->g_max_window;++j){
        read_alleles_mat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+j] = parser->matrix_alleles[matrow*MAX_MATRIX_WIDTH+j];
        //read_match_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+j] = 0;
        //read_mismatch_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+j] = 0;
      }
      // first get the depth of main read
      int depth = parser->depth[matrow];
      //for(uint readindex=0;readindex<depth;++readindex){
        for(uint baseindex=0;baseindex<parser->read_len[matrow];
        ++baseindex){
          if (offset+baseindex<mendelgpu->g_max_window){
            //int phredscore = parser->base_quality[matrow][readindex][baseindex];
            //if(debug) cerr<<"Phred score is "<<phredscore<<endl;
            //read_match_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+offset+baseindex] += parser->logprob_match_lookup[phredscore];
            //read_mismatch_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+offset+baseindex] += parser->logprob_mismatch_lookup[phredscore];
            read_match_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+offset+baseindex] = parser->logmatch_quality[matrow][baseindex];
            read_mismatch_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+offset+baseindex] = parser->logmismatch_quality[matrow][baseindex];
          }
        }
      //}   
      for(uint baseindex=0;baseindex<parser->read_len[matrow];
        ++baseindex){
        if (offset+baseindex<mendelgpu->g_max_window){
          //if (debug)cerr<<"POPULATE for row "<<matrow<<" col "<<offset+baseindex<<" match,mismatch: "<<read_match_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+offset+baseindex]<<","<<read_mismatch_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+offset+baseindex]<<endl;
        }
      }
    } 
    mat_rows_by_subject[subject] = parser->matrix_rows;
  }
  parser->finalize(mendelgpu->g_left_marker);
  cerr<<"Elapsed time for CPU fetch read data on "<<mendelgpu->g_markers<<" markers: "<<(clock()-start)/CLOCKS_PER_SEC<<"\n";
  return;
}

float ReadPenetrance::get_bam_loglikelihood(int len,float *  hap1logprob,float *  hap2logprob){
  float logh1 = 0,logh2 = 0;
  for(int i=0;i<len;++i){
    //cerr<<"debug"<<i<<":"<<hap1logprob[i]<<","<<hap2logprob[i]<<endl;
    logh1+=hap1logprob[i]!=0?(hap1logprob[i]):0;
    logh2+=hap2logprob[i]!=0?(hap2logprob[i]):0;
  }
  float mean = (logh1+logh2)/2.;
  //cerr<<"Mean offset is "<<mean<<endl;
  float loglike = LOG_HALF + (log(exp(logh1-mean)+exp(logh2-mean))+mean);
  return loglike;
}


void ReadPenetrance::process_read_matrices(){
  int penetrance_matrix_size = mendelgpu->g_max_haplotypes*mendelgpu->g_max_haplotypes;
  if (mendelgpu->run_gpu){
    process_read_matrices_opencl();
  }
  if (mendelgpu->run_cpu){
    bool debug_penmat = false;
    int debug_penmat_person = 11;
    double start = clock();
    for(int subject=0;subject<mendelgpu->g_people;++subject){
      bool debug = subject==-10||subject==-11;
      float max_log_pen = -1e10;
      for(int hap1=0;hap1<mendelgpu->g_max_haplotypes;++hap1){
        for(int hap2=hap1;hap2<mendelgpu->g_max_haplotypes;++hap2){
          if (mendelgpu->g_active_haplotype[hap1] && mendelgpu->g_active_haplotype[hap2]){
            if(debug) cerr<<"Computing reads penetrance for subject "<<subject<<" haps "<<hap1<<","<<hap2<<endl;
            if(debug) cerr<<hap1<<":";
            for(int i=0;i<mendelgpu->g_markers;++i) if(debug) cerr<<mendelgpu->g_haplotype[hap1*mendelgpu->g_max_window+i];
            if(debug) cerr<<endl;
            if(debug) cerr<<hap2<<":";
            for(int i=0;i<mendelgpu->g_markers;++i) if(debug) cerr<<mendelgpu->g_haplotype[hap2*mendelgpu->g_max_window+i];
            if(debug) cerr<<endl;
            // GPU friendly implementation
            // get total depth
            int total_depth = mat_rows_by_subject[subject];
            int total_width = mendelgpu->g_max_window;
            if(debug)cerr<<"Total depth is "<<total_depth<<"\n";
            //explore the entire matrix 
            float log_like = 0;
            for(int row=0;row<total_depth;++row){
              float hap1_logprob[mendelgpu->g_max_window];
              float hap2_logprob[mendelgpu->g_max_window];
              for(int col=0;col<total_width;++col){
                int allele = read_alleles_mat[subject*read_compact_matrix_size+row*mendelgpu->g_max_window+col] ;
                //if (debug) cerr<<"Allele is "<<allele<<endl;
                if(allele!=9 && col<mendelgpu->g_markers){
                  hap1_logprob[col] = mendelgpu->g_haplotype[hap1*mendelgpu->g_max_window+col]==allele?
                  read_match_logmat[subject*read_compact_matrix_size+row*mendelgpu->g_max_window+col]:
                  read_mismatch_logmat[subject*read_compact_matrix_size+row*mendelgpu->g_max_window+col];
                  hap2_logprob[col] = mendelgpu->g_haplotype[hap2*mendelgpu->g_max_window+col]==allele?
                  read_match_logmat[subject*read_compact_matrix_size+row*mendelgpu->g_max_window+col]:
                  read_mismatch_logmat[subject*read_compact_matrix_size+row*mendelgpu->g_max_window+col];
                  if (debug)cerr<<"for row,col "<<row<<","<<col<<" GOT HERE with vals "<<hap1_logprob[col]<<","<<hap2_logprob[col]<<"\n";
               
                }else{
                  hap1_logprob[col] = 0;
                  hap2_logprob[col] = 0;
                }
              }
              float ll = get_bam_loglikelihood(mendelgpu->g_max_window,hap1_logprob,hap2_logprob);
              //if (debug) cerr<<"current log like "<<ll<<endl; 
              log_like+=ll;
            }
            if(debug) cerr<<"Log likelihood of all reads for haps "<<hap1<<","<<hap2<<" is "<<log_like<<endl;
            if (log_like>max_log_pen) max_log_pen = log_like;
            mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap1*mendelgpu->g_max_haplotypes+hap2] = log_like;
          } // active hap condition
        } // second hap loop
      }// first hap loop
      //if(debug) cerr<<"Max log penetrance for subject "<<subject<<" is "<<max_log_pen<<endl;
      for(int hap1=0;hap1<mendelgpu->g_max_haplotypes;++hap1){
        for(int hap2=hap1;hap2<mendelgpu->g_max_haplotypes;++hap2){
          if (mendelgpu->g_active_haplotype[hap1] && mendelgpu->g_active_haplotype[hap2]){
            
            float val = mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap1*
            mendelgpu->g_max_haplotypes+hap2]-max_log_pen;
            mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap1*
            mendelgpu->g_max_haplotypes+hap2] = val;
            mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap2*
            mendelgpu->g_max_haplotypes+hap1] = val;
            mendelgpu->penetrance_cache[subject*penetrance_matrix_size+hap1*
            mendelgpu->g_max_haplotypes+hap2] = val>=mendelgpu->gf_logpen_threshold?exp(val):0;
            mendelgpu->penetrance_cache[subject*penetrance_matrix_size+hap2*
            mendelgpu->g_max_haplotypes+hap1] = mendelgpu->penetrance_cache[subject*
            penetrance_matrix_size+hap1*mendelgpu->g_max_haplotypes+hap2];
          } // active hap condition
        } // second hap loop
      }// first hap loop
      if(debug_penmat && (subject==debug_penmat_person)) {
        cerr<<"SUBJECT "<<subject<<endl;
        for(int j=0;j<mendelgpu->g_max_haplotypes;++j){
          if (mendelgpu->g_active_haplotype[j]){
            cerr<<j<<":";
            for(int k=0;k<mendelgpu->g_max_haplotypes;++k){
              if (mendelgpu->g_active_haplotype[k]){
                cerr<<" "<<mendelgpu->penetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k];
                //cerr<<" "<<mendelgpu->logmendelgpu->penetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k]<<","<<mendelgpu->penetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k];
              }
            }
            cerr<<endl;
          }
        }
      }
    }
    cerr<<"Exiting compute penetrance in "<<(clock()-start)/CLOCKS_PER_SEC<<"\n";
    if(debug_penmat) exit(0);
  }
  return;
}
  
  
