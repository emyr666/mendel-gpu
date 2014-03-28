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
    if(mendelgpu->g_markers==mendelgpu->g_snps){
      //int depth = parser->get_total_depth();
      //cerr<<"TOTAL_DEPTH:\t"<<depth<<endl;
    }
    // GPU friendly implementation
    // get total depth
    int non_zero = 0;
    int vec_index = 0;
    for(uint matrow=0;matrow<parser->matrix_rows;++matrow){
      int offset = parser->offset[matrow];
      //if(debug) cerr<<"Working on main read "<<matrow<<" with offset "<<offset<<endl;
      vector_offsets[subject*compact_rows+matrow] = vec_index;
      haplotype_offsets[subject*compact_rows+matrow] = 0;
      //if(debug) cerr<<"vector offset at row "<<matrow<<" is "<<vec_index<<endl;
      int k=0;
      bool hap_found = false;
      for(uint j=0;j<mendelgpu->g_max_window;++j){
        //read_alleles_mat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+j] = parser->matrix_alleles[matrow*MAX_MATRIX_WIDTH+j];
        if (parser->matrix_alleles[matrow*MAX_MATRIX_WIDTH+j]!=9) {
          if(vec_index<max_vector_length){
            read_alleles_vec[subject*max_vector_length+vec_index] =  parser->matrix_alleles[matrow*MAX_MATRIX_WIDTH+j];
            read_match_logvec[subject*max_vector_length+vec_index] =  parser->logmatch_quality[matrow][k];
            read_mismatch_logvec[subject*max_vector_length+vec_index] =  parser->logmismatch_quality[matrow][k];
            //if(debug) cerr<< "at vecindex "<<vec_index<<" alleles,match,mismatch are "<<read_alleles_vec[subject*max_vector_length+vec_index]<<","<<read_match_logvec[subject*max_vector_length+vec_index]<<","<<read_mismatch_logvec[subject*max_vector_length+vec_index]<<endl;
          }
          ++vec_index;
          ++k;
          ++non_zero;
          hap_found = true;
        }else{
          if(!hap_found) ++haplotype_offsets[subject*compact_rows+matrow];
        }
      }
      read_lengths[subject*compact_rows+matrow] = vec_index-vector_offsets[subject*compact_rows+matrow];
      //if(debug) cerr<<"Hap offset at row "<<matrow<<" is "<<haplotype_offsets[subject*compact_rows+matrow]<<endl;
      //if(debug) cerr<<"Read length at row "<<matrow<<" is "<<read_lengths[subject*compact_rows+matrow]<<endl;
      // first get the depth of main read
      int depth = parser->depth[matrow];
        for(uint baseindex=0;baseindex<parser->read_len[matrow];
        ++baseindex){
          if (offset+baseindex<mendelgpu->g_max_window){
            //read_match_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+offset+baseindex] = parser->logmatch_quality[matrow][baseindex];
            //read_mismatch_logmat[subject*read_compact_matrix_size+matrow*mendelgpu->g_max_window+offset+baseindex] = parser->logmismatch_quality[matrow][baseindex];
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
    if (non_zero>max_nonzero_elements) max_nonzero_elements = non_zero;
    mat_rows_by_subject[subject] = parser->matrix_rows;
  }
  parser->finalize(mendelgpu->g_left_marker);
  cerr<<"Elapsed time for CPU fetch read data on "<<mendelgpu->g_markers<<" markers: "<<(clock()-start)/CLOCKS_PER_SEC<<"\n";
  return;
}

float ReadPenetrance::get_bam_loglikelihood(int len,float *  hap1logprob){
  float logh1 = 0;
  for(int i=0;i<len;++i){
    logh1+=hap1logprob[i]!=0?(hap1logprob[i]):0;
  }
  return logh1;
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
  //assert(!isnan(logh1) && !isinf(logh1));
  //assert(!isnan(logh2) && !isinf(logh2));
  if (isnan(loglike) || isinf(loglike)){
    cerr<<"Warning: bamlikelihood failed with overflow: logh1,logh2: "<<logh1<<","<<logh2<<endl;
    // just choose the bigger of the two
    loglike = logh1;
    if (logh2>logh1) loglike = logh2; 
    cerr<<"Assuming loglike of "<<loglike<<endl;
    //throw "Overflow for BAM likelihood";
  }
  //assert(!isnan(loglike) && !isinf(loglike));
  return loglike;
}

struct order_hap_t{
  float logprob;
  float prob;
  int index1,index2;
  order_hap_t(float l, int i){
    logprob = l;
    index1 = i;
  }
  order_hap_t(float l, float p,int i1,int i2){
    logprob = l;
    prob = p;
    index1 = i1;
    index2 = i2;
  }
};
struct byprob{
  bool operator()(const order_hap_t & a,const order_hap_t & b){
return a.logprob>b.logprob;
  }
};

typedef multiset<order_hap_t,byprob> order_hap_set_t;

void ReadPenetrance::find_best_haplotypes(){
  cerr<<"Beginning find best haplotypes\n";
  double start = clock();
  for(int subject=0;subject<mendelgpu->g_people;++subject){
    order_hap_set_t ohs;
    bool debug = true;
    int total_depth = mat_rows_by_subject[subject];
    //if(debug)cerr<<"Total depth is "<<total_depth<<"\n";
    float max_log_pen = -1e10;
    float haploid_logprob[mendelgpu->g_max_haplotypes];
    for(int hap1=0;hap1<mendelgpu->g_max_haplotypes;++hap1){
      best_haplotypes[subject*mendelgpu->g_max_haplotypes+hap1] = 0;
      if (mendelgpu->g_active_haplotype[hap1]){
        // SCREENING STEP
        haploid_logprob[hap1] = 0;
        for(int row=0;row<total_depth;++row){
          int vector_offset = vector_offsets[subject*compact_rows+row];
          int haplotype_offset = haplotype_offsets[subject*compact_rows+row];
          int read_length = read_lengths[subject*compact_rows+row];
          for(int j=0;j<read_length;++j){
            if(haplotype_offset+j<mendelgpu->g_markers){
              if(vector_offset+j < max_vector_length){
                int read_allele = read_alleles_vec[subject*max_vector_length+vector_offset+j];
                float match_logprob = read_match_logvec[subject*max_vector_length+vector_offset+j];
                float mismatch_logprob = read_mismatch_logvec[subject*max_vector_length+vector_offset+j];
                haploid_logprob[hap1] += mendelgpu->g_haplotype[hap1*mendelgpu->g_max_window+haplotype_offset+j]==read_allele?match_logprob:mismatch_logprob;
              }
            }
          }
        }
        order_hap_t oh(haploid_logprob[hap1],hap1);
        ohs.insert(oh);
        mendelgpu->subject_haplotype_screen[subject*mendelgpu->g_max_haplotypes+hap1] = haploid_logprob[hap1];
        // END SCREEN
      }
    }
    int debug_penmat_person_start = 0;
    int debug_penmat_person_end = 0;
    int debug_penmat_person = TEST_SUBJECT;
    bool debug_screen = subject==debug_penmat_person;
    //bool debug_screen = subject>=debug_penmat_person_start && subject<=debug_penmat_person_end;
    int i = 0;
    if(debug_screen)  cerr<<"SCREEN FOR SUBJECT "<<subject<<endl;
    for(order_hap_set_t::iterator it = ohs.begin();it!=ohs.end();it++){
      order_hap_t oh = *it;
      if(i<mendelgpu->g_total_best_haps){
        if(debug_screen) cerr<<oh.index1<<": "<<oh.logprob<<endl;
        best_haplotypes[subject*mendelgpu->g_max_haplotypes+oh.index1] = 1;
      }
      ++i;
    }
  }
  cerr<<"Elapsed time for CPU screen haplotypes: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
}


void ReadPenetrance::process_read_matrices(){
  mendelgpu->g_left_flanking_snps = mendelgpu->g_flanking_snps;
  if (mendelgpu->g_center_snp_start<mendelgpu->g_flanking_snps){
    mendelgpu->g_left_flanking_snps = mendelgpu->g_center_snp_start;
  }
  find_best_haplotypes();
  int penetrance_matrix_size = mendelgpu->g_max_haplotypes*mendelgpu->g_max_haplotypes;
  if (mendelgpu->run_gpu){
    process_read_matrices_opencl();
  }
  if (mendelgpu->run_cpu){
    bool debug_penmat = true;
    int debug_penmat_person = TEST_SUBJECT;
    double start = clock();
    float max_log_pen_arr[mendelgpu->g_people];
    for(int i=0;i<mendelgpu->g_people;++i) max_log_pen_arr[i] = -1e10;
    //float min_log_pen = 0;
    for(int subject=0;subject<mendelgpu->g_people;++subject){
      bool debug = subject==TEST_SUBJECT;
      int total_depth = mat_rows_by_subject[subject];
      if(debug)cerr<<"Total depth is "<<total_depth<<"\n";
      bool old_debug = debug;
      for(int hap1=0;hap1<mendelgpu->g_max_haplotypes;++hap1){
        for(int hap2=hap1;hap2<mendelgpu->g_max_haplotypes;++hap2){
          if (mendelgpu->g_active_haplotype[hap1] && best_haplotypes[subject*mendelgpu->g_max_haplotypes+hap1] && mendelgpu->g_active_haplotype[hap2]&& best_haplotypes[subject*mendelgpu->g_max_haplotypes+hap2]){
            //debug = (hap1==2||hap2==11||hap2==2||hap2==11);
            if(debug) cerr<<"Computing reads penetrance for subject "<<subject<<" haps "<<hap1<<","<<hap2<<endl;
            if(debug) cerr<<hap1<<":";
            for(int i=0;i<mendelgpu->g_markers;++i) if(debug) cerr<<mendelgpu->g_haplotype[hap1*mendelgpu->g_max_window+i];
            if(debug) cerr<<endl;
            if(debug) cerr<<hap2<<":";
            for(int i=0;i<mendelgpu->g_markers;++i) if(debug) cerr<<mendelgpu->g_haplotype[hap2*mendelgpu->g_max_window+i];
            if(debug) cerr<<endl;
            // GPU friendly implementation
            // get total depth
            //explore the entire matrix 
            float log_like = 0;
            if (debug){
              //cerr<<"read_match_vec:";
              for(int i=0;i<max_vector_length;++i){
                //cerr<<" "<<read_mismatch_logvec[subject*max_vector_length+i];
              }
              //cerr<<endl;
            }
            for(int row=0;row<total_depth;++row){
              float hap1_logprob_new[mendelgpu->g_max_window];
              float hap2_logprob_new[mendelgpu->g_max_window];
                int vector_offset = vector_offsets[subject*compact_rows+row];
              int haplotype_offset = haplotype_offsets[subject*compact_rows+row];
              int read_length = read_lengths[subject*compact_rows+row];
              if(read_length>0 && vector_offset<max_vector_length){
                if(debug) {
                  //cerr<<"For read "<<row<<": "<<endl;
                  //cerr<<" Hap offset is "<<haplotype_offset<<endl;
                  //cerr<<" Vec offset is "<<vector_offset<<endl;
                  //cerr<<" Readlen "<<read_length<<endl;
                }
                // BEGIN GPU STYLE ALGORITHM
                int local_temp_read_alleles[mendelgpu->g_max_window];
                float local_temp_match[mendelgpu->g_max_window];
                float local_temp_mismatch[mendelgpu->g_max_window];
                for(int col=0;col<mendelgpu->g_max_window;++col){
                  hap1_logprob_new[col] = 0;
                  hap2_logprob_new[col] = 0;
                  local_temp_read_alleles[col] = 9;
                  local_temp_match[col] = 1;
                  local_temp_mismatch[col] = 1;
                }
                for(int col=0;col<read_length;++col){
                  if (haplotype_offset+col<mendelgpu->g_max_window && vector_offset+col<max_vector_length){
                    local_temp_read_alleles[haplotype_offset+col] = read_alleles_vec[subject*max_vector_length+vector_offset+col];
                    local_temp_match[haplotype_offset+col] = read_match_logvec[subject*max_vector_length+vector_offset+col];
                    local_temp_mismatch[haplotype_offset+col] = read_mismatch_logvec[subject*max_vector_length+vector_offset+col];
                  }
                }
                //if (debug) cerr<<"local_temp_read_alleles: ";
                //for(int col=0;col<mendelgpu->g_max_window;++col){
                //  if(debug) cerr<<local_temp_read_alleles[col];
               // }
               // if (debug) cerr<<endl;
               // if (debug) cerr<<"local_temp_match,mismatch: ";
               // for(int col=0;col<mendelgpu->g_max_window;++col){
               //   if(debug) cerr<<" "<<local_temp_match[col]<<","<<local_temp_mismatch[col];
               // }
               // if (debug) cerr<<endl;
                //for(int col=0;col<read_length;++col){
                for(int col=0;col<mendelgpu->g_markers;++col){
                  int read_allele = local_temp_read_alleles[col];
                  float match_logprob = local_temp_match[col];
                  float mismatch_logprob = local_temp_mismatch[col];
                  if(read_allele!=9){
                    hap1_logprob_new[col] = mendelgpu->g_haplotype[hap1*mendelgpu->g_max_window+col]==read_allele?match_logprob:mismatch_logprob;
                    hap2_logprob_new[col] = mendelgpu->g_haplotype[hap2*mendelgpu->g_max_window+col]==read_allele?match_logprob:mismatch_logprob;
                    //if(debug)cerr<<"NEW hap1prob at "<<col<<" is "<<hap1_logprob_new[col]<<endl;
                    //if(debug)cerr<<"NEW hap2prob at "<<col<<" is "<<hap2_logprob_new[col]<<endl;
                  }

                }
                // END GPU STYLE ALGORITHM
                float ll2 = get_bam_loglikelihood(mendelgpu->g_markers,hap1_logprob_new,hap2_logprob_new);
                //if (debug) cerr<<"current log like "<<ll2<<endl; 
                log_like+=ll2;
              }
            }
            if(debug ) cerr<<"Log likelihood of all reads for haps "<<hap1<<","<<hap2<<" is "<<log_like<<endl;
            if (log_like>max_log_pen_arr[subject]) max_log_pen_arr[subject] = log_like;
            //if (log_like<min_log_pen) min_log_pen = log_like;
            mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap1*mendelgpu->g_max_haplotypes+hap2] = log_like;
          }else{
            mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap1*mendelgpu->g_max_haplotypes+hap2] = -999;
          } // active hap condition
        } // second hap loop
      }// first hap loop
      debug = old_debug;
    } // loop subjects
    float max_log_pen = -1e10;
    for(int subject=0;subject<mendelgpu->g_people;++subject){
      // get the least upper bound on the log penetrances
      if (max_log_pen_arr[subject]>max_log_pen)max_log_pen = max_log_pen_arr[subject];
    }
    for(int subject=0;subject<mendelgpu->g_people;++subject){
      order_hap_set_t ohs_unphased;
      //max_log_pen = 0;
      cerr<<"Max log penetrance across subjects is "<<max_log_pen<<endl;
      if(mendelgpu->g_flanking_snps){
        cerr<<"Subject "<<subject<<" Contents of g_left_hap_imputed:\n";
        for(int j=0;j<mendelgpu->g_left_flanking_snps;++j){
          cerr<<mendelgpu->g_left_hap_imputed[(subject*2)*mendelgpu->g_flanking_snps+j];
        }
        cerr<<endl;
        for(int j=0;j<mendelgpu->g_left_flanking_snps;++j){
          cerr<<mendelgpu->g_left_hap_imputed[(subject*2+1)*mendelgpu->g_flanking_snps+j];
        }
        cerr<<endl;
      }
      for(int hap1=0;hap1<mendelgpu->g_max_haplotypes;++hap1){
        for(int hap2=hap1;hap2<mendelgpu->g_max_haplotypes;++hap2){
          if (mendelgpu->g_active_haplotype[hap1] && mendelgpu->g_active_haplotype[hap2]){
            //cerr<<"Subject "<<subject<<" Contents of current hap1 hap2:\n";
            //mendelgpu->debug_haplotype(cerr,hap1);
            //mendelgpu->debug_haplotype(cerr,hap2);
            float unscaledval = mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap1* mendelgpu->g_max_haplotypes+hap2] ;
            //mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap1*
            //mendelgpu->g_max_haplotypes+hap2] = unscaledval; 
            mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap2*
            mendelgpu->g_max_haplotypes+hap1] = unscaledval;
            float global_scaledval = unscaledval - max_log_pen;
            float scaledval = unscaledval - max_log_pen_arr[subject];
//cerr<<"MAXLOGPEN: "<<max_log_pen<<" MAXLOGPENARR[0] "<<max_log_pen_arr[subject]<<endl;

            mendelgpu->assign_penetrance_element(subject,hap1,hap2,global_scaledval,mendelgpu->penetrance_cache, false );
            mendelgpu->assign_penetrance_element(subject,hap1,hap2,scaledval,mendelgpu->scaled_penetrance_cache, true);
            if (mendelgpu->scaled_penetrance_cache[subject*penetrance_matrix_size+hap1* mendelgpu->g_max_haplotypes+hap2]>0){
              order_hap_t oh(mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap1* mendelgpu->g_max_haplotypes+hap2],mendelgpu->scaled_penetrance_cache[subject*penetrance_matrix_size+hap1* mendelgpu->g_max_haplotypes+hap2],hap1,hap2);
              ohs_unphased.insert(oh);
            }
            if (mendelgpu->scaled_penetrance_cache[subject*penetrance_matrix_size+hap2* mendelgpu->g_max_haplotypes+hap1]>0){
              order_hap_t oh(mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+hap2* mendelgpu->g_max_haplotypes+hap1],mendelgpu->scaled_penetrance_cache[subject*penetrance_matrix_size+hap2* mendelgpu->g_max_haplotypes+hap1],hap2,hap1);
              ohs_unphased.insert(oh);
            }
          }else{
            //cerr<<"Setting hap 1,2: "<<hap1<<","<<hap2<<endl;
            mendelgpu->penetrance_cache[subject*penetrance_matrix_size+hap1* mendelgpu->g_max_haplotypes+hap2] = mendelgpu->penetrance_cache[subject*penetrance_matrix_size+hap2* mendelgpu->g_max_haplotypes+hap1] = 0; 
            mendelgpu->scaled_penetrance_cache[subject*penetrance_matrix_size+hap1* mendelgpu->g_max_haplotypes+hap2] = mendelgpu->scaled_penetrance_cache[subject*penetrance_matrix_size+hap2* mendelgpu->g_max_haplotypes+hap1] = 0; 
          } // active hap condition
        } // second hap loop
      }// first hap loop
      int debug_penmat_person_start = 0;
      int debug_penmat_person_end = mendelgpu->g_people-1;
      if(debug_penmat && (subject==debug_penmat_person)){
        if (mendelgpu->true_haps!=NULL){
          cerr<<"TRUE HAPS FOR SUBJECT "<<subject<<endl;
          for(int k=0;k<2;++k){
            for(int j=mendelgpu->g_left_marker;j<mendelgpu->g_left_marker+mendelgpu->g_markers;++j){
              cerr<<mendelgpu->true_haps[(2*subject+k)*mendelgpu->g_snps+j];
            }
            cerr<<endl;
          }
        }
        cerr<<"PEN RANKING FOR SUBJECT "<<subject<<" (LogLike, Scaled Penetrance)"<<endl;
        int k=0;
        for(order_hap_set_t::iterator it = ohs_unphased.begin();it!=ohs_unphased.end();it++){
          order_hap_t oh = *it;
          //if(k<10 && oh.logprob>mendelgpu->gf_logpen_threshold){
            cerr<<oh.index1<<","<<oh.index2<<": "<<oh.logprob<<","<<oh.prob<<endl;
            mendelgpu->debug_haplotype(cerr,oh.index1);
            mendelgpu->debug_haplotype(cerr,oh.index2);
            ++k;
          //}
        }
      }
// works here
//return;
      debug_penmat_person_start = TEST_SUBJECT;
      debug_penmat_person_end = TEST_SUBJECT;
      if(debug_penmat && (subject>=debug_penmat_person_start && subject<=debug_penmat_person_end)) {
        cerr<<"SUBJECT "<<subject<<endl;
        for(int j=0;j<mendelgpu->g_max_haplotypes;++j){
          if (mendelgpu->g_active_haplotype[j]){
            //cerr<<j<<":";
            for(int k=0;k<mendelgpu->g_max_haplotypes;++k){
              if (mendelgpu->g_active_haplotype[k]){
                if (mendelgpu->penetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k]>0){
                  //cerr<<" "<<j<<","<<k<<":"<<mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k];
                }
                cerr<<" "<<j<<","<<k<<":"<<mendelgpu->logpenetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k];
                cerr<<","<<mendelgpu->penetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k];
              //cerr<<" "<<mendelgpu->logmendelgpu->penetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k]<<","<<mendelgpu->penetrance_cache[subject*penetrance_matrix_size+j*mendelgpu->g_max_haplotypes+k];
              }
            }
            cerr<<endl;
          }
        }
      }
    }
    cerr<<"Exiting compute penetrance in "<<(clock()-start)/CLOCKS_PER_SEC<<"\n";
    //if(debug_penmat) exit(0);
  }
  return;
}

void ReadPenetrance::process_read_matrices_rnaseq(){
  mendelgpu->g_left_flanking_snps = mendelgpu->g_flanking_snps;
  if (mendelgpu->g_center_snp_start<mendelgpu->g_flanking_snps){
    mendelgpu->g_left_flanking_snps = mendelgpu->g_center_snp_start;
  }
  //find_best_haplotypes();
  int penetrance_matrix_size = mendelgpu->g_max_haplotypes*mendelgpu->g_max_haplotypes;
  if (mendelgpu->run_gpu){
    //process_read_matrices_opencl();
  }
  if (mendelgpu->run_cpu){
    bool debug_penmat = true;
    int debug_penmat_person = TEST_SUBJECT;
    double start = clock();
    float max_log_pen_arr[mendelgpu->g_people];
    for(int i=0;i<mendelgpu->g_people;++i) max_log_pen_arr[i] = -1e10;
    //float min_log_pen = 0;
    bool read_aligned = false;
    int * depth_by_hap = new int[mendelgpu->g_people * mendelgpu->g_max_haplotypes];
    for(int h=0;h<mendelgpu->g_people * mendelgpu->g_max_haplotypes;++h) depth_by_hap[h] = 0;
    for(int subject=0;subject<mendelgpu->g_people;++subject){
      bool debug = subject==TEST_SUBJECT;
      int total_depth = mat_rows_by_subject[subject];
      if(debug)cerr<<"Total depth is "<<total_depth<<"\n";
      bool old_debug = debug;
      for(int hap1=0;hap1<mendelgpu->g_max_haplotypes;++hap1){
        if (mendelgpu->g_active_haplotype[hap1]){
          //debug = (hap1==2||hap2==11||hap2==2||hap2==11);
          if(debug) cerr<<"Computing reads penetrance for subject "<<subject<<" hap "<<hap1<<endl;
          if(debug) cerr<<hap1<<":";
          for(int i=0;i<mendelgpu->g_markers;++i) if(debug) cerr<<mendelgpu->g_haplotype[hap1*mendelgpu->g_max_window+i];
          if(debug) cerr<<endl;
          //explore the entire matrix 
          //float log_like = 0;
          double hapfreq = 0;
          int max_read_col = -1;
          int min_read_col = 1000;
          int high_confidence_depth = 0;
          for(int row=0;row<total_depth;++row){
            float hap1_logprob_new[mendelgpu->g_max_window];
            int vector_offset = vector_offsets[subject*compact_rows+row];
            int haplotype_offset = haplotype_offsets[subject*compact_rows+row];
            int read_length = read_lengths[subject*compact_rows+row];
            if(read_length>0 && vector_offset<max_vector_length){
              int local_temp_read_alleles[mendelgpu->g_max_window];
              float local_temp_match[mendelgpu->g_max_window];
              float local_temp_mismatch[mendelgpu->g_max_window];
              for(int col=0;col<mendelgpu->g_max_window;++col){
                  hap1_logprob_new[col] = 0;
                  local_temp_read_alleles[col] = 9;
                  local_temp_match[col] = 1;
                  local_temp_mismatch[col] = 1;
              }
              for(int col=0;col<read_length;++col){
                if (haplotype_offset+col<mendelgpu->g_max_window && vector_offset+col<max_vector_length){
                  local_temp_read_alleles[haplotype_offset+col] = read_alleles_vec[subject*max_vector_length+vector_offset+col];
                  local_temp_match[haplotype_offset+col] = read_match_logvec[subject*max_vector_length+vector_offset+col];
                  local_temp_mismatch[haplotype_offset+col] = read_mismatch_logvec[subject*max_vector_length+vector_offset+col];
                }
              }
              // Josh requests including only reads exceeding read len 
              // specific thresholds
              float log_thresholds[]={0.0000000,-0.1053605,-0.2107210,-0.3160815,-0.4214421,-0.5268026};
              int curr_depth = parser->depth[row];
              for(int col=0;col<mendelgpu->g_markers;++col){
                int read_allele = local_temp_read_alleles[col];
                float match_logprob = local_temp_match[col];
                float mismatch_logprob = local_temp_mismatch[col];
                if(read_allele!=9 && read_length>1){
                  hap1_logprob_new[col] = mendelgpu->g_haplotype[hap1*mendelgpu->g_max_window+col]==read_allele?match_logprob:mismatch_logprob;
                  if (debug) cerr<<"Read col "<<col<<" logprob "<<hap1_logprob_new[col]<<endl;
                  if (col>max_read_col) max_read_col = col;
                  if (col<min_read_col) min_read_col = col;
                }
              }
              float ll2 = get_bam_loglikelihood(mendelgpu->g_markers,hap1_logprob_new);
              if(ll2>log_thresholds[read_length]*curr_depth && read_length>1){
                hapfreq+=exp(ll2);
                high_confidence_depth+=curr_depth;
              }
              //log_like+=ll2;
            }
          }
          depth_by_hap[subject*mendelgpu->g_max_haplotypes+hap1] = high_confidence_depth;
          read_aligned = (min_read_col==0 && max_read_col==(mendelgpu->g_markers-1));
          if(debug) cerr<<"Read col range is ("<<min_read_col<<"-"<<max_read_col<<"). Read aligned is "<<read_aligned<<"\n";
          if(debug) cerr<<"Frequency of all reads for hap "<<hap1<<" is "<<hapfreq<<endl;
          mendelgpu->subject_haplotype_screen[subject*mendelgpu->g_max_haplotypes+hap1] = hapfreq;
        } // active hap condition
      }// first hap loop
      debug = old_debug;
    } // loop subjects
    if (read_aligned){
      mendelgpu->io_manager->writeSubjectHapFreq(mendelgpu->g_people,mendelgpu->g_left_marker,mendelgpu->g_markers,mendelgpu->g_max_window,mendelgpu->g_max_haplotypes,mendelgpu->g_active_haplotype,mendelgpu->g_haplotype,mendelgpu->subject_haplotype_screen, depth_by_hap);
    }
    cerr<<"Exiting compute hap freq in "<<(clock()-start)/CLOCKS_PER_SEC<<"\n";
  }
  return;
}
