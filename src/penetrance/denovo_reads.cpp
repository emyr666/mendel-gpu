#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_penetrance.hpp"
#ifdef USE_GPU
#include"../clsafe.h"
#endif

float DenovoReadsMendelGPU::get_bam_loglikelihood(int len,float *  hap1prob,float *  hap2prob){
  float logh1 = 0,logh2 = 0;
  for(int i=0;i<len;++i){
    //cerr<<"debug"<<i<<":"<<hap1prob[i]<<","<<hap2prob[i]<<endl;
    logh1+=hap1prob[i]!=1?log(hap1prob[i]):0;
    logh2+=hap2prob[i]!=1?log(hap2prob[i]):0;
  }
  float mean = (logh1+logh2)/2.;
  //cerr<<"Mean offset is "<<mean<<endl;
  float loglike = log_half + (log(exp(logh1-mean)+exp(logh2-mean))+mean);
  return loglike;
}

void DenovoReadsMendelGPU::compute_penetrance(){
  for(int subject=0;subject<g_people;++subject){
    float max_log_pen = -1e10;
    for(int hap1=0;hap1<g_max_haplotypes;++hap1){
      for(int hap2=hap1;hap2<g_max_haplotypes;++hap2){
        if (g_active_haplotype[hap1] && g_active_haplotype[hap2]){
          cerr<<"Computing reads penetrance for subject "<<subject<<" haps "<<hap1<<","<<hap2<<endl;
          ReadParser & parser = parsers[subject];
          cerr<<"Extracting region for left marker "<<g_left_marker<<" of markers length "<<g_markers<<"\n";
          parser.extract_region(config->chromosome.data(),g_left_marker,g_markers);
          //cerr<<"Printing data\n";
          //parser.print_data();
          // GPU friendly implementation
          // get total depth
          int total_depth = 0;
          int total_width = parser.MATRIX_WIDTH;
          for(uint matrow=0;matrow<parser.matrix_rows;++matrow){
            total_depth += parser.depth[matrow];
          }
          float prob_matrix1[total_depth*total_width];
          float prob_matrix2[total_depth*total_width];
          for(int i=0;i<total_depth*total_width;++i){
            prob_matrix1[i] = prob_matrix2[i] = 1;
          }
          int expanded_row = 0;
          
        
          // original implementation 
          float loglike = 0;
  
          for(uint matrow=0;matrow<parser.matrix_rows;++matrow){
            int offset = parser.offset[matrow];
            cerr<<"Working on main read "<<matrow<<" with offset "<<offset<<endl;
            // first get the depth of main read
            int depth = parser.depth[matrow];
            for(uint readindex=0;readindex<depth;++readindex){
              //float hap1_prob[parser.read_len[matrow]];
              //float hap2_prob[parser.read_len[matrow]];
              //float hap1_like=1,hap2_like=1;
              for(uint baseindex=0;baseindex<parser.read_len[matrow];
              ++baseindex){
                // APPROACH 1: GPU friendly
                prob_matrix1[expanded_row*total_width+offset+baseindex] =g_haplotype[hap1*g_max_window+baseindex]==
                parser.matrix_alleles[matrow*parser.MATRIX_WIDTH+baseindex]? 
                (parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex])):
                (1.-parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex]));
                prob_matrix2[expanded_row*total_width+offset+baseindex] =g_haplotype[hap2*g_max_window+baseindex]==
                parser.matrix_alleles[matrow*parser.MATRIX_WIDTH+baseindex]? 
                (parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex])):
                (1.-parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex]));

                // APPROACH 2: will not underflow
                //hap1_prob[baseindex] =g_haplotype[hap1*g_max_window+baseindex]==
                //parser.matrix_alleles[matrow*parser.MATRIX_WIDTH+baseindex]? 
                //(parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex])):
                //(1.-parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex]));
                //hap2_prob[baseindex] =g_haplotype[hap2*g_max_window+baseindex]==
                //parser.matrix_alleles[matrow*parser.MATRIX_WIDTH+baseindex]? 
                //(parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex])):
                //(1.-parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex]));

                // APPROACH 3: this approach below may underflow with longer haplotypes 
                //hap1_like*=g_haplotype[hap1*g_max_window+baseindex]==
                //parser.matrix_alleles[matrow*parser.MATRIX_WIDTH+baseindex]?
                //(parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex])):
                //(1.-parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex]));
                //hap2_like*=g_haplotype[hap2*g_max_window+baseindex]==
                //parser.matrix_alleles[matrow*parser.MATRIX_WIDTH+baseindex]?
                //(parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex])):
                //(1.-parser.phred2prob(parser.base_quality[matrow][readindex*parser.read_len[matrow]+offset+baseindex]));
              }
              //float like1 = log(.5*(hap1_like+hap2_like)); 
              //float like2 = get_bam_loglikelihood(parser.read_len[matrow],hap1_prob,hap2_prob);
              //cerr<<"For sub read "<<readindex<<" hap_probs are "<<hap1_like<<" and "<<hap2_like<<" and logprobs are "<<like1<<" and "<<like2<<endl;
              //loglike+=like2;
             
              
              //cout<<matrix_alleles[i*MATRIX_WIDTH+j];
              ++expanded_row;
            }   
          } 
          float loglike2 = 0;
          for(expanded_row=0;expanded_row<total_depth;++expanded_row){
            loglike2+=get_bam_loglikelihood(total_width,
            prob_matrix1+expanded_row*total_width,
            prob_matrix2+expanded_row*total_width);
          }
          cerr<<"Log likelihood of all reads for haps "<<hap1<<","<<hap2<<" is "<<loglike2<<endl;
          if (loglike2>max_log_pen) max_log_pen = loglike2;
          logpenetrance_cache[subject*penetrance_matrix_size+hap1*g_max_haplotypes+hap2] = loglike2;
        } // active hap condition
      } // second hap loop
    }// first hap loop
    //cerr<<"Max log penetrance for subject "<<subject<<" is "<<max_log_pen<<endl;
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
          g_max_haplotypes+hap2] = val>=logpenetrance_threshold?exp(val):0;
          penetrance_cache[subject*penetrance_matrix_size+hap2*
          g_max_haplotypes+hap1] = penetrance_cache[subject*
          penetrance_matrix_size+hap1*g_max_haplotypes+hap2];
          //cerr<<"Penetrance for subject "<<subject<<" haps "<<hap1<<","<<hap2<<": "<<penetrance_cache[subject*penetrance_matrix_size+hap1*g_max_haplotypes+hap2]<<endl;
        } // active hap condition
      } // second hap loop
    }// first hap loop
  }
  return;

  // if we are in sequencing mode, assume non polymorphic sites are ref allele 
  if (g_likelihood_mode == Config::LIKELIHOOD_MODE_READS && !polymorphic_window){
    cerr<<"Skipping non polymorphic site\n";
    return;
  }
  cerr<<"Exiting compute penetrance\n";
}


