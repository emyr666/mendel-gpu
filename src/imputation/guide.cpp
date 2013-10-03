using namespace std;
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<queue>
#include<tr1/unordered_set>
#include<tr1/unordered_map>
#include<set>
#include<map>
#include<list>
#include<math.h>
#include"../cl_constants.h"
#include<cstdlib>
#include<string.h>
#ifdef USE_GPU
#include<CL/cl.hpp>
#include"../clsafe.h"
#endif
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"


void GuidedMendelGPU::impute_genotypes(){
  prep_impute_genotypes_guide();
  impute_diploid_genotypes_guide(1,1,1);
}

void GuidedMendelGPU::prep_impute_genotypes_guide(){
  int left_marker = g_left_marker;
  cerr<<"begin prep impute_genotypes at left SNP "<<left_marker<<endl;
  int max_geno=g_genotype_imputation?3:4;
  int curhapindex = 0;
  if (debug_mode) ofs_debug_haplotype_file<<"Extended haplotypes for left marker "<<left_marker<<endl;
  for(int i=0;i<g_haplotypes;++i){
    string hapstring = hapstr(g_haplotype+i*g_max_window,g_markers);
    if (full_hap_window.find(hapstring)==full_hap_window.end()){
      cerr<<"Orphaned hapstr "<<hapstring<<endl;
      exit(1);
    }
    set<string> extendedhaps = full_hap_window[hapstring];
    float freq = 1.*g_frequency[i]/(extendedhaps.size());
    if (debug_mode) ofs_debug_haplotype_file<<"Compact haplotype: "<<hapstring<<endl;
    float totalcounts = 0;
    int j = 0;
    float extendedfreq[extendedhaps.size()];
    for(set<string>::iterator it  = extendedhaps.begin();
    it!=extendedhaps.end();it++){
       string str = *it;
       if (ref_hap_counts.find(str)==ref_hap_counts.end()){
         cerr<<"Exception: couldn't find refhap counts for "<<str<<endl;
       }
       extendedfreq[j] = ref_hap_counts[str];
       totalcounts+=1.*extendedfreq[j];
       ++j;
    }
    for(int j=0;j<extendedhaps.size();++j) extendedfreq[j] = extendedfreq[j]/totalcounts*g_frequency[i]*1.;
    j = 0;
    for(set<string>::iterator it  = extendedhaps.begin();
    it!=extendedhaps.end();it++){
       string str = *it;
       for(int k=0;k<extended_markers;++k){
         //cerr<<"Making extended hap for "<< curhapindex*g_max_window+k<<"\n";
         extended_haplotype[curhapindex*g_max_window+k] = hapchar2int(str[k]);
       }
       extended_frequency[curhapindex] = extendedfreq[j];
       if (debug_mode) ofs_debug_haplotype_file<<"Extended hap: "<<str<<" base freq "<<g_frequency[i]<<","<<ref_hap_counts[str]<<", extendedfreq: "<<extendedfreq[j]<<endl;
       extended_root_mapping[curhapindex] = i;
       ++j;
       ++curhapindex;
    }
  }
  if (curhapindex!=extended_haplotypes) {
    cerr<<"Assertion failed curhapindex ("<<curhapindex<<")!=extended_haplotypes("<<extended_haplotypes<<")\n ";
    exit(1);
  }
  if (run_gpu){
    #ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotypes, CL_TRUE, 0,sizeof(int), &g_haplotypes, NULL, NULL );
    clSafe(err, "write extended_haplotypes");
    err = commandQueue->enqueueWriteBuffer(*buffer_extended_haplotypes, CL_TRUE, 0,sizeof(int), &extended_haplotypes, NULL, NULL );
    clSafe(err, "write extended_haplotypes");
    err = commandQueue->enqueueWriteBuffer(*buffer_extended_root_mapping, CL_TRUE, 0,sizeof(int)*ref_haplotypes, extended_root_mapping, NULL, NULL );
    clSafe(err, "write extended root mapping");
    err = commandQueue->enqueueWriteBuffer(*buffer_extended_frequency, CL_TRUE, 0,sizeof(float)*ref_haplotypes, extended_frequency, NULL, NULL );
    clSafe(err, "write extended_freq");
    #endif
  }
  cerr<<"done prep impute_genotypes at left SNP "<<left_marker<<endl;
}

void GuidedMendelGPU::impute_diploid_genotypes_guide(int  center_snp_start, 
int  center_snp_end, int  center_snp_offset){
  int max_geno=g_genotype_imputation?3:4;
  int c_snp_offset = g_center_snp_start;
  int c_snp_start = g_center_snp_start - g_left_marker;
  int c_snp_end = c_snp_start+g_flanking_snps;
  bool debug_geno = false;
  bool debug_dosage = c_snp_start==-50;
  bool debug_posterior = c_snp_start==-50;
  bool debug_pen = false;
  cerr<<"Imputing from "<<c_snp_start<<" to "<<(c_snp_end-1)<<" with offset "<<c_snp_offset<<endl;

  for(int hapindex=0;hapindex<extended_haplotypes;++hapindex){
    compress_extendedhap(hapindex,c_snp_start,c_snp_end,extended_haplotype);
    //compress(hapindex,extended_markers,extended_haplotype);
  }
  cerr<<"Compressed haplotypes\n";

  //int current_snp = *d_snp-1;
  //int c_snp = *center_snp-1;
  //cerr<<"begin impute_genotypes_guide at SNP "<<current_snp<<", center: "<<c_snp<<endl;
  //return ;
  for(int i=0;i<extended_haplotypes;++i){
    //extended_center_dosage[i] = extended_haplotype[i*g_max_window+c_snp];
  }
  if (run_gpu){
#ifdef USE_GPU
    double start = clock();
    int last_site = c_snp_end-c_snp_start; 
    err = commandQueue->enqueueWriteBuffer(*buffer_center_snp_end, CL_TRUE, 0,sizeof(int), &last_site, NULL, NULL );
    clSafe(err, "write centersnp end");
    err = commandQueue->enqueueWriteBuffer(*buffer_packedextendedhap, CL_TRUE, 0,sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len, packedextendedhap, NULL, NULL );
    clSafe(err, "write packedhap");
    cerr<<"Launching impute geno with max site "<<last_site<<endl;
    err = commandQueue->enqueueNDRangeKernel(*kernel_impute_genotype_guide,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH_IMPUTE_GUIDE,last_site),cl::NDRange(BLOCK_WIDTH_IMPUTE_GUIDE,1),NULL,NULL);
    clSafe(err,"launch impute_genotype guide");
    // READ GENOTYPE POSTERIORS
    float subject_posterior_prob_block[g_people * g_flanking_snps * 4];
    err = commandQueue->enqueueReadBuffer(*buffer_subject_posterior_prob_block, CL_TRUE, 0, sizeof(float)*g_people*g_flanking_snps*4,subject_posterior_prob_block);
    clSafe(err, "read subject_posterior");
    int len = c_snp_end-c_snp_start;
    for(int j=0;j<len;++j){ // loop over SNPs j
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_posterior_file<<j;
        for(int i=0;i<g_people;++i){
          for(int k=0;k<max_geno;++k){
            ofs_posterior_file<<"\t"<<
            subject_posterior_prob_block[i*g_flanking_snps*4+j*4+k];
          }
        }
        ofs_posterior_file<<endl;
      }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
        for(int i=0;i<g_people;++i){
          ofs_posterior_file<<"GPU_POSTERIOR\t"<<i<<"\t"<<j;
          for(int k=0;k<max_geno;++k){
            ofs_posterior_file<<"\t"<<subject_posterior_prob_block[i*g_flanking_snps*4+j*4+k];
          }
          ofs_posterior_file<<endl;
        }
      }
    }
    int subject_geno_block[g_people*g_flanking_snps];
    err = commandQueue->enqueueReadBuffer(*buffer_subject_genotype_block, CL_TRUE, 0, sizeof(int)*g_people*g_flanking_snps,subject_geno_block);
    clSafe(err, "read subject_genotype");
    cerr<<"Personloop: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
      len = c_snp_end-c_snp_start;
      for(int j=0;j<len;++j){
        int current_snp = c_snp_offset+c_snp_start+j;
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_genotype_file<<current_snp<<"\t";
        }
        for(int i=0;i<g_people;++i){
          if(debug_geno) cout<<"GPUGENO: "<<(c_snp_offset+c_snp_start+j)<<" "<<i<<" "<<subject_geno_block[i*g_flanking_snps+j]<<endl;
          if (outfile_format.compare(FORMAT_MEC)==0){
            ofs_genotype_file<<subject_geno_block[i*g_flanking_snps+j];
          }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
            ofs_genotype_file<<"GPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_geno_block[i*g_flanking_snps+j]<<endl;
          }else{
       //     ofs_genotype_file<<" "<<outfile_format;
          }
        }
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_genotype_file<<endl;
        }
      }
//    for(int i=0;i<g_people;++i){
//      ofs_genotype_file<<"GPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_genotypes[i]<<endl;
//    }
    if (g_genotype_imputation){
      float subject_dosage_block[g_people*g_flanking_snps];
      err = commandQueue->enqueueReadBuffer(*buffer_subject_dosage_block, CL_TRUE, 0, sizeof(float)*g_people*g_flanking_snps,subject_dosage_block);
      clSafe(err, "read subject_dosage");
      int len = c_snp_end-c_snp_start;
      for(int j=0;j<len;++j){
        int current_snp = c_snp_offset+c_snp_start+j;
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_dosage_file<<current_snp<<"\t";
        }
        for(int i=0;i<g_people;++i){
          if (debug_dosage) cout<<"GPUDOSE: "<<current_snp<<" "<<i<<" "<<subject_dosage_block[i*g_flanking_snps+j]<<endl;
          if (outfile_format.compare(FORMAT_MEC)==0){
            char chardose = 
            (int)(subject_dosage_block[i*g_flanking_snps+j]*10)+'A';
            ofs_dosage_file<<chardose;
            //ofs_dosage_file<<subject_dosage_block[i*g_flanking_snps+j]<<endl;
          }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
            ofs_dosage_file<<"GPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosage_block[i*g_flanking_snps+j]<<endl;
          }
        }
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_dosage_file<<endl;
        }
        float rsq = compute_rsq(subject_dosage_block, g_flanking_snps,j);
        ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
      }
    }
//      for(int i=0;i<g_people;++i){
//        //float dose = 0;
//        //for(int j=0;j<max_geno;++j){
//        //  dose+=j*subject_posterior[i*4+j]/denom;
//        //}
//        ofs_dosage_file<<"GPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosages[i]<<endl;
//        if (debug_dosage) cout<<"GPU DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosages[i]<<endl;
//      }
//      float rsq = compute_rsq(subject_dosages);
//      ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
//    }
//    cerr<<"Elapsed time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
//    #endif
#endif
  }
  if(run_cpu){
    float subject_dosages[g_people * g_flanking_snps];
    float subject_posterior_prob[g_people * g_flanking_snps * 4];
    memset(subject_posterior_prob,0,sizeof(float)*g_people*g_flanking_snps * 4);
    double start = clock();
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
      for(int j=0;j<extended_haplotypes;++j){
        //cerr<<"subject hap weight for hap "<<extended_root_mapping[j]<<" is "<<subject_haplotype_weight[i*g_max_haplotypes+extended_root_mapping[j]] <<" for subject "<<i<<endl;
        if(subject_haplotype_weight[i*g_max_haplotypes+
          extended_root_mapping[j]]>0){ 
          //decompress(j,extended_markers,extended_haplotype);
          for(int k=j;k<extended_haplotypes;++k){
            if((!haploid ||  extended_root_mapping[j] == 
            extended_root_mapping[k]) && subject_haplotype_weight[i*
            g_max_haplotypes+extended_root_mapping[k]]>0){ 
              //cerr<<"subject hap weight for hap "<<k<<" is pos for subject "<<i<<endl;
              //decompress(k,extended_markers,extended_haplotype);
              float penetrance =
              penetrance_cache[i*penetrance_matrix_size+extended_root_mapping[j]*
              g_max_haplotypes+extended_root_mapping[k]];
              if (penetrance>0){
                for(int c_snp = c_snp_start;c_snp<c_snp_end; ++c_snp){
                  int packedsite = c_snp - c_snp_start;
                  int current_snp = c_snp_offset + c_snp;
//cerr<<"packed site and current snp: "<<packedsite<<","<<current_snp<<endl;
                  int m;
                  if(g_genotype_imputation){
                    m = (((int)packedextendedhap[j*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) + (((int)packedextendedhap[k*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) ;
                  }else{
                    m = 2*(((int)packedextendedhap[j*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) + (((int)packedextendedhap[k*packedextendedhap_len+(packedsite/32)].octet[(packedsite%32)/8]) >> ((packedsite%32)%8) & 1) ;
                  }
                  float freq = extended_frequency[j]*extended_frequency[k];
                  if (j!=k) freq*=2;
                  float p = freq*penetrance;
                  //p = 1;
                  subject_posterior_prob[i*g_flanking_snps*4+(packedsite)*4+m]+=p;
                  if (debug_pen) cerr<<"i,j,k,m,pen,freq:"<<i<<","<<j<<","<<k<<","<<m<<","<<penetrance<<","<<freq<<endl;
                  //if (debug_pen) cerr<<"packedsite: "<<packedsite<<endl;
                }
              }
            }
          }
        }
      }
      //cerr<<" "<<i;
    }
    cerr<<"Personloop: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    //cerr<<endl;
      //cout<<"CPU POSTERIOR "<<current_snp<<","<<i;
      //for(int j=0;j<max_geno;++j){
       // cout<<" "<<posterior_prob[j];
      //}
      //cout<<endl;
    if (debug_posterior){
      int len = c_snp_end-c_snp_start;
      for(int j=0;j<len;++j){
         for(int i=0;i<g_people;++i){
           cout<<"CPUPOSTERIOR: "<<j<<" "<<i;
           for(int k=0;k<max_geno;++k){
             cout<<" "<<subject_posterior_prob[i*g_flanking_snps*4+j*4+k];
           }
           cout<<endl;
         }
      }
    }
    for(int c_snp = c_snp_start;c_snp<c_snp_end; ++c_snp){
      int current_snp = c_snp_offset + c_snp;
//      if (outfile_format.compare(FORMAT_MEC)==0){
//        ofs_genotype_file<<current_snp<<"\t";
//        ofs_dosage_file<<current_snp<<"\t";
//        ofs_posterior_file<<current_snp;
//      }
      float dosages[g_people];
      float posteriors[g_people*max_geno];
      int genotypes[g_people];
      for(int i=0;i<g_people;++i){
        // compute normalizing constant
        float denom = epsilon;
        float maxval = 0;
        int genotype = 0;
        for(int m=0;m<max_geno;++m){
          float p = subject_posterior_prob[i*g_flanking_snps*4+
          (c_snp-c_snp_start)*4+m];
          denom+=p;
          if (p>maxval) {
            maxval = p;
            genotype = m;
          }
        }
        //cerr<<"Person "<<i<<" denom is "<<denom<<endl;
        if (debug_geno) cout<<"CPUGENO:\t"<<i<<"\t"<<current_snp<<"\t"<<genotype<<endl;
        genotypes[i] = genotype;
//        if (outfile_format.compare(FORMAT_MEC)==0){
//          ofs_genotype_file<<genotype;
//        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
//          ofs_genotype_file<<"CPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<genotype<<endl;
//        }
        dosages[i] = 0;
        for(int j=0;j<max_geno;++j){
          int snpindex = c_snp-c_snp_start;
          float p = subject_posterior_prob[i*g_flanking_snps*4+
          (c_snp-c_snp_start)*4+j]/denom;
          posteriors[i*max_geno+j] = p;
          dosages[i]+=j*p;
        }
        subject_dosages[i*g_flanking_snps+(c_snp-c_snp_start)] = dosages[i];
 

//        float dose=0;
//        if (outfile_format.compare(FORMAT_MEC)==0){
//          for(int j=0;j<max_geno;++j){
//            float p = subject_posterior_prob[i*g_flanking_snps*4+
//            (c_snp-c_snp_start)*4+j]/denom;
//            //ofs_posterior_file<<"\t"<<p;
//            dose+=j*p;
//          }
//          char chardose =  (int)(dose*10)+'A';
//          //ofs_dosage_file<<chardose;
//        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
//          //ofs_posterior_file<<"CPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
//          for(int j=0;j<max_geno;++j){
//            float p = subject_posterior_prob[i*g_flanking_snps*4+
//            (c_snp-c_snp_start)*4+j]/denom;
//            if(debug_mode)ofs_posterior_file<<"\t"<< p;
//            dose+=j*p;
//          }
//          //ofs_dosage_file<<"CPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<dose<<endl;
//          //ofs_posterior_file<<endl;
//        }
      }
      //if (outfile_format.compare(FORMAT_MEC)==0){
        //ofs_genotype_file<<endl;
        //ofs_dosage_file<<endl;
      //}
      float rsq = compute_rsq(subject_dosages,
      g_flanking_snps,(c_snp-c_snp_start));
      if (max_geno==3){
        io_manager->writeDosage(current_snp,dosages,g_people);
        io_manager->writeGenotype(current_snp,genotypes,g_people);
      }
      io_manager->writePosterior(max_geno,current_snp,posteriors,g_people);

      //ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
    }
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_posterior|| debug_dosage||debug_geno) exit(1);
  if (debug_pen) {
    if (debug_mode) ofs_debug_haplotype_file.close();
    exit(1);
  }
}

void GuidedMendelGPU::impute_haploid_genotypes_guide(int  center_snp_start, 
int  center_snp_end, int  center_snp_offset){
  int max_geno=g_genotype_imputation?3:4;
  int c_snp_start = center_snp_start;
  int c_snp_end = center_snp_end-1;
  int c_snp_offset = center_snp_offset;
  bool debug_geno = false;
  bool debug_dosage = c_snp_start==-50;
  bool debug_posterior = false;
  bool debug_pen = false;
  cerr<<"Imputing from "<<c_snp_start<<" to "<<(c_snp_end-1)<<" with offset "<<c_snp_offset<<endl;
  for(int hapindex=0;hapindex<extended_haplotypes;++hapindex){
    compress_extendedhap(hapindex,c_snp_start,c_snp_end,extended_haplotype);
  }
  cerr<<"Compressed haplotypes\n";
  if (run_gpu){
#ifdef USE_GPU
    double start = clock();
    int last_site = c_snp_end-c_snp_start; 
    err = commandQueue->enqueueWriteBuffer(*buffer_center_snp_end, CL_TRUE, 0,sizeof(int), &last_site, NULL, NULL );
    clSafe(err, "write centersnp end");
    err = commandQueue->enqueueWriteBuffer(*buffer_packedextendedhap, CL_TRUE, 0,sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len, packedextendedhap, NULL, NULL );
    clSafe(err, "write packedhap");
    cerr<<"Launching impute geno with max site "<<last_site<<endl;
    err = commandQueue->enqueueNDRangeKernel(*kernel_impute_genotype_guide_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH_IMPUTE_GUIDE,last_site),cl::NDRange(BLOCK_WIDTH_IMPUTE_GUIDE,1),NULL,NULL);
    clSafe(err,"launch impute_genotype guide haploid");
    // READ GENOTYPE POSTERIORS
    float subject_posterior_prob_block[g_people * g_flanking_snps * 4];
    err = commandQueue->enqueueReadBuffer(*buffer_subject_posterior_prob_block, CL_TRUE, 0, sizeof(float)*g_people*g_flanking_snps*4,subject_posterior_prob_block);
    clSafe(err, "read unnormalized subject_posterior");
    float subject_dosages[g_people * g_flanking_snps];
    for(int c_snp = c_snp_start;c_snp<c_snp_end;++c_snp){
      int current_snp = c_snp_offset + c_snp;
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_genotype_file<<current_snp<<"\t";
        ofs_dosage_file<<current_snp<<"\t";
        ofs_posterior_file<<current_snp;
      }
      for(int i=0;i<g_people;++i){
        float * posterior_ptr = subject_posterior_prob_block + 
        (i*g_flanking_snps*4+(c_snp-c_snp_start)*4);
        for(int j=0;j<4;++j) posterior_ptr[j] = 
        posterior_ptr[j]<epsilon?epsilon:posterior_ptr[j];
        // compute normalizing constant
        float denom[2] = {0,0};
        float dose = 0;
        for(int parent=0;parent<2;++parent){
          for(int allele=0;allele<2;++allele){
            denom[parent] += posterior_ptr[2*parent+allele];
          }
          for(int allele=0;allele<2;++allele){
            posterior_ptr[parent*2+allele]/=denom[parent];
          }
        }
        if (debug_posterior){
          cout<<"GPU_POSTERIOR:\t"<<c_snp;
          for(int j=0;j<4;++j) cout<<"\t"<<posterior_ptr[j];
          cout<<endl;
        }
        if (outfile_format.compare(FORMAT_MEC)==0){
          for(int parent=0;parent<2;++parent){
            for(int allele=0;allele<2;++allele){
              ofs_posterior_file<<"\t"<<posterior_ptr[2*parent+allele];
              dose+=allele*posterior_ptr[2*parent+allele];
            }
          }
          char chardose =  (int)(dose*10)+'A';
          ofs_dosage_file<<chardose;
        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
          ofs_posterior_file<<"GPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
          for(int parent=0;parent<2;++parent){
            for(int allele=0;allele<2;++allele){
              ofs_posterior_file<<"\t"<<posterior_ptr[2*parent+allele];
              dose+=allele*posterior_ptr[2*parent+allele];
            }
          }
          ofs_posterior_file<<endl;
          ofs_dosage_file<<"GPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<dose<<endl;
        }
        subject_dosages[i*g_flanking_snps+(c_snp-c_snp_start)] = dose;
        float genoprobs[3];
        genoprobs[0] = posterior_ptr[0]*posterior_ptr[2];
        genoprobs[1] = posterior_ptr[1]*posterior_ptr[2] + 
                       posterior_ptr[0]*posterior_ptr[3];
        genoprobs[2] = posterior_ptr[1]*posterior_ptr[3];
        //cerr<<"genoprobs: "<<genoprobs[0]<<" "<<genoprobs[1]<<" "<<genoprobs[2]<<endl;
        float maxval = genoprobs[0];
        int geno2=0;
        for(int j=1;j<3;++j){
          if (genoprobs[j]>maxval){
            maxval = genoprobs[j];
            geno2 = j;
          }
        }
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_genotype_file<<geno2;
        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
          ofs_genotype_file<<"GPU_GENO:\t"<<i<<"\t"<<current_snp<<
          "\t"<<geno2<<endl;
        }
      }
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_genotype_file<<endl;
        ofs_dosage_file<<endl;
        ofs_posterior_file<<endl;
      }
      float rsq = compute_rsq(subject_dosages,
      g_flanking_snps,(c_snp-c_snp_start));
      ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
    }
#endif
  }
  if(run_cpu){
    float subject_dosages[g_people * g_flanking_snps];
    float subject_posterior_prob[g_people * g_flanking_snps * 4];
    memset(subject_posterior_prob,0,sizeof(float)*g_people*g_flanking_snps * 4);
    double start = clock();
//debug_haplotypes(cerr);
    for(int i=0;i<g_people;++i){
      for(int j=0;j<extended_haplotypes;++j){
        if(subject_haplotype_weight[i*g_max_haplotypes+
        extended_root_mapping[j]]>0){ 
          //cerr<<"extendedrootmapping for ext hap "<<j<<" is "<<extended_root_mapping[j]<<endl;
          for(int parent=0;parent<2;++parent){
            //float penetrance = 1;
            //penetrance_cache[i*2*g_max_haplotypes+2*extended_root_mapping[j] +
            //parent];
            float penetrance =
            penetrance_cache[i*2*g_max_haplotypes+2*extended_root_mapping[j] +
            parent];
            if (penetrance>0){
              float logpenetrance = logpenetrance_cache[i*2*g_max_haplotypes+2*extended_root_mapping[j] +
              parent];
              for(int c_snp = c_snp_start;c_snp<c_snp_end; ++c_snp){
                int packedsite = c_snp - c_snp_start;
                int current_snp = c_snp_offset + c_snp;
                int m;
                if(g_genotype_imputation){
                  m = (((int)packedextendedhap[j*packedextendedhap_len+
                  (packedsite/32)].octet[(packedsite%32)/8]) >> 
                  ((packedsite%32)%8) & 1) ;
                }
                float freq = extended_frequency[j];
                //float p = freq*exp(logpenetrance);
                float p = freq*penetrance;
                subject_posterior_prob[i*g_flanking_snps*4+(packedsite)*4+
                2*parent+m]+=p;

if (i==-10 ){
                cerr<<"person "<<i<<",hap "<<j<<",packedsite "<<packedsite<<",dose "<<m<<",parent "<<parent<<",p "<<p<<",post prob "<<
                subject_posterior_prob[i*g_flanking_snps*4+(packedsite)*4+
                2*parent+m]<<","<<
                subject_posterior_prob[i*g_flanking_snps*4+(packedsite)*4+
                2*parent+0]<<","<<
                subject_posterior_prob[i*g_flanking_snps*4+(packedsite)*4+
                2*parent+1]<<",logpen "<<logpenetrance
                <<endl;
}

              }
            }
          }
        }
      }
    }
    cerr<<"Personloop: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    for(int c_snp = c_snp_start;c_snp<c_snp_end; ++c_snp){
      int current_snp = c_snp_offset + c_snp;
  //cerr<<"currentsnp: "<<current_snp<<endl;
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_genotype_file<<current_snp<<"\t";
        ofs_dosage_file<<current_snp<<"\t";
        ofs_posterior_file<<current_snp;
      }
      for(int i=0;i<g_people;++i){
        float * posterior_ptr = subject_posterior_prob + (i*g_flanking_snps*4+
        (c_snp-c_snp_start)*4);
        for(int j=0;j<4;++j){
          if (i==-10) cerr<<"site: "<<(c_snp-c_snp_start)<<", final person "<<i<<",j "<<j<<",post "<<posterior_ptr[j]<<endl;
          posterior_ptr[j] = 
          posterior_ptr[j]<epsilon?epsilon:posterior_ptr[j];
        }
        // compute normalizing constant
        float denom[2] = {0,0};
        float dose = 0;
        for(int parent=0;parent<2;++parent){
          for(int allele=0;allele<2;++allele){
            denom[parent] += posterior_ptr[2*parent+allele];
          }
          for(int allele=0;allele<2;++allele){
            posterior_ptr[parent*2+allele]/=denom[parent];
          }
        }
        if (debug_posterior){
          cout<<"CPU_POSTERIOR:\t"<<c_snp;
          for(int j=0;j<4;++j) cout<<"\t"<<posterior_ptr[j];
          cout<<endl;
        }
        //cerr<<"Posteriors: ";
        //for(int j=0;j<4;++j) cerr<<" "<<posterior_ptr[j];
        //cerr<<endl;
        if (outfile_format.compare(FORMAT_MEC)==0){
          for(int parent=0;parent<2;++parent){
            for(int allele=0;allele<2;++allele){
              //float p = subject_posterior_prob[i*g_flanking_snps*4+
              //(c_snp-c_snp_start)*4+2*parent+allele]/denom[parent];
              ofs_posterior_file<<"\t"<<posterior_ptr[2*parent+allele];
              dose+=allele*posterior_ptr[2*parent+allele];
            }
          }
          char chardose =  (int)(dose*10)+'A';
          ofs_dosage_file<<chardose;
        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
          ofs_posterior_file<<"CPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
          for(int parent=0;parent<2;++parent){
            for(int allele=0;allele<2;++allele){
              ofs_posterior_file<<"\t"<<posterior_ptr[2*parent+allele];
              dose+=allele*posterior_ptr[2*parent+allele];
            }
          }
          ofs_posterior_file<<endl;
          ofs_dosage_file<<"CPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<dose<<endl;
        }
        subject_dosages[i*g_flanking_snps+(c_snp-c_snp_start)] = dose;
        float genoprobs[3];
        genoprobs[0] = posterior_ptr[0]*posterior_ptr[2];
        genoprobs[1] = posterior_ptr[1]*posterior_ptr[2] + 
                       posterior_ptr[0]*posterior_ptr[3];
        genoprobs[2] = posterior_ptr[1]*posterior_ptr[3];
        //cerr<<"genoprobs: "<<genoprobs[0]<<" "<<genoprobs[1]<<" "<<genoprobs[2]<<endl;
        float maxval = genoprobs[0];
        int geno2=0;
        for(int j=1;j<3;++j){
          if (genoprobs[j]>maxval){
            maxval = genoprobs[j];
            geno2 = j;
          }
        }
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_genotype_file<<geno2;
        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
          ofs_genotype_file<<"CPU_GENO:\t"<<i<<"\t"<<current_snp<<
          "\t"<<geno2<<endl;
        }
      }
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_genotype_file<<endl;
        ofs_dosage_file<<endl;
        ofs_posterior_file<<endl;
      }
      float rsq = compute_rsq(subject_dosages,
      g_flanking_snps,(c_snp-c_snp_start));
      ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
    }
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_posterior|| debug_dosage||debug_geno) exit(1);
}


