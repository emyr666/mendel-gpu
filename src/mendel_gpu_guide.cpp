using namespace std;
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<tr1/unordered_set>
#include<tr1/unordered_map>
#include<set>
#include<map>
#include<list>
#include<math.h>
#include"cl_constants.h"
#include<cstdlib>
#include<string.h>
#ifdef USE_GPU
#include<CL/cl.hpp>
#include"clsafe.h"
#endif
#include"mendel_gpu.hpp"


void MendelGPU::parse_ref_haplotype(){
  string line;
  ifstream ifs(infile_refhap.data());
  if(!ifs.is_open()){
    cerr <<"Cannot open ref haplotypes file "<<infile_refhap<<endl;
    exit(1);
  }
  getline(ifs,line);
  ref_haplotypes = line.length();
  ifs.close();
  cerr<<"Parsing "<<ref_haplotypes<<" ref haplotypes across "<<g_snps<<" SNPs.\n";
  ref_haplotype = new char[ref_haplotypes * g_snps];
  ifs.open(infile_refhap.data());
  for(int i=0;i<g_snps;++i){
    getline(ifs,line);
    for(int j=0;j<ref_haplotypes;++j){
      ref_haplotype[j*g_snps+i] = line[j];
    }
  }
  ifs.close();
}

void MendelGPU::copy_ref_haplotypes_(){
  bool debug_haplotype = false;
  int left_marker = g_left_marker[0]-1;
  int marker_len = g_markers[0];
  // save the original # of markers from MENDEL
  extended_markers = marker_len;
  cerr<<"Including guide haplotypes of length: "<<marker_len<<" at left marker: "<<left_marker<<endl;
  // generate a hash table of unique template haplotypes
  std::tr1::unordered_map<string,set<string> > temp_hap_window;
  //map<string,set<string> > temp_hap_window;
  for(int i=0;i<g_max_haplotypes;++i) g_active_haplotype[i] = 0;
  g_informative_markers = 0;
  for(int j=0;j<marker_len;++j){
    if (informative_snp[left_marker+j]){
      extended_snp_mapping[g_informative_markers] = j;
      ++g_informative_markers;
    }
  }
  // notify MENDEL that the # of markers is temporarily reduced
  g_markers[0] = g_informative_markers;
  // tally up the counts of each reference haplotype
  ref_hap_counts.clear();
  for(int i=0;i<ref_haplotypes;++i){
    string curhap_long(ref_haplotype+i*g_snps+left_marker,marker_len);
    if (ref_hap_counts.find(curhap_long)==ref_hap_counts.end()){
      ref_hap_counts[curhap_long] = 0;
    }
    ++ref_hap_counts[curhap_long];
  }
  // rank the template haplotypes by their frequencies
  multiset<hapobj_t,byHapCountDesc> sorted_guides_extended;
  for(std::tr1::unordered_map<string,int >::iterator it =
  ref_hap_counts.begin();it!=ref_hap_counts.end();it++){
    hapobj_t hapobj;
    hapobj.hapstr = it->first; 
    hapobj.hapcount = it->second;
    sorted_guides_extended.insert(hapobj);
  }
  int extended_count = 0;
  for(multiset<hapobj_t,byHapCountDesc>::iterator it = 
  sorted_guides_extended.begin();it!=sorted_guides_extended.end();it++){
    if (extended_count<g_max_extended_haplotypes){
      hapobj_t hapobj = *it;
      string curhap_long = it->hapstr;
      char shorthap[g_informative_markers];
      for(int j=0;j<g_informative_markers;++j){
        shorthap[j] = curhap_long[extended_snp_mapping[j]];
      }
      string curhap_short(shorthap,g_informative_markers);
      if (temp_hap_window.find(curhap_short)==temp_hap_window.end()){
        set<string> temp;
        temp_hap_window[curhap_short] = temp;
      }
      temp_hap_window[curhap_short].insert(curhap_long);
      ++extended_count;
    }
  }
  // rank the template haplotypes by their frequencies
  multiset<hapobj_t,byHapCountDesc> sorted_guides;
  for(std::tr1::unordered_map<string,set<string> >::iterator it = temp_hap_window.begin();it!=temp_hap_window.end();it++){
    hapobj_t hapobj;
    hapobj.hapstr = it->first; 
    hapobj.hapcount = it->second.size();
    hapobj.extended_set = it->second;
    sorted_guides.insert(hapobj);
  }
  // ADD TO HAP WINDOW THE DENOVO HAPLOTYPES
  int totalcounts = 0;
  string hapstrarr[g_max_haplotypes];
  int countarr[g_max_haplotypes];
  full_hap_window.clear();
  // we need to know how many full haplotypes there are
  int compact_haplotypes = 0;
  extended_haplotypes = 0;
  int max_extended_haplotypes = g_max_haplotypes;
  for(multiset<hapobj_t,byHapCountDesc>::iterator it = sorted_guides.begin();it!=sorted_guides.end();it++){
    if (compact_haplotypes<g_max_haplotypes){
      hapobj_t hapobj = *it;
      for(int j=0;j<g_informative_markers;++j){
        g_haplotype[compact_haplotypes *g_max_window+j] = hapchar2int(hapobj.hapstr[j]);
      }
      g_active_haplotype[compact_haplotypes] = 1;
      hapstrarr[compact_haplotypes] = hapobj.hapstr;
      countarr[compact_haplotypes] = hapobj.hapcount;
      totalcounts+= hapobj.hapcount;
      if (debug_haplotype) cerr<<"Considering "<<hapobj.hapstr<<" count: "<<hapobj.hapcount<<endl;
      for(set<string>::iterator it=hapobj.extended_set.begin();
      it!=hapobj.extended_set.end();it++){
        if(debug_haplotype) cerr<<" extended: "<<*it<<endl;
        ++extended_haplotypes;
      }
      full_hap_window[hapobj.hapstr] = hapobj.extended_set;
      ++compact_haplotypes;
    }
  }
  // UPDATE TOTAL HAPLOTYPE SIZE
  g_haplotypes[0] = compact_haplotypes;
  cerr<<"Working haplotypes after adding templates: "<<compact_haplotypes<<endl;
  cerr<<"Extended haplotypes: "<<extended_haplotypes<<endl;
  if (debug_haplotype) cerr<<"STORING "<<g_haplotypes[0]<<
  " TRUE HAPLOTYPES AND FREQS:\n";
  for(int i=0;i<g_max_haplotypes;++i){
    if(g_active_haplotype[i]){
       g_frequency[i] = 1.*countarr[i]/totalcounts;
       if (debug_haplotype){
         cerr<<i<<":"<<hapstrarr[i]<<","<<g_frequency[i]<<endl;
       }
    }
  }
  if (debug_haplotype) exit(0);
}

void MendelGPU::precompute_penetrance_(){
  int left_marker = g_left_marker[0]-1;
  bool debug_penmat = left_marker<-1;
  cerr<<"Entering precompute_penetrance at left marker "<<left_marker<<"\n";
  int debug_subject = 200;
  for(int hapindex=0;hapindex<g_haplotypes[0];++hapindex){
    compress(hapindex,g_markers[0],g_haplotype);
  }
  if (run_gpu){
#ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotypes, CL_TRUE, 0, sizeof(int),g_haplotypes,NULL,NULL);
    clSafe(err, "write total haplotypes");
    err = commandQueue->enqueueWriteBuffer(*buffer_packedhap, CL_TRUE, 0,sizeof(packedhap_t)*g_max_haplotypes*packedhap_len, packedhap, NULL, NULL );
    clSafe(err, "write packedhaplotypes");
    err = commandQueue->enqueueWriteBuffer(*buffer_extended_snp_mapping, CL_TRUE, 0, sizeof(int)*g_max_window,extended_snp_mapping,NULL,NULL);
    clSafe(err, "write extended_snp_mapping");
    if (geno_dim==4){
      err = commandQueue->enqueueNDRangeKernel(*kernel_precompute_penetrance_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch precompute penetrance haploid");
    }else{
      err = commandQueue->enqueueNDRangeKernel(*kernel_precompute_penetrance,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch precompute penetrance");
    }
    if(debug_penmat){
      if (geno_dim==4){
        float * debug_cache = new float[g_people*g_max_haplotypes*2];
        //err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*2*g_max_haplotypes,debug_cache);
        //clSafe(err, "read logpenetrance cache");
        err = commandQueue->enqueueReadBuffer(*buffer_penetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*2*g_max_haplotypes,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"GPU:\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<debug_cache[debug_subject*2*g_max_haplotypes+2*j+k];
            }
            cout<<endl;
          }
        }
      }else{
        float * debug_cache = new float[g_people*penetrance_matrix_size];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*penetrance_matrix_size,debug_cache);
        clSafe(err, "read logpenetrance cache");
        //err = commandQueue->enqueueReadBuffer(*buffer_penetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*penetrance_matrix_size,debug_cache);
        //clSafe(err, "read penetrance cache");
        cout<<"GPU:\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<debug_cache[debug_subject*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
        delete[] debug_cache;
      }
    }
#endif
  }
  if(run_cpu){
    double start = clock();
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
      float maxlog = -1000;
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j]){
          //decompress(j,g_markers[0],g_haplotype);
          if (geno_dim==4){
            float logpenetrance[2] = {0,0};
            for(int l=0;l<g_markers[0];++l){
              //int dose = g_haplotype[j*g_max_window+l];
              bool dose = (((int)packedhap[j*packedhap_len+
              (l/32)].octet[(l%32)/8]) >> 
              ((l%32)%8) & 1) ;
              //cerr<<"dose1,2: "<<dose<<","<<dose2<<endl;

              for(int parent=0;parent<2;++parent){
                logpenetrance[parent]+=g_region_snp_penetrance[i*
                geno_dim*g_max_region_size+geno_dim*(left_marker+
                extended_snp_mapping[l]-g_region_snp_offset)+2*parent+dose]; 
                //cerr<<"Hap "<<j<<" Person "<<i<<" parent "<<parent<<" SNP "<<left_marker+extended_snp_mapping[l]<<" dose "<<dose<<" value "<<
                //g_region_snp_penetrance[i*
                //geno_dim*g_max_region_size+geno_dim*(left_marker+
                //extended_snp_mapping[l]-g_region_snp_offset)+2*parent+dose]<<
                //endl; 
                
              }
            }
            for(int parent=0;parent<2;++parent){
              //cerr<<"Hap "<<j<<" person "<<i<<" log penetrance for parent "<<parent<<": "<<logpenetrance[parent]<<endl;
            }
            //exit(0);
            for(int parent=0;parent<2;++parent){
              logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              logpenetrance[parent];
              if (logpenetrance[parent]>maxlog) maxlog = logpenetrance[parent];
            }
          }else{
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k] && (!haploid || j==k)){
                decompress(k,g_markers[0],g_haplotype);
                float logpenetrance = 0;
                for(int l=0;l<g_markers[0];++l){
                  int n=g_haplotype[j*g_max_window+l] + g_haplotype[k*g_max_window+l];
                  logpenetrance+=g_region_snp_penetrance[i*(geno_dim*g_max_region_size)+geno_dim*(left_marker+extended_snp_mapping[l]-g_region_snp_offset) +  n]; 
           
                }
                logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = logpenetrance;
                if (logpenetrance>maxlog) maxlog = logpenetrance;
              }
            }
          } 
        }
      }
//maxlog = 0;
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j]){
          if (geno_dim==4){
            //cerr<<"maxlog: "<<maxlog<<endl;
            float val[2];
            for(int parent=0;parent<2;++parent){
              val[parent] = logpenetrance_cache[i*2*g_max_haplotypes+
              2*j+parent]-maxlog;
              logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              val[parent];
              penetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              val[parent]>=logpenetrance_threshold?exp(val[parent]):0;
              //cerr<<"hap:"<<j<<",person "<<i<<",parent:"<<parent<<",logval:"<<logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent]<<",val:"<<penetrance_cache[i*2*g_max_haplotypes+2*j+parent]<<endl;
            }
          }else{
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k] && (!haploid || j==k)){
                float val = logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]-maxlog;
                logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = val;
                penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = val>=logpenetrance_threshold?exp(val):0;
              }
            }
          }
        }
      }
    }
    if(debug_penmat){
      cout<<"CPU:\n";
      if (geno_dim==4){
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<penetrance_cache[debug_subject*2*g_max_haplotypes
              +2*j+k];
              //cout<<" "<<logpenetrance_cache[debug_subject*2*g_max_haplotypes
              //+2*j+k];
            }
            cout<<endl;
          }
        }
      }else{
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<logpenetrance_cache[debug_subject*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
      }
    }
    cerr<<"Standard compute penetrance time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
  }
  if (debug_penmat) exit(0);
  cerr<<"Leaving precompute_penetrance\n";
}

void MendelGPU::prep_impute_genotypes_guide_(){
  int left_marker = g_left_marker[0]-1;
  cerr<<"begin prep impute_genotypes at left SNP "<<left_marker<<endl;
  int max_geno=g_genotype_imputation?3:4;
  int curhapindex = 0;
  if (debug_mode) ofs_debug_haplotype_file<<"Extended haplotypes for left marker "<<left_marker<<endl;
  for(int i=0;i<g_haplotypes[0];++i){
    string hapstring = hapstr(g_haplotype+i*g_max_window,g_markers[0]);
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
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotypes, CL_TRUE, 0,sizeof(int), g_haplotypes, NULL, NULL );
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

void MendelGPU::impute_diploid_genotypes_guide_(int * center_snp_start, 
int * center_snp_end, int * center_snp_offset){
  int max_geno=g_genotype_imputation?3:4;
  int c_snp_start = *center_snp_start-1;
  int c_snp_end = *center_snp_end;
  int c_snp_offset = *center_snp_offset-1;
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
        if(subject_haplotype_weight[i*g_max_haplotypes+
          extended_root_mapping[j]]>0){ 
          //decompress(j,extended_markers,extended_haplotype);
          for(int k=j;k<extended_haplotypes;++k){
            if((!haploid ||  extended_root_mapping[j] == 
            extended_root_mapping[k]) && subject_haplotype_weight[i*
            g_max_haplotypes+extended_root_mapping[k]]>0){ 
              //decompress(k,extended_markers,extended_haplotype);
              float penetrance =
              penetrance_cache[i*penetrance_matrix_size+extended_root_mapping[j]*
              g_max_haplotypes+extended_root_mapping[k]];
              if (penetrance>0){
                for(int c_snp = c_snp_start;c_snp<c_snp_end; ++c_snp){
                  int packedsite = c_snp - c_snp_start;
                  int current_snp = c_snp_offset + c_snp;
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
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_genotype_file<<current_snp<<"\t";
        ofs_dosage_file<<current_snp<<"\t";
        ofs_posterior_file<<current_snp;
      }
      for(int i=0;i<g_people;++i){
        // compute normalizing constant
        float denom = 0;
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
        if (debug_geno) cout<<"CPUGENO:\t"<<i<<"\t"<<current_snp<<"\t"<<genotype<<endl;
        if (outfile_format.compare(FORMAT_MEC)==0){
          ofs_genotype_file<<genotype;
        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
          ofs_genotype_file<<"CPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<genotype<<endl;
        }
        float dose = 0;
        if (outfile_format.compare(FORMAT_MEC)==0){
          for(int j=0;j<max_geno;++j){
            float p = subject_posterior_prob[i*g_flanking_snps*4+
            (c_snp-c_snp_start)*4+j]/denom;
            ofs_posterior_file<<"\t"<<p;
            dose+=j*p;
          }
          char chardose =  (int)(dose*10)+'A';
          ofs_dosage_file<<chardose;
        }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
          ofs_posterior_file<<"CPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
          for(int j=0;j<max_geno;++j){
            float p = subject_posterior_prob[i*g_flanking_snps*4+
            (c_snp-c_snp_start)*4+j]/denom;
            if(debug_mode)ofs_posterior_file<<"\t"<< p;
            dose+=j*p;
          }
          ofs_dosage_file<<"CPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<dose<<endl;
          ofs_posterior_file<<endl;
        }
        subject_dosages[i*g_flanking_snps+(c_snp-c_snp_start)] = dose;
      }
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_genotype_file<<endl;
        ofs_dosage_file<<endl;
      }
      float rsq = compute_rsq(subject_dosages,
      g_flanking_snps,(c_snp-c_snp_start));
      ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
    }
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_posterior|| debug_dosage||debug_geno) exit(1);
  if (debug_pen) {
    if (debug_mode) ofs_debug_haplotype_file.close();
    exit(1);
  }
}

void MendelGPU::impute_haploid_genotypes_guide_(int * center_snp_start, 
int * center_snp_end, int * center_snp_offset){
  int max_geno=g_genotype_imputation?3:4;
  int c_snp_start = *center_snp_start-1;
  int c_snp_end = *center_snp_end;
  int c_snp_offset = *center_snp_offset-1;
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


