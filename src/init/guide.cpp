using namespace std;
#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"

GuidedMendelGPU::GuidedMendelGPU(IO_manager * io):MendelGPU(io){
}

int GuidedMendelGPU::get_max_window_size(){
  return 3 * g_flanking_snps;
}

void GuidedMendelGPU::allocate_memory(){
  MendelGPU::allocate_memory();
  cerr<<"Initializing ref haplotype specific variables\n";
  cerr<<"Maximum ref haps: "<<ref_haplotypes<<endl;
  g_informative_haplotype = new int[ref_haplotypes * g_max_window];
  g_max_extended_haplotypes = config->max_reference_haplotypes;
  packedhap_len = g_max_window/32 + (g_max_window%32!=0);
  packedhap = new packedhap_t[g_max_haplotypes * packedhap_len];
  packedextendedhap_len = g_flanking_snps/32 + (g_flanking_snps%32!=0);
  packedextendedhap = new packedhap_t[ref_haplotypes * packedextendedhap_len];
  extended_snp_mapping = new int[g_max_window];
  extended_haplotype = new int[ref_haplotypes * g_max_window]; 
  extended_frequency = new float[ref_haplotypes];
  extended_root_mapping = new int[ref_haplotypes];

  g_left_marker = 0;
  g_center_snp_start = 0;
  g_center_snp_end = g_flanking_snps-1;
  g_right_marker = 2*g_flanking_snps-1;

}

void GuidedMendelGPU::load_datasets(){
  cerr<<"Loading dataset for guided haplotyper\n";
  io_manager->read_input(ref_haplotype,g_snp_penetrance,this->informative_snp,this->g_people,this->g_snps, this->ref_haplotypes);
}

GuidedMendelGPU::~GuidedMendelGPU(){
  cerr<<"Entering destructor guided haplotyper\n";
  delete[] g_haplotype ;
  delete[] packedhap;
  delete[] packedextendedhap;
  delete[] extended_snp_mapping;
  delete[] extended_haplotype;
  delete[] extended_frequency;
  delete[] extended_root_mapping;
  cerr<<"Exiting destructor guided haplotyper\n";
}

void GuidedMendelGPU::finalize_window(){
  if (g_center_snp_end>g_flanking_snps-1 || g_center_snp_end+g_flanking_snps==g_right_marker){
    g_center_snp_start+=g_flanking_snps;
    g_center_snp_end+=g_flanking_snps;
    if (g_center_snp_end>=g_snps) g_center_snp_end = g_snps-1;
  }
  if (g_left_marker>0 || g_right_marker-g_left_marker+1>=g_max_window){
    g_left_marker+=g_flanking_snps;
  }
  if (g_right_marker<g_snps-1) 
  g_right_marker+=g_flanking_snps;
  if (g_right_marker>=g_snps) g_right_marker = g_snps-1;
}

void GuidedMendelGPU::init_window(){
  MendelGPU::init_window();
  copy_ref_haplotypes(g_left_marker);
  if (run_gpu){
    init_window_opencl();
  }
  //  runKernel("kernel_compute_weights",kernel_compute_weights,g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
//      float subject_haplotype_weight[g_people*g_max_haplotypes];
//      readFromBuffer(buffer_subject_haplotype_weight, g_people*g_max_haplotypes,subject_haplotype_weight,"buffer_subject_haplotype_weight");
//      float haplotype_weight[g_max_haplotypes];
//      memset(haplotype_weight,0,sizeof(float)*g_max_haplotypes);
//      for(int i = 0;i<g_people;++i){
//        for(int j=0;j<g_max_haplotypes;++j){
//          if (g_active_haplotype[j]){
//            cout<<"GPU person "<<i<<" hap "<<j<<" weight "<<subject_haplotype_weight[i*g_max_haplotypes+j]<<endl;
//            haplotype_weight[j]+=subject_haplotype_weight[i*g_max_haplotypes+j];
//          }
//        }
//      }
}


void GuidedMendelGPU::decompress(int hapindex,int markers,int * haplotype){
  //cerr<<"Extracting hapindex: "<<hapindex<<endl;
  for(int site = 0;site<markers;++site){
    int majorindex = site/32;
    int minorindex = (site%32)/8;
    int bitindex = (site%32)%8;
    //cerr<<"extracting val for index: "<<site<<": "<<majorindex<<","<<minorindex<<","<<bitindex<<endl;
    haplotype[hapindex*g_max_window+site] = ((int)packedhap[hapindex*packedhap_len+majorindex].octet[minorindex]) >> bitindex & 1;
  }
}

void GuidedMendelGPU::compress(int hapindex,int markers,int * haplotype){
  //cerr<<"Compressing hapindex: "<<hapindex<<endl;
  for(int site = 0;site<markers;++site){
    int majorindex = site/32;
    int minorindex = (site%32)/8;
    int bitindex = (site%32)%8;
    //cerr<<"inserting "<<i<<" for index: "<<index<<": "<<majorindex<<","<<minorindex<<","<<bitindex<<endl;
    if (bitindex==0) packedhap[hapindex*packedhap_len+majorindex].octet[minorindex] = 0;
    packedhap[hapindex*packedhap_len+majorindex].octet[minorindex] |= (haplotype[hapindex*g_max_window+site]<<bitindex);
  }
}

void GuidedMendelGPU::compress_extendedhap(int extendedhapindex,int begin,int end,int * extendedhaplotype){
  for(int site = begin;site<end;++site){
    int packedsite = site-begin;
    int majorindex = packedsite/32;
    int minorindex = (packedsite%32)/8;
    int bitindex = (packedsite%32)%8;
    if (bitindex==0) packedextendedhap[extendedhapindex*packedextendedhap_len+majorindex].octet[minorindex] = 0;
    packedextendedhap[extendedhapindex*packedextendedhap_len+majorindex].octet[minorindex] |= (extendedhaplotype[extendedhapindex*g_max_window+site]<<bitindex);
  }
}


void GuidedMendelGPU::parse_ref_haplotype(){
//  string line;
//  ifstream ifs(infile_refhap.data());
//  if(!ifs.is_open()){
//    cerr <<"Cannot open ref haplotypes file "<<infile_refhap<<endl;
//    exit(1);
//  }
//  getline(ifs,line);
//  ref_haplotypes = line.length();
//  ifs.close();
//  cerr<<"Parsing "<<ref_haplotypes<<" ref haplotypes across "<<g_snps<<" SNPs.\n";
//  ref_haplotype = new char[ref_haplotypes * g_snps];
//  ifs.open(infile_refhap.data());
//  for(int i=0;i<g_snps;++i){
//    getline(ifs,line);
//    for(int j=0;j<ref_haplotypes;++j){
//      ref_haplotype[j*g_snps+i] = line[j];
//    }
//  }
//  ifs.close();
}

void GuidedMendelGPU::copy_ref_haplotypes(int left_marker){
  bool debug_haplotype = false;
  //int left_marker = g_left_marker[0]-1;
  int marker_len = g_markers;
  // save the original # of markers from MENDEL
  extended_markers = marker_len;
  cerr<<"Including guide haplotypes of length: "<<marker_len<<" at left marker: "<<left_marker<<endl;
  //for(int j=0;j<100;++j){
    //cerr<<"debug load "<<j<<","<<informative_snp[j];
  //}
  //cerr<<endl;
  //cerr<<"got here copy 2\n";
  // generate a hash table of unique template haplotypes
  std::tr1::unordered_map<string,set<string> > temp_hap_window;
  //map<string,set<string> > temp_hap_window;
  for(int i=0;i<g_max_haplotypes;++i) g_active_haplotype[i] = 0;
  g_informative_markers = 0;
  for(int j=0;j<marker_len;++j){
    //cerr<<"Querying informative snp at pos "<<left_marker+j<<" and has value "<<informative_snp[left_marker+j]<<endl;
     if (informative_snp[left_marker+j]){
      extended_snp_mapping[g_informative_markers] = j;
      ++g_informative_markers;
    }
  }
  // notify MENDEL that the # of markers is temporarily reduced
  g_markers = g_informative_markers;
  //cerr<<"Informative markers: "<<g_markers<<endl;
  // tally up the counts of each reference haplotype
  ref_hap_counts.clear();
  for(int i=0;i<ref_haplotypes;++i){
    string curhap_long(ref_haplotype+i*g_snps+left_marker,marker_len);
    //cerr<<"ref hap of "<<curhap_long<<endl;
    if (ref_hap_counts.find(curhap_long)==ref_hap_counts.end()){
      ref_hap_counts[curhap_long] = 0;
    }
    ++ref_hap_counts[curhap_long];
  }
  cerr<<"Tallied up counts for each ref haplotype group. There are "<<ref_hap_counts.size()<<" groups.\n";
  // rank the template haplotypes by their frequencies
  multiset<hapobj_t,byHapCountDesc> sorted_guides_extended;
  for(std::tr1::unordered_map<string,int >::iterator it =
  ref_hap_counts.begin();it!=ref_hap_counts.end();it++){
    hapobj_t hapobj;
    hapobj.hapstr = it->first; 
    hapobj.hapcount = it->second;
    sorted_guides_extended.insert(hapobj);
  }
  cerr<<"Sorted haplotype groups by frequency. There are "<<sorted_guides_extended.size()<<" sorted groups. \n";
  int extended_count = 0;
  cerr<<"There are "<<g_informative_markers<<" informative markers\n";
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
  //cerr<<"Sorted haplotype groups by frequency 2\n";
  // rank the template haplotypes by their frequencies
  multiset<hapobj_t,byHapCountDesc> sorted_guides;
  for(std::tr1::unordered_map<string,set<string> >::iterator it = temp_hap_window.begin();it!=temp_hap_window.end();it++){
    hapobj_t hapobj;
    hapobj.hapstr = it->first; 
    hapobj.hapcount = it->second.size();
    hapobj.extended_set = it->second;
    sorted_guides.insert(hapobj);
  }
  //cerr<<"Ranked haplotype groups by frequency\n";
  // ADD TO HAP WINDOW THE DENOVO HAPLOTYPES
  int totalcounts = 0;
  cerr<<"Max haplotypes was set to "<<g_max_haplotypes<<endl;
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
      //cerr<<"Copy into ref haplotype "<<compact_haplotypes<<":"<<endl;
      for(int j=0;j<g_informative_markers;++j){
        g_haplotype[compact_haplotypes *g_max_window+j] = hapchar2int(hapobj.hapstr[j]);
        //cerr<<g_haplotype[compact_haplotypes *g_max_window+j];
      }
      //cerr<<endl;
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
  cerr<<"Full hap window size "<<full_hap_window.size()<<endl;
  // UPDATE TOTAL HAPLOTYPE SIZE
  g_haplotypes = compact_haplotypes;
  cerr<<"Working haplotypes after adding templates: "<<compact_haplotypes<<endl;
  cerr<<"Extended haplotypes: "<<extended_haplotypes<<endl;
  if (debug_haplotype) cerr<<"STORING "<<g_haplotypes<<
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

