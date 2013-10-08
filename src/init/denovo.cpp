#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"

int DenovoMendelGPU::get_max_window_size(){
  return 2 * g_flanking_snps+1;
}

void DenovoMendelGPU::allocate_memory(){
  MendelGPU::allocate_memory();
  cerr<<"Initializing variables for denovo haplotyper\n";
  twin_hap_index = new int[g_max_haplotypes];
  beyond_left_edge_dosage = new int[g_max_haplotypes];
  right_edge_dosage = new int[g_max_haplotypes];
  center_dosage = new int[g_max_haplotypes];

  //g_haplotype_mode = HAPLOTYPE_MODE_DENOVO;
  occupied_hap_indices.clear();
  free_hap_indices.clear();
  occupied_hap_indices.push_back(0);
  free_hap_indices.push_back(1);
  twin_hap_index[0] = 0;
  twin_hap_index[1] = 1;

  g_left_marker = 0;
  g_center_snp_start = 0;
  g_center_snp_end = 0;
  g_right_marker = 0;
}

DenovoMendelGPU::DenovoMendelGPU(IO_manager * io):MendelGPU(io){
}

DenovoMendelGPU::~DenovoMendelGPU(){
  cerr<<"Entering destructor denovo haplotyper\n";
  delete[] g_haplotype ;
  delete[] twin_hap_index;
  delete[] beyond_left_edge_dosage;
  delete[] right_edge_dosage;
  delete[] center_dosage;
  cerr<<"Exiting destructor denovo haplotyper\n";
}


void DenovoMendelGPU::init_window(){
  double_haplotypes();
  polymorphic_window = true;
  char zero = '0';
  int offset = (int)zero;
  //g_haplotype = haplotype;
  bool debug_haplotype = false;
  if (debug_haplotype){
    if (debug_mode) ofs_debug_haplotype_file<<"ASSUMED DOUBLED HAPLOTYPES:\n";
    for(int i=0;i<g_max_haplotypes;++i){
      if (g_active_haplotype[i]){
        if (debug_mode) ofs_debug_haplotype_file<<i<<":";
        string curhap(g_markers,zero);
        for(int j=0;j<g_markers;++j){
          curhap[j] = (char)(g_haplotype[i*g_max_window+j]+offset); 
          if (debug_mode) ofs_debug_haplotype_file<<curhap[j];
        }
        if (debug_mode) ofs_debug_haplotype_file<<endl;
      }
    }
  }
  if (run_gpu){
    init_window_opencl();
  }
  cerr<<"Max_window: "<<g_max_window<<" Haplotypes: "<<g_haplotypes<<" Markers: "<<g_markers<<" left_marker: "<<g_left_marker<<endl;
  return;

}

void DenovoMendelGPU::finalize_window(){
  cerr<<"Entering finalize window for denovo\n";
  prune_haplotypes_();
  g_prev_left_marker = g_left_marker;
  if (g_center_snp_end>0 || g_center_snp_end+g_flanking_snps==g_right_marker){
    ++g_center_snp_start;
    ++g_center_snp_end;
  }
  if (g_left_marker>0 || g_right_marker-g_left_marker+1==g_max_window){
    ++g_left_marker;
  }
  if (g_right_marker<g_snps-1) ++g_right_marker;
}

void DenovoMendelGPU::double_haplotypes(){
  cerr<<"Entering double_haplotypes with "<<g_haplotypes<<" haplotypes\n";
  cerr<<"Prev and current left marker: "<<g_prev_left_marker<<","<<g_left_marker<<endl;
  int last_marker = g_markers-1;
  cerr<<"Last marker is "<<last_marker<<endl;
  list<int>::iterator it_occupied = occupied_hap_indices.begin();
  list<int>::iterator it_free = free_hap_indices.begin();
  cerr<<"Total occupied: "<<occupied_hap_indices.size()<<" and free: "<<free_hap_indices.size()<<endl;
  std::tr1::unordered_set<string> seen_haplo;
  g_haplotypes=0;
  for(int i=0;i<g_max_haplotypes;++i) g_active_haplotype[i] = 0;
  while(it_occupied!=occupied_hap_indices.end()){
    int occupied_index = *it_occupied;
    int free_index = *it_free;
    //cerr<<"Current occupied/free index: "<<occupied_index<<","<<free_index<<endl;
    if (g_prev_left_marker!=g_left_marker){
      // shift hap alleles one over if we just slided
      for(int i=0;i<last_marker;++i){
        g_haplotype[occupied_index*g_max_window+i] = 
        g_haplotype[occupied_index*g_max_window+i+1] ;
      }
    }
    for(int i=0;i<last_marker;++i){
      g_haplotype[free_index*g_max_window+i] = 
      g_haplotype[occupied_index*g_max_window+i];
      g_frequency[free_index] = g_frequency[occupied_index];
    }
    g_haplotype[occupied_index*g_max_window+last_marker] = 0;
    g_haplotype[free_index*g_max_window+last_marker] = 1;
    int len = g_markers>1?g_markers-1:g_markers;
    string h = hapstr(g_haplotype+occupied_index*g_max_window,len);
    //cerr<<"occupied haplotype: "<<h<<endl;
    if(seen_haplo.find(h)==seen_haplo.end()){
      //cerr<<"Inserting\n";
      seen_haplo.insert(h);
      g_active_haplotype[occupied_index] = g_active_haplotype[free_index] = 1;
      g_haplotypes+=2;
    }else{
      g_active_haplotype[occupied_index] = g_active_haplotype[free_index] = 0;
    }
    it_occupied++;
    it_free++;
  }
  cerr<<"Leaving double_haplotypes, total: "<<g_haplotypes<<"\n";
  cerr<<"Active indices:";
  for(int i=0;i<g_max_haplotypes;++i){
    cerr<<" "<<g_active_haplotype[i];
  }
  cerr<<endl;
}

void DenovoMendelGPU::prune_haplotypes_(){
  cerr<<"Entering prune haplotypes\n";
  int half_max = g_max_haplotypes/2;
  int non_polymorphic_min = 2;
  int threshold = 
  (g_likelihood_mode==Config::LIKELIHOOD_MODE_READS && !polymorphic_window)?
  non_polymorphic_min : half_max;
  if (g_haplotypes>threshold){
    int marker_len = g_markers;
    multiset<hapobj_t,byHapFreqAsc> sorted_guides;
    for(int i=0;i<g_max_haplotypes;++i){
      if (g_active_haplotype[i]){
        hapobj_t hapobj;
        string curhap(marker_len,'0');
        for(int j=0;j<marker_len;++j){
          curhap[j] = hapint2char(g_haplotype[i*g_max_window+j]);
        }
        hapobj.hapstr = curhap;
        hapobj.hapfreq = g_frequency[i];
        hapobj.array_index = i;
        sorted_guides.insert(hapobj);
      }
    }
    cerr<<"Inserted haplotypes into sorted set\n"; 
    int haplotypes2delete = g_haplotypes - threshold;
    cerr<<"Need to remove "<<haplotypes2delete<<" haplotypes\n";
    cerr<<"Sorted guides has "<<sorted_guides.size()<<" elements\n";
    multiset<hapobj_t,byHapFreqAsc>::iterator it = sorted_guides.begin();
    for(int i=0;i<haplotypes2delete;++i){
      //cerr<<"Pruning to 1/2 by deleting haplotype index "<<it->array_index<<" with freq "<<it->hapfreq<<endl;
      g_active_haplotype[it->array_index] = 0;
      --g_haplotypes;
      it++;
    }
    while(it!=sorted_guides.end()){
      if (it->hapfreq<1e-8){
        //cerr<<"Deleting rare haplotype index "<<it->array_index<<" with freq "<<it->hapfreq<<endl;
        g_active_haplotype[it->array_index] = 0;
        --g_haplotypes;
      }
      it++;
    }
    cerr<<"Setting total haplotypes to "<<g_haplotypes<<endl;
  }
  // create the linked lists of ordered occupied and available slots
  occupied_hap_indices.clear();
  free_hap_indices.clear();
  for(int i=0;i<g_max_haplotypes;++i){
    if (g_active_haplotype[i]){
      occupied_hap_indices.push_back(i);
    }else{
      free_hap_indices.push_back(i);
    }
  }
  // at this point persists alleles from the left edge 
  // so they can be retrieved later when this site has expired.
  list<int>::iterator it_occupied = occupied_hap_indices.begin();
  list<int>::iterator it_free = free_hap_indices.begin();
  while(it_occupied!=occupied_hap_indices.end()){
    int occupied_index = *it_occupied;
    int free_index = *it_free;
    // these correspond to a pair where free will copy 
    // penetrances from occupied in the next window 
    twin_hap_index[occupied_index] = occupied_index;
    twin_hap_index[free_index] = occupied_index;
    beyond_left_edge_dosage[free_index] = 
    beyond_left_edge_dosage[occupied_index] = g_haplotype[occupied_index*g_max_window];
    //cerr<<"Setting left edge for "<<occupied_index<<" to "<<g_haplotype[occupied_index*g_max_window]<<endl;
    //cerr<<"free index "<<free_index<<" now points to "<<occupied_index<<endl;
    it_occupied++;
    it_free++;
  }
  cerr<<"Leaving prune haplotypes\n";
}


void DenovoMendelGPU::debug_haplotypes(ostream & os){
  for(int i=0;i<g_max_haplotypes;++i){
    if (g_active_haplotype[i]){
      os<<i<<":";
      for(int j=0;j<g_markers;++j){
        os<<g_haplotype[i*g_max_window+j];
      }
      os<<endl;
    }
  }
}
