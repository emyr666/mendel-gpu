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

using namespace std;

void MendelGPU::double_haplotypes_(){
  cerr<<"Entering double_haplotypes\n";
  cerr<<"Prev and current left marker: "<<g_prev_left_marker[0]<<","<<g_left_marker[0]<<endl;
  int last_marker = g_markers[0]-1;
  //cerr<<"Last marker is "<<last_marker<<endl;
  list<int>::iterator it_occupied = occupied_hap_indices.begin();
  list<int>::iterator it_free = free_hap_indices.begin();
  cerr<<"Total occupied: "<<occupied_hap_indices.size()<<" and free: "<<free_hap_indices.size()<<endl;
  std::tr1::unordered_set<string> seen_haplo;
  *g_haplotypes=0;
  for(int i=0;i<g_max_haplotypes;++i) g_active_haplotype[i] = 0;
  while(it_occupied!=occupied_hap_indices.end()){
    int occupied_index = *it_occupied;
    int free_index = *it_free;
    //cerr<<"Current occupied/free index: "<<occupied_index<<","<<free_index<<endl;
    if (g_prev_left_marker[0]!=g_left_marker[0]){
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
    int len = g_markers[0]>1?g_markers[0]-1:g_markers[0];
    string h = hapstr(g_haplotype+occupied_index*g_max_window,len);
    if(seen_haplo.find(h)==seen_haplo.end()){
      seen_haplo.insert(h);
      g_active_haplotype[occupied_index] = g_active_haplotype[free_index] = 1;
      *g_haplotypes+=2;
    }else{
      g_active_haplotype[occupied_index] = g_active_haplotype[free_index] = 0;
    }
    it_occupied++;
    it_free++;
  }
  cerr<<"Active indices:";
  for(int i=0;i<g_max_haplotypes;++i){
    cerr<<" "<<g_active_haplotype[i];
  }
  cerr<<endl;
  cerr<<"Leaving double_haplotypes, total: "<<*g_haplotypes<<"\n";
}

void MendelGPU::prune_haplotypes_(){
  cerr<<"Entering prune haplotypes\n";
  int half_max = g_max_haplotypes/2;
  if (g_haplotypes[0]>half_max){
    int marker_len = g_markers[0];
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
    int haplotypes2delete = *g_haplotypes - half_max;
    cerr<<"Need to remove "<<haplotypes2delete<<" haplotypes\n";
    cerr<<"Sorted guides has "<<sorted_guides.size()<<" elements\n";
    multiset<hapobj_t,byHapFreqAsc>::iterator it = sorted_guides.begin();
    for(int i=0;i<haplotypes2delete;++i){
      g_active_haplotype[it->array_index] = 0;
      --*g_haplotypes;
      //cerr<<"Pruning to 1/2 by deleting haplotype index "<<it->array_index<<" with freq "<<it->hapfreq<<endl;
      it++;
    }
    while(it!=sorted_guides.end()){
      if (it->hapfreq<1e-8){
        g_active_haplotype[it->array_index] = 0;
        --*g_haplotypes;
        //cerr<<"Deleting rare haplotype index "<<it->array_index<<" with freq "<<it->hapfreq<<endl;
      }
      it++;
    }
    cerr<<"Setting total haplotypes to "<<*g_haplotypes<<endl;
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

void MendelGPU::precompute_penetrance_fast_(){
  cerr<<"Entering precompute_penetrance_fast\n";
  int last_marker = g_markers[0]-1;
  int prev_left_marker = g_prev_left_marker[0]-1;
  int left_marker = g_left_marker[0]-1;
  cerr<<"Taking hap doses at site "<<last_marker<<" and prev and current left marker: "<<prev_left_marker<<","<<left_marker<<"\n";
  bool debug_penmat = false;
  //bool debug_penmat = g_markers[0]>4;
  int debug_penmat_person = 250;
  //bool debug_penmat = g_haplotypes[0]<-8;
  if (debug_penmat) cout<<"Debugging precompute_penetrance_fast_()\n";
  for(int j=0;j<g_max_haplotypes;++j){
    if(g_active_haplotype[j] ){
      right_edge_dosage[j] = g_haplotype[j*g_max_window+last_marker]; 
    }
  }
  if (run_gpu){
#ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_right_edge_dosage, CL_TRUE, 0, sizeof(int)*g_max_haplotypes,right_edge_dosage,NULL,NULL);
    clSafe(err, "write right edge dosage");
    if (prev_left_marker!=left_marker){
      err = commandQueue->enqueueWriteBuffer(*buffer_beyond_left_edge_dosage, CL_TRUE, 0, sizeof(int)*g_max_haplotypes,beyond_left_edge_dosage,NULL,NULL);
      clSafe(err, "write beyond left edge dosage");
    }
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotypes, CL_TRUE, 0, sizeof(int),g_haplotypes,NULL,NULL);
    clSafe(err, "write total haplotypes");
    if (geno_dim==4){
      cerr<<"Launching haploid fast precompute\n";
      err = commandQueue->enqueueNDRangeKernel(*kernel_precompute_penetrance_fast_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch precompute penetrance_fast");
    }else if(geno_dim==3){
      err = commandQueue->enqueueNDRangeKernel(*kernel_precompute_penetrance_fast,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch precompute penetrance_fast");
    }
    if(debug_penmat){
      if (geno_dim==4){
        float * debug_cache = new float[g_people*2*g_max_haplotypes];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*2*g_max_haplotypes,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"fast GPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<debug_cache[debug_penmat_person*g_max_haplotypes*2+j*2+k];
            }
            cout<<endl;
          }
        }
        delete[] debug_cache;
      }else if (geno_dim==3){ 
        float * debug_cache = new float[g_people*penetrance_matrix_size];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*penetrance_matrix_size,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"fast GPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<debug_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
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
      float maxlog = -1e10;
      for(int j=0;j<g_max_haplotypes;++j){
        if(g_active_haplotype[j] ){
          if (geno_dim==4){
            int dose= right_edge_dosage[j];
            float logpenetrance_new[2];
            float logpenetrance_old[2]={0,0};
            float logpenetrance[2];
            for(int parent=0;parent<2;++parent){
              logpenetrance_new[parent] = g_region_snp_penetrance[i*(geno_dim*
              g_max_region_size)+ geno_dim*(left_marker+last_marker-
              g_region_snp_offset) + 2*parent+dose];
      //        cerr<<"Person "<<i<<" parent "<<parent<<" SNP "<<left_marker+last_marker<<" dose "<<dose<<" value "<<
       //       logpenetrance_new[parent]<<endl;
            }
            if (prev_left_marker!=left_marker){
              dose=beyond_left_edge_dosage[j];
              for(int parent=0;parent<2;++parent){
                logpenetrance_old[parent] = g_region_snp_penetrance[i*(geno_dim*
                g_max_region_size)+ geno_dim*(prev_left_marker-
                g_region_snp_offset) + 2*parent+dose];
              }
            }
            for(int parent=0;parent<2;++parent){
              logpenetrance[parent] = g_haplotypes[0]>2?logpenetrance_cache[i*
              2*g_max_haplotypes+2*j+parent] + logpenetrance_new[parent]-
              logpenetrance_old[parent]: logpenetrance_new[parent];
              if (logpenetrance[parent] > maxlog) 
              maxlog = logpenetrance[parent];
              logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              logpenetrance[parent]; 
            }
          }else{
            for(int k=j;k<g_max_haplotypes;++k){
              if(g_active_haplotype[k]){
                int dose= right_edge_dosage[j] + right_edge_dosage[k];
                float logpenetrance_new = g_region_snp_penetrance[i*(geno_dim*
                g_max_region_size)+ geno_dim*(left_marker+last_marker-
                g_region_snp_offset) + dose];
                float logpenetrance_old = 0;
                if (prev_left_marker!=left_marker){
                  dose=beyond_left_edge_dosage[j] + beyond_left_edge_dosage[k];
                  logpenetrance_old = g_region_snp_penetrance[i*(geno_dim*
                  g_max_region_size)+ geno_dim*(prev_left_marker-
                  g_region_snp_offset) + dose];
                }
                float logpenetrance = g_haplotypes[0]>2?logpenetrance_cache[i*
                penetrance_matrix_size+j*g_max_haplotypes+k] + 
                logpenetrance_new-logpenetrance_old:logpenetrance_new;
                if (logpenetrance > maxlog) maxlog = logpenetrance;
                logpenetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k] = logpenetrance; 
              }
            }
          }
        } 
      }
//maxlog = 0;
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] ){
          if (geno_dim==4){
            for(int parent=0;parent<2;++parent){
              float val = logpenetrance_cache[i*2*g_max_haplotypes+2*
              j+parent]-maxlog;
              logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent] = val;
              penetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              val>=logpenetrance_threshold?exp(val):0;
              //if (g_haplotypes[0]>6 && i==0) cerr<<"Log penetrance/penetrance for person "<<i<<" hap "<<j<<" parent "<<parent<<": "<<val<<"/"<<penetrance_cache[i*2*g_max_haplotypes+2*j+parent] <<endl;
            }
          }else{
            for(int k=j;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                float val = logpenetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k]-maxlog;
                logpenetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k] = val;
                logpenetrance_cache[i*penetrance_matrix_size+k*
                g_max_haplotypes+j] = val;
                penetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k] = val>=logpenetrance_threshold?exp(val):0;
                penetrance_cache[i*penetrance_matrix_size+k*
                g_max_haplotypes+j] = penetrance_cache[i*
                penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
          }
        }   
      }
    }
    if(debug_penmat){
      if (geno_dim==4){
        cout<<"fast CPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<logpenetrance_cache[debug_penmat_person*g_max_haplotypes*2+j*2+k];
            }
            cout<<endl;
          }
        }
      }else if (geno_dim==3){ 
        cout<<"fast CPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<logpenetrance_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
      }
    }
    cerr<<"Precompute penetrance fast time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
  }
  if (debug_penmat) exit(1);
  cerr<<"Exiting precompute_penetrance_fast\n";
}

void MendelGPU::impute_penetrance_matrix_(){
  int left_marker=g_left_marker[0]-1;
  bool debug_basic = false;
  bool debug_penmat = false;
  int debug_penmat_person = 600;
  cerr<<"Entering impute penetrance matrix\n";
  if (run_gpu){
#ifdef USE_GPU
    if (debug_basic){
      cerr<<"Twin indices:";
      for(int i=0;i<g_max_haplotypes;++i){
        g_active_haplotype[i] = (i%2==0);
        twin_hap_index[i] = i<g_max_haplotypes/2 ? i+g_max_haplotypes/2 : i;
        cerr<<" "<<i<<":"<<twin_hap_index[i];
      }
      cerr<<endl;
      for(int i=0;i<g_people;++i){
        for(int j=0;j<g_max_haplotypes;++j){
          for(int k=0;k<g_max_haplotypes;++k){
            logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = (j>=g_max_haplotypes/2 && k>=g_max_haplotypes/2)? j*g_max_haplotypes+k : 0;      
          }
        }
      }
      err = commandQueue->enqueueWriteBuffer(*buffer_active_haplotype, CL_TRUE, 0, sizeof(int)*g_max_haplotypes,g_active_haplotype,NULL,NULL);
      clSafe(err, "write active haplotype");
      err = commandQueue->enqueueWriteBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people * penetrance_matrix_size ,logpenetrance_cache,NULL,NULL);
      clSafe(err, "write logpenetrance cache");
    }
    err = commandQueue->enqueueWriteBuffer(*buffer_twin_hap_index, CL_TRUE, 0, sizeof(int)*g_max_haplotypes,twin_hap_index,NULL,NULL);
    clSafe(err, "write twin indices");
    if (geno_dim==4){
      err = commandQueue->enqueueNDRangeKernel(*kernel_impute_penetrance_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch impute penetrance_haploid");
    }else if(geno_dim==3){
      err = commandQueue->enqueueNDRangeKernel(*kernel_impute_penetrance,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch impute penetrance");
    }
    if (debug_penmat){
      if (geno_dim==4){
        float * debug_cache = new float[g_people*2*g_max_haplotypes];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*2*g_max_haplotypes,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"GPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<debug_cache[debug_penmat_person*g_max_haplotypes*2+j*2+k];
            }
            cout<<endl;
          }
        }
        delete[] debug_cache;
      }else if (geno_dim==3){
        float * debug_cache = new float[g_people*penetrance_matrix_size];
        err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*penetrance_matrix_size,debug_cache);
        clSafe(err, "read penetrance cache");
        cout<<"GPU for person "<<debug_penmat_person<<":\n";
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<debug_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
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
  if (run_cpu){
    for(int i=0;i<g_people;++i){
      if (geno_dim==3){  // for diplid case, columns matter
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            for(int k=0;k<g_max_haplotypes;++k){ //copy col
              if (g_active_haplotype[k] && twin_hap_index[k]!=k){
                logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+twin_hap_index[k]]; 
              }
            }
          }
        }
      }
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] && twin_hap_index[j]!=j){
          if (geno_dim==4){
            // copy penetrance for maternal and paternal haplotypes
            for(int parent=0;parent<2;++parent){
              logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent] = 
              logpenetrance_cache[i*2*g_max_haplotypes+2*twin_hap_index[j]
              +parent]; 
              //cerr<<"Copied from "<<twin_hap_index[j]<<" to "<<j<<" value "<<logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent]<<endl;
            }
          }else{
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                logpenetrance_cache[i*penetrance_matrix_size+j*
                g_max_haplotypes+k] = logpenetrance_cache[i*
                penetrance_matrix_size+twin_hap_index[j]*g_max_haplotypes+k]; 
              }
            }
          }
        }
      }
    }
    if(debug_penmat){
      cout<<"CPU for person "<<debug_penmat_person<<":\n";
      if (geno_dim==4){
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<2;++k){
              cout<<" "<<logpenetrance_cache[debug_penmat_person*g_max_haplotypes*2+j*2+k];
            }
            cout<<endl;
          }
        }
      }else if (geno_dim==3){
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<j<<":";
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k]){
                cout<<" "<<logpenetrance_cache[debug_penmat_person*penetrance_matrix_size+j*g_max_haplotypes+k];
              }
            }
            cout<<endl;
          }
        }
      }
    }
  }
  cerr<<"Leaving impute penetrance matrix\n";
  if (debug_penmat) exit(1);
}

void MendelGPU::impute_haploid_genotypes_denovo_(int * center_snp, int * d_snp){
  int current_snp = *d_snp-1;
  int c_snp = *center_snp-1;
  cerr<<"begin impute_genotypes at SNP "<<current_snp<<", center: "<<c_snp<<endl;
  bool debug_dosage = false;
  bool debug_geno = false;
  bool debug_posterior = true;
  bool debug_pen = current_snp == -1;
  int max_geno=g_genotype_imputation?3:4;
  //debug_haplotypes(cerr);
  for(int i=0;i<g_max_haplotypes;++i){
    center_dosage[i] = g_haplotype[i*g_max_window+c_snp];
  }
  float subject_dosages[g_people];
  if (run_gpu){
    float subject_posterior[4*g_people];
    int subject_genotypes[g_people];
    double start = clock();
    #ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_center_dosage, CL_TRUE, 0,sizeof(int)*g_max_haplotypes, center_dosage, NULL, NULL );
    clSafe(err, "write center snp");
    err = commandQueue->enqueueNDRangeKernel(*kernel_impute_genotype_denovo_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
    clSafe(err,"launch impute_genotype_denovo");
    cerr<<"Launched impute_genotype_denovo\n";
    err = commandQueue->enqueueReadBuffer(*buffer_subject_posterior_prob, CL_TRUE, 0, sizeof(float)*g_people*4,subject_posterior);
    clSafe(err, "read subject_posterior");

    for(int i=0;i<g_people;++i){
      float * posterior_prob = subject_posterior + i*4;
      for(int k=0;k<4;++k) posterior_prob[k] = 
      posterior_prob[k]<epsilon?epsilon:posterior_prob[k];
      float denom[2] = {0,0};
      for(int parent=0;parent<2;++parent){
        for(int allele=0;allele<2;++allele){
          denom[parent] += posterior_prob[2*parent+allele];
        }
        for(int allele=0;allele<2;++allele){
          posterior_prob[2*parent+allele]/=denom[parent];
        }
      }
    }
    if (debug_posterior){
      for(int i=0;i<g_people;++i){
        cout<<"GPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
        for(int parent=0;parent<2;++parent){
          for(int allele=0;allele<2;++allele){
            cout<<"\t"<<subject_posterior[i*4+parent*2+allele];
          }
        }
        cout<<endl;
      }
    }
    if (outfile_format.compare(FORMAT_MEC)==0){
      ofs_posterior_file<<current_snp;
      for(int i=0;i<g_people;++i){
        for(int j=0;j<max_geno;++j){
          ofs_posterior_file<<"\t"<<subject_posterior[i*4+j];
        }
      }
      ofs_posterior_file<<endl;
    }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
      for(int i=0;i<g_people;++i){
        ofs_posterior_file<<"GPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
        for(int j=0;j<max_geno;++j){
          ofs_posterior_file<<"\t"<<subject_posterior[i*4+j];
        }
        ofs_posterior_file<<endl;
      }
    }
    cerr<<"Elapsed time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    #endif
  }
  if(run_cpu){
    if (outfile_format.compare(FORMAT_MEC)==0){
      ofs_dosage_file<<current_snp<<"\t";
      ofs_genotype_file<<current_snp<<"\t";
      ofs_posterior_file<<current_snp;
    }
    for(int i=0;i<g_people;++i){
      float posterior_prob[4];
      memset(posterior_prob,0,sizeof(float)*4);
      for(int j=0;j<g_max_haplotypes;++j){
        if(g_active_haplotype[j] ){
          int m;
          for(int parent=0;parent<2;++parent){
            float penetrance = penetrance_cache[i*2*g_max_haplotypes+
            2*j+parent];
            if (penetrance>0){
              float p = penetrance * g_frequency[j];
              posterior_prob[2*parent+center_dosage[j]]+=p;
            }
          }
        }
      }
      for(int k=0;k<4;++k) posterior_prob[k] = posterior_prob[k]<epsilon?epsilon:posterior_prob[k];
      float denom[2] = {0,0};
      float dose = 0;
      for(int parent=0;parent<2;++parent){
        for(int allele=0;allele<2;++allele){
          denom[parent] += posterior_prob[2*parent+allele];
        }
        for(int allele=0;allele<2;++allele){
          posterior_prob[2*parent+allele]/=denom[parent];
        }
      }
      if (debug_posterior){
        cout<<"CPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
        for(int parent=0;parent<2;++parent){
          for(int allele=0;allele<2;++allele){
            float p = posterior_prob[2*parent+allele];
            cout<<"\t"<<p;
          }
        }
        cout<<endl;
      }
      //for(int k=0;k<4;++k) cerr<<" "<<posterior_prob[k]<<endl;
      //for(int k=0;k<2;++k) cerr<<"DENOM: "<<denom[k]<<endl;
      if (outfile_format.compare(FORMAT_MEC)==0){
        for(int parent=0;parent<2;++parent){
          for(int allele=0;allele<2;++allele){
            float p = posterior_prob[2*parent+allele];
            ofs_posterior_file<<"\t"<<p;
            dose+=allele*p;
          }
        }
        char chardose =  (int)(dose*10)+'A';
        ofs_dosage_file<<chardose;
        //ofs_dosage_file.close();ofs_posterior_file.close();exit(0);
      }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
        ofs_posterior_file<<"CPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
        for(int parent=0;parent<2;++parent){
          for(int allele=0;allele<2;++allele){
            float p = posterior_prob[2*parent+allele];
            ofs_posterior_file<<"\t"<<p;
            dose+=allele*p;
          }
        }
        ofs_posterior_file<<endl;
        char chardose =  (int)(dose*10)+'A';
        ofs_dosage_file<<"CPU_DOSE:\t"<<i<<"\t"<<current_snp<<
        "\t"<<chardose<<endl;
      }
      subject_dosages[i] = dose;
      float genoprobs[3];
      genoprobs[0] = posterior_prob[0]*posterior_prob[2];
      genoprobs[1] = posterior_prob[1]*posterior_prob[2] + 
                     posterior_prob[0]*posterior_prob[3];
      genoprobs[2] = posterior_prob[1]*posterior_prob[3];
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
      ofs_dosage_file<<endl;
      ofs_genotype_file<<endl;
      ofs_posterior_file<<endl; 
    }
    float rsq = compute_rsq(subject_dosages,1,0);
    ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_posterior) exit(1);
}


void MendelGPU::impute_diploid_genotypes_denovo_(int * center_snp, int * d_snp){
  int current_snp = *d_snp-1;
  int c_snp = *center_snp-1;
  cerr<<"begin impute_genotypes at SNP "<<current_snp<<", center: "<<c_snp<<endl;
  bool debug_dosage = false;
  bool debug_geno = false;
  bool debug_posterior = false;
  bool debug_pen = current_snp == -1;
  int max_geno=g_genotype_imputation?3:4;
  for(int i=0;i<g_max_haplotypes;++i){
    center_dosage[i] = g_haplotype[i*g_max_window+c_snp];
  }
  float subject_dosages[g_people];
  if (run_gpu){
    float subject_posterior[4*g_people];
    int subject_genotypes[g_people];
    double start = clock();
    #ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_center_dosage, CL_TRUE, 0,sizeof(int)*g_max_haplotypes, center_dosage, NULL, NULL );
    clSafe(err, "write center snp");
    err = commandQueue->enqueueNDRangeKernel(*kernel_impute_genotype_denovo,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
    clSafe(err,"launch impute_genotype_denovo");
    cerr<<"Launched impute_genotype_denovo\n";
    err = commandQueue->enqueueReadBuffer(*buffer_subject_posterior_prob, CL_TRUE, 0, sizeof(float)*g_people*4,subject_posterior);
    clSafe(err, "read subject_posterior");
    if (outfile_format.compare(FORMAT_MEC)==0){
      ofs_posterior_file<<current_snp;
      for(int i=0;i<g_people;++i){
        for(int j=0;j<max_geno;++j){
          ofs_posterior_file<<"\t"<<subject_posterior[i*4+j];
        }
      }
      ofs_posterior_file<<endl;
    }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
      for(int i=0;i<g_people;++i){
        ofs_posterior_file<<"GPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
        for(int j=0;j<max_geno;++j){
          ofs_posterior_file<<"\t"<<subject_posterior[i*4+j];
        }
        ofs_posterior_file<<endl;
      }
    }
    if (debug_posterior){
      for(int i=0;i<g_people;++i){
        float denom = 0;
        for(int j=0;j<max_geno;++j){
          denom+=subject_posterior[i*4+j];
        }
        if(debug_mode) ofs_posterior_file<<"GPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
        for(int j=0;j<max_geno;++j){
          if(debug_mode) ofs_posterior_file<<"\t"<< subject_posterior[i*4+j]/denom;
        }
        if(debug_mode) ofs_posterior_file<<endl;
        if (denom<epsilon){
          //cerr<<"DENOM < EPSILON\n";
          if (denom==0){
             cerr<<"DENOM = ZERO\n";
             err = commandQueue->enqueueReadBuffer(*buffer_logpenetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*penetrance_matrix_size,logpenetrance_cache);
             clSafe(err, "read logpenetrance cache");
             err = commandQueue->enqueueReadBuffer(*buffer_penetrance_cache, CL_TRUE, 0, sizeof(float)*g_people*penetrance_matrix_size,penetrance_cache);
             clSafe(err, "read penetrance cache");
             cerr<<"dosage,j,k,logpen,pen,freq_j,freq_k"<<endl;
             for(int j=0;j<g_max_haplotypes;++j){
               for(int k=0;k<g_max_haplotypes;++k){
                 if (g_active_haplotype[j]&& g_active_haplotype[k]){
                   cerr<<"GPU:"<<center_dosage[j]+center_dosage[k]<<","<<j<<","<<k<<":"<<logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]<<","<<penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]<<","<<g_frequency[j]<<","<<g_frequency[k]<<endl;
                 }
               }
             }
             exit(1);
          }
        }
      }
    }
    err = commandQueue->enqueueReadBuffer(*buffer_subject_genotype, CL_TRUE, 0, sizeof(int)*g_people,subject_genotypes);
    clSafe(err, "read subject_genotype");
    if (debug_geno){
      for(int i=0;i<g_people;++i){
        cout<<"GPU GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_genotypes[i]<<endl;
      }
    }
    if (outfile_format.compare(FORMAT_MEC)==0){
      ofs_genotype_file<<current_snp<<"\t";
      for(int i=0;i<g_people;++i){
        ofs_genotype_file<<subject_genotypes[i];
      }
      ofs_genotype_file<<endl;
    }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
      for(int i=0;i<g_people;++i){
        ofs_genotype_file<<"GPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_genotypes[i]<<endl;
      }
    }
    if (g_genotype_imputation){
      err = commandQueue->enqueueReadBuffer(*buffer_subject_dosage, CL_TRUE, 0, sizeof(float)*g_people,subject_dosages);
      clSafe(err, "read subject_dosage");
      if (outfile_format.compare(FORMAT_MEC)==0){
        ofs_dosage_file<<current_snp<<"\t";
        for(int i=0;i<g_people;++i){
          char chardose = 
          (int)(subject_dosages[i]*10)+'A';
          ofs_dosage_file<<chardose;
        }
        ofs_dosage_file<<endl;
      }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
        for(int i=0;i<g_people;++i){
          ofs_dosage_file<<"GPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosages[i]<<endl;
          if (debug_dosage) cout<<"GPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosages[i]<<endl;
        }
      }
      float rsq = compute_rsq(subject_dosages,1,0);
      ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
    }
    cerr<<"Elapsed time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    #endif
  }
  if(run_cpu){
    if (outfile_format.compare(FORMAT_MEC)==0){
      ofs_dosage_file<<current_snp<<"\t";
      ofs_genotype_file<<current_snp<<"\t";
      ofs_posterior_file<<current_snp;
    }
    for(int i=0;i<g_people;++i){
      bool debug_truth = false;
      //bool debug_truth = (i==14) && (current_snp==0);
      float best_prob = 0;
      float posterior_prob[4];
      int best_pair[2];
      memset(posterior_prob,0,sizeof(float)*4);
      for(int j=0;j<g_max_haplotypes;++j){
        for(int k=j;k<g_max_haplotypes;++k){
          if(g_active_haplotype[j] && g_active_haplotype[k]){
            int m;
            float penetrance = penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
            //cerr<<"Penetrance at "<<j<<" and "<<k<<" :" <<penetrance<<endl;
            if (penetrance>0){
              float freq = g_frequency[j]*g_frequency[k];
              if (j!=k) freq*=2;
              float p = freq*penetrance;
              if (p>best_prob){
                best_prob = p;
                best_pair[0] = j;
                best_pair[1] = k;
              }
              if(g_genotype_imputation){
                m = center_dosage[j] + center_dosage[k];
              }else{
                m = 2*center_dosage[j] + center_dosage[k];
              }
              if (debug_pen) cerr<<"i,j,k,m,pen,freq:"<<i<<","<<j<<","<<k<<","<<m<<","<<penetrance<<","<<freq<<endl;
              posterior_prob[m]+=p;
            }
          }
        }
      }
      if (debug_posterior){
        cout<<"CPU POSTERIOR "<<current_snp<<","<<i;
        for(int j=0;j<max_geno;++j){
          cout<<" "<<posterior_prob[j];
        }
        cout<<endl;
      }
      float denom = 0;
      for(int m=0;m<max_geno;++m){
        //if (posterior_prob[m]<epsilon) posterior_prob[m] = epsilon;
        denom+=posterior_prob[m];
      }
      if (denom==0) {
        cerr<<"divide by zero denom\n"; 
        cerr<<"hap1,hap2:pen,freq,genotype\n";
        for(int j=0;j<g_max_haplotypes;++j){
          int start = g_genotype_imputation?j:0;
          for(int k=start;k<g_max_haplotypes;++k){
            if(g_active_haplotype[j] && g_active_haplotype[k]){
              int m = g_haplotype[j*g_max_window+ c_snp] + 
              g_haplotype[k*g_max_window+ c_snp];
              float logpenetrance = logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
              float freq = g_frequency[j]*g_frequency[k];
              cerr<<j<<","<<k<<":"<<logpenetrance<<","<<freq<<","<<m<<endl;
            }
          }
        }
        exit(1);
      }
      float dose = 0;
      if (outfile_format.compare(FORMAT_MEC)==0){
        for(int j=0;j<max_geno;++j){
          float p = posterior_prob[j]/denom;
          ofs_posterior_file<<"\t"<<p;
          dose+=j*p;
        }
        char chardose =  (int)(dose*10)+'A';
        ofs_dosage_file<<chardose;
      }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
        ofs_posterior_file<<"CPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
        for(int j=0;j<max_geno;++j){
          float p = posterior_prob[j]/denom;
          if(debug_mode)ofs_posterior_file<<"\t"<< p;
          dose+=j*p;
        }
        ofs_dosage_file<<"CPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<dose<<endl;
        ofs_posterior_file<<endl;
      }
      subject_dosages[i] = dose;
      float maxval = 0;
      for(int j=0;j<4;++j){
        if (posterior_prob[j]>maxval) maxval = posterior_prob[j];
      }
      if (maxval>0){
        int j = best_pair[0];
        int k = best_pair[1];
        //cerr<<"best prob "<<best_prob<<" "<<c_snp<<" "<<j<<" "<<k<<" "<<g_haplotype[j*g_max_window+ c_snp] <<" "<<g_haplotype[k*g_max_window+ c_snp]<<endl;
        int geno1 = g_genotype_imputation?
        g_haplotype[j*g_max_window+ c_snp] +
        g_haplotype[k*g_max_window+ c_snp]:
        2*g_haplotype[j*g_max_window+ c_snp] +
        g_haplotype[k*g_max_window+ c_snp];
        float max_prob = posterior_prob[0];
        int geno2 = 0;
        for(int h=1;h<4;++h){
          if (posterior_prob[h]-max_prob>epsilon*max_prob){
            max_prob = posterior_prob[h];
            geno2 = h;
          }
        }  
        if (debug_truth){
          int truegeno = true_geno[i*g_snps+current_snp];
          float truemaf = true_maf[current_snp];
          if ((truegeno!=geno2 || truegeno!=geno1) && truemaf<.01){
            cerr<<"Mismatch at person "<<i<<" SNP "<<current_snp<<endl;
            cerr<<"True geno, geno1, geno2: "<<truegeno<<","<<geno1<<","<<geno2<<endl;
            cerr<<"bestpair: "<<best_pair[0]<<","<<best_pair[1]<<endl;
            cerr<<"hap1,hap2:pen,freq,genotype\n";
            for(int j=0;j<g_max_haplotypes;++j){
              int start = g_genotype_imputation?j:0;
              for(int k=start;k<g_max_haplotypes;++k){
                if(g_active_haplotype[j] && g_active_haplotype[k]){
                  int m = g_haplotype[j*g_max_window+ c_snp] + 
                  g_haplotype[k*g_max_window+ c_snp];
                  float logpenetrance = logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
                  float freq = g_frequency[j]*g_frequency[k];
                  cerr<<j<<","<<k<<":"<<logpenetrance<<","<<freq<<","<<m<<endl;
                }
              }
            }
          }
        }
        if (debug_geno){
          cout<<"CPU GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<geno1<<","<<geno2<<endl;
        }
        ofs_genotype_file<<geno2;
        //ofs_genotype_file<<"CPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<geno1<<","<<geno2<<"\t"<<best_pair[0]<<","<<best_pair[1]<<endl;
      }
    }
    if (outfile_format.compare(FORMAT_MEC)==0){
      ofs_dosage_file<<endl;
      ofs_genotype_file<<endl;
      ofs_posterior_file<<endl; 
    }
    float rsq = compute_rsq(subject_dosages,1,0);
    ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_dosage || debug_posterior|| debug_geno) exit(1);
  if (debug_pen) {
    if (debug_mode) ofs_debug_haplotype_file.close();
    exit(1);
  }
}
