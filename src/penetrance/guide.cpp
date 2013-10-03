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


void GuidedMendelGPU::compute_penetrance(){
  precompute_penetrance();
}

void GuidedMendelGPU::precompute_penetrance(){
  int left_marker = g_left_marker;
  bool debug_penmat = left_marker<-10;
  cerr<<"Entering precompute_penetrance at left marker "<<left_marker<<"\n";
  int debug_subject = 200;
  for(int hapindex=0;hapindex<g_haplotypes;++hapindex){
    compress(hapindex,g_markers,g_haplotype);
  }
  if (run_gpu){
#ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotypes, CL_TRUE, 0, sizeof(int),&g_haplotypes,NULL,NULL);
    clSafe(err, "write total haplotypes");
    err = commandQueue->enqueueWriteBuffer(*buffer_packedhap, CL_TRUE, 0,sizeof(packedhap_t)*g_max_haplotypes*packedhap_len, packedhap, NULL, NULL );
    clSafe(err, "write packedhaplotypes");
    err = commandQueue->enqueueWriteBuffer(*buffer_extended_snp_mapping, CL_TRUE, 0, sizeof(int)*g_max_window,extended_snp_mapping,NULL,NULL);
    clSafe(err, "write extended_snp_mapping");
    if (geno_dim==PHASED_INPUT){
      err = commandQueue->enqueueNDRangeKernel(*kernel_precompute_penetrance_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch precompute penetrance haploid");
    }else{
      err = commandQueue->enqueueNDRangeKernel(*kernel_precompute_penetrance,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch precompute penetrance");
    }
    if(debug_penmat){
      if (geno_dim==PHASED_INPUT){
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
      //cerr<<"Haploid status for person "<<i<<" is "<<haploid<<endl;
      float maxlog = -1000;
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j]){
          decompress(j,g_markers,g_haplotype);
          if (geno_dim==PHASED_INPUT){
            float logpenetrance[2] = {0,0};
            for(int l=0;l<g_markers;++l){
              //int dose2 = g_haplotype[j*g_max_window+l];
              bool dose = (((int)packedhap[j*packedhap_len+
              (l/32)].octet[(l%32)/8]) >> 
              ((l%32)%8) & 1) ;
              //cerr<<"dose1,2: "<<dose<<","<<dose2<<endl;
              for(int parent=0;parent<2;++parent){
                logpenetrance[parent]+=g_snp_penetrance[i*
                geno_dim*g_snps+geno_dim*(left_marker+
                extended_snp_mapping[l])+2*parent+dose]; 
                cerr<<"Hap "<<j<<" Person "<<i<<" parent "<<parent<<" SNP "<<left_marker+extended_snp_mapping[l]<<" dose "<<dose<<" value "<<
                g_snp_penetrance[i*
                geno_dim*g_snps+geno_dim*(left_marker+
                extended_snp_mapping[l])+2*parent+dose]<<
                endl; 
                
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
          }else{ // UNPHASED INPUT
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k] && (!haploid || j==k)){
                //cerr<<"max_window: "<<g_max_window<<" got here for "<<i<<","<<j<<","<<k<<endl;
                decompress(k,g_markers,g_haplotype);
                float logpenetrance = 0;
                for(int l=0;l<g_markers;++l){
                  int n=g_haplotype[j*g_max_window+l] + g_haplotype[k*g_max_window+l];
//cerr<<"haplogenotype for person "<<i<<" hap pairs "<<j<<","<<l<<" : "<<n<<endl;
//cerr<<"extendedsnpmapping l :"<<extended_snp_mapping[l]<<endl;
                  logpenetrance+=g_snp_penetrance[i*(geno_dim*g_snps)+geno_dim*(left_marker+extended_snp_mapping[l]) +  n]; 
//cerr<<"log penetrance is now "<<logpenetrance<<endl;
           
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
          if (geno_dim==PHASED_INPUT){
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
      if (geno_dim==PHASED_INPUT){
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

