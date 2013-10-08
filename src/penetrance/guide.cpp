#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"


void GuidedMendelGPU::compute_penetrance(){
  precompute_penetrance();
}

void GuidedMendelGPU::precompute_penetrance(){
  bool debug_penmat = g_left_marker<0;
  cerr<<"Entering precompute_penetrance at left marker "<<g_left_marker<<"\n";
  int debug_subject = 1;
  for(int hapindex=0;hapindex<g_haplotypes;++hapindex){
    compress(hapindex,g_markers,g_haplotype);
  }
  if(run_gpu){
#ifdef USE_GPU
    precompute_penetrance_opencl();
#endif
  }
  if(run_cpu){
    double start = clock();
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
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
                geno_dim*g_snps+geno_dim*(g_left_marker+
                extended_snp_mapping[l])+2*parent+dose]; 
                cerr<<"Hap "<<j<<" Person "<<i<<" parent "<<parent<<" SNP "<<g_left_marker+extended_snp_mapping[l]<<" dose "<<dose<<" value "<<
                g_snp_penetrance[i*
                geno_dim*g_snps+geno_dim*(g_left_marker+
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
                  logpenetrance+=g_snp_penetrance[i*(geno_dim*g_snps)+geno_dim*(g_left_marker+extended_snp_mapping[l]) +  n]; 
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
      if (i==debug_subject) cerr<<"MAX LOG IS "<<maxlog<<endl;
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
              val[parent]>=gf_logpen_threshold?exp(val[parent]):0;
              //cerr<<"hap:"<<j<<",person "<<i<<",parent:"<<parent<<",logval:"<<logpenetrance_cache[i*2*g_max_haplotypes+2*j+parent]<<",val:"<<penetrance_cache[i*2*g_max_haplotypes+2*j+parent]<<endl;
            }
          }else{
            for(int k=0;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k] && (!haploid || j==k)){
                float val = logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]-maxlog;
                logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = val;
                penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k] = val>=gf_logpen_threshold?exp(val):0;
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

