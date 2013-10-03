using namespace std;
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<tr1/unordered_set>
#include<tr1/unordered_map>
#include<set>
#include<queue>
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

void DenovoMendelGPU::impute_genotypes(){
  impute_diploid_genotypes_denovo_((g_center_snp_start-g_left_marker),g_center_snp_start);
}

void DenovoMendelGPU::impute_haploid_genotypes_denovo_(int center_snp, int  d_snp){
  int current_snp = d_snp;
  int c_snp = center_snp;
  cerr<<"debug begin haploid impute_genotypes at SNP "<<current_snp<<", center: "<<c_snp<<endl;
  return;
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


void DenovoMendelGPU::impute_diploid_genotypes_denovo_(int  center_snp, int  d_snp){
  int c_snp = center_snp;
  int current_snp = d_snp;
  cerr<<"begin diploid impute_genotypes at SNP "<<current_snp<<", center: "<<c_snp<<endl;
  //return ;
  bool debug_dosage = false;
  bool debug_geno = false;
  bool debug_posterior = false;
  bool debug_pen = current_snp == -3;
  int max_geno=g_genotype_imputation?3:4;
  cerr<<"Max geno "<<max_geno<<endl;
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
        cout<<"GPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
        for(int j=0;j<max_geno;++j){
          cout<<"\t"<< subject_posterior[i*4+j];
        }
        cout<<endl;
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
    //if (g_likelihood_mode==Config::LIKELIHOOD_MODE_READS && 
    //  window_polymorphisms[c_snp]){
    float dosages[g_people];
    int genotypes[g_people*max_geno];
    float posteriors[g_people*max_geno];
//    if (outfile_format.compare(FORMAT_MEC)==0){
//      ofs_dosage_file<<current_snp<<"\t";
//      ofs_genotype_file<<current_snp<<"\t";
//      ofs_posterior_file<<current_snp;
//    }
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
      float best_prob = 0;
      float posterior_prob[4];
      int best_pair[2];
      best_pair[0] = best_pair[1] = 0;
      memset(posterior_prob,0,sizeof(float)*4);
      for(int j=0;j<g_max_haplotypes;++j){
        if(g_active_haplotype[j] ){
          for(int k=j;k<g_max_haplotypes;++k){
            if(g_active_haplotype[k] && (!haploid || j==k)){
              int m;
              float penetrance = penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
              //float penetrance = 0;
              //cerr<<"Subject "<<i<<" penetrance at haps "<<j<<" and "<<k<<" :" <<penetrance<<endl;
              if (penetrance>0 ){
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
      }
      // assume base genotype for non polymorphisms in sequencing data
      //if (g_likelihood_mode==Config::LIKELIHOOD_MODE_READS && 
      //!window_polymorphisms[c_snp]){
        //cerr<<"Assuming non polymorphism for this sequenced site\n";
        //posterior_prob[0] = 1;
        //for(int j=1;j<max_geno;++j){
         // posterior_prob[j] = 0;
        //}
      //}else{
        //cerr<<"Assuming polymorphism for this sequenced site\n";
        //cerr<<"CPU POSTERIOR "<<current_snp<<","<<i;
        //for(int j=0;j<max_geno;++j){
          //cerr<<" "<<posterior_prob[j];
        //}
       // cerr<<endl;
        
      //}
      if (debug_posterior){
        cout<<"CPU POSTERIOR "<<current_snp<<","<<i;
        for(int j=0;j<max_geno;++j){
          cout<<" "<<posterior_prob[j];
        }
        cout<<endl;
      }
      float denom = epsilon;
      for(int m=0;m<max_geno;++m){
        //if (posterior_prob[m]<epsilon) posterior_prob[m] = epsilon;
        denom+=posterior_prob[m];
      }
      if (denom==0) {
        cerr<<"divide by zero denom\n"; 
        for(int m=0;m<max_geno;++m){
          cerr<<"posterior for geno "<<m<<": "<<posterior_prob[m]<<endl;
        }
        cerr<<"hap1,hap2:pen,logpen,freq,genotype\n";
        for(int j=0;j<g_max_haplotypes;++j){
          int start = g_genotype_imputation?j:0;
          for(int k=start;k<g_max_haplotypes;++k){
            if(g_active_haplotype[j] && g_active_haplotype[k]){
              int m = g_haplotype[j*g_max_window+ c_snp] + 
              g_haplotype[k*g_max_window+ c_snp];
              float logpenetrance = logpenetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
              float penetrance = penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k];
              float freq = g_frequency[j]*g_frequency[k];
              cerr<<j<<","<<k<<":"<<penetrance<<","<<logpenetrance<<","<<freq<<","<<m<<endl;
            }
          }
        }
        exit(1);
      }
//      float dose = 0;
//      if (outfile_format.compare(FORMAT_MEC)==0){
//        for(int j=0;j<max_geno;++j){
//          float p = posterior_prob[j]/denom;
//          ofs_posterior_file<<"\t"<<p;
//          dose+=j*p;
//        }
//        char chardose =  (int)(dose*10)+'A';
//        ofs_dosage_file<<chardose;
//      }else if (outfile_format.compare(FORMAT_DEFAULT)==0){
//        ofs_posterior_file<<"CPU_POSTERIOR\t"<<i<<"\t"<<current_snp;
//        for(int j=0;j<max_geno;++j){
//          float p = posterior_prob[j]/denom;
//          ofs_posterior_file<<"\t"<< p;
//          dose+=j*p;
//        }
//        ofs_dosage_file<<"CPU_DOSE:\t"<<i<<"\t"<<current_snp<<"\t"<<dose<<endl;
//        ofs_posterior_file<<endl;
//      }
      dosages[i] = 0;
      for(int j=0;j<max_geno;++j){
        float p = posterior_prob[j]/denom;
        posteriors[i*max_geno+j] = p;
        dosages[i]+=j*p;
      }
      //subject_dosages[i] = dose;
      float maxval = 0;
      for(int j=0;j<4;++j){
        if (posterior_prob[j]>maxval) maxval = posterior_prob[j];
      }
      if (maxval>0){
        //int j = best_pair[0];
        //int k = best_pair[1];
        //cerr<<"best prob "<<best_prob<<" "<<c_snp<<" "<<j<<" "<<k<<" "<<g_haplotype[j*g_max_window+ c_snp] <<" "<<g_haplotype[k*g_max_window+ c_snp]<<endl;
        int geno1= -1;
        //int geno1 = g_genotype_imputation?
        //g_haplotype[j*g_max_window+ c_snp] +
        //g_haplotype[k*g_max_window+ c_snp]:
        //2*g_haplotype[j*g_max_window+ c_snp] +
        //g_haplotype[k*g_max_window+ c_snp];
        float max_prob = posterior_prob[0];
        int geno2 = 0;
        for(int h=1;h<4;++h){
          if (posterior_prob[h]-max_prob>epsilon*max_prob){
            max_prob = posterior_prob[h];
            geno2 = h;
          }
        }  
        if (debug_geno){
          cout<<"CPU GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<geno1<<","<<geno2<<endl;
        }
        genotypes[i] = geno2;
        //if (outfile_format.compare(FORMAT_MEC)==0){
        //  ofs_genotype_file<<geno2;
        //}else if (outfile_format.compare(FORMAT_DEFAULT)==0){
         // ofs_genotype_file<<"CPU_GENO\t"<<i<<"\t"<<current_snp<<"\t"<<geno2<<endl;;
        //}
        
        //ofs_genotype_file<<"CPU_GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<geno1<<","<<geno2<<"\t"<<best_pair[0]<<","<<best_pair[1]<<endl;
      }
    }
    if (max_geno==3){
      io_manager->writeDosage(current_snp,dosages,g_people);
      io_manager->writeGenotype(current_snp,genotypes,g_people);
    }
    io_manager->writePosterior(max_geno,current_snp,posteriors,g_people);
    //if (outfile_format.compare(FORMAT_MEC)==0){
    //  ofs_dosage_file<<endl;
    //  ofs_genotype_file<<endl;
    //  ofs_posterior_file<<endl; 
    //}
    float rsq = compute_rsq(dosages,1,0);
    //ofs_quality_file<<current_snp<<"\t"<<rsq<<endl;
    io_manager->writeQuality(current_snp,rsq);
  }//END CPU VERSION
  cerr<<"done impute_geno\n";
  if (debug_dosage || debug_posterior|| debug_geno) exit(1);
  if (debug_pen) {
    if (debug_mode) ofs_debug_haplotype_file.close();
    exit(1);
  }
}
