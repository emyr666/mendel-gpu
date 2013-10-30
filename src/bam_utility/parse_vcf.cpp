#include<assert.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<map>
#include<set>
#include<stdlib.h>
#include<cstdio>
#include<stdio.h>
#include<string.h>
#include"faidx.h"
#include "sam.h"
#include"math.h"

using namespace std;
int snps;
int samples;

map<int,int> position_map;
vector<int> snp_vec;
float * prob_matrix;

void parse_snplist(const char * filename){
  ifstream ifs(filename);
  string line; 
  getline(ifs,line);
  snps = 0;
  while(getline(ifs,line)){
    istringstream iss(line);
    string chr,id,ref,alt;
    int pos;
    float freq;
    iss>>chr>>pos>>id>>ref>>alt>>freq;
    position_map[pos]=snps;
    snp_vec.push_back(pos);
    ++snps;
  }
  ifs.close();
}

void parse_vcf(bool convert_phred){
  for(int i=0;i<samples;++i){
    ostringstream oss;
    oss<<"subject"<<i<<".vcf";
    string filename=oss.str();
    cerr<<"Reading from "<<filename<<endl;
    ifstream ifs(filename.data());
    string line;
    for(int j=0;j<35;++j) getline(ifs,line);
    while(getline(ifs,line)){
      istringstream iss(line);
      string token;
      int pos;
      iss>>token;
      iss>>pos;
      if (position_map.find(pos)!=position_map.end()){
        int row = position_map[pos];
        for(int k=0;k<8;++k) iss>>token;
        //cout<<"pos "<<pos<<" row "<<row<<" phred: "<<token<<"\n";
        char phredstr[token.size()+1];
        strcpy(phredstr,token.data());
        char * pch = strtok(phredstr,",");
        int k = 0;
        float denom = 0;
        while(pch!=NULL){
          if (k<3){
            int iphred = atoi(pch);
            if (convert_phred){
              float prob = pow(10,-iphred/10);
              //cout<<" "<<pch<<","<<prob<<endl;
              prob_matrix[row*samples*3+i*3+k] = prob;
              denom+=prob;
            }else{
              prob_matrix[row*samples*3+i*3+k] = iphred;
            }
          }
          pch = strtok(NULL,",");
          ++k;
        }
        if (convert_phred){
          for(int k=0;k<3;++k){
            prob_matrix[row*samples*3+i*3+k] /= denom;
          }
        }
      }
    }
    ifs.close();
  }
}

void print_matrix(){
  for(int j=0;j<snps;++j){
    cout<<snp_vec[j];
    for(int i=0;i<samples;++i){
      cout<<"\t";
      for(int k=0;k<3;++k){
        if (k) cout<<" ";
        cout<<prob_matrix[j*samples*3+i*3+k];
      }
    }
    cout<<endl;
  }
}

int main(int argc, char *argv[])
{
  try{
    if (argc < 4) {
      fprintf(stderr, "Usage: <samples> <snplistfile> <convert phred?>\n");
      return 1;
    }
    int arg = 0;
    samples = atoi(argv[++arg]);
    const char * snp_filename = argv[++arg];
    bool convert_phred = atoi(argv[++arg]);
    parse_snplist(snp_filename);
    prob_matrix = new float[snps*3*samples];
    for(int j=0;j<snps;++j){
      for(int i=0;i<samples;++i){
        float init[3];
        if (convert_phred){
         init[0] = 1;
         init[1] = 0;
         init[2] = 0;
        }else{
         init[0] = 0;
         init[1] = 100;
         init[2] = 100;
        }
        for(int k=0;k<3;++k){
          prob_matrix[j*samples*3+i*3+k] = init[k];
        }
      }
    }
    cerr<<"Initialized\n";
    parse_vcf(convert_phred);
    print_matrix();
    return 0;
  }catch (const char * & mesg){
    cerr<<"Exception caught of message "<<mesg<<endl;
  }
}
