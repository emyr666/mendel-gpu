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

set<int> position_set;

void parse_snplist(const char * filename){
  ifstream ifs(filename);
  string line; 
  getline(ifs,line); // strip header
  snps = 0;
  while(getline(ifs,line)){
    istringstream iss(line);
    string chr,id,ref,alt;
    int pos;
    float freq;
    iss>>chr>>pos>>id>>ref>>alt>>freq;
    position_set.insert(pos);
  }
  ifs.close();
}

void filter_vcf(){
  string line;
  while(getline(cin,line)){
    if (line[0]=='1'){
      istringstream iss(line);
      string token;
      int pos;
      iss>>token;
      iss>>pos;
      if (position_set.find(pos)!=position_set.end()){
        cout<<line<<endl;
      }
    }else{
      cout<<line<<endl;
    }
  }
}

int main(int argc, char *argv[])
{
  try{
    if (argc < 2) {
      fprintf(stderr, "Usage: <snplistfile>\n");
      return 1;
    }
    int arg = 0;
    const char * snp_filename = argv[++arg];
    parse_snplist(snp_filename);
    filter_vcf();
    return 0;
  }catch (const char * & mesg){
    cerr<<"Exception caught of message "<<mesg<<endl;
  }
}
