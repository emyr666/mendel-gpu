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
#include"math.h"

using namespace std;
int snps;
int samples;

map<int,string> position_map;

void parse_vcf(){
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
      if (position_map.find(pos)==position_map.end()){
        string id,ref,alt;
        iss>>id>>ref>>alt;
        string alleles = "NN";
        alleles[0] = ref[0];
        alleles[1] = alt[0];
        position_map[pos] = alleles;
      }
    }
    ifs.close();
  }
}

void print_matrix(){
  int j = 0;
  cout<<"CHROM\tPOS\tID\tREF\tALT\tFREQ\n";
  for(map<int,string>::iterator it = position_map.begin(); it!=position_map.end(); it++){
     int pos = it->first;
     string alleles = it->second;
     cout<<"1\t"<<pos<<"\tSNP"<<j<<"\t"<<alleles[0]<<"\t"<<alleles[1]<<"\t0\n";
     ++j;
  }
}

int main(int argc, char *argv[])
{
  try{
    if (argc < 2) {
      fprintf(stderr, "Usage: <samples> \n");
      return 1;
    }
    int arg = 0;
    samples = atoi(argv[++arg]);
    parse_vcf();
    print_matrix();
    return 0;
  }catch (const char * & mesg){
    cerr<<"Exception caught of message "<<mesg<<endl;
  }
}
