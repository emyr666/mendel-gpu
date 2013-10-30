#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>

using namespace std;

char * haplo1, * haplo2;

template <class T> void print_fasta(string label,T haplo,int sim_len){
  cout<<">"<<label<<"\n";
  int i = 0;
  for(;i<sim_len;++i){
    cout<<haplo[i];
    if ((i+1)%60==0) cout<<endl;
  }
//  cerr<<"i ended at "<<i<<endl;
  if (i%60!=0) cout<<endl;
}



template<class T> void init_haplo(int sim_len, T haplo){
  cerr<<"Generating chromosome bases\n";
  for(int i=0;i<sim_len;++i) {
    float r = rand()/(RAND_MAX+1.);
    char base;
    if (r<.25){
      base = 'A';
    }else if (r<.5){
      base = 'C';
    }else if (r<.75){
      base = 'G';
    }else{
      base = 'T';
    }
    haplo[i] = base;
  }
  cerr<<"Done!\n";
}

int main(int argc,char * argv[]){
  if (argc<2){
    cout<<"Usage: <simulated_len> \n";
    return 1;
  }
  srand(1);
  int arg=0;
  int sim_len = static_cast<int>(atof(argv[++arg]));
  haplo1 = new char[sim_len];
  init_haplo(sim_len,haplo1);
  print_fasta<char *>("1",haplo1,sim_len);
  return 0;
}
