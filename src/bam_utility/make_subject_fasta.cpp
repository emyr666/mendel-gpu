#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>

using namespace std;

const int MAX_SEQUENCE=300000000;
char * ref_sequence1;
char * ref_sequence2;
int sim_len;

void parse_variants(){
  string line;
  getline(cin,line);
  while(getline(cin,line)){
    //cerr<<"line: "<<line<<endl;
    istringstream iss(line);
    int chr,pos;
    string id,variants,ref,alt;
    float freq;
    iss>>chr>>pos>>id>>ref>>alt>>freq>>variants;
    // the variant file is 1-indexed
    //cerr<<"pos,id,ref,alt,freq: "<<pos<<","<<id<<","<<ref<<","<<alt<<","<<freq<<" variants: "<<variants<<endl;
    ref_sequence1[pos-1] = variants[0];
    //cerr<<"Assigning maternal position "<<pos<<" with "<<haplo1[pos-1]<<endl;
    ref_sequence2[pos-1] = variants[1];
    //cerr<<"Assigning paternal position "<<pos<<" with "<<haplo2[pos-1]<<endl;
  }
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

void load_ref_sequence(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Cannot find fasta file "<<filename<<endl;
    exit(1);
  }
  string line;
  getline(ifs,line);
  int cursor = 0;
  while(getline(ifs,line)){
    int strlen = line.length();
    for(int i=0;i<strlen;++i){
      ref_sequence1[cursor+i] = line[i];
      ref_sequence2[cursor+i] = line[i];
    }
    cursor+=strlen;
  }
  ifs.close();

}

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


int main(int argc,char * argv[]){
  if (argc<3){
    cout<<"Usage: <simulated_len> <refseq fasta>\n";
    return 1;
  }
  srand(1);
  int arg=0;
  sim_len = static_cast<int>(atof(argv[++arg]));
  const char *  fasta_file = argv[++arg];
  ref_sequence1 = new char[MAX_SEQUENCE];
  ref_sequence2 = new char[MAX_SEQUENCE];
  cerr<<"Loading ref seq\n";
  load_ref_sequence(fasta_file);
  cerr<<"parsing variants\n";
  parse_variants();
  print_fasta<char *>("maternal",ref_sequence1,sim_len);
  print_fasta<char *>("paternal",ref_sequence2,sim_len);
  return 0;
}
