#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<cstdlib>

using namespace std;
string chrom;
string basepath;
string permfile;
string sexfile;
unsigned int n;
ofstream * ofs_arr;
int chunksize;
int chunks;
int * perm_arr;
char * sexarr;

bool fexists(const char *filename)
{
  bool exists = false;
  ifstream ifile(filename);
  if (ifile.good()){
    ifile.close();
    exists = true;
  }
  return exists;
}

void load_perm(){
  perm_arr = new int[n];
  ostringstream oss;
  oss<<basepath<<"/"<<permfile;
  ifstream ifs(oss.str().data());
  if (!ifs.is_open()) {
    cerr<< "Warning: Cannot open permutation file\n";
    for(int i=0;i<n;++i){
      perm_arr[i] = i;
    }
  }else{
    string line;
    getline(ifs,line); // skip header
    int i =0;
    while(getline(ifs,line)){
      istringstream iss(line);
      int seq,perm;
      iss>>seq>>perm;
      perm_arr[seq] = perm;
      ++i;
    }
    ifs.close();
  }
}

void load_sex(){
  ofstream * ofs_arr = new ofstream[chunks];      
  bool write_chunk = false;
  for(int i=0;i<chunks;++i){
    ostringstream oss;
    oss<<basepath<<"/"<<sexfile<<".chunk."<<i;
    if (!fexists(oss.str().data())){
      write_chunk = true;
      cerr<<"Need to write sex chunk\n";

      ofs_arr[i].open(oss.str().data());
      ofs_arr[i]<<"seq\tsex\n";
    }else{
      cerr<<"No Need to write sex chunk\n";
    }
  }

  sexarr = new char[n];
  ostringstream oss;
  oss<<basepath<<"/"<<sexfile;
  ifstream ifs(oss.str().data());
  if (!ifs.is_open()) {
    cerr<< "Warning: Cannot open sex file"<<sexfile<<"\n";
    for(int i=0;i<n;++i){
      sexarr[i] = 'F';
    }
  }else{
    string line;
    getline(ifs,line); // skip header
    int i =0;
    while(getline(ifs,line)){
      istringstream iss(line);
      int seq;
      char sex;
      iss>>seq>>sex;
      sexarr[i] = sex;
      if (write_chunk) ofs_arr[i/chunksize]<<i<<"\t"<<sexarr[i]<<endl;
      ++i;
    }
    ifs.close();
  }
  if (write_chunk){
    for(int i=0;i<chunks;++i){
      ofs_arr[i].close();
    }
    delete[] ofs_arr;
  }
}


int main(int argc,char * argv[]){
  if (argc<7){
    cerr<<"Usage: <chr> <basepath> <chunksize> <permfile> <chr> <sexfile>\n";
    exit(1);
  }
  int arg=0;
  chrom = argv[++arg];
  basepath = argv[++arg];
  chunksize = atoi(argv[++arg]);
  permfile = argv[++arg];
  string chr = argv[++arg];
  sexfile = argv[++arg];

  string line;
  getline(cin,line);
  getline(cin,line);
  n = atoi(line.data());
  chunks = n/chunksize+1;
  load_perm();
  load_sex();
  //exit(1);
  ofs_arr = new ofstream[chunks];      
  cerr<<"chr "<<chrom<<", chunks "<<chunks<<endl;
  getline(cin,line);
  for(int i=0;i<chunks;++i){
    ostringstream oss;
    oss<<basepath<<"/chr."<<chrom<<".chunk."<<i;
    ofs_arr[i].open(oss.str().data());
    if (!ofs_arr[i].is_open()) throw "Cannot open file\n";
  }
  while(getline(cin,line)){
    istringstream iss(line);
    string chr,pos,rs,oldvec;
    iss>>chr>>pos>>rs>>oldvec;
    char vec[n];
    for(int i=0;i<n;++i){
      vec[i] = oldvec[perm_arr[i]];
    }
    for(int i=0;i<chunks;++i){
      ofs_arr[i]<<chr<<"\t"<<pos<<"\t"<<rs<<"\t";
    }
    for(int i=0;i<n;++i){
      if (chr.compare("X")==0 && sexarr[i]=='M' && vec[i]=='2'){
        vec[i]='3';
      }
      ofs_arr[i/chunksize]<<vec[i];
    }
    //cerr<<"line: "<<line<<endl;
    for(int i=0;i<chunks;++i){
      ofs_arr[i]<<endl;
    }
  }
  for(int i=0;i<chunks;++i){
    ofs_arr[i].close();
  }
}
