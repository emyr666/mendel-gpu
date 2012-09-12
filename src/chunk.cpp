#include<iostream>
#include<sstream>
#include<fstream>

using namespace std;
string chrom;
string basepath;
unsigned int n;
ofstream * ofs_arr;
int chunksize;

int main(int argc,char * argv[]){
  if (argc<4){
    cerr<<"Usage: <chr> <basepath> <chunksize>\n";
    exit(1);
  }
  int arg=0;
  chrom = argv[++arg];
  basepath = argv[++arg];
  chunksize = atoi(argv[++arg]);
  string line;
  getline(cin,line);
  getline(cin,line);
  n = atoi(line.data());
  int chunks = n/chunksize+1;
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
    string chr,pos,rs,vec;
    iss>>chr>>pos>>rs>>vec;
    for(int i=0;i<chunks;++i){
      ofs_arr[i]<<chr<<"\t"<<pos<<"\t"<<rs<<"\t";
    }
    for(int i=0;i<n;++i){
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
