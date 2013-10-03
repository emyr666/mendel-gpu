using namespace std;

#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<vector>

int main(int argc,char * argv[]){
  if (argc<4){
    cerr<<"Usage: <basefile> <partsize> <header=1|noheader=0>\n";
    exit(0);
  }
  int arg=0;
  const char * basefile = argv[++arg];
  int max = atoi(argv[++arg]);
  bool is_header = atoi(argv[++arg]);
  string line;
  int i=0,c=0;
  string header;
  if (is_header){
    getline(cin,header);
  }
  ofstream ofs;
  string outfile;
  while(getline(cin,line)){
    if ((i % max) == 0){
      if (i){
        ofs.close();
      }
      ++c;
      ostringstream oss;
      oss<<basefile<<c;
      outfile = oss.str();
      ofs.open(outfile.data());
      cerr<<"Writing to "<<outfile<<endl;
      if (is_header){
        ofs<<header<<endl;
      }
    }
    ofs<<line<<endl;
    ++i;
  }
  ofs.close();
}
