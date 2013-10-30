#include<cstdlib>
#include<iostream>
using namespace std;
 
int inrows,incols;
char * inmatrix;

int main(int argc,char * argv[]){
  int arg = 0;
  if (argc<3){
    cerr<<"<inputrows><inputcols>\n";
    return 1;
  }
  inrows = atoi(argv[++arg]);  
  incols = atoi(argv[++arg]);  
  inmatrix = new char[inrows*incols];
  string line;
  int j =0;
  while(getline(cin,line)){
    for(uint i=0;i<line.length();++i){
      inmatrix[j*incols+i] = line[i];
    }
    ++j;
  }
  for(int i=0;i<incols;++i){
    for(int j=0;j<inrows;++j){
      cout<<inmatrix[j*incols+i];
    }
    cout<<endl;
  }
  
  return 0;
}
