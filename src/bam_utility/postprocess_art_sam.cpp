#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>

using namespace std;

int subject;

int main(int argc,char * argv[]){
  if (argc<2){
    cerr<<"Usage: <subject>\n";
    exit(1);
  }
  int arg = 0;
  subject = atoi(argv[++arg]);
  string line;
  while(getline(cin,line)){
    if (line[0]=='m' || line[0]=='p'){
      istringstream iss(line);
      string token;
      int j=0;
      while(iss>>token){
        if(j) cout<<"\t";
        if(j==2) cout<<subject;
        else if(j==5) cout<<"100M";
        else cout<<token;
        ++j;
      }
      cout<<endl;
    }
  }
  return 0;
}
