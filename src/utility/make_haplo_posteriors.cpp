#include <iostream>
#include<sstream>
#include<math.h>
#include<cstdlib>

using namespace std;
const float gf_epsilon = .001;
const float log_gf_epsilon=log(gf_epsilon);

float confidence;

float rescale(float raw){
  return (raw-.5) * confidence + .5;
}

float convert(float raw){
  float logval;
  //cerr<<"before "<<raw;
  //cerr<<"\tafter "<<raw<<endl;
  if (raw<gf_epsilon) logval = log_gf_epsilon;
  else logval = log(raw);
  return logval;
}

int main(int argc,char * argv[]){
  if (argc<5){
    cerr<<"<persons><snps><convert to log?[0|1]><confidence (range 0-1)>\n";
    exit(1);
  }
  int arg=0;
  int n=atoi(argv[++arg]);
  int p=atoi(argv[++arg]);
  bool convertlog=atoi(argv[++arg]);
  confidence=atof(argv[++arg]);
  cerr<<"Confidence is "<<confidence<<endl;
  if (convertlog) cerr<<"INFO: converting output to log space\n";
  string line;
  while(getline(cin,line)){
    istringstream iss(line);
    float a,b,c,d; 
    iss>>a;
    for(int j=0;j<n;++j){
      iss>>a>>b>>c>>d;
      float hap1a=rescale(a+b);
      float hap1b=rescale(c+d);
      float hap2a=rescale(a+c);
      float hap2b=rescale(b+d);
      if (convertlog){
        hap1a = convert(hap1a);
        hap1b = convert(hap1b);
        hap2a = convert(hap2a);
        hap2b = convert(hap2b);
      }
      if (j) cout<<"\t";
      cout<<hap1a<<"\t"<<hap1b<<"\t"<<hap2a<<"\t"<<hap2b;
    }
    cout<<endl;
  }
}
