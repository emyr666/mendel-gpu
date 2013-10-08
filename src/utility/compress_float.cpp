#include<iostream> 
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<math.h>
using namespace std;

int dim;
int n,m;
float * inputmatrix;
const float gf_epsilon = .001;
const float log_gf_epsilon=log(gf_epsilon);

float convert(float raw){
  float logval;
  if (raw<gf_epsilon) logval = log_gf_epsilon;
  else logval = log(raw);
  return logval;
}


void stream_out(float * vec, long unsigned int len){
  char * charvec = reinterpret_cast<char *>(vec);
  long unsigned int outlen = sizeof(float) * len;
  cout.write(charvec,outlen);
  //float * newvec = new float[len];
  //newvec = reinterpret_cast<float * >(charvec);
  //for(int i=0;i<len;++i) cerr<<" "<<newvec[i];
  //cerr<<endl;
}

int main(int argc,char * argv[]){
  if (argc<2){
    cerr<<"Usage: <record element len> \n";
    exit(1);
  }
  int arg=0;
  dim = atoi(argv[++arg]);
  cerr<<"NOTICE: geno dimension is "<<dim<<endl;
  //bool convertlog = atoi(argv[++arg]);
  //if (convertlog) cerr<<"NOTICE: Converting to log space\n";
  string line;
  getline(cin,line);
  istringstream iss(line);
  iss>>n>>m;
  inputmatrix = new float[n*dim*m];
  for(int j=0;j<m;++j){
    getline(cin,line);
    istringstream iss(line);
    for(int i=0;i<n;++i){
      for(int k=0;k<dim;++k){
        float val;
        iss>>val;
        inputmatrix[i*dim*m+dim*j+k]=val;
        //inputmatrix[i*dim*m+dim*j+k]=convertlog?convert(val):val;
 //cerr<<"debug:i,j,k"<<i<<","<<j<<","<<k<<":"<<
 //       inputmatrix[i*dim*m+dim*j+k]<<endl;
      }
    }
  }
  for(int i=0;i<n;++i){
    float * subjectpen = inputmatrix+i*dim*m;
    stream_out(subjectpen,m*dim);
  }
  delete[]inputmatrix;
  return 0;
}
