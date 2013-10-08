#include<iostream>
#include<CL/cl.h>
using namespace std;

const char * clError (cl_int rc);
void clSafe (cl_int rc, string functionname);

inline void checkErr(cl_int err, const char * name)
{
    if (err != CL_SUCCESS) {
        std::cerr << "ERROR: " << name
                 << " (" << err << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
}

template<class T> void MendelGPU::createBuffer(int rw, int dim, const char * label,cl::Buffer * & buf){
  ostringstream oss;
  oss<<"Creating buffer "<<label;
  string mesg = oss.str();
  cerr<<mesg<<" of dimension "<<dim<<endl;
  buf = new cl::Buffer(*context, rw, sizeof(T) * dim, NULL, &err);
  clSafe(err,oss.str().data());
}

template<typename T> void MendelGPU::setArg(cl::Kernel * &  kernel,int & index,T arg,const char * label){
  ostringstream oss;
  oss<<"Setting argument "<<index<<" for kernel "<<label;
  string mesg = oss.str();
  cerr<<mesg<<endl;
  err = kernel->setArg(index++, arg);
  clSafe(err,oss.str().data());
}

template<typename T> void MendelGPU::writeToBuffer(cl::Buffer * & buffer,int dim,T hostArr,const char * label){
  ostringstream oss;
  T t;
  oss<<"Writing to buffer "<<label;
  string mesg = oss.str();
  cerr<<mesg<<" of dimension "<<dim<<endl;
  err = commandQueue->enqueueWriteBuffer(*buffer, CL_TRUE, 0, sizeof(*t)*dim,hostArr,NULL,NULL);
  clSafe(err,oss.str().data());
}

template<typename T> void MendelGPU::readFromBuffer(cl::Buffer * & buffer,int dim,T hostArr,const char * label){
  ostringstream oss;
  T t;
  oss<<"Reading from buffer "<<label;
  string mesg = oss.str();
  cerr<<mesg<<" of dimension "<<dim<<endl;
  err = commandQueue->enqueueReadBuffer(*buffer, CL_TRUE, 0, sizeof(*t)*dim,hostArr,NULL,NULL);
  clSafe(err,oss.str().data());
}

