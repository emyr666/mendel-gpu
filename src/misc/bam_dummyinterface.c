using namespace std;

#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include "bam_interface.h"
#include<vector>
#include<map>

namespace gary_library{
  matrix_struct * matrix = NULL;
}
using namespace gary_library;

int max_reads = 10;
int max_width = 10000;
int window_width;
int leftpos;
int subject;
typedef map<string,vector<string> > mapped_matrices_t;
mapped_matrices_t dataset;

inline float phred2prob(char phred){
  return .5;
}

void allocate_datastructures(){
  // this is the time to read the test datasets in
  cerr<<"BAM API loading all files\n";
  for(int i=0;i<100;++i){
    ostringstream oss;
    oss<<"bams/"<<i;
    ifstream ifs(oss.str().data());
    if (!ifs.is_open()){
      cerr<<"Couldn't find "<<oss.str()<<endl;
      exit(1);
    }
    string line;
    vector<string> curmat;
    for(int r=0;r<max_reads;++r){
      getline(ifs,line);
      curmat.push_back(line);
    }
    dataset[oss.str()] = curmat;
    ifs.close();
  }
  cerr<<"BAM API loaded all files\n";
  matrix = new matrix_struct; 
  matrix->geno = new char[max_reads*max_width];
  matrix->qual = new char[max_reads*max_width];
}

int get_read_data(char * bam_file, uint tid, uint start, uint end, matrix_struct * read_matrices){
   if (matrix == NULL) allocate_datastructures();   
   matrix->tid = tid;
   matrix->start = start;
   matrix->end = end;
   // fill in existing implementation of get_read_data
   mapped_matrices_t::iterator it = dataset.find(bam_file);
   if (it==dataset.end()){
     cerr<<"Cannot find the dataset for the file "<<bam_file<<"!\n";
     return 1;
   }
   vector<string> curmat = it->second;
   int width = end-start+1;
   int height = 10;
   for (int i=0;i<height;++i){
     for (int j=0;j<width;++j){
       matrix->geno[i*width+j] = curmat[i][start+j];
       matrix->qual[i*width+j] = 'A';
     }
   }
   read_matrices = matrix;
   return 0;
}

// deprecated

float unif(){
  return rand()/(RAND_MAX+1.);
}

// the following function returns an error code, 0 for success.
// note that read_length,bit_array,error_prob are return parameters
//void get_read_data(int sample, const char * chrom, int read_index, int * bit_array, float * error_prob){
//  string curvec = dataset[sample][read_index];
//  string readstring = curvec.substr(leftpos,window_width);
//  for(int i=0;i<window_width;++i){
//    bit_array[i] = (int)readstring[i]-(int)'0';
//    error_prob[i] = .01;
//  }
//}

// the following function inits the max number of reads 

//int get_max_reads(){
//  return max_reads;
//}

// inform the API of where the windows begin and end (zero indexed)

//void set_current_window(int sample, const char * chrom, uint first_pos, int length)
//{
//  subject = sample;
//  leftpos = first_pos;
//  window_width = length;
//  //cerr<<"Subject is now "<<subject<<endl;
//}
