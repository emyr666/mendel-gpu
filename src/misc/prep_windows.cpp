#include<vector>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include"bam_interface.h"

using namespace std;

  struct read_t{
    int start;
    int max_width;
    vector<int> alleles_vec;
    vector<float> probs_vec;

    read_t(int start_pos,int max_width){
      this->start = start_pos;
      this->max_width = max_width;
    }

    void push(int allele,float prob){
      alleles_vec.push_back(allele);
      probs_vec.push_back(prob);
    }

    void make_padded_array(int * alleles_arr, float * probs_arr){
      int len = alleles_vec.size();
      int end = start+len;
      for(int i=0;i<start;++i){
        alleles_arr[i] = 9;
        probs_arr[i] = .5;
      }
      for(int i=start;i<end;++i){
        alleles_arr[i] = alleles_vec[i];
        probs_arr[i] = probs_vec[i];
      }
      for(int i=end;i<max_width;++i){
        alleles_arr[i] = 9;
        probs_arr[i] = -.5;
      }
    }
  };

class WindowUtility{
public:
  WindowUtility();
  void loop();
  int fetch_reads(int subject, int left_index, int len, int * true_snpindices, int & total_reads, int * alleles, float * probs);
  static const int MAX_READS=60;
private:
  int region_len;
  int subjects;
  int win_width;
  int win_height;
  bool * polymorphic;
  float * error_probs;
  int * hap_alleles;
  ofstream * ofs_haparr;
  int binary_rowlen; 
};

WindowUtility::WindowUtility(){
  region_len = 1000;
  subjects = 100;
  win_width = 100;
  win_height = get_max_reads();
  cerr<<"Window height is "<<win_height<<endl;
  polymorphic = new bool[win_width];
  hap_alleles = new int[subjects*win_height*win_width];
  error_probs = new float[subjects*win_height*win_width];
  ofs_haparr = new ofstream[subjects];
  binary_rowlen = sizeof(int)+sizeof(int)*win_height+sizeof(float)*win_height;
  
}

void WindowUtility::loop(){
  int left_marker = 0;
  int right_marker = win_width>region_len?region_len:win_width;
  for(int subject=0;subject<subjects;++subject){
    ostringstream oss1,oss2;
    oss1<<"bams/polymorphic."<<subject<<".bin";
    ofs_haparr[subject].open(oss1.str().data(),ios::out|ios::binary);
  }
  while(left_marker<right_marker){
    // do stuff here
    int len = right_marker-left_marker;
    cerr<<"Left and right markers: "<<left_marker<<","<<right_marker<<endl;
    for(int j=0;j<win_width;++j) polymorphic[j] = false;
    // load reads across all subjects
    for(int subject = 0;subject<subjects;++subject){
      //cerr<<"Subject "<<subject<<" at marker "<<left_marker<<endl;
      set_current_window(subject,"1",left_marker,len);
      for(int row=0;row<win_height;++row){
        //cerr<<"Read "<<row<<" for left marker "<<marker<<endl;
        get_read_data(subject,"1",row,hap_alleles+subject*win_height*win_width+row*win_width,error_probs+subject*win_height*win_width+row*win_width);
        for(int j=0;j<win_width;++j){
          // scan for polymorphisms
          //cerr<<hap_alleles[subject*win_height*win_width+row*win_width+j];
          if (hap_alleles[subject*win_height*win_width+row*win_width+j] == 1)
            polymorphic[j] = true;
        }
        //cerr<<endl;
      }
    }
    for(int j=0;j<win_width;++j){
      if (polymorphic[j]){
        int true_pos = left_marker+j;
        cerr<<"Site "<<true_pos<<" is polymorphic"<<endl;
        for(int subject = 0;subject<subjects;++subject){
          char * cvec1 = reinterpret_cast<char *>(&true_pos);
          int temp_allele[win_height];
          float temp_error[win_height];
          for(int r =0;r<win_height;++r){
            temp_allele[r] = hap_alleles[subject*win_height*win_width+r*win_width+j];
            temp_error[r] = error_probs[subject*win_height*win_width+r*win_width+j];
            //if (subject==3) cerr<<"true pos "<<true_pos<<" data: "<<temp_allele[r]<<","<<temp_error[r]<<endl;
          }
          char * cvec2 = reinterpret_cast<char *>(temp_allele);
          char * cvec3 = reinterpret_cast<char *>(temp_error);
          ofs_haparr[subject].write(cvec1,sizeof(int));
          ofs_haparr[subject].write(cvec2,sizeof(int)*win_height);
          ofs_haparr[subject].write(cvec3,sizeof(float)*win_height);
          //ofs_haparr[subject]<<endl;
        }
      }
    }
    
    // end
    left_marker = right_marker;
    right_marker = (left_marker+win_width)>region_len?region_len:(left_marker+win_width);
  }
  for(int subject=0;subject<subjects;++subject){
    ofs_haparr[subject].close();
  }
}

int WindowUtility::fetch_reads(int subject, int left_index, int window_width,int * true_snpindices, int & total_reads, int * alleles, float * probs){
  ostringstream oss1;
  oss1<<"bams/polymorphic."<<subject<<".bin";
  ifstream ifs(oss1.str().data(),ios::in|ios::binary);
  unsigned long seekpos = left_index * binary_rowlen;
  ifs.seekg(seekpos);
  int temp_allelewindow[window_width*win_height];
  float temp_probwindow[window_width*win_height];
  for(int snpindex = 0;snpindex<window_width;++snpindex){
    char * cvec1 = reinterpret_cast<char *>(true_snpindices+snpindex);
    char * cvec2 = reinterpret_cast<char *>(temp_allelewindow+snpindex*win_height);
    char * cvec3 = reinterpret_cast<char *>(temp_probwindow+snpindex*win_height);
    ifs.read(cvec1,sizeof(int));
    ifs.read(cvec2,sizeof(int)*win_height);
    ifs.read(cvec3,sizeof(float)*win_height);
    //cerr<<"Index "<<snpindex<<", true SNP pos is "<<true_snpindices[snpindex]<<endl;
    //for(int j=0;j<win_height;++j){
      //cerr<<" Alleles,errors: "<<temp_allelewindow[snpindex*win_height+j]<<","<<temp_probwindow[snpindex*win_height+j]<<endl;
    //}
  }
  ifs.close();
  vector<read_t> valid_reads;
  for(int row=0;row<win_height;++row){
    int prev_allele;
    int snpindex = 0;
    //cerr<<"At row "<<row<<endl;
    while(snpindex<window_width){
      int curr_allele = temp_allelewindow[snpindex*win_height+row];
      while(snpindex<window_width && 
      temp_allelewindow[snpindex*win_height+row]==9){
        // skip! 
        ++snpindex;
      }
      read_t current_read(snpindex,window_width);
      while(snpindex<window_width && 
      temp_allelewindow[snpindex*win_height+row]!=9){
        current_read.push(temp_allelewindow[snpindex*win_height+row],
        temp_probwindow[snpindex*win_height+row]);
        ++snpindex;
      }
      if (current_read.alleles_vec.size()>1){
        //cerr<<"Found a read of length "<<current_read.alleles_vec.size()<<" with offset "<<current_read.start<<"!"<<endl;
        valid_reads.push_back(current_read);
      }
    }
  }
  total_reads = valid_reads.size();
  vector<read_t>::iterator it;
  int read_index = 0;
  for(it=valid_reads.begin();it!=valid_reads.end();it++){
    read_t curr_read = *it;
    if (read_index<WindowUtility::MAX_READS){
      curr_read.make_padded_array(alleles+read_index*window_width,probs+read_index*window_width);
    }else{
      cerr<<"WARNING! Read index "<<read_index<<" exceed MAX_READS in window of "<<WindowUtility::MAX_READS<<endl;
    }
    ++read_index;
  }
}

int main(int argc, char * argv[]){
  WindowUtility * prepper = new WindowUtility();
  //prepper->loop();
  // testing purposes
  int subject = 3;
  int leftindex = 0;
  int len = 3;
  int true_indices[len];
  int total_reads;
  int alleles[WindowUtility::MAX_READS * len];
  float probs[WindowUtility::MAX_READS * len];
  int errcode =  prepper->fetch_reads(subject, leftindex, len, true_indices, total_reads, alleles,  probs);
  cerr<<"True SNP positions:";
  for(int i = 0;i<len;++i){
    cerr<<" "<<true_indices[i];
  }
  cerr<<endl;
  for(int r=0;r<total_reads;++r){
    cerr<<"Read "<<r<<" alleles\n";
    for(int i = 0;i<len;++i){
      cerr<<alleles[r*len+i];
    }
    cerr<<endl;
    cerr<<"Read "<<r<<" probs\n";
    for(int i = 0;i<len;++i){
      if (i) cerr<<",";
      cerr<<probs[r*len+i];
    }
    cerr<<endl;
  }
  return 0;
}
