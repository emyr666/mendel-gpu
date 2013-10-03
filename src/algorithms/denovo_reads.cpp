#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_penetrance.hpp"

struct read_t{
  int len;
  int * alleles;
  float * errors;
  string allelestr;

  read_t(int len,int * alleles, float * errors);
  read_t(const read_t & dup);
  ~read_t();
  void str_rep();
  void allocate(int len);
};

struct byUniqueReads{
  bool operator()(const read_t & p1,const read_t & p2) const{
    // sort in ascenting order
    return(p1.allelestr<p2.allelestr);
  }
};

DenovoReadsMendelGPU::DenovoReadsMendelGPU(IO_manager * io):DenovoMendelGPU(io){
}

DenovoReadsMendelGPU::~DenovoReadsMendelGPU(){
  cerr<<"Entering destructor denovo reads haplotyper\n";
  delete[] parsers;
  cerr<<"Parsers deleted\n";
  cerr<<"Exiting destructor denovo reads haplotyper\n";

}

void DenovoReadsMendelGPU::allocate_memory(){
  DenovoMendelGPU::allocate_memory();
  cerr<<"Initializing variables for denovo reads haplotyper\n";
  cerr<<"Analysis on "<<g_snps<<" SNPs and "<<g_people<<" persons\n";
  log_half = log(.5);
  string variants = config->bimfile;
  ReadParser::parse_positions(variants.data());
  parsers = new ReadParser[g_people];
  for(int i=0;i<g_people;++i){
    string bam = config->bamfilevec[i];
    //cerr<<"bam file "<<bam<<endl;
    parsers[i].init(bam.data());
    //parsers[i].extract_region("1",0,1);
  }
  cerr<<"Initialized variables for denovo reads haplotyper\n";
}

void DenovoReadsMendelGPU::init_window(){
  DenovoMendelGPU::init_window();
  // Load the reads for this current window
  return;
  if (g_likelihood_mode==Config::LIKELIHOOD_MODE_READS){
    for(int i=0;i<g_max_window;++i) window_polymorphisms[i] = false;
    // generate the window_informative_snp_indices data structure;
    double start = clock();
    float min_pen = 1;
    for(int i=0;i<g_people;++i){
      int temp_alleles[g_max_window];
      float temp_errors[g_max_window];
      set<read_t,byUniqueReads> read_set;
      map<string,int> read_counts;
      set<read_t,byUniqueReads>::iterator read_it;
      g_unique_reads[i] = 0;
      //set_current_window(i,chromosome,left_marker,g_markers);
      for(int r=0;r<MAX_READS;++r){
        // these are haplotype alelles
        // fix this for the new API
        //get_read_data(i, chromosome, r, temp_alleles,temp_errors);
        read_t read(g_markers,temp_alleles,temp_errors);
        read_it = read_set.find(read); 
        if (read_it==read_set.end()){
          read_counts[read.allelestr] = 1;
          for(int j=0;j<g_max_window;++j){
            g_unique_reads_alleles[i*MAX_READS*g_max_window+g_unique_reads[i]*g_max_window+j] = temp_alleles[j];
            g_unique_reads_errors[i*MAX_READS*g_max_window+g_unique_reads[i]*g_max_window+j] = temp_errors[j];
            if (temp_alleles[j]==1) window_polymorphisms[j] = true;
          }
          read_set.insert(read);
          ++g_unique_reads[i];
        }else{
          read_t read_found = *read_it;
          ++read_counts[read_found.allelestr];
        }
      }
      int read_index=0;
      for(read_it = read_set.begin();read_it!=read_set.end();read_it++){
        read_t read = *read_it;
        g_unique_reads_counts[i*MAX_READS+read_index] = read_counts[read.allelestr]; 
        ++read_index;
      }
    } // end loop over individuals
    cerr<<"Polymorphisms:\n";
    polymorphic_window = false; 
    for(int j=0;j<g_markers;++j){
      if (window_polymorphisms[j]) polymorphic_window = true;
      cerr<<window_polymorphisms[j];
    }
    cerr<<endl;
    for(int i=0;i<g_people;++i){
      // let us now cache the penetrance of each haplotype 
      // conditional on each read
      int * read_block =  g_unique_reads_alleles + i*MAX_READS*g_max_window;
      float * error_block = g_unique_reads_errors+i*MAX_READS*g_max_window;
      for(int read_index = 0;read_index<g_unique_reads[i];++read_index){
        for(int hap = 0;hap<g_max_haplotypes;++hap){
          if (g_active_haplotype[hap]){
            float logpen = 0;
            for(int marker = 0;marker<g_markers;++marker){
              if (read_block[read_index*g_max_window+marker]==0||
              read_block[read_index*g_max_window+marker]==1){
                float error = error_block[read_index*g_max_window+marker];
                float non_error = (1.-error);
                if (g_haplotype[hap*g_max_window+marker]==
                read_block[read_index*g_max_window+marker]){
                  logpen+=log(non_error);
                  //cerr<<"M";
                }else{
                  logpen+=log(error);
                  //cerr<<"U";
                }
                //pen*=g_haplotype[marker]==
                //read_block[read_index*g_max_window+marker]? non_error:error;
              }else{
                //cerr<<"N";
              }
            }
            //cerr<<endl;
            //if (pen<epsilon) pen = epsilon;
            //if (pen<min_pen) min_pen = pen;
            g_unique_reads_logpenetrance[i*MAX_READS*g_max_haplotypes+read_index*g_max_haplotypes+hap] = logpen;
            //cerr<<"Person "<<i<<" penetrance for read "<<read.allelestr<<" of hap "<<hap<<" is "<<g_unique_reads_penetrance[i*MAX_READS*g_max_haplotypes+read_index*g_max_haplotypes+hap]<<endl;
          } // if active hap
        } // loop haps
      }
      //cerr<<"Min penetrance across subjects: "<<min_pen<<endl; 
    }
    cerr<<"Total time for precomputing penetrance: "<<
    (clock()-start)/CLOCKS_PER_SEC<<endl;
  }
  cerr<<"Max_window: "<<g_max_window<<" Haplotypes: "<<g_haplotypes<<" Markers: "<<g_markers<<" left_marker: "<<g_left_marker<<endl;
  return;
  //if (g_markers>4) exit(0);
  if (run_gpu){
    init_window_opencl();
  }
}

void read_t::allocate(int len){
  alleles = new int[len];
  errors = new float[len];
} 

read_t::read_t(int len,int * alleles, float * errors){
  //len = 0;
  //alleles = NULL;
  //errors = NULL;
  this->len = len ;
  allocate(len);
  //cerr<<"Calling default constructor\n";
  for(int i=0;i<len;++i){
    this->alleles[i] = alleles[i];
    this->errors[i] = errors[i];
  }
  str_rep();
}

//void init(int * alleles, 

read_t:: read_t(const read_t & dup){
  //cerr<<"Calling copy constructor\n";
  len = dup.len;
  allocate(len);
  for(int i=0;i<len;++i){
    alleles[i] = dup.alleles[i];
    errors[i] = dup.errors[i];
  }
  str_rep();
}

void read_t::str_rep(){
  ostringstream oss;
  //cerr<<"debug_str_rep: ";
  for(int i=0;i<len;++i){
    //cerr<<alleles[i];
    if (alleles[i]==0||alleles[i]==1){
      oss<<alleles[i];
    }else{
      oss<<9;
    }
  }
  //cerr<<endl;
  allelestr = oss.str();
}

read_t::~read_t(){
  if (alleles!=NULL) { delete[] alleles;}
  if (errors!=NULL) delete[] errors;
}
