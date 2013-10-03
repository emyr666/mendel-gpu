#include<math.h>
#include<sstream>
#include<fstream>
#include<iostream>
#include<list>
#include<queue>
#include<set>
#include<tr1/unordered_map>
#include<tr1/unordered_set>
#ifdef USE_GPU
#include<CL/cl.hpp>
#endif

using namespace std;


class MendelGPU{
public:
  MendelGPU(IO_manager *io);
  void init();
  virtual ~MendelGPU();
  void run_sliding_window();
protected:
  virtual void allocate_memory();
  virtual void init_opencl();
  virtual void load_datasets()=0;
  virtual void init_window();
  virtual void compute_penetrance()=0;
  virtual void impute_genotypes()=0;
  virtual void finalize_window();
  // custom data structures
  struct hapobj_t{
    string hapstr;
    int hapcount;
    float hapfreq;
    int array_index;
    set<string> extended_set;
  };
  static const int HAPLOTYPE_MODE_DENOVO = 0;
  static const int HAPLOTYPE_MODE_GUIDE = 1;
  static const int UNPHASED_INPUT = 3;
  static const int PHASED_INPUT = 4;

  IO_manager * io_manager;
  Config  * config;

  // constants
  float epsilon;
  float convergence_criterion;
  float log_half;
  string FORMAT_DEFAULT;
  string FORMAT_MEC;
  float logpenetrance_threshold;
  bool debug_mode;
  bool debug_truth;
  bool run_gpu;
  bool run_cpu;
  int g_snps;
  int g_flanking_snps;
  int g_people;
  int g_max_window;
  int g_max_haplotypes;
  int g_haplotype_mode;
  int g_likelihood_mode;
  float g_delta;
  float g_lambda;
  int geno_dim;
  int penetrance_matrix_size;
  //int g_max_region_size;
  // variables
  int g_markers;
  int g_haplotypes;
  int * g_active_haplotype;
  int * g_haplotype;
  float * g_old_frequency;
  float * g_frequency;
  int g_left_marker;
  int g_right_marker;
  int g_center_snp_start;
  int g_center_snp_end;
  float * penetrance_cache;
  float * logpenetrance_cache;
  int * haploid_arr;
  float * g_snp_penetrance;
  //float * g_region_snp_penetrance;
  //int g_region_snp_offset;
  int g_genotype_imputation;
  ofstream ofs_genotype_file;
  ofstream ofs_dosage_file;
  ofstream ofs_quality_file;
  ofstream ofs_debug_haplotype_file;
  ofstream ofs_posterior_file;
  string outfile_format;
  string outfile_geno;
  string outfile_dosages;
  string outfile_quality;
  string outfile_posterior;
  float * subject_haplotype_weight;
  float * frequency_cache;
  float * g_current_weight;
  float * g_weight;
  float * max_penetrance;
  #ifdef USE_GPU
  cl_int err;
  // commandQueue is used to carry out memory transfers and kernel launches
  // need to be initialized before any thing else
  cl::CommandQueue * commandQueue;
  cl::Context * context;
  vector<cl::Device> devices;
  cl::Program * program;
  cl::Kernel * kernel_simple;
  cl::Kernel * kernel_compute_weights;
  cl::Kernel * kernel_compute_weights_haploid;
  cl::Kernel * kernel_reduce_weights2;
  template<class T> cl::Buffer * createBuffer(int rw,int dim,const char * label);
  cl::Buffer * buffer_simple_in;
  cl::Buffer * buffer_simple_out;
  cl::Buffer * buffer_markers;
  cl::Buffer * buffer_haplotypes;
  cl::Buffer * buffer_haplotype;
  cl::Buffer * buffer_active_haplotype;
  cl::Buffer * buffer_subject_haplotype_weight;
  cl::Buffer * buffer_haplotype_weight;
  cl::Buffer * buffer_left_marker;
  cl::Buffer * buffer_frequency;
  cl::Buffer * buffer_penetrance_cache;
  cl::Buffer * buffer_logpenetrance_cache;
  cl::Buffer * buffer_frequency_cache;

  //cl::Buffer * buffer_region_snp_penetrance;
  //cl::Buffer * buffer_region_snp_offset;
  cl::Buffer * buffer_genotype_imputation;
  cl::Buffer * buffer_max_penetrance;
  cl::Buffer * buffer_haploid_arr;
  cl::Buffer * buffer_iteration;

  #endif

  string hapstr(int * haplotype, int marker_len);
  float compute_rsq(float * dosages,int stride,int index);
  float addLog(float log1, float log2);
private:
  // these set of C++ functions expose the OpenCL API to 
  // the FORTRAN program
  const char * imputation_software;

  //float * genotype_posterior;
  int * true_geno;
  float * true_maf;
  const char * chromosome;
  bool is_sex_chr;
  string infile_refhap;
  string infile_sex;
  string infile_polymorphisms;
  string infile_geno;
  string kernel_path;
  //int left_marker;
  #ifdef USE_GPU
  //cl_int err;
  // handles to the various kernels go here
  #endif
  
  int extractBit(packedhap_t * packed,int index);
  void injectBit(packedhap_t * packed,int index,int i);
  void testpacked();
  void unmarshall(ifstream & ifs_cache, float * outputvec,int len);

  // return true if converged
  bool check_mm_converged();
    
  // this functions are first called by the fortran wrapper
  void init_gpu_(int * platform_id, int * device_id);
  // just a simple test to see that the GPU is working
  void run_simple_(int * scaler, int * return_vec);
  // this function loads settings from the configuration file
  // and stores the values in global variables that are visible
  // to the fortran wrapper. It also loads the datasets in to memory.
  void load_constants_(int * model, int * people,int * snps, int * total_regions,int * haplotype_mode, int * flanking_snps,int * max_haplotypes,int * max_extended_haplotypes,int * platform_id, int * device_id, float * delta, float * lambda, int * i_geno_dim);

  // allocates data structures for the host and the GPU
  void init_buffers_(int * active_haplotype,int * max_window,int * max_haplotypes,int * max_extended_haplotypes,int * max_region_size,int * people, int * snps,  int * genotype_imputation,int * haplotypes, int * markers, int * window1,int * prev_left_marker, int * window2,int * prev_right_marker,float * frequency, float * weight, int * haplotype,int * flanking_snps);
  // allocates region specific data structures for the host and the GPU
  //void init_region_buffers_(int * regionstart,int * regionend);
  // allocates window specific data structures for the host and the GPU
  void init_window_buffers_();
  // allocates iteration specific data structures for the host and the GPU
  void init_iteration_buffers_();
  // functions for computing haplotype frequencies
  void compute_haplotype_weights_(int * iteration);

};

class GuidedMendelGPU:public MendelGPU{
public:
  GuidedMendelGPU(IO_manager *io);
  ~GuidedMendelGPU();
private:
  struct byHapCountDesc{
    bool operator()(const hapobj_t & p1,const hapobj_t & p2) const{
      // sort in descending order
      return(p1.hapcount>p2.hapcount );
    }
  };
  // class variables
  int g_informative_markers;
  int extended_markers;
  bool * informative_snp;
  queue<int> window_informative_snp_indices;
  int extended_haplotypes;
  int * extended_haplotype;
  int * extended_snp_mapping;
  //int * extended_center_dosage;
  float * extended_frequency;
  int * extended_root_mapping;
  std::tr1::unordered_map<string,int > ref_hap_counts;
  std::tr1::unordered_map<string,set<string> > full_hap_window;
  int ref_haplotypes;
  char * ref_haplotype;
  int g_max_extended_haplotypes;
  int packedhap_len;
  packedhap_t * packedhap;
  int packedextendedhap_len;
  packedhap_t * packedextendedhap;
  int * g_informative_haplotype;
  #ifdef USE_GPU
  cl::Buffer * buffer_packedhap;
  cl::Buffer * buffer_extended_snp_mapping;
  cl::Buffer * buffer_extended_haplotypes;
  cl::Buffer * buffer_extended_root_mapping;
  cl::Buffer * buffer_extended_frequency;
  cl::Buffer * buffer_center_snp_end;
  cl::Buffer * buffer_packedextendedhap;
  cl::Buffer * buffer_subject_posterior_prob_block;
  cl::Buffer * buffer_subject_genotype_block;
  cl::Buffer * buffer_subject_dosage_block;
  cl::Kernel * kernel_precompute_penetrance;
  cl::Kernel * kernel_precompute_penetrance_haploid;
  cl::Kernel * kernel_impute_genotype_guide;
  cl::Kernel * kernel_impute_genotype_guide_haploid;
  #endif
  void init_window();
  void finalize_window();
  void impute_genotypes();
  void load_datasets();
  void allocate_memory();
  void init_opencl();
  void compute_penetrance();
  // allocates data structures for the host and the GPU
  void init_buffers_(int * active_haplotype,int * max_window,int * max_haplotypes,int * max_extended_haplotypes,int * max_region_size,int * people, int * snps,  int * genotype_imputation,int * haplotypes, int * markers, int * window1,int * prev_left_marker, int * window2,int * prev_right_marker,float * frequency, float * weight, int * haplotype,int * flanking_snps);
  // functions for reference based haplotyping
  void parse_ref_haplotype();
  void copy_ref_haplotypes(int left_marker);
  // functions for penetrance calculations of SNP data
  // with ref haps
  void compress(int hapindex,int markers,int * haplotype);
  void decompress(int hapindex,int markers,int * haplotype);
  void precompute_penetrance();
  // genotype imputation with ref haplotypes
  void prep_impute_genotypes_guide();
  void impute_diploid_genotypes_guide(int  center_snp_start, 
  int  center_snp_end, int  center_snp_offset);
  void impute_haploid_genotypes_guide(int  center_snp_start, 
  int  center_snp_end, int  center_snp_offset);
  //utility functions
  void compress_extendedhap(int extendedhapindex,int begin,int end,int * extendedhaplotype);
};

class DenovoMendelGPU:public MendelGPU{
public:
  DenovoMendelGPU(IO_manager *io);
  virtual ~DenovoMendelGPU();
protected:
  struct byHapFreqAsc{
    bool operator()(const hapobj_t & p1,const hapobj_t & p2) const{
      // sort in ascenting order
      return(p1.hapfreq<p2.hapfreq);
    }
  };
  bool polymorphic_window;
  int g_prev_left_marker;
  int g_prev_right_marker;
  int * twin_hap_index;
  int * beyond_left_edge_dosage;
  list<int> occupied_hap_indices;
  list<int> free_hap_indices;
  int * right_edge_dosage;
  int * center_dosage;
  bool * window_polymorphisms;

  virtual void compute_penetrance()=0;
  virtual void allocate_memory();
  virtual void init_opencl();
  virtual void init_window();

  void finalize_window();
  void load_datasets();
  // functions for de novo haplotyping
  void impute_genotypes();
  void double_haplotypes();
  void prune_haplotypes_();
  // init opencl
  void init_window_opencl();
  // genotype imputation with no ref haplotypes
  void impute_haploid_genotypes_denovo_(int  center_snp, int  d_snp);
  void impute_diploid_genotypes_denovo_(int  center_snp, int  d_snp);
  // utility functions
  void debug_haplotypes(ostream & os);

  cl::Kernel * kernel_impute_genotype_denovo;
  cl::Kernel * kernel_impute_genotype_denovo_haploid;
  cl::Buffer * buffer_subject_posterior_prob;
  cl::Buffer * buffer_subject_genotype;
  cl::Buffer * buffer_subject_dosage;
  cl::Buffer * buffer_beyond_left_edge_dosage;
  cl::Buffer * buffer_center_snp;
  cl::Buffer * buffer_center_dosage;
  cl::Buffer * buffer_right_edge_dosage;
  cl::Buffer * buffer_twin_hap_index;
  cl::Buffer * buffer_prev_left_marker;
private:
};

class ReadParser;
class DenovoReadsMendelGPU:public DenovoMendelGPU{
public:
  DenovoReadsMendelGPU(IO_manager *io);
  ~DenovoReadsMendelGPU();
private:
  ReadParser *  parsers;
  // for BAM files
  static const int MAX_READS = 100;
  //int g_max_reads;
  int * g_reads_alleles;
  float * g_reads_errors;
  // for the BAM files, there may be redundant reads 
  int * g_unique_reads;
  int * g_unique_reads_counts;
  int * g_unique_reads_alleles;
  float * g_unique_reads_errors;
  float * g_unique_reads_logpenetrance;

  //void impute_genotypes();
  void allocate_memory();
  void init_opencl();
  void compute_penetrance();
  void init_window();

private:
  // functions for penetrance calculations of reads
  float get_bam_loglikelihood(int len,float *  hap1prob,float *  hap2prob);
  float log_half;
};

class DenovoGlfMendelGPU:public DenovoMendelGPU{
public:
  DenovoGlfMendelGPU(IO_manager *io);
  ~DenovoGlfMendelGPU();
private:
  // class variables
  #ifdef USE_GPU
  cl::Kernel * kernel_precompute_penetrance_fast;
  cl::Kernel * kernel_precompute_penetrance_fast_haploid;
  cl::Kernel * kernel_impute_penetrance;
  cl::Kernel * kernel_impute_penetrance_haploid;
  #endif
  // allocates data structures for the host and the GPU
  //void impute_genotypes();
  void allocate_memory();
  void init_opencl();
  void compute_penetrance();

  //void load_datasets();
  //void init_window();
  //void finalize_window();
  //void init_buffers_(int * active_haplotype,int * max_window,int * max_haplotypes,int * max_extended_haplotypes,int * max_region_size,int * people, int * snps,  int * genotype_imputation,int * haplotypes, int * markers, int * window1,int * prev_left_marker, int * window2,int * prev_right_marker,float * frequency, float * weight, int * haplotype,int * flanking_snps);
  // allocates region specific data structures for the host and the GPU
  //void init_region_buffers_(int * regionstart,int * regionend);
  // allocates window specific data structures for the host and the GPU
  void init_window_buffers_();
  // functions for penetrance calculations of SNP data
  // with no ref haps
  void precompute_penetrance_fast(int center_snp);
  void impute_penetrance_matrix();

};


// PUT INLINE FUNCTIONS HERE
inline char hapint2char(int i){
  return (char)(i+(int)'0'); 
}
  
inline int hapchar2int(char c){
  return (int)c-(int)'0';
}
