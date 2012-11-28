
class MendelGPU{
public:
  MendelGPU();
  ~MendelGPU();
  void load_constants_(int * model, int * people,int * snps, int * total_regions,int * haplotype_mode, int * flanking_snps,int * max_haplotypes,int * max_extended_haplotypes,int * platform_id, int * device_id, float * delta, float * lambda, int * i_geno_dim);
  void init_gpu_(int * platform_id, int * device_id);
  void run_simple_(int * scaler, int * return_vec);
  void read_stream_(int * people, int * snps, float * snp_penetrance);
  void init_buffers_(int * active_haplotype,int * max_window,int * max_haplotypes,int * max_extended_haplotypes,int * max_region_size,int * people, int * snps, float * snp_penetrance, int * genotype_imputation,int * haplotypes, int * markers, int * window1,int * prev_left_marker, float * frequency, float * weight, int * haplotype,int * flanking_snps);
  void init_region_buffers_(int * regionstart,int * regionend);
  void double_haplotypes_();
  void prune_haplotypes_();
  void copy_ref_haplotypes_();
  void init_window_buffers_();
  void init_iteration_buffers_();
  void precompute_penetrance_fast_();
  void impute_penetrance_matrix_();
  void precompute_penetrance_();
  void compute_haplotype_weights_(int * iteration);
  void prep_impute_genotypes_guide_();
  void impute_diploid_genotypes_guide_(int * center_snp_start, 
  int * center_snp_end, int * center_snp_offset);
  void impute_haploid_genotypes_denovo_(int * center_snp, int * d_snp);
  void impute_haploid_genotypes_guide_(int * center_snp_start, 
  int * center_snp_end, int * center_snp_offset);
  void impute_diploid_genotypes_denovo_(int * center_snp, int * d_snp);
private:
  // these set of C++ functions expose the OpenCL API to 
  // the FORTRAN program
  float epsilon;
  static const int HAPLOTYPE_MODE_DENOVO = 0;
  static const int HAPLOTYPE_MODE_GUIDE = 1;
  static const int LIKELIHOOD_MODE_GENOTYPES = 0;
  static const int LIKELIHOOD_MODE_READS = 1;
  float logpenetrance_threshold;
  const char * imputation_software;
  int g_haplotype_mode;
  int g_likelihood_mode;
  bool debug_mode;
  bool run_gpu;
  bool run_cpu;
  bool debug_truth;
  int * g_haplotype;
  int * g_informative_haplotype;
  int extended_markers;
  int extended_haplotypes;
  int * extended_haplotype;
  int * extended_snp_mapping;
  //int * extended_center_dosage;
  float * extended_frequency;
  int * extended_root_mapping;
  bool * informative_snp;
  int packedhap_len;
  int packedextendedhap_len;
  packedhap_t * packedhap;
  packedhap_t * packedextendedhap;
  int g_max_region_size;
  int g_flanking_snps;
  ofstream ofs_genotype_file;
  ofstream ofs_dosage_file;
  ofstream ofs_quality_file;
  ofstream ofs_debug_haplotype_file;
  ofstream ofs_posterior_file;
  float * g_snp_penetrance;
  float * g_region_snp_penetrance;
  int g_region_snp_offset;
  //float * genotype_posterior;
  float * g_frequency;
  float * g_current_weight;
  float * g_weight;
  float * subject_haplotype_weight;
  float * logpenetrance_cache;
  float * penetrance_cache;
  float * frequency_cache;
  float * max_penetrance;
  int penetrance_matrix_size;
  int g_genotype_imputation;
  char * ref_haplotype;
  int ref_haplotypes;
  list<int> occupied_hap_indices;
  list<int> free_hap_indices;
  int * twin_hap_index;
  int * beyond_left_edge_dosage;
  int * right_edge_dosage;
  int * center_dosage;
  int * g_active_haplotype;
  int * true_geno;
  int geno_dim;
  float * true_maf;
  int * haploid_arr;
  bool is_sex_chr;
  string FORMAT_DEFAULT;
  string FORMAT_MEC;
  string infile_refhap;
  string infile_sex;
  string infile_geno;
  string outfile_format;
  string outfile_geno;
  string outfile_dosages;
  string outfile_quality;
  string outfile_posterior;
  string kernel_path;
  //int left_marker;
  int g_max_haplotypes,g_max_extended_haplotypes,g_max_window,g_people,g_snps;
  int *g_prev_left_marker,*g_left_marker,*g_markers,*g_haplotypes;
  int g_informative_markers;

  std::tr1::unordered_map<string,set<string> > full_hap_window;
  std::tr1::unordered_map<string,int > ref_hap_counts;
  
  #ifdef USE_GPU
  cl_int err;
  // commandQueue is used to carry out memory transfers and kernel launches
  // need to be initialized before any thing else
  cl::CommandQueue * commandQueue;
  cl::Context * context;
  vector<cl::Device> devices;
  // handles to the various kernels go here
  cl::Program * program;
  cl::Kernel * kernel_simple;
  cl::Kernel * kernel_compute_weights;
  cl::Kernel * kernel_compute_weights_haploid;
  cl::Kernel * kernel_reduce_weights2;
  cl::Kernel * kernel_impute_genotype_denovo;
  cl::Kernel * kernel_impute_genotype_denovo_haploid;
  cl::Kernel * kernel_impute_genotype_guide;
  cl::Kernel * kernel_impute_genotype_guide_haploid;
  cl::Kernel * kernel_precompute_penetrance;
  cl::Kernel * kernel_precompute_penetrance_haploid;
  cl::Kernel * kernel_impute_penetrance;
  cl::Kernel * kernel_impute_penetrance_haploid;
  cl::Kernel * kernel_precompute_penetrance_fast;
  cl::Kernel * kernel_precompute_penetrance_fast_haploid;
  cl::Buffer * buffer_simple_in;
  cl::Buffer * buffer_simple_out;
  cl::Buffer * buffer_markers;
  cl::Buffer * buffer_haplotypes;
  cl::Buffer * buffer_extended_haplotypes;
  cl::Buffer * buffer_prev_left_marker;
  cl::Buffer * buffer_left_marker;
  cl::Buffer * buffer_haplotype;
  cl::Buffer * buffer_region_snp_penetrance;
  cl::Buffer * buffer_region_snp_offset;
  cl::Buffer * buffer_frequency;
  cl::Buffer * buffer_extended_frequency;
  cl::Buffer * buffer_subject_haplotype_weight;
  cl::Buffer * buffer_haplotype_weight;
  cl::Buffer * buffer_center_snp;
  cl::Buffer * buffer_subject_genotype;
  cl::Buffer * buffer_subject_genotype_block;
  cl::Buffer * buffer_subject_dosage;
  cl::Buffer * buffer_subject_dosage_block;
  cl::Buffer * buffer_subject_posterior_prob;
  cl::Buffer * buffer_subject_posterior_prob_block;
  cl::Buffer * buffer_logpenetrance_cache;
  cl::Buffer * buffer_max_penetrance;
  cl::Buffer * buffer_genotype_imputation;
  cl::Buffer * buffer_active_haplotype;
  cl::Buffer * buffer_twin_hap_index;
  cl::Buffer * buffer_penetrance_cache;
  cl::Buffer * buffer_frequency_cache;
  cl::Buffer * buffer_beyond_left_edge_dosage;
  cl::Buffer * buffer_right_edge_dosage;
  cl::Buffer * buffer_center_dosage;
  cl::Buffer * buffer_center_snp_start;
  cl::Buffer * buffer_center_snp_end;
  cl::Buffer * buffer_center_snp_offset;
  //cl::Buffer * buffer_extended_center_dosage;
  cl::Buffer * buffer_extended_root_mapping;
  cl::Buffer * buffer_extended_snp_mapping;
  cl::Buffer * buffer_packedhap;
  cl::Buffer * buffer_packedextendedhap;
  cl::Buffer * buffer_iteration;
  cl::Buffer * buffer_haploid_arr;
  #endif
  
  struct hapobj_t{
    string hapstr;
    int hapcount;
    float hapfreq;
    int array_index;
    set<string> extended_set;
  };
  
  struct byHapFreqAsc{
    bool operator()(const hapobj_t & p1,const hapobj_t & p2) const{
      // sort in ascenting order
      return(p1.hapfreq<p2.hapfreq);
    }
  };
  
  struct byHapCountDesc{
    bool operator()(const hapobj_t & p1,const hapobj_t & p2) const{
      // sort in descending order
      return(p1.hapcount>p2.hapcount );
    }
  };
  
  struct hapstr_pair_t {
    string condensed;
    set<string> extended_set;
  };
  
  struct byUniqueCondensed{
    bool operator()(const hapstr_pair_t & p1,const hapstr_pair_t & p2) const{
      // sort in ascenting order
      return(p1.condensed<p2.condensed);
    }
  };
  
  
  
  float compute_rsq(float * dosages,int stride, int index);
  string hapstr(int * haplotype, int marker_len);
  int extractBit(packedhap_t * packed,int index);
  void decompress(int hapindex,int markers,int * haplotype);
  void injectBit(packedhap_t * packed,int index,int i);
  void compress(int hapindex,int markers,int * haplotype);
  void compress_extendedhap(int extendedhapindex,int begin,int end,int * extendedhaplotype);
  void testpacked();
  void debug_haplotypes(ostream & os);
  void parse_ref_haplotype();
  void unmarshall(ifstream & ifs_cache, float * outputvec,int len);
};

// PUT INLINE FUNCTIONS HERE
inline char hapint2char(int i){
  return (char)(i+(int)'0'); 
}
  
inline int hapchar2int(char c){
  return (int)c-(int)'0';
}
