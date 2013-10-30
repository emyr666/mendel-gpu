#include<assert.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<map>
#include<set>
#include<stdlib.h>
#include<cstdio>
#include<stdio.h>
#include"faidx.h"
#include "sam.h"
#include"math.h"
#ifdef USE_GPU
#include<CL/cl.hpp>
#endif

using namespace std;

class ReadParser;

struct read_constants_t{
  static const uint MAX_MATRIX_WIDTH=20;
  static const uint MAX_SUPERREADS = 50;
  static const uint MAX_SUBREADS = 30;
  static const uint MAX_MATRIX_HEIGHT=MAX_SUPERREADS*MAX_SUBREADS;
  static const float EPSILON=1e-2;
};

struct tmpstruct_t{
  int beg, end;
  samfile_t *in;
  int cursor;
  ReadParser * parser;
};

struct pileup_key_t{
  //pileup_key_t(const pileup_key_t & old);
  pileup_key_t();
  ~pileup_key_t();
  string subject;
  int querypos;
};

struct pileup_val_t:public read_constants_t{
//  pileup_val_t(const pileup_val_t & old);
  pileup_val_t();
  ~pileup_val_t();
  int total_reads;
  int total_positions;
  //string read_names[MAX_SUPERREADS]; 
  int positions[MAX_MATRIX_WIDTH];
  int depths[MAX_SUPERREADS];
  int lengths[MAX_SUPERREADS];
  int alleles[MAX_SUPERREADS * MAX_MATRIX_WIDTH];
  float logmatch_qualities[MAX_SUPERREADS * MAX_MATRIX_WIDTH];
  float logmismatch_qualities[MAX_SUPERREADS * MAX_MATRIX_WIDTH];
  //int basequalities[MAX_SUPERREADS * MAX_SUBREADS * MAX_MATRIX_WIDTH];
  //string  mapqual_str[MAX_SUPERREADS];
  //string  basequal_str[MAX_SUPERREADS];
  //string sequence_str[MAX_SUPERREADS];
  //int  total_cigars[MAX_SUPERREADS];
  //int  cigar_ops[MAX_SUPERREADS][MAX_MATRIX_WIDTH];
  //int  cigar_len[MAX_SUPERREADS][MAX_MATRIX_WIDTH];
};

struct by_subject_pos{
  bool operator()(const pileup_key_t & a,const pileup_key_t & b){
    if (a.subject.compare(b.subject)!=0){
      return a.subject<b.subject;
    }else{
      return a.querypos<b.querypos;
    }
  }
};

typedef map<pileup_key_t,pileup_val_t,by_subject_pos> pileup_map_t;
typedef multimap<pileup_key_t,pileup_val_t,by_subject_pos> pileup_multimap_t;
typedef set<pileup_key_t,by_subject_pos> pileup_set_t;

class ReadParser:public read_constants_t{
public:
  ReadParser();
  //ReadParser(const ReadParser & old);
  ~ReadParser();
  void init(const char * bam, const char * pos_file);

  void extract_region(int subject_index,int snp_offset,int snps,bool clear_cache);
  void print_data();
  
  static void init_phred_lookup();
  // utility function
  static inline float phred2prob1(int phred){
    return 1.-(pow(10,(-phred/10.)));
  }
  void finalize(int left_marker);
  vector<int> found_positions;



  static float min_log_phred,max_log_phred;
  static float logprob_match_lookup[];
  static float logprob_mismatch_lookup[];
  uint matrix_rows;
  int * matrix_alleles;
  int * offset;
  int * depth;
  int * read_len;
  int ** map_quality;
  float ** logmatch_quality;// dimensions are superreads, and bases
  float ** logmismatch_quality;// dimensions are superreads, and bases
  //int *** base_quality;// dimensions are superreads, subreads, and bases
private:

  static bool debug; 
  //static bool use_cache;

  void parse_positions(const char * filename);
  void add_to_pileup_maps(pileup_key_t & key,pileup_val_t & val);
  //vector<int> positions;
  set<int> position_set;
  vector<int> position_vec;
  //pileup_set_t untyped_pileup_set;
  pileup_map_t pileup_map;
  pileup_multimap_t partial_pileup_map;
  int last_left_marker;
  bool pileups_found;
  int cursor;
  bam_index_t *idx;
  const char * bam_filename;
  //const char * variants_file;
  void make_partial_pileup(int lastfound_position,int current_position,int offset,pileup_key_t & key, pileup_val_t & val);
  
  
  

  tmpstruct_t * tmp;
  
  static int hex2int(string input);
  
  
  static void convert_hexstr(string mesg_read,int * qual_arr);

  // callback for bam_fetch()
  static int fetch_func(const bam1_t *b, void *data);
  // callback for bam_plbuf_init()
  static void process_region(const pileup_key_t & key,const pileup_val_t & val,ReadParser  *  parser);
  static int pileup_func(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data);
  static int pileup_func_old(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data);
  
  
  
  string get_region(int subject,int snpindex);
};

class MendelGPU;

class ReadPenetrance:public read_constants_t{
public:
  ReadPenetrance(MendelGPU * mendelgpu);
  ~ReadPenetrance();
  void prefetch_reads(int start_snp, int total_snps);
  void populate_read_matrices();
  void process_read_matrices();
  void init_opencl();
  void free_opencl();
protected:
  ReadParser *  parser;
  void process_read_matrices_opencl();
  #ifdef USE_GPU
  cl::Kernel * kernel_reads_compute_penetrance;
  cl::Kernel * kernel_reads_adjust_penetrance;
  cl::Buffer * buffer_read_alleles_mat;
  cl::Buffer * buffer_read_match_logmat;
  cl::Buffer * buffer_read_mismatch_logmat;
  cl::Buffer * buffer_mat_rows_by_subject;
  #endif
  // functions for penetrance calculations of reads
  float get_bam_loglikelihood(int len,float *  hap1prob,float *  hap2prob);
  float log_half;
  int compact_rows;
  int full_rows;
  int read_compact_matrix_size;
  int read_full_matrix_size;
  int * read_alleles_mat; // has superreads rows
  float * read_match_logmat; // has superreads*subreads rows
  float * read_mismatch_logmat;// has superreads*subreads rows
  int * mat_rows_by_subject;
private:
  MendelGPU * mendelgpu;
};
