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

using namespace std;

class ReadParser;

struct tmpstruct_t{
  int beg, end;
  samfile_t *in;
  int cursor;
  ReadParser * parser;
};

class ReadParser{
public:
  ReadParser();
  ReadParser(const ReadParser & old);
  ~ReadParser();
  void init(const char * bam);

  static void parse_positions(const char * filename);
  void extract_region(const char * chrom,int snp_offset,int snps);
  void print_data();
  
  static void init_phred_lookup(int max_width);
  // utility function
  static inline float phred2prob1(int phred){
    return 1-(pow(10,(-phred/10)));
  }
  vector<int> found_positions;
  static uint max_matrix_width;
  static uint max_matrix_height;
  static const uint MAX_SUPERREADS = 10;
  static const uint MAX_SUBREADS = 50;
  static const uint MAX_READ_LEN = 10;
  static const float EPSILON=1e-2;
  static float min_log_phred,max_log_phred;
  static float logprob_match_lookup[];
  static float logprob_mismatch_lookup[];
  uint matrix_rows;
  int * matrix_alleles;
  int * offset;
  int * depth;
  int * read_len;
  int ** map_quality;
  int *** base_quality;// dimensions are superreads, subreads, and bases
private:
  static bool debug; 
  bam_index_t *idx;
  const char * bam_filename;
  //const char * variants_file;
  
  
  static vector<int> positions;
  

  tmpstruct_t * tmp;
  
  static int hex2int(string input);
  
  
  static void convert_hexstr(string mesg_read,int * qual_arr);

  // callback for bam_fetch()
  static int fetch_func(const bam1_t *b, void *data);
  // callback for bam_plbuf_init()
  static int pileup_func(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data);
  
  
  
  string get_region(const char * chrom,int snpindex);
};

