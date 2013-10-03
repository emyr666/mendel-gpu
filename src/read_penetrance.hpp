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
  static float phred2prob(int phred);
  // utility function
  //static inline float phred2prob(int phred){
   // return 1-(pow(10,(-phred/10)));
  //}
  vector<int> found_positions;
  static const float EPSILON=1e-2;
  static const uint MATRIX_WIDTH=20;
  static const uint MATRIX_HEIGHT=200;
  uint matrix_rows;
  int * matrix_alleles;
  int * offset;
  int * depth;
  int * read_len;
  int ** map_quality;
  int ** base_quality;
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

