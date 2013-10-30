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

bam_index_t *idx;

const uint MATRIX_WIDTH=10;
const uint MATRIX_HEIGHT=20;

int * depth;
int * read_len;
int * matrix_alleles;
int * matrix_quality;
int ** map_quality;
int ** base_quality;
uint matrix_rows;
vector<int> positions;
vector<int> found_positions;


struct tmpstruct_t{
  int beg, end;
  samfile_t *in;
  faidx_t * fasta_index;
  int cursor;
};

inline void int2hex(int input, ostream & output){
  output<<hex<<input; 
}

int hex2int(string input){
  int val;
  stringstream converter(input);
  converter>>hex>>val;
  return val;
}

static inline float phred2prob(int phred){
  return 1-(pow(10,(-phred/10)));
}

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
  bam_plbuf_t *buf = (bam_plbuf_t*)data;
  bam_plbuf_push(b, buf);
  return 0;
}

void convert_hexstr(string mesg_read,int * qual_arr){
  int mesg_len = mesg_read.length();
  for(int i=0;i<mesg_len/2;++i){
    string element = mesg_read.substr(i*2,2);
    qual_arr[i] = hex2int(element);
  }
}

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data)
{
  tmpstruct_t *tmp = (tmpstruct_t*)data;
  //cerr<<"QueryPos: "<<querypos<<" "<<tmp->beg<<","<<tmp->end<<endl;
  if ((int)querypos >= tmp->beg && (int)querypos < tmp->end){
    cerr<<"Chr: "<<tmp->in->header->target_name[tid]<<" pos: "<<querypos+1<<" reads: "<< n<<endl;
    int i=0;
    for(i=0;i<n;++i){
      uint genomic_pos = querypos + 1;
      cerr<<" Read: "<<i<<endl;
      //printf("  indel status: %d\n", pl[i].indel); 
      //printf("  level: %d\n", pl[i].level); 
      //printf("  is_del: %d\n", pl[i].is_del); 
      //printf("  is_head: %d\n", pl[i].is_head); 
      //printf("  is_tail: %d\n", pl[i].is_tail); 
      const bam1_t * bam = pl[i].b;
      char * name = bam1_qname(bam); int pos = bam->core.pos+1;
      uint32_t* cigar = bam1_cigar(bam);
      int cigar_len = bam_cigar2qlen(&(bam->core),cigar);
      //tmp->snp.pileup_pos = pl[i].qpos;
      //int pileup_pos = pos+pl[i].qpos;
      cerr<<"  Name "<<name<<" at pos "<<pos<<" with cigar len "<<cigar_len<<endl;
      //cerr<<" auxiliary data len = "<<bam->l_aux<<endl;
      uint8_t * x1 = bam_aux_get(bam,"x1");
      if (x1){
        string mesg_read = bam_aux2Z(x1);
        //cerr<<"Message: "<<mesg_read<<" data: "<<bam->data[0]<<" of length: "<<mesg_read.length()<<endl;
        int mesg_len = mesg_read.length();
        depth[matrix_rows] = mesg_len/2;
        map_quality[matrix_rows] = new int[depth[matrix_rows]];
        //int read_qual_arr[depth[matrix_rows]];
        convert_hexstr(mesg_read,map_quality[matrix_rows]);
        //for(int i=0;i<depth[matrix_rows];++i) cerr<<" "<<map_quality[matrix_rows][i];
        //cerr<<endl;
      }
      uint8_t * y1 = bam_aux_get(bam,"y1");
      if (y1){
        string mesg_read = bam_aux2Z(y1);
        //cerr<<"Message2: "<<mesg_read<<" of length: "<<mesg_read.length()<<endl;
        int mesg_len = mesg_read.length();
        read_len[matrix_rows] = mesg_len/(2*depth[matrix_rows]);
        base_quality[matrix_rows] = new int[mesg_len/2];
        convert_hexstr(mesg_read,base_quality[matrix_rows]);
        //for(int i=0;i<mesg_len/2;++i) cerr<<" "<<base_quality[matrix_rows][i];
        //cerr<<endl;
      }

      string seq_str(bam->core.l_qseq,' ');
      //char * qseq = (char*)malloc(bam->core.l_qseq+1);
      uint8_t * seq = bam1_seq(bam);
      bool debug_seq = true;
      if (debug_seq) cerr<<"  Query sequence:"; 
      int j;
      for(j=0;j<bam->core.l_qseq;++j){
        int v = bam1_seqi(seq,j);
        seq_str[j] =  bam_nt16_rev_table[v];
        if (debug_seq)cerr<<" "<<j<<":"<<seq_str[j];
      }
      //qseq[j] = 0;
      if (debug_seq) cerr<<endl;
      uint offset = tmp->cursor;
      uint cursor = 0;
      for(int k=0;k<bam->core.n_cigar;++k){
        int cop = cigar[k] & BAM_CIGAR_MASK;
        int c1 = cigar[k] >> BAM_CIGAR_SHIFT;
        cerr<<"  Current cigar len "<<c1<<" of type ";
        if (cop==BAM_CMATCH ){
          string read_substr(seq_str,cursor,1);
          assert(matrix_rows<MATRIX_HEIGHT);
          assert(offset+cursor<MATRIX_WIDTH);
          if (read_substr.compare("=")==0){
            cerr<<" MATCH\n";
            matrix_alleles[matrix_rows * MATRIX_WIDTH+offset+cursor] = 0;
          }else if (read_substr.compare("N")==0){
            cerr<<" MISMATCH\n";
            matrix_alleles[matrix_rows * MATRIX_WIDTH+offset+cursor] = 1;
          }else{
            cerr<<"Unknown read_substr is "<<read_substr<<endl;
            exit(0);
          }
          cerr<<"  Inserted "<< matrix_alleles[matrix_rows * MATRIX_WIDTH+offset+cursor]<<" into row "<<matrix_rows<<" and col "<<offset+cursor<<endl;
          genomic_pos+=c1;
        }else if (cop==BAM_CPAD){
          cerr<<"  CPAD\n";
          genomic_pos+=c1;
          ++cursor;
        }
        cerr<<"  Genomic position is now "<<genomic_pos<<endl;
      }
      ++matrix_rows;
    }
    // the offset is moved by one column for next SNP
    ++tmp->cursor;
    found_positions.push_back(querypos+1); 
    cerr<<"Cursor is now "<<tmp->cursor<<endl;
  }else{
     //cerr<<"No reads found in this interval\n";
  }
  return 0;
}

void extract_region(const char * region,tmpstruct_t & tmp){
  //const char * region = "20:60749-60749";
  //int ref;
  // extract the region next from BAM file
  bam_plbuf_t *buf;
  int ref;
  cerr<<"Querying SNP position: "<<region<<endl;
  bam_parse_region(tmp.in->header, region, &ref,
                  &tmp.beg, &tmp.end); // parse the region
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s\n", region);
    exit(1);
  }
  buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup
  bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);
  bam_plbuf_push(0, buf); // finalize pileup
  bam_plbuf_destroy(buf);
}

void parse_positions(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Cannot open variant file "<<filename<<endl;
    exit(1);
  }
  string line;
  getline(ifs,line);
  while(getline(ifs,line)){
    istringstream iss(line);
    string chr;
    int pos;
    iss>>chr>>pos;
    positions.push_back(pos);
  }
  ifs.close();
}

string get_region(const char * chrom,int snpindex){
  ostringstream oss;
  oss<<chrom<<":"<<positions[snpindex]<<"-"<<positions[snpindex];
  string ossstr = oss.str();
  cerr<<"Regions is "<<ossstr<<endl;
  return ossstr;
}

int main(int argc, char *argv[])
{
  string test="1A";
  //int x=hex2int(test);
  //cout<<"int: "<<x<<endl;
  //cout<<"hex: ";
  //int2hex(x,cout);
  //cout<<endl;
  //exit(0);
  try{
    tmpstruct_t tmp;
    if (argc < 6) {
      fprintf(stderr, "Usage: <in.bam> <variants> <chr> <snp_start> <total_snps> \n");
      return 1;
    }
    depth = new int[MATRIX_HEIGHT];
    read_len = new int[MATRIX_HEIGHT];
    map_quality = new int * [MATRIX_HEIGHT];
    base_quality = new int * [MATRIX_HEIGHT];
    matrix_alleles = new int[MATRIX_WIDTH*MATRIX_HEIGHT];
    matrix_quality = new int[MATRIX_WIDTH*MATRIX_HEIGHT];
    matrix_rows = 0;
    for(uint i=0;i<MATRIX_WIDTH*MATRIX_HEIGHT;++i){
      matrix_alleles[i] = matrix_quality[i] = 9;
    }
    int arg = 0;
    const char * bam_filename = argv[++arg];
    const char * variants_file = argv[++arg];
    parse_positions(variants_file);
    const char * chrom = argv[++arg];
    tmp.beg = 0; tmp.end = 0x7fffffff;
    tmp.in = samopen(bam_filename, "rb", 0);
    if (tmp.in == 0) {
      fprintf(stderr, "Fail to open BAM file %s\n", bam_filename);
      return 1;
    }
    idx = bam_index_load(bam_filename); // load BAM index
    if (idx == 0) {
      fprintf(stderr, "BAM indexing file is not available.\n");
      return 1;
    }
    int offset = atoi(argv[++arg]);
    int total_regions = atoi(argv[++arg]);
    tmp.cursor = 0;
    for(int i=0;i<total_regions;++i){
      cerr<<"Getting region\n";
      string region = get_region(chrom,offset+i);
      //const char * region = "20:60522-60522";
      cerr<<"Region is "<<region<<endl;
      extract_region(region.data(),tmp);
      cerr<<"Matrix row is now "<<matrix_rows<<endl;
    }
    cerr<<"SNP POSITIONS:\n";
    for(uint i=0;i<found_positions.size();++i){
      cout<<found_positions[i]<<endl;
    }
    cerr<<"HAPLOTYPE MATRIX:\n";
    for(uint i=0;i<matrix_rows;++i){
      for(uint j=0;j<MATRIX_WIDTH;++j){
        if (j) cout<<" ";
        cout<<matrix_alleles[i*MATRIX_WIDTH+j];
      }
      cout<<endl;
    }
    cout<<"DEPTH\tREAD_LENGTH\tMAP_QUAL\tBASE_QUAL\n";
    for(uint i=0;i<matrix_rows;++i){
      cout<<depth[i];
      cout<<"\t";
      cout<<read_len[i];
      cout<<"\t";
      for(int j=0;j<depth[i];++j){
        if (j) cout<<",";
        cout<<map_quality[i][j];
      }
      cout<<"\t";
      for(int j=0;j<depth[i];++j){
        if (j) cout<<";";
        for(int k=0;k<read_len[i];++k){
          if (k) cout<<",";
          cout<<base_quality[i][j*read_len[i]+k];
        }
      }
      cout<<endl;

    }
    bam_index_destroy(idx);
    samclose(tmp.in);
    return 0;
  }catch (const char * & mesg){
    cerr<<"Exception caught of message "<<mesg<<endl;
  }
}

