#include"../read_penetrance.hpp"

bool ReadParser::debug = false;
vector<int> ReadParser::positions;

int ReadParser::hex2int(string input){
  int val;
  stringstream converter(input);
  converter>>hex>>val;
  return val;
}

float ReadParser::phred2prob(int phred){
  float prob = 1-(pow(10,(-phred/10)));
  if (prob<ReadParser::EPSILON) {
    prob = ReadParser::EPSILON;
    //cerr<<"A Retuning "<<prob<<endl;
  }
  else if ((1.-prob)<(1.-ReadParser::EPSILON)) {
    prob = 1.-ReadParser::EPSILON;
    //cerr<<"B Retuning "<<prob<<endl;
  }
  return prob;
}


// callback for bam_fetch()
int ReadParser::fetch_func(const bam1_t *b, void *data)
{
  bam_plbuf_t *buf = (bam_plbuf_t*)data;
  bam_plbuf_push(b, buf);
  return 0;
}

void ReadParser::convert_hexstr(string mesg_read,int * qual_arr){
  int mesg_len = mesg_read.length();
  for(int i=0;i<mesg_len/2;++i){
    string element = mesg_read.substr(i*2,2);
    qual_arr[i] = ReadParser::hex2int(element);
  }
}

// callback for bam_plbuf_init()
int ReadParser::pileup_func(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data)
{
  ReadParser *parser = (ReadParser *)data;
  tmpstruct_t *tmp = parser->tmp;
  
  //cerr<<"QueryPos: "<<querypos<<" "<<tmp->beg<<","<<tmp->end<<endl;
  if ((int)querypos >= tmp->beg && (int)querypos < tmp->end){
    if (debug) cerr<<"Chr: "<<tmp->in->header->target_name[tid]<<" pos: "<<querypos+1<<" reads: "<< n<<endl;
    int i=0;
    for(i=0;i<n;++i){
      uint genomic_pos = querypos + 1;
      if (debug) cerr<<" Read: "<<i<<endl;
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
      if (debug) cerr<<"  Name "<<name<<" at pos "<<pos<<" with cigar len "<<cigar_len<<endl;
      //cerr<<" auxiliary data len = "<<bam->l_aux<<endl;
      uint8_t * x1 = bam_aux_get(bam,"x1");
      if (x1){
        string mesg_read = bam_aux2Z(x1);
        //cerr<<"Message: "<<mesg_read<<" data: "<<bam->data[0]<<" of length: "<<mesg_read.length()<<endl;
        int mesg_len = mesg_read.length();
        parser->depth[parser->matrix_rows] = mesg_len/2;
        parser->map_quality[parser->matrix_rows] = new int[parser->depth[parser->matrix_rows]];
        //int read_qual_arr[parser->depth[parser->matrix_rows]];
        ReadParser::convert_hexstr(mesg_read,parser->map_quality[parser->matrix_rows]);
        //for(int i=0;i<parser->depth[parser->matrix_rows];++i) cerr<<" "<<parser->map_quality[parser->matrix_rows][i];
        //cerr<<endl;
      }
      uint8_t * y1 = bam_aux_get(bam,"y1");
      if (y1){
        string mesg_read = bam_aux2Z(y1);
        //cerr<<"Message2: "<<mesg_read<<" of length: "<<mesg_read.length()<<endl;
        int mesg_len = mesg_read.length();
        parser->read_len[parser->matrix_rows] = mesg_len/(2*parser->depth[parser->matrix_rows]);
        parser->base_quality[parser->matrix_rows] = new int[mesg_len/2];
        convert_hexstr(mesg_read,parser->base_quality[parser->matrix_rows]);
        //for(int i=0;i<mesg_len/2;++i) cerr<<" "<<parser->base_quality[parser->matrix_rows][i];
        //cerr<<endl;
      }

      string seq_str(bam->core.l_qseq,' ');
      //char * qseq = (char*)malloc(bam->core.l_qseq+1);
      uint8_t * seq = bam1_seq(bam);
      bool debug_seq = false;
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
        if (debug) cerr<<"  Current cigar len "<<c1<<" of type ";
        if (cop==BAM_CMATCH ){
          string read_substr(seq_str,cursor,1);
          assert(parser->matrix_rows<MATRIX_HEIGHT);
          assert(offset+cursor<MATRIX_WIDTH);
          if (read_substr.compare("=")==0){
            if (debug) cerr<<" MATCH\n";
            parser->matrix_alleles[parser->matrix_rows * MATRIX_WIDTH+offset+cursor] = 0;
          }else if (read_substr.compare("N")==0){
            if (debug) cerr<<" MISMATCH\n";
            parser->matrix_alleles[parser->matrix_rows * MATRIX_WIDTH+offset+cursor] = 1;
          }else{
            cerr<<"Unknown read_substr is "<<read_substr<<endl;
            exit(0);
          }
          if (debug)  cerr<<"  Inserted "<< parser->matrix_alleles[parser->matrix_rows * MATRIX_WIDTH+offset+cursor]<<" into row "<<parser->matrix_rows<<" and col "<<offset+cursor<<endl;
          genomic_pos+=c1;
        }else if (cop==BAM_CPAD){
          if (debug) cerr<<"  CPAD\n";
          genomic_pos+=c1;
          ++cursor;
        }
        if (debug) cerr<<"  Genomic position is now "<<genomic_pos<<endl;
      }
      parser->offset[parser->matrix_rows] = offset;
      ++parser->matrix_rows;
    }
    // the offset is moved by one column for next SNP
    ++tmp->cursor;
    if (debug) cerr<<"  query position is now "<<querypos<<endl;
    parser->found_positions.push_back(querypos+1); 
    if (debug) cerr<<"Cursor is now "<<tmp->cursor<<endl;
  }else{
     //cerr<<"No reads found in this interval\n";
  }
  return 0;
}

void ReadParser::extract_region(const char * chrom,int offset,int snps){
  if (debug) cerr<<"Selected chromosome "<<chrom<<endl;
  matrix_rows = 0;
  tmp->cursor = 0;
  for(uint i=0;i<MATRIX_WIDTH*MATRIX_HEIGHT;++i){
    matrix_alleles[i] = 9;
  }
  found_positions.clear();
  for(int snpindex = offset;snpindex<offset+snps;++snpindex){
    string regionstr = get_region(chrom,snpindex);
    const char * region = regionstr.data();
    //const char * region = "20:60749-60749";
    //int ref;
    // extract the region next from BAM file
    int ref;
    if (debug) cerr<<"Querying SNP position:  "<<region<<endl;
    bam_parse_region(tmp->in->header, region, &ref,
                    &tmp->beg, &tmp->end); // parse the region
    if (ref < 0) {
      fprintf(stderr, "Invalid region %s\n", region);
      exit(1);
    }
    bam_plbuf_t *buf = bam_plbuf_init(pileup_func, this); // initialize pileup
    bam_fetch(tmp->in->x.bam, idx, ref, tmp->beg, tmp->end, buf, fetch_func);
    bam_plbuf_push(0, buf); // finalize pileup
    bam_plbuf_destroy(buf);
    if (debug) cerr<<"Matrix row is now "<<matrix_rows<<endl;
  }
}

void ReadParser::parse_positions(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Cannot open variant file "<<filename<<endl;
    exit(1);
  }
  string line;
  getline(ifs,line);
  while(getline(ifs,line)){
    istringstream iss(line);
    string chr,rs;
    int pos,morgan;
    iss>>chr>>rs>>morgan>>pos;
    //cerr<<"Position is "<<pos<<endl;
    ReadParser::positions.push_back(pos);
  }
  if (debug) cerr<<"There are "<<ReadParser::positions.size()<<" positions available\n";
  ifs.close();
}

string ReadParser::get_region(const char * chrom,int snpindex){
  ostringstream oss;
  if (debug) cerr<<"Fetching region\n";
  oss<<chrom<<":"<<ReadParser::positions[snpindex]<<"-"<<ReadParser::positions[snpindex];
  string ossstr = oss.str();
  if (debug) cerr<<"Region is "<<ossstr<<endl;
  return ossstr;
}

ReadParser::ReadParser(const ReadParser & old){
  //cerr<<"Copy constructor being called!\n";
  init(old.bam_filename);
}

ReadParser::ReadParser(){
}

void ReadParser::init(const char * bam_filename){
  this->bam_filename = bam_filename;
  offset = new int[MATRIX_HEIGHT];
  depth = new int[MATRIX_HEIGHT];
  read_len = new int[MATRIX_HEIGHT];
  map_quality = new int * [MATRIX_HEIGHT];
  base_quality = new int * [MATRIX_HEIGHT];
  matrix_alleles = new int[MATRIX_WIDTH*MATRIX_HEIGHT];
  tmp = new tmpstruct_t;
  tmp->beg = 0; tmp->end = 0x7fffffff;
  tmp->in = samopen(bam_filename, "rb", 0);
  tmp->parser = this;
  if (tmp->in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", bam_filename);
    exit(1);
  }
  idx = bam_index_load(bam_filename); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM indexing file is not available.\n");
    exit(1);
  }
}

ReadParser::~ReadParser(){
  bam_index_destroy(idx);
  samclose(tmp->in);
  delete tmp;
  //cerr<<"Destructed\n";
}

void ReadParser::print_data(){
  cout<<"SNP POSITIONS:\n";
  for(uint i=0;i<found_positions.size();++i){
    cout<<found_positions[i]<<endl;
  }
  cout<<"HAPLOTYPE MATRIX:\n";
  for(uint i=0;i<matrix_rows;++i){
    for(uint j=0;j<MATRIX_WIDTH;++j){
      //if (j) cout<<" ";
      cout<<matrix_alleles[i*MATRIX_WIDTH+j];
    }
    cout<<endl;
  }
  cout<<"OFFSET\tDEPTH\tREAD_LENGTH\tMAP_QUAL\tBASE_QUAL\tBASE_PROB\n";
  for(uint i=0;i<matrix_rows;++i){
    cout<<offset[i];
    cout<<"\t";
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
    cout<<"\t";
    for(int j=0;j<depth[i];++j){
      if (j) cout<<";";
      for(int k=0;k<read_len[i];++k){
        if (k) cout<<",";
        cout<<phred2prob(base_quality[i][j*read_len[i]+k]);
      }
    }
    cout<<endl;
  }
}

int main_readpen(int argc, char *argv[])
{
  try{
    if (argc < 6) {
      fprintf(stderr, "Usage: <in.bam> <variants> <chr> <snp_start> <total_snps> \n");
      return 1;
    }
    int arg = 0;
    const char * bam_filename = argv[++arg];
    const char * variants_file = argv[++arg];
    const char * chrom = argv[++arg];
    int offset = atoi(argv[++arg]);
    int total_regions = atoi(argv[++arg]);
  
    //ReadParser * parser = new ReadParser(bam_filename,variants_file);
    ReadParser::parse_positions(variants_file);
    ReadParser * parser = new ReadParser();
    parser->init(bam_filename);
    parser->extract_region(chrom,offset,total_regions);
      //string region = get_region(chrom,offset+i);
      //const char * region = "20:60522-60522";
    parser->print_data();
    delete parser;
    return 0;
  }catch (const char * & mesg){
    cerr<<"Exception caught of message "<<mesg<<endl;
  }
}

