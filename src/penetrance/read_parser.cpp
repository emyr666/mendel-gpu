#include"../read_utilities.hpp"
#include"../cl_constants.h"

bool ReadParser::debug = false;
//bool ReadParser::use_cache = true;

//vector<int> ReadParser::positions;
//map<pileup_key_t,pileup_val_t,by_subject_pos> ReadParser::pileup_map;
float ReadParser::min_log_phred = log(ReadParser::EPSILON);
float ReadParser::max_log_phred = log(1.-ReadParser::EPSILON);
float ReadParser::logprob_match_lookup[100];
float ReadParser::logprob_mismatch_lookup[100];
int test_subject = TEST_SUBJECT;
string test_subject_str=TEST_SUBJECT_STR;

int ReadParser::hex2int(string input){
  int val;
  stringstream converter(input);
  converter>>hex>>val;
  return val;
}

inline float get_bounded(float val){
  float logmatch = log(val);
  if (logmatch<ReadParser::min_log_phred) {
    logmatch = ReadParser::min_log_phred;
  } else if (val>ReadParser::max_log_phred){
    logmatch = ReadParser::max_log_phred;
  }
  return logmatch;
}

void ReadParser::init_phred_lookup(){
  for(int phred=0;phred<100;++phred){
    float match_prob = 1.-(pow(10,(-phred/10.)));
    logprob_match_lookup[phred] =get_bounded(match_prob);
    float mismatch_prob = 1.-match_prob;
    logprob_mismatch_lookup[phred] =get_bounded(1.-match_prob);
    //cerr<<"phred "<<phred<<" has match value "<<match_prob<<","<<logprob_match_lookup[phred]
    //<<" and mismatch value "<<mismatch_prob<<","<<logprob_mismatch_lookup[phred]<<endl;
  }
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

void ReadParser::process_region(const pileup_key_t & key,const pileup_val_t & val,ReadParser  *  parser)
{
  uint genomic_pos = key.querypos;
  if ( parser->cursor < parser->MAX_MATRIX_WIDTH && 
       parser->matrix_rows<parser->MAX_SUPERREADS){
    bool debug = key.subject.compare(test_subject_str)==0;
    if (debug) cerr<<"For populating matrix, subject: "<<key.subject<<" pos: "<<genomic_pos<<" reads: "<<val.total_reads<<endl;
    for(int i=0;i<val.total_reads;++i){
      parser->read_len[parser->matrix_rows] = val.lengths[i];
      if (val.lengths[i]<1){
        if (debug) cerr<<" Read: "<<i<<" is skipped for being 0 length\n";
      }else{
        if (debug) cerr<<" Read: "<<i<<" at matrix row: "<<parser->matrix_rows<<endl;
        parser->depth[parser->matrix_rows] = val.depths[i];
        if(debug) cerr<<"  Read depth,length is "<<parser->depth[parser->matrix_rows]<<","<<parser->read_len[parser->matrix_rows]<<endl;
        //if (debug) cerr<<"MAX MATRIX WIDTH is "<<MAX_MATRIX_WIDTH<<endl;
        for(int j=0;j<val.lengths[i];++j){
          int matrix_col = j+parser->cursor;
          if (val.offsets[i]+j < MAX_VECTOR_LENGTH){
            if (matrix_col<MAX_MATRIX_WIDTH){
              parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+matrix_col] = val.alleles_compact[val.offsets[i]+j];
              if (debug) cerr<<"  Allele at row,col "<<parser->matrix_rows<<","<<matrix_col<<" is "<<parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+matrix_col]<<endl;
            }else{
              if (debug) cerr<<"matrix col is "<<matrix_col<<endl;
            }
            parser->logmatch_quality[parser->matrix_rows][j] = val.logmatch_qualities_compact[val.offsets[i]+j];
            parser->logmismatch_quality[parser->matrix_rows][j] = val.logmismatch_qualities_compact[val.offsets[i]+j];
          }
        } 
        //  for(int j=0;j<parser->read_len[parser->matrix_rows];++j){
        //    int col = parser->cursor+j;
        //    if(col<MAX_MATRIX_WIDTH){
        //      parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+col] = val.alleles[i*MAX_MATRIX_WIDTH+j];
        //      if (debug) cerr<<"  Allele at row,col "<<parser->matrix_rows<<","<<col<<" is "<<parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+col]<<endl;
        //    }
        //   }
        for(int k=0;k<parser->read_len[parser->matrix_rows];++k){
          //parser->logmatch_quality[parser->matrix_rows][k] = val.logmatch_qualities[i*MAX_MATRIX_WIDTH+k];
          //parser->logmismatch_quality[parser->matrix_rows][k] = val.logmismatch_qualities[i*MAX_MATRIX_WIDTH+k];
        }
        uint offset = parser->cursor;
        parser->offset[parser->matrix_rows] = offset;
        ++parser->matrix_rows;
        // this is the bottom boundary of the matrix, exit for loop
        if (parser->matrix_rows==parser->MAX_SUPERREADS){
          //cerr<<"Hit bottom boundary breaking\n";
          break;
        }
      }
    }
    // the offset is moved by one column for next SNP
    //++parser->cursor;
    parser->found_positions.push_back(key.querypos); 
    //if (debug) cerr<<"Cursor is now "<<parser->cursor<<endl;
  }
}

void ReadParser::add_to_pileup_maps(pileup_key_t & key,pileup_val_t & val){
  if (debug) cerr<<"Adding to pileup map for subject "<<key.subject<<" at pos "<<key.querypos<<endl;
  pileup_map[key] = val;
  for(int i=1;i<val.total_positions;++i){
    if (debug) cerr<<"At position "<<val.positions[i]<<endl;
    pileup_key_t key_partial = key;
    key_partial.querypos = val.positions[i];
    pileup_val_t val_partial;
    make_partial_pileup(val.positions[0],val.positions[i],i,key_partial,val_partial);
    if (debug)cerr<<"Adding into partial pileup map key subject "<<key_partial.subject<<" and pos "<<key_partial.querypos<<endl;
    partial_pileup_map.insert(pair<pileup_key_t,pileup_val_t>(key_partial,val_partial));
  }
}

// callback for bam_plbuf_init()
int ReadParser::pileup_func(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data)
{
  // query pos is 0 based, while SNP list is 1 based
  ReadParser *parser = (ReadParser *)data;
  tmpstruct_t *tmp = parser->tmp;
  pileup_key_t plkey;
  plkey.subject = tmp->in->header->target_name[tid];
  uint genomic_pos = querypos + 1;
  plkey.querypos = genomic_pos;
  bool debug =   plkey.subject.compare(test_subject_str)==0;
  if ((int)querypos >= tmp->beg && (int)querypos < tmp->end ){
    if(debug)cerr<<"Pileup function: Read found: querypos,beg,end: "<<querypos<<","<<tmp->beg<<","<<tmp->end<<"!\n";
    pileup_val_t plval;
    plval.total_reads = n;
    if (debug) cerr<<"Subject: "<<tmp->in->header->target_name[tid]<<" pos: "<<genomic_pos<<" total reads: "<< n<<endl;
    set<int> position_set;
    int vectorized_col = 0;
    float temp_logmatch[MAX_SUPERREADS*MAX_MATRIX_WIDTH];
    float temp_logmismatch[MAX_SUPERREADS*MAX_MATRIX_WIDTH];
    for(int i=0;i<n;++i){
      //assert(i<MAX_SUPERREADS);
      if (i==MAX_SUPERREADS) {
        cerr<<"WARNING: A total of "<<n<<" super reads were encountered at this position "<<genomic_pos<<", exceeding the max allocation of "<<MAX_SUPERREADS<<". Consider increasing the MAX_SUPERREADS value.\n";
        throw "Read parser error";
        break;
      }
      uint current_genomic_pos = genomic_pos;
      const bam1_t * bam = pl[i].b;
      char * name = bam1_qname(bam); 
      int bam_pos = bam->core.pos;
      if (bam_pos+1!=genomic_pos){
        if (debug) cerr<<" Read "<<i<<" at BAM pos "<<bam_pos+1<<" skipped as it would misalign with matrix\n";
        return 0;
      }else{
        //plval.read_names[i] = name;
        uint32_t* cigar = bam1_cigar(bam);
        int cigar_len = bam_cigar2qlen(&(bam->core),cigar);
        if (debug) cerr<<" Read "<<i<<" with name "<<name<<" at BAM pos "<<bam_pos+1<<" will be processed\n";
        uint8_t * x1 = bam_aux_get(bam,"x1");
        string mesg_x1,mesg_y1;
        if (x1){
          mesg_x1 = bam_aux2Z(x1);
          int mesg_len = mesg_x1.length();
          plval.depths[i] = mesg_len/2;
          //if (plval.depths[i]>MAX_SUBREADS){
            //cerr<<"WARNING: A read depth of "<<plval.depths[i]<<" was encountered for super read "<<i<<", exceeding the max allocation of "<<MAX_SUBREADS<<". Consider increasing the MAX_SUBREADS value.\n";
            //throw "Read parser error";
          //}
        }
        uint8_t * y1 = bam_aux_get(bam,"y1");
        if (y1){
          mesg_y1 = bam_aux2Z(y1);
          int mesg_len = mesg_y1.length();
          plval.lengths[i] = mesg_len/(2*plval.depths[i]);
          //if (plval.lengths[i]>MAX_MATRIX_WIDTH){
            //cerr<<"WARNING: A read length of "<<plval.lengths[i]<<" was encountered for super read "<<i<<", exceeding the max allocation of "<<MAX_MATRIX_WIDTH<<". Consider increasing the MAX_MATRIX_WIDTH value.\n";
            //throw "Read parser error";
          //}
          for(int k=0;k<plval.lengths[i];++k){
            if (k<MAX_MATRIX_WIDTH){
              temp_logmatch[i*MAX_MATRIX_WIDTH+k] = 0;
              temp_logmismatch[i*MAX_MATRIX_WIDTH+k] = 0;
            }
          }
          int stroffset = 0;
          for(int j=0;j<plval.depths[i];++j){
            for(int k=0;k<plval.lengths[i];++k){
              string element = mesg_y1.substr(stroffset,2);
              if (k<MAX_MATRIX_WIDTH){
                temp_logmatch[i*MAX_MATRIX_WIDTH+k] += logprob_match_lookup[hex2int(element)];
                temp_logmismatch[i*MAX_MATRIX_WIDTH+k] += logprob_mismatch_lookup[hex2int(element)];
              }
              stroffset += 2;
            }
          }
        }
        string seq_str(bam->core.l_qseq,' ');
        bool debug_seq = false;
        uint8_t * seq = bam1_seq(bam);
        if (debug_seq) cerr<<"  Query sequence"; 
        int j;
        for(j=0;j<bam->core.l_qseq;++j){
          int v = bam1_seqi(seq,j);
          seq_str[j] =  bam_nt16_rev_table[v];
          if (debug_seq)cerr<<" "<<j<<":"<<seq_str[j];
        }
        //plval.sequence_str[i] = seq_str;
        if (debug_seq) cerr<<endl;
        //uint offset = tmp->cursor;
        int cursor = 0;
        //plval.total_cigars[i] =bam->core.n_cigar; 
        current_genomic_pos = bam_pos ;
        plval.offsets[i] = vectorized_col;
        for(int k=0;k<bam->core.n_cigar;++k){
          int cop = cigar[k] & BAM_CIGAR_MASK;
          int c1 = cigar[k] >> BAM_CIGAR_SHIFT;
          if (debug) cerr<<"  At cursor: "<<cursor<<endl;
          if (debug) cerr<<"  Current cigar len "<<c1<<" of type ";
          if (cop==BAM_CMATCH ){
            //assert(cursor<MAX_MATRIX_WIDTH);
            string read_substr(seq_str,cursor,1); // read character at cursor
            if(cursor<MAX_MATRIX_WIDTH && vectorized_col<MAX_VECTOR_LENGTH){
              if (read_substr.compare("=")==0){
                if (debug) cerr<<" MATCH\n";
                //plval.alleles[i*MAX_MATRIX_WIDTH + cursor] = 0;
                plval.alleles_compact[vectorized_col] = 0;
              }else if (read_substr.compare("N")==0){
                if (debug) cerr<<" MISMATCH\n";
                //plval.alleles[i*MAX_MATRIX_WIDTH+cursor] = 1;
                plval.alleles_compact[vectorized_col] = 1;
              }else{
                cerr<<"Unknown read_substr is "<<read_substr<<" for subject "<<plkey.subject<<" and query pos "<<plkey.querypos<<endl;
                exit(0);
              }
              plval.logmatch_qualities_compact[vectorized_col] = 
              temp_logmatch[i*MAX_MATRIX_WIDTH+cursor];
              plval.logmismatch_qualities_compact[vectorized_col] = 
              temp_logmismatch[i*MAX_MATRIX_WIDTH+cursor];
              if (debug) cerr<<"  Inserted "<< plval.alleles_compact[vectorized_col]<<" into row "<<i<<" and col "<<cursor<<endl;
              ++vectorized_col;
            }
            current_genomic_pos+=c1;
            position_set.insert(current_genomic_pos);
            ++cursor;
          }else if (cop==BAM_CPAD){
            if (debug) cerr<<" CPAD\n";
            current_genomic_pos+=c1;
          }
          if (debug) cerr<<"  Genomic position is now "<<current_genomic_pos<<endl;
        }
        if (debug) cerr<<"  For read "<<i<<" offset and lengths are "<<plval.offsets[i]<<" and "<<plval.lengths[i]<<endl;
      }
      if(debug)cerr<<"  parsed cigars\n"; 
    }
    plval.total_positions = position_set.size();
    int p = 0;    
    if (debug) cerr<<" Found positions: ";
    for(set<int>::iterator it = position_set.begin();it!=position_set.end();
    it++){
      int pos = *it;
      if (debug && p) cerr<<",";
      if (debug) cerr<<*it;
      if (debug && parser->position_set.find(pos)==parser->position_set.end()){
        throw "This position is invalid";
      }
      plval.positions[p++] = *it;
    }
    if (debug) cerr<<endl;
    process_region(plkey,plval,parser);
    parser->add_to_pileup_maps(plkey, plval);
    if(debug)cerr<<"Inserted for subject "<<plkey.subject<<" and pos "<<plkey.querypos<<endl;
    // the offset is moved by one column for next SNP
    //++tmp->cursor;
    //if (debug) cerr<<"  query position is now "<<querypos<<endl;
    //parser->found_positions.push_back(querypos+1); 
    //if (debug) cerr<<"Cursor is now "<<tmp->cursor<<endl;
  }else{
    if (parser->position_set.find(plkey.querypos)!=parser->position_set.end()){
      if(debug)cerr<<"Read not found at "<<plkey.querypos<<"\n";
    }
     //exit(1);
  }
  return 0;
}

pileup_key_t::pileup_key_t(){
}

pileup_key_t::~pileup_key_t(){
}

pileup_val_t::pileup_val_t(){
  for(int i=0;i<MAX_MATRIX_WIDTH;++i){
    positions[i] = -1;
  }
}

//pileup_val_t::pileup_val_t(const pileup_val_t & old){
  // do a deep copy
//  this->total_reads = old.total_reads;
//  for(int i=0;i<this->total_reads;++i){
//    this->depths[i] = old.depths[i];
//    this->lengths[i] = old.lengths[i];
//    for(int j=0;j<this->lengths[i];++j){
//      this->alleles[i*MAX_MATRIX_WIDTH+j] = old.alleles[i*MAX_MATRIX_WIDTH+j];
//    }
//    for(int k=0;k<this->depths[i];++k){
//      for(int j=0;j<this->lengths[i];++j){
//        this->basequalities[i*MAX_SUBREADS*MAX_MATRIX_WIDTH+k*MAX_MATRIX_WIDTH+j] = old.basequalities[i*MAX_SUBREADS*MAX_MATRIX_WIDTH+k*MAX_MATRIX_WIDTH+j];
//      }
//    }
//  }
//}

pileup_val_t::~pileup_val_t(){
  //if(debug)cerr<<"Destroying pileup val of "<<total_reads<<" reads.\n";
  //delete[] read_names;
  //delete[] depths;
  //delete[] mapqual_str;
  //delete[] basequal_str;
  //delete[] sequence_str;
  //delete[] total_cigars;
  //for(int i=0;i<total_reads;++i){
  //  delete[]cigar_ops[i];
  //  delete[]cigar_len[i];
 // }
 // delete[]cigar_ops;
 // delete[]cigar_len;
}

// callback for bam_plbuf_init()
int ReadParser::pileup_func_old(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data)
{
  ReadParser *parser = (ReadParser *)data;
  tmpstruct_t *tmp = parser->tmp;
  if (debug)cerr<<"QueryPos: "<<querypos<<" "<<tmp->beg<<","<<tmp->end<<endl;
  if ((int)querypos >= tmp->beg && (int)querypos < tmp->end && tmp->cursor<parser->MAX_MATRIX_WIDTH && parser->matrix_rows<parser->MAX_SUPERREADS){
    if (debug) cerr<<"Chr: "<<tmp->in->header->target_name[tid]<<" pos: "<<querypos+1<<" reads: "<< n<<endl;
    int i=0;
    for(i=0;i<n;++i){
      uint genomic_pos = querypos + 1;
      if (debug) cerr<<" Read: "<<i<<endl;
      const bam1_t * bam = pl[i].b;
      char * name = bam1_qname(bam); int pos = bam->core.pos+1;
      uint32_t* cigar = bam1_cigar(bam);
      int cigar_len = bam_cigar2qlen(&(bam->core),cigar);
      if (debug) cerr<<"  Name "<<name<<" at pos "<<pos<<" with cigar len "<<cigar_len<<endl;
      uint8_t * x1 = bam_aux_get(bam,"x1");
      string mesg_x1,mesg_y1;
      if (x1){
        mesg_x1 = bam_aux2Z(x1);
        int mesg_len = mesg_x1.length();
        parser->depth[parser->matrix_rows] = mesg_len/2;
        for(int i=0;i<mesg_len/2;){
          string element = mesg_x1.substr(i*2,2);
          //parser->map_quality[parser->matrix_rows][i] = ReadParser::hex2int(element);
          ++i;
          if (i==MAX_SUPERREADS) break;
        }
      }
      uint8_t * y1 = bam_aux_get(bam,"y1");
      if (y1){
        mesg_y1 = bam_aux2Z(y1);
        int mesg_len = mesg_y1.length();
        parser->read_len[parser->matrix_rows] = mesg_len/(2*parser->depth[parser->matrix_rows]);
        int stroffset = 0;
        if (parser->depth[parser->matrix_rows]>MAX_SUBREADS){
          cerr<<"A read depth of "<<parser->depth[parser->matrix_rows]<<" was encountered for super read "<<parser->matrix_rows<<", exceeding the max allocation of "<<MAX_SUBREADS<<". Consider increasing the MAX_SUBREADS value.\n";
          throw "Read parser error";
        }
        for(int j=0;j<parser->read_len[parser->matrix_rows];++j){
          parser->logmatch_quality[i][j] = 0;
          parser->logmismatch_quality[i][j] = 0;
        }
        for(int k=0;k<parser->depth[parser->matrix_rows];++k){
          for(int j=0;j<parser->read_len[parser->matrix_rows];++j){
            string element = mesg_y1.substr(stroffset,2);
            if (j<parser->MAX_MATRIX_WIDTH){
              //parser->base_quality[parser->matrix_rows][i][j] = ReadParser::hex2int(element);
              parser->logmatch_quality[i][j] += logprob_match_lookup[hex2int(element)];
              parser->logmismatch_quality[i][j] += logprob_mismatch_lookup[hex2int(element)];
            }
            stroffset += 2;
          }
        }
      }
      string seq_str(bam->core.l_qseq,' ');
      bool debug_seq = false;
      uint8_t * seq = bam1_seq(bam);
      if (debug_seq) cerr<<"  Query sequence:"; 
      int j;
      for(j=0;j<bam->core.l_qseq;++j){
        int v = bam1_seqi(seq,j);
        seq_str[j] =  bam_nt16_rev_table[v];
        if (debug_seq)cerr<<" "<<j<<":"<<seq_str[j];
      }
      if (debug_seq) cerr<<endl;
      uint offset = tmp->cursor;
      uint cursor = 0;
      for(int k=0;k<bam->core.n_cigar;++k){
        int cop = cigar[k] & BAM_CIGAR_MASK;
        int c1 = cigar[k] >> BAM_CIGAR_SHIFT;
        if (debug) cerr<<"  Current cigar len "<<c1<<" of type ";
        if (cop==BAM_CMATCH ){
          string read_substr(seq_str,cursor,1);
          if (offset+cursor<MAX_MATRIX_WIDTH){
            if (read_substr.compare("=")==0){
              if (debug) cerr<<" MATCH\n";
               parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+offset+cursor] = 0;
            }else if (read_substr.compare("N")==0){
              if (debug) cerr<<" MISMATCH\n";
              parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+offset+cursor] = 1;
            }else{
              cerr<<"Unknown read_substr is "<<read_substr<<endl;
              exit(0);
            }
            //if (debug)  cerr<<"  Inserted "<< parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+offset+cursor]<<" into row "<<parser->matrix_rows<<" and col "<<offset+cursor<<endl;
          }
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
      // this is the bottom boundary of the matrix, exit for loop
      if (parser->matrix_rows==parser->MAX_SUPERREADS){
        //cerr<<"Hit bottom boundary breaking\n";
        break;
      }
    }
    // the offset is moved by one column for next SNP
    //++tmp->cursor;
    if (debug) cerr<<"  query position is now "<<querypos<<endl;
    parser->found_positions.push_back(querypos+1); 
    //if(debug)  cerr<<"Cursor is now "<<tmp->cursor<<endl;
  }else{
     //cerr<<"No reads found in this interval\n";
  }
  return 0;
}

void ReadParser::make_partial_pileup(int lastpos,int currentpos,int offset,pileup_key_t & key, pileup_val_t & val){
  bool debug = key.subject.compare(test_subject_str)==0;
  pileup_key_t orig_key;
  orig_key.subject = key.subject;
  orig_key.querypos = lastpos;
  if (debug)cerr<<"Making partial pileup for subject "<<key.subject<<" from "<<lastpos<<" to "<<currentpos<<endl;
  pileup_map_t::iterator it = pileup_map.find(orig_key);
  if (it==pileup_map.end()){
    throw "Couldn't find original pileup!";
  }else{
    //int offset = current_snpindex - lastfound_snpindex;
    pileup_val_t orig_val = it->second;
    val.total_reads = orig_val.total_reads;
//    cerr<<"Original pileup had "<<val.total_reads<<" total reads.\n";
    //for(int i=0;i<val.total_reads;++i){
      //val.depths[i] = orig_val.depths[i];
      //val.lengths[i] = orig_val.lengths[i]-offset;
      //for(int j=0;j<val.lengths[i];++j){
       //if (offset+j<MAX_MATRIX_WIDTH){
         //val.alleles[i*MAX_MATRIX_WIDTH+j] = orig_val.alleles[i*MAX_MATRIX_WIDTH+offset+j];
         //val.logmatch_qualities[i*MAX_MATRIX_WIDTH+j] = orig_val.logmatch_qualities[i*MAX_MATRIX_WIDTH+offset+j];
         //val.logmismatch_qualities[i*MAX_MATRIX_WIDTH+j] = orig_val.logmismatch_qualities[i*MAX_MATRIX_WIDTH+offset+j];
       //}
      //}
    //}
    // new vectorized matrix
    int k=0;
    for(int i=0;i<val.total_reads;++i){
      for(int j=0;j<orig_val.lengths[i];++j){
        //if(debug)cerr<<"read "<<i<<" col "<<j<<" is "<<orig_val.alleles_compact[k]<<endl;
        ++k;
      }
    }
    
    int output_vectorized_col = 0;
    for(int i=0;i<val.total_reads;++i){
      
      val.depths[i] = orig_val.depths[i];
      val.lengths[i] = orig_val.lengths[i]-offset;
      val.offsets[i] = output_vectorized_col;
      //if(debug) cerr<<"lengths,offsets at read "<<i<<": "<<val.lengths[i]<<","<<orig_val.offsets[i]<<endl;
      for(int j=0;j<orig_val.lengths[i];++j){
        if (j>=offset){
          if(output_vectorized_col<MAX_VECTOR_LENGTH){
            val.alleles_compact[output_vectorized_col] = 
            orig_val.alleles_compact[orig_val.offsets[i]+j];
            //if(debug) cerr<<"At val.offset "<<orig_val.offsets[i]<<" j "<<j<<" alleles is "<<val.alleles_compact[output_vectorized_col]<<endl;
            val.logmatch_qualities_compact[output_vectorized_col] = 
            orig_val.logmatch_qualities_compact[orig_val.offsets[i]+j];
            val.logmismatch_qualities_compact[output_vectorized_col] = 
            orig_val.logmismatch_qualities_compact[orig_val.offsets[i]+j];
          }
          ++output_vectorized_col;
        }
      }
    }
    for(int j=0;j<output_vectorized_col;++j){
      //if (debug)cerr<<"Allele at vectorized col "<<j<<" is now "<<val.alleles_compact[j]<<endl;
    }
    
  }
}

void ReadParser::extract_region(int subject,int offset,int snps,bool populate_matrix){
  bool debug = subject==test_subject;
  ostringstream oss_subject;
  oss_subject<<subject;
  string subject_str = oss_subject.str();
  matrix_rows = 0;
  tmp->cursor = cursor = 0;
  for(int i=0;i<MAX_MATRIX_WIDTH*MAX_SUPERREADS;++i){
    matrix_alleles[i] = 9;
  }
  found_positions.clear();
  //int test_subject=269;
  int lastfound_snpindex = -1;
  for(int snpindex = offset;snpindex<offset+snps;++snpindex){
    if (matrix_rows==MAX_SUPERREADS) break;  // we're done with this matrix!
    ostringstream oss_subject;
    oss_subject<<subject;
    // set up a search key for any lookups to pile up tables.
    pileup_key_t key;
    key.subject = subject_str;
    key.querypos = position_vec[snpindex];

    if (snpindex==offset && populate_matrix){
      // If at the first column, we need to extract reads that begin
      // left of the matrix
      if(subject==test_subject) cerr<<"Checking for pileups overlapping left edge\n";
      int p=0;
      for(pileup_multimap_t::iterator it = partial_pileup_map.lower_bound(key);
      it!=partial_pileup_map.upper_bound(key);it++){
       if(subject==test_subject) cerr<<"Found an overlapping pileup: "<<p<<"\n";
       pileup_val_t partial_val= it->second;
       process_region(it->first,partial_val,this);
       ++p;
      }
    }

    if(subject==test_subject)cerr<<"Querying internal pileups in cache for subject "<<key.subject<<" and snpindex "<<snpindex<<" which has pos "<<key.querypos<<endl;
    pileup_map_t::iterator it = pileup_map.find(key);
    if (it!=pileup_map.end()){
      // first attempt to fetch from memory
      pileups_found = true;
      pileup_key_t foundkey = it->first;
      if(subject==test_subject )cerr<<"Found cached pileup\n";
      pileup_val_t foundval = it->second;
      process_region(foundkey,foundval,this);
      lastfound_snpindex = snpindex;
    }else if(untyped_pileup_set.find(key)!=untyped_pileup_set.end()){
      if(debug) cerr<<"Was not found in samtools, will not attempt again\n";
    }else{
      //if memory fetch fails fetch from disk
      //if (untyped_pileup_set.find(key)==untyped_pileup_set.end()){
      string regionstr = get_region(subject,snpindex);
      const char * region = regionstr.data();
      //int ref;
      // extract the region next from BAM file
      int ref;
      if(subject==test_subject)cerr<<"Cache miss, samtool querying region: "<<region<<endl;
      bam_parse_region(tmp->in->header, region, &ref,
        &tmp->beg, &tmp->end); // parse the region
      if (ref < 0) {
        fprintf(stderr, "Invalid region %s\n", region);
        exit(1);
      }
      //bam_plbuf_t *buf = bam_plbuf_init(pileup_func, this);
      bam_plbuf_reset(buf);
      bam_fetch(tmp->in->x.bam, idx, ref, tmp->beg, tmp->end, buf, fetch_func);
      bam_plbuf_push(0, buf); // finalize pileup
      //bam_plbuf_destroy(buf);
      //if (use_cache){
      if(pileup_map.find(key)==pileup_map.end()){
        if(subject==test_subject)cerr<<"No pileups were loaded for subject "<<key.subject<<", position "<<key.querypos<<".\n";
        untyped_pileup_set.insert(key);
      }else{
        // pileups were found for snpindex
        pileups_found = true;
        if(subject==test_subject)cerr<<"Pileups successfully loaded into cache for subject "<<key.subject<<", position "<<key.querypos<<". Marking last loaded snp index as "<<snpindex<<endl;
        lastfound_snpindex = snpindex;
      }
    }
    ++cursor;
    ++tmp->cursor;
    if(subject==test_subject ) cerr<<"Matrix row,col is now "<<matrix_rows<<","<<cursor<<endl;
  }
  if (populate_matrix){
    cerr<<" clear out things from "<<last_left_marker<<" onwards to "<<offset<<"\n";
    for(int snpindex=last_left_marker;snpindex<offset;++snpindex){
      //cerr<<"In deleter of subject "<<subject<<" and snpindex "<<snpindex<<"\n";
      pileup_key_t key;
      key.subject = subject_str;
      key.querypos = position_vec[snpindex];
      if (pileup_map.find(key)!=pileup_map.end()){
        if(subject==test_subject)cerr<<"For internal pileups deleting in cache subject "<<key.subject<<" and pos "<<key.querypos<<endl;
        pileup_map.erase(pileup_map.find(key));
      }
      if (untyped_pileup_set.find(key)!=untyped_pileup_set.end()){
        if(subject==test_subject)cerr<<"For untyped pileups deleting in cache subject "<<key.subject<<" and pos "<<key.querypos<<endl;
        untyped_pileup_set.erase(untyped_pileup_set.find(key));
      }
      //pileup_map_t::iterator it = pileup_map.find(key);
      //if (it!=pileup_map.end()){
       // if(subject==test_subject)cerr<<"For full pileups deleting subject "<<key.subject<<" and pos "<<key.querypos<<endl;
        //pileup_map.erase(it);
      //}
      if (partial_pileup_map.lower_bound(key)!=partial_pileup_map.upper_bound(key)){
        if (subject==test_subject) cerr<<"For overlapping pileups deleting in cache subject "<<key.subject<<" and pos "<<key.querypos<<endl;
        partial_pileup_map.erase(partial_pileup_map.lower_bound(key),partial_pileup_map.upper_bound(key));
      }

      // sanity check
      if (snpindex==-10){
        if (subject==test_subject) cerr<<"At "<<snpindex<<" SANITY CHECK for subject "<<subject<<endl;
        for(int i=0;i<snpindex;++i){
          pileup_key_t testkey = key;
          testkey.querypos = position_vec[i];
          if (pileup_map.find(testkey)!=pileup_map.end()){
            if (subject==test_subject) cerr<<"Existing pileup found at snp, position "<<snpindex<<","<<position_vec[i]<<endl;
          }
          if (partial_pileup_map.find(testkey)!=partial_pileup_map.end()){
            if (subject==test_subject) cerr<<"Using find, Existing partial pileup found at snp, position "<<snpindex<<","<<position_vec[i]<<endl;
          }
          if (partial_pileup_map.lower_bound(testkey)!=partial_pileup_map.upper_bound(key)){
            pileup_key_t found_key = partial_pileup_map.lower_bound(testkey)->first;
            if (subject==test_subject) cerr<<"Using lowerupper bound, Existing partial pileup found at snp, position "<<i<<","<<position_vec[i]<<" key subject,querypos: "<<found_key.subject<<","<<found_key.querypos<<endl;
          }
        }
      }
    }
//    cerr<<"Data structure sizes are untyped_pileup_set: "<<untyped_pileup_set.size()<<" pileup_map: "<<pileup_map.size()<<" partial_pileup_map: "<<partial_pileup_map.size()<<endl;
  }
}

void ReadParser::finalize(int offset){
  last_left_marker = offset;
}
  
void ReadParser::parse_positions(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Cannot open variant file "<<filename<<endl;
    exit(1);
  }
  string line;
  //getline(ifs,line);
  int j=0;
  while(getline(ifs,line)){
    istringstream iss(line);
    string chr,rs;
    int pos,morgan;
    iss>>chr>>rs>>morgan>>pos;
    //cerr<<"Position "<<j<<" is "<<pos<<endl;
    position_vec.push_back(pos);
    position_set.insert(pos);
           ++j;
  }
  if (debug) cerr<<"There are "<<position_vec.size()<<" positions available\n";
  ifs.close();
}

string ReadParser::get_region(int subject,int snpindex){
  ostringstream oss;
  if (debug) cerr<<"Fetching region\n";
  oss<<subject<<":"<<position_vec[snpindex]<<"-"<<position_vec[snpindex];
  string ossstr = oss.str();
  if (debug) cerr<<"Region is "<<ossstr<<endl;
  return ossstr;
}

//ReadParser::ReadParser(const ReadParser & old){
  //cerr<<"Copy constructor being called!\n";
//  init(old.bam_filename,old.);
//}

ReadParser::ReadParser(){
}

void ReadParser::init(const char * bam_filename,const char * pos_filename){
         parse_positions(pos_filename);
  this->pileups_found = false;
  this->bam_filename = bam_filename;
  idx = bam_index_load(bam_filename); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM indexing file is not available.\n");
    exit(1);
  }
  tmp = new tmpstruct_t;
  tmp->beg = 0; tmp->end = 0x7fffffff;
  //tmp->parser = this;
  tmp->in = samopen(bam_filename, "rb", 0);
  if (tmp->in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", bam_filename);
    exit(1);
  }
  buf = bam_plbuf_init(pileup_func, this);
  offset = new int[MAX_SUPERREADS];
  depth = new int[MAX_SUPERREADS];
  read_len = new int[MAX_SUPERREADS];
  map_quality = new int * [MAX_SUPERREADS];
  logmatch_quality = new float * [MAX_SUPERREADS];
  logmismatch_quality = new float * [MAX_SUPERREADS];
  for(int i=0;i<MAX_SUPERREADS;++i){
    logmatch_quality[i] = new float[MAX_MATRIX_WIDTH];
    logmismatch_quality[i] = new float[MAX_MATRIX_WIDTH];
    map_quality[i] = new int[MAX_SUBREADS];
  }
  matrix_alleles = new int[MAX_MATRIX_WIDTH*MAX_SUPERREADS];
  last_left_marker = 0;
  return;
}

ReadParser::~ReadParser(){
  if (idx!=0){
    cerr<<"Freeing idx object\n";
    bam_index_destroy(idx);
  }
  cerr<<"Freeing buf object\n";
  bam_plbuf_destroy(buf);
  cerr<<"Closing tmp->in\n";
  samclose(tmp->in);
  cerr<<"Freeing tmp\n";
  delete tmp;
  delete[] offset;
  delete[] depth;
  delete[] read_len;
  for(int i=0;i<MAX_SUPERREADS;++i){
    delete[] logmatch_quality[i];
    delete[] logmismatch_quality[i];
    delete[] map_quality[i];
  }
  delete[] logmatch_quality;
  delete[] logmismatch_quality;
  delete[] map_quality;
  delete[] matrix_alleles;
}

int ReadParser::get_total_depth(){
  int totaldepth = 0;
  for(uint i=0;i<matrix_rows;++i){
    totaldepth+=depth[i];
  }
  return totaldepth;
}

void ReadParser::print_data(){
  cerr<<"SNP POSITIONS:\n";
  for(uint i=0;i<found_positions.size();++i){
    cerr<<found_positions[i]<<endl;
  }
  cerr<<"ROW\tHAPLOTYPE MATRIX:\n";
  for(uint i=0;i<matrix_rows;++i){
    cerr<<i<<"\t";
    for(uint j=0;j<MAX_MATRIX_WIDTH;++j){
      //if (j) cerr<<" ";
      cerr<<matrix_alleles[i*MAX_MATRIX_WIDTH+j];
    }
    cerr<<endl;
  }
  cerr<<"ROW\tOFFSET\tDEPTH\tREAD_LENGTH\tBASE_LOGPROB\n";
  for(uint i=0;i<matrix_rows;++i){
    cerr<<i<<"\t";
    cerr<<offset[i];
    cerr<<"\t";
    cerr<<depth[i];
    cerr<<"\t";
    cerr<<read_len[i];
//    cerr<<"\t";
//    for(int j=0;j<depth[i];++j){
//      if (j) cerr<<",";
  //    cerr<<map_quality[i][j];
 //   }
    int maxcol = read_len[i]<MAX_MATRIX_WIDTH?read_len[i]:MAX_MATRIX_WIDTH;
    cerr<<"\t";
    for(int k=0;k<maxcol;++k){
      if (k) cerr<<";";
      cerr<<logmatch_quality[i][k]<<","<<logmismatch_quality[i][k];
    }
    cerr<<endl;
  }
}

int main_readpen(int argc, char *argv[])
{
  try{
    if (argc < 6) {
      fprintf(stderr, "Usage: <in.bam> <variants> <subject> <snp_start> <total_snps> \n");
      return 1;
    }
    int arg = 0;
    const char * bam_filename = argv[++arg];
    const char * variants_file = argv[++arg];
    int subject = atoi(argv[++arg]);
    int offset = atoi(argv[++arg]);
    int total_regions = atoi(argv[++arg]);
  
    //ReadParser * parser = new ReadParser(bam_filename,variants_file);
    //ReadParser::parse_positions(variants_file);
    ReadParser * parser = new ReadParser();
    parser->init(bam_filename,variants_file);
    parser->extract_region(subject,offset,total_regions,false);
    parser->print_data();
    delete parser;
    return 0;
  }catch (const char * & mesg){
    cerr<<"Exception caught of message "<<mesg<<endl;
  }
}

