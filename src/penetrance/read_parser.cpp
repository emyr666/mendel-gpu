#include"../read_utilities.hpp"

bool ReadParser::debug = false;
bool ReadParser::use_cache = true;

//vector<int> ReadParser::positions;
//map<pileup_key_t,pileup_val_t,by_subject_pos> ReadParser::pileup_map;
float ReadParser::min_log_phred = log(ReadParser::EPSILON);
float ReadParser::max_log_phred = log(1.-ReadParser::EPSILON);
float ReadParser::logprob_match_lookup[100];
float ReadParser::logprob_mismatch_lookup[100];

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
    if (debug) cerr<<"Subject: "<<key.subject<<" pos: "<<genomic_pos<<" reads: "<<val.total_reads<<endl;
    for(int i=0;i<val.total_reads;++i){
      if (debug) cerr<<" Read: "<<i<<" matrix rows: "<<parser->matrix_rows<<endl;
      //string name = val.read_names[i];
      //if (debug) cerr<<"  Name "<<name<<endl;
      //string mesg_x1 = val.mapqual_str[i];
      //string mesg_y1 = val.basequal_str[i];
      //if (debug) cerr<<"X1 is "<<mesg_x1<<" and Y1 is "<<mesg_y1<<endl;
      //int mesg_len = mesg_x1.length();
      parser->depth[parser->matrix_rows] = val.depths[i];
      if(debug) cerr<<"read depth is "<<parser->depth[parser->matrix_rows]<<endl;
      //parser->depth[parser->matrix_rows] = mesg_len/2;
      //if (parser->depth[parser->matrix_rows]>MAX_SUBREADS){
       //   cerr<<"A read depth of "<<parser->depth[parser->matrix_rows]<<" was encountered for super read "<<parser->matrix_rows<<", exceeding the max allocation of "<<MAX_SUBREADS<<". Consider increasing the MAX_SUBREADS value.\n";
       //   throw "Read parser error";
      //}
      //for(int j=0;j<mesg_len/2;){
       // string element = mesg_x1.substr(j*2,2);
        //parser->map_quality[parser->matrix_rows][j] = ReadParser::hex2int(element);
        //++j;
      //}
      //mesg_len = mesg_y1.length();
      //parser->read_len[parser->matrix_rows] = mesg_len/(2*parser->depth[parser->matrix_rows]);
      parser->read_len[parser->matrix_rows] = val.lengths[i];
      if(debug) cerr<<"Read length is "<<parser->read_len[parser->matrix_rows]<<endl;
      //int stroffset = 0;
      for(int j=0;j<parser->read_len[parser->matrix_rows];++j){
        int col = parser->cursor+j;
        if(col<MAX_MATRIX_WIDTH){
          parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+col] = val.alleles[i*MAX_MATRIX_WIDTH+j];
          if (debug) cerr<<"Allele at col "<<col<<" is "<<parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+col]<<endl;
        }
      }
      for(int k=0;k<parser->read_len[parser->matrix_rows];++k){
        parser->logmatch_quality[parser->matrix_rows][k] = val.logmatch_qualities[i*MAX_MATRIX_WIDTH+k];
        parser->logmismatch_quality[parser->matrix_rows][k] = val.logmismatch_qualities[i*MAX_MATRIX_WIDTH+k];
      }
      //for(int j=0;j<parser->depth[parser->matrix_rows];++j){
        //for(int k=0;k<parser->read_len[parser->matrix_rows];++k){
          //int col = parser->cursor+k;
          //if(col<MAX_MATRIX_WIDTH){
            //parser->base_quality[parser->matrix_rows][j][k] = val.basequalities[i*MAX_SUBREADS*MAX_MATRIX_WIDTH+j*MAX_MATRIX_WIDTH+k];
            //if (debug) cerr<<"Quality at matrix row  "<<parser->matrix_rows<<" subread "<<j<<" col "<<k<<" is "<<parser->base_quality[parser->matrix_rows][j][k]<<endl;
            
          //}
          //string element = mesg_y1.substr(stroffset,2);
          //if (j<parser->MAX_MATRIX_WIDTH){
           // parser->base_quality[parser->matrix_rows][l][j] = ReadParser::hex2int(element);
          //}
          //stroffset += 2;
        //}
      //}
      //string seq_str = val.sequence_str[i];
      //bool debug_seq = false;
      //if (debug_seq) cerr<<"  Query sequence:"; 
      //int j;
      //for(j=0;j<seq_str.length();++j){
        //if (debug_seq)cerr<<" "<<j<<":"<<seq_str[j];
      //}
      //if (debug_seq) cerr<<endl;
      uint offset = parser->cursor;
      //uint cursor = 0;
      //for(int k=0;k<val.total_cigars[i];++k){
      //  int cop = val.cigar_ops[i][k];
      //  int c1 = val.cigar_len[i][k];
      //  if (debug) cerr<<"  Current cigar len "<<c1<<" of type ";
      //  assert(offset+cursor<MAX_MATRIX_WIDTH);
      //  if (cop==BAM_CMATCH ){
      //    string read_substr(seq_str,cursor,1);
      //    if (offset+cursor<MAX_MATRIX_WIDTH){
      //      if (read_substr.compare("=")==0){
      //        if (debug) cerr<<" MATCH\n";
      //         parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+offset+cursor] = 0;
      //      }else if (read_substr.compare("N")==0){
      //        if (debug) cerr<<" MISMATCH\n";
      //        parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+offset+cursor] = 1;
      //      }else{
      //        cerr<<"Unknown read_substr is "<<read_substr<<endl;
      //        exit(0);
      //      }
      //      if (debug)  cerr<<"  Inserted "<< parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+offset+cursor]<<" into row "<<parser->matrix_rows<<" and col "<<offset+cursor<<endl;
      //    }
      //  }else if (cop==BAM_CPAD){
      //    if (debug) cerr<<"  CPAD\n";
      //    ++cursor;
      //  }
     // }
      parser->offset[parser->matrix_rows] = offset;
      ++parser->matrix_rows;
      // this is the bottom boundary of the matrix, exit for loop
      if (parser->matrix_rows==parser->MAX_SUPERREADS){
        //cerr<<"Hit bottom boundary breaking\n";
        break;
      }
    }
    // the offset is moved by one column for next SNP
    //++parser->cursor;
    parser->found_positions.push_back(key.querypos); 
    //if (debug) cerr<<"Cursor is now "<<parser->cursor<<endl;
  }
}

// callback for bam_plbuf_init()
int ReadParser::pileup_func(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data)
{
  ReadParser *parser = (ReadParser *)data;
  tmpstruct_t *tmp = parser->tmp;
  pileup_key_t plkey;
  plkey.subject = tmp->in->header->target_name[tid];
  plkey.querypos = querypos+1; 
  if(debug)cerr<<"Subject: "<<tmp->in->header->target_name[tid]<<endl;
  if(debug)cerr<<"QueryPos: "<<querypos+1<<" "<<tmp->beg<<","<<tmp->end<<endl;
  if ((int)querypos >= tmp->beg && (int)querypos < tmp->end ){
    if(debug)cerr<<"Read found at "<<plkey.querypos<<"\n";
    pileup_val_t plval;
    plval.total_reads = n;
    if (debug) cerr<<"Chr: "<<tmp->in->header->target_name[tid]<<" pos: "<<querypos+1<<" reads: "<< n<<endl;
    if (n>MAX_SUPERREADS){
      cerr<<"A total of "<<n<<" super reads were encountered at this position "<<querypos+1<<", exceeding the max allocation of "<<MAX_SUPERREADS<<". Consider increasing the MAX_SUPERREADS value.\n";
      throw "Read parser error";
     
    }
    for(int i=0;i<n;++i){
      uint genomic_pos = querypos + 1;
      if (debug) cerr<<" Read: "<<i<<endl;
      const bam1_t * bam = pl[i].b;
      char * name = bam1_qname(bam); int pos = bam->core.pos+1;
      //plval.read_names[i] = name;
      uint32_t* cigar = bam1_cigar(bam);
      int cigar_len = bam_cigar2qlen(&(bam->core),cigar);
      if (debug) cerr<<"  Name "<<name<<" at pos "<<pos<<" with cigar len "<<cigar_len<<endl;
      uint8_t * x1 = bam_aux_get(bam,"x1");
      string mesg_x1,mesg_y1;
      if (x1){
        mesg_x1 = bam_aux2Z(x1);
        int mesg_len = mesg_x1.length();
        plval.depths[i] = mesg_len/2;
        if (plval.depths[i]>MAX_SUBREADS){
//        if (parser->depth[parser->matrix_rows]>MAX_SUBREADS){
          cerr<<"A read depth of "<<plval.depths[i]<<" was encountered for super read "<<i<<", exceeding the max allocation of "<<MAX_SUBREADS<<". Consider increasing the MAX_SUBREADS value.\n";
          throw "Read parser error";
        }
        //parser->depth[parser->matrix_rows] = mesg_len/2;
        //for(int i=0;i<mesg_len/2;){
          //string element = mesg_x1.substr(i*2,2);
          //parser->map_quality[parser->matrix_rows][i] = ReadParser::hex2int(element);
          //++i;
          //if (i==MAX_SUPERREADS) break;
        //}
        //plval.mapqual_str[i] = mesg_x1;
      }
      uint8_t * y1 = bam_aux_get(bam,"y1");
      if (y1){
        mesg_y1 = bam_aux2Z(y1);
        //plval.basequal_str[i] = mesg_y1;
        int mesg_len = mesg_y1.length();
        plval.lengths[i] = mesg_len/(2*plval.depths[i]);
        if (plval.lengths[i]>MAX_MATRIX_WIDTH){
//        if (parser->depth[parser->matrix_rows]>MAX_SUBREADS){
          cerr<<"A read length of "<<plval.lengths[i]<<" was encountered for super read "<<i<<", exceeding the max allocation of "<<MAX_MATRIX_WIDTH<<". Consider increasing the MAX_MATRIX_WIDTH value.\n";
          throw "Read parser error";
        }
        ////parser->read_len[parser->matrix_rows] = mesg_len/(2*parser->depth[parser->matrix_rows]);
        for(int k=0;k<plval.lengths[i];++k){
          plval.logmatch_qualities[i*MAX_MATRIX_WIDTH+k] = 0;
          plval.logmismatch_qualities[i*MAX_MATRIX_WIDTH+k] = 0;
        }
        int stroffset = 0;
        for(int j=0;j<plval.depths[i];++j){
        //for(int i=0;i<parser->depth[parser->matrix_rows];++i){
          for(int k=0;k<plval.lengths[i];++k){
        //  for(int j=0;j<parser->read_len[parser->matrix_rows];++j){
            string element = mesg_y1.substr(stroffset,2);
        //    if (j<parser->MAX_MATRIX_WIDTH){
              //plval.basequalities[i*MAX_SUBREADS*MAX_MATRIX_WIDTH+j*MAX_MATRIX_WIDTH+k] = ReadParser::hex2int(element);
              //if (debug) cerr<<"Setting base quality at subread "<<j<<" and col "<<k<<" to "<<plval.basequalities[i*MAX_SUBREADS*MAX_MATRIX_WIDTH+j*MAX_MATRIX_WIDTH+k]<<endl;
              plval.logmatch_qualities[i*MAX_MATRIX_WIDTH+k] += logprob_match_lookup[hex2int(element)];
              plval.logmismatch_qualities[i*MAX_MATRIX_WIDTH+k] += logprob_mismatch_lookup[hex2int(element)];
        //      parser->base_quality[parser->matrix_rows][i][j] = ReadParser::hex2int(element);
        //    }
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
      uint cursor = 0;
      //plval.total_cigars[i] =bam->core.n_cigar; 
      for(int k=0;k<bam->core.n_cigar;++k){
        int cop = cigar[k] & BAM_CIGAR_MASK;
        int c1 = cigar[k] >> BAM_CIGAR_SHIFT;
        //plval.cigar_ops[i][k] = cop;
        //plval.cigar_len[i][k] = c1;
        if (debug) cerr<<"  Current cigar len "<<c1<<" of type ";
        if (cop==BAM_CMATCH ){
          string read_substr(seq_str,cursor,1);
          //if (offset+cursor<MAX_MATRIX_WIDTH){
            if (read_substr.compare("=")==0){
              if (debug) cerr<<" MATCH\n";
                plval.alleles[i*MAX_MATRIX_WIDTH + cursor] = 0;
        //       parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+offset+cursor] = 0;
            }else if (read_substr.compare("N")==0){
              if (debug) cerr<<" MISMATCH\n";
                plval.alleles[i*MAX_MATRIX_WIDTH+cursor] = 1;
        //      parser->matrix_alleles[parser->matrix_rows * MAX_MATRIX_WIDTH+offset+cursor] = 1;
            }else{
              cerr<<"Unknown read_substr is "<<read_substr<<endl;
              exit(0);
            }
            if (debug)  cerr<<"  Inserted "<< plval.alleles[i*MAX_MATRIX_WIDTH+cursor]<<" into row "<<i<<" and col "<<cursor<<endl;
          //}
        //  genomic_pos+=c1;
        }else if (cop==BAM_CPAD){
          if (debug) cerr<<"  CPAD\n";
        //  genomic_pos+=c1;
          ++cursor;
        }
       // if (debug) cerr<<"  Genomic position is now "<<genomic_pos<<endl;
      }
      if(debug)cerr<<"parsed cigars\n"; 
      //parser->offset[parser->matrix_rows] = offset;
      //++parser->matrix_rows;
      // this is the bottom boundary of the matrix, exit for loop
      //if (parser->matrix_rows==parser->MAX_SUPERREADS){
        //cerr<<"Hit bottom boundary breaking\n";
        //break;
      //}
    }
    process_region(plkey,plval,parser);
    parser->pileup_map[plkey] = plval;
    if(debug)cerr<<"Inserted for subject "<<plkey.subject<<" and pos "<<plkey.querypos<<endl;
    // the offset is moved by one column for next SNP
    //++tmp->cursor;
    //if (debug) cerr<<"  query position is now "<<querypos<<endl;
    //parser->found_positions.push_back(querypos+1); 
    //if (debug) cerr<<"Cursor is now "<<tmp->cursor<<endl;
  }else{
    if (parser->position_set.find(plkey.querypos)!=parser->position_set.end()){
      //parser->untyped_pileup_set.insert(plkey);
      if(debug)cerr<<"Read not found at "<<plkey.querypos<<"\n";
    }
     //exit(1);
  }
  return 0;
}

pileup_key_t::pileup_key_t(){
}

//pileup_key_t::pileup_key_t(const pileup_key_t & old){
 // this->subject = old.subject;
 // this->querypos = old.querypos;
//}

pileup_key_t::~pileup_key_t(){
}

pileup_val_t::pileup_val_t(){
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


void ReadParser::extract_region(int subject,int offset,int snps){
  if (debug) cerr<<"Selected subject "<<subject<<endl;
  ostringstream oss_subject;
  oss_subject<<subject;
  string subject_str = oss_subject.str();
  matrix_rows = 0;
  tmp->cursor = cursor = 0;
  for(int i=0;i<MAX_MATRIX_WIDTH*MAX_SUPERREADS;++i){
    matrix_alleles[i] = 9;
  }
  found_positions.clear();
  int test_subject=-10;
  for(int snpindex = offset;snpindex<offset+snps;++snpindex){
    if (matrix_rows==MAX_SUPERREADS) break;  // we're done with this matrix!
    ostringstream oss_subject;
    oss_subject<<subject;
    pileup_key_t key;
    key.subject = subject_str;
    key.querypos = position_vec[snpindex];
    if(subject==test_subject)cerr<<"Querying in CACHE subject "<<key.subject<<" and snpindex "<<snpindex<<" which is pos "<<key.querypos<<endl;
    pileup_map_t::iterator it = pileup_map.find(key);
    if (use_cache && it!=pileup_map.end()){
      pileup_key_t foundkey = it->first;
      if(subject==test_subject )cerr<<"FOUND CACHED old sub pos "<<key.subject<<","<<key.querypos<<" and new "<<foundkey.subject<<","<<foundkey.querypos<<endl;
      //if(debug && subject==test_subject )cerr<<"FOUND CACHED old sub pos "<<key.subject<<","<<key.querypos<<" and new "<<foundkey.subject<<","<<foundkey.querypos<<endl;
      pileup_val_t foundval = it->second;
      process_region(foundkey,foundval,this);
    }else{
      if (untyped_pileup_set.find(key)==untyped_pileup_set.end()){
        string regionstr = get_region(subject,snpindex);
        const char * region = regionstr.data();
        //const char * region = "20:60749-60749";
        //int ref;
        // extract the region next from BAM file
        int ref;
        if(debug)cerr<<"Querying SNP position:  "<<region<<endl;
        bam_parse_region(tmp->in->header, region, &ref,
          &tmp->beg, &tmp->end); // parse the region
        if (ref < 0) {
  		fprintf(stderr, "Invalid region %s\n", region);
  		exit(1);
        }
        bam_plbuf_t *buf = use_cache?bam_plbuf_init(pileup_func, this):bam_plbuf_init(pileup_func_old, this); // initialize pileup
        bam_fetch(tmp->in->x.bam, idx, ref, tmp->beg, tmp->end, buf, fetch_func);
        bam_plbuf_push(0, buf); // finalize pileup
        bam_plbuf_destroy(buf);
        if (use_cache && pileup_map.find(key)==pileup_map.end()){
          if(subject==test_subject)cerr<<"Would add in "<<key.subject<<","<<key.querypos<<endl;
          untyped_pileup_set.insert(key);
        }
      }else{
        if(subject==test_subject )cerr<<"SKIPPING CACHED old sub pos "<<key.subject<<","<<key.querypos<<endl;
      }
    }
    ++cursor;
    ++tmp->cursor;
    if(subject==test_subject ) cerr<<"Matrix row,col is now "<<matrix_rows<<","<<cursor<<endl;
  }
  //cerr<<" clear out things from "<<last_left_marker<<" onwards\n";
  for(int snpindex=last_left_marker;snpindex<offset;++snpindex){
    //cerr<<"In deleter of subject "<<subject<<" and snpindex "<<snpindex<<"\n";
    pileup_key_t key;
    key.subject = subject_str;
    key.querypos = position_vec[snpindex];
    pileup_map_t::iterator it = pileup_map.find(key);
    if (it!=pileup_map.end()){
      if(subject==test_subject)cerr<<"Deleting subject "<<key.subject<<" and pos "<<key.querypos<<endl;
      pileup_map.erase(it);
    }
    pileup_set_t::iterator it2 = untyped_pileup_set.find(key);
    if (it2!=untyped_pileup_set.end()){
      if(subject==test_subject)cerr<<"Deleting untyped subject "<<key.subject<<" and pos "<<key.querypos<<endl;
      untyped_pileup_set.erase(it2);
    }
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
  offset = new int[MAX_SUPERREADS];
  depth = new int[MAX_SUPERREADS];
  read_len = new int[MAX_SUPERREADS];
  map_quality = new int * [MAX_SUPERREADS];
  //base_quality = new int ** [MAX_SUPERREADS];
  logmatch_quality = new float * [MAX_SUPERREADS];
  logmismatch_quality = new float * [MAX_SUPERREADS];
  for(int i=0;i<MAX_SUPERREADS;++i){
    logmatch_quality[i] = new float[MAX_MATRIX_WIDTH];
    logmismatch_quality[i] = new float[MAX_MATRIX_WIDTH];
    map_quality[i] = new int[MAX_SUBREADS];
    //base_quality[i] = new int * [MAX_SUBREADS];
    //for(int j=0;j<MAX_SUBREADS;++j){
      //base_quality[i][j] = new int[MAX_MATRIX_WIDTH];
    //}
  }
  matrix_alleles = new int[MAX_MATRIX_WIDTH*MAX_SUPERREADS];
  return;
}

ReadParser::~ReadParser(){
  bam_index_destroy(idx);
  samclose(tmp->in);
  delete tmp;
  //cerr<<"Destructed\n";
}

void ReadParser::print_data(){
  cerr<<"SNP POSITIONS:\n";
  for(uint i=0;i<found_positions.size();++i){
    cerr<<found_positions[i]<<endl;
  }
  cerr<<"HAPLOTYPE MATRIX:\n";
  for(uint i=0;i<matrix_rows;++i){
    for(uint j=0;j<MAX_MATRIX_WIDTH;++j){
      //if (j) cerr<<" ";
      cerr<<matrix_alleles[i*MAX_MATRIX_WIDTH+j];
    }
    cerr<<endl;
  }
  cerr<<"OFFSET\tDEPTH\tREAD_LENGTH\tBASE_LOGPROB\n";
  for(uint i=0;i<matrix_rows;++i){
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
    parser->extract_region(subject,offset,total_regions);
      //const char * region = "20:60522-60522";
    parser->print_data();
    delete parser;
    return 0;
  }catch (const char * & mesg){
    cerr<<"Exception caught of message "<<mesg<<endl;
  }
}

