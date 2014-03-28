#include<iostream>
#include<istream>
#include<ostream>
#include<iterator>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<vector>
#include<map>
#include<set>
#include<cstdlib>
#include<stdio.h>
#include"faidx.h"
#include "sam.h"
#include"math.h"

using namespace std;

bam_index_t *idx;

int subject;
const uint MAX_READ_LEN = 1000;
const uint SNP_POINT = 0;
const uint SNP_INSERTION = 1;
const uint SNP_DELETION = 2;
int srb_read_index;

void print_output(faidx_t * fasta_index,uint current_position);

struct snp_t{
  uint index;
  uint position;
  uint pileup_pos;
  uint snp_type; //0 for regular, 1 for insertion, 2 for deletion
  string ref_allele;
  string alt_allele;
  string read_allele;
  uint8_t qual;
  bool is_poly;
};

struct by_snp_position{
  bool operator()(const snp_t & snp1, const snp_t & snp2){
    bool returnval = snp1.index<snp2.index;
    //bool returnval = snp1.position<snp2.position;
    return returnval;
  }
};

typedef set<snp_t,by_snp_position> snp_set_t;
//typedef vector<snp_t> snp_vector_t;

struct tmpstruct_t{
  int beg, end;
  samfile_t *in;
  faidx_t * fasta_index;
  snp_t snp;
};


struct reduced_bam_t{
  // mandatory BAM fields
  string qname,rname;
  //string cigar,rnext,seq,qual;
  int pos,  mapq,flag;
  //int tlen;
  // custom fields
  //int cigar_len;
  //snp_vector_t snp_vector;
};

struct small_reduced_bam_t{
  string rname;
  int pos;
  string short_cigar;
  string short_seq;
  vector<int> mapq_vector;
  vector<int> qual_vector;
};

struct by_readname{
  bool operator()(const reduced_bam_t & bam1, const reduced_bam_t & bam2){
    bool returnval = bam1.qname<bam2.qname;
    //bool returnval = bam1.qname.compare(bam2.qname)==0?bam1.flag<bam2.flag:bam1.qname<bam2.qname;
    return returnval;
  }
};

struct by_shortseq{
  bool operator()(const small_reduced_bam_t & bam1, const small_reduced_bam_t & bam2){
    //cerr<<"Debug: "<<bam1.rname<<","<<bam2.rname<<";"<<bam1.pos<<","<<bam2.pos<<";"<<bam1.short_cigar<<","<<bam2.short_cigar<<";"<<bam1.short_seq<<","<<bam2.short_seq<<endl;
    if (bam1.rname.compare(bam2.rname)!=0){
      //cerr<<"condition A\n";
      return bam1.rname<bam1.rname;
    }
    else if(bam1.pos!=bam2.pos){
      //cerr<<"condition B\n";
      return bam1.pos<bam2.pos;
    }
    else if (bam1.short_cigar.compare(bam2.short_cigar)!=0){
      //cerr<<"condition C\n";
      return bam1.short_cigar<bam2.short_cigar;
    }
    else if (bam1.short_seq.compare(bam2.short_seq)!=0){
      //cerr<<"condition D\n";
      return bam1.short_seq<bam2.short_seq;
    }
    else return false;
  }
};

typedef map<reduced_bam_t,snp_set_t, by_readname> bam_map_t;
bam_map_t bam_map;
typedef vector<reduced_bam_t> bam_vector_t;
bam_vector_t bam_vector;
typedef set<small_reduced_bam_t, by_shortseq> small_reduced_bam_set_t;
small_reduced_bam_set_t small_reduced_bam_set;


static inline float phred2prob(int phred){
  return 1-(pow(10,(-phred/10)));
}

inline void int2hex(int input, ostream & output){
  output.fill('0');
  output.width (2);
  output.flags( ios::uppercase);
  output<<hex<<input;
  output.unsetf(ios::uppercase|ios::hex);
}

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
  bam_plbuf_t *buf = (bam_plbuf_t*)data;
  bam_plbuf_push(b, buf);
  return 0;
}

bool variant_in_cigar_match(const bam1_t * & bam,tmpstruct_t *tmp, string seq_str, bool & is_poly){
  uint genomic_pos = bam->core.pos+1;
  uint32_t* cigar = bam1_cigar(bam);
  uint8_t * qual = bam1_qual(bam);
  //cerr<<"debug qual:";
  //for(int i=0;i<100;++i) cerr<<(char)(qual[i]+33);
  //cerr<<endl;
  //int cigar_leng = bam_cigar2qlen(&(bam->core),cigar);
  //cerr<<" Total length of cigars: "<<cigar_leng<<endl;
  //int last_genomic_pos = genomic_pos+cigar_leng-1;
  //cerr<<" Range for genomic position: "<<genomic_pos<<"-"<<last_genomic_pos<<endl;
  bool snp_match = false;
  uint cursor = 0;
  is_poly = false;
  cerr<<"  Number of cigars: "<<bam->core.n_cigar<<endl; 
  for(int k=0;k<bam->core.n_cigar;++k){
    int cop = cigar[k] & BAM_CIGAR_MASK;
    int c1 = cigar[k] >> BAM_CIGAR_SHIFT;
    cerr<<"  Current cigar len "<<c1<<" of type ";
    uint start,last;
    if (cop==BAM_CMATCH || cop==BAM_CEQUAL || cop==BAM_CDIFF){
    //switch (cop) {
     //case BAM_CMATCH:
       cerr<<" M";
       start = genomic_pos;
       last = start+c1-1;
       cerr<<"["<<start<<"-"<<last<<"] ";
       for(int i=0;i<c1;++i){
         //cerr<<"cursor: "<<cursor<<" pileup pos "<<tmp->snp.pileup_pos<<endl;
         if (cursor==tmp->snp.pileup_pos){
           snp_match = true;
           tmp->snp.qual = qual[tmp->snp.pileup_pos];
           //char cqual = static_cast<char>(tmp->snp.qual);
           //cerr<<"quality at "<<tmp->snp.pileup_pos<<" is "<<cqual<<endl;
           if (tmp->snp.snp_type==SNP_POINT){
             string read_substr(seq_str,cursor,1);
             cerr<<"Read variant at "<<genomic_pos<<" for point mutation is "<<read_substr<<endl;
             if (read_substr.compare(tmp->snp.ref_allele)==0)
               is_poly = false;
             else if (read_substr.compare(tmp->snp.alt_allele)==0)
               is_poly = true;
             else {
               cerr<<" Read allele found was "<<read_substr<<", mismatch between read variant and alleles for point mutation.\n";
               is_poly = true;
             }
             return snp_match; 
           }else if (tmp->snp.snp_type==SNP_INSERTION){
             if (genomic_pos==last){
               cerr<<"Found a putative insertion!\n";
               if (k < (bam->core.n_cigar-1)){
                 int next_cop = cigar[k+1] & BAM_CIGAR_MASK;
                 if (next_cop==BAM_CINS){
                   cerr<<"Confirmed insertion!\n";
                   is_poly = true;
                   return snp_match;
                 }
               }
             }
           }else if (tmp->snp.snp_type==SNP_DELETION){
             if (last<genomic_pos+tmp->snp.ref_allele.length()-1){
               cerr<<"Found a putative deletion!\n";
               if (k < (bam->core.n_cigar-1)){
                 int next_cop = cigar[k+1] & BAM_CIGAR_MASK;
                 if (next_cop==BAM_CDEL){
                   cerr<<"Confirmed deletion!\n";
                   is_poly = true;
                   return snp_match;
                 }
               }
             }
             //if (genomic_pos+tmp->snp.ref_allele.length()-1>last && genomic_pos+tmp->snp.ref_allele.length()-1<last_genomic_pos){
              // cerr<<"Found a putative deletion!\n";
               //is_poly = true;
             //}else{
               //cerr<<"Found a putative wild type no deletion!\n";
               //is_poly = false;
             //}
           }
         }
         ++cursor;
         ++genomic_pos;
       }
    //   break;
    }else if (cop==BAM_CHARD_CLIP){
//     case BAM_CHARD_CLIP:
          cerr<<" CHARD_CLIP CHECK";
          /* printf("[%d]", pos);  // No coverage */
          /* pos is not advanced by this operation */
 //         break;
    
    }else if(cop==BAM_CSOFT_CLIP){
    //   case BAM_CSOFT_CLIP:
          cerr<<" CSOFT_CLIP";
          /* printf("[%d]", pos);  // No coverage */
          /* pos is not advanced by this operation */
          cursor+=c1;
     //     break;
    
    }else if (cop== BAM_CDEL){
     //case BAM_CDEL:
          cerr<<" CDEL";
          start = genomic_pos;
          last = start+c1-1;
          cerr<<"["<<start<<"-"<<last<<"] ";
          /* printf("[%d-%d]", pos, pos + cl - 1);  // Spans positions, No Coverage */
          genomic_pos+=c1;
      //    break;
    }else if (cop==BAM_CPAD){ 
     //case BAM_CPAD:
          cerr<<" CPAD CHECK";
          /* printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage */
          genomic_pos+=c1;
      //    break;
    
    }else if (cop==BAM_CINS){
     //case BAM_CINS:
          cerr<<" CINS";
          /* printf("[%d]", pos); // Special case - adds <cl> bp "throughput", but not genomic position "coverage" */
          start = genomic_pos;
          cerr<<"["<<start<<"] ";
          /* How you handle this is application dependent */
          /* pos is not advanced by this operation */
          cursor+=c1;
      //    break;
    }else if (cop==BAM_CREF_SKIP){
     //case BAM_CREF_SKIP:
          cerr<<" CREF_SKIP CHECK";
          /* printf("[%d-%d]", pos, pos + cl - 1);  Spans positions, No Coverage */
          genomic_pos+=c1;
      //    break;
    }else{
     //default:
          cerr<< "Unhandled cigar_op "<<(char)cop<<","<<cop<<": "<< c1<<endl;
    }
    cerr<<endl;
  }
  cerr<<"  Pileup match: "<<snp_match<<" Is polymorphism: "<<is_poly<<endl;
  return snp_match;
}

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t querypos, int n, const bam_pileup1_t *pl, void *data)
{
  tmpstruct_t *tmp = (tmpstruct_t*)data;
  //cerr<<"QueryPos: "<<querypos<<" "<<tmp->beg<<","<<tmp->end<<endl;
  if ((int)querypos >= tmp->beg && (int)querypos < tmp->end){
    cerr<<"Chr "<<tmp->in->header->target_name[tid]<<" pos "<<querypos + 1<<" reads "<< n<<" with alleles: "<<tmp->snp.ref_allele<<","<<tmp->snp.alt_allele<<" snp type: "<<tmp->snp.snp_type<<endl;
    int i=0;
    for(i=0;i<n;++i){
      cerr<<" Read: "<<i<<endl;
      //printf("  indel status: %d\n", pl[i].indel); 
      //printf("  level: %d\n", pl[i].level); 
      //printf("  is_del: %d\n", pl[i].is_del); 
      //printf("  is_head: %d\n", pl[i].is_head); 
      //printf("  is_tail: %d\n", pl[i].is_tail); 
      const bam1_t * bam = pl[i].b;
      char * name = bam1_qname(bam); int pos = bam->core.pos+1;
      reduced_bam_t newbam;
      newbam.qname = name;
      newbam.flag = bam->core.flag;
      newbam.rname = tmp->in->header->target_name[tid];
      newbam.pos = pos;
      newbam.mapq = bam->core.qual;
      //uint32_t* cigar = bam1_cigar(bam);
      //int cigar_len = bam_cigar2qlen(&(bam->core),cigar);
      tmp->snp.pileup_pos = pl[i].qpos;
      //int pileup_pos = pos+pl[i].qpos;
      cerr<<"  Name "<<name<<" at pos "<<pos<<" with flag "<<bam->core.flag<<" pileup base pos "<<tmp->snp.pileup_pos<<endl; 
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
      bool is_poly;
      if (variant_in_cigar_match(bam,tmp, seq_str,is_poly)){
        snp_set_t snp_set;
        if (bam_map.find(newbam)!=bam_map.end()){
          snp_set = bam_map[newbam];
        } 
        snp_t snp;
        snp.index = tmp->snp.index;
        snp.position = tmp->snp.position;
        snp.pileup_pos = pl[i].qpos;
        snp.ref_allele = tmp->snp.ref_allele;
        snp.alt_allele = tmp->snp.alt_allele;
        snp.read_allele = seq_str[pl[i].qpos];
        snp.is_poly = is_poly;
        snp.qual = tmp->snp.qual;
        snp_set.insert(snp);
        bam_map[newbam] = snp_set;
      }else{
        cerr<<"Variant not in cigar match for this read.\n";
      }
      //free(qseq);
      //printf("  l_aux,data_len,m_data: %d,%d,%d\n", bam->l_aux,bam->data_len,bam->m_data); 
      //printf("  mapping qual: %f\n", phred2prob(bam->core.qual)); 
      //printf("  l_qseq: %d\n", bam->core.l_qseq); 
      //char * qqual = (char*)malloc(bam->core.l_qseq+1);
      //uint8_t * qual = bam1_qual(bam);
      //int j;
      //printf("  queryqual:");
      //for(j=0;j<bam->core.l_qseq;++j){
        //char v = bam1_seqi(seq,j);
        //qseq[j] =  bam_nt16_rev_table[v];
        //printf(" %f", phred2prob(qual[j])); 
      //}
      //printf("\n");
      //qseq[j] = 0;
      
    }
  }else{
     //cerr<<"No reads found in this interval\n";
  }
  return 0;
}

void parse_variants(const char * filename, tmpstruct_t & tmp){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Cannot open "<<filename<<endl;
    exit(1);
  }
  string line;
  getline(ifs,line);
  srb_read_index = 0;
  tmp.snp.index = 0;
  while(getline(ifs,line)){
    string chr,id;
    istringstream iss(line);
    iss>>chr>>tmp.snp.position>>id>>tmp.snp.ref_allele>>tmp.snp.alt_allele;
    if (tmp.snp.ref_allele.length()==tmp.snp.alt_allele.length()){
      tmp.snp.snp_type = SNP_POINT;
    }else if (tmp.snp.ref_allele.length()<tmp.snp.alt_allele.length()){
      tmp.snp.snp_type = SNP_INSERTION; 
    }else if (tmp.snp.ref_allele.length()>tmp.snp.alt_allele.length()){
      tmp.snp.snp_type = SNP_DELETION; 
    }
    // generate pileup for this region
    ostringstream oss;
    bool use_subject = subject>=0;
    if (use_subject){
      oss<<subject<<":"<<tmp.snp.position<<"-"<<tmp.snp.position;
    }else{
      oss<<chr<<":"<<tmp.snp.position<<"-"<<tmp.snp.position;
    }
    oss.flush();
    //const char * region = oss.str().data();
    string regionstr = oss.str();
    const char * region = regionstr.data();
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
    print_output(tmp.fasta_index,tmp.snp.position);
    ++tmp.snp.index;
  }
  print_output(tmp.fasta_index,1e9);
  ifs.close();
  //exit(0);
}

void print_output(faidx_t * fasta_index,uint current_position){
  int read_index = 0;
  //for (bam_vector_t::iterator it = bam_vector.begin();it!=bam_vector.end();it++){
  bool debug_rawreads = false;
  ofstream ofs_debug;
  if (debug_rawreads) ofs_debug.open("debug");
  vector<reduced_bam_t> reduced_bam_delete_vector;
  for (bam_map_t::iterator it = bam_map.begin();it!=bam_map.end();it++){
    reduced_bam_t reduced_bam = it->first;
    //cerr<<"Bam pos is "<<reduced_bam.pos<<" and current pos is "<<current_position<<endl;
    if (reduced_bam.pos+MAX_READ_LEN<current_position){
      //reduced_bam_t reduced_bam = *it;
      ostringstream oss_new_cigar,oss_new_seq;
      vector<int> base_qual_vec;
      cerr<<reduced_bam.qname<<":"<<reduced_bam.flag<<endl;
      snp_set_t snp_set = it->second;
      //snp_set_t snp_set = reduced_bam.snp_set;
      int first_position = 0;
      int last_pos = -1;
      int last_index = -1;
      int read_snp_index = 0;
      bool valid_srb = true;
      for (snp_set_t::iterator it2 = snp_set.begin();it2!=snp_set.end();it2++){
        snp_t snp = *it2;
        if (last_index>-1 && snp.index-last_index > 1){
          valid_srb = false;
          //cerr<<"This SRB should be invalid.\n";
        }
        if (read_snp_index==0) first_position = snp.position;
        if (last_pos>=0 && (snp.position-(last_pos+1)))
          oss_new_cigar<<(snp.position-(last_pos+1))<<"P";
        oss_new_cigar<<"1M";      
        char seqchar = snp.is_poly ? 'X' : '=';
        oss_new_seq<<seqchar;
        //if (read_snp_index) oss_new_qual<<",";
        //oss_new_qual<<static_cast<int>(snp.qual);
        base_qual_vec.push_back(snp.qual);
        //oss_new_qual<<static_cast<char>(snp.qual+33);
        cerr<<" snp at "<<snp.position<<" of index "<<snp.index<<" has alleles "<<snp.ref_allele<<" and "<<snp.alt_allele<<" and read allele "<<snp.read_allele<<" poly status: "<<snp.is_poly<<endl;
        int fasta_len;
        ostringstream oss_region;
        oss_region<<reduced_bam.rname<<":"<<snp.position<<"-"<<snp.position<<endl;
        string oss_region_str = oss_region.str();
        bool debug_refseq = true;
        if (debug_refseq){
          const char * fasta_rawseq = fai_fetch(fasta_index,oss_region_str.data(),&fasta_len);
          cerr<<" ref seq reports:";
          for(int i=0;i<fasta_len;++i){
            cerr<<fasta_rawseq[i];
          }
        cerr<<endl;
        }
        last_pos = snp.position;
        last_index = snp.index;
        ++read_snp_index;
      }
      //if (reduced_bam.cigar_len-(last_pileup_pos+1))
       //oss_new_cigar<<(reduced_bam.cigar_len-(last_pileup_pos+1))<<"P";
      // final output
      if (debug_rawreads) ofs_debug<<read_index<<"\t"<<reduced_bam.flag<<"\t"<<reduced_bam.rname<<
      "\t"<<first_position<<"\t"<<reduced_bam.mapq<<"\t"<<oss_new_cigar.str()<<
      "\t*\t0\t0\t"<<oss_new_seq.str()<<"\t";
       
      cerr<<endl;
      // coalesce same reads into a super read
      small_reduced_bam_t srb;
      srb.rname = reduced_bam.rname;
      srb.pos = first_position;
      srb.short_cigar = oss_new_cigar.str();
      string new_seq = oss_new_seq.str();
      srb.short_seq = new_seq;
      small_reduced_bam_set_t::iterator srb_it = small_reduced_bam_set.find(srb);
      if (srb_it!=small_reduced_bam_set.end()){
        cerr<<"Found srb!\n";
        srb = *srb_it;
        small_reduced_bam_set.erase(srb_it);
      }
      srb.mapq_vector.push_back(reduced_bam.mapq);
      for(vector<int>::iterator it = base_qual_vec.begin();it!=base_qual_vec.end();it++){
        srb.qual_vector.push_back(*it);
      }
      
      if (valid_srb) small_reduced_bam_set.insert(srb);
      ++read_index;
      reduced_bam_delete_vector.push_back(reduced_bam);
    }
  }
  if (debug_rawreads) ofs_debug.close();
  //cerr<<"SRB set is "<<small_reduced_bam_set.size()<<endl;
  vector<small_reduced_bam_t> srb_delete_vector;
  for(small_reduced_bam_set_t::iterator it = small_reduced_bam_set.begin();
  it!=small_reduced_bam_set.end();it++){
    small_reduced_bam_t srb = *it;
    if (srb.pos+2*MAX_READ_LEN<current_position){
      cout<<srb_read_index<<"\t0\t"<<srb.rname<<"\t"<<srb.pos<<"\t0\t"<<srb.short_cigar<<"\t*\t0\t0\t"<<srb.short_seq<<"\t*";
      cout<<"\tx1:H:";
      for(vector<int>::iterator it2 = srb.mapq_vector.begin();
      it2!=srb.mapq_vector.end();it2++){
        int2hex(*it2,cout);
        //cout<<","<<*it2;
      }
      cout<<"\ty1:H:";
      for(vector<int>::iterator it3 = srb.qual_vector.begin();
      it3!=srb.qual_vector.end();it3++){
        int2hex(*it3,cout);
        //cout<<","<<*it3;
      }
      cout<<endl;
      ++srb_read_index;
      srb_delete_vector.push_back(srb);
    }else{
      //cerr<<"SRB last pos is "<<srb.pos<<" and current_position is "<<current_position<<endl;
    }
  }
  for(vector<reduced_bam_t>::iterator it = reduced_bam_delete_vector.begin();
  it!=reduced_bam_delete_vector.end();it++){
    bam_map_t::iterator it2 = bam_map.find(*it);
    if (it2==bam_map.end()){
      throw "Should not happen that bam map cannot find deletion item";
    }else{
      bam_map.erase(it2);
    }
  }
  for(vector<small_reduced_bam_t>::iterator it = srb_delete_vector.begin();
  it!=srb_delete_vector.end();it++){
    small_reduced_bam_set_t::iterator it2 = small_reduced_bam_set.find(*it);
    if (it2==small_reduced_bam_set.end()){
      throw "Should not happen that srb set cannot find deletion item";
    }else{
      small_reduced_bam_set.erase(it2);
    }
    
  }
}

int main(int argc, char *argv[])
{
  try{
    tmpstruct_t tmp;
    if (argc < 5) {
      fprintf(stderr, "Usage: compress <subject> <variant_file> <in.bam> <in.fasta>\n");
      return 1;
    }
    int arg = 0;
    subject = atoi(argv[++arg]);   
    const char * variant_filename = argv[++arg];
    const char * bam_filename = argv[++arg];
    const char * fasta_filename = argv[++arg];
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
    tmp.fasta_index = fai_load(fasta_filename);
    if (tmp.fasta_index == 0) {
      fprintf(stderr, "Fail to open FASTA file %s\n", fasta_filename);
      return 1;
    } 
    parse_variants(variant_filename,tmp);
    bam_index_destroy(idx);
    samclose(tmp.in);
    return 0;
  }catch (const char * & mesg){
    cerr<<"Exception caught of message "<<mesg<<endl;
  }
}

