#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>

using namespace std;

char * ref_sequence;
char * snp_alleles;
ofstream ofs_snpinfo;
ofstream ofs_variants;
vector<int> position_vec;
vector<string> haplo_vec;
int subjects;
int snps = 0;
const int MAX_SEQUENCE=300000000;
const int MAX_SNPS = 3000000;

float get_freq(string hapstr){
  float freq = 0;
  int n = hapstr.length();
  for(int i=0;i<n;++i){
    freq+=hapstr[i]=='1';
  }
  freq/=n;
  return freq;
}

char get_alt_base(char ref_base){
  int r = 3*(rand()/(RAND_MAX+1.));
  //cerr<<r;
  char alt_base;
  char a_alts[]={'C','G','T'};
  char c_alts[]={'A','G','T'};
  char g_alts[]={'A','C','T'};
  char t_alts[]={'A','C','G'};
  switch (ref_base){
    case 'A':
      alt_base = a_alts[r];
      break;
    case 'C':
      alt_base = c_alts[r];
      break;
    case 'G':
      alt_base = g_alts[r];
      break;
    case 'T':
      alt_base = t_alts[r];
      break;
    default:
      alt_base = 'N';
  }
  return alt_base;
}

void parse_kgp(){
  string line;
  ofs_snpinfo<<"CHROM\tPOS\tID\tREF\tALT\tFREQ\n";
  snps = 0;
  int last_basepos = -1;
  int offset = 1000;
  while(getline(cin,line)){
    istringstream iss(line);
    string id,haplo;
    int basepos;
    // positions are 1-based
    iss>>id>>basepos>>haplo;
    //cerr<<"orig base "<<basepos<<endl;
    if(snps==0) offset = basepos - offset; 
    basepos = basepos-offset+1;
    //cerr<<"offset "<<offset<<" snp "<<snps<<" basepos "<<basepos<<endl;
    // don't allow duplicate positions
    if (basepos<=last_basepos) basepos=last_basepos+1;
    position_vec.push_back(basepos);
    haplo_vec.push_back(haplo);
    ostringstream oss_common;
    //cerr<<"Fetching base at position "<<basepos-1<<endl;
    char ref_base = ref_sequence[basepos-1];
    char alt_base = get_alt_base(ref_base);
    snp_alleles[snps*2] = ref_base;
    snp_alleles[snps*2+1] = alt_base;
    oss_common<<"1\t"<<basepos<<"\t"<<id<<"\t"<<ref_base<<
    "\t"<<alt_base<<"\t"<<get_freq(haplo);
    ofs_snpinfo<<oss_common.str()<<endl;
    if (snps % 10000 == 0) cerr<<"site "<<snps<<" complete\n";
    last_basepos = basepos;
    ++snps;
  }
}

void write_variants(const char * variants_dir){
  for(int i=0;i<subjects;++i){
    ostringstream oss_outpath;
    oss_outpath<<variants_dir<<"/variant."<<i;
    ofs_variants.open(oss_outpath.str().data());
    ofs_variants<<"GENO\n";
    for(int j=0;j<snps;++j){
      int haplo_allele1 = (int)haplo_vec[j][i*2]-(int)'0';
      int haplo_allele2 = (int)haplo_vec[j][i*2+1]-(int)'0';
      ofs_variants<<snp_alleles[j*2+haplo_allele1]<<
      snp_alleles[j*2+haplo_allele2]<<endl;
    }
    ofs_variants.close();
  }
}

void load_ref_sequence(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Cannot find fasta file "<<filename<<endl;
    exit(1);
  }
  string line;
  getline(ifs,line);
  int cursor = 0;
  while(getline(ifs,line)){
    int strlen = line.length();
    for(int i=0;i<strlen;++i){
      ref_sequence[cursor+i] = line[i];
    }
    cursor+=strlen;
  }
  ifs.close();
 
}

int main(int argc,char * argv[]){
  if (argc<5){
    cout<<"Usage: <fasta_file> <snpinfo_file> <variants_dir> <subjects> \n";
    return 1;
  }
  srand(1);
  int arg=0;
  const char * fasta_file = argv[++arg];
  ref_sequence = new char[MAX_SEQUENCE];
  snp_alleles = new char[2*MAX_SNPS];
  load_ref_sequence(fasta_file);
  const char * snpinfo_file = argv[++arg];
  const char * variants_dir = argv[++arg];
  subjects = atoi(argv[++arg]);
  ofs_snpinfo.open(snpinfo_file);
  parse_kgp();
  ofs_snpinfo.close();
  write_variants(variants_dir);
  return 0;
}
