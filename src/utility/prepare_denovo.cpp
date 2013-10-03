#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<map>
#include<set>
#include<vector>
#include<math.h>
using namespace std;

const float epsilon = .001;
const float log_epsilon = log(epsilon);
const float log_third = log(.3);
//const char * filename_reflegend;
//const char * filename_refhap;
const char * filename_plink_fam;
const char * filename_plink_bim;
const char * filename_plink_bed;
const char * filename_glf_info;
const char * filename_glf_data;
const char * filename_mec_data;
int snps,persons;
//int selectchr,snps,persons;
int * genotypes_snpmajor;
float * penetrances_snpmajor;
const int PLINK_FORMAT = 0;
const int GLF_FORMAT = 1;
const int MEC_FORMAT = 2;
int inputformat;

//set<string> study_snps;
//map<int,int> haplo_positions;


struct study_record_t{
  int chr,genpos,phypos;
  string id,a1,a2;
  int * genotypes;
  float * penetrances;
  int newindex;
  study_record_t(int chr,string id,int genpos, int phypos,string a1,string a2){
    this->chr = chr;
    this->id = id;
    this->genpos = genpos;
    this->phypos = phypos;
    this->a1 = a1;
    this->a2 = a2;
    genotypes = NULL;
    penetrances = NULL;
    newindex = -1;
  }
  study_record_t(int chr,string id,int phypos,string a1,string a2){
    this->chr = chr;
    this->id = id;
    this->genpos = 0;
    this->phypos = phypos;
    this->a1 = a1;
    this->a2 = a2;
    genotypes = NULL;
    penetrances = NULL;
    newindex = -1;
  }
  study_record_t(int chr,int pos){
    this->chr = chr;
    phypos = pos;
  }
  void insert_genotypes(int * geno){
    this->genotypes = geno;
  }
  void insert_penetrances(float * pen){
    this->penetrances = pen;
  }
};

struct bypos{
  bool operator()(const study_record_t & r1,const study_record_t & r2){
    //return (r1.phypos<r2.phypos);
    return (r1.chr<r2.chr || r1.chr==r2.chr && r1.phypos<r2.phypos);
  }
};

vector<study_record_t> study_records;

void load_mec_study_positions(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
     cerr<<"Cannot open "<<filename<<endl;
     exit(1);
  }
  string line;
  int snp=0;
  while(getline(ifs,line)){
    istringstream iss(line);
    string id,a1="a1",a2="a2",chrstr;
    int chr,genpos=0,phypos;
    iss>>chrstr>>phypos>>id;
    if (chrstr.compare("X")==0){
      chr=23;
    }else{
      chr = atoi(chrstr.data());
    }
    study_record_t record(chr, id, genpos, phypos, a1,a2);
    record.insert_genotypes(genotypes_snpmajor+snp*persons);
    study_records.push_back(record);
    ++snp;
  }
  ifs.close();
  //cerr<<"Records loaded\n";
}

void load_plink_study_positions(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
     cerr<<"Cannot open "<<filename<<endl;
     exit(1);
  }
  string line;
  int snp=0;
  while(getline(ifs,line)){
    istringstream iss(line);
    string id,a1,a2;
    int chr,genpos,phypos;
    iss>>chr>>id>>genpos>>phypos>>a1>>a2;
    study_record_t record(chr, id, genpos, phypos, a1,a2);
    record.insert_genotypes(genotypes_snpmajor+snp*persons);
    study_records.push_back(record);
    ++snp;
  }
  ifs.close();
  cerr<<"Records loaded\n";
}

void load_glf_positions(){
  ifstream ifs(filename_glf_info);
  if (!ifs.is_open()){
     cerr<<"Cannot open "<<filename_glf_info<<endl;
     exit(1);
  }
  string line;
  int snp=0;
  while(getline(ifs,line)){
    istringstream iss(line);
    string id,a1,a2;
    int chr,genpos,phypos;
    iss>>chr>>phypos>>id>>a1>>a2;
    study_record_t record(chr, id,  phypos, a1,a2);
    record.insert_penetrances(penetrances_snpmajor+snp*3*persons);
    study_records.push_back(record);
    ++snp;
  }
  ifs.close();
  cerr<<"Records loaded\n";
}

//void load_haplo_data(){
//  ifstream ifs_legend(filename_reflegend);
//  if (!ifs_legend.is_open()){
//     cerr<<"Cannot open "<<filename_reflegend<<endl;
//     exit(1);
//  }
//  ifstream ifs_hap(filename_refhap);
//  if (!ifs_hap.is_open()){
//     cerr<<"Cannot open "<<filename_refhap<<endl;
//     exit(1);
//  }
//  string line,hapline;
//  getline(ifs_legend,line);
//  int i=0;
//  while(getline(ifs_legend,line)){
//    getline(ifs_hap,hapline);
//    istringstream iss(line);
//    istringstream iss2(hapline);
//    string id,a1,a2;
//    int pos;
//    iss>>id>>pos>>a1>>a2;
//    hap_record_t hap_record(id,pos,a1,a2);
//    int allele;
//    while(iss2>>allele){
//      hap_record.haplotypes.push_back(allele);
//    }
//    hap_records.push_back(hap_record);
//    //haplo_positions[pos] = i;
//    study_record_t record(selectchr,pos);
//    set<study_record_t>::iterator it = study_records.find(record);
//    study_record_t newrecord(selectchr, id, 0, pos, a1,a2);
//    if (it!=study_records.end()){
//      newrecord = *it;
//      newrecord.newindex = i;
//      newrecord.a1 = a1;
//      newrecord.a2 = a2;
//      study_records.erase(it);
//    }else{
//      //cerr<<"NOTE: Couldn't find study record at pos "<<pos<<endl;
//    }
//    study_records.insert(newrecord);
//    ++i;
//  }
//  ifs_legend.close();
//}

int colcount(const char * & filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
     cerr<<"Cannot open "<<filename<<endl;
     exit(1);
  }
  string line;
  int cols = 0;
  getline(ifs,line);
  ifs.close();
  istringstream iss(line);
  string token;
  while(iss>>token) ++cols;
  return cols;
}

int linecount(const char * & filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
     cerr<<"Cannot open "<<filename<<endl;
     exit(1);
  }
  int lines=0;
  string line;
  while(getline(ifs,line)){
    ++lines;
  }
  ifs.close();
  return lines;
}

void parse_mec(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Cannot open "<<filename<<endl;
    exit(1);
  }
  for (int snp=0;snp<snps;++snp){
    string line,chr,pos,rs,vec;
    getline(ifs,line);
    istringstream iss(line);
    iss>>chr>>pos>>rs>>vec;
    for(int i=0;i<persons;++i){
      genotypes_snpmajor[snp*persons+i] = (int)vec[i]-(int)'0';
      //cerr<<"snp,person "<<snp<<","<<i<<": "<< genotypes_snpmajor[snp*persons+i];
    }
  }
  ifs.close();
}

void parse_bed(const char * filename){
  int totalrows,totalcols;
  ifstream ifs(filename,ios::in|ios::binary);
  char header[3];
  ifs.read(header,3);
  bool snpmajor = (bool)header[2];
  cerr<<"SNP_MAJOR_MODE:\t"<<snpmajor<<endl;
  if (snpmajor){
    totalrows = snps;
    totalcols = persons;
  }else{
    totalrows = persons;
    totalcols = snps;
  }
  bool remainder = (totalcols % 4!=0)?true:false;
  int veclen = totalcols/4+remainder;
  char vector_read[veclen];
  int masks[4] = {3,12,48,192};
  char ** outmatrix;
  int mapping[] = {1,0,2,3};
  for (int row=0;row<totalrows;++row){
    ifs.read(vector_read,veclen);
    int col = 0;
    for (int byteindex =0;byteindex<veclen;++byteindex){
      for(int pair=0;pair<4;++pair){
        int val = (masks[pair] & vector_read[byteindex]) >> (2*pair);
        int geno = mapping[val];
        if (col<totalcols){
           int matindex = snpmajor? row*totalcols+col:  col*totalrows+row;
           genotypes_snpmajor[matindex] = geno;
           //cout<<geno;
        }
        ++col;
      }
    }
    //cout<<endl;
    //cerr<<"Row: "<<row<<" completed.\n";
  }
  ifs.close();
}

void parse_glf(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Could not open "<<filename<<endl;
    exit(1);
  }
  string line;
  for(int j=0;j<snps;++j){
    getline(ifs,line);
    istringstream iss(line);
    for(int i=0;i<persons;++i){
      for(int k=0;k<3;++k){
        float phred;
        iss>>phred;
        float pen = pow(10,phred/10);
        penetrances_snpmajor[j*3*persons+3*i+k] = pen<epsilon?
        log_epsilon:log(pen);
      }
    }
  }
  ifs.close();
}

void write_output(){
  ofstream ofs_geno1("out.genotype.info");
  ofstream ofs_geno2("out.genotype.data");
  ofs_geno2<<persons<<"\t"<<snps<<endl;
  ofs_geno1<<"chr\tid\tgenpos\tphypos\ta1\ta2"<<endl;
  for(vector<study_record_t>::iterator it = study_records.begin();
  it!=study_records.end(); it++){
    study_record_t record2 = *it;
    ofs_geno1<<record2.chr<<"\t"<<record2.id<<"\t"<<record2.genpos<<"\t"<<record2.phypos<<"\t"<<record2.a1<<"\t"<<record2.a2<<endl;
    if(inputformat==PLINK_FORMAT||inputformat==MEC_FORMAT){
      if (record2.genotypes==NULL){
        //cerr<<"NOTE: SNP "<<record2.id<<" is missing in study\n";
        for(int i=0;i<persons;++i){
          if (i) ofs_geno2<<"\t";
          ofs_geno2<<log_third<<"\t"<<log_third<<"\t"<<log_third;
        }
      }else{
        for(int i=0;i<persons;++i){
          if (i) ofs_geno2<<"\t";
          switch(record2.genotypes[i]){
            case 1:
              ofs_geno2<<0<<"\t"<<log_epsilon<<"\t"<<log_epsilon;
              break;
            case 2:
              ofs_geno2<<log_epsilon<<"\t"<<0<<"\t"<<log_epsilon;
              break;
            case 3:
              ofs_geno2<<log_epsilon<<"\t"<<log_epsilon<<"\t"<<0;
              break;
            default :
              ofs_geno2<<log_third<<"\t"<<log_third<<"\t"<<log_third;
              break;
          }
        }
      }      
    }else if(inputformat==GLF_FORMAT){
      if (record2.penetrances==NULL){
          for(int i=0;i<persons;++i){
            if (i) ofs_geno2<<"\t";
            ofs_geno2<<log_third<<"\t"<<log_third<<"\t"<<log_third;
          }
      }else{
        for(int i=0;i<persons;++i){
          for(int k=0;k<3;++k){
            if (i||k) ofs_geno2<<"\t";
            ofs_geno2<<record2.penetrances[i*3+k];
          }
        }      
      }
    }
    ofs_geno2<<endl;
  }
  ofs_geno1.close();
  ofs_geno2.close();
}

int mec_persons(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
    throw "Cannot find filename for mec data\n";
  }
  string line,chr,pos,rs,vec;
  getline(ifs,line);
  istringstream iss(line);
  iss>>chr>>pos>>rs>>vec;
  ifs.close();
  return vec.length();
}

int main(int argc,char * argv[]){
  if (argc<2){
    cout<<"Inputfiles: [0=plink|1=glf|2=mec] [<plink fam> <plink bim> <plink bed> | <glf info> <glf data>]\n";
    return 1;
  }
  int arg=1;
  //selectchr = atoi(argv[arg++]);
  //filename_reflegend = argv[arg++];
  //filename_refhap = argv[arg++];
  inputformat = atoi(argv[arg++]);  
  if (inputformat==PLINK_FORMAT){
    filename_plink_fam = argv[arg++];
    filename_plink_bim = argv[arg++];
    filename_plink_bed = argv[arg++];
    snps = linecount(filename_plink_bim);
    persons = linecount(filename_plink_fam);
  }else if (inputformat==GLF_FORMAT){
    filename_glf_info = argv[arg++];
    filename_glf_data = argv[arg++];
    snps = linecount(filename_glf_info);
    int cols = colcount(filename_glf_data);
    if (cols % 3 !=0 ) {
      cerr<<"The number of cols in the glf data file needs to be multiple of 3\n";
     return 1;
    }
    persons = cols/3;
  }else if (inputformat==MEC_FORMAT){
    filename_mec_data = argv[arg++];
    //cerr<<"MEC format filename: "<<filename_mec_data<<endl;
    snps = linecount(filename_mec_data);
    persons = mec_persons(filename_mec_data);
    //cerr<<"SNPS/PERSON: "<<snps<<","<<persons<<endl;
  }else{
    cout<<"Unrecognized format "<<inputformat<<endl;
    return 1;
  }
  cerr<<"SNPS: "<<snps<<" persons: "<<persons<<endl;
  if (inputformat==PLINK_FORMAT){
    genotypes_snpmajor = new int[snps*persons];
    parse_bed(filename_plink_bed);
    load_plink_study_positions(filename_plink_bim);
  }else if (inputformat==MEC_FORMAT){
    genotypes_snpmajor = new int[snps*persons];
    parse_mec(filename_mec_data);
    load_mec_study_positions(filename_mec_data);
  }else if (inputformat==GLF_FORMAT){
    penetrances_snpmajor = new float[snps*3*persons];
    parse_glf(filename_glf_data);
    load_glf_positions();
  }
  write_output();
  return 0;
}
