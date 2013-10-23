#include<map>
#include<vector>
#include<fstream>
#include<iostream>

using namespace std;

struct Config{
  static const int LIKELIHOOD_MODE_GENOTYPES = 0;
  static const int LIKELIHOOD_MODE_READS = 1;
  int model;
  bool use_reference_panel;
  string chromosome;
  int max_reference_haplotypes;
  string legendfile;
  string refhapfile;
  string inputformat;
  int g_likelihood_mode;
  string famfile;
  string bimfile;
  string bedfile;
  string sexfile;
  string bamfile;
  bool is_phased;
  bool is_sex_chr;
  string glf;
  bool use_gpu;
  bool use_cpu;
  int total_regions;
  int flanking_snps;
  int max_haplotypes;
  float delta;
  float lambda;
  bool debug;
  string format;
  string file_posterior;
  string file_genotype;
  string file_dosage;
  string file_quality;
  int platform_id;
  int device_id;
};

class IO_manager{
public:
  IO_manager();
  ~IO_manager();
  bool load_config(const char * config_filename);
  // for guide haplotypes glf 
  bool read_input(char * & hap_arr, float * & snp_penetrance, 
  bool * & informative_snp,int & persons, int & snps, int & total_ref_haps);
  // for guide haplotypes reads 
  bool read_input(char * & hap_arr, bool * & informative_snp,
  int & persons, int & snps, int & total_ref_haps);
  // for denovo haplotypes glf
  bool read_input(float * & snp_penetrance,int & persons, int & snps);
  // for denovo haplotypes reads
  bool read_input(int & persons, int & snps);
  //int get_total_persons();
  //int get_total_snps();
  //int get_total_refhaplotypes();
  vector<int> get_positions();
  Config * config;
  void writePosterior(int geno_len,int snp_index,float * val,int len);
  void writeGenotype(int snp_index,int * val,int len);
  void writeDosage(int snp_index,float * val,int len);
  void writeQuality(int snp_index,float val);

private:
  int snps,persons,refhaplotypes;
  //bool use_refhap;
  map<int,string> haplotype_map;
  map<int,string> genotype_map;
  map<int,vector<float> > penetrance_map;
  vector<int> position_vector;
  template<class ofs_T,class array_T> void write_output(ofs_T & ofs,int snp_index,int geno_len,int array_len, array_T * vals);
  bool read_input(char * & hap_arr, float * & snp_penetrance, 
  bool * & informative_snp,int * & haploid_arr,int & persons, int & snps, int & total_ref_haps);
  void parse_ref_haplotypes(const char * legend, const char * hapfile);
  void parse_plink(const char * famfile,const char * bimfile, const char * bedfile, bool * & informative_snp);
  void parse_glf(const char * famfile,const char * bimfile, bool is_phased, const char * glf,bool * & informative_snp);
  void load_bam_settings(const char * famfile, const char * polymorphisms, bool * & informative_snp);
  ofstream ofs_posterior_file;
  ofstream ofs_genotype_file;
  ofstream ofs_dosage_file;
  ofstream ofs_quality_file;
};
