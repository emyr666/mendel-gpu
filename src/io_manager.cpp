using namespace std;
#include<fstream>
#include<sstream>
#include<iostream>
#include<vector>
#include<map>
#include<bitset>
#include<boost/property_tree/ptree.hpp>
#include<boost/property_tree/xml_parser.hpp>
#include<boost/foreach.hpp>
#include<math.h>

#include"io_manager.hpp"

void checkOpen(const char * filename, ifstream & ifs){
  if (!ifs.is_open()){
    cerr<<"Cannot open "<<filename<<endl;
    throw "Program aborting";
  }
}


int linecount(const char * filename){
  int count=0;
  ifstream ifs(filename);
  if (!ifs.is_open()){
    cerr<<"Cannot open "<<filename<<endl;
    throw "Program aborting";
  }
  string line;
  while(getline(ifs,line)){
    ++count;
  }
  ifs.close();
  return count;
}

IO_manager::IO_manager(){
  config = new Config();
}

IO_manager::~IO_manager(){
  ofs_posterior_file.close();
  ofs_genotype_file.close();
  ofs_dosage_file.close();
  ofs_quality_file.close();
  delete config;
}

template <class ofs_T, class array_T> void IO_manager::write_output(ofs_T & ofs,int snp_index, int geno_len, int array_len,array_T * vals){
  ofs<<snp_index;
  for(int i=0;i<array_len;++i){
    for(int j=0;j<geno_len;++j){
      ofs<<"\t";
      ofs<<vals[i*geno_len+j];
    }
  }
  ofs<<endl;
}

void IO_manager::writePosterior(int geno_len,int snp_index,float * val,int len){
  write_output(ofs_posterior_file,snp_index,geno_len,len,val);
}

void IO_manager::writeGenotype(int snp_index,int * val,int len){
  write_output(ofs_genotype_file,snp_index,1,len,val);
}

void IO_manager::writeDosage(int snp_index,float * val,int len){
  write_output(ofs_dosage_file,snp_index,1,len,val);
}

void IO_manager::writeQuality(int snp_index,float  val){
  write_output(ofs_quality_file,snp_index,1,1,&val);
}


int IO_manager::get_total_snps(){
  return this->snps;
}

int IO_manager::get_total_persons(){
  return this->persons;
}

int IO_manager::get_total_refhaplotypes(){
  return this->refhaplotypes;
}

bool IO_manager::load_config(const char * xmlfile){
  cerr<<"Initializing configuration from XML file\n";
  ifstream ifs(xmlfile);
  if (!ifs.is_open()){
    cerr<<"Could not open "<<xmlfile<<" configuration file\n";
    return false;
  }
  ifs.close();
  try{
    boost::property_tree::ptree pt;
    boost::property_tree::read_xml(xmlfile,pt);
    config->model = pt.get<int>("model");
    config->chromosome = pt.get<string>("chromosome");
    config->use_reference_panel = pt.get<bool>("use_reference_panel");
    if (config->use_reference_panel){
      config->max_reference_haplotypes = pt.get<int>("reference_panel_settings.max_reference_haplotypes");
      config->legendfile = pt.get<string>("reference_panel_settings.legend");
      config->refhapfile = pt.get<string>("reference_panel_settings.haplotypes");
    }

    config->inputformat = pt.get<string>("input_type");
    cerr<<"Input format is "<<config->inputformat<<endl;
    if (config->inputformat.compare("plink")==0){
      config->g_likelihood_mode = Config::LIKELIHOOD_MODE_GENOTYPES;
      config->famfile = pt.get<string>("plink_settings.famfile");
      config->bimfile = pt.get<string>("plink_settings.bimfile");
      config->bedfile = pt.get<string>("plink_settings.bedfile");
    }else if (config->inputformat.compare("glf")==0){
      config->g_likelihood_mode = Config::LIKELIHOOD_MODE_GENOTYPES;
      config->famfile = pt.get<string>("glf_settings.famfile");
      config->bimfile = pt.get<string>("glf_settings.bimfile");
      config->is_phased = pt.get<bool>("glf_settings.phased_input");
      config->glf = pt.get<string>("glf_settings.glf");
    }else if (config->inputformat.compare("bam")==0){
      config->g_likelihood_mode = Config::LIKELIHOOD_MODE_READS;
      config->famfile = pt.get<string>("bam_settings.famfile");
      config->bimfile = pt.get<string>("bam_settings.bimfile");
      config->bam_manifest_file = pt.get<string>("bam_settings.bam_manifest_file");
      ifstream ifs_bamlist(config->bam_manifest_file.data());
      if (!ifs_bamlist.is_open()) throw "BAM manifest file not found";
      string bamfilepath; 
      while(getline(ifs_bamlist,bamfilepath)){
        config->bamfilevec.push_back(bamfilepath);
      }
      ifs_bamlist.close();
    }else{
      throw "Invalid input type.  Valid values are plink, glf, bam.";
    }
    config->use_gpu = pt.get<bool>("use_gpu");
    config->use_cpu = pt.get<bool>("use_cpu");
    config->total_regions = pt.get<int>("tuning_parameters.total_regions");
    config->flanking_snps = pt.get<int>("tuning_parameters.flanking_snps");
    config->max_haplotypes = pt.get<int>("tuning_parameters.max_haplotypes");
    config->delta = pt.get<float>("tuning_parameters.delta");
    config->lambda = pt.get<float>("tuning_parameters.lambda");
    config->debug = pt.get<bool>("tuning_parameters.debug");
    config->format = pt.get<string>("output_settings.format");
    config->file_posterior = pt.get<string>("output_settings.posterior");
    config->file_genotype = pt.get<string>("output_settings.genotype");
    config->file_dosage = pt.get<string>("output_settings.dosage");
    config->file_quality = pt.get<string>("output_settings.quality");
    config->platform_id = pt.get<int>("opencl_settings.platform_id");
    config->device_id = pt.get<int>("opencl_settings.device_id");
  }catch (const exception & e){
    cerr<<"Caught an exception attempting to parse XML: "<<e.what()<<"\n";
    return false;
  }
  cerr<<"Configuration loaded\n";
  try{
    ofs_posterior_file.open(config->file_posterior.data());
    ofs_genotype_file.open(config->file_genotype.data());
    ofs_dosage_file.open(config->file_dosage.data());
    ofs_quality_file.open(config->file_quality.data());
  }catch (const exception & e){
    cerr<<"Caught an exception attempting to open outfiles: "<<e.what()<<"\n";
    return false;
  }
  
  return true;
}

bool IO_manager::read_input(char * & haplotype_array, float * & snp_penetrance, bool * & informative_snp, int * & haploid_arr){
  try{
    if (config->use_reference_panel){
      haplotype_map.clear();
      cerr<<"Cleared haplotype map\n";
      parse_ref_haplotypes(config->legendfile.data(),config->refhapfile.data());
      for(uint j=0;j<position_vector.size();++j){
        int position = position_vector[j];
        string hapstr = haplotype_map[position];
        if (j==0){
          haplotype_array = new char[this->refhaplotypes * snps];
        }
        for(int h=0;h<this->refhaplotypes;++h){
          haplotype_array[h*snps+j] = hapstr[h];
        }
      }
    }
  }catch (const exception & e){
    cerr<<"Caught an exception attempting to parse ref haplotypes: "<<e.what()<<"\n";
    return false;
  }
  try{
    if (config->inputformat.compare("plink")==0){
      parse_plink(config->famfile.data(), config->bimfile.data(),
      config->bedfile.data(), informative_snp);
    }else if (config->inputformat.compare("glf")==0){
      parse_glf(config->famfile.data(), config->bimfile.data(),
      config->is_phased,config->glf.data(), informative_snp);
      //int geno_dim = config->is_phased?4:3;
    }else if (config->inputformat.compare("bam")==0){
      load_bam_settings(config->famfile.data(), 
      config->bimfile.data(), config->bamfilevec, informative_snp);
    }
  }catch (const exception & e){
    cerr<<"Caught an exception attempting to parse study data: "<<e.what()<<"\n";
    return false;
  }
  try{
    cerr<<"Loading penetrances...\n";
    float epsilon = 1e-10;
    float callprob = 1.-epsilon;
    float log_zero = log(epsilon);
    float log_one = log(callprob);
    if(config->inputformat.compare("plink")==0 || 
      config->inputformat.compare("glf")==0){
      int geno_dim = (config->inputformat.compare("glf")==0 && 
      config->is_phased)?4:3;
      float nocallprob = 1./geno_dim;
      float log_nocall = log(nocallprob);
      snp_penetrance = new float[persons * snps * geno_dim];
      vector<int>::iterator it;
      cerr<<"Records transposed to subject major:";
      for(uint j = 0;j<position_vector.size();++j){
        int position = position_vector[j];
        // In this case geno_dim is always 3
        if (config->inputformat.compare("plink")==0){
          string genostr = genotype_map[position];
          //cerr<<position;
          for(int i=0;i<persons;++i){
            int g = (int)genostr[i]-(int)'0';
            float g1=log_zero,g2=log_zero,g3=log_zero;
            switch(g){
            case 0:
              g1=log_nocall;
              g2=log_nocall;
              g3=log_nocall;
              break;
            case 1:
              g1=log_one;
              break;
            case 2:
              g2=log_one;
              break;
            case 3:
              g3=log_one;
              break;
            }
            snp_penetrance[i*snps*geno_dim+j*3] = g1;
            snp_penetrance[i*snps*geno_dim+j*3+1] = g2;
            snp_penetrance[i*snps*geno_dim+j*3+2] = g3;
          }
          //cerr<<endl;
        }else if(config->inputformat.compare("glf")==0){
          vector<float> penetrance_vec = penetrance_map[position];
          for(int i=0;i<persons;++i){
            for(int k=0;k<geno_dim;++k){
              snp_penetrance[i*snps*geno_dim+j*geno_dim+k] = log(penetrance_vec[i*geno_dim+k]);
            }
          }
        }
        if (j && j%10000==0) cerr<<" "<<j;
      }
      cerr<<endl;
      if (config->inputformat.compare("plink")==0){
        genotype_map.clear();
        cerr<<"Genotype map cleared\n";
      }else if(config->inputformat.compare("glf")==0){
        penetrance_map.clear();
        cerr<<"Penetrance map cleared\n";
      }
    }
  }catch (const exception & e){
    cerr<<"Caught an exception attempting to load into array: "<<e.what()<<"\n";
    return false;
  } 
  // attempt to read in list of haploid/diploids
  try{
    ifstream ifs_person;
    ifs_person.open(config->sexfile.data());
    haploid_arr = new int[persons];
    if (!ifs_person.is_open()){
      cerr <<"WARNING: Cannot find the person sex file "<<config->sexfile<<". Assuming all females\n";
      for(int i=0;i<persons;++i) haploid_arr[i] = 0;
    }else{
      string line;
      getline(ifs_person,line);
      int person=0;
      int haploids = 0;
      for(int person = 0; person<persons;++person){
        getline(ifs_person,line);
        istringstream iss(line);
        int seq;
        char sex;
        iss>>seq>>sex;
        haploid_arr[person] = (config->is_sex_chr && sex=='M')?1:0;
        haploids+=haploid_arr[person];
      }
      cerr<<"Total haploids in the analysis: "<<haploids<<endl;
      ifs_person.close();
    }
  }catch (const exception & e){
    cerr<<"Caught an exception attempting to read sex file: "<<e.what()<<"\n";
    return false;
  } 
  return true;
}

void IO_manager::parse_ref_haplotypes(const char * legendfile, const char * hapfile){
  ifstream ifs_legend(legendfile);
  ifstream ifs_hap(hapfile);
  checkOpen(legendfile,ifs_legend);
  checkOpen(hapfile,ifs_hap);
  cerr<<"Found legend "<<legendfile<<" and haplotypes "<<hapfile<<endl;
  string line_legend,line_hap;
  getline(ifs_legend,line_legend); // skip header
  int j = 0;
  cerr<<"Records loaded:";
  while(getline(ifs_legend,line_legend)){
    istringstream iss_legend(line_legend);
    string rs;
    int position;
    iss_legend>>rs>>position;
    getline(ifs_hap,line_hap);
    if (haplotype_map.find(position)==haplotype_map.end()){
      ostringstream oss_hap;
      istringstream iss_hap(line_hap);
      string allele;
      int h = 0; 
      while(iss_hap>>allele){
        oss_hap<<allele;
        ++h;
      }
      if (j==0){
        this->refhaplotypes = h;      
      }
      haplotype_map[position] = oss_hap.str();
      //cerr<<"Added hap pos "<<position<<endl;
      position_vector.push_back(position);
      ++j;
      if (j%10000==0) cerr<<" "<<j;
    }
  }
  cerr<<endl;
  this->snps = j; 
  ifs_legend.close();
  ifs_hap.close();
  cerr<<j<<" SNPs will be used for imputation\n";
}


void IO_manager::parse_plink(const char * famfile,const char * bimfile, const char * bedfile, bool * & informative_snp){
  this->persons = linecount(famfile);
  int study_snps = linecount(bimfile);

  ifstream ifs_bim(bimfile);
  string bim_line;
  vector<int> study_pos_vec;
  bool include_snp[study_snps];
  string chr_arr[study_snps];
  int tagsnp_indices[study_snps];
  int j=0,tagsnps=0;
  while(getline(ifs_bim,bim_line)){
    string snpid,morgan,a1,a2;
    int position;
    istringstream iss_bim(bim_line);
    iss_bim>>chr_arr[j]>>snpid>>morgan>>position>>a1>>a2;
    //cerr<<"Userefhap:"<<config->use_reference_panel<<" config->chr: "<<config->chromosome<<" chrarr "<<chr_arr[j]<<endl;
    //if (chr_arr[j].compare(config->chromosome)==0) cerr<<"Finding "<<position<<endl;
    if (!config->use_reference_panel || (chr_arr[j].compare(config->chromosome)==0 && haplotype_map.find(position)!=haplotype_map.end())){
      include_snp[j] = true;
      tagsnp_indices[j] = tagsnps;
      study_pos_vec.push_back(position);
      //cerr<<"Tag SNP at "<<j<<" at pos "<<position<<endl;
      ++tagsnps;
    }else{
      //if (chr_arr[j].compare(config->chromosome)==0 ) cerr<<"Skipping position "<<chr_arr[j]<<"-"<<position<<" in study as it is not present on ref panel\n";
      include_snp[j] = false;
    }
    ++j;
  }
  ifs_bim.close();
  cerr<<"Total study SNPs: "<<study_snps<<endl;
  cerr<<"Total intersecting SNPs: "<<tagsnps<<endl;
  int * genotypes_snpmajor = new int[tagsnps*persons];
  int totalrows,totalcols;
  ifstream ifs_bed(bedfile,ios::in|ios::binary);
  char header[3];
  ifs_bed.read(header,3);
  bool snpmajor = (bool)header[2];
  cerr<<"PLINK file is ";
  if (snpmajor){
    cerr<<"SNP major\n";
  }else{
    cerr<<"Subject major\n";
  }
  if (snpmajor){
    totalrows = study_snps;
    totalcols = persons;
  }else{
    totalrows = persons;
    totalcols = study_snps;
  }
  bool remainder = (totalcols % 4!=0)?true:false;
  int veclen = totalcols/4+remainder;
  char vector_read[veclen];
  int masks[4] = {3,12,48,192};
  //char ** outmatrix;
  int mapping[] = {1,0,2,3};
  cerr<<"Records read:";
  for (int row=0;row<totalrows;++row){
    if (snpmajor && !include_snp[row]){
      ifs_bed.seekg(veclen,ios_base::cur);
    }else{
      ifs_bed.read(vector_read,veclen);
      int col = 0;
      for (int byteindex =0;byteindex<veclen;++byteindex){
        for(int pair=0;pair<4;++pair){
          int val = (masks[pair] & vector_read[byteindex]) >> (2*pair);
          int geno = mapping[val];
          if (col<totalcols &&  (snpmajor || include_snp[col])){ 
             int matindex = snpmajor? tagsnp_indices[row]*totalcols+col:
             tagsnp_indices[col]*totalrows+row;
             genotypes_snpmajor[matindex] = geno;
             //cout<<geno;
          }
          ++col;
        }
      }
      //cout<<endl;
    }
    if (row && row%10000==0) cerr<<" "<<row;
  }
  cerr<<endl;
  ifs_bed.close();
  // first go through and add to genotype_map things that are present in the ref panel
  for(uint j=0;j<study_pos_vec.size();++j){
    int position = study_pos_vec[j];
    ostringstream oss_geno;
    for(int i=0;i<persons;++i){
      oss_geno<<genotypes_snpmajor[j*persons+i];
    }
    genotype_map[position] = oss_geno.str();
    // this line emulates what would have been done if
    // we parse referenced haplotypes
    if(!config->use_reference_panel) position_vector.push_back(position);
  }
  // this line emulates what would have been done if
  // we parse referenced haplotypes
  if (!config->use_reference_panel) this->snps = position_vector.size();
  // allocate memory for informative_snp
  if (config->use_reference_panel) {
    informative_snp = new bool[this->snps];
    cerr<<"Allocated informative SNP of size "<<this->snps<<endl;
  }
  // next, pad genotypes for those that are present in ref hap but not genotype_map
  ostringstream oss_geno;
  for(int i=0;i<persons;++i){
    oss_geno<<'0';
  }
  j = 0;
  for(vector<int>::iterator it = position_vector.begin();
  it!=position_vector.end();it++){
    int position = *it;
    if (config->use_reference_panel)informative_snp[j] = true;
    if (genotype_map.find(position)==genotype_map.end()){
      cerr<<"Will impute position "<<position<<" at index "<<j<<endl;
      genotype_map[position] = oss_geno.str();
      if (config->use_reference_panel) informative_snp[j] = false;
    }
    ++j;
  }
  delete[] genotypes_snpmajor;
  cerr<<"SNPs for imputation: "<<this->snps<<endl;
}

void IO_manager::parse_glf(const char * famfile,const char * bimfile, bool is_phased, const char * glf, bool * & informative_snp){
  int geno_dim = is_phased?4:3;
  this->persons = linecount(famfile);
  int study_snps = linecount(bimfile);
  //penetrance_t * penetrances_snpmajor = new int[snps*persons];
  ifstream ifs_bim(bimfile);
  ifstream ifs_glf(glf);
  string line_bim,line_glf;
  string chr_arr[study_snps];
  for(int j=0;j<study_snps;++j){
    getline(ifs_bim,line_bim);
    getline(ifs_glf,line_glf); 
    string snpid,morgan,a1,a2;
    int position;
    istringstream iss_bim(line_bim);
    iss_bim>>chr_arr[j]>>snpid>>morgan>>position>>a1>>a2;
    // if we are in denovo haplotyping mode or this is a tag snp store this
    // data, else ignore it.
    if(!config->use_reference_panel || (config->chromosome.compare(chr_arr[j])==0 && haplotype_map.find(position)!=haplotype_map.end())){
      istringstream iss_glf(line_glf);
      vector<float> pen_vec;
      for(int i=0;i<persons;++i){
        for(int k=0;k<geno_dim;++k){
          float f;
          iss_glf>>f;
          pen_vec.push_back(f);
        }
      }
      penetrance_map[position] = pen_vec;
    }
    // this line is needed since we don't call parse_ref_haplotypes
    if (!config->use_reference_panel) position_vector.push_back(position);
  }
  ifs_bim.close();
  ifs_glf.close();
  // this line is needed since we don't call parse_ref_haplotypes
  if (!config->use_reference_panel) this->snps = position_vector.size();
  // allocate memory for informative_snp
  if (config->use_reference_panel) informative_snp = new bool[this->snps];
  // next, pad penetrances for those that are present in ref hap but not penetrance_map
  float missing = 1./geno_dim;
  int j = 0;
  for(vector<int>::iterator it = position_vector.begin();
  it!=position_vector.end();it++){
    int position = *it;
    if (config->use_reference_panel) informative_snp[j] = true;
    if (penetrance_map.find(position)==penetrance_map.end()){
      cerr<<"Will impute position "<<position<<endl;
      vector<float> pen_vec;
      for(int i=0;i<persons;++i){
        for(int k=0;k<geno_dim;++k){
          pen_vec.push_back(missing);
        }
      }
      penetrance_map[position] = pen_vec;
      if (config->use_reference_panel) informative_snp[j] = false;
    }
    ++j;
  }
}

void IO_manager::load_bam_settings(const char * famfile, const char * polymorphisms, vector<string> bam_filelist,bool * & informative_snp){
  this->persons = linecount(famfile);
  ifstream ifs_poly(polymorphisms);
  if (!ifs_poly.is_open()) {
    cerr<<"Can't open "<<polymorphisms<<". ";
    throw "File not found.";
  }
  if (config->use_reference_panel){
    throw "Reference panels are only supported for non-read data!";
  }
  string line;
  while(getline(ifs_poly,line)){
    istringstream iss_bim(line);
    string token;
    int position;
    iss_bim>>token>>token>>token>>position>>token>>token;
    position_vector.push_back(position);
  }
  ifs_poly.close();
  this->snps = position_vector.size();
  // allocate memory for informative_snp
  if (config->use_reference_panel){
    informative_snp = new bool[this->snps];
    for(int j=0;j<this->snps;++j) informative_snp[j] = true;
  }
}

int main2(int argc,char * argv[]){
  IO_manager * manager = new IO_manager();
  manager->load_config("settings.xml");
  char * haplotypes;
  float * snp_penetrance;
  bool * info_snp;
  int * hap_arr;
  manager->read_input(haplotypes,snp_penetrance,info_snp,hap_arr);
  int snps = manager->get_total_snps();
  int persons = manager->get_total_persons();
  for(int i=0;i<persons;++i){
    for(int j=0;j<snps;++j){
      cerr<<" ";
      for(int k=0;k<3;++k){
        if(k) cerr<<",";
        cerr<<snp_penetrance[i*3*snps+3*j+k];
      }
    }
    cerr<<endl;
  }
  delete manager;
  return 0;
}
