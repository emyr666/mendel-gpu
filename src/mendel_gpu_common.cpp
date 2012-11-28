using namespace std;
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<tr1/unordered_set>
#include<tr1/unordered_map>
#include<set>
#include<map>
#include<list>
#include<math.h>
#include"cl_constants.h"
#include<cstdlib>
#include<string.h>
#ifdef USE_GPU
#include<CL/cl.hpp>
#include"clsafe.h"
#endif
#include"mendel_gpu.hpp"


MendelGPU::MendelGPU(){
  epsilon = 1e-10;
  FORMAT_DEFAULT="default";
  FORMAT_MEC="mec";
  logpenetrance_threshold = -10;
  debug_mode = false;
  run_gpu = true;
  run_cpu = false;
  debug_truth = false;
#ifdef USE_GPU
  commandQueue = NULL;
  context = NULL;
  program = NULL;
  kernel_simple = NULL;
  kernel_compute_weights = NULL;
  kernel_compute_weights_haploid = NULL;
  kernel_reduce_weights2 = NULL;
  kernel_impute_genotype_denovo = NULL;
  kernel_impute_genotype_denovo_haploid = NULL;
  kernel_impute_genotype_guide = NULL;
  kernel_impute_genotype_guide_haploid = NULL;
  kernel_precompute_penetrance = NULL;
  kernel_precompute_penetrance_haploid = NULL;
  kernel_impute_penetrance = NULL;
  kernel_impute_penetrance_haploid = NULL;
  kernel_precompute_penetrance_fast = NULL;
  kernel_precompute_penetrance_fast_haploid = NULL;
  buffer_simple_in = NULL;
  buffer_simple_out = NULL;
  buffer_markers = NULL;
  buffer_haplotypes = NULL;
  buffer_extended_haplotypes = NULL;
  buffer_prev_left_marker = NULL;
  buffer_left_marker = NULL;
  buffer_haplotype = NULL;
  buffer_region_snp_penetrance = NULL;
  buffer_region_snp_offset = NULL;
  buffer_frequency = NULL;
  buffer_extended_frequency = NULL;
  buffer_subject_haplotype_weight = NULL;
  buffer_haplotype_weight = NULL;
  buffer_center_snp = NULL;
  buffer_subject_genotype = NULL;
  buffer_subject_genotype_block = NULL;
  buffer_subject_dosage = NULL;
  buffer_subject_dosage_block = NULL;
  buffer_subject_posterior_prob = NULL;
  buffer_subject_posterior_prob_block = NULL;
  buffer_logpenetrance_cache = NULL;
  buffer_max_penetrance = NULL;
  buffer_genotype_imputation = NULL;
  buffer_active_haplotype = NULL;
  buffer_twin_hap_index = NULL;
  buffer_penetrance_cache = NULL;
  buffer_frequency_cache = NULL;
  buffer_beyond_left_edge_dosage = NULL;
  buffer_right_edge_dosage = NULL;
  buffer_center_dosage = NULL;
#endif
}

MendelGPU::~MendelGPU(){
  ofs_genotype_file.close();
  ofs_dosage_file.close();
  ofs_quality_file.close();
  if (debug_mode) ofs_debug_haplotype_file.close();;
  if(debug_mode) ofs_posterior_file.close();
  delete[] logpenetrance_cache;
  delete[] g_current_weight;
}

string MendelGPU::hapstr(int * haplotype, int marker_len){
  string curhap(marker_len,'0');
  for(int j=0;j<marker_len;++j){
    curhap[j] = hapint2char(haplotype[j]);
  }
  return curhap;
}

int MendelGPU::extractBit(packedhap_t * packed,int index){
  int majorindex = index/32;
  int minorindex = (index%32)/8;
  int bitindex = (index%32)%8;
  cerr<<"extracting val for index: "<<index<<": "<<majorindex<<","<<minorindex<<","<<bitindex<<endl;
  int val = ((int)packed[majorindex].octet[minorindex]) >> bitindex & 1;
  return val;
}

void MendelGPU::decompress(int hapindex,int markers,int * haplotype){
  //cerr<<"Extracting hapindex: "<<hapindex<<endl;
  for(int site = 0;site<markers;++site){
    int majorindex = site/32;
    int minorindex = (site%32)/8;
    int bitindex = (site%32)%8;
    //cerr<<"extracting val for index: "<<site<<": "<<majorindex<<","<<minorindex<<","<<bitindex<<endl;
    haplotype[hapindex*g_max_window+site] = ((int)packedhap[hapindex*packedhap_len+majorindex].octet[minorindex]) >> bitindex & 1;
  }
}

void MendelGPU::injectBit(packedhap_t * packed,int index,int i){
  int majorindex = index/32;
  int minorindex = (index%32)/8;
  int bitindex = (index%32)%8;
  cerr<<"inserting "<<i<<" for index: "<<index<<": "<<majorindex<<","<<minorindex<<","<<bitindex<<endl;
  if (bitindex==0) packed[majorindex].octet[minorindex] = 0;
  //cerr<<"before: "<<i<<endl;
  packed[majorindex].octet[minorindex] |= (i<<bitindex);
  //cerr<<"after: "<<i<<endl;
}

void MendelGPU::compress(int hapindex,int markers,int * haplotype){
  //cerr<<"Compressing hapindex: "<<hapindex<<endl;
  for(int site = 0;site<markers;++site){
    int majorindex = site/32;
    int minorindex = (site%32)/8;
    int bitindex = (site%32)%8;
    //cerr<<"inserting "<<i<<" for index: "<<index<<": "<<majorindex<<","<<minorindex<<","<<bitindex<<endl;
    if (bitindex==0) packedhap[hapindex*packedhap_len+majorindex].octet[minorindex] = 0;
    packedhap[hapindex*packedhap_len+majorindex].octet[minorindex] |= (haplotype[hapindex*g_max_window+site]<<bitindex);
  }
}

void MendelGPU::compress_extendedhap(int extendedhapindex,int begin,int end,int * extendedhaplotype){
  for(int site = begin;site<end;++site){
    int packedsite = site-begin;
    int majorindex = packedsite/32;
    int minorindex = (packedsite%32)/8;
    int bitindex = (packedsite%32)%8;
    if (bitindex==0) packedextendedhap[extendedhapindex*packedextendedhap_len+majorindex].octet[minorindex] = 0;
    packedextendedhap[extendedhapindex*packedextendedhap_len+majorindex].octet[minorindex] |= (extendedhaplotype[extendedhapindex*g_max_window+site]<<bitindex);
  }
}

void MendelGPU::testpacked(){
  int testhap[60];
  packedhap_t testpack[3];
  for(int i=0;i<60;++i){
    testhap[i] = i%3==0;
    //cerr<<testhap[i];
    injectBit(testpack,i,testhap[i]);
  }
  cerr<<endl;
  ofstream ofs_debug_bin("debugbin");
  for(int i=0;i<3;++i){
    for(int j=0;j<4;++j){
      ofs_debug_bin.write(testpack[i].octet,4);
    }
  }
  ofs_debug_bin.close();
  for(int i=0;i<60;++i){
    testhap[i] = extractBit(testpack,i);
  }
  for(int i=0;i<60;++i){
    cerr<<testhap[i];
  }
  cerr<<endl;
  exit(0);
}

void MendelGPU::debug_haplotypes(ostream & os){
  for(int i=0;i<g_max_haplotypes;++i){
    if (g_active_haplotype[i]){
      os<<i<<":";
      for(int j=0;j<*g_markers;++j){
        os<<g_haplotype[i*g_max_window+j];
      }
      os<<endl;
    }
  }
}


void MendelGPU::load_constants_(int * model, int * people,int * snps, int * total_regions,int * haplotype_mode, int * flanking_snps,int * max_haplotypes,int * max_extended_haplotypes,int * platform_id, int * device_id, float * delta, float * lambda, int * i_geno_dim){
  ifstream ifs("settings.txt");
  if (!ifs.is_open()) {
    cout<<"Cannot load settings.txt\n";
    exit(1);
  }
  string line;
  while(getline(ifs,line)){
    istringstream iss(line);
    string tag;
    string strval;
    iss>>tag;
    g_likelihood_mode = LIKELIHOOD_MODE_GENOTYPES;
    if (tag.compare("PEOPLE")==0){
      iss>>people[0];
    }else if (tag.compare("BINARY_SNPS")==0){
      iss>>snps[0];
    }else if (tag.compare("TOTAL_REGIONS")==0){
      iss>>total_regions[0];
    }else if (tag.compare("HAPLOTYPE_MODE")==0){
      iss>>haplotype_mode[0];
      g_haplotype_mode = haplotype_mode[0];
    }else if (tag.compare("PLATFORM_ID")==0){
      iss>>platform_id[0];
    }else if (tag.compare("DEVICE_ID")==0){
      iss>>device_id[0];
    }else if (tag.compare("FLANKING_SNPS")==0){
      iss>>flanking_snps[0];
    }else if (tag.compare("MAX_HAPLOTYPES")==0){
      iss>>max_haplotypes[0];
    }else if (tag.compare("MAX_EXTENDED_HAPLOTYPES")==0){
      iss>>max_extended_haplotypes[0];
    }else if (tag.compare("MODEL")==0){
      iss>>model[0];
    }else if (tag.compare("DELTA")==0){
      iss>>delta[0];
    }else if (tag.compare("LAMBDA")==0){
      iss>>lambda[0];
    }else if (tag.compare("DEBUG_MODE")==0){
      iss>>debug_mode;
    }else if (tag.compare("RUN_GPU")==0){
      iss>>run_gpu;
    }else if (tag.compare("RUN_CPU")==0){
      iss>>run_cpu;
//    }else if (tag.compare("KERNEL_PATH")==0){
//      iss>>kernel_path;
    }else if (tag.compare("INFILE_REFHAP")==0){
      iss>>infile_refhap;
    }else if (tag.compare("INFILE_SEX")==0){
      iss>>infile_sex;
    }else if (tag.compare("IS_SEX_CHR")==0){
      iss>>is_sex_chr;
    }else if (tag.compare("INFILE_GENO_FILEFORMAT")==0){
      string format;
      iss>>format;
      if (format.compare("bam")==0){
        g_likelihood_mode = LIKELIHOOD_MODE_READS;
      }
    }else if (tag.compare("INFILE_GENO")==0){
      iss>>infile_geno;
    }else if (tag.compare("INFILE_GENO_DIM")==0){
      iss>>geno_dim;
      i_geno_dim[0] = geno_dim;
    }else if (tag.compare("OUTFILE_FORMAT")==0){
      iss>>outfile_format;
    }else if (tag.compare("OUTFILE_GENO")==0){
      iss>>outfile_geno;
    }else if (tag.compare("OUTFILE_DOSAGES")==0){
      iss>>outfile_dosages;
    }else if (tag.compare("OUTFILE_QUALITY")==0){
      iss>>outfile_quality;
    }else if (tag.compare("OUTFILE_POSTERIOR")==0){
      iss>>outfile_posterior;
    }else{
      iss>>strval;
    }
  }
  ifs.close();
}

// function called from main fortran program
//
void MendelGPU::init_gpu_(int * platform_id, int * device_id){
  imputation_software = getenv("IMPUTATION_SOFTWARE");
  
  if (debug_mode){ 
    ofs_debug_haplotype_file.open("DEBUG_HAPLOTYPE");
  }
  ofs_posterior_file.open(outfile_posterior.data());
  ofs_genotype_file.open(outfile_geno.data());
  ofs_dosage_file.open(outfile_dosages.data());
  ofs_quality_file.open(outfile_quality.data());
  if(run_gpu){
#ifdef USE_GPU
    if (commandQueue!=NULL) return;
    cerr<<"Initializing GPU\n";
    cerr<<"Fortran program initializing GPU with platform id "<<*platform_id<<" and device id "<<*device_id<<".\n";
    vector<cl::Platform> platforms;
    err = cl::Platform::get(&platforms);
    cerr<<"Platform ID "<<*platform_id<<" has name "<<platforms[*platform_id].getInfo<CL_PLATFORM_NAME>().c_str()<<endl;
    cl_context_properties cps[3] = {CL_CONTEXT_PLATFORM,(cl_context_properties)(platforms[*platform_id])(),0};
    context = new cl::Context(CL_DEVICE_TYPE_GPU,cps);
    devices = context->getInfo<CL_CONTEXT_DEVICES>();
    cerr<<"There are "<<devices.size()<<" devices\n";
    cerr<<"Device ID "<<*device_id<<" is of type: ";
    cl_device_type dtype = devices[*device_id].getInfo<CL_DEVICE_TYPE>();
    switch(dtype){
      case CL_DEVICE_TYPE_GPU:
        cerr<<"GPU\n";
        break;
      case CL_DEVICE_TYPE_CPU:
        cerr<<"CPU\n";
        break;
    } 
    commandQueue = new cl::CommandQueue(*context,devices[*device_id],0,&err);
    string sources[]={"cl_constants.h","kernels_common.c","kernels_denovo.c","kernels_guide.c"};
    ostringstream full_source_str;
    for(int j=0;j<4;++j){
      ostringstream oss;
      oss<<imputation_software<<"/src/"<<sources[j];
      string full_path = oss.str();
      cerr<<"Opening "<<full_path<<endl;
      ifstream ifs(full_path.data());
      clSafe(ifs.is_open()?CL_SUCCESS:-1,"kernel_path not found");
      string source_str(istreambuf_iterator<char>(ifs),(istreambuf_iterator<char>()));
      full_source_str<<source_str;
    }
    string source_str = full_source_str.str();
    // create a program object from kernel source
    cl::Program::Sources source(1,make_pair(source_str.c_str(),source_str.length()+1));
    program = new cl::Program(*context,source);
    err = program->build(devices);
    if(err!=CL_SUCCESS){
      cerr<<"Build failed:\n";
      string buffer;
      program->getBuildInfo(devices[0],CL_PROGRAM_BUILD_LOG,&buffer);
      cerr<<buffer<<endl;
    }
    // CREATE KERNELS
    kernel_simple = new cl::Kernel(*program,"simple",&err);
    clSafe(err,"creating kernel simple");
    // CREATE BUFFERS
    buffer_simple_in = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer in");
    buffer_simple_out = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * BLOCK_WIDTH, NULL, &err);
    clSafe(err,"creating buffer out");
    // SET KERNEL ARGUMENTS
    int arg;
    int kernelWorkGroupSize;
    arg = 0;
    err = kernel_simple->setArg(arg++, *buffer_simple_in);
    clSafe(err,"set kernel arg for kernel_simple");
    err = kernel_simple->setArg(arg++, *buffer_simple_out);
    clSafe(err,"set kernel arg for kernel_simple");
    kernelWorkGroupSize = kernel_simple->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel simple");
    cerr<<"simple kernel work group size is "<<kernelWorkGroupSize<<endl;
    cerr<<"GPU initialized\n";
#endif
  }
}

void MendelGPU::run_simple_(int * scaler, int * return_vec){
#ifdef USE_GPU
  cl_int err;
  //cerr<<"Scaler: "<<scaler[0];
  for(int i=0;i<256;++i){
    //return_vec[i]*=scaler[0];
    //cerr<<i<<":"<<return_vec[i]<<" ";
  }
  cerr<<endl;
  err = commandQueue->enqueueWriteBuffer(*buffer_simple_in, CL_TRUE, sizeof(int)*0,  sizeof(int)*1, scaler, NULL, NULL );
  clSafe(err, "write scaler");
  err = commandQueue->enqueueWriteBuffer(*buffer_simple_out, CL_TRUE, sizeof(int)*0,  sizeof(int)*BLOCK_WIDTH, return_vec, NULL, NULL );
  clSafe(err, "write vector");
  err = commandQueue->enqueueNDRangeKernel(*kernel_simple,cl::NullRange,cl::NDRange(BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
  clSafe(err,"launch simple kernel");
  cerr<<"launched kernel\n";
  err = commandQueue->enqueueReadBuffer(*buffer_simple_out, CL_TRUE, 0, BLOCK_WIDTH*sizeof(int),return_vec);
  clSafe(err, "read test vec");
  cerr<<"Read buffer\n";
#endif
}

void MendelGPU::unmarshall(ifstream & ifs_cache, float * outputvec,int len){
  char * charvec = reinterpret_cast<char * >(outputvec);
  ifs_cache.read(charvec,sizeof(float) * len);
  outputvec = reinterpret_cast<float * > (charvec);
}

void MendelGPU::read_stream_(int * people, int * snps, float * snp_penetrance){
  int total_persons = people[0];
  int total_snps = snps[0];
  ifstream ifs_cache;
  ifs_cache.open(infile_geno.data(),ios::in|ios::binary);
  if (!ifs_cache.is_open()){
    cerr <<"Cannot find the penetrance stream file\n";
    exit(1);
  }
  cerr<<"Opened stream inputfile of geno_dim "<<geno_dim<<" for "<<total_persons<<" persons and "<<total_snps <<" snps."<<"\n";
  int floatlen = geno_dim * total_snps;
  float logpen_threshold;
  if (geno_dim==3) logpen_threshold = log(.3);
  else if (geno_dim==4) logpen_threshold = log(.5);
  else {
    cerr<<"Invalid geno_dim, could not establish log pen threshold\n";
    exit(1);
  }
  cerr<<"Log penetrance threshold for informativeness is "<<logpen_threshold<<endl;
  for(int i=0;i< total_persons;++i){
    unmarshall(ifs_cache,snp_penetrance + i * geno_dim * total_snps,floatlen);
    for(int j=0;j<total_snps;++j){
      for(int k=0;k<geno_dim;++k){
        int index = i*floatlen + j*geno_dim + k;
        if(snp_penetrance[index] < -7) 
           snp_penetrance[index] = -7 ;
        if(snp_penetrance[index] == 0) 
           snp_penetrance[index] = -epsilon ;
      }
    }
  }
  ifs_cache.close();  
  cerr<<"Closed stream inputfile\n";
  informative_snp = new bool[total_snps];
  cerr<<"Informative SNPs from reference haplotypes:\n";
  for(int j=0;j<total_snps;++j){
    float max_logpen = -7;
    for(int i=0;i< total_persons;++i){
      for(int k=0;k<geno_dim;++k){
        float logpen = snp_penetrance[i * floatlen + j * geno_dim + k];
        if (logpen>max_logpen) max_logpen = logpen;
      }
    }
    informative_snp[j] = (max_logpen-logpen_threshold)>.01;
    if (informative_snp[j]) cerr<<" "<<j;
  }
  cerr<<endl;
  // Also take this opportunity to read in the person meta data
  cerr<<"Total persons "<<total_persons<<endl;
  haploid_arr = new int[total_persons];
  for(int i=0;i<total_persons;++i) haploid_arr[i] = 0;
  ifstream ifs_person;
  ifs_person.open(infile_sex.data());
  if (!ifs_person.is_open()){
    cerr <<"WARNING: Cannot find the person sex file "<<infile_sex<<". Assuming all females\n";
  }else{
    string line;
    getline(ifs_person,line);
    int person=0;
    int haploids = 0;
    for(int person = 0; person<total_persons;++person){
      getline(ifs_person,line);
      istringstream iss(line);
      int seq;
      char sex;
      iss>>seq>>sex;
      if (is_sex_chr && sex=='M') {
        haploid_arr[person] = 1;
        ++haploids;
      }
    }
    cerr<<"Total haploids in the analysis: "<<haploids<<endl;
    ifs_person.close();
  }
}

void MendelGPU::init_buffers_(int * active_haplotype,int * max_window,int * max_haplotypes,int * max_extended_haplotypes,int * max_region_size,int * people, int * snps, float * snp_penetrance, int * genotype_imputation,int * haplotypes, int * markers, int * window1,int * prev_left_marker, float * frequency, float * weight, int * haplotype,int * flanking_snps){
 cerr<<"Buffers initializing...\n";
 g_flanking_snps = *flanking_snps;
 g_active_haplotype = active_haplotype;
 g_prev_left_marker = prev_left_marker;
 g_left_marker = window1;
 g_markers = markers;
 g_haplotypes = haplotypes;
 g_haplotype = haplotype;
 g_frequency = frequency;
 g_weight = weight;
 g_genotype_imputation = *genotype_imputation;
 if (geno_dim==4 && !g_genotype_imputation){
   cerr<<"g_genotype_imputation: "<<g_genotype_imputation<<endl;
   cerr<<"Haplotyping not supported for haplotype input\n";
   exit(1);
 }
 g_max_region_size = *max_region_size;
 g_max_window = *max_window;
 g_max_haplotypes = *max_haplotypes;
 g_max_extended_haplotypes = *max_extended_haplotypes;
 packedhap_len = g_max_window/32 + (g_max_window%32!=0);
 packedhap = new packedhap_t[g_max_haplotypes * packedhap_len];
 g_people = *people;
 g_snps = *snps;
 g_informative_haplotype = new int[ref_haplotypes * g_max_window];
 extended_snp_mapping = new int[g_max_window];
 cerr<<"snps: "<<g_snps<<" people: "<<g_people<<" maxwindow: "<<g_max_window<<" max_haplotypes: "<<g_max_haplotypes<<" flanking: "<<g_flanking_snps<<endl;
 if (g_haplotype_mode == HAPLOTYPE_MODE_GUIDE){
   parse_ref_haplotype();
   extended_haplotype = new int[ref_haplotypes * g_max_window]; 
   //extended_center_dosage = new int[ref_haplotypes];
   extended_frequency = new float[ref_haplotypes];
   extended_root_mapping = new int[ref_haplotypes];
 }
 packedextendedhap_len = g_flanking_snps/32 + (g_flanking_snps%32!=0);
 packedextendedhap = new packedhap_t[ref_haplotypes * packedextendedhap_len];
 g_snp_penetrance = snp_penetrance;
 g_region_snp_penetrance = new float[g_people * geno_dim * g_max_region_size];
 cerr<<"Allocated SNP penetrance matrix of geno_dim "<<geno_dim<<endl;
 //genotype_posterior = new float[g_people * geno_dim * g_snps];
 penetrance_matrix_size = 0;
 if (geno_dim==4){ // for phased input
   logpenetrance_cache = new float[g_people*2*g_max_haplotypes];
   penetrance_cache = new float[g_people*2*g_max_haplotypes];
 }else if (geno_dim==3){
   penetrance_matrix_size = g_max_haplotypes*g_max_haplotypes;
   logpenetrance_cache = new float[g_people*penetrance_matrix_size];
   penetrance_cache = new float[g_people*penetrance_matrix_size];
   frequency_cache = new float[penetrance_matrix_size];
 }else{
   cerr<<"Invalid dim of "<<geno_dim<<endl;
   exit(1);
 }
 g_current_weight = new float[g_max_haplotypes];
 max_penetrance = new float[g_people];
 twin_hap_index = new int[g_max_haplotypes];
 beyond_left_edge_dosage = new int[g_max_haplotypes];
 right_edge_dosage = new int[g_max_haplotypes];
 center_dosage = new int[g_max_haplotypes];
 subject_haplotype_weight = new float[g_people*g_max_haplotypes];
 if(run_gpu){
#ifdef USE_GPU
    // CREATE ALL KERNELS HERE
    kernel_impute_penetrance = new cl::Kernel(*program,"impute_penetrance",&err);
    clSafe(err,"creating kernel impute penetrance");
    kernel_impute_penetrance_haploid = new cl::Kernel(*program,"impute_penetrance_haploid",&err);
    clSafe(err,"creating kernel impute penetrance_haploid");
    kernel_precompute_penetrance_fast = new cl::Kernel(*program,"precompute_penetrance_fast",&err);
    clSafe(err,"creating kernel precompute_fast penetrance");
    kernel_precompute_penetrance_fast_haploid = new cl::Kernel(*program,"precompute_penetrance_fast_haploid",&err);
    clSafe(err,"creating kernel precompute_fast penetrance haploid");
    kernel_precompute_penetrance = new cl::Kernel(*program,"precompute_penetrance",&err);
    clSafe(err,"creating kernel precompute");
    kernel_precompute_penetrance_haploid = new cl::Kernel(*program,"precompute_penetrance_haploid",&err);
    clSafe(err,"creating kernel precompute");
    kernel_compute_weights = new cl::Kernel(*program,"compute_weights",&err);
    clSafe(err,"creating kernel compute_weights");
    kernel_compute_weights_haploid = new cl::Kernel(*program,"compute_weights_haploid",&err);
    clSafe(err,"creating kernel compute_weights");
    kernel_reduce_weights2 = new cl::Kernel(*program,"reduce_weights2",&err);
    clSafe(err,"creating kernel reduce_weights2");
    kernel_impute_genotype_denovo = new cl::Kernel(*program,"impute_genotype_denovo",&err);
    clSafe(err,"creating kernel impute_genotype_denovo");
    kernel_impute_genotype_denovo_haploid = new cl::Kernel(*program,"impute_genotype_denovo_haploid",&err);
    clSafe(err,"creating kernel impute_genotype_denovo_haplod");
    kernel_impute_genotype_guide = new cl::Kernel(*program,"impute_genotype_guide",&err);
    clSafe(err,"creating kernel impute_genotype_guide");
    kernel_impute_genotype_guide_haploid = new cl::Kernel(*program,"impute_genotype_guide_haploid",&err);
    clSafe(err,"creating kernel impute_genotype_guide_haploid");
    // CREATE ALL BUFFERS HERE
    buffer_markers = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer markers");
    buffer_haplotypes = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer haplotypes");
    buffer_extended_haplotypes = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer extended haplotypes");
    buffer_left_marker = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer left marker");
    buffer_prev_left_marker = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer prev left marker");
    buffer_haplotype = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_max_window * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype buffer");
    buffer_region_snp_offset = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating snp offset buffer");
    buffer_region_snp_penetrance = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people * geno_dim * g_max_region_size , NULL, &err);
    clSafe(err,"creating snp penetrance buffer");
    buffer_frequency = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * *max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype frequency buffer");
    buffer_subject_haplotype_weight = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating subject haplotype weights buffer");
    buffer_haplotype_weight = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype weights buffer");
    buffer_center_snp = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer center marker");
    buffer_center_snp_start = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer center start marker");
    buffer_center_snp_end = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer center end marker");
    buffer_center_snp_offset = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer center offset marker");
    buffer_subject_genotype = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_people, NULL, &err);
    clSafe(err,"creating buffer subject genotype");
    buffer_subject_genotype_block = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_people*g_flanking_snps, NULL, &err);
    clSafe(err,"creating buffer subject genotype block");
    buffer_subject_dosage = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people, NULL, &err);
    clSafe(err,"creating buffer subject dosage");
    buffer_subject_dosage_block = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people*g_flanking_snps, NULL, &err);
    clSafe(err,"creating buffer subject dosage block");
    buffer_subject_posterior_prob = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*4, NULL, &err);
    clSafe(err,"creating buffer subject_posterior_prob");
    buffer_subject_posterior_prob_block = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*g_flanking_snps*4, NULL, &err);
    clSafe(err,"creating buffer subject_posterior_prob block");
    buffer_haploid_arr = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_people, NULL, &err);
    clSafe(err,"creating buffer haploid_arr");
    if (geno_dim==4){
      buffer_logpenetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*2*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer logpenetrance_cache");
      buffer_penetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*2*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer penetrance_cache");
      buffer_frequency_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer frequency_cache");
    }else if (geno_dim==3){
      buffer_logpenetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer logpenetrance_cache");
      buffer_penetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer penetrance_cache");
      buffer_frequency_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer frequency_cache");
    }
    buffer_max_penetrance = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people, NULL, &err);
    clSafe(err,"creating buffer max penetrance");
    buffer_genotype_imputation = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int), NULL, &err);
    clSafe(err,"creating buffer geno impute");
    buffer_twin_hap_index = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    clSafe(err,"creating buffer twin index");
    buffer_active_haplotype = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    clSafe(err,"creating buffer active hap");
    buffer_iteration = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*1, NULL, &err);
    clSafe(err,"creating buffer iteration");
    if (g_haplotype_mode==HAPLOTYPE_MODE_GUIDE){
      buffer_packedhap = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(packedhap_t) *g_max_haplotypes * packedhap_len, NULL, &err);
      clSafe(err,"creating packedhaplotype buffer");
      buffer_packedextendedhap = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(packedhap_t) *ref_haplotypes * packedextendedhap_len, NULL, &err);
      clSafe(err,"creating packedextendedhaplotype buffer");
      buffer_extended_root_mapping = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*ref_haplotypes, NULL, &err);
      clSafe(err,"creating buffer extended root mapping");
      buffer_extended_snp_mapping = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_window, NULL, &err);
      clSafe(err,"creating buffer extended snp mapping");
      buffer_extended_frequency = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * ref_haplotypes, NULL, &err);
      clSafe(err,"creating extended frequency buffer");
    }else if (g_haplotype_mode==HAPLOTYPE_MODE_DENOVO){
      buffer_center_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer center dosage");
      buffer_beyond_left_edge_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer beyond left edge dosage");
      buffer_right_edge_dosage = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer right edge dosage");
    }
    cerr<<"GPU Buffers created\n";
    // SET KERNEL ARGUMENTS HERE
    int arg;
    int kernelWorkGroupSize;
    // IMPUTE PENETRANCE MATRIX
    arg = 0;
    err = kernel_impute_penetrance->setArg(arg++, g_max_haplotypes);
    clSafe(err,"set kernel arg for kernel_impute_penetrance");
    err = kernel_impute_penetrance->setArg(arg++, penetrance_matrix_size);
    clSafe(err,"set kernel arg for kernel_impute_penetrance");
    err = kernel_impute_penetrance->setArg(arg++, *buffer_active_haplotype);
    clSafe(err,"set kernel arg for kernel_impute_penetrance");
    err = kernel_impute_penetrance->setArg(arg++, *buffer_twin_hap_index);
    clSafe(err,"set kernel arg for kernel_impute_penetrance");
    err = kernel_impute_penetrance->setArg(arg++, *buffer_logpenetrance_cache);
    clSafe(err,"set kernel arg for kernel_impute_penetrance");
    err = kernel_impute_penetrance->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_impute_penetrance");
    err = kernel_impute_penetrance->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_impute_penetrance");
    err = kernel_impute_penetrance->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_impute_penetrance");
    cerr<<"Total args for impute penetrance: "<<arg<<endl;
    kernelWorkGroupSize = kernel_impute_penetrance->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel impute_penetrance");
    cerr<<"impute penetrance kernel work group size is "<<kernelWorkGroupSize<<endl;
    // IMPUTE PENETRANCE MATRIX
    arg = 0;
    err = kernel_impute_penetrance_haploid->setArg(arg++, g_max_haplotypes);
    clSafe(err,"set kernel arg for kernel_impute_penetrance_haploid");
    err = kernel_impute_penetrance_haploid->setArg(arg++, penetrance_matrix_size);
    clSafe(err,"set kernel arg for kernel_impute_penetrance_haploid");
    err = kernel_impute_penetrance_haploid->setArg(arg++, *buffer_active_haplotype);
    clSafe(err,"set kernel arg for kernel_impute_penetrance_haploid");
    err = kernel_impute_penetrance_haploid->setArg(arg++, *buffer_twin_hap_index);
    clSafe(err,"set kernel arg for kernel_impute_penetrance_haploid");
    err = kernel_impute_penetrance_haploid->setArg(arg++, *buffer_logpenetrance_cache);
    clSafe(err,"set kernel arg for kernel_impute_penetrance_haploid");
    err = kernel_impute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_impute_penetrance_haploid");
    err = kernel_impute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_impute_penetrance_haploid");
    cerr<<"Total args for impute penetrance: "<<arg<<endl;
    kernelWorkGroupSize = kernel_impute_penetrance_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel impute_penetrance_haploid");
    cerr<<"impute penetrance kernel work group size is "<<kernelWorkGroupSize<<endl;
    // COMPUTE SUBJECT SPECIFIC HAPLOTYPE WEIGHTS
    arg = 0;
    err = kernel_compute_weights->setArg(arg++, g_max_haplotypes);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, g_max_window);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, penetrance_matrix_size);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_iteration);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_genotype_imputation);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_markers);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_haplotypes);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_left_marker);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_haploid_arr);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_haplotype);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_region_snp_penetrance);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    //err = kernel_compute_weights->setArg(arg++, *buffer_frequency_cache);
    //clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_frequency);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_penetrance_cache);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_subject_haplotype_weight);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, *buffer_active_haplotype);
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes)); // active haplotypes
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // known haplotype marginals
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // haplotype freqs
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
    clSafe(err,"set kernel arg for kernel_compute_weights");
    err = kernel_compute_weights->setArg(arg++, cl::__local(sizeof(float)*1));
    clSafe(err,"set kernel arg for kernel_compute_weights");
    cerr<<"Total args for kernel_compute_weights: "<<arg<<endl;
    kernelWorkGroupSize = kernel_compute_weights->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel compute_weights");
    cerr<<"compute_weights kernel work group size is "<<kernelWorkGroupSize<<endl;
    // COMPUTE SUBJECT SPECIFIC HAPLOTYPE WEIGHTS
    arg = 0;
    err = kernel_compute_weights_haploid->setArg(arg++, g_max_haplotypes);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, g_max_window);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, penetrance_matrix_size);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_iteration);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_genotype_imputation);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_markers);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_haplotypes);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_left_marker);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_haplotype);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_region_snp_penetrance);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    //err = kernel_compute_weights_haploid->setArg(arg++, *buffer_frequency_cache);
    //clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_frequency);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_penetrance_cache);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_subject_haplotype_weight);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, *buffer_active_haplotype);
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes)); // active haplotypes
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // known haplotype marginals
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // haplotype freqs
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    err = kernel_compute_weights_haploid->setArg(arg++, cl::__local(sizeof(float)*1));
    clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
    kernelWorkGroupSize = kernel_compute_weights_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel compute_weights_haploid");
    cerr<<"compute_weights_haploid kernel work group size is "<<kernelWorkGroupSize<<endl;

    // COMPUTE HAPLOTYPE WEIGHTS BY REDUCTION
    arg = 0;
    err = kernel_reduce_weights2->setArg(arg++, g_people);
    clSafe(err,"set kernel arg for kernel_reduce_weights2");
    err = kernel_reduce_weights2->setArg(arg++, g_max_haplotypes);
    clSafe(err,"set kernel arg for kernel_reduce_weights2");
    err = kernel_reduce_weights2->setArg(arg++, *buffer_haplotypes);
    clSafe(err,"set kernel arg for kernel_reduce_weights2");
    err = kernel_reduce_weights2->setArg(arg++, *buffer_subject_haplotype_weight);
    clSafe(err,"set kernel arg for kernel_reduce_weights2");
    err = kernel_reduce_weights2->setArg(arg++, *buffer_haplotype_weight);
    clSafe(err,"set kernel arg for kernel_reduce_weights2");
    err = kernel_reduce_weights2->setArg(arg++, *buffer_active_haplotype);
    clSafe(err,"set kernel arg for kernel_reduce_weights2");
    err = kernel_reduce_weights2->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_reduce_weights2");
    err = kernel_reduce_weights2->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes));
    clSafe(err,"set kernel arg for kernel_reduce_weights2");
    kernelWorkGroupSize = kernel_reduce_weights2->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel reduce_weights2");
    cerr<<"reduce_weights2 kernel work group size is "<<kernelWorkGroupSize<<endl;
    if (g_haplotype_mode==HAPLOTYPE_MODE_DENOVO){
      // PRECOMPUTE PENETRANCE FAST
      arg = 0;
      err = kernel_precompute_penetrance_fast->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, g_max_region_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, logpenetrance_threshold);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_region_snp_offset);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_markers);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_prev_left_marker);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_left_marker);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_haploid_arr);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_beyond_left_edge_dosage);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_right_edge_dosage);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_region_snp_penetrance);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_logpenetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, *buffer_active_haplotype);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, cl::__local(sizeof(float)*geno_dim));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, cl::__local(sizeof(float)*geno_dim));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      err = kernel_precompute_penetrance_fast->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast");
      cerr<<"Total args for precompute pen: "<<arg<<endl;
      kernelWorkGroupSize = kernel_precompute_penetrance_fast->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel precompute_penetrance_fast");
      cerr<<"precompute_penetrance_fast kernel work group size is "<<kernelWorkGroupSize<<endl;
      // PRECOMPUTE PENETRANCE FAST HAPLOID
      arg = 0;
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, g_max_region_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, logpenetrance_threshold);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_region_snp_offset);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_markers);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_prev_left_marker);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_left_marker);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_beyond_left_edge_dosage);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_right_edge_dosage);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_region_snp_penetrance);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_logpenetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, *buffer_active_haplotype);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, cl::__local(sizeof(float)*geno_dim));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, cl::__local(sizeof(float)*geno_dim));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      err = kernel_precompute_penetrance_fast_haploid->setArg(arg++, cl::__local(sizeof(float)*2*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_fast_haploid");
      cerr<<"Total args for precompute pen: "<<arg<<endl;
      kernelWorkGroupSize = kernel_precompute_penetrance_fast_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel precompute_penetrance_fast_haploid");
      cerr<<"precompute_penetrance_fast_haploid kernel work group size is "<<kernelWorkGroupSize<<endl;
      // IMPUTE UNPHASED GENOTYPE
      arg = 0;
      err = kernel_impute_genotype_denovo->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, g_max_window);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, epsilon);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_genotype_imputation);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_markers);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_left_marker);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_haploid_arr);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_frequency);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_subject_genotype);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_subject_dosage);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_subject_posterior_prob);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_active_haplotype);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_center_dosage);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, *buffer_subject_haplotype_weight);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes)); 
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // for loading haplotype frequencies
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(bestpair_t)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(float)*4));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      err = kernel_impute_genotype_denovo->setArg(arg++, cl::__local(sizeof(float)*1)); // normalizing constant
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo");
      kernelWorkGroupSize = kernel_impute_genotype_denovo->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_denovo");
      cerr<<"impute_genotype_denovo kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;

      // IMPUTE UNPHASED GENOTYPE HAPLOID
      arg = 0;
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, g_max_window);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, epsilon);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_genotype_imputation);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_markers);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_left_marker);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_frequency);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_subject_genotype);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_subject_dosage);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_subject_posterior_prob);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_active_haplotype);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_center_dosage);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, *buffer_subject_haplotype_weight);
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes)); 
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // for loading haplotype frequencies
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(float)*4));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes));
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      err = kernel_impute_genotype_denovo_haploid->setArg(arg++, cl::__local(sizeof(float)*1)); // normalizing constant
      clSafe(err,"set kernel arg for kernel_impute_genotype_denovo_haploid");
      kernelWorkGroupSize = kernel_impute_genotype_denovo_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_denovo_haploid");
      cerr<<"impute_genotype_denovo_haploid kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;
    } else if (g_haplotype_mode == HAPLOTYPE_MODE_GUIDE){
      // PRECOMPUTE PENETRANCE FOR THE GUIDE HAP APPROACH
      arg = 0;
      err = kernel_precompute_penetrance->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, g_max_window);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, g_max_region_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, logpenetrance_threshold);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, packedhap_len);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_region_snp_offset);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_markers);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_left_marker);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_haploid_arr);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_packedhap);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_haplotype);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_region_snp_penetrance);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_logpenetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, *buffer_extended_snp_mapping);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(packedhap_t)*g_max_haplotypes*packedhap_len));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(float)*g_max_window*geno_dim));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      err = kernel_precompute_penetrance->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance");
      cerr<<"Total args for precompute pen: "<<arg<<endl;
      kernelWorkGroupSize = kernel_precompute_penetrance->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel precompute_penetrance");
      cerr<<"precompute_penetrance kernel work group size is "<<kernelWorkGroupSize<<endl;
      arg = 0;
      err = kernel_precompute_penetrance_haploid->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, g_max_window);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, g_max_region_size);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, logpenetrance_threshold);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, packedhap_len);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_region_snp_offset);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_markers);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_left_marker);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_packedhap);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_haplotype);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_region_snp_penetrance);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_logpenetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, *buffer_extended_snp_mapping);
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(packedhap_t)*g_max_haplotypes*packedhap_len));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(int)*g_max_window));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_window*geno_dim));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      err = kernel_precompute_penetrance_haploid->setArg(arg++, cl::__local(sizeof(float)*2*BLOCK_WIDTH));
      clSafe(err,"set kernel arg for kernel_precompute_penetrance_haploid");
      cerr<<"Total args for precompute pen: "<<arg<<endl;
      kernelWorkGroupSize = kernel_precompute_penetrance_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel precompute_penetrance_haploid");
      cerr<<"precompute_penetrance_haploid kernel work group size is "<<kernelWorkGroupSize<<endl;
      arg = 0;
      err = kernel_impute_genotype_guide->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, epsilon);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, packedextendedhap_len);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, g_flanking_snps);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_genotype_imputation);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_extended_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_center_snp_end);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_haploid_arr);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_extended_root_mapping);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_extended_frequency);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_subject_genotype_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_subject_dosage_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_subject_posterior_prob_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_packedextendedhap);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, *buffer_subject_haplotype_weight);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(int)*ref_haplotypes)); // for root mapping
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*ref_haplotypes)); // for loading haplotype frequencies
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); // prob0
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob1
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob2
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob3
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*4));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // hap marginals
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      err = kernel_impute_genotype_guide->setArg(arg++, cl::__local(sizeof(float)*1)); // normalizing constant
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide");
      kernelWorkGroupSize = kernel_impute_genotype_guide->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_guide");
      cerr<<"impute_genotype_guide kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;

      arg = 0;
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, epsilon);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, packedextendedhap_len);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, g_flanking_snps);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_genotype_imputation);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_center_snp_end);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_root_mapping);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_frequency);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_genotype_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_dosage_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_posterior_prob_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_packedextendedhap);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_haplotype_weight);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(int)*ref_haplotypes)); // for root mapping
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*ref_haplotypes)); // for loading haplotype frequencies
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); // prob0
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob1
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob2
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob3
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*4));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // hap marginals
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*1)); // normalizing constant
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      kernelWorkGroupSize = kernel_impute_genotype_guide_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_guide_haploid");
      cerr<<"impute_genotype_guide_haploid kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;

      arg = 0;
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, g_max_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, epsilon);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, penetrance_matrix_size);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, packedextendedhap_len);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, g_flanking_snps);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_genotype_imputation);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_haplotypes);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_center_snp_end);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_root_mapping);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_extended_frequency);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_penetrance_cache);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_genotype_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_dosage_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_posterior_prob_block);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_packedextendedhap);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, *buffer_subject_haplotype_weight);
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(int)*ref_haplotypes)); // for root mapping
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*ref_haplotypes)); // for loading haplotype frequencies
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(packedhap_t)*ref_haplotypes*packedextendedhap_len));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); // prob0
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob1
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob2
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*BLOCK_WIDTH_IMPUTE_GUIDE)); //prob3
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*4));
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*g_max_haplotypes)); // hap marginals
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      err = kernel_impute_genotype_guide_haploid->setArg(arg++, cl::__local(sizeof(float)*1)); // normalizing constant
      clSafe(err,"set kernel arg for kernel_impute_genotype_guide_haploid");
      kernelWorkGroupSize = kernel_impute_genotype_guide_haploid->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
      clSafe(err,"get workgroup size kernel impute_genotype_guide_haploid");
      cerr<<"impute_genotype_guide_haploid kernel total args, "<<arg<<", work group size is "<<kernelWorkGroupSize<<endl;
    }


    cerr<<"GPU kernel arguments assigned.\n";

    // TRANSFER ANY DATA HERE
    err = commandQueue->enqueueWriteBuffer(*buffer_genotype_imputation, CL_TRUE, 0,sizeof(int), &g_genotype_imputation, NULL, NULL );
    clSafe(err, "write imputation mode");
    err = commandQueue->enqueueWriteBuffer(*buffer_haploid_arr, CL_TRUE, 0,sizeof(int)*g_people, haploid_arr, NULL, NULL );
    clSafe(err, "write haploid arr");
#endif
  }
  cerr<<"Buffers initialized\n";
}

void MendelGPU::init_region_buffers_(int * regionstart,int * regionend){
  int begin = g_region_snp_offset = *regionstart-1;
  int end = *regionend-1;
  int region_len = end-begin+1;
  for(int i=0;i<g_people;++i){
    for(int j=0;j<region_len;++j){
      for(int g=0;g<geno_dim;++g){
        g_region_snp_penetrance[i*(geno_dim*g_max_region_size)+j*geno_dim+g]
        = g_snp_penetrance[i*(geno_dim*g_snps)+(begin+j)*geno_dim+g];
      }
    }
  }
  cerr<<"Copied SNP penetrance sub buffer of region length "<<region_len<<endl;
  if(run_gpu){
#ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_region_snp_offset, CL_TRUE, 0,sizeof(int), &begin, NULL, NULL );
    clSafe(err, "write snp offset for region");
    err = commandQueue->enqueueWriteBuffer(*buffer_region_snp_penetrance, CL_TRUE, 0,sizeof(int)*g_people*geno_dim*g_max_region_size, g_region_snp_penetrance, NULL, NULL );
    clSafe(err, "write snp penetrance for region");
#endif
  }
  occupied_hap_indices.clear();
  free_hap_indices.clear();
  occupied_hap_indices.push_back(0);
  free_hap_indices.push_back(1);
  twin_hap_index[0] = 0;
  twin_hap_index[1] = 1;
}

void MendelGPU::init_window_buffers_(){
  char zero = '0';
  int offset = (int)zero;
  //g_haplotype = haplotype;
  int left_marker = g_left_marker[0]-1;
  int prev_left_marker = g_prev_left_marker[0]-1;

  bool debug_haplotype = false;
  if (debug_haplotype){
    
    if (debug_mode) ofs_debug_haplotype_file<<"ASSUMED DOUBLED HAPLOTYPES:\n";
    for(int i=0;i<g_max_haplotypes;++i){
      if (g_active_haplotype[i]){
        if (debug_mode) ofs_debug_haplotype_file<<i<<":";
        string curhap(*g_markers,zero);
        for(int j=0;j<*g_markers;++j){
          curhap[j] = (char)(g_haplotype[i*g_max_window+j]+offset); 
          if (debug_mode) ofs_debug_haplotype_file<<curhap[j];
        }
        if (debug_mode) ofs_debug_haplotype_file<<endl;
      }
    }
  }

  cerr<<"Max_window: "<<g_max_window<<" Haplotypes: "<<*g_haplotypes<<" Markers: "<<*g_markers<<" left_marker: "<<left_marker<<endl;
  if (run_gpu){
    #ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_active_haplotype, CL_TRUE, 0,  sizeof(int) * g_max_haplotypes, g_active_haplotype, NULL, NULL );
    clSafe(err, "write active haplotype indices");
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotypes, CL_TRUE, 0,  sizeof(int), g_haplotypes, NULL, NULL );
    clSafe(err, "write haplotypes");
    err = commandQueue->enqueueWriteBuffer(*buffer_markers, CL_TRUE, 0,  sizeof(int), g_markers, NULL, NULL );
    clSafe(err, "write markers");
    err = commandQueue->enqueueWriteBuffer(*buffer_left_marker, CL_TRUE, 0,  sizeof(int), &left_marker, NULL, NULL );
    clSafe(err, "write left marker");
    err = commandQueue->enqueueWriteBuffer(*buffer_prev_left_marker, CL_TRUE, 0,  sizeof(int), &prev_left_marker, NULL, NULL );
    clSafe(err, "write prev left marker");
    err = commandQueue->enqueueWriteBuffer(*buffer_haplotype, CL_TRUE, sizeof(int)*0,  sizeof(int)*g_max_window*g_max_haplotypes, g_haplotype, NULL, NULL );
    clSafe(err, "write haplotype");
    #endif
    cerr<<"Buffers sent to GPU\n";
  }
}

void MendelGPU::init_iteration_buffers_(){
  if (run_gpu){
#ifdef USE_GPU
  err = commandQueue->enqueueWriteBuffer(*buffer_frequency, CL_TRUE, 0,  sizeof(float)*g_max_haplotypes, g_frequency, NULL, NULL );
  clSafe(err, "write haplotype frequencies");
#endif
  cerr<<"Iteration Buffers sent to GPU\n";
  }
}


void MendelGPU::compute_haplotype_weights_(int * iteration){
  bool b_impute = g_genotype_imputation;
  cerr<<"begin compute_weights at iteration: "<<iteration[0]<<"\n";
  bool debug_personhap = false;
  //bool debug_personhap = iteration[0]==2;
  //bool debug_personhap = g_haplotypes[0]<-4;
  int debug_person = 10;
  bool debug_weights = g_haplotypes[0]<-4;
  bool debug_freq = g_haplotypes[0]<-4;
  if (run_gpu){
    double start = clock();
    #ifdef USE_GPU
    err = commandQueue->enqueueWriteBuffer(*buffer_iteration, CL_TRUE, 0,sizeof(int), iteration, NULL, NULL );
    clSafe(err, "write iteration");
    if (geno_dim==4){
      // for phased input
      err = commandQueue->enqueueNDRangeKernel(*kernel_compute_weights_haploid,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch compute_weights");
    }else if (geno_dim==3){
      // for unphased input
      err = commandQueue->enqueueNDRangeKernel(*kernel_compute_weights,cl::NullRange,cl::NDRange(g_people*BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
      clSafe(err,"launch compute_weights");
    }
    cerr<<"Launched compute_haplotype_weights\n";
    if (debug_personhap){
      float subject_haplotype_weight[g_people*g_max_haplotypes];
      err = commandQueue->enqueueReadBuffer(*buffer_subject_haplotype_weight, CL_TRUE, 0, sizeof(float)*g_people*g_max_haplotypes,subject_haplotype_weight);
      clSafe(err, "read subject_hap weight");
      float haplotype_weight[g_max_haplotypes];
      memset(haplotype_weight,0,sizeof(float)*g_max_haplotypes);
      for(int i = 0;i<g_people;++i){
        for(int j=0;j<g_max_haplotypes;++j){
          if (g_active_haplotype[j]){
            cout<<"GPU person "<<i<<" hap "<<j<<" weight "<<subject_haplotype_weight[i*g_max_haplotypes+j]<<endl;
            haplotype_weight[j]+=subject_haplotype_weight[i*g_max_haplotypes+j];
          }
        }
      }
    }
    err = commandQueue->enqueueNDRangeKernel(*kernel_reduce_weights2,cl::NullRange,cl::NDRange(BLOCK_WIDTH,1),cl::NDRange(BLOCK_WIDTH,1),NULL,NULL);
    clSafe(err,"launch reduce_weights");
    cerr<<"Launched reduce_weights\n";
    err = commandQueue->enqueueReadBuffer(*buffer_haplotype_weight, CL_TRUE, 0, sizeof(float)*g_max_haplotypes,g_weight);
    clSafe(err, "read hap weight");
    if (debug_weights){
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j]){
          cout<<"GPU hap "<<j<<" weight "<<g_weight[j]<<endl;
        }
      }
    }
    cerr<<"Elapsed time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    #endif
  }
  if(run_cpu){
    memset(g_weight,0,sizeof(float)*g_max_haplotypes);
    if (iteration[0]==1){
      memset(subject_haplotype_weight,0,sizeof(float)*g_people*g_max_haplotypes);
    }
    //for(int i=0;i<0;++i){
    double start = clock();
    if (debug_freq) cout<<"CPU freq: "<<endl;
    if (geno_dim==3){
      for(int j=0;j<g_max_haplotypes;++j){
        int start = b_impute?j:0;
        if (g_active_haplotype[j]){
          for(int k=start;k<g_max_haplotypes;++k){
            if (g_active_haplotype[k]){
              frequency_cache[j*g_max_haplotypes+k] = g_frequency[j]*g_frequency[k];
              if (j!=k) frequency_cache[j*g_max_haplotypes+k]*=2;
              if (debug_freq) cout<<" "<<frequency_cache[j*g_max_haplotypes+k];
            }
          }
        }
        if (debug_freq) cout<<endl;
      }
    }
    for(int i=0;i<g_people;++i){
      int haploid = haploid_arr[i];
      float likelihood = 0;
      memset(g_current_weight,0,sizeof(float)*g_max_haplotypes);
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] && (iteration[0]==1 || subject_haplotype_weight[i*g_max_haplotypes+j]>0)){
          if(geno_dim==4){
            for(int parent=0;parent<2;++parent){
              float penetrance = penetrance_cache[i*2*g_max_haplotypes+
              2*j+parent];
              if (penetrance>0){
                float freq = g_frequency[j];
                float p = freq * penetrance;
                likelihood+=p;
                g_current_weight[j]+=p;
              }
            }
          }else if (geno_dim==3){
            for(int k=j;k<g_max_haplotypes;++k){
              if (g_active_haplotype[k] &&
              (iteration[0]==1 || 
              subject_haplotype_weight[i*g_max_haplotypes+k]>0)){
                float penetrance = (haploid_arr[i] || j==k) ?
                penetrance_cache[i*penetrance_matrix_size+j*g_max_haplotypes+k]
                : 0 ;
                if (penetrance>0){
                  float freq = frequency_cache[j*g_max_haplotypes+k];
                  float p = freq*penetrance;
                  likelihood+=p;
                  g_current_weight[j]+=p;
                  g_current_weight[k]+=p;
                }
              }
            } 
          }
        }
      }
      for(int j=0;j<g_max_haplotypes;++j){
        if (g_active_haplotype[j] && (iteration[0]==1 || 
        subject_haplotype_weight[i*g_max_haplotypes+j]>0)){
            subject_haplotype_weight[i*g_max_haplotypes+j] = 
            g_current_weight[j];
        }
      }
      char zero='0';
      int offset=(int)zero;
      if (likelihood>0){
        for(int j=0;j<g_max_haplotypes;++j){
          if(g_active_haplotype[j]){
            g_current_weight[j]/=likelihood;   
            g_weight[j]+=g_current_weight[j];
            if (debug_personhap ) {
             cout<<"CPU person "<<i<<" hap "<<j<<" weight "<<g_current_weight[j]<<" marginal "<<subject_haplotype_weight[i*g_max_haplotypes+j]<<endl;
            }
          }
        }
      }else{
        //cerr<<"Likelihood for person "<<i<<" is zero\n";
        //exit(1);
      }
    }
    cerr<<"Compute haplotype weights: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    if (debug_weights){
      for(int j=0;j<g_max_haplotypes;++j){
        if(g_active_haplotype[j]){
          cout<<"CPU hap "<<j<<" ";
          for(int m2=0;m2<g_markers[0];++m2){
            cout<<g_haplotype[j*g_max_window+m2];
          } 
          cout<<" weight "<<g_weight[j]<<endl;
        }
      }
    }
  }//END CPU VERSION
  cerr<<"done do_iteration\n";
  if (debug_personhap||debug_weights) exit(1);
  if (debug_freq) exit(1);
}

float MendelGPU::compute_rsq(float * dosages,int stride, int index){
  float mean = 0;
  for(int i=0;i<g_people;++i){
    mean+=dosages[i*stride+index];
  }
  if (mean==0) return 0;
  mean/=(1.*g_people);
  float p = mean/2;
  float exp_var = 2*p*(1-p);
  if (exp_var==0) return 0;
  float obs_var = 0;
  for(int i=0;i<g_people;++i){
    obs_var+=pow(dosages[i*stride+index]-mean,2);
  }
  obs_var/=(g_people);
  float rsq = obs_var/exp_var;
  if (rsq>1) rsq = 1;
  return rsq;
}

//int main(){
int main2(){
#ifdef USE_GPU
  int platform_id = 0;
  int device_id = 0;
  const char * kernel_path="kernels.c";
  MendelGPU mendelGPU;
  mendelGPU.init_gpu_(&platform_id,&device_id);
  int scaler = 2;
  int return_vec[256];
  for(int i=0;i<256;++i){
    return_vec[i] = i;
  }
  mendelGPU.run_simple_(&scaler,return_vec);
  for(int i=0;i<256;++i){
    cout<<i<<":"<<return_vec[i]<<" ";
  }
  cout<<endl;
#endif
}

