#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"

MendelGPU::MendelGPU(IO_manager * io){
  this->io_manager = io;
  this->config = io->config;
}

void MendelGPU::run_sliding_window(){
  cerr<<"Running sliding window on "<<g_snps<<" SNPs with "<<g_flanking_snps<<" flanking SNPs\n";

  g_markers = 0;
  cerr<<"center snp start,right_marker,left_marker: "<<g_center_snp_start<<","<<g_right_marker<<","<<g_left_marker<<endl;
  //while(g_center_snp_start<10){
  while(g_center_snp_start<g_snps){
    int window_size = g_right_marker-g_left_marker+1;
    cerr<<"*** New window "<<g_left_marker<<" to "<<g_right_marker<<", snps "<<g_center_snp_start<<" to "<<g_center_snp_end<<endl;
    g_markers = window_size;
    init_window();
    bool converged = false;
    gi_iteration = 0;
    int max_iter = 10;
    do{
      if (gi_iteration==0) compute_penetrance(); 
      compute_haplotype_weights();
      converged = check_mm_converged();// computes haplotype frequencies
      init_iteration_buffers();  
      cerr<<"MM iteration "<<gi_iteration<<" and convergence "<<converged<<endl;
      ++gi_iteration;
    }while(gi_iteration<max_iter && !converged);
    if (converged) cerr<<"Converged at iteration: "<<gi_iteration<<endl;
    else cerr<<"Frequencies did not converge\n";
    if (g_center_snp_start>0 || g_center_snp_end+g_flanking_snps==g_right_marker)
    impute_genotypes();
    finalize_window();
    cerr<<"Finalized window\n";

  }
  cerr<<"Done!\n";

}

void MendelGPU::init(){
  allocate_memory();
  if(run_gpu){
  // initialize the GPU if necessary
#ifdef USE_GPU
    init_opencl();
#endif
  }
}

void MendelGPU::allocate_memory(){
  cerr<<"initialize all the variables in base class\n";
  gf_epsilon = 1e-10;
  gf_convergence_criterion = 1e-4;
  gf_logpen_threshold = -10;
  debug_opencl = false;
  // this will provide the base path for the kernels
  imputation_software = getenv("IMPUTATION_SOFTWARE");
  cerr<<"Software base path is at "<<imputation_software<<endl;
  run_gpu = config->use_gpu;
  run_cpu = config->use_cpu;
  g_flanking_snps = config->flanking_snps;
  g_max_window = get_max_window_size();
  cerr<<"max window size is "<<g_max_window<<endl;
  g_genotype_imputation = config->model==1?1:0;
  g_likelihood_mode = config->g_likelihood_mode;
  g_delta = config->delta;
  geno_dim = config->is_phased?4:3;
  g_max_haplotypes = config->max_haplotypes;
  // load data sets
  load_datasets();
  // allocation begins here
  g_haplotypes = 0;
  g_haplotype = new int[g_max_haplotypes*g_max_window];  
  g_frequency = new float[g_max_haplotypes];
  for(int i=0;i<g_max_haplotypes;++i) g_frequency[i] = 0;
  g_old_frequency = new float[g_max_haplotypes];
  g_active_haplotype = new int[g_max_haplotypes];
  haploid_arr = new int[g_people];
  for(int i=0;i<g_people;++i)haploid_arr[i] = 0;
  penetrance_matrix_size = 0;
  if (geno_dim==PHASED_INPUT){ // for phased input
    cerr<<"Allocating penetrances for phased input\n";
    penetrance_matrix_size = g_people*2*g_max_haplotypes;
    logpenetrance_cache = new float[penetrance_matrix_size];
    penetrance_cache = new float[penetrance_matrix_size];
    frequency_cache = new float[penetrance_matrix_size];
  }else if (geno_dim==UNPHASED_INPUT){
    penetrance_matrix_size = g_max_haplotypes*g_max_haplotypes;
    cerr<<"Allocating penetrances for unphased input for "<<g_people<<" subjects and matrix size "<<penetrance_matrix_size<<"\n";
    logpenetrance_cache = new float[g_people*penetrance_matrix_size];
    penetrance_cache = new float[g_people*penetrance_matrix_size];
    frequency_cache = new float[penetrance_matrix_size];
  }else{
    cerr<<"Invalid dim of "<<geno_dim<<endl;
    exit(1);
  }
  g_current_weight = new float[g_max_haplotypes];
  g_weight = new float[g_max_haplotypes];
  subject_haplotype_weight = new float[g_people*g_max_haplotypes];
  max_penetrance = new float[g_people];
  // initialize variables
  for(int j=0;j<g_max_haplotypes;++j){
    g_frequency[j] = 1./g_haplotypes;
  }
  cerr<<"Allocated base class memory\n";
  return;
}

MendelGPU::~MendelGPU(){
  cerr<<"Entering destructor MendelGPU\n";
  delete[] g_frequency ;
  delete[] g_old_frequency ;
  delete[] g_active_haplotype ;
  delete[] logpenetrance_cache ;
  delete[] penetrance_cache ;
  delete[] frequency_cache ;
  delete[] g_current_weight ;
  delete[] g_weight ;
  delete[] subject_haplotype_weight ;
  delete[] max_penetrance ;
  cerr<<"Exiting destructor MendelGPU\n";
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

void MendelGPU::unmarshall(ifstream & ifs_cache, float * outputvec,int len){
  char * charvec = reinterpret_cast<char * >(outputvec);
  ifs_cache.read(charvec,sizeof(float) * len);
  outputvec = reinterpret_cast<float * > (charvec);
}

void MendelGPU::finalize_window(){
}

void MendelGPU::init_window(){
  for(int j=0;j<g_max_haplotypes;++j) g_old_frequency[j] = g_frequency[j];
}

void MendelGPU::init_iteration_buffers(){
  if (run_gpu){
    init_iteration_buffers_opencl();
  }
}

float MendelGPU::addLog(float log1,float log2){
  float smaller,larger,scale;
  if (log1<log2){
    smaller = log1;
    larger = log2;
  }else{
    smaller = log2;
    larger = log1;
  }
  //cerr<<"Smaller,larger: "<<smaller<<","<<larger<<endl;
  scale = smaller+(larger-smaller)/2.;
  //cerr<<"scale: "<<scale<<endl;
  float sum = exp(log1-scale)+exp(log2-scale); 
  //cerr<<"sum: "<<sum<<endl;
  //cerr<<"logsum: "<<log(sum)+scale<<endl;
  return log(sum)+scale;
}

float MendelGPU::compute_rsq(float * dosages,int stride,int index){
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
  //MendelGPU mendelGPU(NULL);
  //mendelGPU.init_gpu_(&platform_id,&device_id);
  int scaler = 2;
  int return_vec[256];
  for(int i=0;i<256;++i){
    return_vec[i] = i;
  }
  //mendelGPU.run_simple_(&scaler,return_vec);
  for(int i=0;i<256;++i){
    cout<<i<<":"<<return_vec[i]<<" ";
  }
  cout<<endl;
#endif
  return 0;
}

