#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_penetrance.hpp"

DenovoReadsMendelGPU::DenovoReadsMendelGPU(IO_manager * io):DenovoMendelGPU(io){
}


void DenovoReadsMendelGPU::load_datasets(){
  io_manager->read_input(this->g_people,this->g_snps);
  cerr<<"Loaded input in DenovoReadsMendelGPU\n";
}


DenovoReadsMendelGPU::~DenovoReadsMendelGPU(){
  cerr<<"Entering destructor denovo reads haplotyper\n";
  delete[] parsers;
  cerr<<"Parsers deleted\n";
  cerr<<"Exiting destructor denovo reads haplotyper\n";

}

void DenovoReadsMendelGPU::allocate_memory(){
  DenovoMendelGPU::allocate_memory();
  cerr<<"Initializing variables for denovo reads haplotyper\n";
  cerr<<"Analysis on "<<g_snps<<" SNPs and "<<g_people<<" persons\n";
  log_half = log(.5);

  string variants = config->bimfile;
  ReadParser::parse_positions(variants.data());
  ReadParser::init_phred_lookup(g_max_window);
  parsers = new ReadParser[g_people];
  for(int i=0;i<g_people;++i){
    string bam = config->bamfilevec[i];
    //cerr<<i<<": bam file: "<<bam<<endl;
    parsers[i].init(bam.data());
    //parsers[i].extract_region("1",0,1);
  }
  compact_rows = ReadParser::MAX_SUPERREADS;
  full_rows = ReadParser::MAX_SUPERREADS*ReadParser::MAX_SUBREADS;
  read_full_matrix_size = full_rows*g_max_window;
  read_compact_matrix_size = compact_rows*g_max_window;
  cerr<<"Read matrix sizes are "<<read_full_matrix_size<<" and "<<read_compact_matrix_size<<endl;
  superread_indices = new int[g_people*full_rows];
  read_alleles_mat = new int[g_people*read_compact_matrix_size];
  read_match_logmat = new float[g_people*read_full_matrix_size];
  read_mismatch_logmat = new float[g_people*read_full_matrix_size];
  mat_rows_by_subject = new int[g_people];
  cerr<<"Initialized variables for denovo reads haplotyper\n";
}

void DenovoReadsMendelGPU::init_window(){
  DenovoMendelGPU::init_window();
  // Load the reads for this current window
  if (run_gpu){
    init_window_opencl();
  }
  return;
}
