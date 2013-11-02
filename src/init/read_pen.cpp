#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"


ReadPenetrance::~ReadPenetrance(){
  cerr<<"Entering destructor Read penetrance\n";
  cerr<<"The maximum non-zero elements found was "<<max_nonzero_elements<<endl;
  delete parser;
  cerr<<"Parser deleted\n";
  delete[] read_alleles_mat;
  delete[] read_match_logmat;
  delete[] read_mismatch_logmat;
  delete[] mat_rows_by_subject;
  delete[] read_alleles_vec;
  delete[] read_match_logvec;
  delete[] read_mismatch_logvec;
  delete[] vector_offsets;
  delete[] haplotype_offsets;
  delete[] read_lengths;
  cerr<<"Exiting destructor Read penetrance\n";

}

ReadPenetrance::ReadPenetrance(MendelGPU * mendelgpu){
  this->mendelgpu = mendelgpu;
  if (mendelgpu->g_max_window>MAX_MATRIX_WIDTH){
    cerr<<"The reads parser constant MAX_MATRIX_WIDTH is currently set to "<<MAX_MATRIX_WIDTH<<" but Mendel GPU window size is "<<mendelgpu->g_max_window<<". Please ensure the MAX_MATRIX_WIDTH is at least as large as "<<mendelgpu->g_max_window<<".\n";
    throw "Error in init read penetrance";
  }
  cerr<<"Initializing variables for Read Penetrance\n";
  log_half = log(.5);
  max_nonzero_elements = 0;

  Config * config = mendelgpu->config;
  string variants = config->bimfile;
  //ReadParser::parse_positions(variants.data());
  ReadParser::init_phred_lookup();
  parser = new ReadParser();
  parser->init(config->bamfile.data(),variants.data());
  compact_rows = ReadParser::MAX_SUPERREADS;
  full_rows = ReadParser::MAX_SUPERREADS*ReadParser::MAX_SUBREADS;
  read_full_matrix_size = full_rows*mendelgpu->g_max_window;
  read_compact_matrix_size = compact_rows*mendelgpu->g_max_window;
  cerr<<"Read matrix sizes are "<<read_full_matrix_size<<" and "<<read_compact_matrix_size<<endl;
  read_alleles_mat = new int[mendelgpu->g_people*read_compact_matrix_size];
  read_match_logmat = new float[mendelgpu->g_people*read_compact_matrix_size];
  read_mismatch_logmat = new float[mendelgpu->g_people*read_compact_matrix_size];
  mat_rows_by_subject = new int[mendelgpu->g_people];

  read_alleles_vec = new int[mendelgpu->g_people*read_compact_matrix_size];
  read_match_logvec = new float[mendelgpu->g_people*read_compact_matrix_size];
  read_mismatch_logvec = new float[mendelgpu->g_people*read_compact_matrix_size];
  vector_offsets = new int[mendelgpu->g_people*compact_rows];
  haplotype_offsets = new int[mendelgpu->g_people*compact_rows];
  read_lengths = new int[mendelgpu->g_people*compact_rows];
}
