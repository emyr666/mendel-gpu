#include"io_manager.hpp"
#include"cl_constants.h"
#include"mendel_gpu.hpp"
#ifdef USE_GPU
#include"cl_templates.hpp"
#endif


int main(int argc,char * argv[]){
  if (argc<2){
    cerr<<"Usage: mendel_gpu <config xml filename>\n";
    return 1;
  }
  int arg = 0;
  const char * configfile = argv[++arg];
  const char * sourcepath = getenv("IMPUTATION_SOFTWARE");
  MendelGPU * haplotyper = NULL;
  try{
    IO_manager * io = new IO_manager();
    if (!io->load_config(configfile)){
      throw "Cannot load configuration file";
    }
    Config *  config = io->config;
    if (config->use_reference_panel){
      haplotyper = new GuidedMendelGPU(io);
    }else{
      if (config->g_likelihood_mode==Config::LIKELIHOOD_MODE_READS){
        cerr<<"Shot gun reads\n";
        haplotyper = new DenovoReadsMendelGPU(io);
      }else if (config->g_likelihood_mode==Config::LIKELIHOOD_MODE_GENOTYPES){
        cerr<<"GLF input\n";
        haplotyper = new DenovoGlfMendelGPU(io);
      }
    }
    haplotyper->init();
    haplotyper->run_sliding_window();
  }catch(const char * & mesg){
    cerr<<"Caught an exception with message "<<mesg<<endl;
    return 1;
  }catch(exception e){
    cerr<<"Caught an exception with message "<<e.what()<<endl;
    return 1;
  }
  delete haplotyper;
  return 0;
}
