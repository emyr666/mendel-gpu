using namespace std;
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<queue>
#include<tr1/unordered_set>
#include<tr1/unordered_map>
#include<set>
#include<map>
#include<list>
#include<math.h>
#include"../cl_constants.h"
#include<cstdlib>
#include<string.h>
#ifdef USE_GPU
#include<CL/cl.hpp>
#include"../clsafe.h"
#endif
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"

template<class T> cl::Buffer * MendelGPU::createBuffer(int rw, int dim, const char * label){
  cl::Buffer * buf = new cl::Buffer(*context, rw, sizeof(T) * dim, NULL, &err);
  clSafe(err,label);
  return buf;
}

void MendelGPU::init_opencl(){
  if(run_gpu){
  // initialize the GPU if necessary
#ifdef USE_GPU
    int platform_id = config->platform_id;
    int device_id = config->device_id;
    if (commandQueue!=NULL) return;
    cerr<<"Initializing GPU with platform id "<<platform_id<<
    " and device id "<<device_id<<".\n";
    vector<cl::Platform> platforms;
    err = cl::Platform::get(&platforms);
    cerr<<"Platform ID "<<platform_id<<" has name "<<
    platforms[platform_id].getInfo<CL_PLATFORM_NAME>().c_str()<<endl;
    cl_context_properties cps[3] = {CL_CONTEXT_PLATFORM,
    (cl_context_properties)(platforms[platform_id])(),0};
    context = new cl::Context(CL_DEVICE_TYPE_GPU,cps);
    devices = context->getInfo<CL_CONTEXT_DEVICES>();
    cerr<<"There are "<<devices.size()<<" devices\n";
    cerr<<"Device ID "<<device_id<<" is a ";
    cl_device_type dtype = devices[device_id].getInfo<CL_DEVICE_TYPE>();
    switch(dtype){
      case CL_DEVICE_TYPE_GPU:
        cerr<<"GPU\n";
        break;
      case CL_DEVICE_TYPE_CPU:
        cerr<<"CPU\n";
        break;
    } 
    commandQueue = new cl::CommandQueue(*context,devices[device_id],0,&err);
    string sources[]={"cl_constants.h","kernels/common.c","kernels/denovo.c",
    "kernels/guide.c"};
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
    kernel_compute_weights = new cl::Kernel(*program,"compute_weights",&err);
    clSafe(err,"creating kernel compute_weights");
    kernel_compute_weights_haploid = new cl::Kernel(*program,"compute_weights_haploid",&err);
    clSafe(err,"creating kernel compute_weights haploid");
    kernel_reduce_weights2 = new cl::Kernel(*program,"reduce_weights2",&err);
    clSafe(err,"creating kernel reduce_weights2");
    // CREATE BUFFERS
    buffer_simple_in = createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_simple");
    //buffer_simple_in = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating buffer in");
    buffer_simple_out = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * BLOCK_WIDTH, NULL, &err);
    clSafe(err,"creating buffer out");
    buffer_markers = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer markers");
    buffer_haplotypes = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer haplotypes");
    buffer_left_marker = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer left marker");
    buffer_haplotype = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * g_max_window * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype buffer");
    buffer_active_haplotype = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_max_haplotypes, NULL, &err);
    clSafe(err,"creating buffer active hap");
    buffer_subject_haplotype_weight = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating subject haplotype weights buffer");
    buffer_haplotype_weight = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype weights buffer");
    buffer_left_marker = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * 1, NULL, &err);
    clSafe(err,"creating buffer left marker");
    buffer_frequency = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_max_haplotypes, NULL, &err);
    clSafe(err,"creating haplotype frequency buffer");
    if (geno_dim==PHASED_INPUT){
      buffer_penetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*2*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer penetrance_cache");
      buffer_logpenetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*2*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer logpenetrance_cache");
      buffer_frequency_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_max_haplotypes, NULL, &err);
      clSafe(err,"creating buffer frequency_cache");
    }else if (geno_dim==UNPHASED_INPUT){
      buffer_penetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer penetrance_cache");
      buffer_logpenetrance_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer logpenetrance_cache");
      buffer_frequency_cache = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*penetrance_matrix_size, NULL, &err);
      clSafe(err,"creating buffer frequency_cache");
    }
    //buffer_region_snp_penetrance = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people * geno_dim * g_max_region_size , NULL, &err);
    //clSafe(err,"creating snp penetrance buffer");
    //buffer_region_snp_offset = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating snp offset buffer");
    buffer_genotype_imputation = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int), NULL, &err);
    clSafe(err,"creating buffer geno impute");
    buffer_max_penetrance = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people, NULL, &err);
    clSafe(err,"creating buffer max penetrance");
    buffer_haploid_arr = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*g_people, NULL, &err);
    clSafe(err,"creating buffer haploid_arr");
    buffer_iteration = new cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int)*1, NULL, &err);
    clSafe(err,"creating buffer iteration");
    cerr<<"GPU Buffers created\n";
    // SET KERNEL ARGUMENTS HERE
    int arg;
    int kernelWorkGroupSize;
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
    //err = kernel_compute_weights->setArg(arg++, *buffer_region_snp_penetrance);
    //clSafe(err,"set kernel arg for kernel_compute_weights");
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
    //err = kernel_compute_weights_haploid->setArg(arg++, *buffer_region_snp_penetrance);
    //clSafe(err,"set kernel arg for kernel_compute_weights_haploid");
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
    cerr<<"GPU kernel arguments assigned.\n";

    // TRANSFER ANY DATA HERE
    // SET KERNEL ARGUMENTS
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

