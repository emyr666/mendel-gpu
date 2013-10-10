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
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif


void MendelGPU::init_opencl(){
  cerr<<"Initializing OpenCL for MendelGPU\n";
  if(run_gpu){
  // initialize the GPU if necessary
#ifdef USE_GPU
    int platform_id = config->platform_id;
    int device_id = config->device_id;
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
    vector<string> sources;
    sources.push_back("cl_constants.h");
    sources.push_back("kernels/common.c");
    sources.push_back("kernels/denovo.c");
    sources.push_back("kernels/denovo_glf.c");
    sources.push_back("kernels/denovo_reads.c");
    sources.push_back("kernels/guide.c");
    ostringstream full_source_str;
    for(int j=0;j<sources.size();++j){
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
    const char * preproc;
    if (geno_dim==PHASED_INPUT){
      preproc = "-Dphased";
    } else if(geno_dim==UNPHASED_INPUT){
      preproc = "-Dunphased";
    }
    err = program->build(devices,preproc);
    if(err!=CL_SUCCESS){
      cerr<<"Build failed:\n";
      string buffer;
      program->getBuildInfo(devices[0],CL_PROGRAM_BUILD_LOG,&buffer);
      cerr<<buffer<<endl;
      throw "Aborted from OpenCL build fail.";
    }
    // CREATE KERNELS
    createKernel("simple",kernel_simple);
    createKernel("compute_weights",kernel_compute_weights);
    createKernel("reduce_weights2",kernel_reduce_weights2);
    cerr<<"Kernels created\n";
    // CREATE BUFFERS
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_simple",buffer_simple_in);
    createBuffer<int>(CL_MEM_READ_WRITE,BLOCK_WIDTH,"buffer_simple_out",buffer_simple_out);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_iteration",buffer_iteration);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_genotype_imputation",buffer_genotype_imputation);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_markers",buffer_markers);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_haplotypes",buffer_haplotypes);
    createBuffer<int>(CL_MEM_READ_ONLY,1,"buffer_left_marker",buffer_left_marker);
    createBuffer<int>(CL_MEM_READ_ONLY,g_people,"buffer_haploid_arr",buffer_haploid_arr);
    createBuffer<int>(CL_MEM_READ_WRITE,g_max_window * g_max_haplotypes,"buffer_haplotype",buffer_haplotype);
    createBuffer<float>(CL_MEM_READ_WRITE,g_max_haplotypes,"buffer_frequency",buffer_frequency);
    if (geno_dim==PHASED_INPUT){
      createBuffer<float>(CL_MEM_READ_WRITE,g_people*2*g_max_haplotypes,"buffer_penetrance_cache",buffer_penetrance_cache);
      createBuffer<float>(CL_MEM_READ_WRITE,g_people*2*g_max_haplotypes,"buffer_logpenetrance_cache",buffer_logpenetrance_cache);
    }else if (geno_dim==UNPHASED_INPUT){
      createBuffer<float>(CL_MEM_READ_WRITE,g_people*penetrance_matrix_size,"buffer_penetrance_cache",buffer_penetrance_cache);
      createBuffer<float>(CL_MEM_READ_WRITE,g_people*penetrance_matrix_size,"buffer_logpenetrance_cache",buffer_logpenetrance_cache);
    }
    //createBuffer<>(CL_MEM_READ_WRITE,,"buffer_",buffer_);
    createBuffer<float>(CL_MEM_READ_WRITE,g_people * g_max_haplotypes,"buffer_subject_haplotype_weight",buffer_subject_haplotype_weight);
    createBuffer<int>(CL_MEM_READ_ONLY,g_max_haplotypes,"buffer_active_haplotype",buffer_active_haplotype);
    createBuffer<float>(CL_MEM_READ_WRITE,g_max_haplotypes,"buffer_haplotype_weight",buffer_haplotype_weight);
    //buffer_region_snp_penetrance = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float) * g_people * geno_dim * g_max_region_size , NULL, &err);
    //clSafe(err,"creating snp penetrance buffer");
    //buffer_region_snp_offset = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * 1, NULL, &err);
    //clSafe(err,"creating snp offset buffer");
    //buffer_max_penetrance = new cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*g_people, NULL, &err);
    //clSafe(err,"creating buffer max penetrance");
    cerr<<"GPU Buffers created\n";
    // initialize anything here
    writeToBuffer(buffer_genotype_imputation,1,&g_genotype_imputation,"buffer_genotype_imputation");
    writeToBuffer(buffer_haploid_arr,g_people,haploid_arr,"buffer_haploid_arr");
    cerr<<"GPU Buffers initialized\n";
    // SET KERNEL ARGUMENTS HERE
    int arg;
    int kernelWorkGroupSize;
    arg = 0;
    setArg(kernel_simple,arg,*buffer_simple_in,"kernel_simple");
    setArg(kernel_simple,arg,*buffer_simple_out,"kernel_simple");
    cerr<<"Set arguments for simple\n";
    //kernelWorkGroupSize = kernel_simple->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    //clSafe(err,"get workgroup size kernel simple");
    // COMPUTE SUBJECT SPECIFIC HAPLOTYPE WEIGHTS
    arg = 0;
    setArg(kernel_compute_weights,arg,g_max_haplotypes,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,g_max_window,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,penetrance_matrix_size,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_iteration,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_genotype_imputation,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_markers,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_haplotypes,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_left_marker,"kernel_compute_weights");
    if(geno_dim==UNPHASED_INPUT){
      setArg(kernel_compute_weights,arg,*buffer_haploid_arr,"kernel_compute_weights");
    }
    setArg(kernel_compute_weights,arg,*buffer_haplotype,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_frequency,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_penetrance_cache,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_subject_haplotype_weight,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,*buffer_active_haplotype,"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,cl::__local(sizeof(float)*g_max_haplotypes),"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,cl::__local(sizeof(float)*g_max_haplotypes),"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,cl::__local(sizeof(float)*g_max_haplotypes),"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,cl::__local(sizeof(float)*BLOCK_WIDTH),"kernel_compute_weights");
    setArg(kernel_compute_weights,arg,cl::__local(sizeof(float)*1),"kernel_compute_weights");
    cerr<<"Total args for kernel_compute_weights: "<<arg<<endl;
    kernelWorkGroupSize = kernel_compute_weights->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel compute_weights");
    cerr<<"compute_weights kernel work group size is "<<kernelWorkGroupSize<<endl;
    // COMPUTE HAPLOTYPE WEIGHTS BY REDUCTION
    arg = 0;
    setArg(kernel_reduce_weights2,arg,g_people,"kernel_reduce_weights2");
    setArg(kernel_reduce_weights2,arg,g_max_haplotypes,"kernel_reduce_weights2");
    setArg(kernel_reduce_weights2,arg,*buffer_haplotypes,"kernel_reduce_weights2");
    setArg(kernel_reduce_weights2,arg,*buffer_subject_haplotype_weight,"kernel_reduce_weights2");
    setArg(kernel_reduce_weights2,arg,*buffer_haplotype_weight,"kernel_reduce_weights2");
    setArg(kernel_reduce_weights2,arg,*buffer_active_haplotype,"kernel_reduce_weights2");
    setArg(kernel_reduce_weights2,arg,cl::__local(sizeof(int)*g_max_haplotypes),"kernel_reduce_weights2");
    setArg(kernel_reduce_weights2,arg,cl::__local(sizeof(float)*g_max_haplotypes),"kernel_reduce_weights2");
    kernelWorkGroupSize = kernel_reduce_weights2->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(devices[0], &err);
    clSafe(err,"get workgroup size kernel reduce_weights2");
    cerr<<"reduce_weights2 kernel work group size is "<<kernelWorkGroupSize<<endl;
    cerr<<"GPU kernel arguments assigned.\n";
#endif
  }
}

void MendelGPU::init_window_opencl(){
  if (run_gpu){
    #ifdef USE_GPU
    writeToBuffer(buffer_active_haplotype, g_max_haplotypes, g_active_haplotype,"buffer_active_haplotype");
    writeToBuffer(buffer_haplotypes, 1,&g_haplotypes, "buffer_haplotypes" );
    writeToBuffer(buffer_markers, 1, &g_markers, "buffer_markers" );
    writeToBuffer(buffer_left_marker, 1, &g_left_marker, "buffer_left_marker");
    writeToBuffer(buffer_haplotype, g_max_window*g_max_haplotypes, g_haplotype, "buffer_haplotype");
    writeToBuffer(buffer_frequency, g_max_haplotypes, g_frequency, "buffer_frequency");
    #endif
    cerr<<"Buffers sent to GPU for current window\n";
  }
}

void MendelGPU::init_iteration_buffers_opencl(){
#ifdef USE_GPU
  writeToBuffer(buffer_frequency, g_max_haplotypes, g_frequency, "buffer_frequency");
#endif
  cerr<<"Iteration Buffers sent to GPU\n";
}
