
typedef struct {
public:
//VALUES THAT WILL BE WRITTEN BY get_read_data: 
  int n,m;// n dimension, i.e. number of reads in window, m dimension, i.e. end-start+1
  char *geno, // an n*m matrix of genotypes, '0' for missing, '1' for ref allele, '2' for non ref allele
  *qual; // an n*m matrix of phred quality scores
//VALUES USED FOR INTERNAL HANDLING:
  uint n_i, tid, start, end;
  char **qname;
} matrix_struct;


//VALUES TO WRITE TO BEFORE LAUNCHING get_read_data: uint tid, the tid, i.e. chromosome, 0-indexed.  start, start position, 1-indexed end,  end position, 1-indexed.  Because these fields in matrix need to be set each time anyway before get_read_data is called, I suggest integrating into one step, setting the matrix_struct fields within the body of this function so the function signature's intended behavior is more intuitive.

int get_read_data(char *bam_file, uint tid, uint start, uint end, matrix_struct * read_matrices);

// self-explanatory utility function.

inline float phred2prob(char phred);

// possible implementation in pseudo code
// #include "bam_interface.h"
//
// // if you MUST use global variables, at least qualify them to avoid
// // clashes with the rest of the software. this will work with g++
// 
// namespace alex_library{
//   matrix_struct * matrix = NULL;
// }
//
// using namespace alex_library;
//
// inline float phred2prob(char phred){
//   // do simple conversion
//   return myfloat;
// }
//
// void allocate_datastructures(){
//  // allocate matrix as appropriate
// }
//
// int get_read_data(char * bam_file, uint tid, uint start, uint end, matrix_struct * read_matrices){
//   
//   if (matrix == NULL) allocate_datastructures();   
//   matrix.tid = tid;
//   matrix.start = start;
//   matrix.end = end;
//   // fill in existing implementation of get_read_data
//   read_matrices = matrix;
//   if (exception_found) return error_code;
//   else return 0;
// }
