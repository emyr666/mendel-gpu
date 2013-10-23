include ../include.mk

OBJS=common.o read_pen.o read_pen_opencl.o common_opencl.o guide_glf.o guide_glf_opencl.o guide_reads.o denovo_glf.o denovo_glf_opencl.o denovo_reads.o read_parser.o
compile: $(OBJS)

read_parser: read_parser.o
	$(LINKER) -o read_penetrance $(LINK_FLAGS) $< ../../lib/libbam.a

read_parser.o: read_parser.cpp ../read_utilities.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

read_pen.o: read_pen.cpp ../mendel_gpu.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

read_pen_opencl.o: read_pen_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

common.o: common.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

common_opencl.o: common_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_glf.o: denovo_glf.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_reads.o: denovo_reads.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_glf.o: guide_glf.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_glf_opencl.o: guide_glf_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_reads.o: guide_reads.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<


clean:
	rm -f $(OBJS)
