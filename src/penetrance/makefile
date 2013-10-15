include ../include.mk

OBJS=common.o common_opencl.o guide.o guide_opencl.o denovo_glf.o denovo_glf_opencl.o denovo_reads.o denovo_reads_opencl.o read_penetrance.o
compile: $(OBJS)

read_penetrance: read_penetrance.o
	$(LINKER) -o read_penetrance $(LINK_FLAGS) $< ../../lib/libbam.a

read_penetrance.o: read_penetrance.cpp ../read_penetrance.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

common.o: common.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

common_opencl.o: common_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_glf.o: denovo_glf.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_glf_opencl.o: denovo_glf_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_reads.o: denovo_reads.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_reads_opencl.o: denovo_reads_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide.o: guide.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_opencl.o: guide_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

clean:
	rm -f $(OBJS)
