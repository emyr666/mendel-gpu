include ../include.mk

OBJS=common.o common_opencl.o guide.o guide_opencl.o denovo.o denovo_opencl.o denovo_glf.o denovo_glf_opencl.o denovo_reads.o denovo_reads_opencl.o
compile: $(OBJS)

common.o: common.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

common_opencl.o: common_opencl.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide.o: guide.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_opencl.o: guide_opencl.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo.o: denovo.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_opencl.o: denovo_opencl.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_glf.o: denovo_glf.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_glf_opencl.o: denovo_glf_opencl.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_reads.o: denovo_reads.cpp ../read_penetrance.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_reads_opencl.o: denovo_reads_opencl.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

clean:
	rm -f $(OBJS)
