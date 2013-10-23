include ../include.mk

OBJS=read_pen.o read_pen_opencl.o common.o common_opencl.o guide.o guide_opencl.o guide_glf.o guide_glf_opencl.o guide_reads.o guide_reads_opencl.o denovo.o denovo_opencl.o denovo_glf.o denovo_glf_opencl.o denovo_reads.o denovo_reads_opencl.o
compile: $(OBJS)

common.o: common.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

common_opencl.o: common_opencl.cpp ../cl_templates.hpp ../cl_constants.h ../cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

read_pen.o: read_pen.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

read_pen_opencl.o: read_pen_opencl.cpp ../cl_templates.hpp ../cl_constants.h ../cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide.o: guide.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_opencl.o: guide_opencl.cpp ../cl_templates.hpp ../cl_constants.h ../cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_glf.o: guide_glf.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_glf_opencl.o: guide_glf_opencl.cpp ../cl_templates.hpp ../cl_constants.h ../cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_reads.o: guide_reads.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_reads_opencl.o: guide_reads_opencl.cpp ../cl_templates.hpp ../cl_constants.h ../cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo.o: denovo.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_opencl.o: denovo_opencl.cpp ../cl_templates.hpp ../cl_constants.h ../cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_glf.o: denovo_glf.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_glf_opencl.o: denovo_glf_opencl.cpp ../cl_templates.hpp ../cl_constants.h ../cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_reads.o: denovo_reads.cpp ../read_utilities.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_reads_opencl.o: denovo_reads_opencl.cpp ../cl_templates.hpp ../cl_constants.h ../cl_templates.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

clean:
	rm -f $(OBJS)
