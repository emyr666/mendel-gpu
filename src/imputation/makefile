include ../include.mk

OBJS=guide.o guide_opencl.o denovo.o denovo_opencl.o
compile: $(OBJS)

denovo.o: denovo.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_opencl.o: denovo_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide.o: guide.cpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide_opencl.o: guide_opencl.cpp ../cl_templates.hpp ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

clean:
	rm -f $(OBJS)
