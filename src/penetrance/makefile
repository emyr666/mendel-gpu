include ../include.mk

OBJS=common.o guide.o denovo_glf.o denovo_reads.o read_penetrance.o
compile: $(OBJS)

read_penetrance: read_penetrance.o
	$(LINKER) -o read_penetrance $(LINK_FLAGS) $< ../../lib/libbam.a

read_penetrance.o: read_penetrance.cpp ../read_penetrance.hpp
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

common.o: common.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_glf.o: denovo_glf.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

denovo_reads.o: denovo_reads.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide.o: guide.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

clean:
	rm -f $(OBJS)
