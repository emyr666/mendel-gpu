include ../include.mk

OBJS=guide.o denovo.o
compile: $(OBJS)

denovo.o: denovo.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

guide.o: guide.cpp ../clsafe.h ../cl_constants.h
	$(C_COMPILER) $(CINC_FLAGS) $(PREPROC) -c $<

clean:
	rm -f $(OBJS)
