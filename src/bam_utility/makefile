# BELOW ARE EXAMPLES OF POTENTIAL HARDWARE PROFILES IF YOU PLAN TO DEPLOY
# THE PROGRAM ON MULTIPLE MACHINES
include locations.baseline

EXECUTABLES=compress filter_vcf get_consensus make_refseq_fasta make_subject_fasta parse_vcf postprocess_art_sam process_kgp_output process_macs_output read_penetrance transpose

all: $(EXECUTABLES)

compress: compress.cpp
	$(LINKER) -o $(BINDIR)/compress $(LINK_FLAGS)  $(CINC_FLAGS) $(PREPROC) $< $(LIBDIR)/libbam.a

read_penetrance: read_penetrance.cpp 
	$(LINKER) -o $(BINDIR)/read_penetrance $(LINK_FLAGS)  $(CINC_FLAGS) $(PREPROC) $< $(LIBDIR)/libbam.a

make_refseq_fasta: make_refseq_fasta.cpp
	$(LINKER) -o $(BINDIR)/make_refseq_fasta $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<

get_consensus: get_consensus.cpp
	$(LINKER) -o $(BINDIR)/get_consensus $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<

filter_vcf: filter_vcf.cpp
	$(LINKER) -o $(BINDIR)/filter_vcf $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<

parse_vcf: parse_vcf.cpp
	$(LINKER) -o $(BINDIR)/parse_vcf $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<

make_subject_fasta: make_subject_fasta.cpp
	$(LINKER) -o $(BINDIR)/make_subject_fasta $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<

postprocess_art_sam: postprocess_art_sam.cpp
	$(LINKER) -o $(BINDIR)/postprocess_art_sam $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<

process_macs_output: process_macs_output.cpp
	$(LINKER) -o $(BINDIR)/process_macs_output $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<

process_kgp_output: process_kgp_output.cpp
	$(LINKER) -o $(BINDIR)/process_kgp_output $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<

transpose: transpose.cpp
	$(LINKER) -o $(BINDIR)/transpose $(LINK_FLAGS) $(CINC_FLAGS) $(PREPROC) $<
	

clean:
	rm -f $(EXECUTABLES)
