AM_LDFLAGS = -fopenmp -pthread
LIBS = $(top_srcdir)/lib/libassembly.a @LIBS@
AM_CXXFLAGS = -Wall -O3 -fopenmp -pthread
INCLUDES =  -I$(top_srcdir)/src -I$(top_srcdir)/gtest_src 

shuffle_reads_SOURCES = $(top_srcdir)/src/tools/shuffle_reads.cpp
scaffold_SOURCES = $(top_srcdir)/src/tools/scaffold.cpp
validate_component_SOURCES = $(top_srcdir)/src/tools/validate_component.cpp
fa2fq_SOURCES = $(top_srcdir)/src/tools/fa2fq.cpp
validate_rna_SOURCES = $(top_srcdir)/src/tools/validate_rna.cpp
validate_contigs_mummer_SOURCES = $(top_srcdir)/src/tools/validate_contigs_mummer.cpp
print_graph_SOURCES = $(top_srcdir)/src/tools/print_graph.cpp
idba_tran_test_SOURCES = $(top_srcdir)/src/tools/idba_tran_test.cpp
sample_reads_SOURCES = $(top_srcdir)/src/tools/sample_reads.cpp
sort_psl_SOURCES = $(top_srcdir)/src/tools/sort_psl.cpp
sort_reads_SOURCES = $(top_srcdir)/src/tools/sort_reads.cpp
filter_contigs_SOURCES = $(top_srcdir)/src/tools/filter_contigs.cpp
parallel_rna_blat_SOURCES = $(top_srcdir)/src/tools/parallel_rna_blat.cpp
filter_blat_SOURCES = $(top_srcdir)/src/tools/filter_blat.cpp
split_fa_SOURCES = $(top_srcdir)/src/tools/split_fa.cpp
test_SOURCES = $(top_srcdir)/src/tools/test.cpp
validate_reads_blat_SOURCES = $(top_srcdir)/src/tools/validate_reads_blat.cpp
split_fq_SOURCES = $(top_srcdir)/src/tools/split_fq.cpp 

fq2fa_SOURCES = $(top_srcdir)/src/release/fq2fa.cpp
validate_contigs_blat_SOURCES = $(top_srcdir)/src/release/validate_contigs_blat.cpp
idba_ud_SOURCES = $(top_srcdir)/src/release/idba_ud.cpp
raw_n50_SOURCES = $(top_srcdir)/src/release/raw_n50.cpp
idba_tran_SOURCES = $(top_srcdir)/src/release/idba_tran.cpp
idba_hybrid_SOURCES = $(top_srcdir)/src/release/idba_hybrid.cpp
parallel_blat_SOURCES = $(top_srcdir)/src/release/parallel_blat.cpp
sim_reads_tran_SOURCES = $(top_srcdir)/src/release/sim_reads_tran.cpp
filterfa_SOURCES = $(top_srcdir)/src/release/filterfa.cpp
idba_SOURCES = $(top_srcdir)/src/release/idba.cpp
split_scaffold_SOURCES = $(top_srcdir)/src/release/split_scaffold.cpp
sim_reads_SOURCES = $(top_srcdir)/src/release/sim_reads.cpp 

bin_PROGRAMS = \
	idba_hybrid 

noinst_PROGRAMS = \
	shuffle_reads \
	scaffold \
	validate_component \
	fa2fq \
	validate_rna \
	validate_contigs_mummer \
	print_graph \
	idba_tran_test \
	sample_reads \
	sort_psl \
	sort_reads \
	filter_contigs \
	parallel_rna_blat \
	filter_blat \
	split_fa \
	test \
	validate_reads_blat \
	split_fq \
	fq2fa \
	validate_contigs_blat \
	idba_ud \
	raw_n50 \
	idba_tran \
	parallel_blat \
	sim_reads_tran \
	filterfa \
	idba \
	split_scaffold \
	sim_reads 

