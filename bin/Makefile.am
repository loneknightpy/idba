AM_LDFLAGS = -fopenmp -pthread
LIBS = $(top_srcdir)/lib/libassembly.a
AM_CXXFLAGS = -Wall -O3 -fopenmp -pthread
INCLUDES =  -I$(top_srcdir)/src -I$(top_srcdir)/gtest_src 

shuffle_reads_SOURCES = $(top_srcdir)/src/tools/shuffle_reads.cpp
scaffold_SOURCES = $(top_srcdir)/src/tools/scaffold.cpp
sample_reads_SOURCES = $(top_srcdir)/src/tools/sample_reads.cpp
validate_contigs_mummer_SOURCES = $(top_srcdir)/src/tools/validate_contigs_mummer.cpp
sort_psl_SOURCES = $(top_srcdir)/src/tools/sort_psl.cpp
print_graph_SOURCES = $(top_srcdir)/src/tools/print_graph.cpp
validate_component_SOURCES = $(top_srcdir)/src/tools/validate_component.cpp
split_fa_SOURCES = $(top_srcdir)/src/tools/split_fa.cpp
filter_blat_SOURCES = $(top_srcdir)/src/tools/filter_blat.cpp
validate_rna_SOURCES = $(top_srcdir)/src/tools/validate_rna.cpp
fa2fq_SOURCES = $(top_srcdir)/src/tools/fa2fq.cpp
parallel_rna_blat_SOURCES = $(top_srcdir)/src/tools/parallel_rna_blat.cpp
validate_reads_blat_SOURCES = $(top_srcdir)/src/tools/validate_reads_blat.cpp
split_fq_SOURCES = $(top_srcdir)/src/tools/split_fq.cpp
sort_reads_SOURCES = $(top_srcdir)/src/tools/sort_reads.cpp
filter_contigs_SOURCES = $(top_srcdir)/src/tools/filter_contigs.cpp
idba_tran_test_SOURCES = $(top_srcdir)/src/tools/idba_tran_test.cpp
test_SOURCES = $(top_srcdir)/src/tools/test.cpp 

split_scaffold_SOURCES = $(top_srcdir)/src/release/split_scaffold.cpp
parallel_blat_SOURCES = $(top_srcdir)/src/release/parallel_blat.cpp
idba_tran_SOURCES = $(top_srcdir)/src/release/idba_tran.cpp
raw_n50_SOURCES = $(top_srcdir)/src/release/raw_n50.cpp
idba_ud_SOURCES = $(top_srcdir)/src/release/idba_ud.cpp
sim_reads_tran_SOURCES = $(top_srcdir)/src/release/sim_reads_tran.cpp
validate_contigs_blat_SOURCES = $(top_srcdir)/src/release/validate_contigs_blat.cpp
fq2fa_SOURCES = $(top_srcdir)/src/release/fq2fa.cpp
idba_SOURCES = $(top_srcdir)/src/release/idba.cpp
filterfa_SOURCES = $(top_srcdir)/src/release/filterfa.cpp
sim_reads_SOURCES = $(top_srcdir)/src/release/sim_reads.cpp
idba_hybrid_SOURCES = $(top_srcdir)/src/release/idba_hybrid.cpp 

bin_PROGRAMS = \
	idba_hybrid 

noinst_PROGRAMS = \
	shuffle_reads \
	scaffold \
	sample_reads \
	validate_contigs_mummer \
	sort_psl \
	print_graph \
	validate_component \
	split_fa \
	filter_blat \
	validate_rna \
	fa2fq \
	parallel_rna_blat \
	validate_reads_blat \
	split_fq \
	sort_reads \
	filter_contigs \
	idba_tran_test \
	test \
	split_scaffold \
	parallel_blat \
	idba_tran \
	raw_n50 \
	idba_ud \
	sim_reads_tran \
	validate_contigs_blat \
	fq2fa \
	idba \
	filterfa \
	sim_reads 

