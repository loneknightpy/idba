AM_LDFLAGS = -fopenmp -pthread -std=c++0x
LIBS = $(top_srcdir)/lib/libassembly.a @LIBS@
AM_CXXFLAGS = -Wall -O3 -fopenmp -pthread -std=c++0x
INCLUDES =  -I$(top_srcdir)/src -I$(top_srcdir)/gtest_src 

fa2fq_SOURCES = $(top_srcdir)/src/tools/fa2fq.cpp
filter_blat_SOURCES = $(top_srcdir)/src/tools/filter_blat.cpp
filter_contigs_SOURCES = $(top_srcdir)/src/tools/filter_contigs.cpp
idba_tran_test_SOURCES = $(top_srcdir)/src/tools/idba_tran_test.cpp
parallel_rna_blat_SOURCES = $(top_srcdir)/src/tools/parallel_rna_blat.cpp
print_graph_SOURCES = $(top_srcdir)/src/tools/print_graph.cpp
sample_reads_SOURCES = $(top_srcdir)/src/tools/sample_reads.cpp
scaffold_SOURCES = $(top_srcdir)/src/tools/scaffold.cpp
shuffle_reads_SOURCES = $(top_srcdir)/src/tools/shuffle_reads.cpp
sort_psl_SOURCES = $(top_srcdir)/src/tools/sort_psl.cpp
sort_reads_SOURCES = $(top_srcdir)/src/tools/sort_reads.cpp
split_fa_SOURCES = $(top_srcdir)/src/tools/split_fa.cpp
split_fq_SOURCES = $(top_srcdir)/src/tools/split_fq.cpp
test_SOURCES = $(top_srcdir)/src/tools/test.cpp
validate_component_SOURCES = $(top_srcdir)/src/tools/validate_component.cpp
validate_contigs_mummer_SOURCES = $(top_srcdir)/src/tools/validate_contigs_mummer.cpp
validate_reads_blat_SOURCES = $(top_srcdir)/src/tools/validate_reads_blat.cpp
validate_rna_SOURCES = $(top_srcdir)/src/tools/validate_rna.cpp 

filterfa_SOURCES = $(top_srcdir)/src/release/filterfa.cpp
fq2fa_SOURCES = $(top_srcdir)/src/release/fq2fa.cpp
idba_SOURCES = $(top_srcdir)/src/release/idba.cpp
idba_hybrid_SOURCES = $(top_srcdir)/src/release/idba_hybrid.cpp
idba_tran_SOURCES = $(top_srcdir)/src/release/idba_tran.cpp
idba_ud_SOURCES = $(top_srcdir)/src/release/idba_ud.cpp
parallel_blat_SOURCES = $(top_srcdir)/src/release/parallel_blat.cpp
raw_n50_SOURCES = $(top_srcdir)/src/release/raw_n50.cpp
sim_reads_SOURCES = $(top_srcdir)/src/release/sim_reads.cpp
sim_reads_tran_SOURCES = $(top_srcdir)/src/release/sim_reads_tran.cpp
split_scaffold_SOURCES = $(top_srcdir)/src/release/split_scaffold.cpp
validate_contigs_blat_SOURCES = $(top_srcdir)/src/release/validate_contigs_blat.cpp 

bin_PROGRAMS = \
	idba_hybrid 

noinst_PROGRAMS = \
	fa2fq \
	filter_blat \
	filter_contigs \
	idba_tran_test \
	parallel_rna_blat \
	print_graph \
	sample_reads \
	scaffold \
	shuffle_reads \
	sort_psl \
	sort_reads \
	split_fa \
	split_fq \
	test \
	validate_component \
	validate_contigs_mummer \
	validate_reads_blat \
	validate_rna \
	filterfa \
	fq2fa \
	idba \
	idba_tran \
	idba_ud \
	parallel_blat \
	raw_n50 \
	sim_reads \
	sim_reads_tran \
	split_scaffold \
	validate_contigs_blat 

