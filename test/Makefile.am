AM_LDFLAGS = -fopenmp -pthread
LIBS = $(top_srcdir)/lib/libassembly.a $(top_srcdir)/lib/libgtest.a
AM_CXXFLAGS = -Wall -O3 -fopenmp -pthread
INCLUDES =  -I$(top_srcdir)/src -I$(top_srcdir)/gtest_src 

compact_sequence_unittest_SOURCES = $(top_srcdir)/src/test/compact_sequence_unittest.cpp
hash_table_unittest_SOURCES = $(top_srcdir)/src/test/hash_table_unittest.cpp
bit_edges_unittest_SOURCES = $(top_srcdir)/src/test/bit_edges_unittest.cpp
kmer_unittest_SOURCES = $(top_srcdir)/src/test/kmer_unittest.cpp
hash_set_unittest_SOURCES = $(top_srcdir)/src/test/hash_set_unittest.cpp
atomic_integer_unittest_SOURCES = $(top_srcdir)/src/test/atomic_integer_unittest.cpp
vertex_status_unittest_SOURCES = $(top_srcdir)/src/test/vertex_status_unittest.cpp
histgram_unittest_SOURCES = $(top_srcdir)/src/test/histgram_unittest.cpp
managed_list_unittest_SOURCES = $(top_srcdir)/src/test/managed_list_unittest.cpp
sequence_unittest_SOURCES = $(top_srcdir)/src/test/sequence_unittest.cpp
short_sequence_unittest_SOURCES = $(top_srcdir)/src/test/short_sequence_unittest.cpp
hash_map_unittest_SOURCES = $(top_srcdir)/src/test/hash_map_unittest.cpp 

noinst_PROGRAMS = \
	compact_sequence_unittest \
	hash_table_unittest \
	bit_edges_unittest \
	kmer_unittest \
	hash_set_unittest \
	atomic_integer_unittest \
	vertex_status_unittest \
	histgram_unittest \
	managed_list_unittest \
	sequence_unittest \
	short_sequence_unittest \
	hash_map_unittest 

