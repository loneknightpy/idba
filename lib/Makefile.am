AM_LDFLAGS = -fopenmp -pthread
AM_CXXFLAGS = -Wall -O3 -fopenmp -pthread
INCLUDES =  -I$(top_srcdir)/src -I$(top_srcdir)/gtest_src 

libassembly_a_SOURCES = \
	$(top_srcdir)/src/basic/hash.h \
	$(top_srcdir)/src/basic/bit_operation.h \
	$(top_srcdir)/src/basic/math.h \
	$(top_srcdir)/src/basic/pool.h \
	$(top_srcdir)/src/basic/math.cpp \
	$(top_srcdir)/src/basic/atomic_integer.h \
	$(top_srcdir)/src/basic/functional.h \
	$(top_srcdir)/src/basic/kmer.h \
	$(top_srcdir)/src/basic/histgram.h \
	$(top_srcdir)/src/container/hash_table.h \
	$(top_srcdir)/src/container/hash_set.h \
	$(top_srcdir)/src/container/hash_map.h \
	$(top_srcdir)/src/container/managed_list.h \
	$(top_srcdir)/src/misc/hash_aligner.h \
	$(top_srcdir)/src/misc/blat_record.cpp \
	$(top_srcdir)/src/misc/blat_record.h \
	$(top_srcdir)/src/misc/options_description.h \
	$(top_srcdir)/src/misc/utils.cpp \
	$(top_srcdir)/src/misc/utils.h \
	$(top_srcdir)/src/misc/options_description.cpp \
	$(top_srcdir)/src/misc/hash_aligner.cpp \
	$(top_srcdir)/src/misc/log.cpp \
	$(top_srcdir)/src/misc/log.h \
	$(top_srcdir)/src/sequence/sequence_reader.h \
	$(top_srcdir)/src/sequence/sequence.h \
	$(top_srcdir)/src/sequence/compact_sequence.cpp \
	$(top_srcdir)/src/sequence/compact_sequence.h \
	$(top_srcdir)/src/sequence/read_library.h \
	$(top_srcdir)/src/sequence/sequence_writer.h \
	$(top_srcdir)/src/sequence/sequence_io.h \
	$(top_srcdir)/src/sequence/sequence.cpp \
	$(top_srcdir)/src/sequence/short_sequence.cpp \
	$(top_srcdir)/src/sequence/sequence_io.cpp \
	$(top_srcdir)/src/sequence/short_sequence.h \
	$(top_srcdir)/src/sequence/sequence_writer.cpp \
	$(top_srcdir)/src/sequence/sequence_reader.cpp \
	$(top_srcdir)/src/graph/contig_graph_path.h \
	$(top_srcdir)/src/graph/contig_info.cpp \
	$(top_srcdir)/src/graph/contig_graph_branch_group.h \
	$(top_srcdir)/src/graph/bit_edges.h \
	$(top_srcdir)/src/graph/contig_graph.cpp \
	$(top_srcdir)/src/graph/contig_graph.h \
	$(top_srcdir)/src/graph/hash_graph.cpp \
	$(top_srcdir)/src/graph/scaffold_graph_vertex.h \
	$(top_srcdir)/src/graph/contig_graph_branch_group.cpp \
	$(top_srcdir)/src/graph/contig_builder.h \
	$(top_srcdir)/src/graph/hash_graph_branch_group.cpp \
	$(top_srcdir)/src/graph/scaffold_graph.h \
	$(top_srcdir)/src/graph/hash_graph_branch_group.h \
	$(top_srcdir)/src/graph/hash_graph_path.h \
	$(top_srcdir)/src/graph/vertex_status.h \
	$(top_srcdir)/src/graph/scaffold_graph_path.h \
	$(top_srcdir)/src/graph/hash_graph_vertex.h \
	$(top_srcdir)/src/graph/hash_graph.h \
	$(top_srcdir)/src/graph/scaffold_graph.cpp \
	$(top_srcdir)/src/graph/contig_graph_vertex.h \
	$(top_srcdir)/src/graph/contig_info.h \
	$(top_srcdir)/src/assembly/assembly_utility.h \
	$(top_srcdir)/src/assembly/local_assembler.cpp \
	$(top_srcdir)/src/assembly/assembly_utility.cpp \
	$(top_srcdir)/src/assembly/local_assembler.h 

noinst_LIBRARIES = \
	libassembly.a 

