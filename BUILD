
cc_library(
  name = "assembly",
  srcs = glob(
    [
      "src/basic/*.cpp",
      "src/container/*.cpp",
      "src/sequence/*.cpp",
      "src/graph/*.cpp",
      "src/misc/*.cpp",
      "src/assembly/*.cpp"
    ]
  ),
  hdrs = glob(
    [
      "src/basic/*.h",
      "src/container/*.h",
      "src/sequence/*.h",
      "src/graph/*.h",
      "src/misc/*.h",
      "src/assembly/*.h"
    ]
  ),
  includes = ["src"],
  copts = ["-Wall", "-O3", "-fopenmp"],
)

cc_test(
  name = "unit_test",
  srcs = glob(["src/test/*.cpp", "gtest_src/gtest/*.cc", "gtest_src/gtest/*.h"]),
  includes = ["gtest_src", "src"],
  copts = ["-Wall", "-O3", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)


cc_binary(
  name = "filterfa",
  srcs = ["src/release/filterfa.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "fq2fa",
  srcs = ["src/release/fq2fa.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "idba",
  srcs = ["src/release/idba.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "idba_hybrid",
  srcs = ["src/release/idba_hybrid.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "idba_tran",
  srcs = ["src/release/idba_tran.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "idba_ud",
  srcs = ["src/release/idba_ud.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "parallel_blat",
  srcs = ["src/release/parallel_blat.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "raw_n50",
  srcs = ["src/release/raw_n50.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "sim_reads",
  srcs = ["src/release/sim_reads.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "sim_reads_tran",
  srcs = ["src/release/sim_reads_tran.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "split_scaffold",
  srcs = ["src/release/split_scaffold.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "validate_contigs_blat",
  srcs = ["src/release/validate_contigs_blat.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "fa2fq",
  srcs = ["src/tools/fa2fq.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "filter_blat",
  srcs = ["src/tools/filter_blat.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "filter_contigs",
  srcs = ["src/tools/filter_contigs.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "idba_tran_test",
  srcs = ["src/tools/idba_tran_test.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "parallel_rna_blat",
  srcs = ["src/tools/parallel_rna_blat.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "print_graph",
  srcs = ["src/tools/print_graph.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "sample_reads",
  srcs = ["src/tools/sample_reads.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "scaffold",
  srcs = ["src/tools/scaffold.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "shuffle_reads",
  srcs = ["src/tools/shuffle_reads.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "sort_psl",
  srcs = ["src/tools/sort_psl.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "sort_reads",
  srcs = ["src/tools/sort_reads.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "split_fa",
  srcs = ["src/tools/split_fa.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "split_fq",
  srcs = ["src/tools/split_fq.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "test",
  srcs = ["src/tools/test.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "validate_component",
  srcs = ["src/tools/validate_component.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "validate_contigs_mummer",
  srcs = ["src/tools/validate_contigs_mummer.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "validate_reads_blat",
  srcs = ["src/tools/validate_reads_blat.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)

cc_binary(
  name = "validate_rna",
  srcs = ["src/tools/validate_rna.cpp"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)
