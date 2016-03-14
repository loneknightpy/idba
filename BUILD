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
    name = "idba",
    srcs = ["src/release/idba.cpp"],
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

