#! /usr/bin/python
##
# @file bazel-scan.py
# @brief Scan src directory and create BUILD file for bazel
# @author Yu Peng (ypeng@cs.hku.hk)
# @version 1.0.0
# @date 2016-03-13

import os
import re

COMMON = """
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

"""

TEMPLATE = """
cc_binary(
  name = "%s",
  srcs = ["%s"],
  copts = ["-Wall", "-O3", "-Isrc", "-fopenmp"],
  linkopts = ["-pthread", "-fopenmp", "-lm"],
  deps = [":assembly"],
)
"""

if __name__ == "__main__":
  fp = open('BUILD', 'w')
  fp.write(COMMON)
  paths = ["src/release", "src/tools"]
  for path in paths:
    for file in os.listdir(path):
      if re.search("\.cpp$", file):
        basename = os.path.basename(file).rpartition('.')[0]
        fp.write(TEMPLATE % (basename, path + "/" + file))




