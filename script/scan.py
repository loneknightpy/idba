#! /usr/bin/python
##
# @file scan.py
# @brief Scan src directory recursively to build Makefile.am automatically.
# @author Yu Peng (ypeng@cs.hku.hk)
# @version 1.0.0
# @date 2011-08-02

import sys
import os
import re

params = {}
install = []
noninstall = []

def Usage():
    print("Usage: scan.py {bin, lib, test}\n")


def ParamsToString(params):
    return "".join(["%s = %s\n" % (k, v) for k, v in list(params.items())])

def Scan(paths, pattern):
    files = []
    for path in paths:
        for f in os.listdir(path):
            if re.search(pattern, f):
                files.append("$(top_srcdir)/" + path + "/" + f)
    return files

def ScanLibrary(name, paths, pattern):
    noninstall.append(name)
    files = Scan(paths, pattern)
    name = name.replace('.', '_');
    return " \\\n\t".join([name + "_SOURCES ="] + files)

def ScanBinary(paths, pattern):
    files = Scan(paths, pattern)
    sources = []
    for f in files:
        name = os.path.basename(f).rpartition('.')[0]
        name = name.replace('.', '_')
        name = name.replace('-', '_')
        noninstall.append(name)
        sources.append(name + "_SOURCES = " + f)
    return "\n".join(sources);

def SetInstall(name):
    noninstall.remove(name)
    install.append(name)

if __name__ == "__main__":
    if (len(sys.argv) < 2):
        Usage()
        exit(1)

    params["AM_CXXFLAGS"] = "-Wall -O3 -fopenmp -pthread"
    params["AM_LDFLAGS"] = "-fopenmp -pthread"
    params["INCLUDES"] = " " \
            + "-I$(top_srcdir)/src " \
            + "-I$(top_srcdir)/gtest_src "
    
    if (sys.argv[1] == "lib"):
        print(ParamsToString(params))
        #print "noinst_LIBRARIES = libcommon.a libassembly.a\n"

        #print ScanLibrary("libcommon.a", ["include/common", "src/common"], "\.cpp$|\.h$"), "\n"
        print(ScanLibrary("libassembly.a", 
                [#"src/common", "src/common", 
                    "src/basic", 
                    "src/container", 
                    "src/misc", 
                    "src/sequence", 
                    "src/graph", 
                    "src/assembly", ], 
                "\.cpp$|\.h$"), "\n")

        print("noinst_LIBRARIES = \\\n\t" + " \\\n\t".join(noninstall), "\n")

    elif (sys.argv[1] == "lib-test"):
        print(ParamsToString(params))
        #print "noinst_LIBRARIES = libcommon.a libassembly.a\n"

        #print ScanLibrary("libcommon.a", ["include/common", "src/common"], "\.cpp$|\.h$"), "\n"
        print(ScanLibrary("libassembly.a", 
                [#"src/common", "src/common", 
                    "src/basic", 
                    "src/container", 
                    "src/misc", 
                    "src/sequence", 
                    "src/graph", 
                    "src/assembly", ], 
                "\.cpp$|\.h$"), "\n")

        print(ScanLibrary("libgtest.a", 
                [ "gtest_src/gtest" ],
                "\.cpp$|\.h$|.cc$"), "\n")

        print("noinst_LIBRARIES = \\\n\t" + " \\\n\t".join(noninstall), "\n")

    elif (sys.argv[1] == "bin"):
        params["LIBS"] = "$(top_srcdir)/lib/libassembly.a @LIBS@";
        print(ParamsToString(params))
        print(ScanBinary(["src/tools"], "\.cpp$"), "\n")
        print(ScanBinary(["src/release"], "\.cpp$"), "\n")

        SetInstall("idba_hybrid")
        print("bin_PROGRAMS = \\\n\t" + " \\\n\t".join(install), "\n")
        print("noinst_PROGRAMS = \\\n\t" + " \\\n\t".join(noninstall), "\n")

    elif (sys.argv[1] == "release"):
        params["LIBS"] = "$(top_srcdir)/lib/libassembly.a @LIBS@";
        print(ParamsToString(params))
        print(ScanBinary(["src/release"], "\.cpp$"), "\n")

        SetInstall("idba_hybrid")
        print("bin_PROGRAMS = \\\n\t" + " \\\n\t".join(install), "\n")
        print("noinst_PROGRAMS = \\\n\t" + " \\\n\t".join(noninstall), "\n")

    elif (sys.argv[1] == "test"):
        params["LIBS"] = "$(top_srcdir)/lib/libassembly.a $(top_srcdir)/lib/libgtest.a @LIBS@";
        print(ParamsToString(params))
        print(ScanBinary(["src/test"], "\.cpp$"), "\n")

        print("noinst_PROGRAMS = \\\n\t" + " \\\n\t".join(noninstall), "\n")




