#! /usr/bin/python
##
# @file run-unittest.py
# @brief Scan test directory and run all unit tests.
# @author Yu Peng (ypeng@cs.hku.hk)
# @version 1.0.0
# @date 2011-08-02


import sys
import os
import re

if __name__ == "__main__":

    path = "test/";
    #test_programs = []
    for file in os.listdir(path):
        if (file[0] != '.' and re.search("unittest", file) and os.access(path + file, os.X_OK)):
            os.system(path + file);
            #test_programs.append(path + file);




