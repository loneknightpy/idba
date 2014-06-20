/**
 * @file split_fq.cpp
 * @brief Split a set of Fastq reads into two files.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.6
 * @date 2011-10-31
 */

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iostream>
#include <stdexcept>

#include "misc/options_description.h"
#include "misc/utils.h"
#include "sequence/sequence.h"
#include "sequence/sequence_io.h"

using namespace std;

bool is_paired = false;
bool is_merged = false;
bool is_filtered = false;

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    desc.AddOption("paired", "", is_paired, "if the reads are paired-end in one file");
    desc.AddOption("merge", "", is_merged, "if the reads are paired-end in two files");
    desc.AddOption("filter", "", is_filtered, "filter out reads containing 'N'");

    try
    {
        desc.Parse(argc, argv);

        if (argc < 4)
            throw logic_error("not enough parameters");

    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "fq2fa - Convert Fastq sequences to Fasta sequences." << endl;
        cerr << "Usage: fq2fa tmp.fq tmp.fa [...] " << endl;
        cerr << "       fq2fa --paired tmp.fq tmp.fa" << endl;
        cerr << "       fq2fa --merge tmp_1.fq tmp_2.fq tmp.fa" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    FastqReader reader(argv[1]);
    FastqWriter writer1(argv[2]);
    FastqWriter writer2(argv[3]);

    Sequence seq1, seq2;
    string comment1, comment2;
    string quality1, quality2;
    while (reader.Read(seq1, comment1, quality1) && reader.Read(seq2, comment2, quality2))
    {
        writer1.Write(seq1, comment1, quality1);
        writer2.Write(seq2, comment2, quality2);
    }

    return 0;
}

