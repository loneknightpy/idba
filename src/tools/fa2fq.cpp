/**
 * @file fa2fq.cpp
 * @brief Convert fasta to fastq.
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

        if (argc < 3)
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

    FastaReader reader(argv[1]);
    FastqWriter writer(argv[2]);

    Sequence seq;
    string comment;
    while (reader.Read(seq, comment))
    {
        string quality;
        quality.append(seq.size(), 33 + 40);
        writer.Write(seq, comment, quality);
    }

    return 0;
}

