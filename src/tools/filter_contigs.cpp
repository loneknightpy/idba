/**
 * @file filter_contgis.cpp
 * @brief Filter out short contigs.
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

int min_contig = 100;

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    desc.AddOption("min_contig", "", min_contig, "filter out reads containing 'N'");

    try
    {
        desc.Parse(argc, argv);

        if (argc < 3)
            throw logic_error("not enough parameters");

    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "filter_contigs - Convert Fastq sequences to Fasta sequences." << endl;
        cerr << "Usage: fq2fa tmp.fq tmp.fa [...] " << endl;
        cerr << "       fq2fa --paired tmp.fq tmp.fa" << endl;
        cerr << "       fq2fa --merge tmp_1.fq tmp_2.fq tmp.fa" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    FastaReader reader(argv[1]);
    FastaWriter writer(argv[2]);

    Sequence seq;
    string comment;
    while (reader.Read(seq, comment))
    {
        if ((int)seq.size() >= min_contig)
        {
            writer.Write(seq, comment);
        }
    }

    return 0;
}

