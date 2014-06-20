/**
 * @file raw_n50.cpp
 * @brief Calculate the raw N50 of a set of contigs.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-10
 */

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iostream>
#include <stdexcept>

#include "misc/options_description.h"
#include "sequence/sequence_io.h"
#include "misc/utils.h"
#include "sequence/sequence.h"

using namespace std;

int min_contig = 100;
int ref_length = 0;

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    desc.AddOption("min_contig", "", min_contig, "");
    desc.AddOption("ref_length", "", ref_length, "");

    try
    {
        desc.Parse(argc, argv);

        if (argc < 2)
            throw logic_error("not enough parameters");

    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "raw_n50 - Compute the raw n50 of a set of contigs." << endl;
        cerr << "Usage: raw_n50 contigs.fa" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    deque<Sequence> contigs;
    ReadSequence(argv[1], contigs);
    PrintN50(contigs, min_contig, ref_length);

    return 0;
}
