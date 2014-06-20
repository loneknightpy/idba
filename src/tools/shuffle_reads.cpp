/**
 * @file shuffle_reads.cpp
 * @brief Shuffle a set of reads.
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

int main(int argc, char *argv[])
{
    OptionsDescription desc;

    try
    {
        desc.Parse(argc, argv);

        if (argc < 3)
            throw logic_error("not enough parameters");

    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "shuffle_reads - shuffle reads." << endl;
        cerr << "Usage: shuffle_reads reads.fq shuffle_reads.fa [...] " << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    deque<Sequence> reads;
    deque<string> names;
    ReadSequence(argv[1], reads, names);

    int n = reads.size()/2;
    deque<int> aux(n);
    for (int i = 0; i < n; ++i)
        aux[i] = i;

    for (int i = 0; i < n; ++i)
        swap(aux[i], aux[rand()%(n-i) + i]);

    FastaWriter writer(argv[2]);
    for (int i = 0; i < n; ++i)
    {
        writer.Write(reads[aux[i]*2], names[aux[i]*2]);
        writer.Write(reads[aux[i]*2+1], names[aux[i]*2+1]);
    }

    return 0;
}

