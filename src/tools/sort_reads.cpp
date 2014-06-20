/**
 * @file fq2fa.cpp
 * @brief 
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
#include "sequence/short_sequence.h"

using namespace std;

bool is_paired = false;
bool is_merged = false;
bool is_filtered = false;
int min_length = 0;
deque<ShortSequence> reads;
vector<int> aux;

bool Compare(int x, int y)
{
    if (reads[x] != reads[y])
        return reads[x] < reads[y];
    else
        return reads[x+1] < reads[y+1];
}

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    desc.AddOption("paired", "", is_paired, "if the reads are paired-end in one file");
    desc.AddOption("merge", "", is_merged, "if the reads are paired-end in two files");
    desc.AddOption("filter", "", is_filtered, "filter out reads containing 'N'");
    desc.AddOption("min_length", "", min_length, "minimum length ");

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

    ReadSequence(argv[1], reads);

    cout << "read" << endl;
    for (int i = 0; i < (int)reads.size(); i += 2)
    {
        if (reads[i+1] < reads[i])
            swap(reads[i], reads[i+1]);

        if ((int)reads[i].size() >= min_length && (int)reads[i+1].size() >= min_length)
            aux.push_back(i);
    }

    sort(aux.begin(), aux.end(), Compare);
    cout << "sort" << endl;

    int index = 0;
    int last = -1;
    FastaWriter writer(argv[2]);
    for (int i = 0; i < (int)aux.size(); ++i)
    {
        int id = aux[i];
        if (last == -1 || reads[id] != reads[last] || reads[id+1] != reads[last+1])
        {
            Sequence seq1(reads[id]);
            Sequence seq2(reads[id+1]);
            writer.Write(seq1, FormatString("reads_%d/1", index));
            writer.Write(seq2, FormatString("reads_%d/2", index));
            index += 2;
        }

        last = id;
    }

    cout << index << " " << reads.size();

    return 0;
}

