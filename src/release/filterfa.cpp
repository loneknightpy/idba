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
#include "sequence/sequence_io.h"
#include "misc/utils.h"
#include "sequence/sequence.h"

using namespace std;

bool is_paired = false;
bool is_merged = false;

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    desc.AddOption("paired", "", is_paired, "if the reads are paired-end");
    desc.AddOption("merge", "", is_merged, "if the reads are paired-end in two files");

    try
    {
        desc.Parse(argc, argv);

        if (argc < 2)
            throw logic_error("not enough parameters");

    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "fq2fa - Filter out fasta sequence containing N." << endl;
        cerr << "Usage: filterfa tmp.fa out.fa " << endl;
        cerr << "       filterfa --paired tmp.fa out.fa" << endl;
        cerr << "       filterfa --merged tmp_1.fa tmp_2.fa out.fa" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    if (!is_paired && !is_merged)
    {
        FastaReader reader(argv[1]);
        FastaWriter writer(argv[2]);

        Sequence seq;
        string comment;
        while (reader.Read(seq, comment))
        {
            if (seq.IsValid())
            {
                writer.Write(seq, comment);
            }
        }
    }
    else if (is_merged)
    {
        FastaReader reader1(argv[1]);
        FastaReader reader2(argv[2]);
        FastaWriter writer(argv[3]);

        Sequence seq1, seq2;
        string comment1, comment2;
        while (reader1.Read(seq1, comment1) && reader2.Read(seq2, comment2))
        {
            if (seq1.IsValid() && seq2.IsValid())
            {
                writer.Write(seq1, comment1);
                writer.Write(seq2, comment2);
            }
        }
    }
    else if (is_paired)
    {
        FastaReader reader1(argv[1]);
        FastaWriter writer(argv[2]);

        Sequence seq1, seq2;
        string comment1, comment2;
        while (reader1.Read(seq1, comment1) && reader1.Read(seq2, comment2))
        {
            if (seq1.IsValid() && seq2.IsValid())
            {
                writer.Write(seq1, comment1);
                writer.Write(seq2, comment2);
            }
        }
    }

    return 0;
}

