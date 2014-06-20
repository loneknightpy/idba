/**
 * @file split_scaffold.cpp
 * @brief Remove the N in scaffolds and split them into contigs.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.3
 * @date 2011-09-09
 */

#include <algorithm>
#include <cctype>
#include <cmath>
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

int main(int argc, char *argv[])
{
    OptionsDescription desc;

    try
    {
        desc.Parse(argc, argv);

        if (argc < 2)
            throw logic_error("not enough parameters");

    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "split_scaffold - Split the scaffols into contigs by removing 'N' " << endl;
        cerr << "Usage: split_scaffold scaffold.fa" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    FastaReader reader(argv[1]);
    FastaWriter writer(string(argv[1]) + ".split");
    Sequence seq;
    string name;
    while (reader.Read(seq, name))
    {
        int last = -1;
        int index = 0;
        for (unsigned j = 0; j < seq.size(); ++j)
        {
            if (seq[j] < 4 && last == -1)
                last = j;
            else if (seq[j] == 4)
            {
                if (last != -1)
                {
                    Sequence sub_seq(seq, last, j - last);
                    last = -1;
                    writer.Write(sub_seq, FormatString("%s_%d", name.c_str(), index++));
                }
            }
        }

        if (last != -1)
        {
            Sequence sub_seq(seq, last);
            writer.Write(sub_seq, FormatString("%s_%d", name.c_str(), index++));
        }
    }

    return 0;
}

