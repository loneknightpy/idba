/**
 * @file sort_psl.cpp
 * @brief Sort the alignment result of blat.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.3
 * @date 2011-09-06
 */

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "misc/blat_record.h"
#include "misc/options_description.h"
#include "misc/utils.h"
#include "sequence/sequence.h"

using namespace std;

const int MaxLine = (1 << 20);

char line[MaxLine];
char buf[MaxLine];
deque<BlatRecord> records;
deque<string> lines;
deque<int> aux;

bool Compare(int x, int y)
{
    BlatRecord &r1 = records[x];
    BlatRecord &r2 = records[y];

    if (r1.ref_name != r2.ref_name)
        return r1.ref_name < r2.ref_name;
    else if (r1.ref_from != r2.ref_from)
        return r1.ref_from < r2.ref_from;
    else
        return r1.ref_to < r2.ref_to;
}

int main(int argc, char *argv[])
{
    OptionsDescription desc;

    try
    {
        desc.Parse(argc, argv);

        if (argc < 3)
            throw invalid_argument("not enough parameters");
    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "validate_contigs_blat - validate contigs by blat." << endl;
        cerr << "Usage: validate_contigs_blat ref.fa contigs.fa." << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    FILE *fin = fopen(argv[1], "rb");
    FILE *fout = fopen(argv[2], "wb");

    while (fgets(line, MaxLine, fin) != NULL)
    {
        BlatRecord record;
        record.Parse(line);

        records.push_back(record);
        lines.push_back(line);
    }

    aux.resize(records.size());
    for (unsigned i = 0; i < records.size(); ++i)
        aux[i] = i;
    sort(aux.begin(), aux.end(), Compare);

    for (unsigned i = 0; i < records.size(); ++i)
        fprintf(fout, "%s", lines[aux[i]].c_str());

    return 0;
}

