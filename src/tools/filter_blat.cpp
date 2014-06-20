/**
 * @file filter_blat.cpp
 * @brief Filter out low quality blat records.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.3
 * @date 2011-09-06
 */

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "misc/blat_record.h"
#include "misc/options_description.h"
#include "misc/utils.h"
#include "sequence/sequence.h"
#include "sequence/sequence_io.h"

using namespace std;

const int MaxLine = (1 << 20);

char line[MaxLine];
char buf[MaxLine];

int main(int argc, char *argv[])
{
    int min_contig = 100;
    double similar = 0.95;
    bool is_local = false;
    
    OptionsDescription desc;
    desc.AddOption("min_contig", "", min_contig, "minimum contigs");
    desc.AddOption("similar", "", similar, "similarity");
    desc.AddOption("is_local", "", is_local, "is local");

    try
    {
        desc.Parse(argc, argv);
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

    FILE *fin = OpenFile(argv[1], "rb");
    FILE *fout = OpenFile(argv[2], "wb");

    while (fgets(line, MaxLine, fin) != NULL)
    {
        BlatRecord record;
        record.Parse(line);

        if ((is_local || record.match_count > similar * record.query_length)
                && record.query_length >= min_contig
                && record.match_count >= min_contig
                && record.match_count > similar * abs(record.ref_to - record.ref_from))
        {
            fprintf(fout, "%s", line);
        }
    }
    
    return 0;
}

