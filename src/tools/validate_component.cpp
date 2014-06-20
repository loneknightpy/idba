/**
 * @file validate_contigs_blat.cpp
 * @brief 
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
#include "sequence/sequence_io.h"
#include "misc/utils.h"
#include "sequence/sequence.h"
#include "basic/histgram.h"

using namespace std;

const int MaxLine = (1 << 20);

char line[MaxLine];
char buf[MaxLine];

int main(int argc, char *argv[])
{
    int min_contig = 0;
    double similar = 0.95;
    double complete_rate = 0.8;
    bool is_local = false;
    
    OptionsDescription desc;
    desc.AddOption("min_contig", "", min_contig, "minimum contigs");
    desc.AddOption("similar", "", similar, "similarity");
    desc.AddOption("complete_rate", "", complete_rate, "completeness");
    desc.AddOption("is_local", "", is_local, "local align");

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

    deque<Sequence> refs;
    deque<string> ref_names;
    ReadSequence(argv[1], refs, ref_names);

    deque<Sequence> contigs;
    deque<string> contig_names;
    ReadSequence(argv[2], contigs, contig_names);

    FILE *ffound_list = OpenFile(argv[3], "rb");
    double cover;
    map<string, int> real_found_transcripts;
    while (fscanf(ffound_list, "%s %lf", buf, &cover) != EOF)
    {
        real_found_transcripts[buf] = 1;
    }

    vector<int> is_found(refs.size());
    vector<vector<int> > flags(refs.size());
    map<string, int> dict;
    map<string, int> gene_size;
    map<string, int> gene_num;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        flags[i].resize(refs[i].size(), false);
        size_t index = ref_names[i].find(' ');
        if (index != string::npos)
            ref_names[i].resize(index);
        dict[ref_names[i]] = i;

        string tmp = ref_names[i];
        Replace(tmp, '.', ' ');
        sscanf(tmp.c_str(), "%s", buf);
        gene_size[buf] += refs[i].size();
        gene_num[buf] += 1;
    }

    //map<string, int> contig_dict;
    map<string, int> component_map;
    map<string, int> contig_table;
    int maximum = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        string tmp = contig_names[i];
        Replace(tmp, '_', ' ');
        int id;
        sscanf(tmp.c_str(), "%s %d", buf, &id);
        component_map[contig_names[i]] = id;
        maximum = id;
        contig_table[contig_names[i]] = i;
    }

    deque<int> component_size(maximum + 1);
    deque<int> component_num(maximum + 1);
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        component_size[component_map[contig_names[i]]] += contigs[i].size();
        component_num[component_map[contig_names[i]]] += 1;
    }

    Histgram<int> histgram;
    for (unsigned i = 0; i < component_num.size(); ++i)
    {
        histgram.insert(component_num[i]);
    }

    cout << component_num.size() << endl;
    for (unsigned i = 0; i < 20; ++i)
        cout << i << " " << histgram.count(i, i+1) << endl;
    cout << histgram.count(10, 1000000000) << endl;
    cout << histgram.maximum() << endl;
    
    string blat_file = string(argv[2]) + ".blat";
    FILE *fblat = OpenFile(blat_file, "rb");
    
    map<string, int> found_gene;

    deque<deque<BlatRecord> > match_table(component_size.size());
    //deque<map<string, set<int> > > match_table(component_size.size());
    while (fgets(line, MaxLine, fblat) != NULL)
    {
        BlatRecord record;
        record.Parse(line);

        int ref_id = dict[record.ref_name];

        //if (record.match_count > similar * record.query_length 
        if (record.match_count > similar * (record.query_to - record.query_from)
                && abs(record.query_to - record.query_from) >= min_contig
                && record.match_count > similar * abs(record.ref_to - record.ref_from)
           )
        {
            match_table[component_map[record.query_name]].push_back(record);
        }
    }

    int total_pure = 0;
    int total_decompse = 0;
    int total = 0;
    map<string, int> found_transcripts_in_partial;
    map<string, int> found_partial_gene;
    map<int, int> distribution;
    map<int, int> real_distribution;
    map<string, int> all_got_gran;
    for (unsigned i = 0; i < component_size.size(); ++i)
    {
        deque<BlatRecord> &table = match_table[i];

        map<string, int> transcript_coverage;
        map<string, set<int> > transcript_match;
        for (unsigned j = 0; j < table.size(); ++j)
        {
            BlatRecord record = table[j];

            if (record.match_count > similar * record.query_length)
                transcript_match[record.ref_name].insert(contig_table[record.query_name]);

            transcript_coverage[record.ref_name] = 0;

            int ref_id = dict[record.ref_name];
            for (unsigned x = 0; x < record.blocks.size(); ++x)
            {
                BlatBlock block = record.blocks[x];
                for (unsigned y = block.ref_from; y < block.ref_from + block.size; ++y)
                {
                    
                    flags[ref_id][y] = 1;
                }
            }
        }

        map<string, int> gene_coverage;
        map<string, set<int> > gene_match;

        deque<string> got_tran;
        set<string> got_gene;
        bool pure = false;
        for (map<string, int>::iterator p = transcript_coverage.begin(); p != transcript_coverage.end(); ++p)
        {
            //int sum = 0;
                //sum += contigs[*q].size();



//            if (sum > 0.95 * component_size[i])
//                pure = true;

            string tmp = p->first;
            Replace(tmp, '.', ' ');
            sscanf(tmp.c_str(), "%s", buf);

            int ref_id = dict[p->first];

            for (set<int>::iterator q = transcript_match[p->first].begin(); q != transcript_match[p->first].end(); ++q)
                gene_match[buf].insert(*q);

            int cover = 0;
            for (unsigned x = 0; x < refs[ref_id].size(); ++x)
            {
                gene_coverage[buf] += flags[ref_id][x];
                cover += flags[ref_id][x];
                flags[ref_id][x] = 0;
            }

            if (cover > 0.75 * refs[ref_id].size()
                    //&& all_got_gran.find(p->first) == all_got_gran.end()
               )
            {
                got_tran.push_back(p->first);
                got_gene.insert(buf);
            }

        }

        for (map<string, set<int> >::iterator p = gene_match.begin(); p != gene_match.end(); ++p)
        {
            int sum = 0;
            for (set<int>::iterator q = p->second.begin(); q != p->second.end(); ++q)
                sum += contigs[*q].size();

            if (sum > 0.95 * component_size[i])
                pure = true;
        }

        bool composed = false;
        for (map<string, int>::iterator p = gene_coverage.begin(); p != gene_coverage.end(); ++p)
        {
            if (p->second > 0.8 * gene_size[p->first])
            {
                if (composed)
                {
                    composed = false;
                    break;
                }
                else
                {
                    composed = true;
                }

                //found_gene[p->first] = 1;
            }
        }

        
        if (component_size[i] > 300)
        {
            distribution[got_tran.size()]++;
            for (unsigned i = 0; i < got_tran.size(); ++i)
            {
                all_got_gran[got_tran[i]] = 1;
                if (real_found_transcripts.find(got_tran[i]) != real_found_transcripts.end())
                    real_distribution[got_tran.size()]++;
            }

            if (pure)
            {
                ++total_pure;

                if (composed)
                {
                    ++total_decompse;
                    for (map<string, int>::iterator p = gene_coverage.begin(); p != gene_coverage.end(); ++p)
                    {
                        if (p->second > 0.8 * gene_size[p->first])
                        {
                            found_gene[p->first] = 1;
                            //++distribution[gene_num[p->first]];
                            //distribution[got_tran.size()]++;
                        }
                    }
                }
                else
                {
                    if (got_tran.size() > 0 && got_gene.size() == 1)
                    {
                        for (unsigned x = 0; x < got_tran.size(); ++x)
                        {
                            found_transcripts_in_partial[got_tran[x]] = 1;
                            string tmp = got_tran[x];
                            Replace(tmp, '.', ' ');
                            sscanf(tmp.c_str(), "%s", buf);
                            found_partial_gene[buf] = 1;

                            //distribution[got_tran.size()]++;
                        }
                    }
                    else if (got_tran.size() > 0)
                    {
                        //distribution[got_tran.size()]++;
                    }
                }
            }
            //else
                //distribution[got_tran.size()]++;
            ++total;
        }
    }

    cout << total_pure << " " << total_decompse << " " << total << endl;


    map<string, int> found_transcripts;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        string tmp = ref_names[i];
        Replace(tmp, '.', ' ');
        sscanf(tmp.c_str(), "%s", buf);

        if (found_gene.find(buf) != found_gene.end())
            found_transcripts[ref_names[i]] = 1;
    }
    cout << "found gene " << found_gene.size() <<  " " << found_transcripts.size() << endl;
    cout << "found partial gene " << found_partial_gene.size() << " " << found_transcripts_in_partial.size() << endl;

    int count = 0;
    for (map<string, int>::iterator p = found_gene.begin(); p != found_gene.end(); ++p)
    {
        if (found_partial_gene.find(p->first) != found_partial_gene.end())
            ++count;
    }
    cout << count << endl;

    for (map<int, int>::iterator p = distribution.begin(); p != distribution.end(); ++p)
        cout << p->first << " " << p->second << " " << p->first * p->second << " " << real_distribution[p->first] << endl;

//    for (map<int, int>::iterator p = real_distribution.begin(); p != real_distribution.end(); ++p)
//        cout << p->first << " " << p->second << endl;

//    FILE *fin = OpenFile(argv[1] + string(".list"), "rb");
//    FastaWriter writer(argv[2] + string(".ref"));
//    FILE *fout_list = OpenFile(argv[2] + string(".ref.list"), "wb");
//    
//    double cover;
//    for (unsigned i = 0; i < refs.size(); ++i)
//    {
//        fscanf(fin, "%s %lf", buf, &cover);
//        if (all_got_gran.find(buf) != all_got_gran.end())
//        {
//            writer.Write(refs[i], ref_names[i]);
//            fprintf(fout_list, "%s %.4f\n", buf, cover);
//        }
//    }
//

    
    return 0;
}
