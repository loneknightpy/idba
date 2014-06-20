/**
 * @file validate_contigs_mummer.cpp
 * @brief Use Mummer alignment result to validate contigs.
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

int min_contig = 100;
double similar = 0.95;
double min_identity = 0.95;
int max_gap = 10000000;
int max_overlap = 1000;

char line[MaxLine];
char buf[MaxLine];
char ref_name_buf[MaxLine];
char query_name_buf[MaxLine];

struct MummerRecord
{
    int ref_id;
    int query_id;
    bool is_reverse;
    int query_from;
    int query_to;
    int query_length;
    int ref_from;
    int ref_to;
    int ref_length;

    void Init(int query_from, int query_to, int query_length, 
            int ref_from, int ref_to, int ref_length)
    {
        this->query_from = query_from;
        this->query_to = query_to;
        this->query_length = query_length;
        this->ref_from = ref_from;
        this->ref_to = ref_to;
        this->ref_length = ref_length;
    }

    void ReverseComplement()
    {
        is_reverse = !is_reverse;

        swap(query_from, query_to);
        query_from = query_length - query_from;
        query_to = query_length - query_to;

        swap(ref_from, ref_to);
        ref_from = ref_length - ref_from;
        ref_to = ref_length - ref_to;
    }

    bool operator <(const MummerRecord &x) const
    {
        if (ref_id != x.ref_id)
            return ref_id < x.ref_id;
        else if (is_reverse != x.is_reverse)
            return is_reverse < x.is_reverse;
        else if (query_from != x.query_from)
            return query_from < x.query_from;
        else
            return query_to < x.query_to;
    }

    bool IsConsisnte(const MummerRecord &x) const
    {
        return query_to < x.query_from + max_overlap && query_to + max_gap > x.query_from
            && ref_to < x.ref_from + max_overlap && ref_to + max_gap > x.ref_from;
    }
};

int LCS(deque<MummerRecord> &v, deque<MummerRecord> &longest)
{
    vector<int> prev(v.size(), -1);
    vector<int> score(v.size(), 0);

    int sum = 0;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        prev[i] = -1;
        score[i] = v[i].query_to - v[i].query_from;

        sum += score[i];
        //cout << v[i].query_from << " " << v[i].query_to << " " << v[i].ref_from << " " << v[i].ref_to << " " << score[i] << endl;
    }
//    cout << endl;
    cout << "sum " << sum << endl;

    for (unsigned i = 0; i < v.size(); ++i)
    {
        for (unsigned j = i+1; j < v.size(); ++j)
        {
            int s = v[j].query_to - v[j].query_from;
            if (v[i].IsConsisnte(v[j]) && score[i] + s > score[j])
            {
                score[j] = score[i] + s;
                prev[j] = i;
            }
            else
            {
//                cout << "no " << endl;
//                cout << v[i].query_from << " " << v[i].query_to << " " << v[i].ref_from << " " << v[i].ref_to << endl;
//                cout << v[j].query_from << " " << v[j].query_to << " " << v[j].ref_from << " " << v[j].ref_to << endl;
            }
        }
    }

    int best = 0;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (score[i] > score[best])
            best = i;
    }
    cout << best << " " << score[best] << " " << v[0].query_length << endl;

    int k = best;
    longest.resize(0);
    while (k != -1)
    {
        longest.push_back(v[k]);
        k = prev[k];
    }
    reverse(longest.begin(), longest.end());

    //cout << longest.size() << endl;
//
    return score[best];
}

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    desc.AddOption("min_contig", "", min_contig, "minimum contigs");
    desc.AddOption("similar", "", similar, "similarity");

    try
    {
        desc.Parse(argc, argv);
    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "validate_contigs_mummer - validate contigs by mummer." << endl;
        cerr << "Usage: validate_contigs_mummer ref.fa contigs.fa." << endl;
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

    vector<vector<bool> > flags(refs.size());
    map<string, int> ref_table;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        flags[i].resize(refs[i].size(), false);
        size_t index = ref_names[i].find(' ');
        if (index != string::npos)
            ref_names[i].resize(index);
        ref_table[ref_names[i]] = i;
    }

    map<string, int> contig_table;
    int num_gaps = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        size_t index = contig_names[i].find(' ');
        if (index != string::npos)
            contig_names[i].resize(index);
        contig_table[contig_names[i]] = i;

        bool is_new_gap = true;
        for (unsigned j = 0; j < contigs[i].size(); ++j)
        {
            if (contigs[i][j] == 4)
            {
                if (is_new_gap)
                {
                    is_new_gap = false;
                    ++num_gaps;
                }
            }
            else
                is_new_gap = true;
        }
    }

    string coords_file = string(argv[2]) + ".mummer.coords";
    FILE *fcoord = OpenFile(coords_file, "rb");

    fgets(line, MaxLine, fcoord);
    fgets(line, MaxLine, fcoord);
    fgets(line, MaxLine, fcoord);
    fgets(line, MaxLine, fcoord);
    fgets(line, MaxLine, fcoord);

    deque<deque<MummerRecord> > records(contigs.size());
    while (fgets(line, MaxLine, fcoord) != NULL)
    {
        int ref_from;
        int ref_to;
        int query_from;
        int query_to;
        double identity;

        sscanf(line, "%d %d %s %d %d %s %s %s %s %lf %s %s %s", 
                &ref_from, &ref_to, buf, &query_from, &query_to, buf, 
                buf, buf, buf, &identity, buf, ref_name_buf, query_name_buf);
        identity /= 100;

        //if (identity < similar)
        if (identity < min_identity)
            continue;
        else
        {
            int query_id = contig_table[query_name_buf];
            int query_length = contigs[query_id].size();
            int ref_id = ref_table[ref_name_buf];
            int ref_length = refs[ref_id].size();
            bool is_reverse = false;

            if (query_from < query_to)
            {
                --ref_from;
                --query_from;
            }
            else
            {
                --ref_from;
                --query_to;

                swap(query_from, query_to);
                is_reverse = true;
                swap(ref_from, ref_to);
                ref_from = ref_length - ref_from;
                ref_to = ref_length - ref_to;

                //cout <<  line << endl;
            }

            MummerRecord record;
            record.query_id = query_id;
            record.ref_id = ref_id;
            record.is_reverse = is_reverse;
            record.query_from = query_from;
            record.query_to = query_to;
            record.query_length = query_length;
            record.ref_from = ref_from;
            record.ref_to = ref_to;
            record.ref_length = ref_length;

            records[query_id].push_back(record);

            if (record.ref_to - record.ref_from >= min_contig)
            {
                if (record.is_reverse)
                    record.ReverseComplement();

                for (int k = record.ref_from; k < record.ref_to; ++k)
                {
                    flags[record.ref_id][k] = 1;
                }
            }
        }
    }

    cout << "hello" << endl;

    map<string, int> valid_contigs;
    deque<int> valid_lengths;
    for (unsigned i = 0; i < records.size(); ++i)
    {
        sort(records[i].begin(), records[i].end());
        int index = 0;
        while (index < (int)records[i].size())
        {
            //cout << index << " " << records[i].size() << endl;
            deque<MummerRecord> v;
            v.push_back(records[i][index]);
            int k = index + 1;
            while (k < (int)records[i].size()
                    && v.back().ref_id == records[i][k].ref_id && v.back().is_reverse == records[i][k].is_reverse)
            {
                v.push_back(records[i][k++]);
            }
            index = k;
            
            deque<MummerRecord> longest;
            int match_length = LCS(v, longest);

            if (match_length > similar * contigs[i].size())
            {
                for (unsigned j = 0; j < longest.size(); ++j)
                {
                    MummerRecord record = longest[j];
                    if (record.is_reverse)
                        record.ReverseComplement();

                    for (int k = record.ref_from; k < record.ref_to; ++k)
                    {
                        flags[record.ref_id][k] = 1;
                    }
                }

                valid_contigs[contig_names[i]] = 1;
                valid_lengths.push_back(match_length);
            }
        }

    }

    long long count = 0;
    long long total = 0;
    for (unsigned k = 0; k < flags.size(); ++k)
    {
        for (unsigned i = 0; i < flags[k].size(); ++i)
        {
            if (flags[k][i])
                ++count;
            ++total;
        }
    }

    sort(valid_lengths.begin(), valid_lengths.end());
    reverse(valid_lengths.begin(), valid_lengths.end());

    long long n50 = 0;
    long long sum = 0;
    long long n80 = 0;

    for (unsigned i = 0; i < valid_lengths.size(); ++i)
    {
        sum += valid_lengths[i];
        if (sum >= 0.5 * total && n50 == 0)
            n50 = valid_lengths[i];
        if (sum >= 0.8 * total && n80 == 0)
            n80 = valid_lengths[i];
    }
    cout << "total " << total << " " << sum << endl;

    long long maximum = 0;
    long long mean = 0;
    if (valid_lengths.size() > 0)
    {
        maximum = valid_lengths[0];
        mean = sum / valid_lengths.size();
    }

    long long sum_wrong = 0;
    long long num_wrong = 0;
    long long corret_contigs = 0;
    long long sum_corret = 0;
    int last_id = 0;
    int last_error = 0;
    FastaWriter error_writer(FormatString("%s.error.fa", argv[2]));
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() < min_contig)
            continue;

        if (valid_contigs.find(contig_names[i]) == valid_contigs.end())
        {
            ++num_wrong;
            sum_wrong += contigs[i].size();
            error_writer.Write(contigs[i], contig_names[i]);
        }
        else
        {
            last_id = i;
            last_error = sum_wrong;

            ++corret_contigs;
            sum_corret += contigs[i].size();
        }
    }

    int num_mismatch = 0;
    printf("last id %d %d total contigs %d gaps %d\n", last_id, last_error, (int)(num_wrong + corret_contigs), num_gaps);
    printf("contigs: %lld N50: %lld coverage: %.2f%% max: %lld mean: %lld total: %lld/%lld N80: %lld\n",
            (long long)valid_contigs.size(), n50, count * 100.0 / total, maximum, mean, count, total, n80);
    printf("substitution error: %.4f%% wrong contigs: %lld %lld correct: %lld %lld %s\n", 
            num_mismatch * 100.0 /sum, num_wrong, sum_wrong, corret_contigs, sum_corret, argv[2]);

    deque<int> lengths;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        int last = 0;
        for (unsigned j = 0; j < refs[i].size(); ++j)
        {
            if (flags[i][j] == 0)
            {
                if (flags[i][last])
                {
                    lengths.push_back(j - last);
                    last = j;
                }
            }
            else
            {
                if (flags[i][last] == 0)
                    last = j;
            }
        }
    }
    sort(lengths.begin(), lengths.end());
    reverse(lengths.begin(), lengths.end());
    
    return 0;
}

