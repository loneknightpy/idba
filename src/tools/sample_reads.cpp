/**
 * @file sample_reads.cpp
 * @brief Sample reads from a set of reads.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-10
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
#include "misc/utils.h"
#include "sequence/sequence.h"
#include "sequence/sequence_io.h"


using namespace std;


enum DistributionType
{
    UNIFORM, LOG_NORMAL,
};
const double Pi = 3.1415926535;

double NormalRand();
int SimulateErrors(Sequence &read);
int64_t SimulateReads(Sequence &ref, int num_reads, FastaWriter &writer);
int64_t SimulateReads(Sequence &ref, int num_reads, FastaWriter &writer, FastaWriter &writer_correct);

double depth = 30;
double error_rate = 0.01;
int read_length = 100;
bool is_paired = false;
int insert_distance = 500;
int dist_type = 0;
int delta = -1;
int read_id = 0;
int ref_id = 0;
bool print_correct = false;

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    desc.AddOption("depth", "", depth, "sequencing depth of the sample");
    desc.AddOption("error_rate", "", error_rate, "error rate");
    desc.AddOption("read_length", "", read_length, "read length");
    desc.AddOption("paired", "", is_paired, "if paired-end");
    desc.AddOption("sd", "", delta, "sd");
    desc.AddOption("insert_distance", "", insert_distance, "insert distance");
    desc.AddOption("print_correct", "", print_correct, "output correct reads");

    try
    {
        desc.Parse(argc, argv);

        if (argc < 3)
            throw logic_error("not enough parameters");

    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "sim_reads - Simulate sequencing reads from references." << endl;
        cerr << "Usage: sim_reads ref.fa reads.fa" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    if (delta == -1)
        delta = insert_distance / 10;

    deque<Sequence> refs;
    ReadSequence(argv[1], refs);

    deque<Sequence> split_refs;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        Sequence &ref = refs[i];
        int last = -1;
        for (unsigned j = 0; j < ref.size(); ++j)
        {
            if (ref[j] < 4 && last == -1)
                last = j;
            else if (ref[j] == 4)
            {
                if (last != -1)
                {
                    Sequence seq(ref, last, j - last);
                    if ((int)seq.size() > read_length && (!is_paired || (int)seq.size() > insert_distance))
                        split_refs.push_back(seq);
                    last = -1;
                }
            }
        }

        if (last != -1)
        {
            Sequence seq(ref, last);
            if ((int)seq.size() > read_length && (!is_paired || (int)seq.size() > insert_distance))
                split_refs.push_back(seq);
        }
    }
    split_refs.swap(refs);

    deque<int> num_reads(refs.size());
    deque<double> weight(refs.size());
    double sum_weight = 0;
    int64_t sum_refs = 0;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        weight[i] = refs[i].size();
        sum_weight += weight[i];
        sum_refs += refs[i].size();
    }

    cout << "num of refs = " << refs.size() << endl;
    cout << "depth = " << depth << endl;
    cout << "error rate = " << error_rate << endl;
    cout << "read length = " << read_length << endl;
    cout << "is paired = " << is_paired << endl;
    cout << "insert distance = " << insert_distance << endl;
    cout << "bases of references = " << sum_refs << endl;
    cout << "number of reads = " << int64_t(sum_refs * depth / read_length) << endl;

    if (print_correct)
    {
        FastaWriter writer(argv[2]);
        FastaWriter writer_correct(argv[2] + string(".correct"));
        for (unsigned i = 0; i < refs.size(); ++i)
        {
            num_reads[i] = double((sum_refs * depth / read_length) * (weight[i] / sum_weight) + 0.5);
            num_reads[i] += num_reads[i] & 1;

            cout << "ref " << i << " " << refs[i].size() << " " << num_reads[i] << endl;
            if (!is_paired || (int)refs[i].size() > insert_distance)
                SimulateReads(refs[i], num_reads[i], writer, writer_correct);
        }
    }
    else
    {
        FastaWriter writer(argv[2]);
        for (unsigned i = 0; i < refs.size(); ++i)
        {
            num_reads[i] = double((sum_refs * depth / read_length) * (weight[i] / sum_weight) + 0.5);
            num_reads[i] += num_reads[i] & 1;

            cout << "ref " << i << " " << refs[i].size() << " " << num_reads[i] << endl;
            if (!is_paired || (int)refs[i].size() > insert_distance)
                SimulateReads(refs[i], num_reads[i], writer);
        }
    }

    return 0;
}

double NormalRand()
{
    double x = rand()*1.0 / RAND_MAX;
    double y = rand()*1.0 / RAND_MAX;
    return sqrt(-2*log(x)) * cos(2*Pi*y);
}

int SimulateErrors(Sequence &seq)
{
    int error = 0;
    for (unsigned i = 0; i < seq.size(); ++i)
    {
        if (rand()*1.0 < error_rate*RAND_MAX)
        {
            ++error;
            int c = seq[i];
            while (c == seq[i])
                c = rand()/93 % 4;
            seq[i] = c;
        }
    }
    return error;
}

int64_t SimulateReads(Sequence &ref, int num_reads, FastaWriter &writer)
{
    if ((int)ref.size() < read_length)
        return 0;

    if (is_paired && (int)ref.size() < insert_distance)
        return 0;

    if (!is_paired)
    {
        Sequence seq;
        for (int i = 0; i < num_reads; ++i)
        {
            int offset = 0;
            //while (true)
            {
                offset = rand() % (ref.size() - read_length + 1);
                seq.Assign(ref, offset, read_length);
//                if (seq.IsValid())
//                    break;
            }

            SimulateErrors(seq);

            if (rand() < RAND_MAX/2)
                seq.ReverseComplement();

            writer.Write(seq, FormatString("read_%d_%d", ref_id, read_id++));
        }
    }
    else
    {
        Sequence seq1;
        Sequence seq2;
        for (int i = 0; i < num_reads; i += 2)
        {
            int d = delta * NormalRand();
            if (d + insert_distance > (int)ref.size())
                d = ref.size() - insert_distance;

            int offset = 0;
            Sequence seq;
            //while (true)
            {
                offset = rand() % (ref.size() + 1 - insert_distance - d);
                seq.Assign(ref, offset, insert_distance + d);
//                if (seq.IsValid())
//                    break;
            }

            seq1.Assign(seq, 0, read_length);
            seq2.Assign(seq, seq.size() - read_length, read_length);
            seq2.ReverseComplement();

            SimulateErrors(seq1);
            SimulateErrors(seq2);

            if (rand() < RAND_MAX/2)
                swap(seq1, seq2);

            writer.Write(seq2, FormatString("read%d_%d/1", ref_id, read_id));
            writer.Write(seq1, FormatString("read%d_%d/2", ref_id, read_id));
            read_id += 2;
        }
    }

    ++ref_id;

    return num_reads;
}

int64_t SimulateReads(Sequence &ref, int num_reads, FastaWriter &writer, FastaWriter &writer_correct)
{
    if ((int)ref.size() < read_length)
        return 0;

    if (is_paired && (int)ref.size() < insert_distance)
        return 0;

    if (!is_paired)
    {
        Sequence seq;
        for (int i = 0; i < num_reads; ++i)
        {
            int offset = 0;
            while (true)
            {
                offset = rand() % (ref.size() - read_length + 1);
                seq.Assign(ref, offset, read_length);
                if (seq.IsValid())
                    break;
            }

            Sequence correct_seq(seq);
            SimulateErrors(seq);

            if (rand() < RAND_MAX/2)
                seq.ReverseComplement();

            writer_correct.Write(correct_seq, FormatString("read_%d_%d", ref_id, read_id));
            writer.Write(seq, FormatString("read_%d_%d", ref_id, read_id++));
        }
    }
    else
    {
        Sequence seq1;
        Sequence seq2;
        for (int i = 0; i < num_reads; i += 2)
        {
            int d = delta * NormalRand();
            if (d + insert_distance > (int)ref.size())
                d = ref.size() - insert_distance;

            int offset = 0;
            Sequence seq;
            while (true)
            {
                offset = rand() % (ref.size() + 1 - insert_distance - d);
                seq.Assign(ref, offset, insert_distance + d);
                if (seq.IsValid())
                    break;
            }

            seq1.Assign(seq, 0, read_length);
            seq2.Assign(seq, seq.size() - read_length, read_length);
            seq2.ReverseComplement();

            Sequence correct_seq1(seq1);
            Sequence correct_seq2(seq2);

            SimulateErrors(seq1);
            SimulateErrors(seq2);

            if (rand() < RAND_MAX/2)
            {
                swap(seq1, seq2);
                swap(correct_seq1, correct_seq2);
            }

            writer_correct.Write(correct_seq2, FormatString("read%d_%d/1", ref_id, read_id));
            writer.Write(seq2, FormatString("read%d_%d/1", ref_id, read_id));
            writer_correct.Write(correct_seq1, FormatString("read%d_%d/2", ref_id, read_id));
            writer.Write(seq1, FormatString("read%d_%d/2", ref_id, read_id));
            read_id += 2;
        }
    }

    ++ref_id;

    return num_reads;
}

