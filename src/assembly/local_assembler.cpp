/**
 * @file local_assembler.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2012-01-05
 */

#include "assembly/local_assembler.h"

#include <pthread.h>

#include <cstdlib>
#include <deque>
#include <iostream>
#include <vector>

#include "assembly/assembly_utility.h"
#include "basic/histgram.h"
#include "graph/contig_graph.h"
#include "graph/hash_graph.h"
#include "misc/hash_aligner.h"
#include "sequence/sequence.h"
#include "sequence/short_sequence.h"

using namespace std;

static int Match(const Sequence &a, const Sequence &b)
{
    int match = 0;
    for (unsigned i = 0, j = 0; i < a.size() && j < b.size(); ++i, ++j)
    {
        if (a[i] == b[j])
            ++match;
    }
    return match;
}

LocalAssembler::~LocalAssembler()
{
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)locks_.size(); ++i)
        omp_destroy_lock(&locks_[i]);
}

void LocalAssembler::Initialize(AssemblyInfo &assembly_info, deque<Sequence> &contigs)
{
    assembly_info_ = &assembly_info;
    contigs_ = &contigs;
    local_contigs_.resize(contigs.size() * 2);
    extended_contigs_.resize(contigs.size());
    extensions_.resize(contigs.size() * 2);

    locks_.resize(contigs.size());
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)locks_.size(); ++i)
        omp_init_lock(&locks_[i]);
    in_reads_.resize(contigs.size());
    out_reads_.resize(contigs.size());

    read_length_ = assembly_info.read_length();
}

int64_t LocalAssembler::Assemble(deque<Sequence> &contigs)
{
    if (num_threads_ == 0)
        num_threads_ = omp_get_max_threads();

    contigs.clear();
    omp_set_num_threads(1);

    vector<pthread_t> threads(num_threads_);
    vector<LocalAssemblyTask> tasks(num_threads_);
    for (int i = 0; i < num_threads_; ++i)
    {
        tasks[i].id = i;
        tasks[i].local_assembler = this;
        pthread_create(&threads[i], NULL, LocalAssembleThread, (void *)&tasks[i]);
    }
    
    //cout << "thread" << endl;
    for (int i = 0; i < num_threads_; ++i)
    {
        pthread_join(threads[i], NULL);
        contigs.insert(contigs.end(), tasks[i].local_contigs.begin(), tasks[i].local_contigs.end());
    }
    //cout << "end" << endl;
    contigs = local_contigs();

    omp_set_num_threads(num_threads_);

    return contigs.size();
}

void * LocalAssembler::LocalAssembleThread(void *p)
{
    omp_set_num_threads(1);

    LocalAssemblyTask &task = *(LocalAssemblyTask *)p;
    LocalAssembler &local_assembler = *task.local_assembler;
    deque<Sequence> &contigs = local_assembler.contigs();
    deque<Sequence> &local_contigs = local_assembler.local_contigs();
    deque<Sequence> &extended_contigs = local_assembler.extended_contigs();
    deque<Sequence> &extensions = local_assembler.extensions();
    vector<deque<uint64_t> > &in_reads = local_assembler.in_reads();
    vector<deque<uint64_t> > &out_reads = local_assembler.out_reads();
    int num_threads = local_assembler.num_threads();
    int min_contig = local_assembler.min_contig();
    int local_range = local_assembler.local_range();

    for (int64_t i = task.id; i < (int64_t)contigs.size(); i += num_threads)
    {
        extended_contigs[i] = contigs[i];

        if ((int)contigs[i].size() < min_contig)
            continue;

        for (int strand = 0; strand < 2; ++strand)
        {
            deque<Sequence> one_round_local_contigs;
            if (strand == 0)
                local_assembler.LocalAssemble(contigs[i], in_reads[i], one_round_local_contigs);
            else
                local_assembler.LocalAssemble(contigs[i], out_reads[i], one_round_local_contigs);

            task.local_contigs.insert(task.local_contigs.end(), 
                    one_round_local_contigs.begin(), one_round_local_contigs.end());

            Sequence prefix = contigs[i];
            if ((int64_t)prefix.size() > local_range)
                prefix.resize(local_range);

            for (unsigned j = 0; j < one_round_local_contigs.size(); ++j)
            {
                if (one_round_local_contigs[j].size() >= prefix.size())
                {
                    Sequence suffix(one_round_local_contigs[j], one_round_local_contigs[j].size() - prefix.size());
                    int match = Match(prefix, suffix);

                    if (match < prefix.size() * 0.95)
                    {
                        one_round_local_contigs[j].ReverseComplement();
                        suffix.Assign(one_round_local_contigs[j], one_round_local_contigs[j].size() - prefix.size());
                        match = Match(prefix, suffix);
                    }

                    if (match >= prefix.size() * 0.95)
                    {
                        Sequence seq(one_round_local_contigs[j], 0, one_round_local_contigs[j].size() - prefix.size());
                        seq.Append(extended_contigs[i]);
                        extended_contigs[i] = seq;
                        extensions[2*i + strand] = one_round_local_contigs[j];
                        break;
                    }
                }
            }

            if (extensions[2*i + strand].empty())
                extensions[2*i + strand] = prefix;

            Sequence local_contig;
            for (unsigned j = 0; j < one_round_local_contigs.size(); ++j)
            {
                local_contig += one_round_local_contigs[j];
                local_contig += uint8_t(4);
            }
            local_contigs[2*i + strand] = local_contig;

            contigs[i].ReverseComplement();
            extended_contigs[i].ReverseComplement();
        }
    }

    return p;
}

void LocalAssembler::LocalAssemble(Sequence &contig, std::deque<uint64_t> &local_reads, std::deque<Sequence> &local_contigs)
{
    deque<ShortSequence> &reads = assembly_info_->reads;

    int min_num_reads = median_ / read_length_;

    if ((int)local_reads.size() > min_num_reads)
    {
        sort(local_reads.begin(), local_reads.end());

        AssemblyInfo assembly_info;
        int last = -1;
        int count = 0;
        for (unsigned j = 0; j < local_reads.size(); ++j)
        {
            int offset = local_reads[j] >> 36;
            int index = local_reads[j] & ((1ULL << 32) - 1);

            count = (offset == last) ? count+1 : 1;
            last = offset;

            if (count <= 3)
                assembly_info.reads.push_back(reads[index]);
        }
        Sequence seq = contig;
        if ((int64_t)seq.size() > local_range_)
            seq.resize(local_range_);
        assembly_info.long_reads.push_back(seq);
        assembly_info.ClearStatus();

        IterativeAssemble(assembly_info, local_contigs);
    }
}

void LocalAssembler::IterativeAssemble(AssemblyInfo &assembly_info, deque<Sequence> &contigs)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    Sequence contig_tail = assembly_info.long_reads[0];
    assembly_info.long_reads.clear();
    assembly_info.ClearStatus();

    HashGraph hash_graph(mink_);

    hash_graph.reserve(local_range_*4);
    deque<ContigInfo> contig_infos;
    ContigGraph contig_graph;
    for (int kmer_size = mink_; kmer_size <= min(maxk_, read_length_); kmer_size += step_)
    {
        int64_t sum = 0;
        hash_graph.set_kmer_size(kmer_size);
        hash_graph.clear();
        for (int64_t i = 0; i < (int64_t)reads.size(); ++i)
        {
            if ((int)reads[i].size() < kmer_size)
                continue;

            Sequence seq(reads[i]);
            hash_graph.InsertKmers(seq);
            sum += seq.size() - kmer_size + 1;
        }

        Histgram<int> histgram = hash_graph.coverage_histgram();
        double mean = histgram.percentile(1 - 1.0 * local_range_ / hash_graph.num_vertices());
        double threshold = mean;

        hash_graph.InsertKmers(contig_tail);
        for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
            hash_graph.InsertUncountKmers(contigs[i]);

        hash_graph.Assemble(contigs, contig_infos);
        contig_graph.set_kmer_size(kmer_size);
        contig_graph.Initialize(contigs, contig_infos);
        contig_graph.RemoveDeadEnd(kmer_size*2);
        
        //contig_graph.Trim(10);
        contig_graph.RemoveBubble();
        //contig_graph.MergeSimilarPath();
        contig_graph.IterateCoverage(kmer_size*2, 1, threshold);
        //contig_graph.Prune(kmer_size);

        contig_graph.Assemble(contigs, contig_infos);

        if (contigs.size() == 1)
            break;
    }

    int index = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() > min_contig_)
            contigs[index++].swap(contigs[i]);
    }
    contigs.resize(index);
}

LocalAssemblerWithReference::~LocalAssemblerWithReference()
{
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)locks_.size(); ++i)
        omp_destroy_lock(&locks_[i]);
}

void LocalAssemblerWithReference::Initialize(AssemblyInfo &assembly_info, deque<Sequence> &contigs)
{
    assembly_info_ = &assembly_info;
    contigs_ = &contigs;
    local_contigs_.resize(contigs.size());

    locks_.resize(contigs.size());
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)locks_.size(); ++i)
        omp_init_lock(&locks_[i]);

    read_offset_.resize(assembly_info.reads.size());
    read_tables_.resize(contigs.size());
    matches_.resize(contigs.size());
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
        matches_[i].resize(contigs[i].size());

    read_length_ = assembly_info.read_length();
}

void LocalAssemblerWithReference::set_insert_distance(double median, double sd)
{ 
    median_ = median; 
    sd_ = sd; 
    local_range_ = median; 
    for (int64_t i = 0; i < (int64_t)contigs().size(); ++i)
        read_tables_[i].resize(contigs()[i].size()/local_range_ + 1);
}

void LocalAssemblerWithReference::GetSegments(deque<Sequence> &contigs)
{
    deque<Sequence> &refs = this->contigs();

    contigs.clear();
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        Sequence seq = refs[i];
        deque<int> counts = matches_[i];

        int last = -1;
        for (int i = 0; i <= (int)seq.size(); ++i)
        {
            if (i == (int)seq.size() || counts[i] >= 1)
            {
                if (i - last > max_gap_ + 1 || last == -1 || i == (int)seq.size())
                {
                    {
                        for (int j = last + 1; j < i; ++j)
                            seq[j] = 4;
                    }
                }

                last = i;
            }
        }

        double sum = 0;
        int len = 0;
        for (int i = 0; i < (int)seq.size(); ++i)
        {
            if (counts[i] > 0)
            {
                sum += counts[i];
                len++;
            }
        }
        
        Sequence subseq;
        double subsum = 0;
        int sublen = 0;
        for (unsigned i = 0; i <= seq.size(); ++i)
        {
            if (i == seq.size() || seq[i] == 4)
            {
                {
                    if ((int)subseq.size() > min_contig())
                        contigs.push_back(subseq);
                }

                subseq.resize(0);
                subsum = 0;
                sublen = 0;
            }
            else
            {
                subseq.Append(seq[i]);
                if (counts[i] > 0)
                {
                    sublen++;
                    subsum += counts[i];
                }
            }
        }
    }
}

int64_t LocalAssemblerWithReference::Assemble(deque<Sequence> &contigs)
{
    long long sum = 0;
    long long count = 0;
    for (unsigned i = 0; i < matches_.size(); ++i)
    {
        for (unsigned j = 0; j < matches_[i].size(); ++j)
        {
            if (matches_[i][j])
            {
                sum += matches_[i][j];
                ++count;
            }
        }
    }

    double coverage = 1.0 * sum / count;
    expeceted_coverage_ = coverage;
    cout << "coverage " << coverage << endl;

    coverage_table_.resize(matches_.size());
    for (unsigned i = 0; i < matches_.size(); ++i)
    {
        coverage_table_[i].resize(matches_[i].size());
        double sum = 0;
        int len = 0;
        for (unsigned j = 0; j < matches_[i].size(); ++j)
        {
            while (j + len < matches_[i].size() && len < read_length_)
            {
                sum += matches_[i][j+len];
                ++len;
            }

            coverage_table_[i][j] = sum / len;
            sum -= matches_[i][j];
            --len;

            double tmp = 0;
            double tmp_len = 0;
            for (int k = 0; k < read_length_ && j + k < matches_[i].size(); ++k)
            {
                tmp += matches_[i][j+k];
                ++tmp_len;
            }

            //cout << matches_[i][j] << " " << coverage_table_[i][j] << " " << tmp/tmp_len<< endl;
        }
    }


    if (num_threads_ == 0)
        num_threads_ = omp_get_max_threads();

    contigs.clear();
    omp_set_num_threads(1);

    vector<pthread_t> threads(num_threads_);
    vector<LocalAssemblyTask> tasks(num_threads_);
    for (int i = 0; i < num_threads_; ++i)
    {
        tasks[i].id = i;
        tasks[i].local_assembler = this;
        pthread_create(&threads[i], NULL, LocalAssembleThread, (void *)&tasks[i]);
    }
    
    for (int i = 0; i < num_threads_; ++i)
    {
        pthread_join(threads[i], NULL);
        contigs.insert(contigs.end(), tasks[i].local_contigs.begin(), tasks[i].local_contigs.end());
    }

    omp_set_num_threads(num_threads_);

    return contigs.size();
}

void * LocalAssemblerWithReference::LocalAssembleThread(void *p)
{
    omp_set_num_threads(1);

    LocalAssemblyTask &task = *(LocalAssemblyTask *)p;
    LocalAssemblerWithReference &local_assembler = *task.local_assembler;
    deque<Sequence> &contigs = local_assembler.contigs();
    //deque<Sequence> &local_contigs = local_assembler.local_contigs();
    deque<deque<deque<int> > > read_tables = local_assembler.read_tables();
    deque<deque<int> > &matches = local_assembler.matches();
    deque<ShortSequence> &reads = local_assembler.assembly_info().reads;
    deque<deque<long long> > &read_offset = local_assembler.read_offset();
    deque<deque<double> > &coverage_table = local_assembler.coverage_table();
    int max_gap = local_assembler.max_gap();
    double expeceted_coverage = local_assembler.expeceted_coverage();
    int num_threads = local_assembler.num_threads();
    int local_range = local_assembler.local_range();
//    int read_length = local_assembler.read_length();
//    int min_contig = local_assembler.min_contig();

    for (unsigned x = 0; x < read_tables.size(); ++x)
    {
        for (unsigned y = task.id; y < read_tables[x].size(); y += num_threads)
        {
            //cout << x << " " << y << " " << task.id << endl;
//            if (read_tables[x][y].size() < local_range * 2.0 / read_length)
//                continue;
//
//            if (read_tables[x][y].size() > local_range * 10.0 / read_length)
//                continue;

            AssemblyInfo local_assembly_info;
            for (unsigned i = 0; i < read_tables[x][y].size(); ++i)
            {
                int offset = local_range * y;
                int read_id = read_tables[x][y][i];
                bool flag = false;
                int aligned[2] = {0};
                for (int j = 0; j < 2; ++j)
                {
                    for (int k = 0; k < (int)read_offset[read_id+j].size(); ++k)
                    {
                        int a = read_offset[read_id+j][k] >> 32;
                        int b = read_offset[read_id+j][k] & ((1LL << 32)-1);

                        if (a == (int)x && abs(b - offset) < 3*local_range)
                        {
                            if (coverage_table[a][b] < 1.5*expeceted_coverage)
                                flag = true;
                            aligned[j] = 1;
                        }
                    }
                }

                //if (flag || aligned[0] + aligned[1] == 2)
                {
                    local_assembly_info.reads.push_back(reads[read_tables[x][y][i]]);
                    local_assembly_info.reads.push_back(reads[read_tables[x][y][i] + 1]);
                }
            }
            local_assembly_info.ClearStatus();

            //cout << "hello1" << endl;

            Sequence seq(contigs[x], local_range * y, local_range*2);
            deque<int> counts(matches[x].begin() + local_range * y, matches[x].begin() + local_range * y + seq.size());
            int last = -1;
            for (int i = 0; i <= (int)seq.size(); ++i)
            {
                if (i == (int)seq.size() || counts[i] >= 1)
                {
                    if (i - last > max_gap + 1 || last == -1 || i == (int)seq.size())
                    {
                        {
                            for (int j = last + 1; j < i; ++j)
                                seq[j] = 4;
                        }
                    }

                    last = i;
                }
            }

            double sum = 0;
            int len = 0;
            for (int i = 0; i < (int)seq.size(); ++i)
            {
                //if (seq[i] < 4)
                if (counts[i] > 0)
                {
                    sum += counts[i];
                    len++;
                }
            }
            
            //cout << len << " " << read_tables[x][y].size() * read_length / len << endl;

            if (len < local_range)// || sum/len < 0.5 * expeceted_coverage)
                continue;

            if (sum / len > 10)
                continue;

//            if (read_tables[x][y].size() < len * 1.0 / read_length)
//                continue;
//
//            if (read_tables[x][y].size() > len * 5.0 / read_length)
//                continue;


            Sequence subseq;
            double subsum = 0;
            int sublen = 0;
            for (unsigned i = 0; i <= seq.size(); ++i)
            {
                if (i == seq.size() || seq[i] == 4)
                {
//                    if (subseq.size() >= min_contig && subsum/sublen > 0.5 * expeceted_coverage
//                            && sublen*0.5 > subseq.size() - sublen)
                    {
//                        if (sublen != 0 && subsum / sublen > 2 && subsum / sublen < 10)
                        if (subseq.size() > 0)
                            local_assembly_info.long_reads.push_back(subseq);
                    }

                    subseq.resize(0);
                    subsum = 0;
                    sublen = 0;
                }
                else
                {
                    subseq.Append(seq[i]);
                    if (counts[i] > 0)
                    {
                        sublen++;
                        subsum += counts[i];
                    }
                }
            }

            //local_assembly_info.long_reads.push_back(seq);
            
//            if (sum / len > 10)
//                local_assembly_info.long_reads.clear();

            deque<Sequence> result_contigs;
            local_assembler.IterativeAssemble(local_assembly_info, result_contigs);
            task.local_contigs.insert(task.local_contigs.end(), result_contigs.begin(), result_contigs.end());

            //task.local_contigs.insert(task.local_contigs.end(), local_assembly_info.long_reads.begin(), local_assembly_info.long_reads.end());

            //task.local_contigs.push_back(seq);
        }
    }

    return p;
}

void LocalAssemblerWithReference::IterativeAssemble(AssemblyInfo &assembly_info, deque<Sequence> &contigs)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    assembly_info.ClearStatus();
    contigs.clear();

    HashGraph hash_graph(mink_);

    deque<ContigInfo> contig_infos;
    ContigGraph contig_graph;
    for (int kmer_size = mink_; kmer_size <= min(maxk_, read_length_); kmer_size += step_)
    {
        int64_t sum = 0;
        hash_graph.set_kmer_size(kmer_size);
        hash_graph.clear();
#pragma omp parallel for
        for (int64_t i = 0; i < (int64_t)reads.size(); ++i)
        {
            if ((int)reads[i].size() < kmer_size)
                continue;

            Sequence seq(reads[i]);
            hash_graph.InsertKmers(seq);
            sum += seq.size() - kmer_size + 1;
        }

        Histgram<int> histgram = hash_graph.coverage_histgram();
        double mean = histgram.median();
        //double mean = histgram.percentile(1 - 1.0 * local_range_ / hash_graph.num_vertices());
        //double threshold = sqrt(max(mean, 1.0));
        double threshold = mean/2;

        if (kmer_size < (mink_ + maxk_)/2)
        {
#pragma omp parallel for
            for (int64_t i = 0; i < (int64_t)assembly_info.long_reads.size(); ++i)
                hash_graph.InsertUncountKmers(assembly_info.long_reads[i]);
        }

#pragma omp parallel for
        for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
            hash_graph.InsertUncountKmers(contigs[i]);

        hash_graph.Assemble(contigs, contig_infos);
        contig_graph.set_kmer_size(kmer_size);
        contig_graph.Initialize(contigs, contig_infos);
        contig_graph.RemoveDeadEnd(kmer_size*2);
        contig_graph.RemoveBubble();
        //contig_graph.MergeSimilarPath();
        contig_graph.IterateCoverage(kmer_size*2, 1, threshold);
        contig_graph.Assemble(contigs, contig_infos);

        if (contigs.size() == 1)
            break;
    }

    int index = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() > min_contig_)
                //&& contig_infos[i].kmer_count() >= contigs[i].size() - maxk_ + 1)
            contigs[index++].swap(contigs[i]);
    }
    contigs.resize(index);
}

