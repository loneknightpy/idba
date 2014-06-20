/**
 * @file hash_graph.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#include "graph/hash_graph.h"

#include <algorithm>
#include <deque>
#include <stdexcept>

#include "basic/bit_operation.h"
#include "basic/histgram.h"
#include "basic/kmer.h"
#include "container/hash_table.h"
#include "graph/contig_builder.h"
#include "graph/contig_info.h"
#include "graph/hash_graph_branch_group.h"
#include "graph/hash_graph_vertex.h"
#include "sequence/sequence.h"

using namespace std;


#include <iostream>

int64_t HashGraph::InsertKmersWithPrefix(const Sequence &seq, uint64_t prefix, uint64_t mask)
{
    if (seq.size() < kmer_size_)
        return 0;

    Kmer kmer(kmer_size_);
    int length = 0;
    int64_t num_kmers = 0;
    for (uint64_t i = 0; i < seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[i] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        Kmer key = kmer.unique_format();
        if ((((key.hash() * 10619863ULL + 17977) % 790738119649411319ULL) & mask) == prefix)
        {
            HashGraphVertex &vertex = vertex_table_.find_or_insert(HashGraphVertex(key));
            vertex.count() += 1;
            HashGraphVertexAdaptor adaptor(&vertex, kmer != key);

            if (length > (int)kmer_size_ && seq[i-kmer_size_] < 4)
                adaptor.in_edges().Add(3 - seq[i-kmer_size_]);
            if (i+1 < seq.size() && seq[i+1] < 4)
                adaptor.out_edges().Add(seq[i+1]);

            ++num_kmers;
        }
    }

    return num_kmers;
}

int64_t HashGraph::InsertUncountKmers(const Sequence &seq)
{
    if (seq.size() < kmer_size_)
        return 0;

    Kmer kmer(kmer_size_);
    int length = 0;
    int64_t num_kmers = 0;
    for (uint64_t i = 0; i < seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[i] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        Kmer key = kmer.unique_format();

        HashGraphVertex &vertex = vertex_table_.find_or_insert(HashGraphVertex(key));
        HashGraphVertexAdaptor adaptor(&vertex, kmer != key);

        if (length > (int)kmer_size_ && seq[i-kmer_size_] < 4)
            adaptor.in_edges().Add(3 - seq[i-kmer_size_]);
        if (i+1 < seq.size() && seq[i+1] < 4)
            adaptor.out_edges().Add(seq[i+1]);

        ++num_kmers;
    }

    return num_kmers;
}

int64_t HashGraph::InsertInternalKmers(const Sequence &seq, int min_count)
{
    if (seq.size() < kmer_size_)
        return 0;

    Kmer kmer(kmer_size_);
    int length = 0;
    int64_t num_kmers = 0;
    deque<int> found_index;
    deque<HashGraphVertexAdaptor> found_kmer;
    for (uint64_t i = 0; i < seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[i] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        HashGraphVertexAdaptor adaptor = FindVertexAdaptor(kmer);
        if (adaptor.is_null())
            continue;

        if (length > (int)kmer_size_ && seq[i-kmer_size_] < 4)
            adaptor.in_edges().Add(3 - seq[i-kmer_size_] + 4);
        if (i+1 < seq.size() && seq[i+1] < 4)
            adaptor.out_edges().Add(seq[i+1] + 4);

        if (adaptor.count() >= min_count)
        {
            found_index.push_back(i);
            found_kmer.push_back(adaptor);
        }
    }

    deque<int> flags(seq.size(), 0);
    for (uint64_t i = 0; i+1 < found_index.size(); ++i)
    {
        HashGraphVertexAdaptor from = found_kmer[i]; //FindVertexAdaptor(found_kmer[i]);
        HashGraphVertexAdaptor to = found_kmer[i+1]; //FindVertexAdaptor(found_kmer[i+1]);

        if (from.is_null() || to.is_null())
        {
            cout << "error" << endl;
            continue;
        }

        if ((from.out_edges() & 15) == 0 && (to.in_edges() & 15) == 0)
        {
            for (int j = found_index[i] + 1; j < found_index[i+1]; ++j)
                flags[j] = 1;
        }
    }

    if (found_index.size() > 0)
    {
        if (found_kmer.front().in_edges() == 0)
        {
            for (int j = kmer_size_ - 1; j < found_index.front(); ++j)
                flags[j] = 1;
        }

        if (found_kmer.back().out_edges() == 0)
        {
            for (int j = found_index.back() + 1; j < (int)seq.size(); ++j)
                flags[j] = 1;
        }
    }

    length = 0;
    for (uint64_t i = 0; i < seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[i] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        if (flags[i])
        {
            Kmer key = kmer.unique_format();

            HashGraphVertex &vertex = vertex_table_.find_or_insert(HashGraphVertex(key));
            vertex.count() += 1;
            HashGraphVertexAdaptor adaptor(&vertex, kmer != key);

            if (length > (int)kmer_size_ && seq[i-kmer_size_] < 4)
                adaptor.in_edges().Add(3 - seq[i-kmer_size_] + 4);
            if (i+1 < seq.size() && seq[i+1] < 4)
                adaptor.out_edges().Add(seq[i+1] + 4);

            ++num_kmers;
        }
    }

    return num_kmers;
}

int64_t HashGraph::InsertEdges(const Sequence &seq)
{
    if (seq.size() < kmer_size_)
        return 0;

    Kmer kmer(kmer_size_);
    int length = 0;
    int64_t num_kmers = 0;
    for (uint64_t i = 0; i < seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[i] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        HashGraphVertexAdaptor adaptor = FindVertexAdaptor(kmer);

        if (adaptor.is_null())
            continue;

        if (length > (int)kmer_size_ && seq[i-kmer_size_] < 4)
            adaptor.in_edges().Add(3 - seq[i-kmer_size_]);
        if (i+1 < seq.size() && seq[i+1] < 4)
            adaptor.out_edges().Add(seq[i+1]);
    }

    return num_kmers;
}

int64_t HashGraph::InsertExistKmers(const Sequence &seq)
{
    if (seq.size() < kmer_size_)
        return 0;

    Kmer kmer(kmer_size_);
    int length = 0;
    int64_t num_kmers = 0;
    for (uint64_t i = 0; i < seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[i] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        HashGraphVertexAdaptor adaptor = FindVertexAdaptor(kmer);

        if (adaptor.is_null())
            continue;

        adaptor.count() += 1;
        if (length > (int)kmer_size_ && seq[i-kmer_size_] < 4)
            adaptor.in_edges().Add(3 - seq[i-kmer_size_]);
        if (i+1 < seq.size() && seq[i+1] < 4)
            adaptor.out_edges().Add(seq[i+1]);
    }

    return num_kmers;
}


int64_t HashGraph::RemoveKmers(const Sequence &seq)
{
    if (seq.size() < kmer_size_)
        return 0;

    Kmer kmer(kmer_size_);
    int length = 0;
    int64_t num_kmers = 0;
    for (uint64_t i = 0; i < seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[i] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        Kmer key = kmer.unique_format();
        HashGraphVertex &vertex = *vertex_table_.find(key);
        vertex.status().SetDeadFlag();

        ++num_kmers;
    }

    return num_kmers;
}

int64_t HashGraph::ErodeEnd(int min_cover)
{
    ErodeFunc func(this, min_cover);
    vertex_table_.for_each(func);
    uint64_t num_eroded_vertice = RefreshVertices();
    RefreshEdges();
    ClearStatus();
    return num_eroded_vertice;
}

int64_t HashGraph::Trim(int min_length)
{
    deque<Sequence> contigs;
    deque<ContigInfo> contig_infos;
    Assemble(contigs, contig_infos);

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        if ((contig_infos[i].out_edges() == 0 || contig_infos[i].in_edges() == 0)
                && (int)contigs[i].size() < min_length + (int)kmer_size_ - 1)
            RemoveKmers(contigs[i]);
    }

    uint64_t old_num_vertices = vertex_table_.size();
    Refresh();

    return old_num_vertices - vertex_table_.size();
}

int64_t HashGraph::RemoveDeadEnd(int min_length)
{
    uint64_t num_deadend = 0;
    int l = 1;
    while (true)
    {
        l = min(2*l, min_length);
        num_deadend += Trim(l);

        if (l == min_length)
            break;
    }
    num_deadend += Trim(min_length);
    return num_deadend;
}

int64_t HashGraph::RemoveLowCoverage(double min_cover, int min_length)
{
    uint64_t old_num_vertices = vertex_table_.size();

    int l = 1;
    while (true)
    {
        l = min(2*l, min_length);

        deque<Sequence> contigs;
        deque<ContigInfo> contig_infos;
        Assemble(contigs, contig_infos);

#pragma omp parallel for
        for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
        {
            if (contig_infos[i].kmer_count() * 1.0 / (contigs[i].size() - kmer_size_ + 1)  < min_cover
                    && (int)contigs[i].size() < l + (int)kmer_size_ - 1)
                RemoveKmers(contigs[i]);
        }
        Refresh();
        Trim(l);

        if (l == min_length)
            break;
    }

    return old_num_vertices - vertex_table_.size();
}

int64_t HashGraph::RemoveBubble()
{
    BubbleFunc func(this);
    vertex_table_.for_each(func);
    deque<HashGraphVertexAdaptor> &candidates = func.candidates();

    int64_t bubble = 0;
    for (unsigned i = 0; i < candidates.size(); ++i)
    {
        HashGraphVertexAdaptor current = candidates[i];

        if (current.out_edges().size() > 1 && current.in_edges().size() == 1)
        {
            HashGraphBranchGroup branch_group(this, current, 4, kmer_size_*2 + 2);

            if (branch_group.Search())
            {
                HashGraphVertexAdaptor begin = branch_group.begin();
                HashGraphVertexAdaptor end = branch_group.end();

                begin.ReverseComplement();
                end.ReverseComplement();
                std::swap(begin, end);
                HashGraphBranchGroup rev_branch_group(this, begin, 4, kmer_size_*2 + 2);

                if (rev_branch_group.Search() && rev_branch_group.end() == end)
                {
                    branch_group.Merge();
                    ++bubble;
                }
            }
        }
    }

//    for (HashGraph::iterator p = begin(); p != end(); ++p)
//    {
//        for (int strand = 0; strand < 2; ++strand)
//        {
//            HashGraphVertexAdaptor current(&*p, strand);
    for (unsigned i = 0; i < candidates.size(); ++i)
    {
        HashGraphVertexAdaptor current = candidates[i];

        if (current.out_edges().size() > 1 && current.in_edges().size() == 1)
        {
            HashGraphBranchGroup branch_group(this, current, 4, kmer_size_ + 2);

            if (branch_group.Search())
            {
                HashGraphVertexAdaptor begin = branch_group.begin();
                HashGraphVertexAdaptor end = branch_group.end();

                begin.ReverseComplement();
                end.ReverseComplement();
                std::swap(begin, end);
                HashGraphBranchGroup rev_branch_group(this, begin, 4, kmer_size_ + 2);

                if (rev_branch_group.Search() && rev_branch_group.end() == end)
                {
                    branch_group.Merge();
                    ++bubble;
                }
            }
        }
        //}
    }

    Refresh();

    return bubble;
}

int64_t HashGraph::Assemble(std::deque<Sequence> &contigs)
{
    deque<ContigInfo> contig_infos;
    return Assemble(contigs, contig_infos); 
}

int64_t HashGraph::Assemble(std::deque<Sequence> &contigs, std::deque<ContigInfo> &contig_infos)
{
    contigs.clear();
    contig_infos.clear();
    AssembleFunc func(this);
    vertex_table_.for_each(func);
    contigs.swap(func.contigs());
    contig_infos.swap(func.contig_infos());
    ClearStatus();
    return contigs.size();
}

void HashGraph::ErodeFunc::operator ()(HashGraphVertex &vertex)
{
    if (vertex.in_edges().size() > 0 && vertex.out_edges().size() > 0)
        return;

    if (vertex.status().IsDead())
        return;

    if (vertex.count() < min_cover_)
    {
        vertex.status().SetDeadFlag();

        for (int strand = 0; strand < 2; ++strand)
        {
            HashGraphVertexAdaptor current(&vertex, strand);
            for (int x = 0; x < 4; ++x)
            {
                if (current.out_edges()[x])
                {
                    current.out_edges().Remove(x);
                    Kmer kmer = current.kmer();
                    kmer.ShiftAppend(x);
                    HashGraphVertexAdaptor next = hash_graph_->FindVertexAdaptor(kmer);
                    if (!next.is_null())
                    {
                        next.in_edges().Remove(3 - current.kmer()[0]);
                        (*this)(next.vertex());
                    }
                }
            }
        }
    }
}

void HashGraph::TrimFunc::operator ()(HashGraphVertex &vertex)
{
    if (vertex.in_edges().size() > 0 && vertex.out_edges().size() > 0)
        return;

    if (vertex.kmer().IsPalindrome())
        return;

    if (!vertex.status().Lock(omp_get_thread_num()))
        return;

    for (int strand = 0; strand < 2; ++strand)
    {
        HashGraphVertexAdaptor current(&vertex, strand);
        
        if (current.in_edges().size() > 0)
            continue;

        deque<HashGraphVertexAdaptor> path;
        path.push_back(current);
        for (int i = 0; i < min_length_; ++i)
        {
            if (current.out_edges().size() != 1)
                return;

            Kmer next_kmer = current.kmer();
            next_kmer.ShiftAppend(bit_operation::BitToIndex(current.out_edges()));
            HashGraphVertexAdaptor next = hash_graph_->FindVertexAdaptor(next_kmer);
            if (next.in_edges().size() != 1)
                break;

            if (!next.status().LockPreempt(omp_get_thread_num()))
                return;

            current = next;
            path.push_back(current);
        }

        if ((int)path.size() < min_length_)
        {
            for (unsigned i = 0; i < path.size(); ++i)
                path[i].status().SetDeadFlag();
        }
    }
}

void HashGraph::BubbleFunc::operator ()(HashGraphVertex &vertex)
{
    for (int strand = 0; strand < 2; ++strand)
    {
        HashGraphVertexAdaptor current(&vertex, strand);

        if (current.out_edges().size() > 1 && current.in_edges().size() == 1)
        {
            HashGraphBranchGroup branch_group(hash_graph_, current, 4, hash_graph_->kmer_size()*2 + 2);

            if (branch_group.Search())
            {
                HashGraphVertexAdaptor begin = branch_group.begin();
                HashGraphVertexAdaptor end = branch_group.end();

                begin.ReverseComplement();
                end.ReverseComplement();
                std::swap(begin, end);
                HashGraphBranchGroup rev_branch_group(hash_graph_, begin, 4, hash_graph_->kmer_size()*2 + 2);

                if (rev_branch_group.Search() && rev_branch_group.end() == end)
                {
                    omp_set_lock(&bubble_lock_);
                    candidates_.push_back(current);
                    omp_unset_lock(&bubble_lock_);
                }
            }
        }
    }
}

void HashGraph::AssembleFunc::operator ()(HashGraphVertex &vertex)
{
    if (!vertex.status().Lock(omp_get_thread_num()))
        return;

    ContigBuilder contig_builder;
    contig_builder.Append(HashGraphVertexAdaptor(&vertex));

    if (!vertex.kmer().IsPalindrome())
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            HashGraphVertexAdaptor current(&vertex, strand);
            
            while (true)
            {
                HashGraphVertexAdaptor next;
                if (!hash_graph_->GetNextVertexAdaptor(current, next))
                    break;

                if (hash_graph_->IsLoop(contig_builder.contig(), next))
                    break;

                if (!next.status().LockPreempt(omp_get_thread_num()))
                    return;

                contig_builder.Append(next);
                current = next;
            }

            contig_builder.ReverseComplement();
        }
    }

    omp_set_lock(&contig_lock_);
    contigs_.push_back(contig_builder.contig());
    contig_infos_.push_back(contig_builder.contig_info());
    omp_unset_lock(&contig_lock_);
}

