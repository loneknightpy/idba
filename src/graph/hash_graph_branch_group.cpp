/**
 * @file hash_graph_branch_group.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.4
 * @date 2011-09-21
 */

#include "hash_graph_branch_group.h"

#include <algorithm>
#include <cstdio>
#include <iostream>

using namespace std;

bool HashGraphBranchGroup::Search()
{
    branches_.reserve(max_branches_);

    HashGraphPath path;
    path.Append(begin_);
    branches_.push_back(path);

    if (begin_.in_edges().size() != 1 || begin_.out_edges().size() <= 1 || begin_.out_edges().size() > max_branches_)
        return false;

    bool is_converge = false;
    for (int k = 1; k < max_length_; ++k)
    {
        int num_branches = branches_.size();
        for (int i = 0; i < num_branches; ++i)
        {
            HashGraphVertexAdaptor current = branches_[i].back();

            if (current.out_edges().size() == 0)
                return false;

            bool is_first = true;
            HashGraphPath path = branches_[i];
            for (int x = 0; x < 4; ++x)
            {
                if (current.out_edges()[x])
                {
                    Kmer kmer = current.kmer();
                    kmer.ShiftAppend(x);
                    HashGraphVertexAdaptor next = hash_graph_->FindVertexAdaptor(kmer);

                    if (next.status().IsDead())
                        return false;

                    if (is_first)
                    {
                        branches_[i].Append(next);
                        is_first = false;
                    }
                    else
                    {
                        if ((int)branches_.size() == max_branches_)
                            return false;

                        path.Append(next);
                        branches_.push_back(path);
                        path.Pop();
                    }
                }
            }
        }

        end_ = branches_[0].back();

        if (end_.out_edges().size() == 1)
        {
            is_converge = true;
            for (unsigned i = 1; i < branches_.size(); ++i)
            {
                if (branches_[i].back() != end_)
                {
                    is_converge = false;
                    break;
                }
            }

            if (is_converge)
                break;
        }
    }

    return is_converge && begin_ != end_;
}


void HashGraphBranchGroup::Merge()
{
    unsigned best = 0;
    for (unsigned i = 1; i < branches_.size(); ++i)
    {
        if (branches_[i].kmer_count() > branches_[best].kmer_count())
            best = i;
    }

    int kmer_size = begin_.kmer_size();
    for (unsigned i = 0; i < branches_.size(); ++i)
    {
        HashGraphPath &path = branches_[i];
        path.front().out_edges() = 0;
        path.back().in_edges() = 0;
        for (unsigned j = 1; j+1 < path.num_nodes(); ++j)
        {
            path[j].in_edges() = 0;
            path[j].out_edges() = 0;
            path[j].status().SetDeadFlag();
        }
    }

    HashGraphPath &path = branches_[best];
    for (unsigned j = 1; j+1 < path.num_nodes(); ++j)
        path[j].status().ResetDeadFlag();

    for (unsigned j = 0; j+1 < path.num_nodes(); ++j)
    {
        hash_graph_->AddEdge(path[j], path[j+1].kmer()[kmer_size-1]);
    }
}

