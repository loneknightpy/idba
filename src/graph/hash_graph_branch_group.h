/**
 * @file hash_graph_branch_group.h
 * @brief HashGraphBranchGroup Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.4
 * @date 2011-09-21
 */

#ifndef __GRAPH_HASH_GRAPH_BRANCH_GROUP_H_

#define __GRAPH_HASH_GRAPH_BRANCH_GROUP_H_

#include "basic/kmer.h"
#include "graph/hash_graph.h"
#include "graph/hash_graph_path.h"

#include <vector>


/**
 * @brief It is used to contain a branch group in de Bruijn graph (HashGraph).
 */
class HashGraphBranchGroup
{
public:
    HashGraphBranchGroup(HashGraph *graph, HashGraphVertexAdaptor begin, 
            int max_branches = 2, int max_length = 0)
    {
        hash_graph_ = graph;
        begin_ = begin;
        max_branches_ = max_branches;
        max_length_ = max_length;

        if (max_length_ == 0)
            max_length_ = begin_.kmer().size() + 2;
    }

    bool Search();
    void Merge();

    HashGraphVertexAdaptor begin() { return begin_; }
    HashGraphVertexAdaptor end() { return end_; }

private:
    HashGraph *hash_graph_;
    HashGraphVertexAdaptor begin_;
    HashGraphVertexAdaptor end_;
    std::vector<HashGraphPath> branches_;
    int max_branches_;
    int max_length_;
};

#endif

