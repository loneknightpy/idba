/**
 * @file scaffold_graph_path.h
 * @brief ScaffoldGraphPath Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.10
 * @date 2011-08-16
 */

#ifndef __GRAPH_SCAFFOLD_GRAPH_PATH_H_

#define __GRAPH_SCAFFOLD_GRAPH_PATH_H_

#include <stdint.h>

#include <deque>

#include "graph/scaffold_graph_vertex.h"


/**
 * @brief It is a path of scaffolds vertices in ScaffoldGraph.
 */
class ScaffoldGraphPath
{
public:
    ScaffoldGraphPath() {}
    ScaffoldGraphPath(const ScaffoldGraphPath &path)
        : vertices_(path.vertices_), distances_(path.distances_) {}

    const ScaffoldGraphPath &operator =(const ScaffoldGraphPath &path)
    { vertices_ = path.vertices_; distances_ = path.distances_; return *this; }

    ScaffoldGraphVertexAdaptor &operator [](uint32_t index) 
    { return vertices_[index]; }
    const ScaffoldGraphVertexAdaptor &operator [](uint32_t index) const 
    { return vertices_[index]; }

    void Append(const ScaffoldGraphVertexAdaptor &vertex, int d)
    {
        vertices_.push_back(vertex);
        if (vertices_.size() > 1)
            distances_.push_back(d);
    }

    void Pop()
    {
        vertices_.pop_back();
        if (!distances_.empty())
            distances_.pop_back();
    }

    const ScaffoldGraphPath &ReverseComplement()
    {
        std::reverse(vertices_.begin(), vertices_.end());
        for (unsigned i = 0; i < vertices_.size(); ++i)
            vertices_[i].ReverseComplement();
        std::reverse(distances_.begin(), distances_.end());
        return *this;
    }

    void Assemble(ContigGraphPath &path)
    {
        path.clear();

        for (unsigned i = 0; i < vertices_.size(); ++i)
        {
            if (i == 0)
                path.Append(vertices_[i].path(), 0);
            else
                path.Append(vertices_[i].path(), distances_[i-1]);
        }
    }

    void swap(ScaffoldGraphPath &path) 
    {
        if (this != &path)
        {
            vertices_.swap(path.vertices_); 
            distances_.swap(path.distances_); 
        }
    }

    ScaffoldGraphVertexAdaptor &front() { return vertices_.front(); }
    const ScaffoldGraphVertexAdaptor &front() const { return vertices_.front(); }

    ScaffoldGraphVertexAdaptor &back() { return vertices_.back(); }
    const ScaffoldGraphVertexAdaptor &back() const { return vertices_.back(); }

//    uint64_t kmer_count() const
//    {
//        uint64_t sum = 0;
//        for (unsigned i = 0; i < vertices_.size(); ++i)
//            sum += vertices_[i].kmer_count();
//        return sum;
//    }

    uint32_t size() const
    {
        uint32_t size = 0;
        for (unsigned i = 0; i < vertices_.size(); ++i)
            size += vertices_[i].path_size();
        for (unsigned i = 0; i < distances_.size(); ++i)
            size += distances_[i];
        return size;
    }

//    uint32_t internal_size(int kmer_size) const
//    {
//        if (vertices_.size() <= 1)
//            return vertices_.size();
//
//        uint32_t size = kmer_size + 1;
//        for (unsigned i = 1; i+1 < vertices_.size(); ++i)
//            size += vertices_[i].contig_size();
//        for (unsigned i = 0; i < distances_.size(); ++i)
//            size += distances_[i];
//        return size;
//    }

    uint32_t num_nodes() const { return vertices_.size(); }
    void clear() { vertices_.clear(); distances_.clear(); }

    std::deque<int> &distances() { return distances_; }
    const std::deque<int> &distances() const { return distances_; }

private:
    std::deque<ScaffoldGraphVertexAdaptor> vertices_;
    std::deque<int> distances_; 
};

namespace std
{
template <> inline void swap(ScaffoldGraphPath &x, ScaffoldGraphPath &y) { x.swap(y); }
}

#endif

