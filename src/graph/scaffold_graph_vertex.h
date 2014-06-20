/**
 * @file scaffold_graph_vertex.h
 * @brief ScaffoldGraphVertex Class and ScaffoldGraphVertexAdaptor Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.10
 * @date 2012-04-16
 */

#ifndef __GRAPH_SCAFFOLD_GRAPH_VERTEX_H_

#define __GRAPH_SCAFFOLD_GRAPH_VERTEX_H_

#include <stdint.h>

#include <algorithm>
#include <string>

#include "basic/kmer.h"
#include "graph/bit_edges.h"
#include "graph/vertex_status.h"
#include "sequence/sequence.h"
#include "graph/contig_info.h"
#include "graph/contig_graph_path.h"


/**
 * @brief It is the vertex class used in ScaffoldGraph class.
 */
class ScaffoldGraphVertex
{
public:
    explicit ScaffoldGraphVertex(const ContigGraphPath &path = ContigGraphPath())
        : path_(path) {}

    ScaffoldGraphVertex(const ScaffoldGraphVertex &x)
        : path_(x.path_), id_(x.id_), status_(x.status_) {}

    const ScaffoldGraphVertex &operator =(const ScaffoldGraphVertex &x)
    {
        if (this != &x)
        {
            path_ = x.path_;
            id_ = x.id_;
        }
        return *this;
    }

    const ContigGraphPath &path() const { return path_; }
    void set_path(const ContigGraphPath &path) { path_ = path; }

    uint32_t path_size() const { return path_.size(); }
//    const Sequence &contig() const { return contig_; }
//    void set_contig(const Sequence &contig) { contig_ = contig; }

//    uint32_t contig_size() const { return contig_.size(); }
//    uint32_t num_kmer() const { return contig_.size() - kmer_size() + 1; }

    uint32_t id() const { return id_; }
    void set_id(uint32_t id) { id_ = id; }

    VertexStatus &status() { return status_; }
    const VertexStatus &status() const { return status_; }

//    Kmer begin_kmer(int kmer_size) const { return contig_.GetKmer(0, kmer_size); }
//    Kmer end_kmer(int kmer_size) const { return contig_.GetKmer(contig_.size() - kmer_size, kmer_size); }
//
//    double coverage() const { return 1.0 * contig_info_.kmer_count() / (contig_size() - kmer_size() + 1); }
    double coverage() const
    {
        double sum = 0;
        int count = 0;
        for (unsigned i = 0; i < path_.num_nodes(); ++i)
        {
            sum += path_[i].kmer_count();
            count += path_[i].contig_size() - path_[i].kmer_size() + 1;
        }
        return sum / count;
    }

//    const SequenceCount &counts() const { return contig_info_.counts(); }
//    void set_counts(const SequenceCount &counts) { contig_info_.set_counts(counts); }

//    char get_base(uint32_t index) const { return contig_[index]; }

//    SequenceCountUnitType get_count(uint32_t index) const { return contig_info_.counts()[index]; }

    void swap(ScaffoldGraphVertex &x)
    { 
        if (this != &x)
        {
            path_.swap(x.path_); 
            std::swap(id_, x.id_); 
            status_.swap(x.status_); 
        }
    }

    void clear() 
    { 
        path_.clear(); 
        id_ = 0; 
        status_.clear(); 
    }

private:
    ContigGraphPath path_;
    uint32_t id_;
    VertexStatus status_;
};

/**
 * @brief It is a adaptor class used to access ScaffoldGraphVertex. Becase 
 * a scaffold vertex and its reverse complement share the same vertex, 
 * using adaptor makes sure that modification to the vertex consistant.
 */
class ScaffoldGraphVertexAdaptor
{
public:
    explicit ScaffoldGraphVertexAdaptor(ScaffoldGraphVertex *vertex = NULL, bool is_reverse = false)
    { vertex_ = vertex; is_reverse_ = is_reverse; }
    ScaffoldGraphVertexAdaptor(const ScaffoldGraphVertexAdaptor &x)
    { vertex_ = x.vertex_, is_reverse_ = x.is_reverse_; }

    const ScaffoldGraphVertexAdaptor &operator =(const ScaffoldGraphVertexAdaptor &x)
    { vertex_ = x.vertex_; is_reverse_ = x.is_reverse_; return *this; }
        
    bool operator <(const ScaffoldGraphVertexAdaptor &x) const
    { return (vertex_ != x.vertex_) ? (vertex_ < x.vertex_) : (is_reverse_ < x.is_reverse_); }
    bool operator >(const ScaffoldGraphVertexAdaptor &x) const
    { return (vertex_ != x.vertex_) ? (vertex_ > x.vertex_) : (is_reverse_ > x.is_reverse_); }

    bool operator ==(const ScaffoldGraphVertexAdaptor &x) const
    { return vertex_ == x.vertex_ && is_reverse_ == x.is_reverse_; }
    bool operator !=(const ScaffoldGraphVertexAdaptor &x) const
    { return vertex_ != x.vertex_ || is_reverse_ != x.is_reverse_; }

    const ScaffoldGraphVertexAdaptor &ReverseComplement() { is_reverse_ = !is_reverse_; return *this; }

    ContigGraphPath path() const
    {
        ContigGraphPath path = vertex_->path();
        return !is_reverse_ ? path : path.ReverseComplement();
    }

    uint32_t path_size() const { return vertex_->path_size(); }
//    Sequence contig() const
//    {
//        Sequence contig = vertex_->contig();
//        return !is_reverse_ ? contig : contig.ReverseComplement();
//    }

//    uint32_t contig_size() const { return vertex_->contig().size(); }
//    uint32_t num_kmer() const { return vertex_->num_kmer(); }

    void set_vertex(ScaffoldGraphVertex *vertex, bool is_reverse)
    { vertex_ = vertex; is_reverse_ = is_reverse; }

//    ContigInfo contig_info() const
//    {
//        ContigInfo contig_info = vertex_->contig_info();
//        return (!is_reverse_ ? contig_info : contig_info.ReverseComplement()); 
//    }
//
//    uint64_t kmer_size() const { return vertex_->kmer_size(); }
//    void set_kmer_size(uint64_t kmer_size) { vertex_->set_kmer_size(kmer_size); }
//
//    uint64_t kmer_count() const { return vertex_->kmer_count(); }
//    void set_kmer_count(uint64_t kmer_count) { vertex_->set_kmer_count(kmer_count); }

    uint32_t id() const { return vertex_->id(); }
    void set_id(uint32_t id) { vertex_->set_id(id); }

    VertexStatus &status() { return vertex_->status(); }
    const VertexStatus &status() const { return vertex_->status(); }

//    BitEdges &in_edges() { return !is_reverse_ ? vertex_->in_edges() : vertex_->out_edges(); }
//    const BitEdges &in_edges() const { return !is_reverse_ ? vertex_->in_edges() : vertex_->out_edges(); }
//
//    BitEdges &out_edges() { return !is_reverse_ ? vertex_->out_edges() : vertex_->in_edges(); }
//    const BitEdges &out_edges() const { return !is_reverse_ ? vertex_->out_edges() : vertex_->in_edges(); }

//    uint32_t &in_kmer_count() { return !is_reverse_ ? vertex_->in_kmer_count() : vertex_->out_kmer_count(); }
//    const uint32_t &in_kmer_count() const { return !is_reverse_ ? vertex_->in_kmer_count() : vertex_->out_kmer_count(); }
//
//    uint32_t &out_kmer_count() { return !is_reverse_ ? vertex_->out_kmer_count() : vertex_->out_kmer_count(); }
//    const uint32_t &out_kmer_count() const { return !is_reverse_ ? vertex_->out_kmer_count() : vertex_->out_kmer_count(); }

//    SequenceCount counts()
//    {
//        if (!is_reverse_) return vertex_->counts();
//        else { SequenceCount counts = vertex_->counts(); std::reverse(counts.begin(), counts.end()); return counts; }
//    }
//
//    char get_base(uint32_t index) const
//    { return (!is_reverse_) ? vertex_->get_base(index) : 3 - vertex_->get_base(contig_size() - 1 - index); }
//
//    SequenceCountUnitType get_count(uint32_t index) const 
//    { return (!is_reverse_) ? vertex_->get_count(index) : vertex_->get_count(vertex_->counts().size() - 1 - index); }
//
//    Kmer begin_kmer(int kmer_size) const 
//    { return !is_reverse_ ? vertex_->begin_kmer(kmer_size) : vertex_->end_kmer(kmer_size).ReverseComplement(); }
//
//    Kmer end_kmer(int kmer_size) const
//    { return !is_reverse_ ? vertex_->end_kmer(kmer_size) : vertex_->begin_kmer(kmer_size).ReverseComplement(); }
//
    double coverage() const 
    { return vertex_->coverage(); }

    bool is_reverse() const { return is_reverse_; }

    void swap(ScaffoldGraphVertexAdaptor &x)
    { 
        if (this != &x)
        {
            std::swap(vertex_, x.vertex_); 
            std::swap(is_reverse_, x.is_reverse_); 
        }
    }

    bool is_null() const { return vertex_ == NULL; }
    void clear() { vertex_->clear(); }

private:
    ScaffoldGraphVertex *vertex_;
    bool is_reverse_;
};

namespace std
{
template <> inline void swap(ScaffoldGraphVertex &x, ScaffoldGraphVertex &y) { x.swap(y); }
template <> inline void swap(ScaffoldGraphVertexAdaptor &x, ScaffoldGraphVertexAdaptor &y) { x.swap(y); }
}

#endif

