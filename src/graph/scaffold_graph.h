/**
 * @file scaffold_graph.h
 * @brief ScaffoldGraph Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.1
 * @date 2011-09-02
 */

#ifndef __GRAPH_SCAFFOLD_GRAPH_H_

#define __GRAPH_SCAFFOLD_GRAPH_H_

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "graph/contig_graph.h"
#include "graph/scaffold_graph_path.h"
#include "graph/scaffold_graph_vertex.h"


/**
 * @brief It is a class for storing position of a pair of reads.
 */
class ScaffoldGraphPair
{
public:
    ScaffoldGraphPair() {}
    ScaffoldGraphPair(int level, int from, int to, int distance)
        : level_(level), from_(from), to_(to), distance_(distance) {}

    int level() { return level_; }
    int from() { return from_; }
    int to() { return to_; }
    int distance() { return distance_; }

private:
    int level_;
    int from_;
    int to_;
    int distance_;
};

/**
 * @brief It is the edge in ScaffoldGraph.
 */
class ScaffoldGraphEdge
{
public:
    ScaffoldGraphEdge() {}

    ScaffoldGraphEdge(int level, ScaffoldGraphVertexAdaptor from, 
            ScaffoldGraphVertexAdaptor to, int d)
    { from_ = from; to_ = to; values_.push_back(d); level_ = level; }

    ScaffoldGraphVertexAdaptor from() const { return from_; }
    ScaffoldGraphVertexAdaptor to() const { return to_; }

    std::vector<int> &values() { return values_; }
    const std::vector<int> &values() const { return values_; }

    int distance() const { return distance_; }
    int level() const { return level_; }

    VertexStatus &status() { return status_; }
    const VertexStatus &status() const { return status_; }

    void Parse()
    {
        std::sort(values_.begin(), values_.end());
        if (values_.empty())
            distance_ = (1 << 30);
        else
            distance_ = values_[values_.size()/2];
    }

private:
    ScaffoldGraphVertexAdaptor from_;
    ScaffoldGraphVertexAdaptor to_;
    std::vector<int> values_;
    int distance_;
    int level_;
    VertexStatus status_;
};

/**
 * @brief It is an adaptor class for access ScaffoldGraphEdge. Because the edge
 * an its reverse complement share the same edge instance, using adaptor makes
 * sure the modification of edge consistant.
 */
class ScaffoldGraphEdgeAdaptor
{
public:
    explicit ScaffoldGraphEdgeAdaptor(ScaffoldGraphEdge *edge = NULL, bool is_reverse = false)
    { edge_ = edge; is_reverse_ = is_reverse; }
    ScaffoldGraphEdgeAdaptor(const ScaffoldGraphEdgeAdaptor &x)
    { edge_ = x.edge_, is_reverse_ = x.is_reverse_; }

    ScaffoldGraphVertexAdaptor from() const 
    { return !is_reverse_ ? edge_->from() : edge_->to().ReverseComplement(); }

    ScaffoldGraphVertexAdaptor to() const 
    { return !is_reverse_ ? edge_->to() : edge_->from().ReverseComplement(); }

    std::vector<int> &values() { return edge_->values(); }
    const std::vector<int> &values() const { return edge_->values(); }

    int distance() const { return edge_->distance(); }
    int level() const { return edge_->level(); }

    VertexStatus &status() { return edge_->status(); }
    const VertexStatus &status() const { return edge_->status(); }

    const ScaffoldGraphEdgeAdaptor &ReverseComplement() 
    { is_reverse_ = !is_reverse_; return *this; }

    void Parse() { edge_->Parse(); }

private:
    ScaffoldGraphEdge *edge_;
    bool is_reverse_;
};

/**
 * @brief It is a contig graph built upon paired-end reads information. Each 
 * vertex is a contg, maybe a path of contigs in ContigGraph, and each edge 
 * between vertex u and vertex v means there are at least min_pairs paired-end
 * reads connecting u and v.
 */
class ScaffoldGraph
{
public:
    explicit ScaffoldGraph(uint32_t kmer_size = 0)
        : contig_graph_(kmer_size), min_pairs_(5)
    { Initialize(); }

    explicit ScaffoldGraph(uint32_t kmer_size, const std::deque<Sequence> &contigs)
        : contig_graph_(kmer_size, contigs), min_pairs_(5)
    { Initialize(); }

    explicit ScaffoldGraph(uint32_t kmer_size, const std::deque<Sequence> &contigs,
            const std::deque<ContigInfo> &contig_infos)
        : contig_graph_(kmer_size, contigs, contig_infos), min_pairs_(5)
    { Initialize(); }

    ~ScaffoldGraph() { clear(); }

    void Initialize();
    void Initialize(std::deque<ContigGraphPath> &paths);

    void BuildContigToScaffoldMap();
    void BuildEdges();
    void RefreshEdges();
    void ParseEdges(bool is_uneven = false);
    void FilterEdges(int min_pairs, int min_length);

    void ClearStatus();

    bool IsConnected(int level, ScaffoldGraphVertexAdaptor from, ScaffoldGraphVertexAdaptor to);
    int64_t RemoveTransitiveConnections(int level);

    bool IsConsistent(int level, ScaffoldGraphVertexAdaptor current);
    bool IsConsistentMulti(int level, ScaffoldGraphVertexAdaptor current);

    bool ExtendPath(int level, ScaffoldGraphPath &scaffold_path);
    bool ExtendPathMulti(int level, ScaffoldGraphPath &scaffold_path);

    int64_t Assemble(int level, std::deque<ContigGraphPath> &paths);
    int64_t Assemble(int level, std::deque<Sequence> &contigs);

    int64_t AssembleMulti(int level, std::deque<ContigGraphPath> &paths);
    int64_t AssembleMulti(int level, std::deque<Sequence> &contigs);

    void AddPair(int level, int from, int to, int d)
    { pairs_.push_back(ScaffoldGraphPair(level, from, to, d)); }

    void AddEdge(int level, ScaffoldGraphVertexAdaptor from, ScaffoldGraphVertexAdaptor to, int d)
    {
        std::deque<ScaffoldGraphEdgeAdaptor> &all_edges = GetEdges(from);
        for (unsigned i = 0; i < all_edges.size(); ++i)
        {
            if (all_edges[i].level() == level && all_edges[i].to() == to)
            {
                all_edges[i].values().push_back(d);
                return;
            }
        }
        AddNewEdge(level, from, to, d);
    }

    void AddNewEdge(int level, ScaffoldGraphVertexAdaptor from, ScaffoldGraphVertexAdaptor to, int d)
    {
        ScaffoldGraphEdge edge(level, from, to, d);
        edge_data_.push_back(edge);

        ScaffoldGraphEdgeAdaptor adp(&edge_data_.back());
        GetEdges(adp.from()).push_back(adp);

        adp.ReverseComplement();
        GetEdges(adp.from()).push_back(adp);
    }

    std::deque<ScaffoldGraphEdgeAdaptor> GetEdges(int level, ScaffoldGraphVertexAdaptor current)
    {
        std::deque<ScaffoldGraphEdgeAdaptor> edges;
        std::deque<ScaffoldGraphEdgeAdaptor> &all_edges = GetEdges(current);
        for (unsigned i = 0; i < all_edges.size(); ++i)
        {
            if (all_edges[i].level() == level
                    && !all_edges[i].status().IsDead())
                edges.push_back(all_edges[i]);
        }
        return edges;
    }

    std::deque<ScaffoldGraphEdgeAdaptor> &GetEdges(ScaffoldGraphVertexAdaptor current)
    { return !current.is_reverse() ? out_edges_[current.id()] : in_edges_[current.id()]; }

    std::deque<ScaffoldGraphVertex> &vertices() { return vertices_; }
    const std::deque<ScaffoldGraphVertex> &vertices() const { return vertices_; }

    ContigGraph &contig_graph() { return contig_graph_; }
    const ContigGraph &contig_graph() const { return contig_graph_; }

    int num_edges(int level) const
    {
        int count = 0;
        for (unsigned i = 0; i < edge_data_.size(); ++i)
        {
            if (edge_data_[i].status().IsDead())
                continue;

            if (edge_data_[i].level() != level)
                continue;

            ++count;
        }

        return count;
    }

    int min_pairs() const { return min_pairs_; }
    void set_min_pairs(int min_pairs) { min_pairs_ = min_pairs; }

    int kmer_size() const { return contig_graph_.kmer_size(); }

    int read_length(int level) const { return read_length_[level]; }
    double expected_coverage(int level) const { return expected_coverage_[level]; }
    double mean(int level) const { return mean_[level]; }
    double sd(int level) const { return sd_[level]; }
    int num_level() const { return read_length_.size(); }

    void set_library_info(int level, int read_length, double coverage, double mean, double sd)
    {
        if ((int)read_length_.size() < level + 1)
        {
            read_length_.resize(level+1);
            expected_coverage_.resize(level+1);
            mean_.resize(level+1);
            sd_.resize(level+1);
        }
        read_length_[level] = read_length;
        expected_coverage_[level] = coverage;
        mean_[level] = mean;
        sd_[level] = sd;
    }

    void clear()
    { contig_graph_.clear(); }

private:
    ScaffoldGraph(const ScaffoldGraph &);
    const ScaffoldGraph &operator =(const ScaffoldGraph &);

    double ExpectedEdges(int level, int len1, int len2, int distance);
    double ExpectedEdges(int level, int len1, int len2, int distance, double expected_coverage);

    static const int kTimeLimit = 500;

    ContigGraph contig_graph_;
    std::deque<ScaffoldGraphVertex> vertices_;

    std::map<ContigGraphVertexAdaptor, ScaffoldGraphVertexAdaptor> contig_to_scaffold_;
    std::map<ContigGraphVertexAdaptor, int> contig_to_scaffold_position_;
    std::deque<std::deque<ScaffoldGraphEdgeAdaptor> > in_edges_;
    std::deque<std::deque<ScaffoldGraphEdgeAdaptor> > out_edges_;
    std::deque<ScaffoldGraphEdge> edge_data_;
    std::deque<ScaffoldGraphPair> pairs_;
    
    int min_pairs_;
    std::vector<int> read_length_;
    std::vector<double> expected_coverage_;
    std::vector<double> mean_;
    std::vector<double> sd_;
};

#endif

