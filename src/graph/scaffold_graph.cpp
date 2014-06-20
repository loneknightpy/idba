/**
 * @file scaffold_graph.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.10
 * @date 2012-04-16
 */

#include "graph/scaffold_graph.h"

#include <cmath>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <sstream>

#include "basic/math.h"

using namespace std;

static bool Compare(const ScaffoldGraphEdgeAdaptor &x,
        const ScaffoldGraphEdgeAdaptor &y)
{ return x.distance() > y.distance(); }

void ScaffoldGraph::Initialize()
{
    deque<ContigGraphVertex> &contig_graph_vertices = contig_graph_.vertices();
    vertices_.resize(contig_graph_vertices.size());
    for (unsigned i = 0; i < contig_graph_vertices.size(); ++i)
    {
        ContigGraphPath path;
        path.Append(ContigGraphVertexAdaptor(&contig_graph_vertices[i]), 0);
        vertices_[i].set_path(path);
        vertices_[i].set_id(i);
        vertices_[i].status().clear();
    }
}

void ScaffoldGraph::Initialize(std::deque<ContigGraphPath> &paths)
{
    vertices_.resize(paths.size());
    for (unsigned i = 0; i < paths.size(); ++i)
    {
        vertices_[i].set_path(paths[i]);
        vertices_[i].set_id(i);
        vertices_[i].status().clear();
    }
}

void ScaffoldGraph::BuildContigToScaffoldMap()
{
    contig_to_scaffold_.clear();
    contig_to_scaffold_position_.clear();
    for (unsigned i = 0; i < vertices_.size(); ++i)
    {
        if (vertices_[i].status().IsDead())
            continue;

        for (int strand = 0; strand < 2; ++strand)
        {
            ScaffoldGraphVertexAdaptor current(&vertices_[i], strand);
            ContigGraphPath path = current.path();
            int d = 0;
            for (unsigned j = 0; j < path.num_nodes(); ++j)
            {
                if (contig_to_scaffold_.find(path[j]) != contig_to_scaffold_.end())
                {
                    contig_to_scaffold_[path[j]] = ScaffoldGraphVertexAdaptor();
                    contig_to_scaffold_position_[path[j]] = -1;
                }
                else
                {
                    contig_to_scaffold_[path[j]] = current;
                    contig_to_scaffold_position_[path[j]] = d;
                }

                if (j+1 < path.num_nodes())
                {
                    d += path[j].contig_size();
                    d += path.distances()[j];
                }
            }
        }
    }
}

void ScaffoldGraph::BuildEdges()
{
    BuildContigToScaffoldMap();

    in_edges_.clear();
    out_edges_.clear();
    in_edges_.resize(vertices_.size());
    out_edges_.resize(vertices_.size());
    edge_data_.clear();

    for (unsigned i = 0; i < pairs_.size(); ++i)
    {
        ContigGraphVertexAdaptor contig_from(&contig_graph_.vertices()[pairs_[i].from()>>1], pairs_[i].from()&1);
        ContigGraphVertexAdaptor contig_to(&contig_graph_.vertices()[pairs_[i].to()>>1], pairs_[i].to()&1);
        ScaffoldGraphVertexAdaptor from = contig_to_scaffold_[contig_from];
        ScaffoldGraphVertexAdaptor to = contig_to_scaffold_[contig_to];

        if (from.is_null() || to.is_null() || from.id() == to.id())
            continue;

        int level = pairs_[i].level();
        int distance = pairs_[i].distance() 
            - (from.path().size() - (contig_to_scaffold_position_[contig_from] + contig_from.contig_size()))
            - contig_to_scaffold_position_[contig_to];

        AddEdge(level, from, to, distance);
    }
}

void ScaffoldGraph::RefreshEdges()
{
    for (unsigned i = 0; i < vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ScaffoldGraphVertexAdaptor current(&vertices_[i], strand);
            deque<ScaffoldGraphEdgeAdaptor> &edges = GetEdges(current);
            int index = 0;
            for (unsigned j = 0; j < edges.size(); ++j)
            {
                if (edges[j].status().IsDead()
                        || edges[j].from().status().IsDead()
                        || edges[j].to().status().IsDead())
                {
                    edges[j].status().SetDeadFlag();
                }
                else
                {
                    edges[index++] = edges[j];
                }
            }
            edges.resize(index);
        }
    }
}

void ScaffoldGraph::ParseEdges(bool is_uneven)
{
    for (unsigned i = 0; i < vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ScaffoldGraphVertexAdaptor current(&vertices_[i], strand);
            deque<ScaffoldGraphEdgeAdaptor> &edges = GetEdges(current);
            for (unsigned j = 0; j < edges.size(); ++j)
            {
                if (edges[j].status().IsDead())
                    continue;

                edges[j].Parse();
                double expected = ExpectedEdges(edges[j].level(), 
                        edges[j].from().path().size(),
                        edges[j].to().path().size(),
                        edges[j].distance()
                        );

                int level = edges[j].level();
                if (is_uneven 
                        && edges[j].from().path().size() > mean_[level]*4
                        && edges[j].to().path().size() > mean_[level]*4
                        )
                {
                    double tmp = ExpectedEdges(edges[j].level(), 
                        edges[j].from().path().size(),
                        edges[j].to().path().size(),
                        edges[j].distance(),
                        //(edges[j].from().coverage() + edges[j].to().coverage()) / 2
                        //(edges[j].from().coverage() + edges[j].to().coverage()) / 2
                        //max(edges[j].from().coverage(), edges[j].to().coverage())
                        min(edges[j].from().coverage(), edges[j].to().coverage())
                        //min(edges[j].from().path().back().coverage(), edges[j].to().path().front().coverage())
                        );

                    //if (tmp > expected)
                        expected = tmp;
                }

                double rate = 0.3;
//                if (mean_[edges[j].level()] > 1000)
//                    rate = 0.1;

                if ((int)edges[j].values().size() < expected * rate
                        || edges[j].distance() > 0.75 * mean(edges[j].level())
                        || edges[j].distance() < -2*sd(edges[j].level()) - kmer_size()
                   )
                {
                    edges[j].status().SetDeadFlag();
                }
            }
        }
    }

    RefreshEdges();
}

void ScaffoldGraph::FilterEdges(int min_pairs, int min_length)
{
    for (unsigned i = 0; i < vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ScaffoldGraphVertexAdaptor current(&vertices_[i], strand);
            deque<ScaffoldGraphEdgeAdaptor> &edges = GetEdges(current);
            for (unsigned j = 0; j < edges.size(); ++j)
            {
                if (edges[j].status().IsDead())
                    continue;

                if ((int)edges[j].values().size() < min_pairs
                        || (int)edges[j].from().path_size() < min_length
                        || (int)edges[j].to().path_size() < min_length
                   )
                {
                    edges[j].status().SetDeadFlag();
                }
            }
        }
    }

    RefreshEdges();
}

void ScaffoldGraph::ClearStatus()
{
    for (unsigned i = 0; i < vertices().size(); ++i)
        vertices()[i].status().clear();
}

bool ScaffoldGraph::IsConnected(int level, ScaffoldGraphVertexAdaptor from, ScaffoldGraphVertexAdaptor to)
{
    deque<ScaffoldGraphVertexAdaptor> qu;
    qu.push_back(from);
    from.status().SetUsedFlag();

    int time = 0;
    int index = 0;
    bool is_found = false;
    while (++time < kTimeLimit && index < (int)qu.size() && !is_found)
    {
        ScaffoldGraphVertexAdaptor current = qu[index++];

        deque<ScaffoldGraphEdgeAdaptor> edges = GetEdges(level, current);
        for (unsigned i = 0; i < edges.size(); ++i)
        {
            if (edges[i].to() == to)
            {
                is_found = true;
                break;
            }
            else if (!edges[i].to().status().IsUsed())
            {
                qu.push_back(edges[i].to());
                edges[i].to().status().SetUsedFlag();
            }
        }
    }

    for (unsigned i = 0; i < qu.size(); ++i)
        qu[i].status().clear();

    return is_found;
}

int64_t ScaffoldGraph::RemoveTransitiveConnections(int level)
{
    int removed = 0;
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ScaffoldGraphVertexAdaptor current(&vertices()[i], strand);

            deque<ScaffoldGraphEdgeAdaptor> edges = GetEdges(level, current);
            if (edges.size() < 2U)
                continue;

            for (unsigned j = 0; j < edges.size(); ++j)
            {
                edges[j].status().SetDeadFlag();
                if (!IsConnected(level, current, edges[j].to()))
                    edges[j].status().ResetDeadFlag();
                else
                    ++removed;
            }
        }
    }

    return removed;
}

bool ScaffoldGraph::IsConsistent(int level, ScaffoldGraphVertexAdaptor current)
{
    deque<ScaffoldGraphEdgeAdaptor> edges = GetEdges(level, current);
    return edges.size() == 1;
}

bool ScaffoldGraph::IsConsistentMulti(int level, ScaffoldGraphVertexAdaptor current)
{
    deque<ScaffoldGraphEdgeAdaptor> edges = GetEdges(level, current);
    return edges.size() == 1;
}

bool ScaffoldGraph::ExtendPath(int level, ScaffoldGraphPath &scaffold_path)
{
    ScaffoldGraphVertexAdaptor current = scaffold_path.back();
    if (!IsConsistent(level, current))
        return false;

    deque<ScaffoldGraphEdgeAdaptor> edges = GetEdges(level, current);
    ScaffoldGraphVertexAdaptor next = edges[0].to();

    ScaffoldGraphVertexAdaptor rev_next = next;
    rev_next.ReverseComplement();
    if (!IsConsistent(level, rev_next))
        return false;

    if (next.status().IsDead())
        return false;

    if (next.status().IsUsed())
        return false;

    next.status().SetUsedFlag();
    scaffold_path.Append(next, edges[0].distance());
    return true;
}

bool ScaffoldGraph::ExtendPathMulti(int level, ScaffoldGraphPath &scaffold_path)
{
    ScaffoldGraphVertexAdaptor current = scaffold_path.back();

    deque<ScaffoldGraphEdgeAdaptor> edges = GetEdges(level, current);
    if (edges.size() == 0)
        return false;
    else
    {
        ScaffoldGraphVertexAdaptor next = edges[0].to();
        ScaffoldGraphVertexAdaptor rev_next = next;
        rev_next.ReverseComplement();

        if (IsConsistentMulti(level, current))
        {
            if (!IsConsistentMulti(level, rev_next))
                return false;

            if (next.status().IsDead())
                return false;

            if (next.status().IsUsed())
                return false;

            next.status().SetUsedFlag();
            scaffold_path.Append(next, edges[0].distance());
            return true;
        }
        else
        {
            set<ScaffoldGraphVertexAdaptor> table;
            for (unsigned i = 0; i < scaffold_path.num_nodes(); ++i)
            {
                ScaffoldGraphVertexAdaptor x = scaffold_path[i];
                for (int j = level+1; j < num_level(); ++j)
                {
                    deque<ScaffoldGraphEdgeAdaptor> long_edges = GetEdges(j, x);
                    for (unsigned k = 0; k < long_edges.size(); ++k)
                        table.insert(long_edges[k].to());
                }
            }

            int index = 0;
            for (unsigned i = 0; i < edges.size(); ++i)
            {
                if (table.find(edges[i].to()) != table.end())
                    edges[index++] = edges[i];
            }
            edges.resize(index);

            if (edges.size() == 1)
            {
                cout << "succeed" << endl;
                ScaffoldGraphVertexAdaptor next = edges[0].to();

                if (next.status().IsDead())
                    return false;

                if (next.status().IsUsed())
                    return false;

                next.status().SetUsedFlag();
                scaffold_path.Append(next, edges[0].distance());
                return true;
            }
            else
            {
                sort(edges.begin(), edges.end(), Compare);

                ScaffoldGraphVertexAdaptor next = edges[0].to();
                bool is_consistent = true;
                for (unsigned i = 1; i < edges.size(); ++i)
                {
                    deque<ScaffoldGraphEdgeAdaptor> middle_edges = GetEdges(level, edges[i].to());
                    is_consistent = false;
                    for (unsigned j = 0; j < middle_edges.size(); ++j)
                    {
                        if (middle_edges[j].to() == next)
                            is_consistent = true;
                    }

                    if (!is_consistent)
                        break;
                }

                if (!is_consistent)
                {
                    cout << "failed" << endl;
                    return false;
                }

                cout << "succeed2" << endl;

                if (next.status().IsDead())
                    return false;

                if (next.status().IsUsed())
                    return false;

                next.status().SetUsedFlag();
                scaffold_path.Append(next, edges[0].distance());
                return true;
            }

            return false;
        }
    }
}

int64_t ScaffoldGraph::Assemble(int level, deque<ContigGraphPath> &paths)
{
    paths.clear();

    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        if (vertices()[i].status().IsDead())
            continue;

        if (vertices()[i].status().IsUsed())
            continue;

        ScaffoldGraphVertexAdaptor current(&vertices()[i]);
        ScaffoldGraphPath scaffold_path;
        scaffold_path.Append(current, 0);
        current.status().SetUsedFlag();

        for (int strand = 0; strand < 2; ++strand)
        {
            while (ExtendPath(level, scaffold_path))
                ;
            scaffold_path.ReverseComplement();
        }

        ContigGraphPath path;
        scaffold_path.Assemble(path);
        paths.push_back(path);
    }

    ClearStatus();
    return paths.size();
}

int64_t ScaffoldGraph::AssembleMulti(int level, deque<ContigGraphPath> &paths)
{
    paths.clear();

    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        if (vertices()[i].status().IsDead())
            continue;

        if (vertices()[i].status().IsUsed())
            continue;

        ScaffoldGraphVertexAdaptor current(&vertices()[i]);
        ScaffoldGraphPath scaffold_path;
        scaffold_path.Append(current, 0);
        current.status().SetUsedFlag();

        for (int strand = 0; strand < 2; ++strand)
        {
            while (ExtendPathMulti(level, scaffold_path))
                ;
            scaffold_path.ReverseComplement();
        }

        ContigGraphPath path;
        scaffold_path.Assemble(path);
        paths.push_back(path);
    }

    ClearStatus();
    return paths.size();
}

int64_t ScaffoldGraph::Assemble(int level, deque<Sequence> &contigs)
{
    deque<ContigGraphPath> paths;
    Assemble(level, paths);
    contigs.resize(paths.size());
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        ContigInfo contig_info;
        paths[i].Assemble(contigs[i], contig_info);
    }
    return contigs.size();
}

int64_t ScaffoldGraph::AssembleMulti(int level, deque<Sequence> &contigs)
{
    deque<ContigGraphPath> paths;
    AssembleMulti(level, paths);
    contigs.resize(paths.size());
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        ContigInfo contig_info;
        paths[i].Assemble(contigs[i], contig_info);
    }
    return contigs.size();
}

double ScaffoldGraph::ExpectedEdges(int level, int len1, int len2, int distance)
{
    //cout << len1 << " " << len2 << " " << distance << endl;

    int read_length = read_length_[level];
    double mean = mean_[level];
    double sd = sd_[level];
    double expected_coverage = expected_coverage_[level];

    int from = max(0, int(len1 - mean - 2*sd));
    int to = max(0, len1 - read_length);

    from = from - len1 - distance;
    to = to - len1 - distance;

    double sum = 0;
    for (int i = from; i < to; ++i)
    {
        double expected = i + mean;
        double left = read_length - expected;
        double right = min(len2, int(mean + 2*sd)) - expected;
        sum += NormalCDF(right/sd) - NormalCDF(left/sd);
    }

    double p = sum / (to - from);

    return p * (to - from) * expected_coverage/2;
}

double ScaffoldGraph::ExpectedEdges(int level, int len1, int len2, int distance, double expected_coverage)
{
    //cout << len1 << " " << len2 << " " << distance << endl;

    int read_length = read_length_[level];
    double mean = mean_[level];
    double sd = sd_[level];
    //expected_coverage = min(expected_coverage_[level], expected_coverage);
    //double expected_coverage = expected_coverage_[level];

    int from = max(0, int(len1 - mean - 2*sd));
    int to = max(0, len1 - read_length);

    from = from - len1 - distance;
    to = to - len1 - distance;

    double sum = 0;
    for (int i = from; i < to; ++i)
    {
        double expected = i + mean;
        double left = read_length - expected;
        double right = min(len2, int(mean + 2*sd)) - expected;
        sum += NormalCDF(right/sd) - NormalCDF(left/sd);
    }

    double p = sum / (to - from);

    return p * (to - from) * expected_coverage/2;
}



