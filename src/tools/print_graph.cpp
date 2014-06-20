/**
 * @file print_graph.cpp
 * @brief Print a contig graph.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.4
 * @date 2011-09-21
 */

#include "graph/contig_graph.h"
#include "graph/hash_graph.h"
#include "misc/options_description.h"
#include "sequence/sequence.h"
#include "sequence/sequence_io.h"

#include <algorithm>
#include <cstdio>
#include <deque>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    int kmer_size = 50;
    int max_length = 1000000;
    
    OptionsDescription desc;
    desc.AddOption("kmer", "k", kmer_size, "k value");
    desc.AddOption("max_length", "", max_length, "max length");

    desc.Parse(argc, argv);

    deque<Sequence> refs;
    ReadSequence(argv[1], refs);

    HashGraph hash_graph(kmer_size);
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        if ((int)refs[i].size() > max_length)
            refs[i].resize(max_length);
        hash_graph.InsertKmers(refs[i]);
    }

    hash_graph.Refresh();
    hash_graph.AddAllEdges();

    deque<Sequence> contigs;
    deque<ContigInfo> contig_infos;
    hash_graph.Assemble(contigs, contig_infos);

    cerr << "build" << endl;

    ContigGraph contig_graph(kmer_size);
    contig_graph.Initialize(contigs, contig_infos);

    cerr << "kmer " << hash_graph.num_vertices() << " branches " << contigs.size()<< endl;

    deque<deque<ContigGraphVertexAdaptor> > components;
    deque<string> component_strings;
    contig_graph.GetComponents(components, component_strings);

    for (unsigned i = 0; i < component_strings.size(); ++i)
        cout << component_strings[i] << endl;

    //FastaWriter writer(argv[2]);
    WriteSequence(argv[2], contigs, "conitg");

    return 0;
}

