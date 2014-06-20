#include <iostream>
#include <algorithm>
#include <cstdio>
#include <deque>

#include "graph/contig_graph.h";
#include "sequence/sequence.h";

using namespace std;

int main(int argc, char *argv[])
{
    int kmer_size = 25;

    ContigGraph contig_graph(kmer_size);
    deque<Sequence> contigs;

    contig_graph.Initialize(contigs);

    return 0;
}


