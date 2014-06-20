/**
 * @file idba_tran.cpp
 * @brief An iterative de Bruijn graph assembler for transcriptome sequencing reads.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.13
 * @date 2012-09-14
 */

#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

#include "assembly/assembly_utility.h"
#include "assembly/local_assembler.h"
#include "basic/bit_operation.h"
#include "basic/histgram.h"
#include "graph/contig_graph.h"
#include "graph/hash_graph.h"
#include "misc/hash_aligner.h"
#include "misc/log.h"
#include "misc/options_description.h"
#include "misc/utils.h"
#include "sequence/sequence.h"
#include "sequence/sequence_io.h"
#include "sequence/short_sequence.h"


using namespace std;

struct IDBAOption
{
    string directory;
    string read_file;
    string long_read_file;
    int mink;
    int maxk;
    int step;
    int inner_mink;
    int inner_step;
    int prefix_length;
    int min_count;
    int min_support;
    int min_contig;
    double similar;
    int max_mismatch;
    int seed_kmer_size;
    int num_threads;
    int min_pairs;
    int max_gap;
    bool is_no_local;
    bool is_use_all;
    bool is_no_coverage;
    bool is_no_correct;
    bool is_pre_correction;
    bool is_no_internal;
    string reference;

    IDBAOption()
    {
        directory = "out";
        mink = 20;
        maxk = 60;
        step = 10;
        inner_mink = 10;
        inner_step = 5;
        prefix_length = 3;
        min_count = 2;
        min_support = 1;
        min_contig = 200;
        similar = 0.95;
        max_mismatch = 3;
        seed_kmer_size = 30;
        num_threads = 0;
        min_pairs = 3;
        max_gap = 0;
        is_no_local = false;
        is_use_all = false;
        is_no_coverage = false;
        is_no_correct = false;
        is_pre_correction = false;
    }

    string log_file()
    { return directory + "/log"; }

    string kmer_file()
    { return directory + "/kmer"; }

    string align_file(int kmer_size)
    { return directory + FormatString("/align-%d", kmer_size); }

    string align_local_file(int kmer_size)
    { return directory + FormatString("/align-local-%d", kmer_size); }

    string graph_file(int kmer_size)
    { return directory + FormatString("/graph-%d.fa", kmer_size); }

    string contig_file(int kmer_size)
    { return directory + FormatString("/contig-%d.fa", kmer_size); }

    string contig_info_file(int kmer_size)
    { return directory + FormatString("/contig-info-%d.fa", kmer_size); }

    string local_contig_file(int kmer_size)
    { return directory + FormatString("/local-contig-%d.fa", kmer_size); }

    string contig_file()
    { return directory + "/contig.fa"; }

    string scaffold_file()
    { return directory + "/scaffold.fa"; }

    string ref_contig_file()
    { return directory + "/ref_contig.fa"; }

    string transcript_file(int kmer_size)
    { return directory + FormatString("/transcript-%d.fa", kmer_size); }

    string label_transcript_file(int kmer_size)
    { return directory + FormatString("/label-transcript-%d.fa", kmer_size); }

    string component_file(int kmer_size)
    { return directory + FormatString("/component-%d", kmer_size); }
};

AssemblyInfo assembly_info;
IDBAOption option;
double median = 0;
double sd = 0;
int read_length = 0;
int max_isoforms = 3;
int max_component_size = 30;

void BuildHashGraph(int kmer_size);
void Assemble(HashGraph &hash_graph);
void AlignReads(const string &contig_file, const string &align_file);
void CorrectReads(int kmer_size);
void LocalAssembly(int kmer_size, int new_kmer_size);
void Iterate(int kmer_size, int new_kmer_size);
void FindIsoforms(ContigGraph &contig_graph, deque<Sequence> &transcripts);

bool CompareCoverage(const ContigGraphVertexAdaptor &x, const ContigGraphVertexAdaptor &y)
{
    return x.coverage() > y.coverage();
}

bool CompareCoveragePath(const ContigGraphPath &x, const ContigGraphPath &y)
{
    return CompareCoverage(x.back(), y.back());
}

double Similarity(ContigGraphPath &x, ContigGraphPath &y)
{
    deque<ContigGraphVertexAdaptor> a;
    deque<ContigGraphVertexAdaptor> b;

    for (unsigned i = 0; i < x.num_nodes(); ++i)
        a.push_back(x[i]);
    for (unsigned i = 0; i < y.num_nodes(); ++i)
        b.push_back(y[i]);

    sort(a.begin(), a.end());
    sort(b.begin(), b.end());

    unsigned i = 0;
    unsigned j = 0;
    double common = 0;
    while (i < a.size() && j < b.size())
    {
        if (a[i] < b[j])
            ++i;
        else if (a[i] > b[j])
            ++j;
        else
        {
            common += a[i].contig_size();
            ++i;
            ++j;
        }
    }

    return common / max(x.size(), y.size());
}

void FindIsoforms(ContigGraphPath &path, deque<ContigGraphPath> &isoforms, ContigGraph &contig_graph, deque<ContigGraphPath> &component_isoforms)
{
    if ((int)isoforms.size() >= max_isoforms)
        return;

    int kmer_size = contig_graph.kmer_size();

    ContigGraphVertexAdaptor current = path.back();
    if (current.out_edges().size() != 0)
    {
        deque<ContigGraphVertexAdaptor> candidates;
        for (int x = 0; x < 4; ++x)
        {
            if (current.out_edges()[x])
            {
                ContigGraphVertexAdaptor next = contig_graph.GetNeighbor(current, x);

                if (next.status().IsUsed())
                    continue;

                candidates.push_back(next);
            }
        }

        sort(candidates.begin(), candidates.end(), CompareCoverage);

        for (unsigned i = 0; i < candidates.size(); ++i)
        {
            ContigGraphVertexAdaptor next = candidates[i];
            path.Append(next, - kmer_size + 1);
            next.status().SetUsedFlag();
            FindIsoforms(path, isoforms, contig_graph, component_isoforms);
            next.status().ResetUsedFlag();
            path.Pop();
        }
    }
    else
    {
        bool flag = true;

        if (flag && path.size() > 300)
        {
            isoforms.push_back(path);
        }
    }
}

bool Compare(int x, int y)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    if (reads[x] != reads[y])
        return reads[x] < reads[y];
    else
        return reads[x+1] < reads[y+1];
}

void SortReads(AssemblyInfo &assembly_info)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    deque<int> aux;
    for (unsigned i = 0; i < reads.size(); i += 2)
    {
        if (reads[i+1] < reads[i])
            swap(reads[i], reads[i+1]);
        aux.push_back(i);
    }

    sort(aux.begin(), aux.end(), Compare);

    deque<ShortSequence> new_reads;
    int last = -1;
    for (int i = 0; i < (int64_t)aux.size(); ++i)
    {
        int id = aux[i];
        if (last == -1 || reads[id] != reads[last] || reads[id+1] != reads[last+1])
        {
            new_reads.push_back(reads[id]);
            new_reads.push_back(reads[id+1]);
        }

        last = id;
    }

    reads.swap(new_reads);
}

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    
    desc.AddOption("out", "o", option.directory, "output directory");
    desc.AddOption("read", "r", option.read_file, FormatString("fasta read file (<=%d)", ShortSequence::max_size()));
    desc.AddOption("long_read", "l", option.long_read_file, FormatString("fasta long read file (>%d)", ShortSequence::max_size()));
    desc.AddOption("mink", "", option.mink, FormatString("minimum k value (<=%d)", Kmer::max_size()));
    desc.AddOption("maxk", "", option.maxk, FormatString("maximum k value (<=%d)", Kmer::max_size()));
    desc.AddOption("step", "", option.step, "increment of k-mer of each iteration");
    desc.AddOption("inner_mink", "", option.inner_mink, "inner minimum k value");
    desc.AddOption("inner_step", "", option.inner_step, "inner increment of k-mer");
    desc.AddOption("prefix", "", option.prefix_length, "prefix length used to build sub k-mer table");
    desc.AddOption("min_count", "", option.min_count, "minimum multiplicity for filtering k-mer when building the graph");
    desc.AddOption("min_support", "", option.min_support, "minimum supoort in each iteration");
    desc.AddOption("num_threads", "", option.num_threads, "number of threads");
    desc.AddOption("seed_kmer", "", option.seed_kmer_size, "seed kmer size for alignment");
    desc.AddOption("min_contig", "", option.min_contig, "minimum size of contig");
    desc.AddOption("similar", "", option.similar, "similarity for alignment");
    desc.AddOption("max_mismatch", "", option.max_mismatch, "max mismatch of error correction");
    //desc.AddOption("min_pairs", "", option.min_pairs, "minimum number of pairs");
    //desc.AddOption("max_gap", "", option.max_gap, "maximum gap in reference");
    desc.AddOption("no_local", "", option.is_no_local, "do not use local assembly");
    desc.AddOption("no_coverage", "", option.is_no_coverage, "do not iterate on coverage");
    desc.AddOption("no_correct", "", option.is_no_correct, "do not do correction");
    desc.AddOption("pre_correction", "", option.is_pre_correction, "perform pre-correction before assembly");
    desc.AddOption("max_isoforms", "", max_isoforms, "maximum number of isoforms");
    desc.AddOption("max_component_size", "", max_component_size, "maximum size of components");

    try
    {
        desc.Parse(argc, argv);

        if (option.read_file == "" && option.long_read_file == "")
            throw logic_error("not enough parameters");

        if (option.maxk < option.mink)
            throw invalid_argument("mink is larger than maxk");

        if (option.maxk > (int)Kmer::max_size())
            throw invalid_argument("maxk is too large");
    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "IDBA-Tran - Iterative de Bruijn Graph Assembler for next-generation transcriptome sequencing data." << endl;
        cerr << "Usage: idba_tran -r read.fa -o output_dir" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    MakeDir(option.directory);

    LogThread log_thread(option.log_file());

    string begin_file = option.directory + "/begin";
    fclose(OpenFile(begin_file, "wb"));

    if (option.num_threads == 0)
        option.num_threads = omp_get_max_threads();
    else
        omp_set_num_threads(option.num_threads);
    cout << "number of threads " << option.num_threads << endl;

    ReadInput(option.read_file, option.long_read_file, assembly_info);
    cout << "reads " << assembly_info.reads.size() << endl;
    cout << "long reads " << assembly_info.long_reads.size() << endl;

    SortReads(assembly_info);
    cout << "reads " << assembly_info.reads.size() << endl;
    cout << "long reads " << assembly_info.long_reads.size() << endl;

    read_length = assembly_info.read_length();
    cout << "read_length " << read_length << endl;

    if (option.is_pre_correction)
    {
        int kmer_size = (option.maxk + option.mink)/2;
        cout << "kmer " << kmer_size << endl;
        BuildHashGraph(kmer_size);
        AlignReads(option.contig_file(kmer_size), option.align_file(kmer_size));
        CorrectReads(kmer_size);
        assembly_info.ClearStatus();
    }

    int old_kmer_size = 0;
    int kmer_size = option.mink;
    while (true)
    {
        cout << "kmer " << kmer_size << endl;

        if (kmer_size == option.mink)
            BuildHashGraph(kmer_size);
        else
            Iterate(old_kmer_size, kmer_size);

        assembly_info.ClearStatus();
        AlignReads(option.contig_file(kmer_size), option.align_file(kmer_size));

        CorrectReads(kmer_size);
        assembly_info.ClearStatus();

        old_kmer_size = kmer_size;
        kmer_size = min(option.maxk, kmer_size + option.step);
        LocalAssembly(old_kmer_size, kmer_size);

        if (old_kmer_size == option.maxk)
            break;
    }

    kmer_size = option.maxk;

    deque<Sequence> contigs;
    deque<string> names;
    ReadSequence(option.contig_file(kmer_size), contigs, names);
    FastaWriter writer(option.contig_file());
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() >= option.min_contig)
            writer.Write(contigs[i], names[i]);
    }

    string end_file = option.directory + "/end";
    fclose(OpenFile(end_file, "wb"));

    fflush(stdout);

    return 0;
}

void BuildHashGraph(int kmer_size)
{
    BuildKmerFile(assembly_info, kmer_size, option.min_count, option.prefix_length, option.kmer_file());

    HashGraph hash_graph(kmer_size);
    ReadKmerFile(option.kmer_file(), hash_graph);

    hash_graph.RefreshEdges();

    if (!option.is_no_internal)
        InsertInternalKmers(assembly_info, hash_graph, option.min_count);

    if (option.reference != "")
    {
        deque<Sequence> ref_contigs;
        ReadSequence(option.ref_contig_file(), ref_contigs);
#pragma omp parallel for
        for (int64_t i = 0; i < (int64_t)ref_contigs.size(); ++i)
            hash_graph.InsertUncountKmers(ref_contigs[i]);
        hash_graph.RefreshEdges();
    }

    Assemble(hash_graph);
}

void Assemble(HashGraph &hash_graph)
{
    cout << "kmers " << hash_graph.num_vertices() << " "<< hash_graph.num_edges() << endl;

    int kmer_size = hash_graph.kmer_size();
    double min_cover = max(1, (kmer_size == option.mink ? option.min_count : option.min_support));

    Histgram<int> hist = hash_graph.coverage_histgram();
    //double expected_coverage = hist.mean();

    deque<Sequence> contigs;
    deque<ContigInfo> contig_infos;
    hash_graph.Assemble(contigs, contig_infos);
    hash_graph.clear();

    {
        HashGraph tmp_hash_graph;
        tmp_hash_graph.swap(hash_graph);
    }

    ContigGraph contig_graph(kmer_size, contigs, contig_infos);
    contigs.clear();
    contig_infos.clear();

    if (!option.is_no_coverage)
    {
        contig_graph.RemoveStandAlone(kmer_size);

        int bubble = contig_graph.RemoveBubble();
        cout << "merge bubble " << bubble << endl;

        contig_graph.RemoveLocalLowCoverage(min_cover, option.min_contig, 0.1);
    }

    contig_graph.SortVertices();
    contig_graph.GetContigs(contigs, contig_infos);
    WriteSequence(option.graph_file(kmer_size), contigs);
    contigs.clear();
    contig_infos.clear();

    if (!option.is_no_coverage)
    {
        double ratio = 0.25;

        deque<Sequence> multi_contigs;
        deque<ContigInfo> multi_contig_infos;
        contig_graph.GetContigs(multi_contigs, multi_contig_infos);
        PrintN50(multi_contigs);

        contig_graph.Trim(10);
        contig_graph.MergeSimilarPath();
        contig_graph.GetContigs(multi_contigs, multi_contig_infos);

        contig_graph.InitializeTable();
        contig_graph.IterateComponentCoverage2(option.min_contig, ratio, min_cover, 1e100, 1.1, max_component_size);
        contig_graph.GetContigs(multi_contigs, multi_contig_infos);

        contig_graph.Trim(10);
        contig_graph.Prune(kmer_size);
        contig_graph.GetContigs(multi_contigs, multi_contig_infos);

        contig_graph.MergeSimilarPath();
    }

    deque<Sequence> multi_contigs;
    deque<ContigInfo> multi_contig_infos;
    contig_graph.SortVertices();
    contig_graph.GetContigs(multi_contigs, multi_contig_infos);
    PrintN50(multi_contigs);
    WriteSequence(option.contig_file(kmer_size), multi_contigs);
    WriteContigInfo(option.contig_info_file(kmer_size), multi_contig_infos);

    deque<Sequence> transcripts;

    FindIsoforms(contig_graph, transcripts);

    int index = 0;
    for (unsigned i = 0; i < transcripts.size(); ++i)
    {
        if (transcripts[i].size() >= 300)
            transcripts[index++] = transcripts[i];
    }
    transcripts.resize(index);

    PrintN50(transcripts);
    WriteSequence(option.transcript_file(kmer_size), transcripts, FormatString("transcript-%d", kmer_size));
}

void AlignReads(const string &contig_file, const string &align_file)
{
    deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);

    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, align_file, true);
    cout << "aligned " << num_aligned_reads << " reads" << endl;
}

void CorrectReads(int kmer_size)
{
    if (option.is_no_correct)
        return;

    deque<Sequence> contigs;
    deque<string> names;
    deque<ContigInfo> contig_infos;
    ReadSequence(option.contig_file(kmer_size), contigs, names);
    CorrectReads(assembly_info, contigs, contig_infos, option.align_file(kmer_size), option.max_mismatch);
    //WriteSequence(option.contig_file(kmer_size), contigs);
    WriteContig(option.contig_file(kmer_size), contigs, contig_infos, FormatString("contig-%d", kmer_size));
}

void LocalAssembly(int kmer_size, int new_kmer_size)
{
    EstimateDistance(option.align_file(kmer_size), median, sd);
    if (median < 0 || median != median || sd != sd || sd > 2*median)
    {
        cout << "invalid insert distance" << endl;
        deque<Sequence> local_contigs;
        WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
        return;
    }

    deque<ShortSequence> &reads = assembly_info.reads;

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    LocalAssembler local_assembler;
    local_assembler.Initialize(assembly_info, contigs);
    local_assembler.set_num_threads(option.num_threads);
    local_assembler.set_mink(option.inner_mink);
    local_assembler.set_maxk(new_kmer_size);
    local_assembler.set_step(option.inner_step);
    local_assembler.set_min_contig(option.min_contig);
    local_assembler.set_insert_distance(median, sd);

    FILE *falign = OpenFile(option.align_file(kmer_size), "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);

        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0)
                local_assembler.AddReadByHashAlignerRecord(record, offset + i);
        }
    }
    fclose(falign);

    deque<Sequence> local_contigs;

    if (!option.is_no_local)
        local_assembler.Assemble(local_contigs);
    
    int num_seed_contigs = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() > option.min_contig)
            ++num_seed_contigs;
    }

    cout << "seed contigs " << num_seed_contigs <<  " local contigs " << local_contigs.size() << endl;
    WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
}

void Iterate(int kmer_size, int new_kmer_size)
{
    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    deque<Sequence> local_contigs;
    ReadSequence(option.local_contig_file(kmer_size), local_contigs);

    deque<Sequence> multi_contigs;
    ReadSequence(option.graph_file(kmer_size), multi_contigs);

    deque<Sequence> transcripts;
    ReadSequence(option.transcript_file(kmer_size), transcripts);

    uint64_t sum = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
        sum += contigs[i].size();
    HashGraph hash_graph(kmer_size);
    hash_graph.reserve(sum);

    deque<Sequence> old_contigs;
    old_contigs.insert(old_contigs.end(), contigs.begin(), contigs.end());
    old_contigs.insert(old_contigs.end(), local_contigs.begin(), local_contigs.end());
    old_contigs.insert(old_contigs.end(), multi_contigs.begin(), multi_contigs.end());
    old_contigs.insert(old_contigs.end(), transcripts.begin(), transcripts.end());
    contigs.clear();
    local_contigs.clear();
    multi_contigs.clear();

    IterateHashGraph(assembly_info, new_kmer_size, option.min_support, hash_graph, old_contigs);
    kmer_size = new_kmer_size;
    old_contigs.clear();

    if (kmer_size < read_length)
        hash_graph.RefreshEdges();
    else
        hash_graph.AddAllEdges();

    Assemble(hash_graph);
}

void FindIsoforms(ContigGraph &contig_graph, deque<Sequence> &transcripts)
{
    transcripts.clear();

    deque<deque<ContigGraphVertexAdaptor> > components;
    deque<string> component_strings;
    contig_graph.GetComponents(components, component_strings);

    int num[100] = {0};
    for (unsigned i = 0; i < components.size(); ++i)
    {
        if ((int)components[i].size() < max_component_size)
            ++num[components[i].size()];
        else
        {
            int total = 0;
            for (unsigned j = 0; j < components[i].size(); ++j)
                total += components[i][j].contig_size();
        }
    }

    int kmer_size = contig_graph.kmer_size();
    int output_component = 0;

    //FastaWriter transcript_writer(option.label_transcript_file(kmer_size));
    FastaWriter component_writer(option.component_file(kmer_size));
    for (unsigned i = 0; i < components.size(); ++i)
    {
        for (unsigned j = 0; j < components[i].size(); ++j)
        {
            Sequence seq = components[i][j].contig();
            component_writer.Write(seq, FormatString("component_%d_%d", i, j));
        }
    }

    //int tran_index = 0;
    for (unsigned i = 0; i < components.size(); ++i)
    {
        int num_begin = 0;
        int num_end = 0;
        for (unsigned j = 0; j < components[i].size(); ++j)
        {
          if (components[i][j].in_edges().size() == 0)
              ++num_begin;
          if (components[i][j].out_edges().size() == 0)
              ++num_end;
        }

        if (num_end < num_begin)
        {
            for (unsigned j = 0; j < components[i].size(); ++j)
                components[i][j].ReverseComplement();
            swap(num_begin, num_end);
        }

        if ((int)components[i].size() <= max_component_size)//  * 100)
        {
            {
                ++output_component;
            }

            deque<ContigGraphPath> all_isoforms;
            for (unsigned j = 0; j < components[i].size(); ++j)
            {
              if (components[i][j].in_edges().size() == 0)
                {
                    ContigGraphPath path;
                    path.Append(components[i][j], 0);
                    components[i][j].status().SetUsedFlag();
                    deque<ContigGraphPath> isoforms;
                    FindIsoforms(path, isoforms, contig_graph, all_isoforms);
                    components[i][j].status().ResetUsedFlag();

                    all_isoforms.insert(all_isoforms.end(), isoforms.begin(), isoforms.end());

                    for (unsigned k = 0; k < isoforms.size(); ++k)
                    {
                        Sequence contig;
                        ContigInfo contig_info;
                        isoforms[k].Assemble(contig, contig_info);

                        transcripts.push_back(contig);

                        //transcript_writer.Write(contig, FormatString("transcript_%d_%d", i, tran_index++));
                    }
                }
            }

            for (unsigned j = 0; j < all_isoforms.size(); ++j)
            {
                for (unsigned k = 0; k < all_isoforms[j].num_nodes(); ++k)
                    all_isoforms[j][k].status().SetUsedFlag();
            }

            for (unsigned j = 0; j < components[i].size(); ++j)
            {
                if (!components[i][j].status().IsUsed())
                {
                    transcripts.push_back(components[i][j].contig());

                    //transcript_writer.Write(components[i][j].contig(), FormatString("transcript_%d_%d", i, tran_index++));
                }
            }
        }
        else
        {
            for (unsigned j = 0; j < components[i].size(); ++j)
            {
                transcripts.push_back(components[i][j].contig());

                //transcript_writer.Write(components[i][j].contig(), FormatString("transcript_%d_%d", i, tran_index++));
            }
        }
    }

    contig_graph.ClearStatus();
}

