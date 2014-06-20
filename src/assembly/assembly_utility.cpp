/**
 * @file assembly_utility.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-06
 */

#include "assembly/assembly_utility.h"

#include <fstream>
#include <iostream>

#include "graph/contig_graph.h"
#include "graph/hash_graph.h"
#include "misc/hash_aligner.h"
#include "misc/hash_aligner.h"
#include "misc/utils.h"
#include "sequence/sequence.h"
#include "sequence/sequence_io.h"
#include "sequence/short_sequence.h"


using namespace std;

struct Score
{
    Score() { counts[0] = counts[1] = counts[2] = counts[3] = 0; best = -1; }
    int &operator [](int index) { return counts[index]; }
    int counts[4];
    int best;
};

int AssemblyInfo::read_length() const
{
    Histgram<int> length_histgram;
    for (unsigned i = 0; i < reads.size(); ++i)
        length_histgram.insert(reads[i].size());
    return length_histgram.median();
}

bool CompareHashAlignerRecord(const HashAlignerRecord &x, const HashAlignerRecord &y)
{
    return x.match_length > y.match_length;
}

void ReadInput(const std::string &read_file, const std::string &long_read_file, AssemblyInfo &assembly_info)
{
    if (read_file != "")
        ReadSequence(read_file, assembly_info.reads);
    if (long_read_file != "")
        ReadSequence(long_read_file, assembly_info.long_reads);

    assembly_info.read_flags.resize(assembly_info.reads.size());
    assembly_info.long_read_flags.resize(assembly_info.long_reads.size());

    assembly_info.ClearStatus();
}

bool ReadHashAlignerRecords(FILE *fp, deque<HashAlignerRecord> &records)
{
    int size = 0;
    if (fread(&size, sizeof(int), 1, fp) != 1)
        return false;

    records.resize(size);
    for (unsigned j = 0; j < records.size(); ++j)
    {
        if (fread(&records[j], sizeof(HashAlignerRecord), 1, fp) != 1)
        {
            records.resize(j);
            return false;
        }
    }

    return true;
}

bool WriteHashAlignerRecords(FILE *fp, deque<HashAlignerRecord> &records)
{
    int size = records.size();
    if (fwrite(&size, sizeof(int), 1, fp) != 1)
        return false;

    for (unsigned j = 0; j < records.size(); ++j)
    {
        if (fwrite(&records[j], sizeof(HashAlignerRecord), 1, fp) != 1)
            return false;
    }

    return true;
}

void ReadHashAlignerRecordBlock(FILE *fp, vector<HashAlignerRecord> &records)
{
    deque<HashAlignerRecord> tmp_records;
    for (unsigned i = 0; i < records.size(); ++i)
    {
        records[i].match_length = 0;
        ReadHashAlignerRecords(fp, tmp_records);
        if (tmp_records.size() == 1)
            records[i] = tmp_records[0];
    }
}

void WriteHashAlignerRecordBlock(FILE *fp, vector<HashAlignerRecord> &records)
{
    for (unsigned i = 0; i < records.size(); ++i)
    {
        HashAlignerRecord &record = records[i];
        int num = 0;
        if (record.match_length != 0)
        {
            num = 1;
            fwrite(&num, sizeof(int), 1, fp);
            fwrite(&record, sizeof(HashAlignerRecord), 1, fp);
        }
        else
            fwrite(&num, sizeof(int), 1, fp);
    }
}

uint64_t WriteContig(const std::string &filename, const std::deque<Sequence> &contigs, 
        const std::deque<ContigInfo> &contig_infos, const std::string &prefix, int min_contig)
{
    FastaWriter writer(filename);
    return WriteContig(writer, contigs, contig_infos, prefix, min_contig);
}

uint64_t WriteContig(FastaWriter &writer, const std::deque<Sequence> &contigs, 
        const std::deque<ContigInfo> &contig_infos, const std::string &prefix, int min_contig)
{
    int64_t num_contigs = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int64_t)contigs[i].size() >= min_contig)
        {
            ++num_contigs;
            writer.Write(contigs[i], 
                    FormatString("%s_%d length_%d read_count_%d", 
                        prefix.c_str(), i, 
                        (int)contigs[i].size(), 
                        (int)contig_infos[i].kmer_count()));
        }
    }
    return num_contigs;
}

void EstimateDistance(const string &align_file, double &mean, double &sd)
{
    Histgram<int> histgram;
    FILE *falign = OpenFile(align_file, "rb");

    deque<HashAlignerRecord> records1;
    deque<HashAlignerRecord> records2;
    
    while (ReadHashAlignerRecords(falign, records1) 
            && ReadHashAlignerRecords(falign, records2))
    {
        if (records1.size() == 1 && records2.size() == 1)
        {
            HashAlignerRecord record1 = records1[0];
            HashAlignerRecord record2 = records2[0];
            record2.ReverseComplement();

            if (record1.ref_id == record2.ref_id && record1.is_reverse == record2.is_reverse)
            {
                histgram.insert((record2.ref_to + record2.query_length - record2.query_to)
                        - (record1.ref_from - record1.query_from));
            }
        }
    }

    histgram.Trim(0.01);
    cout << "distance mean " << histgram.mean() << " sd " << histgram.sd() << endl;
    mean = histgram.mean();
    sd = histgram.sd();

    fclose(falign);
}

void BuildKmerFile(AssemblyInfo &assembly_info, int kmer_size, 
        int min_count, int prefix_length, const std::string &kmer_file)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    vector<bool> &read_flags = assembly_info.read_flags;
    deque<Sequence> &long_reads = assembly_info.long_reads;
    vector<bool> &long_read_flags = assembly_info.long_read_flags;

    HashGraph hash_graph(kmer_size);
    ofstream fkmer(kmer_file.c_str(), ios_base::out | ios_base::binary);

    uint64_t mask = (1ULL << prefix_length) - 1;
    for (int prefix = 0; prefix < (1 << prefix_length); ++prefix)
    {
#pragma omp parallel for schedule(static, 1)
        for (int64_t i = 0; i < (int64_t)reads.size(); ++i)
        {
            if (!read_flags[i])
                continue;
            Sequence seq(reads[i]);
            hash_graph.InsertKmersWithPrefix(seq, prefix, mask);
        }

#pragma omp parallel for schedule(static, 1)
        for (int64_t i = 0; i < (int64_t)long_reads.size(); ++i)
        {
            if (!long_read_flags[i])
                continue;
            hash_graph.InsertKmersWithPrefix(long_reads[i], prefix, mask);
        }

        hash_graph.RefreshVertices(min_count);
        fkmer << hash_graph;
        hash_graph.clear();
    }
}

void ReadKmerFile(const std::string &kmer_file, HashGraph &hash_graph)
{
    ifstream fkmer(kmer_file.c_str(), ios_base::binary | ios_base::in);
    fkmer.seekg(0, ios_base::end);
    int64_t num_nodes = fkmer.tellg() / sizeof(HashGraphVertex);
    hash_graph.reserve(num_nodes);
    fkmer.seekg(0, ios_base::beg);
    fkmer >> hash_graph;
}

void WriteKmerFile(const std::string &kmer_file, HashGraph &hash_graph)
{
    ofstream fkmer(kmer_file.c_str(), ios_base::binary | ios_base::out);
    fkmer << hash_graph;
}

void InsertInternalKmers(AssemblyInfo &assembly_info, HashGraph &hash_graph, int min_count)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    deque<Sequence> &long_reads = assembly_info.long_reads;

#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)reads.size(); ++i)
    {
        Sequence seq(reads[i]);
        hash_graph.InsertInternalKmers(seq, min_count);
    }

#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)long_reads.size(); ++i)
        hash_graph.InsertInternalKmers(long_reads[i], min_count);

    hash_graph.RestoreAndMergeEdges();
    hash_graph.RefreshEdges();
}

void IterateHashGraph(AssemblyInfo &assembly_info, int new_kmer_size, int min_support, HashGraph &hash_graph, deque<Sequence> &old_contigs)
{
    int old_kmer_size = hash_graph.kmer_size();
    deque<ShortSequence> &reads = assembly_info.reads;
    deque<Sequence> &long_reads = assembly_info.long_reads;
    vector<bool> &read_flags = assembly_info.read_flags;
    vector<bool> &long_read_flags = assembly_info.long_read_flags;

#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)old_contigs.size(); ++i)
        hash_graph.InsertUncountKmers(old_contigs[i]);
    hash_graph.AddAllEdges();

    deque<Sequence> contigs;
    hash_graph.Assemble(contigs);
    hash_graph.clear();

    uint64_t sum = 0;
    int d = new_kmer_size - old_kmer_size;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() - old_kmer_size + 1 >= 2*d + 2)
            sum += 2*d + 2;
        else if ((int)contigs[i].size() >= old_kmer_size)
            sum += contigs[i].size() - old_kmer_size + 1;
    }

    HashGraph old_hash_graph(old_kmer_size);
    old_hash_graph.reserve(sum);
#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        Sequence seq;
        seq.Assign(contigs[i], 0, min(new_kmer_size, (int)contigs[i].size()));
        old_hash_graph.InsertKmers(seq);

        seq.Assign(contigs[i], max(0, (int)contigs[i].size() - new_kmer_size), min(new_kmer_size, (int)contigs[i].size()));
        old_hash_graph.InsertKmers(seq);
    }
    //cout << "old kmer " << old_hash_graph.num_vertices() << endl;

    hash_graph.set_kmer_size(new_kmer_size);
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)reads.size(); ++i)
    {
        if (!read_flags[i])
            continue;

        Sequence seq(reads[i]);
        InsertIterativeKmers(old_hash_graph, seq, hash_graph);
    }

#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)long_reads.size(); ++i)
    {
        if (!long_read_flags[i])
            continue;

        InsertIterativeKmers(old_hash_graph, long_reads[i], hash_graph);
    }
    
#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)assembly_info.ref_contigs.size(); ++i)
        InsertIterativeKmers(old_hash_graph, assembly_info.ref_contigs[i], hash_graph);

    old_hash_graph.clear();
    {
        HashGraph tmp_hash_graph;
        tmp_hash_graph.swap(old_hash_graph);
    }

    hash_graph.RefreshVertices(min_support);

#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)old_contigs.size(); ++i)
        hash_graph.InsertUncountKmers(old_contigs[i]);
    hash_graph.ClearCount();

    InsertExistKmers(assembly_info, hash_graph);
}

void InsertExistKmers(AssemblyInfo &assembly_info, HashGraph &hash_graph)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    deque<Sequence> &long_reads = assembly_info.long_reads;
    vector<bool> &read_flags = assembly_info.read_flags;
    vector<bool> &long_read_flags = assembly_info.long_read_flags;

//#pragma omp parallel for
//    for (int64_t i = 0; i < (int64_t)assembly_info.ref_contigs.size(); ++i)
//        hash_graph.InsertExistKmers(assembly_info.ref_contigs[i]);
    //hash_graph.ClearCount();

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)reads.size(); ++i)
    {
        if (!read_flags[i])
            continue;

        Sequence seq(reads[i]);
        //read_flags[i] = (hash_graph.InsertExistKmers(seq) != 0);
        hash_graph.InsertExistKmers(seq);
    }

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)long_reads.size(); ++i)
    {
        if (!long_read_flags[i])
            continue;

        //long_read_flags[i] = (hash_graph.InsertExistKmers(long_reads[i]) != 0);
        hash_graph.InsertExistKmers(long_reads[i]);
    }
}

int64_t InsertIterativeKmers(const HashGraph &old_hash_graph, const Sequence &seq, HashGraph &hash_graph, int kmer_count)
{
    int old_kmer_size = old_hash_graph.kmer_size();
    int new_kmer_size = hash_graph.kmer_size();

    Kmer old_kmer(old_kmer_size);
    Kmer new_kmer(new_kmer_size);
    int length = 0;
    int count = 0;
    int num_iterative_kmers = 0;
    for (uint32_t j = 0; j < seq.size(); ++j)
    {
        old_kmer.ShiftAppend(seq[j]);
        new_kmer.ShiftAppend(seq[j]);

        length = (seq[j] < 4) ? length + 1 : 0;

        count = (length >= old_kmer_size && old_hash_graph.FindVertex(old_kmer) != NULL) ? count+1 : 0;
        if (count >= new_kmer_size - old_kmer_size + 1)
        {
            ++num_iterative_kmers;
            HashGraphVertex *vertex = hash_graph.InsertVertex(new_kmer, kmer_count);
            HashGraphVertexAdaptor adaptor(vertex, new_kmer != vertex->kmer());
            if (length > new_kmer_size && seq[j-new_kmer_size] < 4)
                adaptor.in_edges().Add(3 - seq[j-new_kmer_size]);
            if (j+1 < seq.size() && seq[j+1] < 4)
                adaptor.out_edges().Add(seq[j+1]);
        }
    }

    return num_iterative_kmers;
}

int64_t AlignReads(AssemblyInfo &assembly_info, HashAligner &hash_aligner, double similar, const std::string &align_file, bool is_all)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    vector<bool> &read_flags = assembly_info.read_flags;

    FILE *falign = OpenFile(align_file, "wb");

    int64_t num_aligned_reads = 0;
    int buffer_size = (1 << 20) * omp_get_max_threads();
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<deque<HashAlignerRecord> > buffer_records(omp_get_max_threads());
        vector<HashAlignerRecord> all_records(size);

#pragma omp parallel for schedule(static, 1) reduction(+: num_aligned_reads)
        for (int64_t i = 0; i < size; ++i)
        {
            all_records[i].match_length = 0;
            if (read_flags[offset + i] || is_all)
            {
                deque<HashAlignerRecord> &records = buffer_records[omp_get_thread_num()];
                Sequence seq(reads[offset + i]);
                hash_aligner.AlignRead(seq, records, seq.size() * similar, 1);

                if (records.size() == 1)
                {
                    ++num_aligned_reads;
                    all_records[i] = records[0];
                }
            }
        }

        WriteHashAlignerRecordBlock(falign, all_records);
    }

    fclose(falign);

    return num_aligned_reads;
}

int64_t AlignReadsLocal(AssemblyInfo &assembly_info, HashAligner &hash_aligner, int min_match, int max_mismatch, const std::string &align_file)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    vector<bool> &read_flags = assembly_info.read_flags;

    FILE *falign = OpenFile(align_file, "wb");

    int64_t num_aligned_reads = 0;
    int buffer_size = (1 << 16) * omp_get_max_threads();
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<deque<HashAlignerRecord> > buffer_records(omp_get_max_threads());
        vector<deque<HashAlignerRecord> > all_records(size);

#pragma omp parallel for schedule(static, 1) reduction(+: num_aligned_reads)
        for (int64_t i = 0; i < size; ++i)
        {
            all_records[i].clear();
            if (read_flags[offset + i])
            {
                deque<HashAlignerRecord> &records = buffer_records[omp_get_thread_num()];
                Sequence seq(reads[offset + i]);
                hash_aligner.AlignReadLocal(seq, records, min_match, max_mismatch, 1000);
                all_records[i] = records;
            }
        }

        for (int64_t i = 0; i < size; ++i)
        {
            if (all_records[i].size() > 0)
                ++num_aligned_reads;
            WriteHashAlignerRecords(falign, all_records[i]);
        }
    }

    fclose(falign);

    return num_aligned_reads;
}

int64_t AlignReadsMultiple(AssemblyInfo &assembly_info, HashAligner &hash_aligner, double similar, const std::string &align_file, bool is_all)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    vector<bool> &read_flags = assembly_info.read_flags;

    FILE *falign = OpenFile(align_file, "wb");

    int64_t num_aligned_reads = 0;
    int buffer_size = (1 << 20) * omp_get_max_threads();
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<deque<HashAlignerRecord> > buffer_records(omp_get_max_threads());
        vector<deque<HashAlignerRecord> > all_records(size);

#pragma omp parallel for schedule(static, 1) reduction(+: num_aligned_reads)
        for (int64_t i = 0; i < size; ++i)
        {
            all_records[i].resize(0);
            if (read_flags[offset + i] || is_all)
            {
                deque<HashAlignerRecord> &records = buffer_records[omp_get_thread_num()];
                Sequence seq(reads[offset + i]);
                hash_aligner.AlignReadLocal(seq, records, seq.size() * similar, 10);
                all_records[i] = records;
            }
        }

        for (int64_t i = 0; i < size; ++i)
        {
            if (all_records[i].size() > 0)
                ++num_aligned_reads;

            WriteHashAlignerRecords(falign, all_records[i]);
        }
    }

    fclose(falign);

    return num_aligned_reads;
}

void CorrectReads(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs, 
        std::deque<ContigInfo> &contig_infos, const std::string &align_file, int max_mismatch, int min_match)
{
    int read_length = assembly_info.read_length();

    vector<vector<uint32_t> > read_counts(contigs.size());
    deque<vector<Score> > scores_table(contigs.size());
    deque<vector<short> > overlap(contigs.size());

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        scores_table[i].resize(contigs[i].size());
        overlap[i].resize(contigs[i].size(), 0);
        read_counts[i].resize(omp_get_max_threads(), 0);
    }

    FILE *falign = OpenFile(align_file, "rb");
    deque<ShortSequence> &reads = assembly_info.reads;

    int64_t buffer_size = (1 << 20) * omp_get_max_threads();
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);
        ReadHashAlignerRecordBlock(falign, all_records);

#pragma omp parallel for
        for (int64_t i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];
            if (record.match_length != 0)
            //if (record.match_length != 0 && record.query_length - record.match_length <= max_mismatch)
            {
//#pragma omp atomic
                ++read_counts[record.ref_id][omp_get_thread_num()];

                Sequence seq(reads[offset + i]);
                bool is_reverse = record.is_reverse;
                if (is_reverse)
                {
                    record.ReverseComplement();
                    seq.ReverseComplement();
                }

                vector<Score> &scores = scores_table[record.ref_id];
                int id = record.ref_id;
                for (unsigned j = 0; j < seq.size(); ++j)
                {
                    int index = record.ref_from + j;
#pragma omp atomic
                    scores[index][seq[j]]++;
                    if (min(int(j), int(seq.size() - j)) >= read_length/3)
                    {
#pragma omp atomic
                        overlap[id][index]++;
                    }
                }

                if (is_reverse)
                {
                    record.ReverseComplement();
                    seq.ReverseComplement();
                }
            }
        }
    }
    fclose(falign);

    int64_t conformed_bases = 0;
#pragma omp parallel for reduction(+: conformed_bases)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        for (unsigned j = 0; j < contigs[i].size(); ++j)
        {
            int aux[4];
            for (int k = 0; k < 4; ++k)
                aux[k] = scores_table[i][j][k];
            sort(aux, aux + 4);
            if (aux[3] > aux[2] * 4 && aux[3] >= min_match && overlap[i][j] >= 1)
            {
                ++conformed_bases;
                int best = -1;
                for (int k = 0; k < 4; ++k)
                {
                    if (scores_table[i][j][k] == aux[3])
                        best = k;
                }
                contigs[i][j] = best;
                scores_table[i][j].best = best;
            }
        }
    }

    falign = OpenFile(align_file, "rb");
    int64_t corrected_reads = 0;
    int64_t corrected_bases = 0;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);
        vector<int> flags(size, 1);
        ReadHashAlignerRecordBlock(falign, all_records);

#pragma omp parallel for reduction(+: corrected_reads, corrected_bases)
        for (int64_t i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];
            if (record.match_length != 0)
            //if (record.match_length != 0 && record.query_length - record.match_length <= max_mismatch)
            {
                Sequence seq(reads[offset + i]);

                bool is_reverse = record.is_reverse;
                if (is_reverse)
                {
                    record.ReverseComplement();
                    seq.ReverseComplement();
                }

                int match_count = 0;
                bool is_all_confirmed = true;
                vector<Score> &scores = scores_table[record.ref_id];
                Sequence &contig = contigs[record.ref_id];
                for (unsigned j = 0; j < seq.size(); ++j)
                {
                    if (seq[j] == contig[record.ref_from + j])
                        ++match_count;
                    if (scores[record.ref_from + j].best == -1)
                        is_all_confirmed = false;
                }

                if (match_count >= (int)seq.size() - max_mismatch && is_all_confirmed)
                {
                    ++corrected_reads;

                    flags[i] = false;
                    for (unsigned j = 0; j < seq.size(); ++j)
                    {
                        if (scores[record.ref_from + j].best != -1
                                && seq[j] != contig[record.ref_from + j])
                        {
                            seq[j] = scores[record.ref_from + j].best;
                            ++corrected_bases;
                        }
                    }
                }

                if (is_reverse)
                {
                    record.ReverseComplement();
                    seq.ReverseComplement();
                }

                reads[offset + i] = seq;
            }
        }

        copy(flags.begin(), flags.end(), assembly_info.read_flags.begin() + offset);
    }
    fclose(falign);

    contig_infos.resize(contigs.size());
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        uint32_t sum = 0;
        for (unsigned j = 0; j < read_counts[i].size(); ++j)
            sum += read_counts[i][j];
        contig_infos[i].set_kmer_count(sum);
    }

    cout << "confirmed bases: " << conformed_bases << " correct reads: " << corrected_reads << " bases: " << corrected_bases << endl;
    //cout << "t1 t2 t3: " << t1 << " " << t2 << " " << t3 << endl;
}

void CorrectPairedReads(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs, 
        std::deque<ContigInfo> &contig_infos, const std::string &align_file, double mean, double sd, int max_mismatch, int min_match)
{
    int read_length = assembly_info.read_length();

    vector<vector<uint32_t> > read_counts(contigs.size());
    deque<vector<Score> > scores_table(contigs.size());
    deque<vector<short> > overlap(contigs.size());

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        scores_table[i].resize(contigs[i].size());
        overlap[i].resize(contigs[i].size(), 0);
        read_counts[i].resize(omp_get_max_threads(), 0);
    }

    FILE *falign = OpenFile(align_file, "rb");
    deque<ShortSequence> &reads = assembly_info.reads;

    int64_t buffer_size = (1 << 20) * omp_get_max_threads();
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);
        ReadHashAlignerRecordBlock(falign, all_records);

#pragma omp parallel for
        for (int64_t i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];
            if (record.match_length != 0)
            //if (record.match_length != 0 && record.query_length - record.match_length <= max_mismatch)
            {
//#pragma omp atomic
                ++read_counts[record.ref_id][omp_get_thread_num()];

                Sequence seq(reads[offset + i]);
                bool is_reverse = record.is_reverse;
                if (is_reverse)
                {
                    record.ReverseComplement();
                    seq.ReverseComplement();
                }

                vector<Score> &scores = scores_table[record.ref_id];
                int id = record.ref_id;
                for (unsigned j = 0; j < seq.size(); ++j)
                {
                    int index = record.ref_from + j;
#pragma omp atomic
                    scores[index][seq[j]]++;
                    if (min(int(j), int(seq.size() - j)) >= read_length/3)
                    {
#pragma omp atomic
                        overlap[id][index]++;
                    }
                }

                if (is_reverse)
                {
                    record.ReverseComplement();
                    seq.ReverseComplement();
                }
            }
        }
    }
    fclose(falign);

    int64_t conformed_bases = 0;
//    int64_t t1 = 0;
//    int64_t t2 = 0;
//    int64_t t3 = 0;
#pragma omp parallel for reduction(+: conformed_bases)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        for (unsigned j = 0; j < contigs[i].size(); ++j)
        {
            int aux[4];
            for (int k = 0; k < 4; ++k)
                aux[k] = scores_table[i][j][k];
            sort(aux, aux + 4);
            if (aux[3] > aux[2] * 4 && aux[3] >= min_match && overlap[i][j] >= 1)
            {
                ++conformed_bases;
                int best = -1;
                for (int k = 0; k < 4; ++k)
                {
                    if (scores_table[i][j][k] == aux[3])
                        best = k;
                }
                contigs[i][j] = best;
                scores_table[i][j].best = best;
            }
        }
    }

    falign = OpenFile(align_file, "rb");
    int64_t corrected_reads = 0;
    int64_t corrected_bases = 0;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);
        vector<int> flags(size, 1);
        ReadHashAlignerRecordBlock(falign, all_records);

#pragma omp parallel for reduction(+: corrected_reads, corrected_bases)
        for (int64_t i = 0; i < size; i += 2)
        {
            HashAlignerRecord r1 = all_records[i];
            HashAlignerRecord r2 = all_records[i+1];
            HashAlignerRecord tmp = r2;

            tmp.ReverseComplement();
            double d = fabs(tmp.ref_to - r1.ref_from - mean);

            if (r1.match_length != 0 && r2.match_length != 0
                    && r1.ref_id == r2.ref_id
                    && r1.is_reverse != r2.is_reverse
                    && d < 2*sd)
            {
                Sequence seq1(reads[offset + i]);
                Sequence seq2(reads[offset + i+1]);

                seq2.ReverseComplement();
                r2.ReverseComplement();

                bool is_reverse = r1.is_reverse;
                if (is_reverse)
                {
                    seq1.ReverseComplement();
                    r1.ReverseComplement();
                    seq2.ReverseComplement();
                    r2.ReverseComplement();
                }

                vector<Score> &scores = scores_table[r1.ref_id];
                Sequence &contig = contigs[r1.ref_id];

                bool is_all_confirmed = true;
                for (unsigned j = 0; j < seq1.size(); ++j)
                {
                    if (scores[r1.ref_from + j].best == -1)
                        is_all_confirmed = false;
                }
                for (unsigned j = 0; j < seq2.size(); ++j)
                {
                    if (scores[r2.ref_from + j].best == -1)
                        is_all_confirmed = false;
                }

                if (!is_all_confirmed)
                {
                    goto SINGLE_CORRECT;
                }

                corrected_reads += 2;
                for (unsigned j = 0; j < seq1.size(); ++j)
                {
                    if (seq1[j] != scores[r1.ref_from+j].best)
                    {
                        seq1[j] = contig[r1.ref_from+j];
                        ++corrected_bases;
                    }
                }

                for (unsigned j = 0; j < seq2.size(); ++j)
                {
                    if (seq2[j] != scores[r2.ref_from+j].best)
                    {
                        seq2[j] = contig[r2.ref_from+j];
                        ++corrected_bases;
                    }
                }

                if (is_reverse)
                {
                    seq1.ReverseComplement();
                    r1.ReverseComplement();
                    seq2.ReverseComplement();
                    r2.ReverseComplement();
                }

                seq2.ReverseComplement();
                r2.ReverseComplement();

                reads[offset + i] = seq1;
                reads[offset + i+1] = seq2;
            }
            else 
            {
SINGLE_CORRECT:
                r1 = all_records[i];
                r2 = all_records[i+1];

                if (r1.match_length != 0)
                {
                    int id = i;
                    Sequence seq(reads[offset + i]);
                    HashAlignerRecord record = r1;

                    bool is_reverse = record.is_reverse;
                    if (is_reverse)
                    {
                        record.ReverseComplement();
                        seq.ReverseComplement();
                    }

                    int match_count = 0;
                    bool is_all_confirmed = true;
                    vector<Score> &scores = scores_table[record.ref_id];
                    Sequence &contig = contigs[record.ref_id];
                    for (unsigned j = 0; j < seq.size(); ++j)
                    {
                        if (seq[j] == contig[record.ref_from + j])
                            ++match_count;
                        if (scores[record.ref_from + j].best == -1)
                            is_all_confirmed = false;
                    }

                    if (match_count >= (int)seq.size() - max_mismatch && is_all_confirmed)
                    {
                        ++corrected_reads;

                        for (unsigned j = 0; j < seq.size(); ++j)
                        {
                            if (scores[record.ref_from + j].best != -1
                                    && seq[j] != contig[record.ref_from + j])
                            {
                                seq[j] = scores[record.ref_from + j].best;
                                ++corrected_bases;
                            }
                        }
                    }

                    if (is_reverse)
                    {
                        record.ReverseComplement();
                        seq.ReverseComplement();
                    }

                    reads[offset + id] = seq;
                }
                if (r2.match_length != 0)
                {
                    int id = i+1;
                    Sequence seq(reads[offset + i+1]);
                    HashAlignerRecord record = r2;

                    bool is_reverse = record.is_reverse;
                    if (is_reverse)
                    {
                        record.ReverseComplement();
                        seq.ReverseComplement();
                    }

                    int match_count = 0;
                    bool is_all_confirmed = true;
                    vector<Score> &scores = scores_table[record.ref_id];
                    Sequence &contig = contigs[record.ref_id];
                    for (unsigned j = 0; j < seq.size(); ++j)
                    {
                        if (seq[j] == contig[record.ref_from + j])
                            ++match_count;
                        if (scores[record.ref_from + j].best == -1)
                            is_all_confirmed = false;
                    }

                    if (match_count >= (int)seq.size() - max_mismatch && is_all_confirmed)
                    {
                        ++corrected_reads;

                        for (unsigned j = 0; j < seq.size(); ++j)
                        {
                            if (scores[record.ref_from + j].best != -1
                                    && seq[j] != contig[record.ref_from + j])
                            {
                                seq[j] = scores[record.ref_from + j].best;
                                ++corrected_bases;
                            }
                        }
                    }

                    if (is_reverse)
                    {
                        record.ReverseComplement();
                        seq.ReverseComplement();
                    }

                    reads[offset + id] = seq;
                }
            }
        }

        copy(flags.begin(), flags.end(), assembly_info.read_flags.begin() + offset);
    }
    fclose(falign);

    contig_infos.resize(contigs.size());
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        uint32_t sum = 0;
        for (unsigned j = 0; j < read_counts[i].size(); ++j)
            sum += read_counts[i][j];
        contig_infos[i].set_kmer_count(sum);
    }

    cout << "confirmed bases: " << conformed_bases << " correct reads: " << corrected_reads << " bases: " << corrected_bases << endl;
    //cout << "t1 t2 t3: " << t1 << " " << t2 << " " << t3 << endl;
}

void CorrectContigs(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs, 
        std::deque<ContigInfo> &contig_infos, const std::string &align_file, int max_mismatch, int min_match)
{
    int read_length = assembly_info.read_length();

    vector<vector<uint32_t> > read_counts(contigs.size());
    deque<vector<Score> > scores_table(contigs.size());
    deque<vector<short> > overlap(contigs.size());

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        scores_table[i].resize(contigs[i].size());
        overlap[i].resize(contigs[i].size(), 0);
        read_counts[i].resize(omp_get_max_threads(), 0);
    }

    FILE *falign = OpenFile(align_file, "rb");
    deque<ShortSequence> &reads = assembly_info.reads;

    int64_t buffer_size = (1 << 20) * omp_get_max_threads();
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);
        ReadHashAlignerRecordBlock(falign, all_records);

#pragma omp parallel for
        for (int64_t i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];
            if (record.match_length != 0)
            //if (record.match_length != 0 && record.query_length - record.match_length <= max_mismatch)
            {
//#pragma omp atomic
                ++read_counts[record.ref_id][omp_get_thread_num()];

                Sequence seq(reads[offset + i]);
                bool is_reverse = record.is_reverse;
                if (is_reverse)
                {
                    record.ReverseComplement();
                    seq.ReverseComplement();
                }

                vector<Score> &scores = scores_table[record.ref_id];
                int id = record.ref_id;
                for (unsigned j = 0; j < seq.size(); ++j)
                {
                    int index = record.ref_from + j;
#pragma omp atomic
                    scores[index][seq[j]]++;
                    if (min(int(j), int(seq.size() - j)) >= read_length/3)
                    {
#pragma omp atomic
                        overlap[id][index]++;
                    }
                }

                if (is_reverse)
                {
                    record.ReverseComplement();
                    seq.ReverseComplement();
                }
            }
        }
    }
    fclose(falign);

    int64_t conformed_bases = 0;
#pragma omp parallel for reduction(+: conformed_bases)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        for (unsigned j = 0; j < contigs[i].size(); ++j)
        {
            int aux[4];
            for (int k = 0; k < 4; ++k)
                aux[k] = scores_table[i][j][k];
            sort(aux, aux + 4);
            if (aux[3] > aux[2] * 4 && aux[3] >= min_match && overlap[i][j] >= 1)
            {
                ++conformed_bases;
                int best = -1;
                for (int k = 0; k < 4; ++k)
                {
                    if (scores_table[i][j][k] == aux[3])
                        best = k;
                }
                contigs[i][j] = best;
                scores_table[i][j].best = best;
            }
        }
    }

    cout << "confirmed bases: " << conformed_bases << endl;
}

