/**
 * @file local_assembler.h
 * @brief Local Assembler which assembles paired-end reads locally.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2012-01-05
 */

#ifndef __ASSEMBLY_LOCAL_ASSEMBLER_H_

#define __ASSEMBLY_LOCAL_ASSEMBLER_H_

#include "assembly/assembly_utility.h"
#include "sequence/sequence.h"
#include "graph/hash_graph.h"
#include "misc/hash_aligner.h"

#include <omp.h>

#include <deque>
#include <vector>


/**
 * @brief It is a local assembler which assemble the ends of each contigs 
 * and reads with another end aligned to it. The input is the reads, contigs
 * and alignment result between reads and contigs.
 */
class LocalAssembler
{
public:
    LocalAssembler() {}
    ~LocalAssembler();

    void Initialize(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs);

    AssemblyInfo &assembly_info() { return *assembly_info_; }
    std::deque<Sequence> &contigs() { return *contigs_; }
    std::deque<Sequence> &local_contigs() { return local_contigs_; }
    std::deque<Sequence> &extended_contigs() { return extended_contigs_; }
    std::deque<Sequence> &extensions() { return extensions_; }
    std::vector<std::deque<uint64_t> > &in_reads() { return in_reads_; }
    std::vector<std::deque<uint64_t> > &out_reads() { return out_reads_; }

    int num_threads() const { return num_threads_; }
    void set_num_threads(int num_threads) { num_threads_ = num_threads; }

    int mink() const { return mink_; }
    void set_mink(int mink) { mink_ = mink; }

    int maxk() const { return maxk_; }
    void set_maxk(int maxk) { maxk_ = maxk; }
    
    int step() const { return step_; }
    void set_step(int step) { step_ = step; }

    int min_contig() const { return min_contig_; }
    void set_min_contig(int min_contig) { min_contig_ = min_contig; }

    void set_insert_distance(double median, double sd)
    //{ median_ = median; sd_ = sd; local_range_ = median; }
    { median_ = median; sd_ = sd; local_range_ = std::min(median*2, median + 3*sd); }
    double median() const { return median_; }
    double sd() const { return sd_; }
    double local_range() const { return local_range_; }
    int read_length() const { return read_length_; }

    int64_t Assemble(std::deque<Sequence> &contigs);

    void AddReadByHashAlignerRecord(HashAlignerRecord &record, int read_id)
    {
        if (record.ref_from >= record.ref_length - local_range_)
        {
            int offset = ((record.ref_length - record.ref_from) << 4) | (record.query_length - record.match_length);
            if (record.is_reverse == false)
                AddOutRead(record.ref_id, offset, read_id ^ 1);
            else
                AddInRead(record.ref_id, offset, read_id ^ 1);
        }
    }

    void AddInRead(int contig_id, int offset, int read_id)
    { 
        omp_set_lock(&locks_[contig_id]); 
        in_reads_[contig_id].push_back(((uint64_t)offset << 32) | read_id); 
        omp_unset_lock(&locks_[contig_id]); 
    }

    void AddOutRead(int contig_id, int offset, int read_id)
    { 
        omp_set_lock(&locks_[contig_id]); 
        out_reads_[contig_id].push_back(((uint64_t)offset << 32) | read_id); 
        omp_unset_lock(&locks_[contig_id]); 
    }

    static void *LocalAssembleThread(void *p);

    void LocalAssemble(Sequence &contig, std::deque<uint64_t> &local_reads, std::deque<Sequence> &local_contigs);
    void IterativeAssemble(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs);

private:
    struct LocalAssemblyTask
    {
        int id;
        LocalAssembler *local_assembler;
        std::deque<Sequence> local_contigs;
    };

    AssemblyInfo *assembly_info_;
    std::deque<Sequence> *contigs_;
    std::deque<Sequence> local_contigs_;
    std::deque<Sequence> extended_contigs_;
    std::deque<Sequence> extensions_;
    std::vector<omp_lock_t> locks_;
    std::vector<std::deque<uint64_t> > in_reads_;
    std::vector<std::deque<uint64_t> > out_reads_;
    int read_length_;
    int num_threads_;
    int mink_;
    int maxk_;
    int step_;
    int min_contig_;
    double median_;
    double sd_;
    double local_range_;
};

/**
 * @brief It is a assembler which get contigs from reference. If a region is
 * supported by many reads, it will be extracted as contigs. 
 */
class LocalAssemblerWithReference
{
public:
    LocalAssemblerWithReference() {}
    ~LocalAssemblerWithReference();

    void Initialize(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs);

    AssemblyInfo &assembly_info() { return *assembly_info_; }
    std::deque<Sequence> &contigs() { return *contigs_; }
    std::deque<Sequence> &local_contigs() { return local_contigs_; }
    std::deque<std::deque<std::deque<int> > > &read_tables() { return read_tables_; }
    std::deque<std::deque<int> > &matches() { return matches_; }
    std::deque<std::deque<long long> > &read_offset() { return read_offset_; }
    std::deque<std::deque<double> > &coverage_table() { return coverage_table_; }

    int num_threads() const { return num_threads_; }
    void set_num_threads(int num_threads) { num_threads_ = num_threads; }

    int mink() const { return mink_; }
    void set_mink(int mink) { mink_ = mink; }

    int maxk() const { return maxk_; }
    void set_maxk(int maxk) { maxk_ = maxk; }
    
    int step() const { return step_; }
    void set_step(int step) { step_ = step; }

    int min_contig() const { return min_contig_; }
    void set_min_contig(int min_contig) { min_contig_ = min_contig; }

    void set_insert_distance(double median, double sd);
    double median() const { return median_; }
    double sd() const { return sd_; }
    double local_range() const { return local_range_; }
    int read_length() const { return read_length_; }
    double expeceted_coverage() const { return expeceted_coverage_; }
    void set_max_gap(int max_gap) { max_gap_ = max_gap; }
    int max_gap() const { return max_gap_; }

    void GetSegments(std::deque<Sequence> &contigs);
    int64_t Assemble(std::deque<Sequence> &contigs);

    void AddReadByHashAlignerRecord(HashAlignerRecord &record, int read_id)
    {
        if (record.is_reverse)
            record.ReverseComplement();

        int x = record.ref_id;
//        int y = record.ref_from / local_range();
//
//        read_offset_[read_id].push_back((int64_t(record.ref_id) << 32) | record.ref_from);
        for (int z = record.ref_from; z < record.ref_to; ++z)
            matches_[x][z]++;

//        if (read_tables_[x][y].empty() || read_tables_[x][y].back() != (int)(read_id & ~1))
//            read_tables_[x][y].push_back(read_id & ~1);
//        if (y > 0 && (read_tables_[x][y-1].empty() || read_tables_[x][y-1].back() != (int)(read_id & ~1)))
//            read_tables_[x][y-1].push_back(read_id & ~1);
    }

    static void *LocalAssembleThread(void *p);

    void IterativeAssemble(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs);

private:
    struct LocalAssemblyTask
    {
        int id;
        LocalAssemblerWithReference *local_assembler;
        std::deque<Sequence> local_contigs;
    };
    
    AssemblyInfo *assembly_info_;
    std::deque<Sequence> *contigs_;
    std::deque<Sequence> local_contigs_;
    std::vector<omp_lock_t> locks_;
    std::deque<std::deque<int> > matches_;
    std::deque<std::deque<std::deque<int> > > read_tables_;
    std::deque<std::deque<long long> > read_offset_;
    std::deque<std::deque<double> > coverage_table_;
    int read_length_;
    int num_threads_;
    int mink_;
    int maxk_;
    int step_;
    int min_contig_;
    int max_gap_;
    double median_;
    double sd_;
    double local_range_;
    double expeceted_coverage_;
};

#endif

