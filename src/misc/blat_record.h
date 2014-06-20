/**
 * @file blat_record.h
 * @brief BlatRecord Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.3
 * @date 2011-09-06
 */

#ifndef __MISC_BLAT_RECORD_H_

#define __MISC_BLAT_RECORD_H_

#include <stdint.h>

#include <string>
#include <vector>


struct BlatBlock
{
    int64_t query_from;
    int64_t ref_from;
    int64_t size;
};

/**
 * @brief It is an alignment record in BLAT output.
 */
struct BlatRecord
{
    std::string query_name;
    std::string ref_name;
    int64_t match_count;
    int64_t mismatch_count;
    int64_t query_from;
    int64_t query_to;
    int64_t query_length;
    int64_t ref_from;
    int64_t ref_to;
    int64_t ref_length;
    bool is_reverse;

    bool operator <(const BlatRecord &r) const
    { return match_count > r.match_count; }

    void ReverseComplement()
    {
        int to = ref_length - ref_from;
        int from = ref_length - ref_to;
        ref_from = from;
        ref_to = to;
        is_reverse = !is_reverse;
    }

    std::vector<BlatBlock> blocks;

    void Parse(const std::string &record);
};

#endif

