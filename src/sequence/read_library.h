/**
 * @file read_library.h
 * @brief ShortReadLibrary and LongReadLibrary Class for holding a set of short or long reads.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.10
 * @date 2012-02-14
 */

#ifndef __SEQUENCE_READ_LIBRARY_H_

#define __SEQUENCE_READ_LIBRARY_H_

#include <stdint.h>

#include <deque>
#include <string>

#include "sequence/sequence.h"
#include "sequence/short_sequence.h"


/**
 * @brief It is class for storing short reads.
 */
class ShortReadLibrary
{
public:
    friend uint64_t ReadLibrary(const std::string &filename, ShortReadLibrary &library);

    std::deque<ShortSequence> &reads() { return reads_; }
    const std::deque<ShortSequence> &reads() const { return reads_; }

    ShortSequence &operator [](int64_t index) 
    { return reads_[index]; }

    const ShortSequence &operator [](int64_t index) const
    { return reads_[index]; }

    uint64_t size() const { return reads_.size(); }
    void clear() { reads_.clear(); }

    bool is_paired() const { return is_paired_; }
    void set_paired(bool is_paired) { is_paired_ = is_paired; }

    double insert_distance() const { return insert_distance_; }
    void set_insert_distance(double insert_distance) { insert_distance_ = insert_distance; }

    double sd() const { return sd_; }
    void set_sd(double sd) { sd_ = sd; }

private:
    std::deque<ShortSequence> reads_;
    bool is_paired_;
    double insert_distance_;
    double sd_;
};

/**
 * @brief It is class for storing long reads.
 */
class LongReadLibrary
{
public:
    friend uint64_t ReadLibrary(const std::string &filename, LongReadLibrary &library);

    std::deque<Sequence> &reads() { return reads_; }
    const std::deque<Sequence> &reads() const { return reads_; }

    Sequence &operator [](int64_t index) 
    { return reads_[index]; }

    const Sequence &operator [](int64_t index) const
    { return reads_[index]; }

    uint64_t size() const { return reads_.size(); }
    void clear() { reads_.clear(); }

private:
    std::deque<Sequence> reads_;
};

#endif

