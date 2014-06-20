/**
 * @file short_sequence.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#include "sequence/short_sequence.h"

#include <iostream>
#include <stdexcept>

#include "basic/kmer.h"
#include "sequence/compact_sequence.h"
#include "sequence/sequence.h"


using namespace std;

const ShortSequence &ShortSequence::Assign(const Sequence &seq, int offset, size_t length)
{
    if (length == std::string::npos || length > seq.size() - offset)
        length = seq.size() - offset;

    if (length > max_size())
        throw logic_error("ShortSequence: Sequence is too long.");

    for (unsigned i = 0; i < length; ++i)
        set_base(i, seq[i + offset]);
    resize(length);
    return *this;
}

const ShortSequence &ShortSequence::Assign(const CompactSequence &compact_seq, int offset, size_t length)
{
    if (length == std::string::npos || length > compact_seq.size() - offset)
        length = compact_seq.size() - offset;

    if (length > max_size())
        throw logic_error("ShortSequence: Sequence is too long.");

    for (unsigned i = 0; i < length; ++i)
        set_base(i, compact_seq[i + offset]);
    resize(length);
    return *this;
}

const ShortSequence &ShortSequence::Assign(const Kmer &kmer)
{
    if (kmer.size() > max_size())
        throw logic_error("ShortSequence: Sequence is too long.");

    for (unsigned i = 0; i < kmer.size(); ++i)
        set_base(i, kmer[i]);
    resize(kmer.size());
    return *this;
}
        

