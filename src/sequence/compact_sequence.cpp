/**
 * @file compact_sequence.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#include "sequence/compact_sequence.h"

#include <string>

#include "basic/bit_operation.h"
#include "basic/kmer.h"
#include "sequence/sequence.h"
#include "sequence/short_sequence.h"


using namespace std;

std::istream &operator >>(std::istream &stream, CompactSequence &compact_seq)
{
    Sequence seq;
    stream >> seq;
    compact_seq = seq;
    return stream;
}

std::ostream &operator <<(std::ostream &stream, const CompactSequence &compact_seq)
{
    Sequence seq(compact_seq);
    return stream << seq;
}

const CompactSequence &CompactSequence::Assign(const ShortSequence &short_seq)
{
    resize(short_seq.size());
    for (unsigned i = 0; i < size(); ++i)
        set_base(i, short_seq[i]);
    return *this;
}

const CompactSequence &CompactSequence::Assign(const Kmer &kmer)
{
    resize(kmer.size());
    for (unsigned i = 0; i < size(); ++i)
        set_base(i, kmer[i]);
    return *this;
}

const CompactSequence &CompactSequence::Append(const CompactSequence &compact_seq, int offset, size_t length)
{
    if (length == std::string::npos || length > compact_seq.size() - offset)
        length = compact_seq.size() - offset;

    uint32_t old_size = size();
    resize(old_size + length);

    for (unsigned i = 0; i < length; ++i)
        set_base(i + old_size, compact_seq[i + offset]);
    return *this;
}

const CompactSequence &CompactSequence::Append(const Sequence &seq, int offset, size_t length)
{
    if (length == std::string::npos || length > seq.size() - offset)
        length = seq.size() - offset;

    uint32_t old_size = size();
    resize(old_size + length);

    for (unsigned i = 0; i < length; ++i)
        set_base(i + old_size, seq[i + offset]);
    return *this;
}

const CompactSequence &CompactSequence::Append(uint8_t ch)
{
    resize(size() + 1);
    set_base(size()-1, ch);
    return *this;
}

const CompactSequence &CompactSequence::ReverseComplement()
{
    uint8_t remainder = data_[data_.size()-1];
    reverse(data_.begin(), data_.begin() + data_.size() - 1);
    for (unsigned i = 0; i+1 < data_.size(); ++i)
        bit_operation::ReverseComplement((uint8_t &)data_[i]);
    for (unsigned i = 0; i+2 < data_.size(); ++i)
        data_[i] = (data_[i] >> (remainder<<1)) | (data_[i+1] << ((3-remainder)<<1));
    data_[data_.size()-2] >>= (remainder << 1);
    return *this;
}

Kmer CompactSequence::GetKmer(uint32_t offset, uint32_t kmer_size) const
{
    Kmer kmer(kmer_size);
    for (unsigned i = 0; i < kmer_size; ++i)
        kmer.set_base(i, get_base(i + offset));
    return kmer;
}

