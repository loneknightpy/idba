/**
 * @file short_sequence.h
 * @brief Short Sequence Class.  
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#ifndef __SEQUENCE_SHORT_SEQUENCE_H_

#define __SEQUENCE_SHORT_SEQUENCE_H_

#include <stdint.h>

#include <algorithm>
#include <cstring>
#include <string>


class Sequence;
class CompactSequence;
class Kmer;

/**
 * @brief It store a short sequence in a very space effient way. The length 
 * limit can be found by max_size().
 */
class ShortSequence
{
public:
    ShortSequence() { data_[kNumBytes-2] = data_[kNumBytes-1] = 0; }
    ShortSequence(const ShortSequence &short_seq) 
    { Assign(short_seq); }
    explicit ShortSequence(const Sequence &seq, int offset = 0, size_t length = std::string::npos) 
    { Assign(seq, offset, length); }
    explicit ShortSequence(const CompactSequence &compact_seq, int offset = 0, size_t length = std::string::npos)
    { Assign(compact_seq, offset, length); }
    explicit ShortSequence(const Kmer &kmer)
    { Assign(kmer); }

    const ShortSequence &operator =(const ShortSequence &short_seq) { Assign(short_seq); return *this; }
    const ShortSequence &operator =(const Sequence &seq) { Assign(seq); return *this; }
    const ShortSequence &operator =(const CompactSequence &compact_seq) { Assign(compact_seq); return *this; }
    const ShortSequence &operator =(const Kmer &kmer) { Assign(kmer); return *this; }

    const ShortSequence &Assign(const ShortSequence &short_seq)
    { std::memcpy(data_, short_seq.data_, kNumBytes); return *this; }
    const ShortSequence &Assign(const Sequence &seq, int offset = 0, size_t length = std::string::npos);
    const ShortSequence &Assign(const CompactSequence &seq, int offset = 0, size_t length = std::string::npos);
    const ShortSequence &Assign(const Kmer &kmer);

    uint8_t operator [](uint32_t index) const
    { return (data_[index>>2] >> ((index&3) << 1)) & 3; }
    uint8_t get_base(uint32_t index) const
    { return (data_[index>>2] >> ((index&3) << 1)) & 3; }
    void set_base(uint32_t index, uint8_t ch)
    { data_[index>>2] = (data_[index>>2] & ~(3ULL << ((index&3) << 1))) | ((ch&3) << ((index&3) << 1)); }

    void swap(ShortSequence &short_seq)
    { 
        if (this != &short_seq)
        {
            for (unsigned i = 0; i < kNumBytes; ++i)
                std::swap(data_[i], short_seq.data_[i]);
        }
    }

    uint32_t size() const 
    { return data_[kNumBytes-2] + (data_[kNumBytes-1] << 8); }
    void resize(uint32_t new_size) 
    { data_[kNumBytes-2] = (new_size & ((1<<8)-1)); data_[kNumBytes-1] = (new_size >> 8); }
    bool empty() const 
    { return data_[kNumBytes-2] == 0 && data_[kNumBytes-1] == 0; }
    static uint32_t max_size() { return kMaxShortSequence; }

    bool operator ==(const ShortSequence &seq) const
    {
        if (size() != seq.size())
            return false;
        for (unsigned i = 0; i < size(); ++i)
        {
            if ((*this)[i] != seq[i])
                return false;
        }
        return true;
    }

    bool operator !=(const ShortSequence &seq) const
    { return !(*this == seq); }

    bool operator <(const ShortSequence &seq)
    {
        int len = std::min(size(), seq.size());
        for (int i = 0; i < len; ++i)
        {
            if ((*this)[i] != seq[i])
                return (*this)[i] < seq[i];
        }
        return size() < seq.size();
    }

    static const uint32_t kMaxShortSequence = 128;
    static const uint32_t kNumBytes = (kMaxShortSequence + 3) / 4 + 2;

private:
    uint8_t data_[kNumBytes];
};

namespace std
{
template <> inline void swap(ShortSequence &short_seq1, ShortSequence &short_seq2)
{ short_seq1.swap(short_seq2); }
}

#endif

