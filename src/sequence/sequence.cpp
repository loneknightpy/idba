/**
 * @file sequence.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#include "sequence/sequence.h"

#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

#include "basic/kmer.h"
#include "sequence/compact_sequence.h"
#include "sequence/short_sequence.h"


using namespace std;

istream &operator >>(istream &is, Sequence &seq)
{
    getline(is, seq.bases_);
    if (!is)
        return is;

    string line;
    while (is && (isalnum(is.peek()) || is.peek() == '\n') && getline(is, line))
        seq.bases_ += line;
    is.clear();

    seq.Encode();

    return is;
}

ostream &operator <<(ostream &os, const Sequence &seq)
{
    Sequence tmp = seq;
    tmp.Decode();
    return os << tmp.bases_;
}

const Sequence &Sequence::Assign(const CompactSequence &compact_seq, int offset, size_t length)
{
    if (length == string::npos || length > compact_seq.size() - offset)
        length = compact_seq.size() - offset;

    bases_.resize(length);
    for (unsigned i = 0; i < bases_.size(); ++i)
        bases_[i] = compact_seq[i + offset];
    return *this;
}

const Sequence &Sequence::Assign(const ShortSequence &short_seq)
{
    bases_.resize(short_seq.size());
    for (unsigned i = 0; i < bases_.size(); ++i)
        bases_[i] = short_seq[i];
    return *this;
}

const Sequence &Sequence::Assign(const Kmer &kmer)
{
    bases_.resize(kmer.size());
    for (unsigned i = 0; i < bases_.size(); ++i)
        bases_[i] = kmer[i];
    return *this;
}

const Sequence &Sequence::ReverseComplement()
{
    reverse(bases_.begin(), bases_.end());
    for (unsigned i = 0; i < bases_.size(); ++i)
    {
        switch(bases_[i])
        {
            case 'a':
            case 'A':
                bases_[i] = 'T';
                break;

            case 'c':
            case 'C':
                bases_[i] = 'G';
                break;

            case 'g':
            case 'G':
                bases_[i] = 'C';
                break;

            case 't':
            case 'T':
                bases_[i] = 'A';
                break;

            case 0:
            case 1:
            case 2:
            case 3:
                bases_[i] = 3 - bases_[i];
                break;

            case 'N':
            case 4:
                break;

            default:
                throw logic_error("reverse complement error: unkown character.");
        }
    }

    return *this;
}

bool Sequence::IsValid() const
{
    for (unsigned i = 0; i < bases_.size(); ++i)
    {
        if (!IsValid(bases_[i]))
            return false;
    }
    return true;
}

bool Sequence::IsPalindrome() const
{
    unsigned half = (bases_.size() + 1)/2;
    for (unsigned i = 0; i < half; ++i)
    {
        if (bases_[i] + bases_[bases_.size()-1 - i] != 3)
            return false;
    }
    return true;
}

void Sequence::TrimN()
{
    int len = size();
    while (len > 0 && !IsValid(bases_[len-1]))
        --len;
    bases_.resize(len);
}

Kmer Sequence::GetKmer(uint32_t offset, uint32_t kmer_size) const
{
    Kmer kmer(kmer_size);
    for (unsigned i = 0; i < kmer_size; ++i)
        kmer.set_base(i, bases_[offset + i]);
    return kmer;
}

void Sequence::Encode()
{
    for (unsigned i = 0; i < bases_.size(); ++i)
    {
        switch(bases_[i])
        {
            case 'a':
            case 'A':
                bases_[i] = 0;
                break;

            case 'c':
            case 'C':
                bases_[i] = 1;
                break;

            case 'g':
            case 'G':
                bases_[i] = 2;
                break;

            case 't':
            case 'T':
                bases_[i] = 3;
                break;

            case 'n':
            case 'N':
                bases_[i] = 4;
                break;

            default:
                bases_[i] = 4;
                break;
                //throw logic_error("encode error: unkown character.");
        }
    }
}

void Sequence::Decode()
{
    for (unsigned i = 0; i < bases_.size(); ++i)
    {
        switch(bases_[i])
        {
            case 0:
                bases_[i] = 'A';
                break;

            case 1:
                bases_[i] = 'C';
                break;

            case 2:
                bases_[i] = 'G';
                break;

            case 3:
                bases_[i] = 'T';
                break;

            case 4:
                bases_[i] = 'N';
                break;

            default:
                throw logic_error("decode error: unkown character.");
        }
    }
}

istream &ReadFasta(istream &is, Sequence &seq, string &comment)
{
    string line;
    getline(is, line);

    if (!is)
        return is;

    comment = line.substr(1);

    return is >> seq;
}

ostream &WriteFasta(ostream &os, const Sequence &seq, const string &comment)
{
    return os << ">" << comment << "\n" << seq << "\n";
}

istream &ReadFastq(istream &is, Sequence &seq, string &comment, string &quality)
{
    string line;
    getline(is, line);
    if (!is)
        return is;

    comment = line.substr(1);

    comment = line.substr(1);
    is >> seq;
    getline(is, line);

    quality = "";
    getline(is, line);
    quality = line;
//    while (is && is.peek() != '@' && getline(is, line))
//        quality += line;
    is.clear();

    return is;
}

ostream &WriteFastq(ostream &os, const Sequence &seq, const string &comment, const string &quality)
{
    return os << "@" << comment << "\n" << seq << "\n" << "+" << "\n" << quality << "\n";
}

