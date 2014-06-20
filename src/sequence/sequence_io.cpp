/**
 * @file sequence_io.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-06
 */

#include "sequence/sequence_io.h"

#include <string>

#include "misc/utils.h"
#include "sequence/read_library.h"
#include "sequence/sequence.h"
#include "sequence/sequence_reader.h"
#include "sequence/sequence_writer.h"
#include "sequence/short_sequence.h"

using namespace std;

uint64_t ReadLibrary(const std::string &filename, ShortReadLibrary &library)
{
    library.clear();
    FastaReader reader(filename);
    Sequence seq;
    while (reader.Read(seq))
    {
        seq.TrimN();
        ShortSequence short_seq(seq);
        library.reads_.push_back(short_seq);
    }
    return library.size();
}

uint64_t ReadLibrary(const std::string &filename, LongReadLibrary &library)
{
    library.clear();
    FastaReader reader(filename);
    Sequence seq;
    while (reader.Read(seq))
        library.reads_.push_back(seq);
    return library.size();
}

uint64_t ReadSequence(const string &filename, deque<Sequence> &sequences)
{
    sequences.clear();
    FastaReader reader(filename);
    Sequence seq;
    while (reader.Read(seq))
        sequences.push_back(seq);
    return sequences.size();
}

uint64_t ReadSequence(const string &filename, deque<Sequence> &sequences, deque<string> &names)
{
    sequences.clear();
    FastaReader reader(filename);
    Sequence seq;
    string name;
    while (reader.Read(seq, name))
    {
        sequences.push_back(seq);
        names.push_back(name);
    }
    return sequences.size();
}

uint64_t ReadSequence(const string &filename, deque<ShortSequence> &sequences)
{
    sequences.clear();
    FastaReader reader(filename);
    Sequence seq;
    while (reader.Read(seq))
    {
        seq.TrimN();
        ShortSequence short_seq(seq);
        sequences.push_back(short_seq);
    }
    return sequences.size();
}

uint64_t WriteSequence(const string &filename, const deque<Sequence> &sequences, const string &prefix)
{
    FastaWriter writer(filename);
    for (uint64_t i = 0; i < sequences.size(); ++i)
        writer.Write(sequences[i], FormatString("%s_%d", prefix.c_str(), i));
    return sequences.size();
}

uint64_t WriteSequence(FastaWriter &writer, const deque<Sequence> &sequences, const string &prefix)
{
    for (uint64_t i = 0; i < sequences.size(); ++i)
        writer.Write(sequences[i], FormatString("%s_%d", prefix.c_str(), i));
    return sequences.size();
}

uint64_t WriteSequence(const string &filename, const deque<ShortSequence> &sequences, const string &prefix)
{
    FastaWriter writer(filename);
    return WriteSequence(writer, sequences, prefix);
}

uint64_t WriteSequence(FastaWriter &writer, const deque<ShortSequence> &sequences, const string &prefix)
{
    for (uint64_t i = 0; i < sequences.size(); ++i)
    {
        Sequence seq(sequences[i]);
        writer.Write(seq, FormatString("%s_%d", prefix.c_str(), i));
    }
    return sequences.size();
}

