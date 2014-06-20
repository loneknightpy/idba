/**
 * @file sequence_writer.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#include "sequence/sequence_writer.h"

#include <fstream>
#include <istream>
#include <string>

#include "sequence/sequence.h"


using namespace std;

SequenceWriter::SequenceWriter(ostream &os, bool is_own)
{
    if (os.fail())
        throw logic_error("SequenceWriter::SequenceWriter() ostream is invalid");

    os_ = &os;
    is_own_ = is_own;
}

bool SequenceWriter::Write(const Sequence &seq)
{
    string comment;
    string quality;
    return WriteRecord(seq, comment, quality);
}

bool SequenceWriter::Write(const Sequence &seq, const string &comment)
{
    string quality;
    return WriteRecord(seq, comment, quality);
}

bool SequenceWriter::Write(const Sequence &seq, const string &comment, const string &quality)
{
    return WriteRecord(seq, comment, quality);
}

bool FastaWriter::WriteRecord(const Sequence &seq, const string &comment, const string &quality)
{
    return WriteFasta(*os_, seq, comment);
}

bool FastqWriter::WriteRecord(const Sequence &seq, const string &comment, const string &quality)
{
    return WriteFastq(*os_, seq, comment, quality);
}

