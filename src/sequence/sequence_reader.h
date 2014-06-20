/**
 * @file sequence_reader.h
 * @brief Sequence Reader Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#ifndef __SEQUENCE_SEQUENCE_READER_H_

#define __SEQUENCE_SEQUENCE_READER_H_

#include <fstream>
#include <istream>
#include <string>


class Sequence;

/**
 * @brief It is a abstract reader class which provides simple interface to 
 * read sequence from different format of files.
 */
class SequenceReader
{
public:
    SequenceReader() { is_ = NULL, is_own_ = false; }
    explicit SequenceReader(std::istream &is, bool is_own = false);
    virtual ~SequenceReader() { Dispose(); }

    virtual void Dispose() { if (is_own_) delete is_; is_ = NULL; }

    bool Read(Sequence &seq);
    bool Read(Sequence &seq, std::string &comment);
    bool Read(Sequence &seq, std::string &comment, std::string &quality);

protected:
    SequenceReader(const SequenceReader &);
    const SequenceReader &operator =(const SequenceReader &);

    virtual bool ReadRecord(Sequence &seq, std::string &comment, std::string &quality) = 0;

    std::istream *is_; 
    bool is_own_;
};

/**
 * @brief It is the reader class for Fasta file.
 */
class FastaReader: public SequenceReader
{
public:
    FastaReader() {}
    explicit FastaReader(std::istream &is, bool is_own = false)
        : SequenceReader(is, false) {}
    explicit FastaReader(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary)
        : SequenceReader(*(new std::ifstream(filename.c_str(), mode)), true) {}

protected:
    void Open(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary)
    {
        Dispose();
        is_ = new std::ifstream(filename.c_str(), mode);
        is_own_ = true;
    }

    virtual bool ReadRecord(Sequence &seq, std::string &comment, std::string &quality);

private:
    FastaReader(const FastaReader &);
    const FastaReader &operator = (const FastaReader &);
};

/**
 * @brief It is the reader class for Fastq file.
 */
class FastqReader: public SequenceReader
{
public:
    explicit FastqReader(std::istream &is, bool is_own = false)
        : SequenceReader(is, false) {}
    explicit FastqReader(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary)
        : SequenceReader(*(new std::ifstream(filename.c_str(), mode)), true) {}

protected:
    virtual bool ReadRecord(Sequence &seq, std::string &comment, std::string &quality);

private:
    FastqReader(const FastqReader &);
    const FastqReader &operator =(const FastqReader &);
};

#endif

