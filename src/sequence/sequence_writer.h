/**
 * @file sequence_writer.h
 * @brief Sequence Writer Class. 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#ifndef __MISC_SEQUENCE_WRITER_H_

#define __MISC_SEQUENCE_WRITER_H_

#include <fstream>
#include <istream>
#include <stdexcept>
#include <string>


class Sequence;

/**
 * @brief It is a abstract writer class which provides simple interface to 
 * write sequence into different format of files.
 */
class SequenceWriter
{
public:
    SequenceWriter() { os_ = NULL; is_own_ = false; }
    SequenceWriter(std::ostream &os, bool is_own = false);
    virtual ~SequenceWriter() { Dispose(); }

    virtual void Dispose() { if (is_own_) delete os_; os_ = NULL; }

    bool Write(const Sequence &seq);
    bool Write(const Sequence &seq, const std::string &comment);
    bool Write(const Sequence &seq, const std::string &comment, const std::string &quality);

protected:
    SequenceWriter(const SequenceWriter &);
    const SequenceWriter &operator =(const SequenceWriter &);

    virtual bool WriteRecord(const Sequence &seq, const std::string &comment, const std::string &quality) = 0;

    std::ostream *os_; 
    bool is_own_;
};

/**
 * @brief It is the writer class for Fasta file.
 */
class FastaWriter: public SequenceWriter
{
public:
    FastaWriter() {}
    explicit FastaWriter(std::ostream &os, bool is_own = false)
        : SequenceWriter(os, is_own) {}
    explicit FastaWriter(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary)
        : SequenceWriter(*(new std::ofstream(filename.c_str(), mode)), true) {}

    void Open(const std::string &filename,
            std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary)
    {
        Dispose();
        os_ = new std::ofstream(filename.c_str(), mode);
        is_own_ = true;
    }
protected:

    virtual bool WriteRecord(const Sequence &seq, const std::string &comment, const std::string &quality);

private:
    FastaWriter(const FastaWriter &);
    const FastaWriter &operator =(const FastaWriter &);
};

/**
 * @brief It is the writer class for Fastq file.
 */
class FastqWriter: public SequenceWriter
{
public:
    explicit FastqWriter(std::ostream &os, bool is_own = false)
        : SequenceWriter(os, is_own) {}
    explicit FastqWriter(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary)
        : SequenceWriter(*(new std::ofstream(filename.c_str(), mode)), true) {}

protected:
    virtual bool WriteRecord(const Sequence &seq, const std::string &comment, const std::string &quality);

private:
    FastqWriter(const FastqWriter &);
    const FastqWriter &operator =(const FastqWriter &);
};

#endif

