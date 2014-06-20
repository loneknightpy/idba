/**
 * @file sequence_io.h
 * @brief Useful functions for sequence I/O.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-06
 */

#ifndef __SEQUENCE_SEQUENCE_IO_H_

#define __SEQUENCE_SEQUENCE_IO_H_

#include <stdint.h>

#include <string>
#include <deque>

#include "sequence/sequence_reader.h"
#include "sequence/sequence_writer.h"

class Sequence;
class ShortSequence;
class FastaWriter;
class ShortReadLibrary;
class LongReadLibrary;

uint64_t ReadLibrary(const std::string &filename, ShortReadLibrary &library);
uint64_t ReadLibrary(const std::string &filename, LongReadLibrary &library);

uint64_t ReadSequence(const std::string &filename, std::deque<Sequence> &sequences);
uint64_t ReadSequence(const std::string &filename, std::deque<Sequence> &sequences, std::deque<std::string> &names);

uint64_t ReadSequence(const std::string &filename, std::deque<ShortSequence> &sequences);
uint64_t ReadSequence(const std::string &filename, std::deque<ShortSequence> &sequences, std::deque<std::string> &names);

uint64_t WriteSequence(const std::string &filename, const std::deque<Sequence> &sequences, const std::string &prefix = "");
uint64_t WriteSequence(FastaWriter &writer, const std::deque<Sequence> &sequences, const std::string &prefix = "");

uint64_t WriteSequence(const std::string &filename, const std::deque<ShortSequence> &sequences, const std::string &prefix = "");
uint64_t WriteSequence(FastaWriter &writer, const std::deque<ShortSequence> &sequences, const std::string &prefix = "");

#endif

