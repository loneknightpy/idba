/**
 * @file utils.h
 * @brief Utility functions.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-06
 */

#ifndef __MISC_UTILS_H_

#define __MISC_UTILS_H_

#include <cstdio>
#include <deque>
#include <string>

class Sequence;
class Sequence;

std::string FormatString(const char *fmt, ...);

bool IsExist(const std::string &filename);
FILE *OpenFile(const std::string &filename, const std::string &mode);
void MakeDir(const std::string &diretory);
void Replace(std::string &s, int a, int b);
void ToLower(std::string &s);
void PrintN50(const std::deque<Sequence> &contigs, int min_contig = 0, int ref_length = 0);

#endif

