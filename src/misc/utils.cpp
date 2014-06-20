/**
 * @file utils.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-06
 */

#include "misc/utils.h"

#include <algorithm>
#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <dirent.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "sequence/sequence.h"

using namespace std;

static const uint32_t kMaxLine = (1 << 20);
static char buffer[kMaxLine];

string FormatString(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsprintf(buffer, fmt, ap);
    va_end(ap);

    return buffer;
}

bool IsExist(const string &filename)
{
    FILE *fp = fopen(filename.c_str(), "rb");
    if (fp == NULL)
    {
        return false;
    }
    else
    {
        fclose(fp);
        return true;
    }
}

FILE *OpenFile(const string &filename, const string &mode)
{
    FILE *fp = fopen(filename.c_str(), mode.c_str());

    if (fp == NULL)
    {
        fprintf(stderr, "open %s failed\n", filename.c_str());
        exit(1);
    }

    return fp;
}

void MakeDir(const string &directory)
{
    DIR *dir = opendir(directory.c_str());
    if (dir == NULL)
    {
        if (mkdir(directory.c_str(), 0777) == -1)
            throw invalid_argument("can't open directory");
    }
}

void Replace(string &s, int a, int b)
{
    for (unsigned i = 0; i < s.size(); ++i)
    {
        if (s[i] == a)
            s[i] = b;
    }
}

void ToLower(string &s)
{
    for (unsigned i = 0; i < s.size(); ++i)
        s[i] = tolower(s[i]);
}

void PrintN50(const deque<Sequence> &contigs, int min_contig, int ref_length)
{
    int64_t sum = 0;
    deque<int> lengths;
    for (uint32_t i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() >= min_contig)
        {
            lengths.push_back(contigs[i].size());
            sum += contigs[i].size();
        }
    }

    if (ref_length != 0)
        sum = ref_length;

    if (lengths.size() == 0)
        lengths.push_back(0);

    sort(lengths.begin(), lengths.end());
    reverse(lengths.begin(), lengths.end());

    int64_t partial_sum = 0;
    int n50 = 0;
    int n80 = 0;
    for (uint32_t i = 0; i < lengths.size(); ++i)
    {
        partial_sum += lengths[i];
        if (n50 == 0 && partial_sum > sum * 0.5)
            n50 = lengths[i];
        if (n80 == 0 && partial_sum > sum * 0.8)
            n80 = lengths[i];
    }

    cout << FormatString("contigs: %d n50: %d max: %d mean: %lld total length: %lld n80: %d",
           int(lengths.size()), n50, (int)lengths[0], (long long)(sum/lengths.size()), (long long)sum, n80) << endl;
    cout.flush();
}

