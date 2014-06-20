/**
 * @file blat_record.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.3
 * @date 2011-09-06
 */

#include "misc/blat_record.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "misc/utils.h"

using namespace std;

void BlatRecord::Parse(const string &record)
{
    stringstream ss(record);
    string tmp;
    string strand;
    string block_size_buf;
    string query_from_buf;
    string ref_from_buf;
    int block_count = 0;
    
    ss >> match_count;
    ss >> mismatch_count;
    ss >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;
    ss >> strand;
    ss >> query_name;
    ss >> query_length;
    ss >> query_from;
    ss >> query_to;
    ss >> ref_name;
    ss >> ref_length;
    ss >> ref_from;
    ss >> ref_to;
    ss >> block_count;
    ss >> block_size_buf;
    ss >> query_from_buf;
    ss >> ref_from_buf;

    if (strand == "+")
        is_reverse = false;
    else
        is_reverse = true;

    Replace(block_size_buf, ',', ' ');
    Replace(query_from_buf, ',', ' ');
    Replace(ref_from_buf, ',', ' ');

    stringstream block_size_stream(block_size_buf);
    stringstream query_from_stream(query_from_buf);
    stringstream ref_from_stream(ref_from_buf);

    blocks.resize(block_count);
    for (unsigned i = 0; i < blocks.size(); ++i)
    {
        block_size_stream >> blocks[i].size;
        query_from_stream >> blocks[i].query_from;
        ref_from_stream >> blocks[i].ref_from;
    }
}

