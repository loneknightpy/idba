/**
 * @file hash_map_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-24
 */

#include "gtest/gtest.h"

#include "container/hash_map.h"

#include <map>

#include "basic/kmer.h"
#include "sequence/sequence.h"


using namespace std;

class HashMapTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        for (int i = 0; i < 10000; ++i)
            seq.Append(rand() & 3);

        Kmer kmer(25);
        for (uint32_t i = 0; i < seq.size(); ++i)
        {
            kmer.ShiftAppend(seq[i]);
            simple_hash_map[kmer]++;
            simple_map[kmer]++;
        }
    }

    HashMap<Kmer, int> simple_hash_map;
    map<Kmer, int> simple_map;
    Sequence seq;
};

TEST_F(HashMapTest, DefaultConstructor)
{
    HashMap<Kmer, int> hash_map;
    EXPECT_EQ(0, int(hash_map.size()));
}

TEST_F(HashMapTest, CopyConstructor)
{
    HashMap<Kmer, int> hash_map(simple_hash_map);
    HashMap<Kmer, int>::iterator p = simple_hash_map.begin();
    HashMap<Kmer, int>::iterator q = hash_map.begin();

    while (p != simple_hash_map.end())
    {
        EXPECT_EQ(p->first, q->first);
        EXPECT_EQ(p->second, q->second);
        ++p;
        ++q;
    }

    EXPECT_EQ(hash_map.end(), q);
}

TEST_F(HashMapTest, Iterator)
{
    HashMap<Kmer, int>::iterator iter = simple_hash_map.begin();
    while (iter != simple_hash_map.end())
    {
        EXPECT_NE(simple_map.end(), simple_map.find(iter->first));
        EXPECT_EQ(simple_map[iter->first], iter->second);
        ++iter;
    }
}

TEST_F(HashMapTest, Insert)
{
    HashMap<Kmer, int> hash_map;
    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        pair<HashMap<Kmer, int>::iterator, bool> result = hash_map.insert(*iter);
        EXPECT_EQ(true, result.second);
        EXPECT_EQ(iter->first, result.first->first);
    }

    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        pair<HashMap<Kmer, int>::iterator, bool> result = hash_map.insert(*iter);
        EXPECT_EQ(false, result.second);
        EXPECT_EQ(iter->first, result.first->first);
    }

    EXPECT_EQ(simple_map.size(), hash_map.size());
}

TEST_F(HashMapTest, Find)
{
    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        EXPECT_NE(simple_hash_map.end(), simple_hash_map.find(iter->first));
        EXPECT_EQ(simple_hash_map[iter->first], iter->second);
    }

    for (HashMap<Kmer, int>::iterator iter = simple_hash_map.begin(); iter != simple_hash_map.end(); ++iter)
    {
        EXPECT_NE(simple_map.end(), simple_map.find(iter->first));
        EXPECT_EQ(simple_map[iter->first], iter->second);
    }

    Kmer kmer;
    EXPECT_EQ(simple_hash_map.end(), simple_hash_map.find(kmer));
}

TEST_F(HashMapTest, Remove)
{
    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        uint32_t num_removed_nodes = simple_hash_map.remove(iter->first);
        EXPECT_EQ(1, int(num_removed_nodes));
    }

    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        uint32_t num_removed_nodes = simple_hash_map.remove(iter->first);
        EXPECT_EQ(0, int(num_removed_nodes));
    }

    EXPECT_EQ(0, int(simple_hash_map.size()));
}

TEST_F(HashMapTest, Size)
{
    EXPECT_EQ(simple_map.size(), simple_hash_map.size());
}

TEST_F(HashMapTest, Empty)
{
    EXPECT_EQ(false, simple_hash_map.empty());
    simple_hash_map.clear();
    EXPECT_EQ(true, simple_hash_map.empty());
}

TEST_F(HashMapTest, Swap)
{
    HashMap<Kmer, int> hash_map;
    Kmer kmer;
    hash_map[kmer]++;

    hash_map.swap(simple_hash_map);
    EXPECT_EQ(hash_map.end(), hash_map.find(kmer));
    EXPECT_NE(simple_hash_map.end(), simple_hash_map.find(kmer));
}

TEST_F(HashMapTest, Clear)
{
    simple_hash_map.clear();
    EXPECT_EQ(0, int(simple_hash_map.size()));
}

