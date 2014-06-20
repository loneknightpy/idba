/**
 * @file hash_set_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-24
 */


#include "gtest/gtest.h"

#include "container/hash_set.h"

#include <set>

#include "basic/kmer.h"
#include "sequence/sequence.h"


using namespace std;

class HashSetTest: public testing::Test
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
            simple_hash_set.insert(kmer);
            simple_set.insert(kmer);
        }
    }

    HashSet<Kmer> simple_hash_set;
    set<Kmer> simple_set;
    Sequence seq;
};

TEST_F(HashSetTest, DefaultConstructor)
{
    HashSet<Kmer> hash_set;
    EXPECT_EQ(0, int(hash_set.size()));
}

TEST_F(HashSetTest, CopyConstructor)
{
    HashSet<Kmer> hash_set(simple_hash_set);
    HashSet<Kmer>::iterator p = simple_hash_set.begin();
    HashSet<Kmer>::iterator q = hash_set.begin();

    while (p != simple_hash_set.end())
    {
        EXPECT_EQ(*p, *q);
        ++p;
        ++q;
    }

    EXPECT_EQ(hash_set.end(), q);
}

TEST_F(HashSetTest, Iterator)
{
    HashSet<Kmer>::iterator iter = simple_hash_set.begin();
    while (iter != simple_hash_set.end())
    {
        EXPECT_NE(simple_set.end(), simple_set.find(*iter));
        ++iter;
    }
}

TEST_F(HashSetTest, Insert)
{
    HashSet<Kmer> hash_set;
    for (set<Kmer>::iterator iter = simple_set.begin(); iter != simple_set.end(); ++iter)
    {
        pair<HashSet<Kmer>::iterator, bool> result = hash_set.insert(*iter);
        EXPECT_EQ(true, result.second);
        EXPECT_EQ(*iter, *result.first);
    }

    for (set<Kmer>::iterator iter = simple_set.begin(); iter != simple_set.end(); ++iter)
    {
        pair<HashSet<Kmer>::iterator, bool> result = hash_set.insert(*iter);
        EXPECT_EQ(false, result.second);
        EXPECT_EQ(*iter, *result.first);
    }

    EXPECT_EQ(simple_set.size(), hash_set.size());
}

TEST_F(HashSetTest, Find)
{
    for (set<Kmer>::iterator iter = simple_set.begin(); iter != simple_set.end(); ++iter)
        EXPECT_NE(simple_hash_set.end(), simple_hash_set.find(*iter));

    for (HashSet<Kmer>::iterator iter = simple_hash_set.begin(); iter != simple_hash_set.end(); ++iter)
        EXPECT_NE(simple_set.end(), simple_set.find(*iter));

    Kmer kmer;
    EXPECT_EQ(simple_hash_set.end(), simple_hash_set.find(kmer));
}

TEST_F(HashSetTest, Remove)
{
    for (set<Kmer>::iterator iter = simple_set.begin(); iter != simple_set.end(); ++iter)
    {
        uint32_t num_removed_nodes = simple_hash_set.remove(*iter);
        EXPECT_EQ(1, int(num_removed_nodes));
    }

    for (set<Kmer>::iterator iter = simple_set.begin(); iter != simple_set.end(); ++iter)
    {
        uint32_t num_removed_nodes = simple_hash_set.remove(*iter);
        EXPECT_EQ(0, int(num_removed_nodes));
    }

    EXPECT_EQ(0, int(simple_hash_set.size()));
}

TEST_F(HashSetTest, Size)
{
    EXPECT_EQ(simple_set.size(), simple_hash_set.size());
}

TEST_F(HashSetTest, Empty)
{
    EXPECT_EQ(false, simple_hash_set.empty());
    simple_hash_set.clear();
    EXPECT_EQ(true, simple_hash_set.empty());
}

TEST_F(HashSetTest, Swap)
{
    HashSet<Kmer> hash_set;
    Kmer kmer;
    hash_set.insert(kmer);

    hash_set.swap(simple_hash_set);
    EXPECT_EQ(hash_set.end(), hash_set.find(kmer));
    EXPECT_NE(simple_hash_set.end(), simple_hash_set.find(kmer));
}

TEST_F(HashSetTest, Clear)
{
    simple_hash_set.clear();
    EXPECT_EQ(0, int(simple_hash_set.size()));
}

