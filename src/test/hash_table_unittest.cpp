/**
 * @file hash_table_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-04
 */

#include "gtest/gtest.h"

#include "container/hash_table.h"

#include <sstream>
#include <map>

#include "basic/kmer.h"
#include "sequence/sequence.h"


using namespace std;

class HashTableTest: public testing::Test
{
protected:
    struct Node
    {
        typedef Kmer key_type;

        Node(const Kmer &kmer = Kmer())
        { this->kmer = kmer; count = 0; }

        const key_type &key() const { return kmer; }
        void set_key(const key_type &key) { kmer = key; }

        Kmer kmer;
        int count;
    };

    struct Counter
    {
        Counter() { count = 0; }

        void operator ()(Node &node)
        {
#pragma omp atomic
            count += 1;
        }

        uint64_t count;
    };

    struct Equal
    {
        Equal(const Kmer &kmer) { this->kmer = kmer; }

        bool operator() (const Node &node) const { return kmer == node.kmer; }

        Kmer kmer;
    };

    virtual void SetUp()
    {
        for (int i = 0; i < 10000; ++i)
            seq.Append(rand() & 3);

        Kmer kmer(25);
        for (uint32_t i = 0; i < seq.size(); ++i)
        {
            kmer.ShiftAppend(seq[i]);
            Node &node = simple_hash_table.find_or_insert(Node(kmer));
            node.count++;
            ++simple_map[kmer];
        }
    }

    HashTable<Node, Kmer> simple_hash_table;
    map<Kmer, int> simple_map;
    Sequence seq;
};

TEST_F(HashTableTest, DefaultConstructor)
{
    HashTable<Node, Kmer> hash_table;
    EXPECT_EQ(0, int(hash_table.size()));
    uint64_t default_num_buckets = HashTable<Node, Kmer>::kDefaultNumBuckets;
    EXPECT_EQ(default_num_buckets, hash_table.bucket_count());
}

TEST_F(HashTableTest, CopyConstructor)
{
    HashTable<Node, Kmer> hash_table(simple_hash_table);
    HashTable<Node, Kmer>::iterator p = simple_hash_table.begin();
    HashTable<Node, Kmer>::iterator q = hash_table.begin();

    while (p != simple_hash_table.end())
    {
        EXPECT_EQ(p->kmer, q->kmer);
        EXPECT_EQ(p->count, q->count);
        ++p;
        ++q;
    }

    EXPECT_EQ(hash_table.end(), q);
}

TEST_F(HashTableTest, Iterator)
{
    HashTable<Node, Kmer>::iterator iter = simple_hash_table.begin();
    while (iter != simple_hash_table.end())
    {
        EXPECT_EQ(iter->count, simple_map[iter->kmer]);
        ++iter;
    }
}

TEST_F(HashTableTest, InsertUnique)
{
    HashTable<Node, Kmer> hash_table;
    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        pair<HashTable<Node, Kmer>::iterator, bool> result = hash_table.insert_unique(iter->first);
        EXPECT_EQ(true, result.second);
        EXPECT_EQ(iter->first, result.first->kmer);
    }

    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        pair<HashTable<Node, Kmer>::iterator, bool> result = hash_table.insert_unique(iter->first);
        EXPECT_EQ(false, result.second);
        EXPECT_EQ(iter->first, result.first->kmer);
    }

    EXPECT_EQ(simple_map.size(), hash_table.size());
}

TEST_F(HashTableTest, Find)
{
    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        HashTable<Node, Kmer>::iterator p = simple_hash_table.find(iter->first);
        EXPECT_EQ(iter->second, p->count);
    }

    for (HashTable<Node, Kmer>::iterator iter = simple_hash_table.begin(); iter != simple_hash_table.end(); ++iter)
        EXPECT_EQ(simple_map[iter->kmer], iter->count);

    Kmer kmer;
    EXPECT_EQ(simple_hash_table.end(), simple_hash_table.find(kmer));
}

TEST_F(HashTableTest, Remove)
{
    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        uint32_t num_removed_nodes = simple_hash_table.remove(iter->first);
        EXPECT_EQ(1, int(num_removed_nodes));
    }

    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        uint32_t num_removed_nodes = simple_hash_table.remove(iter->first);
        EXPECT_EQ(0, int(num_removed_nodes));
    }

    EXPECT_EQ(0, int(simple_hash_table.size()));
}

TEST_F(HashTableTest, RemoveIf)
{
    Kmer kmer = simple_map.begin()->first;
    simple_hash_table.remove_if(Equal(kmer));
    HashTable<Node, Kmer>::iterator p = simple_hash_table.find(kmer);
    EXPECT_EQ(simple_hash_table.end(), p);
}

TEST_F(HashTableTest, ForEach)
{
    Counter counter;
    simple_hash_table.for_each(counter);
    EXPECT_EQ(simple_hash_table.size(), counter.count);
}

TEST_F(HashTableTest, Size)
{
    EXPECT_EQ(simple_map.size(), simple_hash_table.size());
}

TEST_F(HashTableTest, Empty)
{
    EXPECT_EQ(false, simple_hash_table.empty());
    simple_hash_table.clear();
    EXPECT_EQ(true, simple_hash_table.empty());
}

TEST_F(HashTableTest, Swap)
{
    HashTable<Node, Kmer> hash_table;
    Kmer kmer;
    Node node(kmer);
    hash_table.insert_unique(node);

    hash_table.swap(simple_hash_table);
    EXPECT_EQ(hash_table.end(), hash_table.find(kmer));
    EXPECT_NE(simple_hash_table.end(), simple_hash_table.find(kmer));
}

TEST_F(HashTableTest, Clear)
{
    simple_hash_table.clear();
    EXPECT_EQ(0, int(simple_hash_table.size()));
}

TEST_F(HashTableTest, WriteAndRead)
{
    stringstream ss;
    HashTable<Node, Kmer> hash_table;

    ss << simple_hash_table;
    ss >> hash_table;

    for (map<Kmer, int>::iterator iter = simple_map.begin(); iter != simple_map.end(); ++iter)
    {
        HashTable<Node, Kmer>::iterator p = hash_table.find(iter->first);
        EXPECT_EQ(iter->second, p->count);
    }

    for (HashTable<Node, Kmer>::iterator iter = simple_hash_table.begin(); iter != simple_hash_table.end(); ++iter)
        EXPECT_EQ(iter->count, simple_map[iter->kmer]);
}

