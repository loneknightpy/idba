/**
 * @file managed_list_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-25
 */

#include "gtest/gtest.h"

#include "container/managed_list.h"

#include <list>

#include "basic/kmer.h"


using namespace std;

class ManagedListTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        simple_list.set_pool(pool);

        for (int i = 0; i < 100; ++i)
        {
            simple_list.push_front(i);
            stl_list.push_front(i);
        }
    }

    ManagedList<int>::node_pool_type pool;
    ManagedList<int> simple_list;
    list<int> stl_list;
};

TEST_F(ManagedListTest, DefaultConstructor)
{
    ManagedList<int> list;
    EXPECT_EQ(0, (int)list.size());
    EXPECT_EQ(list.begin(), list.end());
}

TEST_F(ManagedListTest, ConstructorWithPool)
{
    ManagedList<int> list(pool);
    EXPECT_EQ(0, (int)list.size());
    EXPECT_EQ(list.begin(), list.end());
}

TEST_F(ManagedListTest, CopyConstructor)
{
    ManagedList<int> list(simple_list);
    EXPECT_EQ(simple_list.size(), list.size());
    ManagedList<int>::iterator p = simple_list.begin();
    ManagedList<int>::iterator q = list.begin();

    while (p != simple_list.end())
    {
        EXPECT_EQ(*p, *q);
        ++p;
        ++q;
    }

    EXPECT_EQ(q, list.end());
}

TEST_F(ManagedListTest, PushFront)
{
    ManagedList<int> list(pool);
    for (int i = 0; i < 100; ++i)
    {
        list.push_front(i);
        EXPECT_EQ(i, list.front());
    }
    EXPECT_EQ(100, int(list.size()));
}

TEST_F(ManagedListTest, PopFront)
{
    for (int i = 0; i < 100; ++i)
    {
        EXPECT_EQ(stl_list.front(), simple_list.front());
        stl_list.pop_front();
        simple_list.pop_front();
    }
}

TEST_F(ManagedListTest, Remove)
{
    simple_list.remove(10);
//    for (ManagedList<int>::iterator iter = simple_list.begin(); iter != simple_list.end(); ++iter)
//        EXPECT_NE(10, *iter);
//    EXPECT_EQ(99, int(simple_list.size()));
}

TEST_F(ManagedListTest, Erase)
{
    simple_list.erase(simple_list.find(10));
    for (ManagedList<int>::iterator iter = simple_list.begin(); iter != simple_list.end(); ++iter)
        EXPECT_NE(10, *iter);
    EXPECT_EQ(99, int(simple_list.size()));
}

TEST_F(ManagedListTest, Find)
{
    for (int i = 0; i < 100; ++i)
        EXPECT_NE(simple_list.end(), simple_list.find(i));
    EXPECT_EQ(simple_list.end(), simple_list.find(100));
}

TEST_F(ManagedListTest, Swap)
{
    ManagedList<int> list(pool);
    list.push_front(100);
    list.swap(simple_list);
    EXPECT_EQ(1, int(simple_list.size()));
    EXPECT_EQ(simple_list.begin(), simple_list.find(100));
    EXPECT_EQ(list.end(), list.find(100));
    EXPECT_EQ(100, int(list.size()));
}

TEST_F(ManagedListTest, Size)
{
    ManagedList<int> list(pool);
    for (int i = 0; i < 100; ++i)
    {
        list.push_front(i);
        EXPECT_EQ(i+1, (int)list.size());
    }
}

TEST_F(ManagedListTest, Empty)
{
    ManagedList<int> list(pool);
    EXPECT_EQ(true, list.empty());
    list.push_front(100);
    EXPECT_EQ(false, list.empty());
}

TEST_F(ManagedListTest, Clear)
{
    simple_list.clear();
    EXPECT_EQ(0, (int)simple_list.size());
    EXPECT_EQ(true, simple_list.empty());
    EXPECT_EQ(simple_list.end(), simple_list.end());
}

