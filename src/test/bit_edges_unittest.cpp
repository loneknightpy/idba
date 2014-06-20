/**
 * @file bit_edges_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-26
 */

#include "gtest/gtest.h"

#include "graph/bit_edges.h"

#include <algorithm>

using namespace std;


class BitEdgesTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        x = 7;
        y = 3;
    }

    BitEdges x;
    BitEdges y;
};


TEST_F(BitEdgesTest, DefaultConstructor)
{
    BitEdges edges;
    EXPECT_EQ(0, (int)edges);
    EXPECT_EQ(0, edges.size());
}

TEST_F(BitEdgesTest, CopyConstructor)
{
    BitEdges edges(x);
    EXPECT_EQ(x, edges);
}

TEST_F(BitEdgesTest, ConstructorWithUint8)
{
    uint8_t e = 2;
    BitEdges edges(e);
    EXPECT_EQ(e, edges);
}

TEST_F(BitEdgesTest, Assignment)
{
    BitEdges edges;
    edges = x;
    EXPECT_EQ(x, edges);
}

TEST_F(BitEdgesTest, AssignmentWithUint8)
{
    uint8_t e = 5;
    BitEdges edges;
    edges = e;
    EXPECT_EQ(e, edges);
}

TEST_F(BitEdgesTest, Add)
{
    BitEdges edges;
    EXPECT_EQ(false, edges[0]);
    edges.Add(0);
    EXPECT_EQ(true, edges[0]);
}

TEST_F(BitEdgesTest, Remove)
{
    BitEdges edges;
    EXPECT_EQ(false, edges[1]);
    edges.Add(1);
    EXPECT_EQ(true, edges[1]);
    edges.Remove(1);
    EXPECT_EQ(false, edges[1]);
}

TEST_F(BitEdgesTest, Swap)
{
    BitEdges e1 = x;
    BitEdges e2 = y;
    e1.swap(e2);

    EXPECT_EQ(x, e2);
    EXPECT_EQ(y, e1);
}

TEST_F(BitEdgesTest, Size)
{
    BitEdges edges;
    for (int i = 0; i < 4; ++i)
    {
        EXPECT_EQ(i, edges.size());
        edges.Add(i);
        EXPECT_EQ(i+1, edges.size());
    }
}

TEST_F(BitEdgesTest, Clear)
{
    x.clear();
    EXPECT_EQ(0, x.size());

    y.clear();
    EXPECT_EQ(0, y.size());
}

