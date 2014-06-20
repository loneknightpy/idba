/**
 * @file vertex_status_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-07
 */

#include "gtest/gtest.h"

#include "graph/vertex_status.h"

#include <algorithm>

using namespace std;


class VertexStatusTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        s1.SetUsedFlag();
        s2.SetDeadFlag();
        s1.LockPreempt(1);
        s2.LockPreempt(2);
    }

    VertexStatus s1;
    VertexStatus s2;
};

TEST_F(VertexStatusTest, DefaultConstructor)
{
    VertexStatus s;
    EXPECT_EQ(false, s.IsUsed());
    EXPECT_EQ(false, s.IsDead());
    EXPECT_EQ(-1, s.GetLockID());
}

TEST_F(VertexStatusTest, CopyConstructor)
{
    VertexStatus s(s1);
    EXPECT_EQ(s1.IsUsed(), s.IsUsed());
    EXPECT_EQ(s1.IsDead(), s.IsDead());
    EXPECT_EQ(s1.GetLockID(), s.GetLockID());
}

TEST_F(VertexStatusTest, Assignment)
{
    VertexStatus s;
    s = s1;
    EXPECT_EQ(s1.IsUsed(), s.IsUsed());
    EXPECT_EQ(s1.IsDead(), s.IsDead());
    EXPECT_EQ(s1.GetLockID(), s.GetLockID());
}

TEST_F(VertexStatusTest, UsedFlag)
{
    VertexStatus s;
    EXPECT_EQ(false, s.IsUsed());
    s.SetUsedFlag();
    EXPECT_EQ(true, s.IsUsed());
    s.SetUsedFlag();
    EXPECT_EQ(true, s.IsUsed());
    s.ResetUsedFlag();
    EXPECT_EQ(false, s.IsUsed());
    s.ResetUsedFlag();
    EXPECT_EQ(false, s.IsUsed());
}

TEST_F(VertexStatusTest, DeadFlag)
{
    VertexStatus s;
    EXPECT_EQ(false, s.IsDead());
    s.SetDeadFlag();
    EXPECT_EQ(true, s.IsDead());
    s.SetDeadFlag();
    EXPECT_EQ(true, s.IsDead());
    s.ResetDeadFlag();
    EXPECT_EQ(false, s.IsDead());
    s.ResetDeadFlag();
    EXPECT_EQ(false, s.IsDead());
}

TEST_F(VertexStatusTest, LockPreempt)
{
    VertexStatus s;

    EXPECT_EQ(true, s.LockPreempt(0));

    EXPECT_EQ(false, s1.LockPreempt(0));
    EXPECT_EQ(false, s1.LockPreempt(1));
    EXPECT_EQ(true, s1.LockPreempt(2));

    EXPECT_EQ(false, s2.LockPreempt(0));
    EXPECT_EQ(false, s2.LockPreempt(1));
    EXPECT_EQ(false, s2.LockPreempt(2));
    EXPECT_EQ(true, s2.LockPreempt(3));

}

TEST_F(VertexStatusTest, Swap)
{
    VertexStatus x = s1;
    VertexStatus y = s2;
    swap(x, y);
    EXPECT_EQ(s2.IsDead(), x.IsDead());
    EXPECT_EQ(s2.IsUsed(), x.IsUsed());
    EXPECT_EQ(s2.GetLockID(), x.GetLockID());
    EXPECT_EQ(s1.IsDead(), y.IsDead());
    EXPECT_EQ(s1.IsUsed(), y.IsUsed());
    EXPECT_EQ(s1.GetLockID(), y.GetLockID());
}


