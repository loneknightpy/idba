/**
 * @file atomic_integer_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-07
 */

#include "gtest/gtest.h"

#include "basic/atomic_integer.h"

#include <algorithm>

using namespace std;


class AtomicIntegerTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        a = 17;
        b = 22;
    }

    AtomicInteger<int> a;
    AtomicInteger<int> b;
};

TEST_F(AtomicIntegerTest, DefaultConstructor)
{
    AtomicInteger<int> x;
    EXPECT_EQ(0, x);
}

TEST_F(AtomicIntegerTest, ContructorWithValue)
{
    AtomicInteger<int> x(10);
    EXPECT_EQ(10, x);
}

TEST_F(AtomicIntegerTest, CopyConstructor)
{
    AtomicInteger<int> x(a);
    EXPECT_EQ(a, x);
}

TEST_F(AtomicIntegerTest, Assignment)
{
    AtomicInteger<int> x;
    x = a;
    EXPECT_EQ(a, x);

    int z = a + b;

    a = b = (a + b);
    EXPECT_EQ(z, a);
    EXPECT_EQ(z, b);
}

TEST_F(AtomicIntegerTest, AssignmentWithValue)
{
    AtomicInteger<int> x;
    x = 1;
    EXPECT_EQ(1, x);
}

TEST_F(AtomicIntegerTest, Greater)
{
    EXPECT_EQ(0, int(a > b));
    EXPECT_EQ(1, int(b > a));
}

TEST_F(AtomicIntegerTest, Less)
{
    EXPECT_EQ(1, int(a < b));
    EXPECT_EQ(0, int(b < a));
}

TEST_F(AtomicIntegerTest, Equal)
{
    AtomicInteger<int> x(a);
    EXPECT_EQ(0, int(a == b));
    EXPECT_EQ(1, int(a == x));
}

TEST_F(AtomicIntegerTest, NotEqual)
{
    AtomicInteger<int> x(a);
    EXPECT_EQ(1, int(a != b));
    EXPECT_EQ(0, int(a != x));
}

TEST_F(AtomicIntegerTest, AddAssignment)
{
    AtomicInteger<int> x = a;
    x += b;
    EXPECT_EQ(a + b, x);
}

TEST_F(AtomicIntegerTest, SubAssignment)
{
    AtomicInteger<int> x = a;
    x -= b;
    EXPECT_EQ(a - b, x);
}

TEST_F(AtomicIntegerTest, OrAssignment)
{
    AtomicInteger<int> x = a;
    x |= b;
    EXPECT_EQ(a | b, x);
}

TEST_F(AtomicIntegerTest, AndAssignment)
{
    AtomicInteger<int> x = a;
    x &= b;
    EXPECT_EQ(a & b, x);
}

TEST_F(AtomicIntegerTest, XorAssignment)
{
    AtomicInteger<int> x = a;
    x ^= b;
    EXPECT_EQ(a ^ b, x);
}

TEST_F(AtomicIntegerTest, Increment)
{
    AtomicInteger<int> x = a;
    EXPECT_EQ(a, x++);
    x = a;
    EXPECT_EQ(a+1, ++x);
}

TEST_F(AtomicIntegerTest, Decrement)
{
    AtomicInteger<int> x = a;
    EXPECT_EQ(a, x--);
    x = a;
    EXPECT_EQ(a-1, --x);
}

TEST_F(AtomicIntegerTest, CompareAndSet)
{
    int old_a = a;
    int new_a = a + 10;
    EXPECT_EQ(true, a.CompareAndSet(old_a, new_a));

    old_a = a;
    new_a = a + 10;
    a += 5;
    EXPECT_EQ(false, a.CompareAndSet(old_a, new_a));
}

TEST_F(AtomicIntegerTest, Swap)
{
    AtomicInteger<int> x = a;
    AtomicInteger<int> y = b;
    swap(x, y);
    EXPECT_EQ(b, x);
    EXPECT_EQ(a, y);
}



