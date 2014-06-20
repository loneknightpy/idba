/**
 * @file kmer_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#include "gtest/gtest.h"

#include "basic/kmer.h"


using namespace std;

class KmerTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        kmer0123.resize(4);
        for (int i = 0; i < 4; ++i)
            kmer0123.ShiftAppend(i);

        kmer3210.resize(4);
        for (int i = 0; i < 4; ++i)
            kmer3210.ShiftPreappend(i);
    }

    Kmer kmer0123;
    Kmer kmer3210;
};

TEST_F(KmerTest, DefaultConstructor)
{
    Kmer kmer;
    EXPECT_EQ(0, int(kmer.size()));
}

TEST_F(KmerTest, CopyConstructor)
{
    Kmer kmer(kmer0123);
    EXPECT_EQ(1, int(kmer == kmer0123));
}

TEST_F(KmerTest, ConstructorWithSize)
{
    Kmer kmer(10);
    EXPECT_EQ(10, int(kmer.size()));

    for (int i = 0; i < 10; ++i)
        EXPECT_EQ(0, kmer[i]);
}

TEST_F(KmerTest, Assignment)
{
    Kmer kmer;
    kmer = kmer0123;
    EXPECT_EQ(kmer.size(), kmer0123.size());
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(kmer[i], kmer0123[i]);
}

TEST_F(KmerTest, Less)
{
    EXPECT_EQ(0, int(kmer0123 < kmer3210));
    Kmer kmer = kmer0123;
    EXPECT_EQ(0, int(kmer < kmer0123));
}

TEST_F(KmerTest, Greater)
{
    EXPECT_EQ(1, int(kmer0123 > kmer3210));
    Kmer kmer = kmer0123;
    EXPECT_EQ(0, int(kmer > kmer0123));
}

TEST_F(KmerTest, Equal)
{
    EXPECT_EQ(0, int(kmer0123 == kmer3210));
    Kmer kmer = kmer0123;
    EXPECT_EQ(1, int(kmer == kmer0123));
}

TEST_F(KmerTest, NotEqual)
{
    EXPECT_EQ(1, int(kmer0123 != kmer3210));
    Kmer kmer = kmer0123;
    EXPECT_EQ(0, int(kmer != kmer0123));
}

TEST_F(KmerTest, ReverseComplement)
{
    Kmer kmer = kmer0123;
    EXPECT_EQ(1, int(kmer.ReverseComplement() == kmer0123));

    kmer.set_base(0, 3);
    Kmer rev_comp = kmer;
    rev_comp.ReverseComplement();
    
    EXPECT_EQ(0, int(rev_comp == kmer));
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(3, kmer[i] + rev_comp[3-i]);
}

TEST_F(KmerTest, ShiftAppend)
{
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(i, kmer0123[i]);

    Kmer kmer = kmer3210;
    for (int i = 0; i < 4; ++i)
        kmer.ShiftAppend(i);

    EXPECT_EQ(1, int(kmer == kmer0123));
}

TEST_F(KmerTest, ShiftPreappend)
{
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(3-i, kmer3210[i]);

    Kmer kmer = kmer0123;
    for (int i = 0; i < 4; ++i)
        kmer.ShiftPreappend(i);

    EXPECT_EQ(1, int(kmer == kmer3210));
}

TEST_F(KmerTest, set_base)
{
    kmer0123.set_base(1, 3);
    EXPECT_EQ(3, kmer0123.get_base(1));
}

TEST_F(KmerTest, IsPalindrome)
{
    EXPECT_EQ(1, int(kmer0123.IsPalindrome()));
    Kmer kmer = kmer0123;
    kmer.set_base(0, 1);
    EXPECT_EQ(0, int(kmer.IsPalindrome()));
}

TEST_F(KmerTest, Swap)
{
    Kmer kmer1 = kmer0123;
    Kmer kmer2 = kmer3210;
    swap(kmer1, kmer2);

    for (int i = 0; i < 4; ++i)
    {
        EXPECT_EQ(kmer0123[i], kmer2[i]);
        EXPECT_EQ(kmer3210[i], kmer1[i]);
    }
}

TEST_F(KmerTest, Resize)
{
    Kmer kmer;
    kmer.resize(10);
    EXPECT_EQ(10, int(kmer.size()));
}

TEST_F(KmerTest, Clear)
{
    kmer0123.clear();
    EXPECT_EQ(4, int(kmer0123.size()));
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(0, kmer0123[i]);
}

