/**
 * @file short_sequence_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-03
 */

#include "gtest/gtest.h"

#include "sequence/short_sequence.h"

#include <string>

#include "basic/kmer.h"
#include "sequence/sequence.h"
#include "sequence/compact_sequence.h"


using namespace std;

class ShortSequenceTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        simple_sequence.Assign("ACgttgca");
        simple_short_sequence = simple_sequence;
        simple_string = "ACgttgca";
        encoded_string = "01233210";

        for (unsigned i = 0; i < encoded_string.size(); ++i)
            encoded_string[i] -= '0';

        kmer0123.resize(4);
        for (int i = 0; i < 4; ++i)
            kmer0123.ShiftAppend(i);

        kmer3210.resize(4);
        for (int i = 0; i < 4; ++i)
            kmer3210.ShiftPreappend(i);
    }

    ShortSequence simple_short_sequence;
    Sequence simple_sequence;
    string simple_string;
    string encoded_string;

    Kmer kmer0123;
    Kmer kmer3210;
};

TEST_F(ShortSequenceTest, DefaultConstructor)
{
    ShortSequence short_seq;
    EXPECT_EQ(0, int(short_seq.size()));
}

TEST_F(ShortSequenceTest, CopyConstructor)
{
    ShortSequence short_seq(simple_short_sequence);
    EXPECT_EQ(simple_short_sequence.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_short_sequence[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, ConstructorWithSequence)
{
    ShortSequence short_seq(simple_sequence);
    EXPECT_EQ(simple_sequence.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, ConstructorWithCompactString)
{
    CompactSequence compact_seq(simple_sequence);
    ShortSequence short_seq(compact_seq);
    EXPECT_EQ(simple_sequence.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, ConstructorWithKmer)
{
    ShortSequence short_seq(kmer0123);
    EXPECT_EQ(kmer0123.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(kmer0123[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, Index)
{
    EXPECT_EQ(encoded_string.size(), simple_short_sequence.size());
    for (unsigned i = 0; i < simple_short_sequence.size(); ++i)
        EXPECT_EQ(encoded_string[i], simple_short_sequence[i]);

    simple_short_sequence.set_base(0, 3);
    EXPECT_EQ(3, simple_short_sequence[0]);
}

TEST_F(ShortSequenceTest, Assignment)
{
    ShortSequence short_seq;
    short_seq = simple_short_sequence;
    EXPECT_EQ(simple_short_sequence.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_short_sequence[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, AssignmentWithSequence)
{
    ShortSequence short_seq;
    short_seq = simple_sequence;
    EXPECT_EQ(simple_short_sequence.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_short_sequence[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, AssignmentWithCompactSequence)
{
    CompactSequence compact_seq(simple_sequence);
    ShortSequence short_seq;
    short_seq = compact_seq;
    EXPECT_EQ(simple_short_sequence.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_short_sequence[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, AssignmentWithKmer)
{
    ShortSequence short_seq;
    short_seq = kmer0123;
    EXPECT_EQ(kmer0123.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(kmer0123[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, Assign)
{
    ShortSequence short_seq;
    short_seq.Assign(simple_short_sequence);
    EXPECT_EQ(simple_short_sequence.size(), short_seq.size());
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_short_sequence[i], short_seq[i]);
}

TEST_F(ShortSequenceTest, AssignWithSequence)
{
    ShortSequence short_seq;
    short_seq.Assign(simple_sequence, 2, 4);
    EXPECT_EQ(4, int(short_seq.size()));
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[2+i], short_seq[i]);
}

TEST_F(ShortSequenceTest, AssignWithCompactSequence)
{
    CompactSequence compact_seq(simple_sequence);
    ShortSequence short_seq;
    short_seq.Assign(compact_seq, 2, 4);
    EXPECT_EQ(4, int(short_seq.size()));
    for (unsigned i = 0; i < short_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[2+i], short_seq[i]);
}

TEST_F(ShortSequenceTest, Swap)
{
    ShortSequence short_seq1;
    ShortSequence short_seq2;
    short_seq1 = simple_sequence;
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
        EXPECT_EQ(simple_sequence[i], short_seq1[i]);

    swap(short_seq1, short_seq2);
    EXPECT_EQ(0U, short_seq1.size());
    EXPECT_EQ(simple_sequence.size(), short_seq2.size());
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
        EXPECT_EQ(simple_sequence[i], short_seq2[i]);
}

TEST_F(ShortSequenceTest, Size)
{
    EXPECT_EQ(simple_string.size(), simple_short_sequence.size());
}

TEST_F(ShortSequenceTest, Resize)
{
    simple_short_sequence.resize(5);
    EXPECT_EQ(5, int(simple_short_sequence.size()));
}

TEST_F(ShortSequenceTest, Empty)
{
    EXPECT_EQ(false, simple_short_sequence.empty());
    ShortSequence short_seq;
    EXPECT_EQ(true, short_seq.empty());
}
