/**
 * @file compact_sequence_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-03
 */

#include "gtest/gtest.h"

#include "sequence/compact_sequence.h"

#include <string>

#include "basic/kmer.h"
#include "sequence/sequence.h"
#include "sequence/short_sequence.h"


using namespace std;

class CompactSequenceTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        simple_sequence.Assign("ACgttgca");
        simple_compact_sequence.Assign(simple_sequence);
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

    Sequence simple_sequence;
    CompactSequence simple_compact_sequence;
    string simple_string;
    string encoded_string;

    Kmer kmer0123;
    Kmer kmer3210;
};

TEST_F(CompactSequenceTest, StreamIn)
{
    CompactSequence compact_seq;
    stringstream ss(simple_string);
    ss >> compact_seq;
    EXPECT_EQ(simple_sequence.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, StreamOut)
{
    stringstream ss;
    ss << simple_compact_sequence;
    EXPECT_STRCASEEQ(simple_string.c_str(), ss.str().c_str());
}

TEST_F(CompactSequenceTest, DefaultConstructor)
{
    CompactSequence compact_seq;
    EXPECT_EQ(0, int(compact_seq.size()));
}

TEST_F(CompactSequenceTest, CopyConstructor)
{
    CompactSequence compact_seq(simple_compact_sequence);
    EXPECT_EQ(simple_compact_sequence.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_compact_sequence[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, ConstructorWithSequence)
{
    CompactSequence compact_seq(simple_sequence);
    EXPECT_EQ(simple_sequence.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, ConstructorWithShortSequence)
{
    ShortSequence short_seq(simple_sequence);
    CompactSequence compact_seq(short_seq);
    EXPECT_EQ(simple_sequence.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, ConstructorWithKmer)
{
    CompactSequence compact_seq(kmer0123);
    EXPECT_EQ(kmer0123.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(kmer0123[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, Index)
{
    EXPECT_EQ(encoded_string.size(), simple_compact_sequence.size());
    for (unsigned i = 0; i < simple_compact_sequence.size(); ++i)
        EXPECT_EQ(encoded_string[i], simple_compact_sequence[i]);

    simple_compact_sequence.set_base(0, 3);
    EXPECT_EQ(3, simple_compact_sequence[0]);
}

TEST_F(CompactSequenceTest, Assignment)
{
    CompactSequence compact_seq;
    compact_seq = simple_compact_sequence;
    EXPECT_EQ(simple_compact_sequence.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_compact_sequence[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, AssignmentWithSequence)
{
    CompactSequence compact_seq;
    compact_seq = simple_sequence;
    EXPECT_EQ(simple_compact_sequence.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_compact_sequence[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, AssignmentWithShortSequence)
{
    ShortSequence short_seq(simple_sequence);
    CompactSequence compact_seq;
    compact_seq = short_seq;
    EXPECT_EQ(simple_compact_sequence.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_compact_sequence[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, AssignmentWithKmer)
{
    CompactSequence compact_seq;
    compact_seq = kmer0123;
    EXPECT_EQ(kmer0123.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(kmer0123[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, AddAssignment)
{
    CompactSequence compact_seq;
    compact_seq += simple_compact_sequence;
    compact_seq += simple_compact_sequence;
    EXPECT_EQ(simple_compact_sequence.size() * 2, compact_seq.size());
    for (unsigned i = 0; i < simple_compact_sequence.size(); ++i)
    {
        EXPECT_EQ(simple_compact_sequence[i], compact_seq[i]);
        EXPECT_EQ(simple_compact_sequence[i], compact_seq[i + simple_compact_sequence.size()]);
    }
}

TEST_F(CompactSequenceTest, AddAssignmentWithSequence)
{
    CompactSequence compact_seq;
    compact_seq += simple_sequence;
    compact_seq += simple_sequence;
    EXPECT_EQ(simple_sequence.size() * 2, compact_seq.size());
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
    {
        EXPECT_EQ(simple_sequence[i], compact_seq[i]);
        EXPECT_EQ(simple_sequence[i], compact_seq[i + simple_compact_sequence.size()]);
    }
}

TEST_F(CompactSequenceTest, AddAssignmentWithChar)
{
    uint32_t old_size = simple_compact_sequence.size();
    simple_compact_sequence += 3;
    EXPECT_EQ(old_size + 1, simple_compact_sequence.size());
    EXPECT_EQ(3, simple_compact_sequence[simple_compact_sequence.size()-1]);
}

TEST_F(CompactSequenceTest, Assign)
{
    CompactSequence compact_seq;
    compact_seq.Assign(simple_compact_sequence, 2, 4);
    EXPECT_EQ(4, int(compact_seq.size()));
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[2+i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, AssignWithSequence)
{
    CompactSequence compact_seq;
    compact_seq.Assign(simple_sequence, 2, 4);
    EXPECT_EQ(4, int(compact_seq.size()));
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[2+i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, AssignWithShortSequence)
{
    ShortSequence short_seq(simple_sequence);
    CompactSequence compact_seq;
    compact_seq.Assign(short_seq);
    EXPECT_EQ(simple_sequence.size(), compact_seq.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], compact_seq[i]);
}

TEST_F(CompactSequenceTest, Compare)
{
    CompactSequence compact_seq = simple_compact_sequence;
    EXPECT_EQ(1, int(simple_compact_sequence == compact_seq));
    EXPECT_EQ(0, int(simple_compact_sequence != compact_seq));
    EXPECT_EQ(0, int(simple_compact_sequence < compact_seq));
    EXPECT_EQ(0, int(simple_compact_sequence > compact_seq));

    compact_seq.set_base(0, 1);
    EXPECT_EQ(0, int(simple_compact_sequence == compact_seq));
    EXPECT_EQ(1, int(simple_compact_sequence != compact_seq));

    EXPECT_EQ(1, int(simple_compact_sequence < compact_seq));
    EXPECT_EQ(0, int(simple_compact_sequence > compact_seq));
}

TEST_F(CompactSequenceTest, ReverseComplement)
{
    CompactSequence compact_seq = simple_compact_sequence;
    compact_seq.set_base(0, 1);
    CompactSequence rev_comp = compact_seq;
    rev_comp.ReverseComplement();

    EXPECT_EQ(compact_seq.size(), rev_comp.size());
    for (unsigned i = 0; i < compact_seq.size(); ++i)
        EXPECT_EQ(3, compact_seq[i] + rev_comp[compact_seq.size()-1-i]);
}

TEST_F(CompactSequenceTest, GetKmer)
{
    Kmer kmer = simple_compact_sequence.GetKmer(0, 4);
    EXPECT_EQ(1, int(kmer0123 == kmer));

    kmer = simple_compact_sequence.GetKmer(4, 4);
    EXPECT_EQ(1, int(kmer3210 == kmer));
}

TEST_F(CompactSequenceTest, Swap)
{
    CompactSequence compact_seq1;
    CompactSequence compact_seq2;
    compact_seq1 = simple_sequence;
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
        EXPECT_EQ(simple_sequence[i], compact_seq1[i]);

    swap(compact_seq1, compact_seq2);
    EXPECT_EQ(0U, compact_seq1.size());
    EXPECT_EQ(simple_sequence.size(), compact_seq2.size());
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
        EXPECT_EQ(simple_sequence[i], compact_seq2[i]);
}

TEST_F(CompactSequenceTest, Size)
{
    EXPECT_EQ(simple_string.size(), simple_compact_sequence.size());
}

TEST_F(CompactSequenceTest, Resize)
{
    simple_compact_sequence.resize(5);
    EXPECT_EQ(5, int(simple_compact_sequence.size()));
}

TEST_F(CompactSequenceTest, Empty)
{
    EXPECT_EQ(false, simple_compact_sequence.empty());
    CompactSequence compact_seq;
    EXPECT_EQ(true, compact_seq.empty());
}


