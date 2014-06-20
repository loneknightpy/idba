/**
 * @file sequence_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#include "gtest/gtest.h"

#include "sequence/sequence.h"

#include <sstream>
#include <string>

#include "basic/kmer.h"
#include "sequence/compact_sequence.h"
#include "sequence/short_sequence.h"


using namespace std;

class SequenceTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        simple_sequence.Assign("ACgttgca");
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
    string simple_string;
    string encoded_string;

    Kmer kmer0123;
    Kmer kmer3210;
};

TEST_F(SequenceTest, StreamIn)
{
    Sequence seq;
    stringstream ss(simple_string);
    ss >> seq;
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, StreamOut)
{
    stringstream ss;
    ss << simple_sequence;
    EXPECT_STRCASEEQ(simple_string.c_str(), ss.str().c_str());
}

TEST_F(SequenceTest, DefaultConstructor)
{
    Sequence seq;
    EXPECT_EQ(0, int(seq.size()));
}

TEST_F(SequenceTest, CopyConstructor)
{
    Sequence seq(simple_sequence);
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, ConstructorWithString)
{
    Sequence seq(simple_string);
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, ConstructorWithNumAndChar)
{
    Sequence seq(10, 3);
    EXPECT_EQ(10, int(seq.size()));
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(3, int(seq[i]));
}

TEST_F(SequenceTest, ConstructorWithCompactString)
{
    CompactSequence compact_seq(simple_sequence);
    Sequence seq(compact_seq);
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, ConstructorWithShortSequence)
{
    ShortSequence short_seq(simple_sequence);
    Sequence seq(short_seq);
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, ConstructorWithKmer)
{
    Sequence seq(kmer0123);
    EXPECT_EQ(kmer0123.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(kmer0123[i], seq[i]);
}

TEST_F(SequenceTest, Index)
{
    EXPECT_EQ(encoded_string.size(), simple_sequence.size());
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
        EXPECT_EQ(encoded_string[i], simple_sequence[i]);

    simple_sequence.set_base(0, 3);
    EXPECT_EQ(3, simple_sequence.get_base(0));
}

TEST_F(SequenceTest, Assignment)
{
    Sequence seq;
    seq = simple_sequence;
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, AssignmentWithString)
{
    Sequence seq;
    seq = simple_string;
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, AssignmentWithNumAndChar)
{
    Sequence seq;
    seq.Assign(10, 3);
    EXPECT_EQ(10, int(seq.size()));
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(3, int(seq[i]));
}

TEST_F(SequenceTest, AssignmentWithCompactSequence)
{
    CompactSequence compact_seq(simple_sequence);
    Sequence seq;
    seq = compact_seq;
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, AssignmentWithShortSequence)
{
    ShortSequence short_seq(simple_sequence);
    Sequence seq;
    seq = short_seq;
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, AssignmentWithKmer)
{
    Sequence seq;
    seq = kmer0123;
    EXPECT_EQ(kmer0123.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(kmer0123[i], seq[i]);
}

TEST_F(SequenceTest, AddAssignment)
{
    Sequence seq = simple_sequence;
    seq += simple_sequence;
    EXPECT_EQ(simple_sequence.size() * 2, seq.size());
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
    {
        EXPECT_EQ(simple_sequence[i], seq[i]);
        EXPECT_EQ(simple_sequence[i], seq[i + simple_sequence.size()]);
    }
}

TEST_F(SequenceTest, AddAssignmentWithChar)
{
    uint32_t old_size = simple_sequence.size();
    simple_sequence += 3;
    EXPECT_EQ(old_size + 1, simple_sequence.size());
    EXPECT_EQ(3, simple_sequence[simple_sequence.size()-1]);
}

TEST_F(SequenceTest, Assign)
{
    Sequence seq;
    seq.Assign(simple_sequence, 2, 4);
    EXPECT_EQ(4, int(seq.size()));
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[2+i], seq[i]);
}

TEST_F(SequenceTest, AssignWithString)
{
    Sequence seq;
    seq.Assign(simple_string, 2, 4);
    EXPECT_EQ(4, int(seq.size()));
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[2+i], seq[i]);
}

TEST_F(SequenceTest, AssignWithCompactSequence)
{
    CompactSequence compact_seq(simple_sequence);
    Sequence seq;
    seq.Assign(compact_seq, 2, 4);
    EXPECT_EQ(4, int(seq.size()));
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[2+i], seq[i]);
}

TEST_F(SequenceTest, AssignWithShortSequence)
{
    ShortSequence short_seq(simple_sequence);
    Sequence seq;
    seq.Assign(short_seq);
    EXPECT_EQ(simple_sequence.size(), seq.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq[i]);
}

TEST_F(SequenceTest, Compare)
{
    Sequence seq = simple_sequence;
    EXPECT_EQ(1, int(simple_sequence == seq));
    EXPECT_EQ(0, int(simple_sequence != seq));
    EXPECT_EQ(0, int(simple_sequence < seq));
    EXPECT_EQ(0, int(simple_sequence > seq));

    seq[0] = 1;
    EXPECT_EQ(0, int(simple_sequence == seq));
    EXPECT_EQ(1, int(simple_sequence != seq));

    EXPECT_EQ(1, int(simple_sequence < seq));
    EXPECT_EQ(0, int(simple_sequence > seq));
}

TEST_F(SequenceTest, IsValid)
{
    EXPECT_EQ(1, int(simple_sequence.IsValid()));
    simple_sequence[0] = 4;
    EXPECT_EQ(0, int(simple_sequence.IsValid()));
}

TEST_F(SequenceTest, ReverseComplement)
{
    Sequence seq = simple_sequence;
    seq[0] = 1;
    Sequence rev_comp = seq;
    rev_comp.ReverseComplement();

    EXPECT_EQ(seq.size(), rev_comp.size());
    for (unsigned i = 0; i < seq.size(); ++i)
        EXPECT_EQ(3, seq[i] + rev_comp[seq.size()-1-i]);
}

TEST_F(SequenceTest, GetKmer)
{
    Kmer kmer = simple_sequence.GetKmer(0, 4);
    EXPECT_EQ(1, int(kmer0123 == kmer));

    kmer = simple_sequence.GetKmer(4, 4);
    EXPECT_EQ(1, int(kmer3210 == kmer));
}

TEST_F(SequenceTest, Swap)
{
    Sequence seq1;
    Sequence seq2;
    seq1 = simple_sequence;
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq1[i]);

    swap(seq1, seq2);
    EXPECT_EQ(0U, seq1.size());
    EXPECT_EQ(simple_sequence.size(), seq2.size());
    for (unsigned i = 0; i < simple_sequence.size(); ++i)
        EXPECT_EQ(simple_sequence[i], seq2[i]);
}

TEST_F(SequenceTest, Size)
{
    EXPECT_EQ(simple_string.size(), simple_sequence.size());
}

TEST_F(SequenceTest, Resize)
{
    simple_sequence.resize(5);
    EXPECT_EQ(5, int(simple_sequence.size()));
}

TEST_F(SequenceTest, Empty)
{
    EXPECT_EQ(false, simple_sequence.empty());
    Sequence seq;
    EXPECT_EQ(true, seq.empty());
}


