/**
 * @file histgram_unittest.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-13
 */

#include "gtest/gtest.h"

#include "basic/histgram.h"

#include <vector>

using namespace std;


class HistgramTest: public testing::Test
{
protected:
    virtual void SetUp()
    {
        for (int i = 0; i <= 100; ++i)
        {
            x.push_back(i);
            y.push_back(1.0 * rand() / RAND_MAX);
        }

        for (int i = 0; i <= 100; ++i)
        {
            histgram_x.insert(x[i]);
            histgram_y.insert(y[i]);
        }
    }

    vector<int> x;
    vector<double> y;
    Histgram<int> histgram_x;
    Histgram<double> histgram_y;
};

TEST_F(HistgramTest, DefaultConstructor)
{
    Histgram<int> hist;
    EXPECT_EQ(0U, hist.size());
    EXPECT_EQ(0, hist.minimum());
    EXPECT_EQ(0, hist.maximum());
}

TEST_F(HistgramTest, CopyConstructor)
{
    Histgram<int> hist(histgram_x);
    EXPECT_EQ(histgram_x.size(), hist.size());
    EXPECT_EQ(histgram_x.minimum(), hist.minimum());
    EXPECT_EQ(histgram_x.maximum(), hist.maximum());
}

TEST_F(HistgramTest, Assignment)
{
    Histgram<int> hist;
    hist = histgram_x;
    EXPECT_EQ(histgram_x.size(), hist.size());
    EXPECT_EQ(histgram_x.minimum(), hist.minimum());
    EXPECT_EQ(histgram_x.maximum(), hist.maximum());
}

TEST_F(HistgramTest, Insert)
{
    Histgram<int> hist;
    EXPECT_EQ(0, hist.count(1));
    hist.insert(1);
    EXPECT_EQ(1, hist.count(1));
    hist.insert(1);
    EXPECT_EQ(2, hist.count(1));
}

TEST_F(HistgramTest, Count)
{
    EXPECT_EQ(4-2, histgram_x.count(2, 4));
    EXPECT_EQ(100-3, histgram_x.count(3, 100));
}

TEST_F(HistgramTest, Minimum)
{
    EXPECT_EQ(*min_element(x.begin(), x.end()), histgram_x.minimum());
    EXPECT_EQ(*min_element(y.begin(), y.end()), histgram_y.minimum());
}

TEST_F(HistgramTest, Maximum)
{
    EXPECT_EQ(*max_element(x.begin(), x.end()), histgram_x.maximum());
    EXPECT_EQ(*max_element(y.begin(), y.end()), histgram_y.maximum());
}

TEST_F(HistgramTest, Mean)
{
    double sum = 0;
    for (unsigned i = 0; i < x.size(); ++i)
        sum += x[i];
    EXPECT_EQ(true, (sum / x.size()) - histgram_x.mean() < 1e-6);

    sum = 0;
    for (unsigned i = 0; i < y.size(); ++i)
        sum += y[i];
    EXPECT_EQ(true, (sum / y.size()) - histgram_y.mean() < 1e-6);
}

TEST_F(HistgramTest, Variance)
{
    double sum = 0;
    double square_sum = 0;
    for (unsigned i = 0; i < x.size(); ++i)
    {
        sum += x[i];
        square_sum += x[i] * x[i];
    }
    EXPECT_EQ(true, (square_sum/x.size() - (sum/x.size()) * (sum/x.size())) - histgram_x.variance() < 1e-6);

    sum = 0;
    square_sum = 0;
    for (unsigned i = 0; i < y.size(); ++i)
    {
        sum += y[i];
        square_sum += y[i] * y[i];
    }
    EXPECT_EQ(true, (square_sum/y.size() - (sum/y.size()) * (sum/y.size())) - histgram_y.variance() < 1e-6);
}

TEST_F(HistgramTest, Median)
{
    sort(x.begin(), x.end());
    sort(y.begin(), y.end());

    EXPECT_EQ(x[x.size()/2], histgram_x.median());
    EXPECT_EQ(y[y.size()/2], histgram_y.median());
}

TEST_F(HistgramTest, Trim)
{
    histgram_x.Trim(0.1);
    EXPECT_EQ(int(x.size() * 0.9 + 0.5), int(histgram_x.size()));
    histgram_y.Trim(0.1);
    EXPECT_EQ(int(y.size() * 0.9 + 0.5), int(histgram_y.size()));
}

TEST_F(HistgramTest, Swap)
{
    Histgram<int> histgram;
    swap(histgram, histgram_x);

    EXPECT_EQ(0U, histgram_x.size());
    EXPECT_EQ(x.size(), histgram.size());
}

TEST_F(HistgramTest, Empty)
{
    Histgram<int> hist;
    EXPECT_EQ(true, hist.empty());
    EXPECT_EQ(false, histgram_x.empty());
}

TEST_F(HistgramTest, Size)
{
    EXPECT_EQ(x.size(), histgram_x.size());
    EXPECT_EQ(y.size(), histgram_y.size());
}


