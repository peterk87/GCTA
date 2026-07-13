#include <gtest/gtest.h>

#include "CommFunc.h"

TEST(CommFunc, MedianOdd) {
    EXPECT_DOUBLE_EQ(CommFunc::median(std::vector<double>{3, 1, 2}), 2.0);
}

TEST(CommFunc, MedianEven) {
    EXPECT_DOUBLE_EQ(CommFunc::median(std::vector<double>{4, 1, 3, 2}), 2.5);
}

TEST(CommFunc, MedianSingle) {
    EXPECT_DOUBLE_EQ(CommFunc::median(std::vector<double>{7.5}), 7.5);
}

TEST(CommFunc, FloatEqual) {
    EXPECT_TRUE(CommFunc::FloatEqual(1.0, 1.0));
    EXPECT_TRUE(CommFunc::FloatEqual(1.0, 1.0 + CommFunc::FloatErr / 2.0));
    EXPECT_FALSE(CommFunc::FloatEqual(1.0, 1.1));
}

TEST(CommFunc, MeanAndSum) {
    std::vector<double> x{1.0, 2.0, 3.0};
    EXPECT_DOUBLE_EQ(CommFunc::sum(x), 6.0);
    EXPECT_DOUBLE_EQ(CommFunc::mean(x), 2.0);
}
