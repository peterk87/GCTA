#include <gtest/gtest.h>

#include "StrFunc.h"

TEST(StrFunc, SplitWhitespace) {
    std::vector<std::string> parts;
    EXPECT_EQ(StrFunc::split_string("a b\tc", parts), 3);
    ASSERT_EQ(parts.size(), 3u);
    EXPECT_EQ(parts[0], "a");
    EXPECT_EQ(parts[1], "b");
    EXPECT_EQ(parts[2], "c");
}

TEST(StrFunc, SplitEmpty) {
    std::vector<std::string> parts;
    EXPECT_EQ(StrFunc::split_string("", parts), 0);
    EXPECT_TRUE(parts.empty());
}

TEST(StrFunc, SplitCustomSeparator) {
    std::vector<std::string> parts;
    EXPECT_EQ(StrFunc::split_string("rs1,rs2,rs3", parts, ","), 3);
    ASSERT_EQ(parts.size(), 3u);
    EXPECT_EQ(parts[0], "rs1");
    EXPECT_EQ(parts[2], "rs3");
}
