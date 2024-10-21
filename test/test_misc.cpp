//
// Created by andreas on 16.10.24.
//
#include "gtest/gtest.h"
#include "misc.h"
#include <cmath>


template<int X, int N>
struct TestNthRoot {
    static constexpr int computed_root = nth_root(X, N);
};

TEST(ComputeNthRoot, simple1)
{
    constexpr int x = 1000; // value to find the root of
    constexpr int n = 2;       // root degree
    constexpr auto root = nth_root(x, n);
    EXPECT_EQ(root, std::floor(std::pow(x, 1.0/n)));
}
