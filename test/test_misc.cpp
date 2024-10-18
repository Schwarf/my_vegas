//
// Created by andreas on 16.10.24.
//
#include "gtest/gtest.h"
#include "misc.h"
#include <cmath>



template<int X, int N>
struct TestNthRoot {
    static void run() {
        constexpr int computed_root = nth_root(X, N);
        int expected_root = static_cast<int>(std::pow(X, 1.0 / N));
        EXPECT_EQ(computed_root, expected_root);
    }
};

TEST(ComputeNthRoot, simple1)
{
    long long unsigned x{1000};

    for (int i{0}; i < 15; ++i) {
        x *= 10;
        for (int n{2}; n < 16; ++n) {
            // Compute NumberOfStratifications
            auto NumberOfStratifications = std::floor(std::pow(x, 1.0/n));
            // Compute NumberOfHyperCubes
            auto NumberOfHyperCubes = std::pow(NumberOfStratifications, n);

            // Set the cout to print fixed float point numbers with no digits after the decimal
            std::cout << std::fixed << std::setprecision(0);

            // Output the templated struct with both NumberOfStratifications and NumberOfHyperCubes
            std::cout << "template<> struct ValueMapper<" << x << "," << n << "> {\n"
                      << "    static constexpr int NumberOfStratifications = " << NumberOfStratifications << ";\n"
                      << "    static constexpr long long NumberOfHyperCubes = " << NumberOfHyperCubes << ";\n"
                      << "};" << std::endl;
        }
    }
}
