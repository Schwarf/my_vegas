//
// Created by andreas on 20.10.24.
//

#ifndef ABS_VEGAS_TEST_INTEGRANDS_H
#define ABS_VEGAS_TEST_INTEGRANDS_H
#include <array>
#include <numbers>
#include <cmath>

constexpr double sigma_range{3.0};

// This seems a hard function
template<int dimension>
double sin_cos_tan(std::array<double, dimension> x, void *param) {
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0] * pi;
    auto x1 = x[1];
    auto x2 = 2.0 * x[2] - 1.0;
    constexpr auto jacobi = 2.0 * pi;
    return jacobi * std::sin(1.0 - x0 * x0) * std::cos(x1) * std::tan(x2 * x2);
}

template<int dimension>
double sinus(std::array<double, dimension> x, void *param) {
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];

    return std::sin(x0 * x1 * x2);
}

template<int dimension>
double sinus_10dim(std::array<double, dimension> x, void *param) {

    return std::sin(x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*x[8]*x[9]);
}

template<int dimension>
double polynom(std::array<double, dimension> x, void *param) {
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];

    return (1.0 + x0 + 2.0 * x1 + x2 * x2) * (1.0 + x0 + 2.0 * x1 + x2 * x2) / (1.0 + x2) / (1.0 + x2) / (1.0 + x2);
}

template<int dimension>
double polynom2(std::array<double, dimension> x, void *param) {
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];
    auto x3 = x[3];
    auto x4 = x[4];
    return x3 * x4 * (1.0 + x0 + 2.0 * x1 + x2 * x2) * (1.0 + x0 + 2.0 * x1 + x2 * x2) / (1.0 + x2) / (1.0 + x2) /
           (1.0 + x2);
}

template<int dimension>
double log_exp(std::array<double, dimension> x, void *param) {
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0];
    auto x1 = 2.0 * pi * x[1] - pi;
    auto x2 = 2.0 * x[2] - 1.0;
    constexpr auto jacobi = 2.0 * pi * 2.0;
    return jacobi * std::log(x0) * std::exp(x0) * x1 * x1 / (2.0 + x2 * x2 * x2);
}

#endif //ABS_VEGAS_TEST_INTEGRANDS_H
