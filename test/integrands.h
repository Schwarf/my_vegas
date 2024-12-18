//
// Created by andreas on 20.10.24.
//

#ifndef ABS_VEGAS_INTEGRANDS_H
#define ABS_VEGAS_INTEGRANDS_H
#include <array>
#include <numbers>
#include <cmath>
#include <iostream>


template <int dimension>
double sinus_3dim(std::array<double, dimension> x, void* param)
{
    (void)param;
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];

    return std::sin(x0 * x1 * x2);
}

template <int dimension>
double sinus_10dim(std::array<double, dimension> x, void* param)
{
    (void)param;
    return std::sin(x[0] * x[1] * x[2] * x[3] * x[4] * x[5] * x[6] * x[7] * x[8] * x[9]);
}


template <int dimension>
double polynom1(std::array<double, dimension> x, void* param)
{
    (void)param;
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];

    return (1.0 + x0 + 2.0 * x1 + x2 * x2) * (1.0 + x0 + 2.0 * x1 + x2 * x2) / (1.0 + x2) / (1.0 + x2) / (1.0 + x2);
}

template <int dimension>
double polynom2(std::array<double, dimension> x, void* param)
{
    (void)param;
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];
    auto x3 = x[3];
    auto x4 = x[4];
    return (1.0 + x3) * (1.0 + x3) * x4 * x4 * (5.0 + x0 + 2.0 * x1 + x2 * x2) * (-1.0 + x0 + 2.0 * x1 + x2 * x2) / (1.0
            + x2) / (1.0 + x2) /
        (1.0 + x2);
}

template <int dimension>
double polynom3(std::array<double, dimension> x, void* param)
{
    (void)param;
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];
    auto x3 = x[3];
    return (1.0 + x3 * x3 + x3 * x1 + 6.0 * x2 * x1 * x0 - x2 * x2 * x0) / (2 + x1 * x2 - x0 * x3 * x3 * x3 * x3);
}

// This seems a hard function
template <int dimension>
double sin_cos_tan(std::array<double, dimension> x, void* param)
{
    (void)param;
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0] * pi;
    auto x1 = x[1];
    auto x2 = 2.0 * x[2] - 1.0;
    constexpr auto jacobi = 2.0 * pi;
    return jacobi * std::sin(1.0 - x0 * x0) * std::cos(x1) * std::tan(x2 * x2);
}


template <int dimension>
double log_exp(std::array<double, dimension> x, void* param)
{
    (void)param;
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0];
    auto x1 = 2.0 * pi * x[1] - pi;
    auto x2 = 2.0 * x[2] - 1.0;
    constexpr auto jacobi = 2.0 * pi * 2.0;
    return jacobi * std::log(x0) * std::exp(x0) * x1 * x1 / (2.0 + x2 * x2 * x2);
}

template <int dimension>
double almost_singular(std::array<double, dimension> x, void* param)
{
    (void)param;
    constexpr double epsilon = 1.E-16;
    if (x[0] < epsilon && x[1] < epsilon)
        return 0.0;
    return 1.0 / std::sqrt(x[0] * x[0] + x[1] * x[1]);
}

template <int dimension>
double difficult_for_vegas(std::array<double, dimension> x, void* param)
{
    (void)param;
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0] * 2.0 * pi;
    auto x1 = x[1] * 2.0 * pi;
    constexpr auto jacobi = 4.0 * pi * pi;
    return std::sin(x0 * x1) / x0 * x1 * jacobi;
}

template <int dimension>
double dilogarithm(std::array<double, dimension> x, void* param)
{
    (void)param;
    return -std::log(1.0 - x[0]) / x[0];
}

template <int dimension>
double some_logs(std::array<double, dimension> x, void* param)
{
    (void)param;
    return std::log(x[0]) * std::log(x[0]) * std::log(1.0 - x[0]);
}


template <int dimension>
double gaussian_4d(std::array<double, dimension> x, void* param)
{
    (void)param;
    // Integrate from [-1,1] for x0,x1,x2,x3
    auto x0 = (2.0 * x[0] - 1.0);
    auto x1 = (2.0 * x[1] - 1.0);
    auto x2 = (2.0 * x[2] - 1.0);
    auto x3 = (2.0 * x[3] - 1.0);

    constexpr auto jacobi = 16.0;
    return jacobi * std::exp(-x0 * x0 - x1 * x1 - x2 * x2 - x3 * x3);
}

template <int dimension>
double volume_unit_sphere_3d(std::array<double, dimension> x, void* param)
{
    (void)param;
    // Integrate from [-1,1] for x0,x1,x2
    double value{};
    double jacobi{1};
    for(int dim{}; dim < dimension; ++dim)
    {
        x[dim] = 2.0 * x[dim] - 1.0;
        value += x[dim] * x[dim];
        jacobi *= 2.0;
    }
    return value < 1.0 ? jacobi : 0.0;
}

template <int dimension>
double surface_unit_sphere_3d(std::array<double, dimension> x, void* param)
{
    (void)param;
    // Integral ranges from [-1,1] for all dimensions
    constexpr auto pi = std::numbers::pi;
    auto x0 = pi * x[0];
    auto x1 = 2.0 * pi * x[1] ;
    auto x2 = x[2]; // radius
    constexpr auto jacobi = 2.0*pi*pi;
    return jacobi * std::sin(x0);
}

template <int dimension>
double arctan_derivative(std::array<double, dimension> x, void* param)
{
    (void)param;
    // Integrate from [-1,1] for x0,x1,x2
    return 1.0/(1.0 + x[0]*x[0])/(1.0 + x[1]*x[1])/(1.0 + x[2]*x[2]);
}

#endif //ABS_VEGAS_INTEGRANDS_H
