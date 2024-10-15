//
// Created by andreas on 13.10.24.
//

#include "gtest/gtest.h"
#include "VEGAS_Integrator.h"
#include <cmath>
#include <numbers>

constexpr double sigma_range{3.0};

double sin_cos_tan(std::vector<double> x, void *param) {
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0] * pi;
    auto x1 = x[1];
    auto x2 = 2.0 * x[2] - 1.0;
    constexpr auto jacobi = 2.0 * pi;
    return jacobi * std::sin(1.0 - x0 * x0) * x1 * std::tan(x2 * x2) * std::tan(x2 * x2);
}

double sinus(std::vector<double> x, void *param) {
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];

    return std::sin(x0 * x1 * x2);
}

double polynom(std::vector<double> x, void *param) {
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];

    return (1.0 + x0 + 2.0 * x1 + x2 * x2) * (1.0 + x0 + 2.0 * x1 + x2 * x2) / (1.0 + x2) / (1.0 + x2) / (1.0 + x2);
}

double polynom2(std::vector<double> x, void *param) {
    constexpr auto pi = std::numbers::pi;
    auto x0 = x[0];
    auto x1 = x[1];
    auto x2 = x[2];
    auto x3 = x[3];
    auto x4 = x[4];
    return x3 * x4 * (1.0 + x0 + 2.0 * x1 + x2 * x2) * (1.0 + x0 + 2.0 * x1 + x2 * x2) / (1.0 + x2) / (1.0 + x2) /
           (1.0 + x2);
}

TEST(SimpleFunctionTest, polynom2) {
    constexpr double expected_result{0.7186547465398496};
    constexpr int dimensions{5};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_integrand(std::move(polynom2), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}


TEST(SimpleFunctionTest, polynom) {
    constexpr double expected_result{2.874618986159398};
    constexpr int dimensions{3};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_integrand(std::move(polynom), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}

TEST(SimpleFunctionTest, sinus) {
    constexpr double expected_result{0.12243402879673784};
    constexpr int dimensions{3};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_integrand(std::move(sinus), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}
