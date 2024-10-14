//
// Created by andreas on 13.10.24.
//

#include "gtest/gtest.h"
#include "VEGAS_Integrator.h"
#include <cmath>
#include <numbers>

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


TEST(SimpleFunctionTest, polynom) {
    constexpr double expected_result{2.874618986159398};
    constexpr int dimensions{3};
    VEGAS_Integrator integrator;
    integrator.Set_Integrand(std::move(polynom), dimensions, nullptr);
    integrator.Improve_Grid();
    integrator.Integration();
    std::cout << integrator.Get_Result() << " +/- " << integrator.Get_Error() << " with chi-square: "
              << integrator.Get_Chisq() << std::endl;

    EXPECT_NEAR(expected_result, integrator.Get_Result(), 2.0 * integrator.Get_Error());
}

TEST(SimpleFunctionTest, sinus) {
    constexpr double expected_result{0.12243402879673784};
    constexpr int dimensions{3};
    VEGAS_Integrator integrator;
    integrator.Set_Integrand(std::move(sinus), dimensions, nullptr);
    integrator.Improve_Grid();
    integrator.Integration();
    std::cout << integrator.Get_Result() << " +/- " << integrator.Get_Error() << " with chi-square: "
              << integrator.Get_Chisq() << std::endl;

    EXPECT_NEAR(expected_result, integrator.Get_Result(), 2.0 * integrator.Get_Error());
}
