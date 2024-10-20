//
// Created by andreas on 13.10.24.
//

#include "gtest/gtest.h"
#include "VEGAS_Integrator.h"
#include "test_integrands.h"

TEST(SimpleFunctionTest, polynom2) {
    constexpr double expected_result{0.7186547465398496};
    constexpr int dimensions{5};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_integrand(std::move(polynom2<dimensions>), nullptr);
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
    integrator.set_integrand(std::move(polynom<dimensions>), nullptr);
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
    integrator.set_integrand(std::move(sinus<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}

TEST(SimpleFunctionTest, log_exp) {
    constexpr double expected_result{-28.37381254755662};
    constexpr int dimensions{3};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_integrand(std::move(log_exp<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}

TEST(SimpleFunctionTest, sinus_10dim) {
    constexpr double expected_result{0.00097640369191418530169};
    constexpr int dimensions{10};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_integrand(std::move(sinus_10dim<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}

