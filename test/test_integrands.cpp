//
// Created by andreas on 13.10.24.
//

#include "gtest/gtest.h"
#include "VEGAS_Integrator.h"
#include "integrands.h"

constexpr double sigma_range{3.0};


TEST(SimpleFunctionTest, almost_singular) {
    constexpr double expected_result{1.762747174039086};
    constexpr int dimensions{2};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(almost_singular<dimensions>), nullptr);
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
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(polynom1<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;
    const auto result =  integrator.get_result();
    EXPECT_NEAR(expected_result, result , sigma_range * integrator.get_error());
}



TEST(SimpleFunctionTest, sinus_10dim) {
    constexpr double expected_result{0.00097640369191418530169};
    constexpr int dimensions{10};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(sinus_10dim<dimensions>), nullptr);
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
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(log_exp<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}

TEST(SimpleFunctionTest, polynom2) {
    constexpr double expected_result{1.466821492328336};
    constexpr int dimensions{5};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(polynom2<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}



TEST(SimpleFunctionTest, sinus_3dim) {
    constexpr double expected_result{0.12243402879673784};
    constexpr int dimensions{3};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(sinus_3dim<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}


TEST(SimpleFunctionTest, dilogarithm) {
    constexpr double expected_result{std::numbers::pi*std::numbers::pi/6.0};
    constexpr int dimensions{1};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(dilogarithm<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chisquare() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}


//TEST(SimpleFunctionTest, difficult_for_vegas) {
//    constexpr double expected_result{6.679915983895715};
//    constexpr int dimensions{2};
//    VegasNumericalIntegration<dimensions> integrator;
//    integrator.set_verbosity(VegasVerbosity::Info);
//    integrator.set_integrand(std::move(difficult_for_vegas<dimensions>), nullptr);
//    integrator.improve_grid();
//    integrator.integrate();
//    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
//              << integrator.get_chisquare() << std::endl;
//
//    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
//}
