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
              << integrator.get_chi_square() << std::endl;

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
              << integrator.get_chi_square() << std::endl;
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
              << integrator.get_chi_square() << std::endl;

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
              << integrator.get_chi_square() << std::endl;

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
              << integrator.get_chi_square() << std::endl;

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
              << integrator.get_chi_square() << std::endl;

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
              << integrator.get_chi_square() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}


// TEST(SimpleFunctionTest, difficult_for_vegas) {
//     constexpr double expected_result{30.98141031607467};
//     constexpr int dimensions{2};
//     VegasNumericalIntegration<dimensions> integrator;
//     integrator.set_verbosity(VegasVerbosity::Info);
//     integrator.set_integrand(std::move(difficult_for_vegas<dimensions>), nullptr);
//     integrator.improve_grid();
//     integrator.integrate();
//     std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
//               << integrator.get_chisquare() << std::endl;
//
//     EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
// }

TEST(SimpleFunctionTest, some_logs) {
    constexpr double expected_result{-0.3060180599843586};
    constexpr int dimensions{1};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(some_logs<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chi_square() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}


TEST(SimpleFunctionTest, gaussian_4d) {
    constexpr double expected_result{4.977294701166026};
    constexpr int dimensions{4};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(gaussian_4d<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chi_square() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}

double expected_volume_unit_sphere(int dimensions) {
    return std::pow(std::numbers::pi, dimensions / 2.0) / std::tgamma(dimensions / 2.0 + 1);
}

template <typename Dimension>
class VolumeUnitSphereTest : public ::testing::Test {
public:
    static constexpr int dimensions = Dimension::value;
    double expected_result = expected_volume_unit_sphere(dimensions);
};

template <int Dimension>
struct DimensionWrapper {
    static constexpr int value = Dimension;
};

// Checked up to 8 dimensions, but it takes a few minutes
using DimensionsList = ::testing::Types<
    DimensionWrapper<1>,
    DimensionWrapper<2>,
    DimensionWrapper<3>,
    DimensionWrapper<4>,
    DimensionWrapper<5>
    // DimensionWrapper<6>,
    // DimensionWrapper<7>,
    // DimensionWrapper<8>
>;

// Use TYPED_TEST_SUITE to create the type list for Google Test
TYPED_TEST_SUITE(VolumeUnitSphereTest, DimensionsList);

TYPED_TEST(VolumeUnitSphereTest, ComputeVolume) {
    constexpr int dimensions = TestFixture::dimensions;
    double expected_result = TestFixture::expected_result;
    std::cout << "Current dimension:"  << dimensions << std::endl;
    // Create the integrator with the given dimension
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(volume_unit_sphere_3d<dimensions>), nullptr); // Replace with the actual integrand
    integrator.improve_grid();
    integrator.integrate();

    std::cout << integrator.get_result() << " +/- " << integrator.get_error()
              << " with chi-square: " << integrator.get_chi_square() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}

// TEST(SimpleFunctionTest, volume_unit_sphere_3d) {
//     double expected_result{4.0/3.0*std::numbers::pi};
//     constexpr int dimensions{3};
//     VegasNumericalIntegration<dimensions> integrator;
//     integrator.set_verbosity(VegasVerbosity::None);
//     integrator.set_integrand(std::move(volume_unit_sphere_3d<dimensions>), nullptr);
//     integrator.improve_grid();
//     integrator.integrate();
//     std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
//               << integrator.get_chi_square() << std::endl;
//
//     EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
// }

TEST(SimpleFunctionTest, surface_unit_sphere_3d) {
    constexpr double expected_result{4.0*std::numbers::pi};
    constexpr int dimensions{3};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(surface_unit_sphere_3d<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chi_square() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}


TEST(SimpleFunctionTest, arctan_derivative) {
    constexpr double expected_result{std::numbers::pi*std::numbers::pi*std::numbers::pi/64.0};
    constexpr int dimensions{3};
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(arctan_derivative<dimensions>), nullptr);
    integrator.improve_grid();
    integrator.integrate();
    std::cout << integrator.get_result() << " +/- " << integrator.get_error() << " with chi-square: "
              << integrator.get_chi_square() << std::endl;

    EXPECT_NEAR(expected_result, integrator.get_result(), sigma_range * integrator.get_error());
}
