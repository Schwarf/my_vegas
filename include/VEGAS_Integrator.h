#ifndef VEGAS_INTEGRATOR_H
#define VEGAS_INTEGRATOR_H

#include "VEGAS_map.h"
#include "VEGAS_Stratify.h"
#include <random> // The random number generator and distributions. c++11
#include <vector>
#include <functional>
#include <any>

template <int NumberOfDimensions>
using VEGAS_INTEGRAND = std::function<double(const std::array<double, NumberOfDimensions> &, void *param)>;
//typedef double (*INTEGRAND)(std::vector<double> x, void *param);

enum class VegasVerbosity {
    None = 0,
    Info = 1,
    Debug = 2,
    All = 3
};

template<int NumberOfDimensions>
class VegasNumericalIntegration {
public:
    VegasNumericalIntegration() =default;
    ~VegasNumericalIntegration() = default;
    void set_integrand(VEGAS_INTEGRAND<NumberOfDimensions> &&integrand, void *parameters);
    void improve_grid();
    void integrate(double eps_rel = 1e-3, double eps_abs = 1e-9);
    void set_verbosity(const VegasVerbosity & verbose) {verbosity = verbose;}
    double get_result();
    double get_error();
    double get_chi_square();
private:
    void *integrand_parameters{};
    VegasVerbosity verbosity{VegasVerbosity::Info};
    VEGAS_INTEGRAND<NumberOfDimensions> function_integrand;
    VegasMap<NumberOfDimensions> mapping{};
    VEGAS_Stratify<NumberOfDimensions> stratification;
    std::mt19937 random_number_generator; // Mersenne twister random number engine
    std::uniform_real_distribution<double> distribution; // uniform distribution in double in [0.0, 1.0)
    std::vector<double> results;
    std::vector<double> sigma2;
};

#include "VEGAS_Integrator.inl"

#endif //VEGAS_INTEGRATOR_H