#ifndef VEGAS_INTEGRATOR_H
#define VEGAS_INTEGRATOR_H

#include "VEGAS_map.h"
#include "VEGAS_Stratify.h"
#include <random> // The random number generator and distributions. c++11
#include <vector>
#include <functional>
#include <any>

using VEGAS_INTEGRAND = std::function<double(const std::vector<double>&, void* param)>;
//typedef double (*INTEGRAND)(std::vector<double> x, void *param);

enum class VegasVerbosity
{
    None = 0,
    Info = 1,
    Debug = 2,
    All = 3
};

template<int NumberOfDimensions>
class VegasNumericalIntegration
{
private:
    VegasVerbosity verbosity;

    VEGAS_INTEGRAND function_integrand;
    void* userdata;

    VegasMap map;
    VEGAS_Stratify<NumberOfDimensions> strat;

    std::mt19937 random_number_generator; // Mersenne twister random number engine
    std::uniform_real_distribution<double> distribution; // uniform distribution in double in [0.0, 1.0)

    std::vector<double> results;
    std::vector<double> sigma2;


public:
    VegasNumericalIntegration(){ verbosity = VegasVerbosity::Info;};
    ~VegasNumericalIntegration() = default;

    void Set_Verbose(VegasVerbosity level);

    void set_integrand(VEGAS_INTEGRAND && integrand, void* param);
    void improve_grid();
    void integrate(double eps_rel = 1e-3, double eps_abs = 1e-9);
    
    
    double get_result();
    double get_error();
    double get_chisquare();

};

#include "VEGAS_Integrator.inl"
#endif //VEGAS_INTEGRATOR_H