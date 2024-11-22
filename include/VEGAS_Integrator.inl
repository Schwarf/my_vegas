#include "VEGAS_Integrator.h"
#include <iostream>
#include <iomanip>
#include <chrono>


template<int NumberOfDimensions>
void
VegasNumericalIntegration<NumberOfDimensions>::set_integrand(VEGAS_INTEGRAND<NumberOfDimensions> &&integrand,
                                                             void *parameters) {
    function_integrand = integrand;
    integrand_parameters = parameters;
    results.clear();
    sigma2.clear();
    mapping = VegasMap<NumberOfDimensions>();
    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    random_number_generator.seed(seed);
}

template<int NumberOfDimensions>
void VegasNumericalIntegration<NumberOfDimensions>::improve_grid() {
    std::array<double, NumberOfDimensions> random_numbers{};
    std::array<double, NumberOfDimensions> x{}; // The argument for integrand;
    std::array<double, NumberOfDimensions> y{}; // Random number between 0 to 1;
    double evaluated_integrand_value{}; // evaluated integrand value;
    double jacobian{};
    int number_of_iterations{};
    int starting_number_of_evaluations = 10000;
    constexpr double alpha_start{0.5};
    double result{};
    double error{};
    int number_of_evaluations{};
    int NEVAL_REAL{};
    double Jf{};
    double Jf2{};
    double Ih{};
    double Sig2{};
    double accuracy{};

    double dV = stratification.get_V_cubic();

    mapping.set_alpha(alpha_start);
    // Warm Up with just MAP improvement
    if (verbosity >= VegasVerbosity::Info) {
        std::cout << "======================================================================================"
                  << std::endl;
        std::cout << "| Warm Up the VEGAS Map                                                              |"
                  << std::endl;
        std::cout << "======================================================================================"
                  << std::endl;
        std::cout << "|  Iter  |    N_Eval    |     Result     |      Error     |    Acc    |  Map Changes |"
                  << std::endl;
    }
    for (int warm_iter = 0; warm_iter < 5; warm_iter++) {
        results.push_back(0);
        sigma2.push_back(0);
        Jf = 0;
        Jf2 = 0;
        for (int evaluation{}; evaluation < starting_number_of_evaluations; ++evaluation) {
            for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
                random_numbers[dimension] = distribution(random_number_generator);
            }
            x = mapping.get_x(random_numbers);
            evaluated_integrand_value = function_integrand(x, integrand_parameters);
            jacobian = mapping.get_jacobian();
            if (std::isnan(evaluated_integrand_value) || std::isnan(jacobian)) {
                evaluation--;
                continue;
            }
            mapping.accumulate_weights(evaluated_integrand_value);
            Jf += evaluated_integrand_value * jacobian;
            Jf2 += (evaluated_integrand_value * jacobian) * (evaluated_integrand_value * jacobian);
        }
        Ih = Jf / starting_number_of_evaluations;
        Sig2 = Jf2 / starting_number_of_evaluations - Ih * Ih;
        results[results.size() - 1] += Ih;
        sigma2[sigma2.size() - 1] += Sig2 / starting_number_of_evaluations;
        mapping.update_map();
        accuracy = sqrt(sigma2[sigma2.size() - 1]) / results[results.size() - 1];
        if (verbosity >= VegasVerbosity::Info) {
            std::cout << "| " << std::setw(6) << warm_iter << " | " << std::setw(12) << starting_number_of_evaluations
                      << " | " << std::setw(14) << std::scientific << std::setprecision(5)
                      << results[results.size() - 1] << " | " << std::setw(14) << std::scientific
                      << std::setprecision(5) << sqrt(sigma2[sigma2.size() - 1]) << " | "
                      << resetiosflags(std::ios::scientific) << std::fixed << std::setw(8) << std::setprecision(3)
                      << accuracy * 100 << "% | " << resetiosflags(std::ios::fixed) << std::setw(12) << std::scientific
                      << std::setprecision(5) << mapping.checking_map() << " |" << std::endl;
        }
    }
    result = get_result();
    error = get_error();
    accuracy = error / result;
    if (verbosity >= VegasVerbosity::Info) {
        std::cout << "| Summary of Warm up 5 Iter:   result = " << std::setw(11) << std::scientific
                  << std::setprecision(5) << result << "   error = " << std::setw(11) << std::scientific
                  << std::setprecision(5) << error << "   Acc = " << resetiosflags(std::ios::scientific) << std::fixed
                  << std::setw(6) << std::setprecision(3) << accuracy * 100 << "% |" << std::endl;
    }
    results.clear();
    sigma2.clear();

    if (verbosity >= VegasVerbosity::Info) {
        std::cout << "======================================================================================"
                  << std::endl;
        std::cout << "| Improving the mapping grid and stratification grid                                 |"
                  << std::endl;
        std::cout << "======================================================================================"
                  << std::endl;
        std::cout << "|  Iter  |    N_Eval    |     Result     |      Error     |    Acc    |  Map Changes |"
                  << std::endl;
    }
    while (true) {
        // we decide to end the improvement of grid and strata when the accuracy is about 1%
        // Every 5 iteration, we can check the accuracy, and addjust the number of evaluation
        // Map and Strata improves every another iteration.
        number_of_iterations++;
        stratification.Set_NEVAL(starting_number_of_evaluations);
        results.push_back(0);
        sigma2.push_back(0);
        NEVAL_REAL = 0;
        for (int cube_index = 0; cube_index < stratification.Get_NHYPERCUBICS(); cube_index++) {
            Jf = 0.0;
            Jf2 = 0.0;
            number_of_evaluations = stratification.get_expected_events_per_hyper_cube(cube_index);
            NEVAL_REAL += number_of_evaluations;
            for (int evaluation{}; evaluation < number_of_evaluations; ++evaluation) {
                for (int i_dim = 0; i_dim < NumberOfDimensions; i_dim++) {
                    random_numbers[i_dim] = distribution(random_number_generator);
                }
                y = stratification.get_y(cube_index, random_numbers);
                x = mapping.get_x(y);
                evaluated_integrand_value = function_integrand(x, integrand_parameters);
                jacobian = mapping.get_jacobian();
                if (std::isnan(evaluated_integrand_value) || std::isnan(jacobian)) {
                    evaluation--;
                    continue;
                }
                mapping.accumulate_weights(evaluated_integrand_value);
                stratification.accumulate_weights(cube_index, evaluated_integrand_value * jacobian);
                Jf += evaluated_integrand_value * jacobian;
                Jf2 += (evaluated_integrand_value * jacobian) * (evaluated_integrand_value * jacobian);
            }
            Ih = Jf / number_of_evaluations * dV;
            Sig2 = Jf2 / number_of_evaluations * dV * dV - Ih * Ih;
            results[results.size() - 1] += Ih;
            sigma2[sigma2.size() - 1] += Sig2 / number_of_evaluations;
        }
        if (number_of_iterations % 2 != 0) {
            // if (alpha > 0.05)
            // {
            mapping.update_map();
            // alpha = alpha_start*exp(-number_of_iterations/5.0);
            // mapping.set_alpha(alpha);
            // }
        } else {
            stratification.update_hypercubic_weights();
        }
        accuracy = sqrt(sigma2[sigma2.size() - 1]) / results[results.size() - 1];
        if (verbosity >= VegasVerbosity::Info) {
            std::cout << "| " << std::setw(6) << number_of_iterations << " | " << std::setw(12) << NEVAL_REAL << " | "
                      << std::setw(14) << std::scientific << std::setprecision(5) << results[results.size() - 1]
                      << " | " << std::setw(14) << std::scientific << std::setprecision(5)
                      << sqrt(sigma2[sigma2.size() - 1]) << " | " << resetiosflags(std::ios::scientific) << std::fixed
                      << std::setw(8) << std::setprecision(3) << accuracy * 100 << "% | "
                      << resetiosflags(std::ios::fixed) << std::setw(12) << std::scientific << std::setprecision(5)
                      << mapping.checking_map() << " |" << std::endl;
        }
        if (number_of_iterations % 5 == 0) {
            result = get_result();
            error = get_error();
            accuracy = error / result;
            if (verbosity >= VegasVerbosity::Info) {
                std::cout << "| Summary of Last 5 Iter:      result = " << std::setw(11) << std::scientific
                          << std::setprecision(5) << result << "   error = " << std::setw(11) << std::scientific
                          << std::setprecision(5) << error << "   Acc = " << resetiosflags(std::ios::scientific)
                          << std::fixed << std::setw(6) << std::setprecision(3) << accuracy * 100 << "% |" << std::endl;
            }
            if (accuracy < 0.01) {
                break;
            }
            starting_number_of_evaluations *= static_cast<int>(sqrt(accuracy / 0.01));
            results.clear();
            sigma2.clear();
        }
    }
    if (verbosity >= VegasVerbosity::Info) {
        std::cout << "======================================================================================"
                  << std::endl;
    }
}

template<int NumberOfDimensions>
void VegasNumericalIntegration<NumberOfDimensions>::integrate(double eps_rel, double eps_abs) {
    // We try to reach either relative error (eps_rel) or absolute error (eps_abs)
    // But we also need to make sure chi2 is not bigger than the iteration numbers
    std::array<double, NumberOfDimensions> random_numbers{};
    std::array<double, NumberOfDimensions> y{}; // Random number between 0 to 1;
    std::array<double, NumberOfDimensions> x{}; // The argument for integrand;
    double evaluated_integrand_value{}; // evaluated integrand value;
    double jacobian{}; // The Jacobian from y to x;
    int starting_number_of_evaluations{50000};
    double dV = stratification.get_V_cubic();
    int number_of_iterations{};
    double result{};
    double error{};
    double chi_squared{};
    int number_of_evaluations{};
    int number_of_real_evaluations{};
    double Jf{};
    double Jf2{};
    double Ih{};
    double Sig2{};
    double accuracy{};
    if (verbosity >= VegasVerbosity::Info) {
        std::cout << "=======================================================================" << std::endl;
        std::cout << "| Fixing the mapping grid, still improve strata and Integral          |" << std::endl;
        std::cout << "=======================================================================" << std::endl;
        std::cout << "|  Iter  |    N_Eval    |     Result     |      Error     |    Acc    |" << std::endl;
    }
    while (true) {
        number_of_iterations++;
        stratification.Set_NEVAL(starting_number_of_evaluations);
        results.push_back(0);
        sigma2.push_back(0);
        for (int cube_index = 0; cube_index < stratification.Get_NHYPERCUBICS(); cube_index++) {
            Jf = 0.0;
            Jf2 = 0.0;
            number_of_evaluations = stratification.get_expected_events_per_hyper_cube(cube_index);
            number_of_real_evaluations += number_of_evaluations;
            for (int evaluation{}; evaluation < number_of_evaluations; evaluation++) {
                for (int dimension_index = 0; dimension_index < NumberOfDimensions; dimension_index++) {
                    random_numbers[dimension_index] = distribution(random_number_generator);
                }
                y = stratification.get_y(cube_index, random_numbers);
                x = mapping.get_x(y);
                evaluated_integrand_value = function_integrand(x, integrand_parameters);
                jacobian = mapping.get_jacobian();
                if (std::isnan(evaluated_integrand_value) || std::isnan(jacobian)) {
                    evaluation--;
                    continue;
                }
                stratification.accumulate_weights(cube_index, evaluated_integrand_value * jacobian);
                Jf += evaluated_integrand_value * jacobian;
                Jf2 += (evaluated_integrand_value * jacobian) * (evaluated_integrand_value * jacobian);
            }
            Ih = Jf / number_of_evaluations * dV;
            Sig2 = Jf2 / number_of_evaluations * dV * dV - Ih * Ih;
            results[results.size() - 1] += Ih;
            sigma2[sigma2.size() - 1] += Sig2 / number_of_evaluations;
        }
        stratification.update_hypercubic_weights();
        accuracy = sqrt(sigma2[sigma2.size() - 1]) / results[results.size() - 1];
        if (verbosity >= VegasVerbosity::Info) {
            std::cout << "| " << std::setw(6) << number_of_iterations << " | " << std::setw(12)
                      << number_of_real_evaluations << " | " << std::setw(14) << std::scientific << std::setprecision(5)
                      << results[results.size() - 1] << " | " << std::setw(14) << std::scientific
                      << std::setprecision(5) << sqrt(sigma2[sigma2.size() - 1]) << " | "
                      << resetiosflags(std::ios::scientific) << std::fixed << std::setw(8) << std::setprecision(3)
                      << accuracy * 100 << "% |" << std::endl;
        }
        if (number_of_iterations % 5 == 0) {
            // Every 5 iteration, we check whether we fullfil the condition
            result = get_result();
            error = get_error();
            chi_squared = get_chisquare();
            accuracy = error / result;
            if (verbosity >= VegasVerbosity::Info) {
                std::cout << "| Summary of Last 5 Iter: " << std::setw(14) << std::scientific << std::setprecision(5)
                          << result << " | " << std::setw(14) << std::scientific << std::setprecision(5) << error
                          << " | " << resetiosflags(std::ios::scientific) << std::fixed << std::setw(8)
                          << std::setprecision(3) << accuracy * 100 << "% | chi_squared = " << chi_squared << std::endl;
            }
            if ((accuracy < eps_rel || error < eps_abs) && chi_squared / 5.0 < 1.0) {
                break;
            }
            if (chi_squared / 5.0 < 1.0) {
                starting_number_of_evaluations *= static_cast<int>(sqrt(accuracy / eps_rel));
                results.clear();
                sigma2.clear();
                continue;
            }
            if (chi_squared / 5.0 > 1.0) {
                starting_number_of_evaluations += 5000;
                results.clear();
                sigma2.clear();
                continue;
            }
        }
    }
    if (verbosity >= VegasVerbosity::Info) {
        result = get_result();
        error = get_error();
        chi_squared = get_chisquare();
        accuracy = error / result;
        std::cout << "=======================================================================" << std::endl;
        std::cout << "Summary: " << std::endl;
        std::cout << "Result: " << std::setw(12) << std::scientific << std::setprecision(5) << result << "  Error: "
                  << std::setw(12) << std::scientific << std::setprecision(5) << error << "  Acc: "
                  << resetiosflags(std::ios::scientific) << std::fixed << std::setw(6) << std::setprecision(3)
                  << accuracy * 100 << "%  chi_squared: " << chi_squared << std::endl;
        std::cout << "=======================================================================" << std::endl;
        std::cout << resetiosflags(std::ios::fixed) << std::setprecision(8);
    }
}

template<int NumberOfDimensions>
double VegasNumericalIntegration<NumberOfDimensions>::get_result() {
    double res_num = 0;
    double res_den = 0;
    for (size_t i{}; i < results.size(); i++) {
        res_num += results[i] / sigma2[i];
        res_den += 1.0 / sigma2[i];
    }
    return res_num / res_den;
}

template<int NumberOfDimensions>
double VegasNumericalIntegration<NumberOfDimensions>::get_error() {
    double res = 0;
    for (double sigma_element: sigma2) {
        res += 1.0 / sigma_element;
    }
    return 1.0 / sqrt(res);
}

template<int NumberOfDimensions>
double VegasNumericalIntegration<NumberOfDimensions>::get_chisquare() {
    double Ifinal = get_result();
    double chi2 = 0;
    for (size_t i{}; i < results.size(); i++) {
        chi2 += (results[i] - Ifinal) * (results[i] - Ifinal) / sigma2[i];
    }
    return chi2;
}