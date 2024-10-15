#include "VEGAS_Integrator.h"
#include <iostream>
#include <iomanip>
#include <chrono>


void VegasNumericalIntegration::Set_Verbose(VegasVerbosity level)
{
    verbosity = level;
}

void VegasNumericalIntegration::set_integrand(VEGAS_INTEGRAND && integrand, int dimensions, void* param)
{
    function_integrand = integrand;
    number_of_dimensions = dimensions;
    userdata = param;
    results.clear();
    sigma2.clear();
    map = VegasMap(dimensions);
    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    random_number_generator.seed(seed);
}

void VegasNumericalIntegration::improve_grid()
{
    std::vector<double> random_numbers(number_of_dimensions);
    std::vector<double> y(number_of_dimensions); // Random number between 0 to 1;
    std::vector<double> x(number_of_dimensions); // The argument for integrand;
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

    strat.Set_Dimension(number_of_dimensions);
    double dV = strat.Get_V_Cubic();

    map.set_alpha(alpha_start);
    // Warm Up with just MAP improvement
    if (verbosity >= VegasVerbosity::Info)
    {
        std::cout<<"======================================================================================"<<std::endl;
        std::cout<<"| Warm Up the VEGAS Map                                                              |"<<std::endl;
        std::cout<<"======================================================================================"<<std::endl;
        std::cout<<"|  Iter  |    N_Eval    |     Result     |      Error     |    Acc    |  Map Changes |"<<std::endl;
    }
    for (int warm_iter = 0; warm_iter < 5; warm_iter++)
    {
        results.push_back(0);
        sigma2.push_back(0);
        Jf = 0;
        Jf2 = 0;
        for (int evaluation{}; evaluation < starting_number_of_evaluations; ++evaluation)
        {
            for (int dimension_index{}; dimension_index < number_of_dimensions; ++dimension_index)
            {
                random_numbers[dimension_index] = distribution(random_number_generator);
            }
            x = map.get_x(random_numbers);
            evaluated_integrand_value = function_integrand(x, userdata);
            jacobian = map.get_jacobian(random_numbers);
            if (std::isnan(evaluated_integrand_value) || std::isnan(jacobian))
            {
                evaluation--;
                continue;
            }
            map.accumulate_weight(random_numbers, evaluated_integrand_value);
            Jf += evaluated_integrand_value * jacobian;
            Jf2 += (evaluated_integrand_value * jacobian) * (evaluated_integrand_value * jacobian);
        }
        Ih = Jf / starting_number_of_evaluations;
        Sig2 = Jf2 / starting_number_of_evaluations - Ih * Ih;
        results[results.size() - 1] += Ih;
        sigma2[sigma2.size() - 1] += Sig2 / starting_number_of_evaluations;
        map.update_map();
        accuracy = sqrt(sigma2[sigma2.size() - 1]) / results[results.size() - 1];
        if (verbosity >= VegasVerbosity::Info)
        {
            std::cout << "| " << std::setw(6) << warm_iter << " | " << std::setw(12) << starting_number_of_evaluations << " | " << std::setw(14) << std::scientific << std::setprecision(5) << results[results.size() - 1] << " | " << std::setw(14) << std::scientific << std::setprecision(5) << sqrt(sigma2[sigma2.size() - 1]) << " | " << resetiosflags(std::ios::scientific) << std::fixed << std::setw(8) << std::setprecision(3) << accuracy * 100 << "% | " << resetiosflags(std::ios::fixed) << std::setw(12) << std::scientific << std::setprecision(5) << map.checking_map() << " |" << std::endl;
        }
    }
    result = get_result();
    error = get_error();
    accuracy = error / result;
    if (verbosity >= VegasVerbosity::Info)
    {
        std::cout << "| Summary of Warm up 5 Iter:   result = " << std::setw(11) << std::scientific << std::setprecision(5) << result << "   error = " << std::setw(11) << std::scientific << std::setprecision(5) << error << "   Acc = " << resetiosflags(std::ios::scientific) << std::fixed << std::setw(6) << std::setprecision(3) << accuracy * 100 << "% |" << std::endl;
    }
    results.clear();
    sigma2.clear();

    if (verbosity >= VegasVerbosity::Info)
    {
        std::cout<<"======================================================================================"<<std::endl;
        std::cout<<"| Improving the mapping grid and stratification grid                                 |"<<std::endl;
        std::cout<<"======================================================================================"<<std::endl;
        std::cout<<"|  Iter  |    N_Eval    |     Result     |      Error     |    Acc    |  Map Changes |"<<std::endl;
    }
    while (true)
    {
        // we decide to end the improvement of grid and strata when the accuracy is about 1%
        // Every 5 iteration, we can check the accuracy, and addjust the number of evaluation
        // Map and Strata improves every another iteration.
        number_of_iterations++;
        strat.Set_NEVAL(starting_number_of_evaluations);
        results.push_back(0);
        sigma2.push_back(0);
        NEVAL_REAL = 0;
        for (int inc = 0; inc < strat.Get_NHYPERCUBICS(); inc++)
        {
            Jf = 0;
            Jf2 = 0;
            number_of_evaluations = strat.Get_NH(inc);
            NEVAL_REAL += number_of_evaluations;
            for (int evaluation{}; evaluation < number_of_evaluations; ++evaluation)
            {
                for (int i_dim = 0; i_dim < number_of_dimensions; i_dim++)
                {
                    random_numbers[i_dim] = distribution(random_number_generator);
                }
                y = strat.Get_Y(inc, random_numbers);
                x = map.get_x(y);
                evaluated_integrand_value = function_integrand(x, userdata);
                jacobian = map.get_jacobian(y);
                if (std::isnan(evaluated_integrand_value) || std::isnan(jacobian))
                {
                    evaluation--;
                    continue;
                }
                map.accumulate_weight(y, evaluated_integrand_value);
                strat.Accumulate_Weight(inc, evaluated_integrand_value * jacobian);
                Jf += evaluated_integrand_value * jacobian;
                Jf2 += (evaluated_integrand_value * jacobian) * (evaluated_integrand_value * jacobian);
            }
            Ih = Jf / number_of_evaluations * dV;
            Sig2 = Jf2 / number_of_evaluations * dV * dV - Ih * Ih;
            results[results.size() - 1] += Ih;
            sigma2[sigma2.size() - 1] += Sig2 / number_of_evaluations;
        }
        if (number_of_iterations % 2 != 0)
        {
            // if (alpha > 0.05)
            // {
            map.update_map();
                // alpha = alpha_start*exp(-number_of_iterations/5.0);
                // map.set_alpha(alpha);
            // }
        }
        else
        {
            strat.Update_DH();
        }
        accuracy = sqrt(sigma2[sigma2.size() - 1]) / results[results.size() - 1];
        if (verbosity >= VegasVerbosity::Info)
        {
            std::cout << "| " << std::setw(6) << number_of_iterations << " | " << std::setw(12) << NEVAL_REAL << " | " << std::setw(14) << std::scientific << std::setprecision(5) << results[results.size() - 1] << " | " << std::setw(14) << std::scientific << std::setprecision(5) << sqrt(sigma2[sigma2.size() - 1]) << " | " << resetiosflags(std::ios::scientific) << std::fixed << std::setw(8) << std::setprecision(3) << accuracy * 100 << "% | " << resetiosflags(std::ios::fixed) << std::setw(12) << std::scientific << std::setprecision(5) << map.checking_map() << " |" << std::endl;
        }
        if (number_of_iterations % 5 == 0)
        {
            result = get_result();
            error = get_error();
            accuracy = error / result;
            if (verbosity >= VegasVerbosity::Info)
            {
                std::cout << "| Summary of Last 5 Iter:      result = " << std::setw(11) << std::scientific << std::setprecision(5) << result << "   error = " << std::setw(11) << std::scientific << std::setprecision(5) << error << "   Acc = " << resetiosflags(std::ios::scientific) << std::fixed << std::setw(6) << std::setprecision(3) << accuracy * 100 << "% |" << std::endl;
            }
            if (accuracy < 0.01)
            {
                break;
            }
            starting_number_of_evaluations *= static_cast<int>(sqrt(accuracy / 0.01));
            results.clear();
            sigma2.clear();
        }
    }
    if (verbosity >= VegasVerbosity::Info)
    {
        std::cout<<"======================================================================================"<<std::endl;
    }
}
void VegasNumericalIntegration::integrate(double eps_rel, double eps_abs)
{
    // We try to reach either relative error (eps_rel) or absolute error (eps_abs)
    // But we also need to make sure chi2 is not bigger than the iteration numbers
    std::vector<double> random_numbers(number_of_dimensions);
    std::vector<double> y(number_of_dimensions); // Random number between 0 to 1;
    std::vector<double> x(number_of_dimensions); // The argument for integrand;
    double evaluated_integrand_value{}; // evaluated integrand value;
    double jacobian{}; // The Jacobian from y to x;
    int starting_number_of_evaluations{50000};
    double dV = strat.Get_V_Cubic();
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
    if (verbosity >= VegasVerbosity::Info)
    {
        std::cout<<"======================================================================="<<std::endl;
        std::cout<<"| Fixing the mapping grid, still improve strata and Integral          |"<<std::endl;
        std::cout<<"======================================================================="<<std::endl;
        std::cout<<"|  Iter  |    N_Eval    |     Result     |      Error     |    Acc    |"<<std::endl;
    }
    while (true)
    {
        number_of_iterations++;
        strat.Set_NEVAL(starting_number_of_evaluations);
        results.push_back(0);
        sigma2.push_back(0);
        for (int inc = 0; inc < strat.Get_NHYPERCUBICS(); inc++)
        {
            Jf = 0.0;
            Jf2 = 0.0;
            number_of_evaluations = strat.Get_NH(inc);
            number_of_real_evaluations += number_of_evaluations;
            for (int evaluation{}; evaluation < number_of_evaluations; evaluation++)
            {
                for (int dimension_index = 0; dimension_index < number_of_dimensions; dimension_index++)
                {
                    random_numbers[dimension_index] = distribution(random_number_generator);
                }
                y = strat.Get_Y(inc, random_numbers);
                x = map.get_x(y);
                evaluated_integrand_value = function_integrand(x, userdata);
                jacobian = map.get_jacobian(y);
                if (std::isnan(evaluated_integrand_value) || std::isnan(jacobian))
                {
                    evaluation--;
                    continue;
                }
                strat.Accumulate_Weight(inc, evaluated_integrand_value * jacobian);
                Jf += evaluated_integrand_value * jacobian;
                Jf2 += (evaluated_integrand_value * jacobian) * (evaluated_integrand_value * jacobian);
            }
            Ih = Jf / number_of_evaluations * dV;
            Sig2 = Jf2 / number_of_evaluations * dV * dV - Ih*Ih;
            results[results.size() - 1] += Ih;
            sigma2[sigma2.size() - 1] += Sig2 / number_of_evaluations;
        }
        strat.Update_DH();
        accuracy = sqrt(sigma2[sigma2.size() - 1]) / results[results.size() - 1];
        if (verbosity >= VegasVerbosity::Info)
        {
            std::cout << "| " << std::setw(6) << number_of_iterations << " | " << std::setw(12) << number_of_real_evaluations << " | " << std::setw(14) << std::scientific << std::setprecision(5) << results[results.size() - 1] << " | " << std::setw(14) << std::scientific << std::setprecision(5) << sqrt(sigma2[sigma2.size() - 1]) << " | " << resetiosflags(std::ios::scientific) << std::fixed << std::setw(8) << std::setprecision(3) << accuracy * 100 << "% |" << std::endl;
        }
        if (number_of_iterations % 5 == 0)
        {
            // Every 5 iteration, we check whether we fullfil the condition
            result = get_result();
            error = get_error();
            chi_squared = get_chisquare();
            accuracy = error / result;
            if (verbosity >= VegasVerbosity::Info)
            {
                std::cout << "| Summary of Last 5 Iter: " << std::setw(14) << std::scientific << std::setprecision(5) << result << " | " << std::setw(14) << std::scientific << std::setprecision(5) << error << " | " << resetiosflags(std::ios::scientific) << std::fixed << std::setw(8) << std::setprecision(3) << accuracy * 100 << "% | chi_squared = " << chi_squared << std::endl;
            }
            if ((accuracy < eps_rel || error < eps_abs) && chi_squared / 5.0 < 1.0 )
            {
                break;
            }
            if (chi_squared / 5.0 < 1.0)
            {
                starting_number_of_evaluations *= static_cast<int>(sqrt(accuracy / eps_rel));
                results.clear();
                sigma2.clear();
                continue;
            }
            if (chi_squared / 5.0 > 1.0)
            {
                starting_number_of_evaluations += 5000;
                results.clear();
                sigma2.clear();
                continue;
            }
        }
    }
    if (verbosity >= VegasVerbosity::Info)
    {
        result = get_result();
        error = get_error();
        chi_squared = get_chisquare();
        accuracy = error / result;
        std::cout<<"======================================================================="<<std::endl;
        std::cout<<"Summary: "<<std::endl;
        std::cout << "Result: " << std::setw(12) << std::scientific << std::setprecision(5) << result << "  Error: " << std::setw(12) << std::scientific << std::setprecision(5) << error << "  Acc: " << resetiosflags(std::ios::scientific) << std::fixed << std::setw(6) << std::setprecision(3) << accuracy * 100 << "%  chi_squared: " << chi_squared << std::endl;
        std::cout<<"======================================================================="<<std::endl;
        std::cout<<resetiosflags(std::ios::fixed)<<std::setprecision(8);
    }
}
double VegasNumericalIntegration::get_result()
{
    double res_num = 0;
    double res_den = 0;
    for (int i = 0; i < results.size(); i++)
    {
        res_num += results[i] / sigma2[i];
        res_den += 1.0 / sigma2[i];
    }
    return res_num/res_den;
}
double VegasNumericalIntegration::get_error()
{
    double res = 0;
    for (double sigma_element : sigma2)
    {
        res += 1.0 / sigma_element;
    }
    return 1.0/sqrt(res);
}
double VegasNumericalIntegration::get_chisquare()
{
    double Ifinal = get_result();
    double chi2 = 0;
    for (int i = 0; i < results.size(); i++)
    {
        chi2 += (results[i] - Ifinal)*(results[i] - Ifinal) / sigma2[i];
    }
    return chi2;
}