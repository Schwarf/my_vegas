#include "VEGAS_Integrator.h"
#include <iostream>
#include <iomanip>
#include <chrono>


void VEGAS_Integrator::Set_Verbose(VegasVerbosity level)
{
    verbosity = level;
}

void VEGAS_Integrator::Set_Integrand(VEGAS_INTEGRAND && integrand, int dim, void* param)
{
    function_integrand = integrand;
    dimensions = dim;
    userdata = param;
    Results.clear();
    Sigma2.clear();
    map = VEGAS_Map(dim);
    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    random_number_generator.seed(seed);
}

void VEGAS_Integrator::Improve_Grid()
{
    std::vector<double> yrnd(dimensions);
    std::vector<double> y(dimensions); // Random number between 0 to 1;
    std::vector<double> x(dimensions); // The argument for integrand;
    double f_eval; // evaluated integrand value;
    double Jac; 
    int iter = 0;
    int NEVAL_START = 10000;
    double alpha_start = 0.5;
    double alpha = alpha_start;
    double dV;
    double Res;
    double Err;
    int neval;
    int NEVAL_REAL;
    double Jf;
    double Jf2;
    double Ih;
    double Sig2;
    double acc;
    strat.Set_Dimension(dimensions);
    dV = strat.Get_V_Cubic();
    map.Set_alpha(alpha_start);
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
        Results.push_back(0);
        Sigma2.push_back(0);
        Jf = 0;
        Jf2 = 0;
        for (int ne = 0; ne < NEVAL_START; ne++)
        {
            for (int i_dim = 0; i_dim < dimensions; i_dim++)
            {
                yrnd[i_dim] = distribution(random_number_generator);
            }
            x = map.Get_X(yrnd);
            f_eval = function_integrand(x, userdata);
            Jac = map.Get_Jac(yrnd);
            if (std::isnan(f_eval) || std::isnan(Jac))
            {
                ne--;
                continue;
            }
            map.Accumulate_Weight(yrnd,f_eval);
            Jf += f_eval*Jac;
            Jf2 += pow(f_eval*Jac,2);
        }
        Ih = Jf/NEVAL_START;
        Sig2 = Jf2/NEVAL_START - pow(Jf/NEVAL_START,2);
        Results[Results.size()-1] += Ih;
        Sigma2[Sigma2.size()-1] += Sig2/NEVAL_START;
        map.Update_Map();
        acc = sqrt(Sigma2[Sigma2.size()-1])/Results[Results.size()-1];
        if (verbosity >= VegasVerbosity::Info)
        {
            std::cout<<"| "<<std::setw(6)<<warm_iter<<" | "<<std::setw(12)<<NEVAL_START<<" | "<<std::setw(14)<<std::scientific<<std::setprecision(5)<<Results[Results.size()-1]<<" | "<<std::setw(14)<<std::scientific<<std::setprecision(5)<<sqrt(Sigma2[Sigma2.size()-1])<<" | "<<resetiosflags(std::ios::scientific)<<std::fixed<<std::setw(8)<<std::setprecision(3)<<acc*100<<"% | "<<resetiosflags(std::ios::fixed)<<std::setw(12)<<std::scientific<<std::setprecision(5)<<map.Checking_Map()<<" |"<<std::endl;
        }
    }
    Res = Get_Result();
    Err = Get_Error();
    acc = Err/Res;
    if (verbosity >= VegasVerbosity::Info)
    {
        std::cout<<"| Summary of Warm up 5 Iter:   Res = "<<std::setw(11)<<std::scientific<<std::setprecision(5)<< Res <<"   Err = "<<std::setw(11)<<std::scientific<<std::setprecision(5)<< Err <<"   Acc = "<<resetiosflags(std::ios::scientific)<<std::fixed<<std::setw(6)<<std::setprecision(3)<<acc*100<<"% |"<<std::endl;
    }
    Results.clear();
    Sigma2.clear();

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
        iter++;
        strat.Set_NEVAL(NEVAL_START);
        Results.push_back(0);
        Sigma2.push_back(0);
        NEVAL_REAL = 0;
        for (int inc = 0; inc < strat.Get_NHYPERCUBICS(); inc++)
        {
            Jf = 0;
            Jf2 = 0;
            neval = strat.Get_NH(inc);
            NEVAL_REAL += neval;
            for (int ne = 0; ne < neval; ne++)
            {
                for (int i_dim = 0; i_dim < dimensions; i_dim++)
                {
                    yrnd[i_dim] = distribution(random_number_generator);
                }
                y = strat.Get_Y(inc,yrnd);
                x = map.Get_X(y);
                f_eval = function_integrand(x, userdata);
                Jac = map.Get_Jac(y);
                if (std::isnan(f_eval) || std::isnan(Jac))
                {
                    ne--;
                    continue;
                }
                map.Accumulate_Weight(y,f_eval);
                strat.Accumulate_Weight(inc,f_eval*Jac);
                Jf += f_eval*Jac;
                Jf2 += pow(f_eval*Jac,2);
            }
            Ih = Jf/neval*dV;
            Sig2 = Jf2/neval*dV*dV - pow(Jf/neval*dV,2);
            Results[Results.size()-1] += Ih;
            Sigma2[Sigma2.size()-1] += Sig2/neval;
        }
        if (iter % 2 != 0)
        {
            // if (alpha > 0.05)
            // {
                map.Update_Map();
                // alpha = alpha_start*exp(-iter/5.0);
                // map.Set_alpha(alpha);
            // }
        }
        else
        {
            strat.Update_DH();
        }
        acc = sqrt(Sigma2[Sigma2.size()-1])/Results[Results.size()-1];
        if (verbosity >= VegasVerbosity::Info)
        {
            std::cout<<"| "<<std::setw(6)<<iter<<" | "<<std::setw(12)<<NEVAL_REAL<<" | "<<std::setw(14)<<std::scientific<<std::setprecision(5)<<Results[Results.size()-1]<<" | "<<std::setw(14)<<std::scientific<<std::setprecision(5)<<sqrt(Sigma2[Sigma2.size()-1])<<" | "<<resetiosflags(std::ios::scientific)<<std::fixed<<std::setw(8)<<std::setprecision(3)<<acc*100<<"% | "<<resetiosflags(std::ios::fixed)<<std::setw(12)<<std::scientific<<std::setprecision(5)<<map.Checking_Map()<<" |"<<std::endl;
        }
        if (iter % 5 == 0)
        {
            Res = Get_Result();
            Err = Get_Error();
            acc = Err/Res;
            if (verbosity >= VegasVerbosity::Info)
            {
                std::cout<<"| Summary of Last 5 Iter:      Res = "<<std::setw(11)<<std::scientific<<std::setprecision(5)<< Res <<"   Err = "<<std::setw(11)<<std::scientific<<std::setprecision(5)<< Err <<"   Acc = "<<resetiosflags(std::ios::scientific)<<std::fixed<<std::setw(6)<<std::setprecision(3)<<acc*100<<"% |"<<std::endl;
            }
            if (acc < 0.01)
            {
                break;
            }
            NEVAL_START = NEVAL_START * sqrt(acc/0.01);
            Results.clear();
            Sigma2.clear();
        }
    }
    if (verbosity >= VegasVerbosity::Info)
    {
        std::cout<<"======================================================================================"<<std::endl;
    }
}
void VEGAS_Integrator::Integration(double eps_rel, double eps_abs)
{
    // We try to reach either relative error (eps_rel) or absolute error (eps_abs)
    // But we also need to make sure chi2 is not bigger than the iteration numbers
    std::vector<double> yrnd(dimensions);
    std::vector<double> y(dimensions); // Random number between 0 to 1;
    std::vector<double> x(dimensions); // The argument for integrand;
    double f_eval; // evaluated integrand value;
    double Jac; // The Jacobian from y to x;
    int NEVAL_START = 50000;
    double dV = strat.Get_V_Cubic();
    int iter = 0;
    double Res;
    double Err;
    double Chi2;
    int neval;
    int NEVAL_REAL;
    double Jf;
    double Jf2;
    double Ih;
    double Sig2;
    double acc;
    if (verbosity >= VegasVerbosity::Info)
    {
        std::cout<<"======================================================================="<<std::endl;
        std::cout<<"| Fixing the mapping grid, still improve strata and Integral          |"<<std::endl;
        std::cout<<"======================================================================="<<std::endl;
        std::cout<<"|  Iter  |    N_Eval    |     Result     |      Error     |    Acc    |"<<std::endl;
    }
    while (true)
    {
        iter++;
        strat.Set_NEVAL(NEVAL_START);
        Results.push_back(0);
        Sigma2.push_back(0);
        NEVAL_REAL = 0;
        for (int inc = 0; inc < strat.Get_NHYPERCUBICS(); inc++)
        {
            Jf = 0;
            Jf2 = 0;
            neval = strat.Get_NH(inc);
            NEVAL_REAL += neval;
            for (int ne = 0; ne < neval; ne++)
            {
                for (int i_dim = 0; i_dim < dimensions; i_dim++)
                {
                    yrnd[i_dim] = distribution(random_number_generator);
                }
                y = strat.Get_Y(inc,yrnd);
                x = map.Get_X(y);
                f_eval = function_integrand(x, userdata);
                Jac = map.Get_Jac(y);
                if (std::isnan(f_eval) || std::isnan(Jac))
                {
                    ne--;
                    continue;
                }
                strat.Accumulate_Weight(inc,f_eval*Jac);
                Jf += f_eval*Jac;
                Jf2 += pow(f_eval*Jac,2);
            }
            Ih = Jf/neval*dV;
            Sig2 = Jf2/neval*dV*dV - pow(Jf/neval*dV,2);
            Results[Results.size()-1] += Ih;
            Sigma2[Sigma2.size()-1] += Sig2/neval;
        }
        strat.Update_DH();
        acc = sqrt(Sigma2[Sigma2.size()-1])/Results[Results.size()-1];
        if (verbosity >= VegasVerbosity::Info)
        {
            std::cout<<"| "<<std::setw(6)<<iter<<" | "<<std::setw(12)<<NEVAL_REAL<<" | "<<std::setw(14)<<std::scientific<<std::setprecision(5)<<Results[Results.size()-1]<<" | "<<std::setw(14)<<std::scientific<<std::setprecision(5)<<sqrt(Sigma2[Sigma2.size()-1])<<" | "<<resetiosflags(std::ios::scientific)<<std::fixed<<std::setw(8)<<std::setprecision(3)<<acc*100<<"% |"<<std::endl;
        }
        if (iter%5==0)
        {
            // Every 5 iteration, we check whether we fullfil the condition
            Res = Get_Result();
            Err = Get_Error();
            Chi2 = Get_Chisq();
            acc = Err/Res;
            if (verbosity >= VegasVerbosity::Info)
            {
                std::cout<<"| Summary of Last 5 Iter: "<<std::setw(14)<<std::scientific<<std::setprecision(5)<< Res <<" | "<<std::setw(14)<<std::scientific<<std::setprecision(5)<< Err <<" | "<<resetiosflags(std::ios::scientific)<<std::fixed<<std::setw(8)<<std::setprecision(3)<<acc*100<<"% | Chi2 = "<<Chi2<<std::endl;
            }
            if ( (acc < eps_rel || Err < eps_abs) && Chi2/5.0 < 1.0 )
            {
                break;
            }
            if (Chi2/5.0 < 1.0)
            {
                NEVAL_START = NEVAL_START * sqrt(acc/eps_rel);
                Results.clear();
                Sigma2.clear();
                continue;
            }
            if (Chi2/5.0 > 1.0)
            {
                NEVAL_START += 5000;
                Results.clear();
                Sigma2.clear();
                continue;
            }
        }
    }
    if (verbosity >= VegasVerbosity::Info)
    {
        Res = Get_Result();
        Err = Get_Error();
        Chi2 = Get_Chisq();
        acc = Err/Res;
        std::cout<<"======================================================================="<<std::endl;
        std::cout<<"Summary: "<<std::endl;
        std::cout<<"Result: "<<std::setw(12)<<std::scientific<<std::setprecision(5)<<Res<<"  Error: "<<std::setw(12)<<std::scientific<<std::setprecision(5)<<Err<<"  Acc: "<<resetiosflags(std::ios::scientific)<<std::fixed<<std::setw(6)<<std::setprecision(3)<<acc*100<<"%  Chi2: "<<Chi2<<std::endl;
        std::cout<<"======================================================================="<<std::endl;
        std::cout<<resetiosflags(std::ios::fixed)<<std::setprecision(8);
    }
}
double VEGAS_Integrator::Get_Result()
{
    double res_num = 0;
    double res_den = 0;
    for (int i = 0; i < Results.size(); i++)
    {
        res_num += Results[i]/Sigma2[i];
        res_den += 1.0/Sigma2[i];
    }
    return res_num/res_den;
}
double VEGAS_Integrator::Get_Error()
{
    double res = 0;
    for (int i = 0; i < Sigma2.size(); i++)
    {
        res += 1.0/Sigma2[i];
    }
    return 1.0/sqrt(res);
}
double VEGAS_Integrator::Get_Chisq()
{
    double Ifinal = Get_Result();
    double chi2 = 0;
    for (int i = 0; i < Results.size(); i++)
    {
        chi2 += pow(Results[i]-Ifinal,2)/Sigma2[i];
    }
    return chi2;
}