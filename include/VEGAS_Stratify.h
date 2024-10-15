#ifndef VEGAS_STRATIFY_H
#define VEGAS_STRATIFY_H

#include <vector>
#include <cmath>

class VEGAS_Stratify
{
private:
    int number_of_dimensions{};
    int N_STRAT;
    double beta{};
    double V_cubic{};
    std::vector<double> JF2; // size = (N_STRAT)^(number_of_dimensions)
    std::vector<double> JF; // size = (N_STRAT)^(number_of_dimensions)
    std::vector<double> Counts; // size = (N_STRAT)^(number_of_dimensions)
    std::vector<double> dh; // size = (N_STRAT)^(number_of_dimensions)
    // int N_EVALUATES_TRAINED; // The evaluates number used to train the stratification
    int number_of_expected_evaluations{};
    int number_of_hyper_cubes{};
    int maximum_number_of_hyper_cubes{};


public:
    void Set_Dimension(int ndim)
    {
        number_of_dimensions = ndim;
        Reset_Storage();
    }
    void Set_NEVAL(int NEVAL_EXP)
    {
        number_of_expected_evaluations = NEVAL_EXP;
    }
// void Set_Stratification_System(int ndim, int NEVAL_TRAIN)
// {
//     number_of_dimensions = ndim;
//     N_EVALUATES_TRAINED = NEVAL_TRAIN;
//     Reset_Storage();
// }
    void Reset_Storage()
    {
        // N_STRAT = floor(pow(N_EVALUATES_TRAINED/4.0,1.0/number_of_dimensions));

        number_of_hyper_cubes = pow(N_STRAT, number_of_dimensions);
        if (number_of_hyper_cubes > maximum_number_of_hyper_cubes || number_of_dimensions > 9) // if number_of_dimensions too large, number_of_hyper_cubes will exceed the MAXIMUM number an integer can store
        {
            N_STRAT = floor(pow(maximum_number_of_hyper_cubes, 1.0 / number_of_dimensions));
            number_of_hyper_cubes = pow(N_STRAT, number_of_dimensions);
        }


        V_cubic = pow(1.0/N_STRAT, number_of_dimensions);
        JF2 = std::vector<double>(number_of_hyper_cubes, 0);
        JF  = std::vector<double>(number_of_hyper_cubes, 0);
        Counts = std::vector<double>(number_of_hyper_cubes, 0);
        dh  =std:: vector<double>(number_of_hyper_cubes, 1.0 / number_of_hyper_cubes);
    }
    std::vector<int> Get_Indices(int index)
    {
        std::vector<int> res(number_of_dimensions, 0);
        int Quotient;
        int tmp = index;
        int Remainder;
        for (int i = 0; i < number_of_dimensions; i++)
        {
            Quotient = tmp/N_STRAT;
            Remainder = tmp - Quotient*N_STRAT;
            res[i] = Remainder;
            tmp = Quotient;
        }
        return res;
    }

    std::vector<double> Get_Y(int index, std::vector<double> random_uni)
    {
        double dy = 1.0/N_STRAT;
        std::vector<double> res(number_of_dimensions, 0);
        std::vector<int> ID = Get_Indices(index);
        for (int i = 0; i < number_of_dimensions; i++)
        {
            res[i] = random_uni[i]*dy + ID[i]*dy;
        }
        return res;
    }
    void Accumulate_Weight(int index, double weight)
    {
        // This weight is J*f;
        JF2[index] += weight*weight;
        JF[index] += weight;
        Counts[index] += 1;
    }
    void Update_DH()
    {
        double d_sum = 0;
        double d_tmp;
        for (int i = 0; i < number_of_hyper_cubes; i++)
        {
            d_tmp = V_cubic*V_cubic/Counts[i]*JF2[i] - std::pow(V_cubic/Counts[i]*JF[i],2);
            dh[i] = std::pow(d_tmp,beta);
            d_sum += dh[i];
        }
        for (int i = 0; i < number_of_hyper_cubes; i++)
        {
            dh[i] = dh[i]/d_sum;
        }
    }
    int Get_NH(int index) // Get the expected number of events in each hypercubic.
    {
        int nh = dh[index] * number_of_expected_evaluations;
        return nh<2?2:nh;
    }

    VEGAS_Stratify(): number_of_dimensions{1}, N_STRAT{10}, beta{0.75}, maximum_number_of_hyper_cubes{10000}{
        Reset_Storage();
    };
    ~VEGAS_Stratify(){};
    double Get_V_Cubic(){return V_cubic;}
    int Get_NHYPERCUBICS(){return number_of_hyper_cubes;};
    // void Set_Stratification_System(int number_of_dimensions, int NEVAL_TRAIN);
};


#endif //VEGAS_STRATIFY_H