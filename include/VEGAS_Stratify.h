#ifndef VEGAS_STRATIFY_H
#define VEGAS_STRATIFY_H

#include <vector>

#include <cmath>
class VEGAS_Stratify
{
private:
    int N_DIM;
    int N_STRAT;
    double beta;
    double V_cubic;
    std::vector<double> JF2; // size = (N_STRAT)^(number_of_dimensions)
    std::vector<double> JF; // size = (N_STRAT)^(number_of_dimensions)
    std::vector<double> Counts; // size = (N_STRAT)^(number_of_dimensions)
    std::vector<double> dh; // size = (N_STRAT)^(number_of_dimensions)
    // int N_EVALUATES_TRAINED; // The evaluates number used to train the stratification
    int N_EVALUATES_EXPECTED;
    int N_HYPERCUBICS;
    int N_HYPERCUBICS_MAX;


public:
    void Set_Dimension(int ndim)
    {
        N_DIM = ndim;
        Reset_Storage();
    }
    void Set_NEVAL(int NEVAL_EXP)
    {
        N_EVALUATES_EXPECTED = NEVAL_EXP;
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

        N_HYPERCUBICS = pow(N_STRAT,N_DIM);
        if (N_HYPERCUBICS > N_HYPERCUBICS_MAX || N_DIM > 9) // if number_of_dimensions too large, N_HYPERCUBICS will exceed the MAXIMUM number an integer can store
        {
            N_STRAT = floor(pow(N_HYPERCUBICS_MAX,1.0/N_DIM));
            N_HYPERCUBICS = pow(N_STRAT,N_DIM);
        }


        V_cubic = pow(1.0/N_STRAT, N_DIM);
        JF2 = std::vector<double>(N_HYPERCUBICS,0);
        JF  = std::vector<double>(N_HYPERCUBICS,0);
        Counts = std::vector<double>(N_HYPERCUBICS,0);
        dh  =std:: vector<double>(N_HYPERCUBICS,1.0/N_HYPERCUBICS);
    }
    std::vector<int> Get_Indices(int index)
    {
        std::vector<int> res(N_DIM,0);
        int Quotient;
        int tmp = index;
        int Remainder;
        for (int i = 0; i < N_DIM; i++)
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
        std::vector<double> res(N_DIM,0);
        std::vector<int> ID = Get_Indices(index);
        for (int i = 0; i < N_DIM; i++)
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
        for (int i = 0; i < N_HYPERCUBICS; i++)
        {
            d_tmp = V_cubic*V_cubic/Counts[i]*JF2[i] - std::pow(V_cubic/Counts[i]*JF[i],2);
            dh[i] = std::pow(d_tmp,beta);
            d_sum += dh[i];
        }
        for (int i = 0; i < N_HYPERCUBICS; i++)
        {
            dh[i] = dh[i]/d_sum;
        }
    }
    int Get_NH(int index) // Get the expected number of events in each hypercubic.
    {
        int nh = dh[index]*N_EVALUATES_EXPECTED;
        return nh<2?2:nh;
    }

    VEGAS_Stratify(){N_DIM = 1; N_STRAT = 10; beta = 0.75; N_HYPERCUBICS_MAX = 10000;};
    ~VEGAS_Stratify(){};
    double Get_V_Cubic(){return V_cubic;}
    int Get_NHYPERCUBICS(){return N_HYPERCUBICS;};
    // void Set_Stratification_System(int number_of_dimensions, int NEVAL_TRAIN);
};


#endif //VEGAS_STRATIFY_H