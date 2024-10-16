#ifndef VEGAS_STRATIFY_H
#define VEGAS_STRATIFY_H

#include <vector>
#include <cmath>
#include "misc.h"
template<int NumberOfDimensions>
class VEGAS_Stratify {
public:
    VEGAS_Stratify() : N_STRAT{10}, beta{0.75}, maximum_number_of_hyper_cubes{10000} {
        // N_STRAT = floor(pow(N_EVALUATES_TRAINED/4.0,1.0/NumberOfDimensions));

        number_of_hyper_cubes = pow(N_STRAT, NumberOfDimensions);
        // if NumberOfDimensions too large, number_of_hyper_cubes will exceed the MAXIMUM number an integer can store
        if (number_of_hyper_cubes > maximum_number_of_hyper_cubes || NumberOfDimensions > 9)
        {
            N_STRAT = floor(pow(maximum_number_of_hyper_cubes, 1.0 / NumberOfDimensions));
            number_of_hyper_cubes = pow(N_STRAT, NumberOfDimensions);
        }


        V_cubic = pow(1.0 / N_STRAT, NumberOfDimensions);

        squared_accumulated_function_values = std::vector<double>(number_of_hyper_cubes, 0);
        accumulated_function_values = std::vector<double>(number_of_hyper_cubes, 0);
        counts = std::vector<double>(number_of_hyper_cubes, 0);
        hypercubic_weights = std::vector<double>(number_of_hyper_cubes, 1.0 / number_of_hyper_cubes);
    };

    ~VEGAS_Stratify() = default;

private:
    int N_STRAT;
    double beta{};
    double V_cubic{};
    std::vector<double> squared_accumulated_function_values; // size = (N_STRAT)^(NumberOfDimensions)
    std::vector<double> accumulated_function_values; // size = (N_STRAT)^(NumberOfDimensions)
    std::vector<double> counts; // size = (N_STRAT)^(NumberOfDimensions)
    std::vector<double> hypercubic_weights; // size = (N_STRAT)^(NumberOfDimensions)
    // int N_EVALUATES_TRAINED; // The evaluates number used to train the stratification
    int number_of_expected_evaluations{};
    int number_of_hyper_cubes{};
    int maximum_number_of_hyper_cubes{};


public:
    void Set_NEVAL(int NEVAL_EXP) {
        number_of_expected_evaluations = NEVAL_EXP;
    }

    std::vector<int> get_indices(int index) {
        std::vector<int> indices(NumberOfDimensions, 0);
        for (int i = 0; i < NumberOfDimensions; i++) {
            int quotient = index / N_STRAT;
            int remainder = index - quotient * N_STRAT;
            indices[i] = remainder;
            index = quotient;
        }
        return indices;
    }

    std::vector<double> get_y(int index, std::vector<double> random_uni) {
        const double dy = 1.0 / N_STRAT;
        std::vector<double> res(NumberOfDimensions, 0);
        std::vector<int> ID = get_indices(index);
        for (int i = 0; i < NumberOfDimensions; i++) {
            res[i] = (random_uni[i] + ID[i]) * dy;
        }
        return res;
    }

    void accumulate_weights(int cube_index, double function_value_times_jacobian) {
        // This function_value_times_jacobian is J*f;
        squared_accumulated_function_values[cube_index] +=
                function_value_times_jacobian * function_value_times_jacobian;
        accumulated_function_values[cube_index] += function_value_times_jacobian;
        counts[cube_index] += 1;
    }

    void update_hypercubic_weights() {
        double d_sum = 0;
        double d_tmp;
        for (int i = 0; i < number_of_hyper_cubes; i++) {
            d_tmp = V_cubic * V_cubic / counts[i] * squared_accumulated_function_values[i] -
                    std::pow(V_cubic / counts[i] * accumulated_function_values[i], 2);
            hypercubic_weights[i] = std::pow(d_tmp, beta);
            d_sum += hypercubic_weights[i];
        }
        for (int i = 0; i < number_of_hyper_cubes; i++) {
            hypercubic_weights[i] = hypercubic_weights[i] / d_sum;
        }
    }

    int get_expected_events_per_hyper_cube(int index) // Get the expected number of events in each hypercubic.
    {
        int nh = hypercubic_weights[index] * number_of_expected_evaluations;
        return nh < 2 ? 2 : nh;
    }


    double get_V_cubic() { return V_cubic; }

    int Get_NHYPERCUBICS() { return number_of_hyper_cubes; };
    // void Set_Stratification_System(int NumberOfDimensions, int NEVAL_TRAIN);
};


#endif //VEGAS_STRATIFY_H