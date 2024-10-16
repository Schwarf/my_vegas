#ifndef VEGAS_STRATIFY_H
#define VEGAS_STRATIFY_H

#include <vector>
#include <cmath>
#include "misc.h"

template<int NumberOfDimensions>
class VEGAS_Stratify {
public:
    static constexpr int NumberOfStratifications{10};
    static constexpr double MaximumAllowedNumberOfHypercubes{10000.0};
    static constexpr double RealNumberOfHyperCubes = pow_constexpr(static_cast<double>(NumberOfStratifications), NumberOfDimensions);


    VEGAS_Stratify() : beta{0.75} {
        // NumberOfStratifications = floor(pow(N_EVALUATES_TRAINED/4.0,1.0/NumberOfDimensions));

        number_of_hyper_cubes = pow(NumberOfStratifications, NumberOfDimensions);
        // if NumberOfDimensions too large, number_of_hyper_cubes will exceed the MAXIMUM number an integer can store
        if (number_of_hyper_cubes > MaximumAllowedNumberOfHypercubes || NumberOfDimensions > 9) {
            NumberOfStratifications = floor(pow(MaximumAllowedNumberOfHypercubes, 1.0 / NumberOfDimensions));
            number_of_hyper_cubes = pow(NumberOfStratifications, NumberOfDimensions);
        }


        V_cubic = pow(1.0 / NumberOfStratifications, NumberOfDimensions);

        squared_accumulated_function_values = std::vector<double>(number_of_hyper_cubes, 0);
        accumulated_function_values = std::vector<double>(number_of_hyper_cubes, 0);
        counts = std::vector<double>(number_of_hyper_cubes, 0);
        hypercubic_weights = std::vector<double>(number_of_hyper_cubes, 1.0 / number_of_hyper_cubes);
    };

    ~VEGAS_Stratify() = default;

public:
    void Set_NEVAL(int NEVAL_EXP) {
        number_of_expected_evaluations = NEVAL_EXP;
    }

    std::array<int, NumberOfDimensions> get_indices(int index) {
        std::array<int, NumberOfDimensions> indices{};
        for (int i = 0; i < NumberOfDimensions; i++) {
            int quotient = index / NumberOfStratifications;
            int remainder = index - quotient * NumberOfStratifications;
            indices[i] = remainder;
            index = quotient;
        }
        return indices;
    }

    std::array<double, NumberOfDimensions> get_y(int index, const std::array<double, NumberOfDimensions> &random_uni) {
        const double dy = 1.0 / NumberOfStratifications;
        std::array<double, NumberOfDimensions> res{};
        std::array<int, NumberOfDimensions> ID = get_indices(index);
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

private:
    double beta{};
    double V_cubic{};
    std::vector<double> squared_accumulated_function_values; // size = (NumberOfStratifications)^(NumberOfDimensions)
    std::vector<double> accumulated_function_values; // size = (NumberOfStratifications)^(NumberOfDimensions)
    std::vector<double> counts; // size = (NumberOfStratifications)^(NumberOfDimensions)
    std::vector<double> hypercubic_weights; // size = (NumberOfStratifications)^(NumberOfDimensions)
    // int N_EVALUATES_TRAINED; // The evaluates number used to train the stratification
    int number_of_expected_evaluations{};
    int number_of_hyper_cubes{};

};



#endif //VEGAS_STRATIFY_H