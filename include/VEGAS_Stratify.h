#ifndef VEGAS_STRATIFY_H
#define VEGAS_STRATIFY_H

#include <vector>
#include <cmath>
#include "misc.h"


template<int NumberOfDimensions>
class VEGAS_Stratify {
private:
    static constexpr int InitialNumberOfStratifications{10};
    static constexpr int MaximumNumberOfHypercubes{10000};
    static constexpr long long unsigned InitialNumberOfHyperCubes = pow_constexpr(InitialNumberOfStratifications, NumberOfDimensions);
    static constexpr long long NumberOfStratifications = (InitialNumberOfHyperCubes > MaximumNumberOfHypercubes || NumberOfDimensions > 9) ? ValueMapper<MaximumNumberOfHypercubes, NumberOfDimensions>::NumberOfStratifications : InitialNumberOfStratifications;
    static constexpr long long NumberOfHyperCubes = (InitialNumberOfHyperCubes > MaximumNumberOfHypercubes || NumberOfDimensions > 9) ? ValueMapper<MaximumNumberOfHypercubes, NumberOfDimensions>::NumberOfHyperCubes : InitialNumberOfHyperCubes;

    double beta{};
    double V_cubic{};
    std::array<double, NumberOfHyperCubes> squared_accumulated_function_values;
    std::array<double, NumberOfHyperCubes> accumulated_function_values;
    std::array<double, NumberOfHyperCubes> counts;
    std::array<double, NumberOfHyperCubes> hypercubic_weights;
    int number_of_expected_evaluations{};

public:
    VEGAS_Stratify() : beta{0.75}{
        V_cubic = std::pow(1.0 / NumberOfStratifications, NumberOfDimensions);
        squared_accumulated_function_values = {};
        accumulated_function_values = {};
        counts = {};
        hypercubic_weights.fill(1.0 / NumberOfHyperCubes);
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
        for (int i = 0; i < NumberOfHyperCubes; i++) {
            d_tmp = V_cubic * V_cubic / counts[i] * squared_accumulated_function_values[i] -
                    std::pow(V_cubic / counts[i] * accumulated_function_values[i], 2);
            hypercubic_weights[i] = std::pow(d_tmp, beta);
            d_sum += hypercubic_weights[i];
        }
        for (int i = 0; i < NumberOfHyperCubes; i++) {
            hypercubic_weights[i] = hypercubic_weights[i] / d_sum;
        }
    }

    int get_expected_events_per_hyper_cube(int index) // Get the expected number of events in each hypercubic.
    {
        int nh = hypercubic_weights[index] * number_of_expected_evaluations;
        return nh < 2 ? 2 : nh;
    }


    double get_V_cubic() { return V_cubic; }

    int Get_NHYPERCUBICS() { return NumberOfHyperCubes; };
    // void Set_Stratification_System(int NumberOfDimensions, int NEVAL_TRAIN);
};


#endif //VEGAS_STRATIFY_H