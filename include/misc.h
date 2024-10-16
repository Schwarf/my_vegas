//
// Created by andreas on 16.10.24.
//

#ifndef ABS_VEGAS_MISC_H
#define ABS_VEGAS_MISC_H
#include <iostream>
#include <limits>


constexpr double pow_constexpr(double base, unsigned int exp) {
    double result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result *= base;
        }
        base *= base;
        exp /= 2;
    }
    return result;
}

constexpr double nth_root_constexpr(double value, unsigned int root, double tolerance ) {
    double low = 1, high = value, mid, mid_pow;
    while (low < high) {
        mid = low + (high - low) / 2;
        mid_pow = pow_constexpr(mid, root);

        if (mid_pow == value) {
            return mid;
        } else if (mid_pow < value) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return high - 1;
}

template<unsigned long long N_STRAT, unsigned int NumberOfDimensions>
struct HyperCubes {
    static constexpr unsigned long long maximum_number_of_hyper_cubes = 10000;

    static constexpr unsigned long long initial_value = pow_constexpr(N_STRAT, NumberOfDimensions);
    static constexpr bool needs_adjustment = initial_value > maximum_number_of_hyper_cubes || NumberOfDimensions > 9;

    static constexpr unsigned long long adjusted_N_STRAT = needs_adjustment ? nth_root_constexpr(maximum_number_of_hyper_cubes, NumberOfDimensions) : N_STRAT;
    static constexpr unsigned long long number_of_hyper_cubes = pow_constexpr(adjusted_N_STRAT, NumberOfDimensions);

};

//int main() {
//    std::cout << "Number of hyper cubes: " << HyperCubes<2, 10>::number_of_hyper_cubes << std::endl;
//    return 0;
//}


#endif //ABS_VEGAS_MISC_H
