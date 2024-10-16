//
// Created by andreas on 15.10.24.
//
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iostream>
int main() {
//    const int NumberOfDimensions = 1000000;
//    std::vector<double> random_numbers(NumberOfDimensions);
//    std::vector<double> ID(NumberOfDimensions);
//    double NumberOfIntervals = 10.0;
//
//    // Fill random_numbers with random values
//    std::generate(random_numbers.begin(), random_numbers.end(), []() { return rand() / (double)RAND_MAX; });
//
//    auto start = std::chrono::high_resolution_clock::now();
//
//    // for loop version
//    for (int i = 0; i < NumberOfDimensions; i++) {
//        ID[i] = std::floor(random_numbers[i] * NumberOfIntervals);
//    }
//    // std::transform version
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> diff = end - start;
//    std::cout << "Transform: " << diff.count() << " s\n";
//
//    start = std::chrono::high_resolution_clock::now();
//    std::transform(random_numbers.begin(), random_numbers.end(), ID.begin(),
//                   [NumberOfIntervals](const double random_number) { return std::floor(random_number * NumberOfIntervals); });
//
////    }
//
//    end = std::chrono::high_resolution_clock::now();
//    diff = end - start;
//    std::cout << "For loop: " << diff.count() << " s\n";
    constexpr int N_STRAT{10};
    constexpr int NumberOfDimensions{6};
    auto start = std::chrono::high_resolution_clock::now();
    const double inverse{1.0/N_STRAT};
    double V_cubic = 1.0;
    for (int i{}; i < NumberOfDimensions; ++i)
    {
        V_cubic *= inverse;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "for-loop: " << diff.count() << " s\n";

    start = std::chrono::high_resolution_clock::now();
    const double inverse2{1.0/N_STRAT};
    double V_cubic2 = pow(inverse2, NumberOfDimensions);
    end = std::chrono::high_resolution_clock::now();
    diff =  end-start;
    std::cout << "pow: " << diff.count() << " s\n";
    return 0;
}
