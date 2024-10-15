#include "VEGAS_map.h"
#include <cmath>
#include <numeric>
#include <iostream>

template<int NumberOfDimensions, int NumberOfIntervals>
VegasMap<NumberOfDimensions, NumberOfIntervals>::VegasMap() {
    ID = std::vector<int>(NumberOfDimensions);
    alpha = 0.5;
    constexpr double step_tmp = 1.0 / NumberOfIntervals;
    std::vector<double> x_edges_tmp((NumberOfIntervals + 1));
    std::vector<double> dx_steps_tmp((NumberOfIntervals + 1) - 1, step_tmp);
    for (int i = 1; i < (NumberOfIntervals + 1); i++) {
        x_edges_tmp[i] = i * step_tmp;
    }
    x_edges = std::vector<std::vector<double> >(NumberOfDimensions, x_edges_tmp);
    dx_steps = std::vector<std::vector<double> >(NumberOfDimensions, dx_steps_tmp);
    x_edges_last = std::vector<std::vector<double> >(NumberOfDimensions, x_edges_tmp);
    dx_steps_last = std::vector<std::vector<double> >(NumberOfDimensions, dx_steps_tmp);
    weights = std::vector<std::vector<double> >(NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0.0));
    counts = std::vector<std::vector<double> >(NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0.0));
    smoothed_weights = std::vector<std::vector<double> >(NumberOfDimensions,
                                                         std::vector<double>(NumberOfIntervals, 0.0));
    summed_weights = std::vector<double>(NumberOfDimensions, 0.0);
    delta_weights = std::vector<double>(NumberOfDimensions, 0.0);
    average_weight = std::vector<double>(NumberOfDimensions, 0.0);
    std_weight = std::vector<double>(NumberOfDimensions, 0.0);
}



template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::reset_weight() {
    for (auto &row: weights) {
        std::fill(row.begin(), row.end(), 0.0);
    }

    for (auto &row: counts) {
        std::fill(row.begin(), row.end(), 0.0);
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::compute_interval_ID(const std::vector<double> &random_numbers) {
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        ID[dimension] = std::floor(random_numbers[dimension] * NumberOfIntervals);
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
std::vector<double>
VegasMap<NumberOfDimensions, NumberOfIntervals>::get_interval_offset(const std::vector<double> &random_numbers) const {
//    auto ID = compute_interval_ID(random_numbers);
    std::vector<double> res(NumberOfDimensions);
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        res[dimension] = random_numbers[dimension] * NumberOfIntervals - ID[dimension];
    }
    return res; // Profiling spot compare with move semantics
}

template<int NumberOfDimensions, int NumberOfIntervals>
std::vector<double> VegasMap<NumberOfDimensions, NumberOfIntervals>::get_x(const std::vector<double> &random_numbers) {
//    auto ID = compute_interval_ID(random_numbers);
    compute_interval_ID(random_numbers);
    auto offset = get_interval_offset(random_numbers);
    std::vector<double> res(NumberOfDimensions);
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        int id = ID[dimension];
        res[dimension] = x_edges[dimension][id] + dx_steps[dimension][id] * offset[dimension];
    }
    return res; // Profiling spot compare with move semantics
}

template<int NumberOfDimensions, int NumberOfIntervals>
double VegasMap<NumberOfDimensions, NumberOfIntervals>::get_jacobian(const std::vector<double> &random_numbers) {
//    auto ID = compute_interval_ID(random_numbers);
    double jacobian{1.0};
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        int id = ID[dimension];
        jacobian *= NumberOfIntervals * dx_steps[dimension][id];
    }
    return jacobian;
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::accumulate_weight(const std::vector<double> &y, double f) {
    // f is the value of integrand!
//    auto ID = compute_interval_ID(y);
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        int id = ID[dimension];
        weights[dimension][id] += (f * get_jacobian(y)) * (f * get_jacobian(y));
        counts[dimension][id] += 1;
        // std::cout<<"ID: "<<id<<" weight: "<<weights[dimension][id]<<" counts: "<<counts[dimension][id]<<std::endl;
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::smooth_weight() {
    // std::cout<<"Smoothing weight"<<std::endl;
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        for (int interval{}; interval < NumberOfIntervals; ++interval) {
            if (counts[dimension][interval] != 0) {
                weights[dimension][interval] /= counts[dimension][interval];
            }
        }
    }
    // std::cout<<"Count devided!"<<std::endl;
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        double d_tmp;
        double d_sum = std::accumulate(weights[dimension].begin(), weights[dimension].end(), 0.0);
        summed_weights[dimension] = 0;

        // Handle the first interval (i == 0) outside the loop
        d_tmp = (7.0 * weights[dimension][0] + weights[dimension][1]) / (8.0 * d_sum);
        if (d_tmp != 0.0) {
            d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
        }
        smoothed_weights[dimension][0] = d_tmp;
        summed_weights[dimension] += d_tmp;

        // Handle the main loop (1 <= i < NumberOfIntervals - 1)
        for (int interval = 1; interval < NumberOfIntervals - 1; ++interval) {
            d_tmp = (weights[dimension][interval - 1] + 6.0 * weights[dimension][interval] + weights[dimension][interval + 1]) / (8.0 * d_sum);
            if (d_tmp != 0.0) {
                d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
            }
            smoothed_weights[dimension][interval] = d_tmp;
            summed_weights[dimension] += d_tmp;
        }
        // Handle the last interval (i == NumberOfIntervals - 1) outside the loop
        d_tmp = (weights[dimension][NumberOfIntervals - 2] + 7.0 * weights[dimension][NumberOfIntervals - 1]) /
                (8.0 * d_sum);
        if (d_tmp != 0.0) {
            d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
        }

        smoothed_weights[dimension][NumberOfIntervals - 1] = d_tmp;
        summed_weights[dimension] += d_tmp;

        delta_weights[dimension] = summed_weights[dimension] / NumberOfIntervals;

    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::update_map() {
    smooth_weight();
    x_edges_last = x_edges;
    dx_steps_last = dx_steps;
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        int current_old_interval{};
        int current_new_interval{1};
        double d_accu = 0;
        while (true) {
            d_accu += delta_weights[dimension];
            while (d_accu > smoothed_weights[dimension][current_old_interval]) {
                d_accu -= smoothed_weights[dimension][current_old_interval];
                current_old_interval++;
            }
            x_edges[dimension][current_new_interval] = x_edges_last[dimension][current_old_interval] +
                                                       d_accu / smoothed_weights[dimension][current_old_interval] *
                                                       dx_steps_last[dimension][current_old_interval];
            dx_steps[dimension][current_new_interval - 1] =
                    x_edges[dimension][current_new_interval] - x_edges[dimension][current_new_interval - 1];
            current_new_interval++;
            if (current_new_interval >= NumberOfIntervals) {
                break;
            }
        }
        dx_steps[dimension][NumberOfIntervals - 1] =
                x_edges[dimension][NumberOfIntervals] - x_edges[dimension][NumberOfIntervals - 1];
    }
    reset_weight();
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::check_weight() {
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        average_weight[dimension] = 0;
        for (int interval{}; interval < weights[dimension].size(); ++interval) {
            average_weight[dimension] += weights[dimension][interval];
        }
        average_weight[dimension] /= static_cast<double>(weights[dimension].size());
        for (int interval{}; interval < weights[dimension].size(); ++interval) {
            std_weight[dimension] += (weights[dimension][interval] - average_weight[dimension]) *
                    (weights[dimension][interval] - average_weight[dimension]);
        }
        std_weight[dimension] = sqrt(std_weight[dimension]);// /average_weight;
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::print_edges() {
    std::cout << "Grid Map:" << std::endl;
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        std::cout << "\tx_" << dimension << ":";
        for (int interval{}; interval < (NumberOfIntervals + 1); ++interval) {
            std::cout << "\t" << x_edges[dimension][interval];
        }
        std::cout << std::endl;
        std::cout << "\tdx_" << dimension << ":";
        for (int interval{}; interval < NumberOfIntervals; ++interval) {
            std::cout << "\t" << dx_steps[dimension][interval];
        }
        std::cout << std::endl;
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::print_weights() {
    std::cout << "Weights:" << std::endl;
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        std::cout << "\tweight_" << dimension << ":";
        for (int interval{}; interval < NumberOfIntervals; interval++) {
            std::cout << "\t" << weights[dimension][interval];
        }
        std::cout << std::endl;
        check_weight();
        std::cout << "\t\tAverage: " << average_weight[dimension] << "  STD: " << std_weight[dimension] << " std/ave: "
                  << std_weight[dimension] / average_weight[dimension] / NumberOfIntervals << std::endl;
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
double VegasMap<NumberOfDimensions, NumberOfIntervals>::checking_map() {
    double dx_ave = 1.0 / NumberOfIntervals;
    double chi2{};
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        for (int edges{}; edges < (NumberOfIntervals + 1); ++edges) {
            chi2 += (x_edges[dimension][edges] - x_edges_last[dimension][edges]) * (x_edges[dimension][edges] - x_edges_last[dimension][edges]) / dx_ave /
                    dx_ave;
        }
    }
    return chi2 / NumberOfDimensions /  (NumberOfIntervals + 1);
}