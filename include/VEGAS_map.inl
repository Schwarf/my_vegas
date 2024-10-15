#include "VEGAS_map.h"
#include <cmath>
#include <numeric>
#include <iostream>

template<int NumberOfDimensions, int NumberOfIntervals>
VegasMap<NumberOfDimensions, NumberOfIntervals>::VegasMap()
        : ID{NumberOfDimensions, 0}, // Initialize ID with zero
          x_edges{NumberOfDimensions, std::vector<double>(NumberOfIntervals + 1, 0.0)},
          dx_steps{NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0.0)},
          x_edges_last{NumberOfDimensions, std::vector<double>(NumberOfIntervals + 1, 0.0)},
          dx_steps_last{NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0.0)},
          weights{NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0.0)},
          counts{NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0.0)},
          smoothed_weights{NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0.0)},
          summed_weights{NumberOfDimensions, 0.0},
          delta_weights{NumberOfDimensions, 0.0},
          average_weight{NumberOfDimensions, 0.0},
          std_weight{NumberOfDimensions, 0.0},
          alpha{0.5} // Initialize alpha
{
    reset_map();
}


template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::reset_map() {
    double step_tmp = 1.0 / NumberOfIntervals;

    // Recalculate edges and steps
    for (int i = 0; i <= NumberOfIntervals; i++) {
        double current_x_edge = i * step_tmp;
        if (i > 0) {
            double prev_x_edge = (i - 1) * step_tmp;
            for (int d = 0; d < NumberOfDimensions; d++) {
                dx_steps[d][i - 1] = current_x_edge - prev_x_edge; // Set dx_steps
            }
        }
        for (int d = 0; d < NumberOfDimensions; d++) {
            x_edges[d][i] = current_x_edge; // Set x_edges
            x_edges_last[d][i] = current_x_edge; // Set x_edges_last
        }
    }

    // Reset weights and related statistics
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        std::fill(weights[dimension].begin(), weights[dimension].end(), 0.0);
        std::fill(counts[dimension].begin(), counts[dimension].end(), 0.0);
        std::fill(smoothed_weights[dimension].begin(), smoothed_weights[dimension].end(), 0.0);
        summed_weights[dimension] = 0.0;
        delta_weights[dimension] = 0.0;
        average_weight[dimension] = 0.0;
        std_weight[dimension] = 0.0;
    }
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
    for (int i = 0; i < NumberOfDimensions; i++) {
        ID[i] = std::floor(random_numbers[i] * NumberOfIntervals);
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
std::vector<double>
VegasMap<NumberOfDimensions, NumberOfIntervals>::get_interval_offset(const std::vector<double> &random_numbers) const {
//    auto ID = compute_interval_ID(random_numbers);
    std::vector<double> res(NumberOfDimensions);
    for (int i = 0; i < NumberOfDimensions; i++) {
        res[i] = random_numbers[i] * NumberOfIntervals - ID[i];
    }
    return res; // Profiling spot compare with move semantics
}

template<int NumberOfDimensions, int NumberOfIntervals>
std::vector<double> VegasMap<NumberOfDimensions, NumberOfIntervals>::get_x(const std::vector<double> &random_numbers) {
//    auto ID = compute_interval_ID(random_numbers);
    compute_interval_ID(random_numbers);
    auto offset = get_interval_offset(random_numbers);
    std::vector<double> res(NumberOfDimensions);
    for (int i = 0; i < NumberOfDimensions; i++) {
        int id = ID[i];
        res[i] = x_edges[i][id] + dx_steps[i][id] * offset[i];
    }
    return res; // Profiling spot compare with move semantics
}

template<int NumberOfDimensions, int NumberOfIntervals>
double VegasMap<NumberOfDimensions, NumberOfIntervals>::get_jacobian(const std::vector<double> &random_numbers) {
//    auto ID = compute_interval_ID(random_numbers);
    double jacobian{1.0};
    for (int i = 0; i < NumberOfDimensions; i++) {
        int id = ID[i];
        jacobian *= NumberOfIntervals * dx_steps[i][id];
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
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::smooth_weight() {
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        for (int intervals{}; intervals < weights[dimension].size(); intervals++) {
            if (counts[dimension][intervals] != 0) {
                weights[dimension][intervals] /= counts[dimension][intervals];
            }
        }
    }
    // std::cout<<"Count devided!"<<std::endl;
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        double d_tmp{};
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
        for (int i = 1; i < NumberOfIntervals - 1; i++) {
            d_tmp = (weights[dimension][i - 1] + 6.0 * weights[dimension][i] + weights[dimension][i + 1]) / (8.0 * d_sum);
            if (d_tmp != 0.0) {
                d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
            }
            smoothed_weights[dimension][i] = d_tmp;
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
    for (int dimensions{}; dimensions < NumberOfDimensions; dimensions++) {
        int current_old_interval{};
        int current_new_interval{1};
        double d_accu = 0;
        while (true) {
            d_accu += delta_weights[dimensions];
            while (d_accu > smoothed_weights[dimensions][current_old_interval]) {
                d_accu -= smoothed_weights[dimensions][current_old_interval];
                current_old_interval++;
            }
            x_edges[dimensions][current_new_interval] = x_edges_last[dimensions][current_old_interval] +
                                                        d_accu / smoothed_weights[dimensions][current_old_interval] *
                                                        dx_steps_last[dimensions][current_old_interval];
            dx_steps[dimensions][current_new_interval - 1] =
                    x_edges[dimensions][current_new_interval] - x_edges[dimensions][current_new_interval - 1];
            current_new_interval++;
            if (current_new_interval >= NumberOfIntervals) {
                break;
            }
        }
        dx_steps[dimensions][NumberOfIntervals - 1] =
                x_edges[dimensions][(NumberOfIntervals + 1) - 1] - x_edges[dimensions][(NumberOfIntervals + 1) - 2];
    }
    reset_weight();
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::check_weight() {
    for (int dimensions = 0; dimensions < NumberOfDimensions; dimensions++) {
        average_weight[dimensions] = 0;
        for (int interval{}; interval < NumberOfIntervals; ++interval) {
            average_weight[dimensions] += weights[dimensions][interval];
        }
        average_weight[dimensions] /= static_cast<double>(weights[dimensions].size());
        for (int interval{}; interval < NumberOfIntervals; ++interval) {
            std_weight[dimensions] +=
                    (weights[dimensions][interval] - average_weight[dimensions]) * (weights[dimensions][interval] - average_weight[dimensions]);
        }
        std_weight[dimensions] = sqrt(std_weight[dimensions]);// /average_weight;
    }
}

template<int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::print_edges() {
    std::cout << "Grid Map:" << std::endl;
    for (int dimension{}; dimension < NumberOfDimensions; ++dimension) {
        std::cout << "\tx_" << dimension << ":";
        for (int edge{}; edge < (NumberOfIntervals + 1); ++edge) {
            std::cout << "\t" << x_edges[dimension][edge];
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
        for (int interval{}; interval < NumberOfIntervals; ++interval) {
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
        for (int edge{}; edge < (NumberOfIntervals + 1); ++edge) {
            chi2 += (x_edges[dimension][edge] - x_edges_last[dimension][edge]) * (x_edges[dimension][edge] - x_edges_last[dimension][edge]) / dx_ave /
                    dx_ave;
        }
    }
    return chi2 / NumberOfDimensions / (NumberOfIntervals + 1);
}