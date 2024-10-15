#include "VEGAS_map.h"
#include <cmath>
#include <numeric>
#include <iostream>

template <int NumberOfDimensions, int NumberOfIntervals>
VegasMap<NumberOfDimensions, NumberOfIntervals>::VegasMap() {
    number_of_edges = NumberOfIntervals + 1;
    ID = std::vector<int>(NumberOfDimensions);
    alpha = 0.5;
    reset_map();
}


template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::reset_map() {
    x_edges.clear();
    dx_steps.clear();
    double step_tmp = 1.0 / NumberOfIntervals;
    std::vector<double> x_edges_tmp(number_of_edges);
    std::vector<double> dx_steps_tmp(number_of_edges - 1);
    for (int i = 1; i < number_of_edges; i++) {
        x_edges_tmp[i] = i * step_tmp;
        dx_steps_tmp[i - 1] = x_edges_tmp[i] - x_edges_tmp[i - 1];
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

template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::reset_weight() {
    weights = std::vector<std::vector<double> >(NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0));
    counts = std::vector<std::vector<double> >(NumberOfDimensions, std::vector<double>(NumberOfIntervals, 0));
}

template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::compute_interval_ID(const std::vector<double> &random_numbers)  {

    for (int i = 0; i < NumberOfDimensions; i++) {
        ID[i] = std::floor(random_numbers[i] * NumberOfIntervals);
    }
}

template <int NumberOfDimensions, int NumberOfIntervals>
std::vector<double> VegasMap<NumberOfDimensions, NumberOfIntervals>::get_interval_offset(const std::vector<double> &random_numbers) const {
//    auto ID = compute_interval_ID(random_numbers);
    std::vector<double> res(NumberOfDimensions);
    for (int i = 0; i < NumberOfDimensions; i++) {
        res[i] = random_numbers[i] * NumberOfIntervals - ID[i];
    }
    return res; // Profiling spot compare with move semantics
}

template <int NumberOfDimensions, int NumberOfIntervals>
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

template <int NumberOfDimensions, int NumberOfIntervals>
double VegasMap<NumberOfDimensions, NumberOfIntervals>::get_jacobian(const std::vector<double> &random_numbers) {
//    auto ID = compute_interval_ID(random_numbers);
    double jacobian{1.0};
    for (int i = 0; i < NumberOfDimensions; i++) {
        int id = ID[i];
        jacobian *= NumberOfIntervals * dx_steps[i][id];
    }
    return jacobian;
}

template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::accumulate_weight(const std::vector<double> &y, double f) {
    // f is the value of integrand!
//    auto ID = compute_interval_ID(y);
    for (int i = 0; i < NumberOfDimensions; i++) {
        int id = ID[i];
        weights[i][id] += (f * get_jacobian(y))*(f * get_jacobian(y));
        counts[i][id] += 1;
        // std::cout<<"ID: "<<id<<" weight: "<<weights[i][id]<<" counts: "<<counts[i][id]<<std::endl;
    }
}

template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::smooth_weight() {
    // std::cout<<"Smoothing weight"<<std::endl;
    for (int i_dim = 0; i_dim < NumberOfDimensions; i_dim++) {
        for (int i_inter = 0; i_inter < weights[i_dim].size(); i_inter++) {
            if (counts[i_dim][i_inter] != 0) {
                weights[i_dim][i_inter] /= counts[i_dim][i_inter];
            }
        }
    }
    // std::cout<<"Count devided!"<<std::endl;
    for (int i_dim = 0; i_dim < NumberOfDimensions; i_dim++) {
        double d_tmp;
        double d_sum = std::accumulate(weights[i_dim].begin(), weights[i_dim].end(), 0.0);
        summed_weights[i_dim] = 0;

        // Handle the first interval (i == 0) outside the loop
        d_tmp = (7.0 * weights[i_dim][0] + weights[i_dim][1]) / (8.0 * d_sum);
        if (d_tmp != 0.0) {
            d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
        }
        smoothed_weights[i_dim][0] = d_tmp;
        summed_weights[i_dim] += d_tmp;

        // Handle the main loop (1 <= i < NumberOfIntervals - 1)
        for (int i = 1; i < NumberOfIntervals - 1; i++) {
            d_tmp = (weights[i_dim][i - 1] + 6.0 * weights[i_dim][i] + weights[i_dim][i + 1]) / (8.0 * d_sum);
            if (d_tmp != 0.0) {
                d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
            }
            smoothed_weights[i_dim][i] = d_tmp;
            summed_weights[i_dim] += d_tmp;
        }
        // Handle the last interval (i == NumberOfIntervals - 1) outside the loop
        d_tmp = (weights[i_dim][NumberOfIntervals - 2] + 7.0 * weights[i_dim][NumberOfIntervals - 1]) /
                (8.0 * d_sum);
        if (d_tmp != 0.0) {
            d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
        }

        smoothed_weights[i_dim][NumberOfIntervals - 1] = d_tmp;
        summed_weights[i_dim] += d_tmp;

        delta_weights[i_dim] = summed_weights[i_dim] / NumberOfIntervals;

    }
}

template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::update_map() {
    smooth_weight();
    x_edges_last = x_edges;
    dx_steps_last = dx_steps;
    for (int i_dim = 0; i_dim < NumberOfDimensions; i_dim++) {
        int current_old_interval{};
        int current_new_interval{1};
        double d_accu = 0;
        while (true) {
            d_accu += delta_weights[i_dim];
            while (d_accu > smoothed_weights[i_dim][current_old_interval]) {
                d_accu -= smoothed_weights[i_dim][current_old_interval];
                current_old_interval++;
            }
            x_edges[i_dim][current_new_interval] = x_edges_last[i_dim][current_old_interval] +
                                                   d_accu / smoothed_weights[i_dim][current_old_interval] *
                                                   dx_steps_last[i_dim][current_old_interval];
            dx_steps[i_dim][current_new_interval - 1] =
                    x_edges[i_dim][current_new_interval] - x_edges[i_dim][current_new_interval - 1];
            current_new_interval++;
            if (current_new_interval >= NumberOfIntervals) {
                break;
            }
        }
        dx_steps[i_dim][NumberOfIntervals - 1] =
                x_edges[i_dim][number_of_edges - 1] - x_edges[i_dim][number_of_edges - 2];
    }
    reset_weight();
}

template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::check_weight() {
    for (int i_dim = 0; i_dim < NumberOfDimensions; i_dim++) {
        average_weight[i_dim] = 0;
        for (int i = 0; i < weights[i_dim].size(); i++) {
            average_weight[i_dim] += weights[i_dim][i];
        }
        average_weight[i_dim] /= static_cast<double>(weights[i_dim].size());
        for (int i = 0; i < weights[i_dim].size(); i++) {
            std_weight[i_dim] += (weights[i_dim][i] - average_weight[i_dim])*(weights[i_dim][i] - average_weight[i_dim]);
        }
        std_weight[i_dim] = sqrt(std_weight[i_dim]);// /average_weight;
    }
}

template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::print_edges() {
    std::cout << "Grid Map:" << std::endl;
    for (int i_dim = 0; i_dim < NumberOfDimensions; i_dim++) {
        std::cout << "\tx_" << i_dim << ":";
        for (int i = 0; i < number_of_edges; i++) {
            std::cout << "\t" << x_edges[i_dim][i];
        }
        std::cout << std::endl;
        std::cout << "\tdx_" << i_dim << ":";
        for (int i = 0; i < NumberOfIntervals; i++) {
            std::cout << "\t" << dx_steps[i_dim][i];
        }
        std::cout << std::endl;
    }
}

template <int NumberOfDimensions, int NumberOfIntervals>
void VegasMap<NumberOfDimensions, NumberOfIntervals>::print_weights() {
    std::cout << "Weights:" << std::endl;
    for (int i_dim = 0; i_dim < NumberOfDimensions; i_dim++) {
        std::cout << "\tweight_" << i_dim << ":";
        for (int i = 0; i < NumberOfIntervals; i++) {
            std::cout << "\t" << weights[i_dim][i];
        }
        std::cout << std::endl;
        check_weight();
        std::cout << "\t\tAverage: " << average_weight[i_dim] << "  STD: " << std_weight[i_dim] << " std/ave: "
                  << std_weight[i_dim] / average_weight[i_dim] / NumberOfIntervals << std::endl;
    }
}

template <int NumberOfDimensions, int NumberOfIntervals>
double VegasMap<NumberOfDimensions, NumberOfIntervals>::checking_map() {
    double dx_ave = 1.0 / NumberOfIntervals;
    double chi2{};
    for (int idim = 0; idim < NumberOfDimensions; idim++) {
        for (int i = 0; i < number_of_edges; i++) {
            chi2 += (x_edges[idim][i] - x_edges_last[idim][i])*(x_edges[idim][i] - x_edges_last[idim][i]) / dx_ave/ dx_ave;
        }
    }
    return chi2 / NumberOfDimensions / number_of_edges;
}