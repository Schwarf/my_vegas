#include "VEGAS_map.h"
#include <cmath>
#include <numeric>
#include <iostream>


VegasMap::VegasMap() {
    number_of_dimensions = 1;
    number_of_intervals = 1000;
    number_of_edges = number_of_intervals + 1;
    alpha = 0.5;
    reset_map();
    reset_map();
}

VegasMap::VegasMap(int dimensions) {
    number_of_dimensions = dimensions;
    number_of_intervals = 1000;
    number_of_edges = number_of_intervals + 1;
    alpha = 0.5;
    reset_map();
}

VegasMap::VegasMap(int NDIM, int Intervals) {
    number_of_dimensions = NDIM;
    number_of_intervals = Intervals;
    number_of_edges = number_of_intervals + 1;
    alpha = 0.5;
    reset_map();
}

void VegasMap::reset_map() {
    x_edges.clear();
    dx_steps.clear();
    double step_tmp = 1.0 / number_of_intervals;
    std::vector<double> x_edges_tmp(number_of_edges);
    std::vector<double> dx_steps_tmp(number_of_edges - 1);
    for (int i = 1; i < number_of_edges; i++) {
        x_edges_tmp[i] = i * step_tmp;
        dx_steps_tmp[i - 1] = x_edges_tmp[i] - x_edges_tmp[i - 1];
    }
    x_edges = std::vector<std::vector<double> >(number_of_dimensions, x_edges_tmp);
    dx_steps = std::vector<std::vector<double> >(number_of_dimensions, dx_steps_tmp);
    x_edges_last = std::vector<std::vector<double> >(number_of_dimensions, x_edges_tmp);
    dx_steps_last = std::vector<std::vector<double> >(number_of_dimensions, dx_steps_tmp);
    weights = std::vector<std::vector<double> >(number_of_dimensions, std::vector<double>(number_of_intervals, 0.0));
    counts = std::vector<std::vector<double> >(number_of_dimensions, std::vector<double>(number_of_intervals, 0.0));
    smoothed_weights = std::vector<std::vector<double> >(number_of_dimensions,
                                                         std::vector<double>(number_of_intervals, 0.0));
    summed_weights = std::vector<double>(number_of_dimensions, 0.0);
    delta_weights = std::vector<double>(number_of_dimensions, 0.0);
    average_weight = std::vector<double>(number_of_dimensions, 0.0);
    std_weight = std::vector<double>(number_of_dimensions, 0.0);
}

void VegasMap::reset_weight() {
    weights = std::vector<std::vector<double> >(number_of_dimensions, std::vector<double>(number_of_intervals, 0));
    counts = std::vector<std::vector<double> >(number_of_dimensions, std::vector<double>(number_of_intervals, 0));
}

std::vector<int> VegasMap::get_interval_ID(const std::vector<double> &y) const {
    std::vector<int> result(number_of_dimensions);
    for (int i = 0; i < number_of_dimensions; i++) {
        result[i] = std::floor(y[i] * number_of_intervals);
    }
    return result; // Profiling spot compare with move semantics
}

std::vector<double> VegasMap::get_interval_offset(const std::vector<double> &y) const {
    auto ID = get_interval_ID(y);
    std::vector<double> res(number_of_dimensions);
    for (int i = 0; i < number_of_dimensions; i++) {
        res[i] = y[i] * number_of_intervals - ID[i];
    }
    return res; // Profiling spot compare with move semantics
}

std::vector<double> VegasMap::get_x(const std::vector<double> &y) const {
    auto ID = get_interval_ID(y);
    auto offset = get_interval_offset(y);
    std::vector<double> res(number_of_dimensions);
    for (int i = 0; i < number_of_dimensions; i++) {
        int id = ID[i];
        res[i] = x_edges[i][id] + dx_steps[i][id] * offset[i];
    }
    return res; // Profiling spot compare with move semantics
}

double VegasMap::get_jacobian(const std::vector<double> &y) {
    auto ID = get_interval_ID(y);
    double jacobian{1.0};
    for (int i = 0; i < number_of_dimensions; i++) {
        int id = ID[i];
        jacobian *= number_of_intervals * dx_steps[i][id];
    }
    return jacobian;
}

void VegasMap::accumulate_weight(const std::vector<double> &y, double f) {
    // f is the value of integrand!
    auto ID = get_interval_ID(y);
    for (int i = 0; i < number_of_dimensions; i++) {
        int id = ID[i];
        weights[i][id] += pow(f * get_jacobian(y), 2);
        counts[i][id] += 1;
        // std::cout<<"ID: "<<id<<" weight: "<<weights[i][id]<<" counts: "<<counts[i][id]<<std::endl;
    }
}
void VegasMap::smooth_weight() {
    // std::cout<<"Smoothing weight"<<std::endl;
    for (int i_dim = 0; i_dim < number_of_dimensions; i_dim++) {
        for (int i_inter = 0; i_inter < weights[i_dim].size(); i_inter++) {
            if (counts[i_dim][i_inter] != 0) {
                weights[i_dim][i_inter] /= counts[i_dim][i_inter];
            }
        }
    }
    // std::cout<<"Count devided!"<<std::endl;
    for (int i_dim = 0; i_dim < number_of_dimensions; i_dim++) {
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

        // Handle the main loop (1 <= i < number_of_intervals - 1)
        for (int i = 1; i < number_of_intervals - 1; i++) {
            d_tmp = (weights[i_dim][i - 1] + 6.0 * weights[i_dim][i] + weights[i_dim][i + 1]) / (8.0 * d_sum);
            if (d_tmp != 0.0) {
                d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
            }
            smoothed_weights[i_dim][i] = d_tmp;
            summed_weights[i_dim] += d_tmp;
        }
        // Handle the last interval (i == number_of_intervals - 1) outside the loop
        d_tmp = (weights[i_dim][number_of_intervals - 2] + 7.0 * weights[i_dim][number_of_intervals - 1]) /
                (8.0 * d_sum);
        if (d_tmp != 0.0) {
            d_tmp = pow((d_tmp - 1.0) / log(d_tmp), alpha);
        }

        smoothed_weights[i_dim][number_of_intervals - 1] = d_tmp;
        summed_weights[i_dim] += d_tmp;

        delta_weights[i_dim] = summed_weights[i_dim] / number_of_intervals;

    }
}

void VegasMap::update_map() {
    smooth_weight();
    x_edges_last = x_edges;
    dx_steps_last = dx_steps;
    for (int i_dim = 0; i_dim < number_of_dimensions; i_dim++) {
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
            if (current_new_interval >= number_of_intervals) {
                break;
            }
        }
        dx_steps[i_dim][number_of_intervals - 1] =
                x_edges[i_dim][number_of_edges - 1] - x_edges[i_dim][number_of_edges - 2];
    }
    reset_weight();
}

void VegasMap::check_weight() {
    for (int i_dim = 0; i_dim < number_of_dimensions; i_dim++) {
        average_weight[i_dim] = 0;
        for (int i = 0; i < weights[i_dim].size(); i++) {
            average_weight[i_dim] += weights[i_dim][i];
        }
        average_weight[i_dim] /= static_cast<double>(weights[i_dim].size());
        for (int i = 0; i < weights[i_dim].size(); i++) {
            std_weight[i_dim] += pow(weights[i_dim][i] - average_weight[i_dim], 2);
        }
        std_weight[i_dim] = sqrt(std_weight[i_dim]);// /average_weight;
    }
}

void VegasMap::print_edges() {
    std::cout << "Grid Map:" << std::endl;
    for (int i_dim = 0; i_dim < number_of_dimensions; i_dim++) {
        std::cout << "\tx_" << i_dim << ":";
        for (int i = 0; i < number_of_edges; i++) {
            std::cout << "\t" << x_edges[i_dim][i];
        }
        std::cout << std::endl;
        std::cout << "\tdx_" << i_dim << ":";
        for (int i = 0; i < number_of_intervals; i++) {
            std::cout << "\t" << dx_steps[i_dim][i];
        }
        std::cout << std::endl;
    }
}

void VegasMap::print_weights() {
    std::cout << "Weights:" << std::endl;
    for (int i_dim = 0; i_dim < number_of_dimensions; i_dim++) {
        std::cout << "\tweight_" << i_dim << ":";
        for (int i = 0; i < number_of_intervals; i++) {
            std::cout << "\t" << weights[i_dim][i];
        }
        std::cout << std::endl;
        check_weight();
        std::cout << "\t\tAverage: " << average_weight[i_dim] << "  STD: " << std_weight[i_dim] << " std/ave: "
                  << std_weight[i_dim] / average_weight[i_dim] / number_of_intervals << std::endl;
    }
}

double VegasMap::checking_map() {
    double dx_ave = 1.0 / number_of_intervals;
    double chi2{};
    for (int idim = 0; idim < number_of_dimensions; idim++) {
        for (int i = 0; i < number_of_edges; i++) {
            chi2 += pow(x_edges[idim][i] - x_edges_last[idim][i], 2) / pow(dx_ave, 2);
        }
    }
    return chi2 / number_of_dimensions / number_of_edges;
}