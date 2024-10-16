#ifndef VEGAS_MAP_H
#define VEGAS_MAP_H

#include <vector>
#include <array>
#include <string>

// * This is the grid mapping from y to x:
// * y: the variable upon which we generate uniformly distributed random numbers
// * x: the integral variable, the upper and lower limit set to be 0 and 1 (The mapping from true integrate limits to [0,1] should be done by users in the integrand)
// * The function of this class:
// * 1. Keep track the mapping between y to x
// * 2. Keep the Jacobian from y to x
// * 3. Take care of the grid mapping improvements
template<int NumberOfDimensions, int NumberOfIntervals = 1000>
class VegasMap {
public:
    VegasMap();
    ~VegasMap() = default;
    std::array<double, NumberOfDimensions> get_x(const std::array<double, NumberOfDimensions> &random_numbers);
    void accumulate_weights(double evaluated_integrand); // evaluated_integrand is the integrand, no other manupulation
    double get_jacobian();
    void update_map();
    double checking_map();
    void set_alpha(double alp) { alpha = alp; };
    void print_edges();
    void print_weights();

private:
    double alpha{0.5}; // The parameter control the smooth of weight
    std::array<std::array<double, NumberOfIntervals + 1>, NumberOfDimensions> x_edges;
    std::array<std::array<double, NumberOfIntervals>, NumberOfDimensions> dx_steps;
    std::array<std::array<double, NumberOfIntervals + 1>, NumberOfDimensions> x_edges_last;
    std::array<std::array<double, NumberOfIntervals>, NumberOfDimensions> dx_steps_last;
    std::array<std::array<double, NumberOfIntervals>, NumberOfDimensions> weights = {};
    std::array<std::array<double, NumberOfIntervals>, NumberOfDimensions> counts = {};
    std::array<std::array<double, NumberOfIntervals>, NumberOfDimensions> smoothed_weights = {};
    std::array<double, NumberOfDimensions> summed_weights = {};
    std::array<double, NumberOfDimensions> delta_weights = {};
    std::array<double, NumberOfDimensions> average_weight = {};
    std::array<double, NumberOfDimensions> std_weight = {};
    std::array<int, NumberOfDimensions> ID = {};

    void smooth_weights();
    void reset_weights();
    void check_weight();
    void compute_interval_ID(const std::array<double, NumberOfDimensions> &random_numbers);
    std::array<double, NumberOfDimensions>get_interval_offset(const std::array<double, NumberOfDimensions> &random_numbers) const;
};

#include "VEGAS_map.inl"


#endif // VEGAS_MAP_H