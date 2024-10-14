#ifndef VEGAS_MAP_H
#define VEGAS_MAP_H
#include <vector>
#include <string>

// * This is the grid map from y to x:
// * y: the variable upon which we generate uniformly distributed random numbers
// * x: the integral variable, the upper and lower limit set to be 0 and 1 (The mapping from true integrate limits to [0,1] should be done by users in the integrand)
// * The function of this class:
// * 1. Keep track the mapping between y to x
// * 2. Keep the Jacobian from y to x
// * 3. Take care of the grid map improvements
class VegasMap
{
private:
    int number_of_dimensions;
    int number_of_intervals;
    int number_of_edges; // number_of_intervals + 1;
    double alpha; // The parameter control the smooth of weight
    
    std::vector<std::vector<double> > x_edges; // The edges in x, size = dimensions x number_of_edges;
    std::vector<std::vector<double> > dx_steps; // The step for each interval, size = dimensions x number_of_intervals;

    std::vector<std::vector<double> > x_edges_last; // The edges in x, size = dimensions x number_of_edges;
    std::vector<std::vector<double> > dx_steps_last; // The step for each interval, size = dimensions x number_of_intervals;
    
    std::vector<std::vector<double> > weights; // The weight in each interval, used to improve the grid map, size = dimensions x number_of_intervals;
    std::vector<std::vector<double> > counts; // Count the numbers of random numbers in specific interval, size = dimensions x number_of_intervals;
    
    std::vector<std::vector<double> > smoothed_weights; // Smoothed weights, also renormalized, size = dimensions x  number_of_intervals
    std::vector<double> summed_weights; // The all summed smoothed weights, size = dimensions
    std::vector<double> delta_weights; // The step for weights, size = dimensions

    void smooth_weight();
    void reset_weight();
    void check_weight();
    std::vector<double> average_weight; // size = dimensions
    std::vector<double> std_weight; // size = dimensions
    std::vector<int> get_interval_ID(const std::vector<double> & y) const;
    std::vector<double> get_interval_offset(const std::vector<double> & y) const;


public:
    VegasMap();
    explicit VegasMap(int dimensions);
    VegasMap(int NDIM, int Intervals);
    ~VegasMap(){};

    void reset_map();
    void set_alpha(double alp){ alpha = alp;};
    void accumulate_weight(const std::vector<double> & y, double f); // f is the integrand, no other manupulation
    void update_map();

    int Get_N_Interval() const {return number_of_intervals;}

    std::vector<double> get_x(const std::vector<double> &y) const;
    double get_jacobian(const std::vector<double> &y);
    
    void print_edges();
    void print_weights();

    double checking_map();
};



#endif // VEGAS_MAP_H