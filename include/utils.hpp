#ifndef __UTILS_HPP__
#define __UTILS_HPP__


#include <vector>


double linspace(double inf, double sup, int n, int i);
std::vector<double> generate_random_vector(int n, double inf, double sup, int seed = 0);


#endif
