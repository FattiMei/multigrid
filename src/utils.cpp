#include "utils.hpp"
#include <random>


double linspace(double inf, double sup, int n, int i) {
	return inf + static_cast<double>(i) * (sup - inf) / static_cast<double>(n - 1);
}


std::vector<double> generate_random_vector(int n, double inf, double sup, int seed) {
	std::default_random_engine gen(seed);
	std::uniform_real_distribution<double> uniform(inf, sup);
	std::vector<double> result(n);

	for (auto &x : result) {
		x = uniform(gen);
	}

	return result;
}
