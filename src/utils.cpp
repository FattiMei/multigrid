#include "utils.hpp"
#include <cmath>
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


double linfnorm(const std::vector<double> &x, const std::vector<double> &y) {
	double max = 0.0;

	for (size_t i = 0; i < x.size(); ++i) {
		max = std::max(max, std::abs(x[i] - y[i]));
	}

	return max;
}


double l2norm(const std::vector<double> &x, const std::vector<double> &y) {
	double acc = 0.0;

	for (size_t i = 0; i < x.size(); ++i) {
		const double diff = x[i] - y[i];

		acc += diff * diff;
	}

	return std::sqrt(acc);
}
