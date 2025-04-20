#include "utils.hpp"
#include <cassert>


std::vector<double> linspace(double inf, double sup, int n) {
	assert(n > 0);
	std::vector<double> mesh(n);

	for (int i = 0; i < n; ++i) {
		mesh[i] = inf + static_cast<double>(i) * (sup - inf) / static_cast<double>(n - 1);
	}

	return mesh;
}
