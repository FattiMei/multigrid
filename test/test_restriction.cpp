#include <iostream>
#include <array>
#include <cmath>
#include "multigrid.hpp"


void print_2d(const std::vector<double>& X) {
	const int N = std::sqrt(X.size());

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			std::cout << X[i*N+j] << ' ';
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}


int main() {
	std::vector<double> fine {
		0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 2.0, 3.0, 0.0,
		0.0, 4.0, 5.0, 6.0, 0.0,
		0.0, 7.0, 8.0, 9.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0
	};

	std::vector<double> coarse(3*3);
	full_weight_restriction_2d({5,5}, fine.data(), coarse.data());

	print_2d(fine);
	print_2d(coarse);

	linear_prolongation_2d({5,5}, coarse.data(), fine.data());
	print_2d(fine);

	return 0;
}
