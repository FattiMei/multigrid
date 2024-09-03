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
	std::vector<double> fine(9*9);
	for (size_t i = 0; i < fine.size(); ++i) fine[i] = i;

	std::vector<double> coarse(5*5);
	full_weight_restriction_2d({9,9}, fine.data(), coarse.data());

	print_2d(fine);
	print_2d(coarse);

	linear_prolongation_2d({9,9}, coarse.data(), fine.data());
	print_2d(fine);

	return 0;
}
