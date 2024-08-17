#include <iostream>
#include <chrono>
#include "poisson.hpp"
#include "solvers.hpp"


int main() {
	const double inf = 0.0;
	const double sup = 1.0;
	const std::function<double(double)> f = [](double x) {return x;};
	const std::pair<double,double> boundary{0.0, 0.0};

	std::cout << "n,residual,wall_time[ms]" << std::endl;

	for (int n = 10; n <= 1'000'000; n *= 10) {
		const IsotropicPoisson1D problem(inf, sup, n, f, boundary);

		const auto start = std::chrono::high_resolution_clock::now();
		EigenDirectSolver solver(&problem);
		solver.solve();
		const auto stop = std::chrono::high_resolution_clock::now();

		std::cout << n << ',' << solver.get_residual_norm() << ',' << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << std::endl;
	}

	return 0;
}
