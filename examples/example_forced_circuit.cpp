#include <iostream>
#include <cmath>
#include "circuit.hpp"
#include "direct.hpp"


double f(double x) {
	return 3.0 * std::sin(2.0 * x + M_PI / 4.0) + 2.0 * std::cos(3.0 * x - M_PI / 6.0);
}


int main() {
	ForcedRC problem(
		4.0 * M_PI,
		1.0,
		10'000'000,
		f
	);

	DirectSolver<Eigen::SparseLU<Eigen::SparseMatrix<double>>> solver(&problem);
	solver.solve();

	std::cout << "t,input,output" << std::endl;

	const auto& mesh = problem.get_mesh();
	const auto& rhs  = problem.get_rhs();
	const auto& sol  = solver.get_solution();

	for (size_t i = 0; i < mesh.size(); ++i) {
		std::cout 
			<< mesh[i]
			<< ','
			<< rhs[i]
			<< ','
			<< sol[i]
			<< std::endl;
	}

	return 0;
}
