#include <iostream>
#include "symbolic.h"
#include "poisson.hpp"
#include "direct.hpp"
#include "multigrid.hpp"
#include "utils.hpp"
#include <eigen3/Eigen/SparseLU>


int main() {
	const std::pair<double,double> bottom_left_corner{0.0, 0.0};
	const std::pair<double,double> top_right_corner{1.0, 1.0};

	std::cout << "n,err" << std::endl;

	for (int n = 2; n <= 128; n *= 2) {
		const int m = n + 1;
		const IsotropicPoisson2D problem(
			bottom_left_corner,
			top_right_corner,
			m,
			forcing_term_2d,
			solution_2d
		);

		DirectSolver<Eigen::SparseLU<Eigen::SparseMatrix<double>>> direct(&problem);

		MgSolver iterative(
			&problem,
			MgCycle::V(2),
			InitializationStrategy::Zeros,
			UpdateStrategy::GaussSeidel
		);

		direct.solve();

		std::vector<double> exact_solution(m*m);
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < m; ++j) {
				exact_solution[i*m+j] = solution_2d(
					linspace(bottom_left_corner.first, top_right_corner.first, m, j),
					linspace(bottom_left_corner.second, top_right_corner.second, m, i)
				);
			}
		}

		std::cout
			<< n
			<< ','
			<< linfnorm(exact_solution, direct.get_solution())
			<< std::endl;
	}

	return 0;
}
