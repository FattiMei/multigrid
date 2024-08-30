#include <iostream>
#include "symbolic.h"
#include "poisson.hpp"
#include "direct.hpp"
#include "multigrid.hpp"
#include "utils.hpp"
#include <eigen3/Eigen/SparseLU>


int main() {
	const double inf = 0.0;
	const double sup = 1.0;
	const std::function<double(double)> exact = [](double x) { return solution_1d(x); };
	const std::function<double(double)> f     = [](double x) { return forcing_term_1d(x); };
	const std::pair<double,double> boundary{exact(inf), exact(sup)};

	std::cout << "n,direct,naive,precise" << std::endl;

	for (int n = 10; n <= 1'000'000; n *= 10) {
		const int m = n + 1;
		const IsotropicPoisson1D    problem    (inf, sup, m, f, boundary);
		const PoissonPreciseVariant alternative(inf, sup, m, f, boundary);

		DirectSolver<Eigen::SparseLU<Eigen::SparseMatrix<double>>> direct(&problem);
		MgSolver iterative(
			&problem,
			MgCycle::V(2),
			InitializationStrategy::Zeros,
			UpdateStrategy::GaussSeidel
		);
		MgSolver precise(
			&alternative,
			MgCycle::V(2),
			InitializationStrategy::Zeros,
			UpdateStrategy::GaussSeidel
		);

		direct.solve();
		iterative.solve(1e-16, 100);
		precise.solve(1e-16, 100);

		std::vector<double> exact_solution(m);
		for (int i = 0; i < m; ++i) {
			exact_solution[i] = exact(linspace(inf, sup, m, i));
		}

		std::cout
			<< n
			<< ','
			<< l2norm(exact_solution, direct.get_solution())
			<< ','
			<< l2norm(exact_solution, iterative.get_solution())
			<< ','
			<< l2norm(exact_solution, precise.get_solution())
			<< std::endl;
	}

	return 0;
}
