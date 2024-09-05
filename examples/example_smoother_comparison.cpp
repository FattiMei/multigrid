#include <iostream>
#include <memory>
#include "symbolic.h"
#include "poisson.hpp"
#include "solvers.hpp"
#include "multigrid.hpp"


int main() {
	constexpr int n       = (1 << 10) + 1;
	constexpr int maxiter = 15;

	IsotropicPoisson2D problem(
		{0.0, 0.0},
		{1.0, 1.0},
		n,
		forcing_term_2d,
		solution_2d
	);

	const auto restriction  = full_weight_restriction_2d;
	const auto prolongation = linear_prolongation_2d;

	const std::vector<std::pair<std::string,std::shared_ptr<IterativeSolver>>> solvers{
		{"jacobi"         , std::make_shared<MgSolver>(&problem, MgCycle::V(5, 3), InitializationStrategy::Zeros, UpdateStrategy::Jacobi,         restriction, prolongation)},
		{"sor"            , std::make_shared<MgSolver>(&problem, MgCycle::V(5, 3), InitializationStrategy::Zeros, UpdateStrategy::SOR,         restriction, prolongation)},
		{"GS injective"   , std::make_shared<MgSolver>(&problem, MgCycle::V(5, 3), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, injection_restriction_2d, prolongation)},
		{"GS full weight" , std::make_shared<MgSolver>(&problem, MgCycle::V(5, 3), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, full_weight_restriction_2d, prolongation)},
		{"red-black"      , std::make_shared<MgSolver>(&problem, MgCycle::V(5, 3), InitializationStrategy::Zeros, UpdateStrategy::RedBlack,    restriction, prolongation)},
	};

	std::cout << "it";
	for (auto [label, _] : solvers) {
		std::cout << "," << label;
	}
	std::cout << std::endl;


	for (int i = 0; i < maxiter; ++i) {
		std::cout << i;

		for (auto& [_, solver] : solvers) {
			std::cout << "," << solver->get_residual_norm();
			solver->step();
		}

		std::cout << std::endl;
	}
}
