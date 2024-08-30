#include <iostream>
#include <cmath>
#include <chrono>
#include <memory>
#include "circuit.hpp"
#include "direct.hpp"
#include "solvers.hpp"
#include "multigrid.hpp"
#include "utils.hpp"


const double tau = 1.0;


double exact(double x) {
	return 3.0 * std::sin(2.0 * x + M_PI / 4.0) + 2.0 * std::cos(3.0 * x - M_PI / 6.0);
}


double f(double x) {
	// u' + tau * u
	return 6.0 * std::cos(2.0 * x + M_PI / 4.0) - 6.0 * std::sin(3.0 * x - M_PI / 6.0) + tau * exact(x);
}


std::vector<double> compute_exact_solution(ForcedRC &problem) {
	const auto& mesh = problem.get_mesh();

	std::vector<double> solution(mesh.size());
	for (size_t i = 0; i < mesh.size(); ++i) {
		solution[i] = exact(mesh[i]);
	}

	return solution;
}


int main() {
	ForcedRC problem(
		4.0 * M_PI,
		tau,
		1'000'001,
		f
	);
	const auto exact_solution = compute_exact_solution(problem);
	const std::vector<std::pair<std::string,std::shared_ptr<IterativeSolver>>> solvers{
		{"2-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(2, 3), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, full_weight_restriction_1d, linear_prolongation_2d)},
		{"3-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(3, 3), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, full_weight_restriction_1d, linear_prolongation_2d)},
		{"5-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(5, 3), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, full_weight_restriction_1d, linear_prolongation_2d)},
		{"7-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(7, 3), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, full_weight_restriction_1d, linear_prolongation_2d)},
	};
	double tol = 1e-16;

	std::cerr << "solver,err,walltime" << std::endl;

	{
		const auto t1 = std::chrono::high_resolution_clock::now();

		DirectSolver<Eigen::SparseLU<Eigen::SparseMatrix<double>>> solver(&problem);
		solver.solve();

		const auto t2 = std::chrono::high_resolution_clock::now();

		std::cerr
			<< "direct"
			<< ','
			<< linfnorm(exact_solution, solver.get_solution())
			<< ','
			<< std::chrono::duration<double,std::milli>(t2-t1)
			<< std::endl;

		tol = solver.get_residual_norm();
	}

	for (auto& [label, solver] : solvers) {
		const auto t1 = std::chrono::high_resolution_clock::now();

		solver->solve(tol, 10);

		const auto t2 = std::chrono::high_resolution_clock::now();

		std::cerr
			<< label
			<< ','
			<< linfnorm(exact_solution, solver->get_solution())
			<< ','
			<< std::chrono::duration<double,std::milli>(t2-t1)
			<< std::endl;
	}


	const auto& mesh = problem.get_mesh();
	const auto& in   = problem.get_rhs();
	const auto& u    = solvers[0].second->get_solution();

	std::cout << "t,input,output" << std::endl;
	for (size_t i = 0; i < mesh.size(); ++i) {
		std::cout
			<< mesh[i]
			<< ','
			<< in[i]
			<< ','
			<< u[i]
			<< std::endl;
	}

	return 0;
}
