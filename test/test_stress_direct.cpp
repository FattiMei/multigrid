#include <iostream>
#include "poisson.hpp"
#include "direct.hpp"
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/SparseLU>


int main() {
	const double inf = 0.0;
	const double sup = 1.0;
	const std::function<double(double)> f = [](double x) {return x;};
	const std::pair<double,double> boundary{0.0, 0.0};

	std::cout << "n,SimplicialLDLT,SparseLU,Symmetrized" << std::endl;

	for (int n = 10; n <= 1'000'000; n *= 10) {
		const IsotropicPoisson1D problem(inf, sup, n, f, boundary);

		DirectSolver<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> ldlt(&problem);
		DirectSolver<Eigen::SparseLU<Eigen::SparseMatrix<double>>>       splu(&problem);
		SymmetricSolver<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>> symm(&problem);

		ldlt.solve();
		splu.solve();
		symm.solve();

		std::cout
			<< n
			<< ','
			<< ldlt.get_residual_norm()
			<< ','
			<< splu.get_residual_norm()
			<< ','
			<< symm.get_residual_norm()
			<< std::endl;
	}

	return 0;
}
