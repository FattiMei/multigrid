#include <iostream>
#include <eigen3/Eigen/Sparse>
#include "stencil.hpp"


int main() {
	ThreePointStencilOperator S(5, {1, 2, 3});
	std::cout << S.build_sparse_repr();

	return 0;
}
