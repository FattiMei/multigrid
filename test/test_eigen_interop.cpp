#include <eigen3/Eigen/Dense>
#include "utils.hpp"


bool test_pointer_to_eigen_vector(const int n, const double data[]) {
	Eigen::Map<const Eigen::VectorXd> equivalent(data, n);

	for (int i = 0; i < n; ++i) {
		if (equivalent[i] != data[i]) return false;
	}

	return true;
}


bool test_eigen_vector_to_pointer(const int n) {
	std::vector<double> tmp(n);
	Eigen::Map<Eigen::VectorXd> x(tmp.data(), tmp.size());

	x = Eigen::VectorXd::Random(n);

	for (int i = 0; i < n; ++i) {
		if (x[i] != tmp[i]) {
			return false;
		}
	}

	return true;
}


int main() {
	std::vector<double> arr = generate_random_vector(1'000'000, 0.0, 1.0);

	assert(test_pointer_to_eigen_vector(arr.size(), arr.data()));
	assert(test_eigen_vector_to_pointer(arr.size()));

	return 0;
}
