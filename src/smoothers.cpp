#include "smoothers.hpp"


using namespace Smoother;


void Jacobi::operator()(Update formula, const int n, const double b[], double u[]) {
	// @DESIGN: this condition will be met at the first call of the function and then will never happen. The branch predictor will optimize
	if (n > static_cast<int>(local.size())) {
		local.resize(n);
	}

	// Jersey style, does two sweeps
	for (int i = 1; i < n-1; ++i) {
		formula(i, b, u, local.data());
	}

	for (int i = 1; i < n-1; ++i) {
		formula(i, b, local.data(), u);
	}
}


void GSeidel::operator()(Update formula, const int n, const double b[], double u[]) {
	for (int i = 1; i < n-1; ++i) {
		formula(i, b, u, u);
	}
}


void RedBlack::operator()(Update formula, const int n, const double b[], double u[]) {
	for (int i = 1; i < n-1; i += 2) {
		formula(i, b, u, u);
	}

	for (int i = 2; i < n-1; i += 2) {
		formula(i, b, u, u);
	}
}


void BlackRed::operator()(Update formula, const int n, const double b[], double u[]) {
	for (int i = 2; i < n-1; i += 2) {
		formula(i, b, u, u);
	}

	for (int i = 1; i < n-1; i += 2) {
		formula(i, b, u, u);
	}
}
