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


/*
 * TODO: add resize logic
void Jacobi::operator()(const SparseOperator &A, const double b[], double u[]) {
	for (int i = 0; i < A.problem_size; ++i) {
		double acc = 0.0;

		for (int j = A.offset[i]+1; j < A.offset[i+1]; ++j) {
			acc += A.weights[j] * u[A.neighbours[j]];
		}

		local[i] = (b[i] - acc) / A.weights[A.offset[i]];
	}

	for (int i = 0; i < A.problem_size; ++i) {
		double acc = 0.0;

		for (int j = A.offset[i]+1; j < A.offset[i+1]; ++j) {
			acc += A.weights[j] * local[A.neighbours[j]];
		}

		u[i] = (b[i] - acc) / A.weights[A.offset[i]];
	}
}
*/


void Jacobi::operator()(const ThreePointStencilOperator &A, const double b[], double u[]) {
	const int n      = A.problem_size;
	const double left   = A.stencil[0];
	const double center = A.stencil[1];
	const double right  = A.stencil[2];

	if (n > static_cast<int>(local.size())) {
		local.resize(n);
	}

	for (int i = 1; i < n-1; ++i) {
		local[i] = (b[i] - left * u[i-1] - right * u[i+1]) / center;
	}

	for (int i = 1; i < n-1; ++i) {
		u[i] = (b[i] - left * local[i-1] - right * local[i+1]) / center;
	}
}


void GSeidel::operator()(Update formula, const int n, const double b[], double u[]) {
	for (int i = 1; i < n-1; ++i) {
		formula(i, b, u, u);
	}
}


void GSeidel::operator()(const ThreePointStencilOperator &A, const double b[], double u[]) {
	const int    n      = A.problem_size;
	const double left   = A.stencil[0];
	const double center = A.stencil[1];
	const double right  = A.stencil[2];

	for (int i = 1; i < n-1; ++i) {
		u[i] = (b[i] - left * u[i-1] - right * u[i+1]) / center;
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


void RedBlack::operator()(const ThreePointStencilOperator &A, const double b[], double u[]) {
	const int    n      = A.problem_size;
	const double left   = A.stencil[0];
	const double center = A.stencil[1];
	const double right  = A.stencil[2];

	for (int i = 2; i < n-1; i += 2) {
		u[i] = (b[i] - left * u[i-1] - right * u[i+1]) / center;
	}

	for (int i = 1; i < n-1; i += 2) {
		u[i] = (b[i] - left * u[i-1] - right * u[i+1]) / center;
	}
}
