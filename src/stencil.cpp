#include "stencil.hpp"
#include <cmath>
#include <iostream>


NaiveThreePointStencil::NaiveThreePointStencil(
	const int n,
	const std::array<double,3> weights
) :
	DiscreteOperator(n),
	stencil(weights)
{}


// there is no silver bullet
void NaiveThreePointStencil::relax(const double b[], double u[], UpdateStrategy strategy) {
	switch (strategy) {
		// double sweep, Jersey style
		case UpdateStrategy::SOR:
		case UpdateStrategy::Jacobi: {
			if (local.size() < static_cast<size_t>(n)) local.resize(n);

			local[0] = b[0];

			for (int i = 1; i < n-1; ++i) {
				local[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

			local[n-1] = b[n-1];

			for (int i = 0; i < n; ++i) {
				u[i] = local[i];
			}

		} break;

		case UpdateStrategy::GaussSeidel: {
			u[0] = b[0];
			u[n-1] = b[n-1];

			for (int i = 1; i < n-1; ++i) {
				u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

		} break;

		case UpdateStrategy::RedBlack: {
			u[0] = b[0];

			for (int i = 1; i < n-1; i += 2) {
				u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

			for (int i = 2; i < n-1; i += 2) {
				u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

			u[n-1] = b[n-1];
		} break;
	}
}


Eigen::SparseMatrix<double> NaiveThreePointStencil::get_sparse_repr() const {
	Eigen::SparseMatrix<double> A(n,n);
	A.reserve(Eigen::VectorXi::Constant(n, 3));

	A.coeffRef(0,0) = 1.0;
	for (int i = 1; i < n-1; ++i) {
		A.coeffRef(i,i-1) = stencil[0];
		A.coeffRef(i,i)   = stencil[1];
		A.coeffRef(i,i+1) = stencil[2];
	}
	A.coeffRef(n-1,n-1) = 1.0;

	A.makeCompressed();
	return A;
}


void NaiveThreePointStencil::compute_residual(const double b[], const double u[], double r[]) const {
	r[0] = b[0] - u[0];

	for (int i = 1; i < n-1; ++i) {
		r[i] = b[i] - stencil[0] * u[i-1] - stencil[1] * u[i] - stencil[2] * u[i+1];
	}

	r[n-1] = b[n-1] - u[n-1];
}


double NaiveThreePointStencil::compute_residual_norm(const double b[], const double u[]) const {
	double acc = 0.0;
	double diff;

	diff = b[0] - u[0];
	acc += diff * diff;

	for (int i = 1; i < n-1; ++i) {
		diff = b[i] - u[i-1] * stencil[0] - u[i] * stencil[1] - u[i+1] * stencil[2];
		acc += diff * diff;
	}

	diff = b[n-1] - u[n-1];
	acc += diff * diff;

	return std::sqrt(acc);
}


ThreePointStencil::ThreePointStencil(
	const int    n,
	const double step,
	const std::array<double,3> weights
) :
	DiscreteOperator(n),
	stencil(weights),
	h(step)
{}


void ThreePointStencil::relax(const double b[], double u[], UpdateStrategy strategy) {
	switch (strategy) {
		case UpdateStrategy::SOR:
		case UpdateStrategy::Jacobi: {
			if (local.size() < static_cast<size_t>(n)) local.resize(n);

			local[0] = b[0];

			for (int i = 1; i < n-1; ++i) {
				local[i] = (h * h * b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

			local[n-1] = b[n-1];

			for (int i = 0; i < n; ++i) {
				u[i] = local[i];
			}

		} break;

		case UpdateStrategy::GaussSeidel: {
			u[0] = b[0];
			u[n-1] = b[n-1];

			for (int i = 1; i < n-1; ++i) {
				u[i] = (h * h * b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

		} break;

		case UpdateStrategy::RedBlack: {
			u[0] = b[0];

			for (int i = 1; i < n-1; i += 2) {
				u[i] = (h * h * b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

			for (int i = 2; i < n-1; i += 2) {
				u[i] = (h * h * b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

			u[n-1] = b[n-1];
		} break;
	}
}


Eigen::SparseMatrix<double> ThreePointStencil::get_sparse_repr() const {
	Eigen::SparseMatrix<double> A(n,n);
	A.reserve(Eigen::VectorXi::Constant(n, 3));

	for (int i = 1; i < n-1; ++i) {
		A.coeffRef(i,i-1) = stencil[0] / (h*h);
		A.coeffRef(i,i)   = stencil[1] / (h*h);
		A.coeffRef(i,i+1) = stencil[2] / (h*h);
	}

	A.coeffRef(0,0) = 1.0;
	A.coeffRef(n-1,n-1) = 1.0;

	A.makeCompressed();
	return A;
}


void ThreePointStencil::compute_residual(const double b[], const double u[], double r[]) const {
	r[0] = b[0] - u[0];

	for (int i = 1; i < n-1; ++i) {
		r[i] = b[i] - (stencil[0] * u[i-1] + stencil[1] * u[i] + stencil[2] * u[i+1]) / (h * h);
	}

	r[n-1] = b[n-1] - u[n-1];
}


double ThreePointStencil::compute_residual_norm(const double b[], const double u[]) const {
	double acc = 0.0;
	double diff;

	diff = b[0] - u[0];
	acc += diff * diff;

	for (int i = 1; i < n-1; ++i) {
		diff = b[i] - (stencil[0] * u[i-1] + stencil[1] * u[i] + stencil[2] * u[i+1]) / (h * h);
		acc += diff * diff;
	}

	diff = b[n-1] - u[n-1];
	acc += diff * diff;

	return std::sqrt(acc);
}


void ThreePointPeriodicStencil::relax(const double b[], double u[], UpdateStrategy strategy) {
	switch(strategy) {
		// double sweep, Jersey style
		case UpdateStrategy::SOR:
		case UpdateStrategy::Jacobi: {
			if (local.size() < static_cast<size_t>(n)) local.resize(n);

			{
				local[0] = (b[0] - stencil[0] * u[n-2] - stencil[2] * u[1]) / stencil[1];
				local[n-1] = local[0];

				for (int i = 1; i < n-1; ++i) {
					local[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
				}
			}

			{
				u[0] = (b[0] - stencil[0] * local[n-2] - stencil[2] * local[1]) / stencil[1];
				u[n-1] = u[0];

				for (int i = 1; i < n-1; ++i) {
					u[i] = (b[i] - stencil[0] * local[i-1] - stencil[2] * local[i+1]) / stencil[1];
				}
			}
		} break;

		case UpdateStrategy::GaussSeidel: {
			u[0] = (b[0] - stencil[0] * u[n-2] - stencil[2] * u[1]) / stencil[1];
			u[n-1] = u[0];

			for (int i = 1; i < n-1; ++i) {
				u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}
		} break;

		case UpdateStrategy::RedBlack: {
			u[0] = (b[0] - stencil[0] * u[n-2] - stencil[2] * u[1]) / stencil[1];
			u[n-1] = u[0];

			for (int i = 2; i < n-1; i += 2) {
				u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}

			for (int i = 1; i < n-1; i += 2) {
				u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1]) / stencil[1];
			}
		} break;
	}
}


Eigen::SparseMatrix<double> ThreePointPeriodicStencil::get_sparse_repr() const {
	Eigen::SparseMatrix<double> A(n,n);
	A.reserve(Eigen::VectorXi::Constant(n, 3));

	A.coeffRef(0,n-2) = stencil[0];
	A.coeffRef(0,0)   = stencil[1];
	A.coeffRef(0,1)   = stencil[2];

	for (int i = 1; i < n-1; ++i) {
		A.coeffRef(i,i-1) = stencil[0];
		A.coeffRef(i,i)   = stencil[1];
		A.coeffRef(i,i+1) = stencil[2];
	}

	A.coeffRef(n-1,n-2) = stencil[0];
	A.coeffRef(n-1,n-1) = stencil[1];
	A.coeffRef(n-1,1)   = stencil[2];

	A.makeCompressed();
	return A;
}


void ThreePointPeriodicStencil::compute_residual(const double b[], const double u[], double r[]) const {
	r[0] = b[0] - stencil[0] * u[n-1] - stencil[1] * u[0] - stencil[2] * u[1];

	for (int i = 1; i < n-1; ++i) {
		r[i] = b[i] - stencil[0] * u[i-1] - stencil[1] * u[i] - stencil[2] * u[i+1];
	}

	r[n-1] = b[n-1] - stencil[0] * u[n-2] - stencil[1] * u[n-1] - stencil[2] * u[0];
}


double ThreePointPeriodicStencil::compute_residual_norm(const double b[], const double u[]) const {
	double acc = 0.0;
	double diff;

	diff = b[0] - stencil[0] * u[n-1] - stencil[1] * u[0] - stencil[2] * u[1];
	acc += diff * diff;

	for (int i = 1; i < n-1; ++i) {
		diff = b[i] - stencil[0] * u[i-1] - stencil[1] * u[i] - stencil[2] * u[i+1];
		acc += diff * diff;
	}

	diff = b[n-1] - stencil[0] * u[n-2] - stencil[1] * u[n-1] - stencil[2] * u[0];
	acc += diff * diff;

	return std::sqrt(acc);
}


FivePointStencil::FivePointStencil(const int n, const int m, const std::array<double,5> weights) :
	DiscreteOperator(n * m),
	stencil(weights),
	rows(n),
	cols(m),
	local(0)
{}


void FivePointStencil::relax(const double b[], double u[], UpdateStrategy strategy) {
	switch (strategy) {
		case UpdateStrategy::SOR: {
			const double rho = 0.5;

			if (local.size() < static_cast<size_t>(n)) {
				local.resize(n);
			}

			#pragma omp parallel for
			for (int row = 1; row < rows-1; ++row) {
				const int start = cols * row;
				const int end   = start + cols - 1;

				for (int i = start + 1; i < end; ++i) {
					local[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1] - stencil[3] * u[i+cols] - stencil[4] * u[i-cols]) / stencil[1];
				}
			}

#ifdef SOR_SAME_DATA_ACCESS_PATTERN
			#pragma omp parallel for
			for (int row = 1; row < rows-1; ++row) {
				const int start = cols * row;
				const int end   = start + cols - 1;

				for (int i = start + 1; i < end; ++i) {
					u[i] = (1.0 - rho) * u[i] + rho * local[i];
				}
			}
#else
			#pragma omp parallel for
			for (int i = 0; i < n; ++i) {
				u[i] = (1.0 - rho) * u[i] + rho * local[i];
			}
#endif
		} break;

		case UpdateStrategy::Jacobi: {
			if (local.size() < static_cast<size_t>(n)) {
				local.resize(n);
			}

			#pragma omp parallel for
			for (int row = 1; row < rows-1; ++row) {
				const int start = cols * row;
				const int end   = start + cols - 1;

				for (int i = start + 1; i < end; ++i) {
					local[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1] - stencil[3] * u[i+cols] - stencil[4] * u[i-cols]) / stencil[1];
				}
			}

			#pragma omp parallel for
			for (int i = 0; i < n; ++i) {
				u[i] = local[i];
			}

		} break;

		case UpdateStrategy::GaussSeidel: {
			for (int row = 1; row < rows-1; ++row) {
				const int start = cols * row;
				const int end   = start + cols - 1;

				for (int i = start + 1; i < end; ++i) {
					u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1] - stencil[3] * u[i+cols] - stencil[4] * u[i-cols]) / stencil[1];
				}
			}
		} break;

		case UpdateStrategy::RedBlack: {
			#pragma omp parallel for
			for (int row = 1; row < rows-1; ++row) {
				const int start = cols * row;
				const int end   = start + cols - 1;

				for (int i = start + 1 + (row % 2); i < end; i += 2) {
					u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1] - stencil[3] * u[i+cols] - stencil[4] * u[i-cols]) / stencil[1];
				}
			}

			#pragma omp parallel for
			for (int row = 1; row < rows-1; ++row) {
				const int start = cols * row;
				const int end   = start + cols - 1;

				for (int i = start + 1 + (1 - (row % 2)); i < end; i += 2) {
					u[i] = (b[i] - stencil[0] * u[i-1] - stencil[2] * u[i+1] - stencil[3] * u[i+cols] - stencil[4] * u[i-cols]) / stencil[1];
				}
			}

		} break;
	}
}


Eigen::SparseMatrix<double> FivePointStencil::get_sparse_repr() const {
	Eigen::SparseMatrix<double> A(n,n);

	for (int i = 0; i < cols; ++i) {
		A.coeffRef(i,i) = 1.0;
	}

	for (int row = 1; row < rows-1; ++row) {
		const int start = cols * row;
		const int end   = start + cols - 1;

		A.coeffRef(start,start) = 1.0;

		for (int i = start + 1; i < end; ++i) {
			A.coeffRef(i,i-1)    = stencil[0];
			A.coeffRef(i,i)      = stencil[1];
			A.coeffRef(i,i+1)    = stencil[2];
			A.coeffRef(i,i+cols) = stencil[3];
			A.coeffRef(i,i-cols) = stencil[4];
		}

		A.coeffRef(end,end) = 1.0;
	}

	for (int i = cols * (rows-1); i < rows * cols; ++i) {
		A.coeffRef(i,i) = 1.0;
	}

	A.makeCompressed();
	return A;
}


void FivePointStencil::compute_residual(const double b[], const double u[], double r[]) const {
	#pragma omp parallel for
	for (int row = 1; row < rows-1; ++row) {
		const int start = cols * row;
		const int end   = start + cols - 1;

		for (int i = start + 1; i < end; ++i) {
			r[i] = b[i] - stencil[0] * u[i-1] - stencil[1] * u[i] - stencil[2] * u[i+1] - stencil[3] * u[i+cols] - stencil[4] * u[i-cols];
		}
	}
}


void FivePointStencil::compute_residual_and_restrict(const std::pair<int,int> dim, const double b[], const double u[], double dest[]) {
	const int target_cols = 1 + (dim.second - 1) / 2;
	double K[9]{};

	// 6 7 8
	// 3 4 5
	// 0 1 2

	#pragma omp parallel for private(K)
	for (int i = 2; i < dim.first-1; i += 2) {
		{
			const int linear_index = i * dim.second + 2;

			{
				const int center = linear_index - cols - 1;

				K[0] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}
			{
				const int center = linear_index - cols;

				K[1] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}
			{
				const int center = linear_index - 1;

				K[3] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}
			{
				const int center = linear_index;

				K[4] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}
			{
				const int center = linear_index + cols - 1;

				K[6] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}
			{
				const int center = linear_index + cols;

				K[7] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}
		}

		for (int j = 2; j < dim.second-1; j += 2) {
			const int linear_index = i * dim.second + j;

			// fill the third column of K, then compute the mean
			{
				const int center = linear_index - cols + 1;

				K[2] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}
			{
				const int center = linear_index + 1;

				K[5] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}
			{
				const int center = linear_index + cols + 1;

				K[8] = b[center] - stencil[0] * u[center-1] - stencil[1] * u[center] - stencil[2] * u[center+1] - stencil[3] * u[center+cols] - stencil[4] * u[center-cols];
			}

			dest[(i / 2) * target_cols + (j / 2)] =
				  0.25   * K[4]
				+ 0.125  * K[1]
				+ 0.125  * K[3]
				+ 0.125  * K[5]
				+ 0.125  * K[7]
				+ 0.0625 * K[0]
				+ 0.0625 * K[2]
				+ 0.0625 * K[6]
				+ 0.0625 * K[8];

			K[0] = K[1];
			K[1] = K[2];
			K[3] = K[4];
			K[4] = K[5];
			K[6] = K[7];
			K[7] = K[8];
		}
	}
}


double FivePointStencil::compute_residual_norm(const double b[], const double u[]) const {
	double acc = 0.0;
	double diff;

	for (int row = 1; row < rows-1; ++row) {
		const int start = cols * row;
		const int end   = start + cols - 1;

		for (int i = start + 1; i < end; ++i) {
			diff = b[i] - stencil[0] * u[i-1] - stencil[1] * u[i] - stencil[2] * u[i+1] - stencil[3] * u[i+cols] - stencil[4] * u[i-cols];
			acc += diff * diff;
		}
	}

	return std::sqrt(acc);
}
