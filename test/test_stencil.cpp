#include <iostream>
#include "utils.hpp"
#include "stencil.hpp"


int main() {
	const int    n = 10;
	const double h = 1.0 / (n - 1.0);

	const auto b = generate_random_vector(n, 0.0, 10.0);
	std::vector<double> rnaive(n);
	std::vector<double> rgood(n);

	auto u = generate_random_vector(n, 0.0, 10.0);
	u[0] = b[0];
	u[n-1] = b[n-1];

	const auto naive = NaiveThreePointStencil(
		n,
		{-1.0 / (h*h), 2.0 / (h*h), -1.0 / (h*h)}
	);

	const auto good = ThreePointStencil(
		n,
		h,
		{-1.0, 2.0, -1.0}
	);

	std::cout << "||r|| naive formula: " << naive.compute_residual_norm(b.data(), u.data()) << std::endl;
	std::cout << "||r|| good  formula: " << good.compute_residual_norm(b.data(), u.data()) << std::endl;

	naive.compute_residual(b.data(), u.data(), rnaive.data());
	good.compute_residual(b.data(), u.data(), rgood.data());

	std::cout << "l2   diff: " << l2norm  (rnaive, rgood) << std::endl;
	std::cout << "linf diff: " << linfnorm(rnaive, rgood) << std::endl;

	return 0;
}
