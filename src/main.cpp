#include <iostream>
#include "poisson1D.hpp"


int main() {
	Poisson1D problem(
		0.0,
		1.0,
		100,
		[](double x) {return x;},
		{0.0, 0.0}
	);

	std::cout << "Hello, World\n" << std::endl;
}
