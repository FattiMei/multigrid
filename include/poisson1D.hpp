#ifndef __POISSON1D_HPP__
#define __POISSON1D_HPP__


#include <functional>


enum class InitializationStrategy {
	Zeros,
	Lerp
};


class Poisson1D {
	public:
		Poisson1D(double inf, double sup, int n, std::function<double(double)> f, std::pair<double,double> boundary);
		~Poisson1D();

		void set_initial_approximation(double *u, InitializationStrategy strategy);
		int  get_problem_size();
		const double *get_rhs();



	private:
		const int n;
		double *b;
};


#endif
