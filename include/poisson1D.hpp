#ifndef __POISSON1D_HPP__
#define __POISSON1D_HPP__


#include <functional>


enum class InitializationStrategy {
	Zeros,
	Lerp
};


using Update = std::function<void(int, const double *, const double *, double *)>;


class Poisson1D {
	public:
		Poisson1D(double inf, double sup, int n, std::function<double(double)> f, std::pair<double,double> boundary);
		~Poisson1D();

			void	set_initial_approximation(double *u, InitializationStrategy strategy) const;
			int	get_problem_size() const;
			double*	get_rhs() const;
		const	Update	get_iteration_formula(const int level = 0) const;
		const	Update	get_residual_formula (const int level = 0) const;
			double	get_residual_norm(const double u[]) const;


	private:
		const int n;
		const double h;

		// @DESIGN: once b is initialized, it will be constant. How do we solve this problem?
		double *b;
};


#endif
