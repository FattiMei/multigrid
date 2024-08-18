#ifndef __PROBLEM_HPP__
#define __PROBLEM_HPP__


#include "operator.hpp"
#include <vector>


enum class InitializationStrategy {
	Zeros,
	Lerp
};


class Problem {
	public:
		Problem(const int nodes) : n(nodes), rhs(nodes) {};
		virtual ~Problem() = default;

		virtual DiscreteOperator* get_discrete_operator(const int level = 0) const = 0;
		int get_size() const { return n; };
		const double* get_rhs() const { return static_cast<const double*>(rhs.data()); };


	protected:
		const int n;
		std::vector<double> rhs;
};



#endif
