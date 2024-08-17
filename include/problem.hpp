#ifndef __PROBLEM_HPP__
#define __PROBLEM_HPP__


#include "operator.hpp"


class Problem {
	public:
		Problem(const int nodes) : n(nodes) {};
		virtual ~Problem() = default;

		int size() { return n; };
		virtual DiscreteOperator* get_discrete_operator(const int level = 0) = 0;


	protected:
		const int n;
};



#endif
