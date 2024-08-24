#ifndef __PERIODIC_HPP__
#define __PERIODIC_HPP__


#include "problem.hpp"
#include <functional>


class ForcedRC : public Problem {
	public:
		ForcedRC(
			const double period,
			const double tau,
			const int n,
			const std::function<double(double)> signal
		);
		~ForcedRC() = default;

		double get_step() const override { return h; };
		DiscreteOperator* get_discrete_operator(const int level = 0) const override;
		const std::vector<double>& get_mesh() const;


	private:
		const double tau;
		const double h;
		std::vector<double> mesh;
};


#endif
