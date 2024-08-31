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
		int get_dimension(const int dim) const override {
			return dim == 0 ? n : 0;
		};
		DiscreteOperator* get_discrete_operator(const int level = 0) const override;
		void set_initial_approximation(double u[], const InitializationStrategy strategy) const override;

		const std::vector<double>& get_mesh() const;

	private:
		const double tau;
		const double h;
		std::vector<double> mesh;
};


#endif
