#include <iostream>
#include "multigrid.hpp"


bool equal(const std::vector<MgOp> &reference, const std::vector<MgOp> &alternative) {
	if (reference.size() != alternative.size()) {
		std::cerr << "Recipes have different size" << std::endl;
		return false;
	}

	for (size_t i = 0; i < reference.size(); ++i) {
		if (reference[i] != alternative[i]) {
			return false;
		}
	}

	return true;
}


int main() {
	const std::vector<MgOp> vcycle_reference {
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::DirectSolve,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Prolong,
		MgOp::Relax,
	};
	const auto vcycle_alternative = MgCycle::V(3);
	assert(equal(vcycle_reference, vcycle_alternative));


	const std::vector<MgOp> fcycle_reference {
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::DirectSolve,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::DirectSolve,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::DirectSolve,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Prolong,
		MgOp::Relax,
	};
	const auto fcycle_alternative = MgCycle::F(4);
	assert(equal(fcycle_reference, fcycle_alternative));


	const std::vector<MgOp> wcycle_reference {
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::DirectSolve,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::DirectSolve,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::DirectSolve,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::DirectSolve,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Prolong,
		MgOp::Relax,
		MgOp::Prolong,
		MgOp::Relax
	};
	const auto wcycle_alternative = MgCycle::W(4);
	assert(equal(wcycle_reference, wcycle_alternative));

	return 0;
}
