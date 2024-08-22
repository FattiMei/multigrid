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

	for (auto op : fcycle_alternative) {
		std::cout << op << std::endl;
	}

	assert(equal(fcycle_reference, fcycle_alternative));

	return 0;
}
