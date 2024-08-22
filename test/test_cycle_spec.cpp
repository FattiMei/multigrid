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

	const auto alternative = MgCycle::V(2);

	assert(equal(vcycle_reference, alternative));

	return 0;
}
