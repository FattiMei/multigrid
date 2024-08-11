#include "multigrid.hpp"


// all these implementations assume that pointers have sufficient size
void injective_restriction(const int n, const double src[], double dest[]) {
	for (int i = 0, write_index = 0; i < n; i += 2) {
		dest[write_index++] = src[i];
	}
}


void full_weight_restriction(const int n, const double src[], double dest[]) {
	int write_index = 0;
	dest[write_index++] = src[0];

	for (int i = 2; i < n-1; i += 2) {
		dest[write_index++] = 0.25 * src[i-1] + 0.5 * src[i] + 0.25 * src[i+1];
	}

	dest[write_index] = src[n-1];
}


void linear_prolongation(const int m, const double src[], double dest[]) {
	int write_index = 0;

	for (int i = 0; i < m-1; ++i) {
		dest[write_index++] = src[i];
		dest[write_index++] = 0.5 * src[i] + 0.5 * src[i+1];
	}

	dest[write_index] = src[m-1];
}
