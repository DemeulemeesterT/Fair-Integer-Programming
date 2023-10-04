#include "Ufl.h"

void Ufl::generate_instance(Ufl_parameters U, bool print) {
	// This is based on the code to generate Bilde-Krarup benchmark instances, included in UflLib
		// https://resources.mpi-inf.mpg.de/departments/d1/projects/benchmarks/UflLib/BildeKrarup.html

	int i, j;

	m_customers = U.m;
	n_locations = U.n;
	int F_LO = U.f_min;
	int F_HI = U.f_max;
	int C_LO = U.c_min;
	int C_HI = U.c_max;
	c_alloc_max = U.c_max;

	// Initialize vectors f and c
	f = std::vector<double>(n_locations, -1);
	c = std::vector<std::vector<double>>(n_locations, std::vector<double>(m_customers, -1));

	// Initate random generators
	std::mt19937 generator((unsigned)U.seed);
	std::uniform_int_distribution<int> distr_f(F_LO, F_HI);
	std::uniform_int_distribution<int> distr_c(C_LO, C_HI);

	for (i = 0; i < n_locations; i++) {
		f[i] = distr_f(generator);
		for (j = 0; j < m_customers; j++) {
			c[i][j] = distr_c(generator);
		}
	}
	printf("Pause");
}


Ufl::Ufl(Ufl_parameters Ufl_input, bool print) {
	generate_instance(Ufl_input, print);


}


Ufl::~Ufl() {

}

