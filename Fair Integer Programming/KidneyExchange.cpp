#include "KidneyExchange.h"

void KidneyExchange::read_data(std::string name, bool print) {
	std::ifstream StFile;
	std::string full_name = "Kidney exchange instances/" + name + ".INPUT";
	StFile.open(full_name);

	// If instance can't be found, generate it.
	if (!StFile) {
		std::cout << "\n\n This is embarrassing, I can't find that file...\n";
		char stop;
		std::cin >> stop;

		return;
	}

	StFile >> n_pairs;
	StFile >> n_arcs;

	// Initialize 'A'
	A = std::vector<std::vector<bool>>(n_pairs, std::vector<bool>(n_pairs, 0));

	int a1, a2; // The agents that we will read from the file
	int not_needed;
	for (int i = 0; i < n_arcs; i++) {
		StFile >> a1; 
		StFile >> a2; 
		StFile >> not_needed;

		// Add this arc to our graph
		A[a2][a1] = 1;
	}

	data_name = name;

	StFile.close();

}

void KidneyExchange::export_to_dat(inst I, bool print) {
	std::ofstream S;
	std::string export_name = "Data instances/KE/" + I.data_name + ".dat";
	S.open(export_name);
	S << I.m << " ";
	S << I.n << " ";
	S << I.t << " ";
	if (I.Y_bool == true) {
		S << "B" << "\n";
	}
	else {
		S << "I" << "\n";
	}
	for (int i = 0; i < I.n_var; i++) {
		S << I.v[i] << " ";
	}
	S << "\n";
	for (int l = 0; l < I.m; l++) {
		S << I.C[l] << " ";
	}
	S << "\n";
	for (int l = 0; l < I.m; l++) {
		for (int i = 0; i < I.n_var; i++) {
			S << I.A[l][i] << " ";
		}
		S << "\n";
	}

	S.close();
}

void KidneyExchange::find_cycles(bool print) {
	std::vector<bool> empty(n_pairs, 0);

	// Go through all agents
	for (int i = 0; i < n_pairs; i++) {
		for (int j = i; j < n_pairs; j++) {
			if (A[i][j] == 1) { // Edge (i,j) is in the graph
				if (A[j][i] == 1) { // Cycle of length two
					C.push_back(empty);
					C.back()[i] = 1;
					C.back()[j] = 1;
				}
				else {
					for (int k = 0; k < n_pairs; k++) {
						if (A[j][k] == 1) {
							if (A[k][i] == 1) {
								C.push_back(empty);
								C.back()[i] = 1;
								C.back()[j] = 1;
								C.back()[k] = 1;
							}
						}
					}
				}
			}
		}
	}
}

inst KidneyExchange::generate_instance(bool export_inst, bool print) {
	inst I;

	I.n = n_pairs; // One X-variable for each pair
	I.t = C.size(); // One Y-variable for each cycle
	I.n_var = I.n + I.t;
	I.n_var = I.n + I.t;
	I.Y_bool = true;
	I.Y_coeff_zero = true;

	I.m = 2 * I.n; // Twice one constraint for each agent

	// Objective coefficients
	for (int i = 0; i < I.n; i++) {
		I.v.push_back(1);
	}
	for (int i = 0; i < I.t; i++) {
		I.v.push_back(0);
	}

	// Capacities of the constraints
	for (int i = 0; i < I.n; i++) {
		I.C.push_back(1);
	}
	for (int i = 0; i < I.n; i++) {
		I.C.push_back(0);
	}

	// Constraint matrix
		// The constraints are, in this order:
			// \sum_{c:i\in c} y_c \leq 1, \forall i \in A (each pair in at most one selected cycle)
			// x_i - \sum_{c:i\in c}y_c \leq 0, \forall i \in A (fix x_i's based on Y-variables)

	// Initialize
	I.A = std::vector<std::vector<double>>(I.m, std::vector<double>(I.n_var, 0.0));

	// Fill in
	// First go through all cycles:
	int x_index, y_index;
	for (int c = 0; c < I.t; c++) {
		for (int j = 0; j < I.n; j++) {
			if (C[c][j] == 1) { // If agent j is in cycle c
				x_index = I.n + j;
				y_index = I.n + c;
					// This is added to avoid C++ warnings

				// Add in first set of constraints
				I.A[j][y_index] = 1;

				// Add in second set of constraints
				I.A[x_index][y_index] = -1;
			}
		}
	}
	for (int i = 0; i < I.n; i++) {
		x_index = I.n + i;
		I.A[x_index][i] = 1;
	}

	I.data_name = data_name;

	if (export_inst) {
		export_to_dat(I, print);
	}

	return I;
}

KidneyExchange::KidneyExchange(std::string name, bool print) {
	read_data(name, print);
	find_cycles(print);
}

KidneyExchange::~KidneyExchange() {

}
