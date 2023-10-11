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
	filename = U.filename;

	// Initialize vectors f and c
	f = std::vector<double>(n_locations, -1);
	c = std::vector<std::vector<double>>(n_locations, std::vector<double>(m_customers, -1));

	// Initate random generators
	std::mt19937 generator((unsigned)U.seed);
	std::uniform_int_distribution<int> distr_f(F_LO, F_HI);
	std::uniform_int_distribution<int> distr_c(C_LO, C_HI);

	for (i = 0; i < n_locations; i++) {
		f[i] = distr_f(generator);
	}

	// You can not have more indifference classes than the difference
	if (U.nr_indiff_classes > U.c_max - U.c_min) {
		U.nr_indiff_classes = U.c_max - U.c_min;
	}

	// Randomly sample integers
	if (U.nr_indiff_classes == -1) {
		for (i = 0; i < n_locations; i++) {
			for (j = 0; j < m_customers; j++) {
				c[i][j] = distr_c(generator);
			}
		}
	}

	else {
		double step_size = (U.c_max - U.c_min) / U.nr_indiff_classes;

		int generated_c;
		for (i = 0; i < n_locations; i++) {
			for (j = 0; j < m_customers; j++) {
				generated_c = distr_c(generator);

				double test_value = U.c_min;
				while (test_value <= generated_c) {
					test_value = test_value + step_size;
					if (test_value > generated_c) {
						c[i][j] = (int)test_value - step_size;
					}
				}
			}
		}
	}
	
	
}

inst Ufl::generate_formulation(bool export_inst, bool print) {
	inst I;

	I.n = m_customers; // One X-variable for each customer
	I.t = n_locations + n_locations * m_customers; // One Y-variable for location, and one Y-variable for each customer-location pair
		// The Y-variables are ordered as follows:
			// First the n*m variables for each customer-location pair
			// Then the n variables that capture whether a facility is opened in that location
		
	I.n_var = I.n + I.t;
	I.Y_bool = true;
	I.Y_coeff_zero = false;

	I.X_bool = false;
	I.X_integer = false;

	I.m = 4*m_customers + m_customers * n_locations; 
		// Two constraints for each customer (equality constraints, so 2*m_customers constraints for each set of constraints),
		// and one constraint for each customer-location pair

	// ORDER of the constraints
		// \sum_i M_{ij} = 1, for each customer j (each customer allocated to at most one facility)
			// Equality constraints, implemented as two sets of inequalities.
		// M_{ij} <= Y_i, for each customer-location pair (j,i) (to calculate fixed costs in objective function)
		// X_j = c_alloc_max - \sum_i c_{ij} M_{ij}, for each customer j (utility of an agent) 
			// Equality constraints, implemented as two sets of inequalities.

 		

	// Objective coefficients
	for (int i = 0; i < I.n; i++) {
		I.v.push_back(0); // X-variables, utilities of the customers
	}
	for (int i = 0; i < n_locations; i++) {
		for (int j = 0; j < m_customers; j++) {
			I.v.push_back(c[i][j]); // One for each customer-location pair
		}
	}
	for (int i = 0; i < n_locations; i++) {
		I.v.push_back(f[i]); // For each location
	}

	// Capacities of the constraints (order of the constraints, see above)
	for (int i = 0; i < m_customers; i++) {
		I.C.push_back(1);
		I.C.push_back(-1);
	}
	for (int i = 0; i < n_locations; i++) {
		for (int j = 0; j < m_customers; j++) {
			I.C.push_back(0);
		}
	}
	for (int i = 0; i < m_customers; i++) {
		I.C.push_back(c_alloc_max);
		I.C.push_back(-c_alloc_max);
	}

	// Constraint matrix
	// Initialize
	I.A = std::vector<std::vector<double>>(I.m, std::vector<double>(I.n_var, 0.0));

	// 'counter' keeps track of which row in A we are editing
	int counter = 0;
	// Fill in
	// First, go through all agents
	for (int j = 0; j < m_customers; j++) {
		for (int i = 0; i < n_locations; i++) {
			// The variables are ordered: M_11, M_12, M_13, ..., M_1m, M_21, M_22, ...
			int y_index = m_customers + m_customers * i + j;

			I.A[counter][y_index] = 1;
		}
		counter++;

		for (int i = 0; i < n_locations; i++) {
			// The variables are ordered: M_11, M_12, M_13, ..., M_1m, M_21, M_22, ...
			int y_index = m_customers + m_customers * i + j;

			I.A[counter][y_index] = -1;
		}
		counter++;
	}

	// Second, go through all customer-location pairs
	int variable_counter = m_customers;
	int Y_counter;
	for (int i = 0; i < n_locations; i++) {
		Y_counter = m_customers + m_customers * n_locations + i;
		for (int j = 0; j < m_customers; j++) {
			I.A[counter][variable_counter] = 1;
			variable_counter++;

			I.A[counter][Y_counter] = -1;
			counter++;
		}
	}

	// Third, go through all customers again
	variable_counter = m_customers;
	for (int j = 0; j < m_customers; j++) {
		I.A[counter][j] = 1; // X-variable
		for (int i = 0; i < n_locations; i++) {
			int y_index = m_customers + m_customers * i + j;
			I.A[counter][y_index] = 1;
		}
		counter++;

		I.A[counter][j] = -1; // X-variable
		for (int i = 0; i < n_locations; i++) {
			int y_index = m_customers + m_customers * i + j;
			I.A[counter][y_index] = -1;
		}
		counter++;
	}

	I.data_name = filename;

	if (export_inst) {
		export_to_dat(I, print);
	}

	return I;
}

void Ufl::read_data_OR_lib(std::string name, bool print) {
	std::ifstream StFile;
	std::string full_name = "Uncapacitated Facility Location instances/" + name + ".txt";
	StFile.open(full_name);

	// If instance can't be found, generate it.
	if (!StFile) {
		std::cout << "\n\n This is embarrassing, I can't find that file...\n";
		char stop;
		std::cin >> stop;

		return;
	}

	StFile >> n_locations;
	StFile >> m_customers;

	// Initialize vectors
	c = std::vector<std::vector<double>>(n_locations, std::vector<double>(m_customers, 0));
	f = std::vector<double>(n_locations, 0);

	double cap, demand;

	// Read info locations
	for (int j = 0; j < n_locations; j++) {
		StFile >> cap;
		StFile >> f[j];
	}

	// Read info customers
	int counter = 0;
	for (int i = 0; i < m_customers; i++) {
		StFile >> demand;
		for (int j = 0; j < n_locations; j++) {
			StFile >> c[j][i];
		}
		counter++;
	}

	filename = name;

	StFile.close();

}

void Ufl::export_to_dat(inst I, bool print) {
	std::ofstream S;
	std::string export_name = "Data instances/UFL/" + I.data_name + ".dat";
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

Ufl::Ufl(Ufl_parameters Ufl_input, bool print) {
	generate_instance(Ufl_input, print);
}

Ufl::Ufl(std::string filename, bool print) {
	read_data_OR_lib(filename, print);
}



Ufl::~Ufl() {

}

