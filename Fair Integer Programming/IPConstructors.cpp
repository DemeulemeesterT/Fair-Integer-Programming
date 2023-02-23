#include "IPSolver.h"

inst random_generator(bool print) {
	inst I;
	printf("\t Number of constraints:\n\t\t");
	std::cin >> I.m;
	printf("\t Number of agents/X-variables: \n\t\t");
	std::cin >> I.n;
	printf("\t Number of auxiliary Y-variables: \n\t\t");
	std::cin >> I.t;
	I.n_var = I.n + I.t;
	I.Y_bool = false;
	if (I.t > 0) {
		std::string input = "K";
		while ((input != "b") && (input != "i")) {
			printf("\t Press i for integer Y-variables, press b for binary Y-variables:\n\t\t");
			std::cin >> input;
			if ((input != "b") && (input != "i")) {
				printf("\t\t That is not a valid input. Press i or b.\n\n");
			}
		}
		if (input == "b") {
			I.Y_bool = true;
		}
		else {
			I.Y_bool = false;
		}
		while ((input != "z") && (input != "p")) {
			printf("\t Press z for zero-coefficients of the Y-variables in the objective function, press p otherwise: \n\t\t");
			std::cin >> input;
			if ((input != "z") && (input != "p")) {
				printf("\t\t That is not a valid input. Press z or p.\n\n");
			}
		}
		if (input == "z") {
			I.Y_coeff_zero = true;
		}
		else {
			I.Y_coeff_zero = false;
		}
	}
	else {
		I.Y_bool = true;
		I.Y_coeff_zero = true;
	}

	printf("\t Random seed generator:\n\t\t");
	int seed;
	std::cin >> seed;

	// GENERATION (change here if wanted)
	// In default, we picked parameters that are not too far from each other to create multiple optimal solutions
	std::mt19937 generator((unsigned)seed);
	std::uniform_int_distribution<int> distr_v(1, 3);
	std::uniform_int_distribution<int> distr_A(0, 2);

	// Objective coefficients.
	if (I.Y_coeff_zero == false) {
		for (int i = 0; i < I.n_var; i++) {
			I.v.push_back(distr_v(generator));
		}
	}
	else {
		for (int i = 0; i < I.n; i++) {
			I.v.push_back(distr_v(generator));
		}
		for (int i = 0; i < I.t; i++) {
			I.v.push_back(0);
		}
	}

	// Coefficient matrix
	std::vector<double> empty(I.n_var, 0);
	for (int l = 0; l < I.m; l++) {
		I.A.push_back(empty);
		for (int i = 0; i < I.n_var; i++) {
			I.A[l][i] = distr_A(generator);
		}
	}

	// Capacities of the knapsacks
	for (int l = 0; l < I.m; l++) {
		// Make the constraint somewhat limiting
		// By letting it take a value in the interval [1, sum_i(A_li) - 1]
		int sum = 0;
		for (int i = 0; i < I.n_var; i++) {
			sum += I.A[l][i];
		}
		std::uniform_int_distribution<int> distr_C(1, sum - 1);
		I.C.push_back(distr_C(generator));
	}

	// EXPORT
	std::ofstream S;
	S.open("Data instances/Generated.dat");
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

	printf("\n\tYour instance has been generated and exported to 'Generated.dat'. \n");
	printf("\t\t(If you want to change the generation process, modify the 'random_generator' function in 'IPConstructors.cpp')\n\n");

	return I;
}

inst read_data(std::string name) {
	std::ifstream StFile;
	std::string full_name = "Data instances/" + name + ".dat";
	StFile.open(full_name);

	inst I;

	// If instance can't be found, generate it.
	if (!StFile) {
		std::cout << "\n\n This is embarrassing, I can't find that file...\n";
		std::cout << "\t Instead, we will generate a file randomly:\n\n";

		I = random_generator(false);
		return I;
	}

	StFile >> I.m;
	StFile >> I.n;
	StFile >> I.t;
	I.n_var = I.n + I.t;
	std::string bool_int;
	StFile >> bool_int;
	if (bool_int == "B") {
		I.Y_bool = true;
	}
	else if (bool_int == "I") {
		I.Y_bool = false;
	}
	else {
		std::printf("\n\n\tThere is a mistake in the input file.\n\tPlease clearly state whether the Y-variables are integer or boolean by writing I or B on the end of the first line of the .dat file.\n\tNo other entries than I or B are allowed.\n\n");
		std::printf("\tWe assume the use of integer variables.\n\n");
		I.Y_bool = false;
	}
	// weights of the objects
	double dummy;
	for (int i = 0; i < I.n_var; i++) {
		StFile >> dummy;
		I.v.push_back(dummy);
	}
	// capacities of the knapsacks
	for (int i = 0; i < I.m; i++) {
		StFile >> dummy;
		I.C.push_back(dummy);
	}
	// A
	std::vector<double> empty(I.n_var, 0);
	for (int i = 0; i < I.m; i++) {
		I.A.push_back(empty);
		for (int j = 0; j < I.n_var; j++) {
			StFile >> dummy;
			I.A[i][j] = dummy;
		}
	}

	I.data_name = name;

	// Check if all Y-variables have coefficient zero in the objective function
	I.Y_coeff_zero = true;
	for (int i = 0; i < I.t; i++) {
		if (I.v[I.n + i] != 0) {
			I.Y_coeff_zero = false;
			i = I.t;
		}
	}

	StFile.close();

	return I;
}

void IPSolver::read_MIPLIB_info(bool print) {
	// The file has to be located inside the MIPLIB instances folder in order to read the file name correctly!

	std::ifstream StFile;
	StFile.open("MIPLIB instances/miplib2017-v22.solu");
	std::string line;
	MIPLIB_instance M;

	if (!StFile) {
		std::cout << "\n\n This is embarrassing, I can't find that file...\n";
	}

	while (std::getline(StFile, line)) {
		std::istringstream issLine(line);
		issLine >> M.solution_info >> M.name;

		if (M.solution_info != "=unkn=") {
			issLine >> M.obj_val;
		}
		else {
			M.obj_val = NAN;
		}
		MIPLIB_summary.push_back(M);
	}

	StFile.close();
}

void IPSolver::initializeGurobi(bool print) {
	try {
		env = new GRBEnv();
		env->set(GRB_StringParam_LogFile, "Generated Log-files/IPTests.log");
		model = new GRBModel(*env);
		model->set(GRB_IntParam_PoolSolutions, 10); // Limit the number of solutions that will be stored.
		model->set(GRB_DoubleParam_PoolGap, 0.00001); // Limit the search space by setting a gap for the worst possible solution that will be accepted

		X = new GRBVar[I.n];

		// Build
		for (int i = 0; i < I.n; i++) {
			char name_x[13];
			sprintf_s(name_x, "x_%i", i);
			X[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name_x);
		}

		Y_var = new GRBVar[I.t];
		for (int i = 0; i < I.t; i++) {
			char name_y[13];
			sprintf_s(name_y, "y_%i", i);
			if (I.Y_bool == true) {
				Y_var[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name_y);
			}
			else {
				Y_var[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, name_y);
			}
		}

		GRBLinExpr obj;
		obj = 0.0;

		// Objective function
		for (int i = 0; i < I.n; i++) {
			obj += I.v[i] * X[i];
		}
		for (int i = 0; i < I.t; i++) {
			obj += I.v[I.n + i] * Y_var[i];
		}
		model->setObjective(obj, GRB_MAXIMIZE);

		// Constraints for each knapsack
		GRBLinExpr con;
		for (int i = 0; i < I.m; i++) {
			con = 0;
			for (int j = 0; j < I.n; j++) {
				con += I.A[i][j] * X[j];
			}
			for (int j = 0; j < I.t; j++) {
				con += I.A[i][I.n + j] * Y_var[j];
 			}
			model->addConstr(con <= I.C[i]);
		}

		//model->write("Generated Formulations/IPModel.lp");

		//model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

		env_identical = new GRBEnv();
		env_identical->set(GRB_StringParam_LogFile, "Generated Log-files/IPIdentical.log");
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Error during optimization" << std::endl;
	}
}

void IPSolver::initializeVariables(inst I_in, bool print) {
	I = inst();
	I.m = I_in.m;
	I.n = I_in.n;
	I.t = I_in.t;
	I.n_var = I_in.n_var;
	I.Y_bool = I_in.Y_bool;
	I.Y_coeff_zero = I_in.Y_coeff_zero;
	I.data_name = I_in.data_name;
	opt = -112233445566;
	for (int i = 0; i < I.n_var; i++) {
		I.v.push_back(I_in.v[i]);
	}
	std::vector<double> empty(I.n_var, 0);
	for (int j = 0; j < I.m; j++) {
		I.C.push_back(I_in.C[j]);
		I.A.push_back(empty);
		for (int i = 0; i < I.n_var; i++) {
			I.A[j][i] = I_in.A[j][i];
		}
	}

	Y = std::vector<int>(I.n, -1);
	M = std::vector<int>(I.n, -1);
	N = std::vector<int>(I.n, -1);

	M_size = 0;

	// Initialize the matrix 'identical'
	identical = std::vector<std::vector<int>>(I.n, std::vector<int>(I.n, -1));
	identical_clustered = std::vector<int>(I.n, -1);

	solver_times = 0;

	S = std::vector<solution>();
	S_greedy = std::vector<solution>();

	greedy_wrongly_classified = std::vector<int>();

	max_Y_value = 10;
	max_Y_increment = 5;
	done_partition = false;
	done_greedy_partition = false;
	done_find_identical = false;
	done_solve = false;

	time_partition = 0;
	time_greedy_partition = 0;

	number_of_agents_at_once_RSD = 19;

	read_MIPLIB_info(false);
}

IPSolver::IPSolver(bool print) {
	printf("\n\n Which file do you want to read? (Press r + ENTER to enter the random generator module)\n\t ");
	std::string name;
	std::cin >> name;
	printf("\n\n");
	inst I_in;
	if (name == "r") {
		I_in = random_generator(true);
	}
	else {
		I_in = read_data(name);
	}

	initializeVariables(I_in, print);
	initializeGurobi(print);
	
}

IPSolver::IPSolver(inst I_in, bool print) {
	initializeVariables(I_in, print);
	initializeGurobi(print);
}

IPSolver::IPSolver(IPSolver* K_in, bool print) {
	I = inst();
	I.m = K_in->I.m;
	I.n = K_in->I.n;
	I.t = K_in->I.t; 
	I.n_var = K_in->I.n_var;
	I.Y_bool = K_in->I.Y_bool;
	I.Y_coeff_zero = K_in->I.Y_coeff_zero;
	opt = K_in->opt;

	solver_times = K_in->solver_times;
	for (int i = 0; i < I.n_var; i++) {
		I.v.push_back(K_in->I.v[i]);
	}
	std::vector<double> empty(I.n_var, 0);
	for (int j = 0; j < I.m; j++) {
		I.C.push_back(K_in->I.C[j]);
		I.A.push_back(empty);
		for (int i = 0; i < I.n_var; i++) {
			I.A[j][i] = K_in->I.A[j][i];
		}
	}

	Y = K_in->Y;
	M = K_in->M;
	N = K_in->N;

	M_size = K_in->M_size;

	identical = std::vector<std::vector<int>>(I.n, std::vector<int>(I.n, -1));
	for (int i = 0; i < I.n; i++) {
		identical.push_back(K_in->identical[i]);
	}
	identical_clustered = K_in->identical_clustered;

	S = std::vector<solution>();
	for (int i = 0; i < K_in->S.size(); i++) {
		S.push_back(K_in->S[i]);
	}

	S_greedy = std::vector<solution>();
	for (int i = 0; i < K_in->S_greedy.size(); i++) {
		S_greedy.push_back(K_in->S_greedy[i]);
	}

	greedy_wrongly_classified = K_in->greedy_wrongly_classified;

	max_Y_value = K_in->max_Y_value;
	max_Y_increment = K_in->max_Y_increment;

	done_partition = K_in->done_partition;
	done_greedy_partition = K_in->done_greedy_partition;
	done_find_identical = K_in->done_find_identical;
	done_solve = K_in->done_solve;

	time_partition = K_in->time_partition;
	time_greedy_partition = K_in->time_greedy_partition;

	number_of_agents_at_once_RSD = K_in->number_of_agents_at_once_RSD;

	initializeGurobi(print);
}

IPSolver::IPSolver(const std::string& file_name, bool print) {
	// Create a different model to first read the model
	// and then transfer it to the model in which our notation is used
	try {
		env = new GRBEnv();
		model = new GRBModel(*env, file_name);

		if (model->get(GRB_IntAttr_IsMIP) == 0) {
			throw GRBException("Model is not a MIP");
		}
		if (print) {
			//model->write("Generated Formulations/sorrell3.lp");
		}

		// Extract everything from this model, and place it in the correct format
		GRBVar* vars = model->getVars();
		GRBConstr* cons = model->getConstrs();

		int numvars = model->get(GRB_IntAttr_NumVars);
		int numcons = model->get(GRB_IntAttr_NumConstrs);
		int modelsense = model->get(GRB_IntAttr_ModelSense); // 1 for minimization, -1 for maximization
		for (int i = 0; i < numcons; i++) {
			std::string name = cons[i].get(GRB_StringAttr_ConstrName);
			//std::cout << name << "\n";
		}

		// Store which variables belong to X and which to Y
		std::vector<bool> belongs_to_X(numvars, false);

		for (int i = 0; i < numvars; i++) {
			char v_type = vars[i].get(GRB_CharAttr_VType);
			int UB = vars[i].get(GRB_DoubleAttr_UB);
			int LB = vars[i].get(GRB_DoubleAttr_LB);
			double obj_coeff = vars[i].get(GRB_DoubleAttr_Obj);

			// Check if objective coefficient is integer:
			double intpart;
			if (std::modf(obj_coeff, &intpart) != 0.0) {
				throw GRBException("One of the objective coefficients is not integer.");
			}

			int index = vars[i].index();
			// Does this variable belong to X or to Y?
			if (v_type == 'B' || v_type == 'I') { // Boolean or Integer
				if (UB == 1) {
					if (LB == 0) {
						// Is objective coefficient integer and > 0 (for maximization problem)?
						double intpart;
						if (std::modf(obj_coeff, &intpart) == 0.0) {
							if (modelsense == 1) {
								// Minimization, so the value should be negative
								if (obj_coeff < 0) {
									belongs_to_X[i] = true;
								}
							}
							else if (modelsense == -1) {
								if (obj_coeff > 0) {
									belongs_to_X[i] = true;
								}
							}
						}
					}
				}
			}
			else {
				// Then we have a variable that is continuous, this is not modeled in our setting
				throw GRBException("One of the variables is not binary or integer.");
			}
		}

		// Store pointers to the variables in X and Y
		X = new GRBVar[numvars];
		Y_var = new GRBVar[numvars];

		int x_counter = 0;
		int y_counter = 0;
		for (int i = 0; i < numvars; i++) {
			if (belongs_to_X[i] == true) {
				X[x_counter] = vars[i];
				x_counter++;
			}
			else {
				Y_var[y_counter] = vars[i];
				y_counter++;
			}
		}

		// Initialize the instance object
		inst I_in;
		I_in.n = x_counter;
		I_in.t = y_counter;
		I_in.m = numcons;
		I_in.data_name = file_name;
		I_in.n_var = numvars;

		for (int i = 0; i < I_in.n_var; i++) {
			I_in.v.push_back(vars[i].get(GRB_DoubleAttr_Obj));
		}
		std::vector<double> empty(I_in.n_var, 0);
		for (int j = 0; j < I_in.m; j++) {
			char sense = cons[j].get(GRB_CharAttr_Sense);
			I_in.A.push_back(empty);

			// We model all constraints as <=
			// If >=: we change direction
			// If =: we could in theory add an additional Y-variable with a zero coefficient in the objective function to turn this into a <= constraint
			if (sense == '<' || sense == '=') {
				I_in.C.push_back(cons[j].get(GRB_DoubleAttr_RHS));
				// First of the variables in X, then of the variables in Y
				for (int i = 0; i < I_in.n; i++) {
					I_in.A[j][i] = model->getCoeff(cons[j], X[i]);
				}
				for (int i = 0; i < I_in.t; i++) {
					I_in.A[j][i + I_in.n] = model->getCoeff(cons[j], Y_var[i]);
				}
			}
			else if (sense == '>') {
				I_in.C.push_back(-cons[j].get(GRB_DoubleAttr_RHS));
				// First of the variables in X, then of the variables in Y
				for (int i = 0; i < I_in.n; i++) {
					I_in.A[j][i] = - model->getCoeff(cons[j], X[i]);
				}
				for (int i = 0; i < I_in.t; i++) {
					I_in.A[j][i + I_in.n] = - model->getCoeff(cons[j], Y_var[i]);
				}
			}
			else {
				throw GRBException("One of the constraints is not >, <, or =.");
			}
		}

		I_in.Y_bool = true;
		I_in.Y_coeff_zero = true;
		for (int i = 0; i < I_in.t; i++) {
			if (Y_var[i].get(GRB_CharAttr_VType) != 'B') {
				if (Y_var[i].get(GRB_CharAttr_VType) == 'I') {
					I_in.Y_bool = false;
				}
				else {
					throw GRBException("One of the variables is not binary or integer.");
				}
			}
			
			if (Y_var[i].get(GRB_DoubleAttr_Obj) != 0.0) {
				I_in.Y_coeff_zero = false;
			}
		}

		initializeVariables(I_in, print);

		delete[] cons;
		delete[] vars;

		// Check for the MIPLIB instances what the optimal solution is and store it in 'opt'.
		// First, remove the unnecessary parts of the 'file_name'
		if (file_name.rfind("MIPLIB instances/", 0) == 0) { // pos=0 limits the search to the prefix
			std::string t = file_name;
			std::string s = "MIPLIB instances/";
			std::string::size_type i = file_name.find(s);
			

			if (i != std::string::npos)
				t.erase(i, s.length());

			std::string u = ".mps";
			std::string::size_type j = t.find(u);
			if (j != std::string::npos)
				t.erase(j, u.length());

			// Find whether the optimal value is known
			int k = 0;
			while (k < MIPLIB_summary.size()) {
				if (MIPLIB_summary[k].name == t) {
					if (MIPLIB_summary[k].solution_info == "=opt=") {
						opt = MIPLIB_summary[k].obj_val;
					}
					k = MIPLIB_summary.size();
				}
				else {
					k++;
				}
			}
		}

	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Error during optimization" << std::endl;
	}

	
}

IPSolver::~IPSolver() {
	delete[] X;
	delete[] Y_var;
	delete model;
	delete env;

	delete env_identical;

}
