#include "LotteryDesigner.h"

//******************************************
//********* AUXILIARY FUNCTIONS ************
//******************************************

void LotteryDesigner::print_lottery(lottery L) {
	printf("\nThe selection probabilities are:\n");
	for (int i = 0; i < K->I.n; i++) {
		printf("\tAgent %2.i", i);
		if (K->Y[i] == 1) {
			printf(" (Y)");
		}
		else if (K->N[i] == 1) {
			printf(" (N)");
		}
		else {
			printf(" (M)");
		}
		printf(":   %.3f\n", L.p[i]);
	}

	printf("\n LOWEST selection probability:\t %.3f\n", lowest_selection_prob_M(L));
	printf(" HIGHEST selection probability:\t %.3f\n", highest_selection_prob_M(L));
	printf("\n Number of used solutions:\t %i\n", count_non_zero_weights(L));
	printf("\n Time used:\t %.2f seconds\n", L.t);
}

int LotteryDesigner::count_non_zero_weights(lottery L) {
	int counter = 0;
	for (int i = 0; i < L.S.size(); i++) {
		if (L.w[i] > 0.0000001) {
			counter++;
		}
	}

	return counter;
}

double LotteryDesigner::lowest_selection_prob_M(lottery L) {
	double lowest = 1;
	for (int i = 0; i < L.p.size(); i++) {
		if (K->M[i] == 1) {
			if (L.p[i] < lowest) {
				lowest = L.p[i];
			}
		}
	}
	return lowest;
}

double LotteryDesigner::highest_selection_prob_M(lottery L) {
	double highest = 0;
	for (int i = 0; i < L.p.size(); i++) {
		if (K->M[i] == 1) {
			if (L.p[i] > highest) {
				highest = L.p[i];
			}
		}
	}
	return highest;
}

void LotteryDesigner::reset_K(bool print) {
	delete K;
	K = new IPSolver(K_initial, print);
}

std::vector<std::vector<int>> LotteryDesigner::order_agents(std::vector<int> v, int iterations, bool print, unsigned seed) {
	// Define random number generator
	unsigned seed_init = (unsigned)time(NULL);
	std::mt19937 generator(seed_init);
	if (seed != 123456789) {
		// Change the value of the seed if it was submitted.
		generator.seed(seed);
	}
	std::uniform_int_distribution<int> distribution(0, 99999999);

	// Create a vector that will store all orderings
	std::vector<std::vector<int>> orders;

	// if the number of agents in M is smaller than 8, we can simply go over all orderings
	// Otherwise we have to take a sample of all orderings.
	if (v.size() <= 8) {
		iterations = 1;
		for (int i = 0; i < v.size(); i++) {
			iterations = iterations * (i + 1); // Computes the factorial number
		}

		// Find all permutations
		do {
			orders.push_back(v);
		} while (std::next_permutation(v.begin(), v.end()));
	}
	else {
		for (int i = 0; i < iterations; i++) {
			// Randomly order the agents
			// We use Durstenfeld's interpretation of the Fischer-Yates shuffle here (https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle)
			std::vector<int> ordered;
			for (int k = 0; k < v.size(); k++) {
				ordered.push_back(v[k]);
			}

			int next_index, next;

			for (int k = (v.size() - 1); k > 0 - 1; k--) {
				// We will shuffle the agents by in each iteration placing the 'next' agent at the back of the array (worst position)
				std::uniform_int_distribution<int> rand(0, k);
				next_index = rand(generator);
				// These two steps could be made shorter by only calling the random generator once (on interval [0,nstudents])
				// And by then taking the modulus of the number wrt i.
				// Experiment with that if time is really crucial.
				next = ordered[next_index];

				// Swap the positions of 'next' and the currently last element in the array
				ordered[next_index] = ordered[k];
				ordered[k] = next;
			}

			orders.push_back(ordered);
		}
	}

	return orders;
}


//******************************************
//************ MAIN FUNCTIONS **************
//******************************************

lottery LotteryDesigner::uniform(bool print) {
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	// Disable output if not wanted
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   

	// Find all optimal solutions
	// First we exclude all already found solutions to be found again
	GRBLinExpr lin;
	int counter;

	for (int i = 0; i < K->S.size(); i++) {
		K->block_solution(K->S[i]);
	}

	// Enforce the model to find the optimal objective value
	lin = 0;
	for (int j = 0; j < K->I.n; j++) {
		lin += K->I.v[j] * K->X[j];
	}
	for (int j = 0; j < K->I.t; j++) {
		lin += K->I.v[K->I.n + j] * K->Y_var[j];
	}
	K->model->addConstr(lin == K->opt);
	K->model->write("Generated formulations/IPModel.lp");

	clock_t start_time = clock();

	// Now keep solving the model until we can no longer find solutions with the optimal objective value
	bool feasible = true;
	int status;
	while (feasible) {
		K->solve(false);
		status = K->model->get(GRB_IntAttr_Status);
		if (status == 3) {
			// The model is infeasible
			feasible = false;
		}
		else {
			// Make sure the last found solution (which is added at the back of S) can no longer be found
			K->block_solution(K->S.back());
		}

		if (K->S.size() > 1000) {
			// To not keep looking forever, stop the search after 1000 found solutions.
			feasible = false;
		}
	}

	// The final solution assigns an equal weight to each of the optimal solutions
	for (int i = 0; i < K->S.size(); i++) {
		L.S.push_back(K->S[i]);
		L.w.push_back((1 / (double)K->S.size()));
		for (int j = 0; j < K->I.n; j++) {
			L.p[j] += (1 / (double)K->S.size()) * K->S[i].x[j];
		}
	}

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	if (print) {
		print_lottery(L);

		printf("SOLUTIONS:\n");
		for (int i = 0; i < L.S.size(); i++) {
			for (int j = 0; j < L.S[i].x.size(); j++) {
				if (L.S[i].x[j] == true) {
					printf("1 ");
				}
				else {
					printf("0 ");
				}
			}
			printf("\n");
		}
	}

	// Reset 'K', so that we don't "cheat" by making it easier for the other methods.
	reset_K(print);

	return L;
}

lottery LotteryDesigner::leximax(bool print) {
	// We want to build a column generation framework that maximizes the selection probability
	// of the least-selected agent in M
		// First we solve a master problem, based on the optimal solutions we already know.
		// Then we solve our knapack problem with an alternative objective function,
			// and the constraint that the original objective value is equal to the optimal value
		// Add the solution to the master problem again, etc.

	// Disable output if not wanted
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver


	LeximaxMaster* M = new LeximaxMaster(K, print);
	lottery L = M->solve(print);
	delete M;

	if (print) {
		print_lottery(L);
	}

	// Reset 'K', so that we don't "cheat" by making it easier for the other methods.
	reset_K(print);
	
	return L;
}


lottery LotteryDesigner::RSD(int iterations, bool print, unsigned seed) {
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	// Disable output if not wanted
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	// We only want to permute the agents in 'M', because letting the agents in 'Y' or 'N' 
	// choose among the optimal solutions will not have an impact.
	// Create a list of (random) orders of the agents in 'M'
	std::vector<int> M2;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			M2.push_back(i);
		}
	}

	clock_t start_time = clock();

	std::vector<std::vector<int>> orders = order_agents(M2, iterations, print, seed);

	for (int i = 0; i < orders.size(); i++) {
		solution sol = K->RSD_once(orders[i], print);

		// Go through all solutions in S to see if another solution has this ID
			// This method will not work for large numbers
			// Largest number that can be stored in a double is 1.7e308
			// Just add the solution if the ID is too large, don't perform check
		bool found = false;
		if (sol.ID != -1) { // This is the value we give the ID if exponents are too large
			for (int i = 0; i < L.S.size(); i++) {
				if (L.S[i].ID == sol.ID) {
					found = true;
					L.w[i] += 1 / (double)orders.size();
					// Increase the selection weight of this solution
					i = L.S.size();
				}
			}
		}
		// Now add this solution to 'L', with a selection probability equal to 1/orders.size()
		if (found == false) {
			L.S.push_back(sol);
			L.w.push_back(1 / (double)orders.size());
		}


	}

	// Calculate the selection probabilities of the agents, based on the solutions
	for (int i = 0; i < L.S.size(); i++) {
		for (int j = 0; j < K->I.n; j++) {
			L.p[j] += L.w[i] * L.S[i].x[j];
		}
	}

	if (print) {
		print_lottery(L);
	}

	// Reset 'K', so that we don't "cheat" by making it easier for the other methods.
	reset_K(print);

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	return L;
}

lottery LotteryDesigner::RSD_once(bool print, unsigned seed) {
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	// Disable output if not wanted
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	// We only want to permute the agents in 'M', because letting the agents in 'Y' or 'N' 
	// choose among the optimal solutions will not have an impact.
	// Create a list of (random) orders of the agents in 'M'
	std::vector<int> M2;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			M2.push_back(i);
		}
	}

	clock_t start_time = clock();

	std::vector<std::vector<int>> orders = order_agents(M2, 1, print, seed);

	for (int i = 0; i < orders.size(); i++) {
		solution sol = K->RSD_once(orders[i], print);

		// Go through all solutions in S to see if another solution has this ID
			// This method will not work for large numbers
			// Largest number that can be stored in a double is 1.7e308
			// Just add the solution if the ID is too large, don't perform check
		bool found = false;
		if (sol.ID != -1) { // This is the value we give the ID if exponents are too large
			for (int i = 0; i < L.S.size(); i++) {
				if (L.S[i].ID == sol.ID) {
					found = true;
					L.w[i] += 1 / (double)orders.size();
					// Increase the selection weight of this solution
					i = L.S.size();
				}
			}
		}
		// Now add this solution to 'L', with a selection probability equal to 1/orders.size()
		if (found == false) {
			L.S.push_back(sol);
			L.w.push_back(1 / (double)orders.size());
		}


	}

	// Calculate the selection probabilities of the agents, based on the solutions
	for (int i = 0; i < L.S.size(); i++) {
		for (int j = 0; j < K->I.n; j++) {
			L.p[j] += L.w[i] * L.S[i].x[j];
		}
	}

	if (print) {
		print_lottery(L);
	}

	// Reset 'K', so that we don't "cheat" by making it easier for the other methods.
	reset_K(print);

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	return L;
}

lottery LotteryDesigner::rename_variables(int iterations, bool print, unsigned seed) {
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	clock_t start_time = clock();

	// Create a list of (random) orders of the agents
	std::vector<int> A;
	for (int i = 0; i < K->I.n; i++) {
		A.push_back(i);
	}
	std::vector<std::vector<int>> orders = order_agents(A, iterations, print, seed);

	GRBEnv* env = new GRBEnv();
	env->set(GRB_StringParam_LogFile, "Generated Log-files/Rename_indices.log");

	// Now build the knapsack model by adding the variables in these orders, and store the solutions
	for (int j = 0; j < orders.size(); j++) {
		GRBModel* model = new GRBModel(*env);
		//model->set(GRB_IntParam_PoolSolutions, 10); // Limit the number of solutions that will be stored.
		model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver
		//model->getEnv().set(GRB_IntParam_Presolve, 0);

		GRBVar* X = new GRBVar[K->I.n];
		// Build
		for (int i = 0; i < K->I.n; i++) {
			char name_x[13];
			sprintf_s(name_x, "x_%i", i);
			X[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name_x);
		}

		GRBVar* Y_var = new GRBVar[K->I.t];
		// Build
		for (int i = 0; i < K->I.t; i++) {
			char name_y[13];
			sprintf_s(name_y, "y_%i", i);
			if (K->I.Y_bool == 1) {
				Y_var[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name_y);
			}
			else {
				Y_var[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, name_y);
			}
		}

		GRBLinExpr obj;
		obj = 0.0;

		// Objective function
		int agent;
		for (int i = 0; i < K->I.n; i++) {
			agent = orders[j][i];
			obj += K->I.v[agent] * X[i];
		}
		for (int i = 0; i < K->I.t; i++) {
			obj += K->I.v[K->I.n + i] * Y_var[i];
		}
		model->setObjective(obj, GRB_MAXIMIZE);

		// Constraints for each knapsack
		GRBLinExpr con;
		for (int i = 0; i < K->I.m; i++) {
			con = 0;
			for (int k = 0; k < K->I.n; k++) {
				agent = orders[j][k];
				con += K->I.A[i][agent] * X[k];
			}
			for (int k = 0; k < K->I.t; k++) {
				con += K->I.A[i][K->I.n + k] * Y_var[k];
			}
			model->addConstr(con <= K->I.C[i]);
		}

		//model->write("ModelReindexed.lp");
		//model->getEnv().set(GRB_IntParam_Presolve, 0);
		model->optimize();

		// Store solution
		solution sol;
		sol.x = std::vector<bool>(K->I.n, 0);
		for (int i = 0; i < K->I.n; i++) {
			// Now we have to assign this solution value again to the agent that is represented by this variable
			agent = orders[j][i];
			sol.x[agent] = (bool)X[i].get(GRB_DoubleAttr_X);
		}
		sol.y = std::vector<int>(K->I.t, 0);
		for (int i = 0; i < K->I.t; i++) {
			sol.y[i] = Y_var[i].get(GRB_DoubleAttr_X);
		}

		// Compute the binary number represented by the solution.
		int bin = 0;
		int exponent = 0;
		int exponent_sum = 0;
		for (int i = K->I.n - 1; i > -1; i--) {
			if (sol.x[i] == true) {
				bin += pow(2, exponent);
				exponent_sum += exponent;
			}
			exponent++;
		}

		// This method will not work for large numbers
		// Largest number that can be stored in a double is 1.7e308
		// Just add the solution if the ID is too large, don't perform check
		bool found = false;
		if (exponent_sum <= 80) {
			// Go through all solutions in S to see if another solution has this ID
			for (int i = 0; i < L.S.size(); i++) {
				if (L.S[i].ID == bin) {
					found = true;
					L.w[i] += 1 / (double)orders.size();
					// Increase the selection weight of this solution
					i = L.S.size();
				}
			}
		}

		// Now add this solution to 'L', with a selection probability equal to 1/orders.size()
		if (found == false) {
			L.S.push_back(sol);
			L.w.push_back(1 / (double)orders.size());
		}

		delete[] X;
		delete[] Y_var;
		delete model;
	}
	delete env;

	// Calculate the selection probabilities of the agents, based on the solutions
	for (int i = 0; i < L.S.size(); i++) {
		for (int j = 0; j < K->I.n; j++) {
			L.p[j] += L.w[i] * L.S[i].x[j];
		}
	}

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	if (print) {
		print_lottery(L);
	}

	return L;
}

lottery LotteryDesigner::perturb_objective(int iterations, bool print, unsigned seed) {
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	clock_t start_time = clock();

	// NOTE: this only works when the coefficients in the objective function are integer

	// We choose the interval in which the perturbations should lie in order to always still find an optimal solution
	// epsilon \in [-Gamma, Gamma]
		// with Gamma = 1/(|A| * max{v})
	// First, find the maximum weight of the agents in the objective function
	double max_weight = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->I.v[i] > max_weight) {
			max_weight = K->I.v[i];
		}
	}
	double gamma = 1 / ((double)K->I.n * (double)max_weight);

	// Define random number generator
	unsigned seed_init = (unsigned)time(NULL);
	std::mt19937 generator(seed_init);
	if (seed != 123456789) {
		// Change the value of the seed if it was submitted.
		generator.seed(seed);
	}
	std::uniform_real_distribution<double> perturbation(-gamma, gamma);

	GRBLinExpr obj;
	for (int i = 0; i < iterations; i++)
	{
		// Shift the objective function of the problem
		obj = 0.0;

		for (int j = 0; j < K->I.n; j++) {
			obj += K->I.v[j] * (1 + perturbation(generator)) * K->X[j];
		}
		K->model->setObjective(obj, GRB_MAXIMIZE);

		K->model->optimize();

		// Store solution
		solution sol;
		sol.x = std::vector<bool>(K->I.n, 0);
		for (int i = 0; i < K->I.n; i++) {
			sol.x[i] = (bool)K->X[i].get(GRB_DoubleAttr_X);
		}

		// Compute the binary number represented by the solution.
		int bin = 0;
		int exponent = 0;
		for (int i = K->I.n - 1; i > -1; i--) {
			if (sol.x[i] == true) {
				bin += pow(2, exponent);
			}
			exponent++;
		}
		sol.ID = bin;

		// Go through all solutions in S to see if another solution has this ID
		bool found = false;
		for (int i = 0; i < L.S.size(); i++) {
			if (L.S[i].ID == bin) {
				found = true;
				L.w[i] += 1 / (double)iterations;
				// Increase the selection weight of this solution
				i = L.S.size();
			}
		}

		// Now add this solution to 'L', with a selection probability equal to 1/orders.size()
		if (found == false) {
			L.S.push_back(sol);
			L.w.push_back(1 / (double)iterations);
		}
	}

	// Calculate the selection probabilities of the agents, based on the solutions
	for (int i = 0; i < L.S.size(); i++) {
		for (int j = 0; j < K->I.n; j++) {
			L.p[j] += L.w[i] * L.S[i].x[j];
		}
	}

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	if (print) {
		print_lottery(L);
	}

	// Reset 'K', so that we don't "cheat" by making it easier for the other methods.
	reset_K(print);

	return L;
}

lottery LotteryDesigner::Nash(bool print) {
	// Store the result in a lottery object 'L'
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	clock_t start_time = clock();

	// So far, we first find all optimal solutions using 'uniform', and then solve our Nash model for those solutions
	lottery L_uniform = uniform(false);

	try {
		GRBEnv* env = new GRBEnv();
		env->set(GRB_StringParam_LogFile, "Generated Log-files/Nash.log");

		GRBModel* model = new GRBModel(*env);
		//model->set(GRB_IntParam_PoolSolutions, 10); // Limit the number of solutions that will be stored.
		model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

		int size_M = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				size_M++;
			}
		}

		GRBVar* w = new GRBVar[L_uniform.S.size()];
		// The weights of all the optimal solutions
		GRBVar* p = new GRBVar[size_M];
		// The selection probabilities of the agents
		GRBVar* log_p = new GRBVar[size_M];
		// The logarithms of 'p'

		for (int i = 0; i < L_uniform.S.size(); i++) {
			std::ostringstream name_w;
			name_w << "w_" << i;
			w[i] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, name_w.str());
		}
		int counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				char name_p[13];
				sprintf_s(name_p, "p_%i", i);
				p[counter] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, name_p);
				counter++;
			}
		}
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				char name_log_p[13];
				sprintf_s(name_log_p, "log_p_%i", i);
				log_p[counter] = model->addVar(-GRB_INFINITY, 0.0, 1.0, GRB_CONTINUOUS, name_log_p);
				// Note that each of the variables has objective coefficient 1, as we want to maximize the sum of the log_p variables
				counter++;
			}
		}

		GRBConstr* C_p = new GRBConstr[size_M];
		// Constraints that determine 'p' based on the weights 'w'
		GRBGenConstr* C_log_p = new GRBGenConstr[size_M];
		GRBConstr C_sum_1;

		// Declare that p_i is equal to the sum of the weights of the solutions in which i is selected
		GRBLinExpr lin;
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			lin = 0.0;
			if (K->M[i] == 1) {
				for (int s = 0; s < L_uniform.S.size(); s++) {
					if (L_uniform.S[s].x[i] == 1) {
						lin += w[s];
					}
				}
				C_p[counter] = model->addConstr(lin == p[counter]);
				counter++;
			}
		}

		// Declare that the sum of the weights should be equal to one
		lin = 0.0;
		for (int s = 0; s < L_uniform.S.size(); s++) {
			lin += w[s];
		}
		C_sum_1 = model->addConstr(lin == 1);

		// Set objective function
		counter = 0;
		GRBLinExpr obj = 0.0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				obj += log_p[counter];
				counter++;
			}
		}
		model->setObjective(obj, GRB_MAXIMIZE);

		model->write("Generated Formulations/ModelNash.lp");

		// Lastly, let 'log_p[i]' be the log of 'p[i]'
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				C_log_p[counter] = model->addGenConstrLog(p[counter], log_p[counter]);
				counter++;
			}
		}
		model->write("Generated Formulations/ModelNash.lp");

		model->optimize();
		//model->write("ModelNash.lp");

		// Store the result in 'L'
		for (int i = 0; i < L_uniform.S.size(); i++) {
			L.S.push_back(L_uniform.S[i]);
			L.w.push_back(w[i].get(GRB_DoubleAttr_X));
		}

		// Store the selection probabilities of the agents
		for (int i = 0; i < L_uniform.S.size(); i++) {
			for (int j = 0; j < K->I.n; j++) {
				L.p[j] += L.w[i] * L_uniform.S[i].x[j];
			}
		}

		delete[] C_log_p;
		delete[] C_p;
		delete[] log_p;
		delete[] p;
		delete[] w;

		delete model;
		delete env;
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Exception during optimization" << std::endl;

	}

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	if (print) {
		print_lottery(L);
	}

	return L;
}

lottery LotteryDesigner::Nash_LP_relax(bool print) {
	// Store the result in a lottery object 'L'
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	clock_t start_time = clock();

	try {
		GRBEnv* env = new GRBEnv();
		env->set(GRB_StringParam_LogFile, "Generated Log-files/Nash_LP_relax.log");


		GRBModel* model = new GRBModel(*env);
		//model->set(GRB_IntParam_PoolSolutions, 10); // Limit the number of solutions that will be stored.
		//model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver
		//model->set(GRB_DoubleParam_PoolGap, 0.000001); // Limit the search space by setting a gap for the worst possible solution that will be accepted


		int size_M = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				size_M++;
			}
		}

		GRBVar* X = new GRBVar[K->I.n];

		// Build
		for (int i = 0; i < K->I.n; i++) {
			char name_x[13];
			sprintf_s(name_x, "x_%i", i);
			X[i] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, name_x);
		}

		GRBVar* Y_var = new GRBVar[K->I.t];
		for (int i = 0; i < K->I.t; i++) {
			char name_y[13];
			sprintf_s(name_y, "y_%i", i);
			if (K->I.Y_bool == true) {
				Y_var[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name_y);
			}
			else {
				Y_var[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, name_y);
			}
		}

		GRBVar* log_x = new GRBVar[K->I.n];
		// The logarithms of 'x'. These values will only be added to the objective function for agents who are in M

		int counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			char name_log_p[13];
			sprintf_s(name_log_p, "log_x_%i", i);
			log_x[counter] = model->addVar(-6, 0.0, 1.0, GRB_CONTINUOUS, name_log_p);
			// Note that each of the variables has objective coefficient 1, as we want to maximize the sum of the log_p variables
			counter++;
		}

		GRBGenConstr* C_log_x = new GRBGenConstr[size_M];

		GRBLinExpr con;
		for (int i = 0; i < K->I.m; i++) {
			con = 0;
			for (int j = 0; j < K->I.n; j++) {
				con += K->I.A[i][j] * X[j];
			}
			for (int j = 0; j < K->I.t; j++) {
				con += K->I.A[i][K->I.n + j] * Y_var[j];
			}
			model->addConstr(con <= K->I.C[i]);
		}

		// Add constraint to fix the objective value
		GRBLinExpr obj;
		obj = 0.0;

		// Objective function
		for (int i = 0; i < K->I.n; i++) {
			obj += K->I.v[i] * X[i];
		}
		for (int i = 0; i < K->I.t; i++) {
			obj += K->I.v[K->I.n + i] * Y_var[i];
		}
		model->addConstr(obj == K->opt);

		// Lastly, let 'log_x[i]' be the log of 'x[i]' for the agents in M
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				C_log_x[counter] = model->addGenConstrLog(X[i], log_x[i]);
				counter++;
			}
		}

		// Set objective function
		obj = 0.0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				obj += log_x[i];
			}
		}

		model->setObjective(obj, GRB_MAXIMIZE);

		// Lastly, fix the probabilities of the agents in Y and in N
		for (int i = 0; i < K->I.n; i++) {
			if (K->Y[i] == 1) {
				model->addConstr(X[i] == 1);
			}
			else if (K->N[i] == 1) {
				model->addConstr(X[i] == 0);
			}
		}

		model->write("Generated Formulations/ModelNash_LP_relax.lp");

		GRBModel presolved_model = model->presolve();
		presolved_model.write("Generated Formulations/ModelNash_LP_relax_presolved.lp");

		model->optimize();
		//model->write("ModelNash.lp");

		// Store the selection probabilities of the agents
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				L.p[i] = X[i].get(GRB_DoubleAttr_X);
			}
			else if (K->N[i] == 1) {
				L.p[i] = 0;
			}
			else if (K->Y[i] == 1) {
				L.p[i] = 1;
			}
		}

		delete[] C_log_x;
		delete[] log_x;

		delete model;
		delete env;
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Exception during optimization" << std::endl;

	}

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	if (print) {
		print_lottery(L);
	}

	return L;
	
}


std::vector<lottery> LotteryDesigner::compare_methods(std::string s, int iterations, bool interactive_console, unsigned seed) {
	// Solve the named solution concepts and store the lotteries in a vector of lotteries
	std::vector<lottery> L;
	
	// First check if the partitioning of the agents into Y, N, and M is already done
	if (K->done_partition == false) {
		printf("\n\t We first have to find a partitioning of the agents to know which agents are always/never/sometimes selected.\n");
		K->partition(true);
	}
	
	// Not all rules work (yet) for integer Y-variables
	if (K->I.Y_bool == true) {
		// There is only something to compare if the set 'M' is non-empty
		int size_M = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				size_M++;
			}
		}

		if (size_M > 0) {
			if (interactive_console) {
				printf("\n\n ****************************************************************************");
				printf("\n Which methods do you want to compare? Type the corresponding letters + ENTER\n");
				printf("\t U = Uniform\n");
				printf("\t L = Leximax\n");
				printf("\t R = Random Serial Dictatorship\n");
				printf("\t T = RSD (once)\n");
				printf("\t O = Re-index the variables randomly\n");
				printf("\t P = Perturb the objective function randomly\n");
				printf("\t N = Nash solution (maximizes the geometric mean) - HEURISTIC\n");
				printf("\t S = Nash solution (maximizes the geometric mean) - SIMPLICAL DECOMPOSITION\n");
				printf("\t Q = Nash solution (maximizes the geometric mean) - LP-relaxation\n");
				printf("\t C = Nash solution (maximizes the geometric mean) - Column generation (Flanigan et al., 2021)\n");
				printf(" ****************************************************************************\n\n\t ");
				std::cin >> s;				
			}

			int size_M = 0;
			for (int i = 0; i < K->I.n; i++) {
				if (K->M[i] == 1) {
					size_M++;
				}
			}

			std::string letter;

			for (int i = 0; i < s.size(); i++) {
				letter = s[i];
				if (letter == "U") {
					L.push_back(uniform(false));
				}
				else if (letter == "L") {
					L.push_back(leximax(false));
				}
				else if (letter == "R") {
					L.push_back(RSD(iterations, false, seed));
				}
				else if (letter == "T") {
					L.push_back(RSD(iterations, false, seed));
				}
				else if (letter == "O") {
					L.push_back(rename_variables(iterations, false, seed));
				}
				else if (letter == "P") {
					L.push_back(perturb_objective(iterations, false, seed));
				}
				else if (letter == "N") {
					L.push_back(Nash(false));
				}
				else if (letter == "S") {
					SimplicalDecomposition* SD = new SimplicalDecomposition(K, false);
					L.push_back(SD->SD_Nash(true));
					delete SD;
				}
				else if (letter == "Q") {
					L.push_back(Nash_LP_relax(false));
				}
				else if (letter == "C") {
					SimplicalDecomposition* SD = new SimplicalDecomposition(K, false);
					L.push_back(SD->Nash_CG_Flanigan(false));
					delete SD;
				}
			}

			// Now print everything in a nice, fancy table
			printf("\n\n SELECTION PROBABILITIES (of the agents in M):\n\n");
			std::string display;
			printf("----------------");
			for (int i = 0; i < L.size(); i++) {
				printf("--------");
			}
			printf("\n\t");
			for (int i = 0; i < s.size(); i++) {
				letter = s[i];
				if (letter == "U") {
					display = "UNIF*";
				}
				else if (letter == "R") {
					display = " RSD";
				}
				else if (letter == "L") {
					display = "LEXI";
				}
				else if (letter == "O") {
					display = "INDEX*";
				}
				else if (letter == "P") {
					display = "PERT";
				}
				else if (letter == "N") {
					display = "NASH_h*";
				}
				else if (letter == "S") {
					display = "NASH_sd";
				}
				else if (letter == "Q") {
					display = "NASHLP";
				}
				else if (letter == "C") {
					display = "NASHCG";
				}
				printf("\t %s", display.c_str());
			}
			printf("\n----------------");
			for (int i = 0; i < L.size(); i++) {
				printf("--------");
			}
			printf("\n");

			for (int i = 0; i < K->I.n; i++) {
				if (K->M[i] == 1) {
					if (i == 0) {
						printf(" Agent  0:");
					}
					else {
						printf(" Agent %2.i:", i);
					}
					for (int j = 0; j < L.size(); j++) {
						printf("\t %.2f", L[j].p[i]);
					}
					printf("\n");
				}
			}
			printf("----------------");
			for (int i = 0; i < L.size(); i++) {
				printf("--------");
			}

			printf("\n Minimum :");
			double minimum;
			for (int i = 0; i < L.size(); i++) {
				minimum = 1;
				for (int j = 0; j < K->I.n; j++) {
					if (K->M[j] == 1) {
						if (L[i].p[j] < minimum) {
							minimum = L[i].p[j];
						}
					}
				}
				printf("\t %.2f", minimum);
			}

			printf("\n\n Mean :\t");
			double average;
			for (int i = 0; i < L.size(); i++) {
				average = 0;
				for (int j = 0; j < K->I.n; j++) {
					if (K->M[j] == 1) {
						average += L[i].p[j];
					}
				}
				printf("\t %.2f", average / (double)size_M);
			}

			printf("\n Geom. mean :");
			double g_average;
			for (int i = 0; i < L.size(); i++) {
				g_average = 1;
				for (int j = 0; j < K->I.n; j++) {
					if (K->M[j] == 1) {
						g_average = g_average * L[i].p[j];
					}
				}
				printf("\t %.2f", pow(g_average, (double)1 / (double)size_M));
			}

			printf("\n\n Solut. used:");
			int sum;
			for (int i = 0; i < L.size(); i++) {
				sum = 0;
				for (int j = 0; j < L[i].S.size(); j++) {
					if (L[i].w[j] > 0) {
						sum++;
					}
				}
				printf("\t %i", sum);
			}

			printf("\n Time :\t");
			for (int i = 0; i < L.size(); i++) {
				printf("\t %.2f", L[i].t);
			}
			printf("\n----------------");
			for (int i = 0; i < L.size(); i++) {
				printf("--------");
			}

			// DISCLAIMERS:
			printf("\n");
			if (K->I.Y_coeff_zero == false) {
				printf("\n\n *DISCLAIMER: not all objective coefficients of the Y-variables are equal to zero.");
				printf("\n\t This means that UNIFORM and RSD might not output a Pareto-optimal lottery over the solutions.");
			}
			for (int i = 0; i < s.size(); i++) {
				letter = s[i];
				if (letter == "U") {
					printf("\n\n *DISCLAIMER 'UNIFORM': current maximum number of optimal solutions considered is set to 1000.");
				}
				else if (letter == "O") {
					printf("\n\n *DISCLAIMER 'INDEX': presolve is set to off.");
				}
				else if (letter == "N") {
					printf("\n\n *DISCLAIMER 'NASH': current maximum number of optimal solutions considered is set to 1000.");
				}
			}

			printf("\n");
			for (int i = 0; i < K->I.n; i++) {
				for (int s = 0; s < L[0].S.size(); s++) {
					printf("\t%i", (int)L[0].S[s].x[i]);
				}
				printf("\n");
			}

			printf("\n\n");
		}


		else {
			printf("\n\n\t All agents are either always or never selected in this instance.\n");

			// Return an empty lottery;
			lottery L_empty;
			L_empty.p = std::vector<double>(K->I.n, 0);
			L_empty.w = std::vector<double>();
			L_empty.S = std::vector<solution>();
			L.push_back(L_empty);
		}
	}
	else {
		printf("\n\n\t Not all functions work (yet) for integer Y-variables.\n");
		printf("\t (More specifically, think about how to exclude a solution in 'block_solution' for integer Ys.)\n");
	}

	return L;
}



//******************************************
//****** CONSTRUCTOR AND DESTRUCTOR ********
//******************************************

LotteryDesigner::LotteryDesigner(IPSolver* K_in, bool print) {
	K = new IPSolver(K_in, print);
	K_initial = new IPSolver(K_in, print);

}

LotteryDesigner::~LotteryDesigner() {
	delete K;
	delete K_initial;
}
