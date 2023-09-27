#include "IPSolver.h"

double IPSolver::solve(bool print) {

	solver_times++;

	//model->write("Generated Formulations/IPModel.lp");

	model->optimize();
	int status = model->get(GRB_IntAttr_Status);
	double z;
	if (status != 3) { // If feasible
		z = model->get(GRB_DoubleAttr_ObjVal);

		// Store solution
		std::vector<double> x(I.n, 0);
		for (int j = 0; j < I.n; j++) {
			x[j] = (bool)X[j].get(GRB_DoubleAttr_X);
		}

		std::vector<double> y(I.t, 0);
		for (int j = 0; j < I.t; j++) {
			y[j] = Y_var[j].get(GRB_DoubleAttr_X);
		}

		if (print) {
			for (int j = 0; j < I.n; j++) {
				if (x[j] == 1) {
					printf("\tX_%i = %.2f\n", j, (double)x[j]);
				}
			}
			for (int j = 0; j < I.t; j++) {
				if (y[j] > 0) {
					printf("\tY_%i = %.2f\n", j, (double)y[j]);
				}
			}
		}

		// Check: capacities satisfied?
		std::vector<double> sum(I.m, 0);
		for (int k = 0; k < I.m; k++) {
			for (int j = 0; j < I.n; j++) {
				if (x[j] > 0) {
					sum[k] += I.A[k][j];
				}
			}
		}

		if (solver_times == 1) {
			// Initialize 'opt' if solved for the first time
			opt = z;
		}

		// Fill in Y, M and N
		double obj_val = 0;	// We should calculate it ourselves, because when we use this problem as a pricing problem
							// the obtained objective value is different
		for (int i = 0; i < I.n; i++) {
			obj_val += I.v[i] * x[i];
		}
		for (int i = 0; i < I.t; i++) {
			obj_val += I.v[I.n + i] * y[i];
		}
		if (obj_val == opt) {
			for (int i = 0; i < I.n; i++) {
				if (Y[i] == -1) {
					if (x[i] == 0) {
						// i is definitely not selected in all solutions
						Y[i] = 0;
					}
				}
				if (N[i] == -1) {
					if (x[i] == 1) {
						// i is selected in at least one solution
						N[i] = 0;
					}
				}
			}

			// Add the solution to S if it is not yet in there
			solution s;

			// Compute the binary number represented by the solution.
			int bin = 0;
			int exponent = 0;
			int exponent_sum = 0;
			for (int i = I.n - 1; i > -1; i--) {
				if (x[i] == true) {
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
				for (int i = 0; i < S.size(); i++) {
					if (S[i].ID == bin) {
						found = true;
						i = S.size();
					}
				}
			}

			// Add solution to S
			s.x = x;
			s.y = y;
			s.ID = bin;
			if (found == false) {
				S.push_back(s);
			}
		}
	}
	else {
		z = -11223344;
	}

	done_solve = true;

	return z;
}

double IPSolver::solve_partition(bool print) {

	solver_times++;

	//model->write("Generated Formulations/IPModel.lp");

	// Create callback and add it to the model.
	PartitionCallback cb = PartitionCallback(opt, print);
		// 'opt' is the optimal solution of the problem
	model->setCallback(&cb);

	model->optimize();
	int status = model->get(GRB_IntAttr_Status);
	double z;
	if (status == 2) { // If optimal
		z = model->get(GRB_DoubleAttr_ObjVal);

		// Store solution
		std::vector<double> x(I.n, 0);
		for (int j = 0; j < I.n; j++) {
			x[j] = (bool)X[j].get(GRB_DoubleAttr_X);
		}

		std::vector<double> y(I.t, 0);
		for (int j = 0; j < I.t; j++) {
			y[j] = Y_var[j].get(GRB_DoubleAttr_X);
		}

		if (print) {
			for (int j = 0; j < I.n; j++) {
				if (x[j] == 1) {
					printf("\tX_%i = %.2f\n", j, (double)x[j]);
				}
			}
			for (int j = 0; j < I.t; j++) {
				if (y[j] > 0) {
					printf("\tY_%i = %.2f\n", j, (double)y[j]);
				}
			}
		}

		// Check: capacities satisfied?
		/*std::vector<double> sum(I.m, 0);
		for (int k = 0; k < I.m; k++) {
			for (int j = 0; j < I.n; j++) {
				if (x[j] > 0) {
					sum[k] += I.A[k][j];
				}
			}
		}*/

		if (solver_times == 1) {
			// Initialize 'opt' if solved for the first time
			opt = z;
		}

		// Fill in Y, M and N
		double obj_val = 0;	// We should calculate it ourselves, because when we use this problem as a pricing problem
							// the obtained objective value is different
		for (int i = 0; i < I.n; i++) {
			obj_val += I.v[i] * x[i];
		}
		for (int i = 0; i < I.t; i++) {
			obj_val += I.v[I.n + i] * y[i];
		}
		if (obj_val == opt) {
			for (int i = 0; i < I.n; i++) {
				if (Y[i] == -1) {
					if (x[i] == 0) {
						// i is definitely not selected in all solutions
						Y[i] = 0;
					}
				}
				if (N[i] == -1) {
					if (x[i] == 1) {
						// i is selected in at least one solution
						N[i] = 0;
					}
				}
			}

			// Add the solution to S if it is not yet in there
			solution s;

			// Compute the binary number represented by the solution.
			int bin = 0;
			int exponent = 0;
			int exponent_sum = 0;
			for (int i = I.n - 1; i > -1; i--) {
				if (x[i] == true) {
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
				for (int i = 0; i < S.size(); i++) {
					if (S[i].ID == bin) {
						found = true;
						i = S.size();
					}
				}
			}

			// Add solution to S
			s.x = x;
			s.y = y;
			s.ID = bin;
			if (found == false) {
				S.push_back(s);
			}
		}
	}
	else {
		z = -11223344;
	}

	// We should remove the callback 
	model->setCallback(NULL);

	done_solve = true;

	return z;
}

IP_report IPSolver::solve_return_solution(bool print) {

	IP_report IP_R;
	solution s;
	bool found;

	solver_times++;

	//model->write("Generated Formulations/IPModel.lp");

	model->optimize();
	int status = model->get(GRB_IntAttr_Status);
	double z;
	if (status != 3) { // If feasible
		z = model->get(GRB_DoubleAttr_ObjVal);

		// Store solution
		std::vector<double> x(I.n, 0);
		for (int j = 0; j < I.n; j++) {
			x[j] = (bool)X[j].get(GRB_DoubleAttr_X);
		}

		std::vector<double> y(I.t, 0);
		for (int j = 0; j < I.t; j++) {
			y[j] = Y_var[j].get(GRB_DoubleAttr_X);
		}

		if (print) {
			for (int j = 0; j < I.n; j++) {
				if (x[j] == 1) {
					printf("\tX_%i = %.2f\n", j, (double)x[j]);
				}
			}
			for (int j = 0; j < I.t; j++) {
				if (y[j] > 0) {
					printf("\tY_%i = %.2f\n", j, (double)y[j]);
				}
			}
		}

		// Check: capacities satisfied?
		std::vector<double> sum(I.m, 0);
		for (int k = 0; k < I.m; k++) {
			for (int j = 0; j < I.n; j++) {
				if (x[j] > 0) {
					sum[k] += I.A[k][j];
				}
			}
		}

		if (solver_times == 1) {
			// Initialize 'opt' if solved for the first time
			opt = z;
		}

		// Fill in Y, M and N
		double obj_val = 0;	// We should calculate it ourselves, because when we use this problem as a pricing problem
							// the obtained objective value is different
		for (int i = 0; i < I.n; i++) {
			obj_val += I.v[i] * x[i];
		}
		for (int i = 0; i < I.t; i++) {
			obj_val += I.v[I.n + i] * y[i];
		}
		if (obj_val >= opt - 0.00001) {
			for (int i = 0; i < I.n; i++) {
				if (Y[i] == -1) {
					if (x[i] == 0) {
						// i is definitely not selected in all solutions
						Y[i] = 0;
					}
				}
				if (N[i] == -1) {
					if (x[i] == 1) {
						// i is selected in at least one solution
						N[i] = 0;
					}
				}
			}

			// Add the solution to S if it is not yet in there
			

			// Compute the binary number represented by the solution.
			int bin = 0;
			int exponent = 0;
			int exponent_sum = 0;
			for (int i = I.n - 1; i > -1; i--) {
				if (x[i] == true) {
					bin += pow(2, exponent);
					exponent_sum += exponent;
				}
				exponent++;
			}

			// This method will not work for large numbers
			// Largest number that can be stored in a double is 1.7e308
			// Just add the solution if the ID is too large, don't perform check
			found = false;
			if (exponent_sum <= 80) {
				// Go through all solutions in S to see if another solution has this ID
				for (int i = 0; i < S.size(); i++) {
					if (S[i].ID == bin) {
						found = true;
						i = S.size();
					}
				}
			}

			// Add solution to S
			s.x = x;
			s.y = y;
			s.ID = bin;
			if (found == false) {
				S.push_back(s);
			}
		}
		else {
			found = true;
		}
	}
	else {
		z = -11223344;
		found = true;
	}

	done_solve = true;

	IP_R.s = s;
	IP_R.opt_obj_value = z;
	IP_R.solution_already_in_S = found;
	return IP_R;
}


void IPSolver::analyze(bool print) {
	solve(false);
	printf("\n*******************************************\n\t         PARTITION\n*******************************************\n");
	partition(print);
	//printf("\n\n*******************************************\n\t     IDENTICAL AGENTS\n*******************************************\n");
	//find_identical_agents(print);
}

void IPSolver::partition(bool print) {
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	// First, check whether X-variables are binary (dichotomous) or not (cardinal preferences)
	if (I.X_bool) {

		clock_t start_time = clock();

		for (int i = 0; i < I.n; i++) {
			bool unselected = true;
			// Check for every agent:
			// Whether he belongs to Y:
			if (Y[i] == -1) {
				GRBConstr c = model->addConstr(X[i] == 0);
				int z = this->solve_partition(print);
				if (z != opt) {
					Y[i] = 1;
					M[i] = 0;
					N[i] = 0;
					unselected = false;
					Xmin[i] = 1;
					Xmax[i] = 1;
				}
				model->remove(c);
			}
			if (N[i] == -1) {
				GRBConstr c = model->addConstr(X[i] == 1);
				int z = this->solve_partition(print);
				if (z != opt) {
					N[i] = 1;
					M[i] = 0;
					Y[i] = 0;
					Xmin[i] = 0;
					Xmax[i] = 0;
					unselected = false;
				}
				model->remove(c);
			}

			if (unselected) {
				M[i] = 1;
				Y[i] = 0;
				N[i] = 0;
				Xmin[i] = 0;
				Xmax[i] = 1;
			}

		}

		double length = ((double)(clock() - start_time) / CLK_TCK);

		// Determine the nunmber of agents in M
		M_size = 0;
		for (int i = 0; i < I.n; i++) {
			if (M[i] == 1) {
				M_size++;
			}
		}
		if (print) {
			print_partition(length);
		}

		time_partition = length;

		done_partition = true;
	}


	// CARDINAL PREFERENCES
	else {
		partition_cardinal(print);
	}
}

void IPSolver::greedy_partition(bool print) {
	// We keep track of which agents are not yet covered by an optimal solution
	std::vector<int> covered;

	clock_t start_time = clock();

	for (int i = 0; i < I.n; i++) {
		covered.push_back(i);
	}

	bool terminate = false;
	int constraints_added = 0;

	std::vector<GRBConstr> added_constraints;

	while (terminate == false) {
		//model->write("Generated Formulations/ModelGreedy.lp");
		model->optimize();
		int status = model->get(GRB_IntAttr_Status);
		double z;
		if (status != 3) { // If feasible
			z = model->get(GRB_DoubleAttr_ObjVal);

			// Check if this solution is optimal
			if (z == opt) {
				std::vector<double> x(I.n, 0);
				for (int j = 0; j < I.n; j++) {
					x[j] = (bool)X[j].get(GRB_DoubleAttr_X);
				}

				std::vector<double> y(I.t, 0);
				for (int j = 0; j < I.t; j++) {
					y[j] = Y_var[j].get(GRB_DoubleAttr_X);
				}

				if (print) {
					for (int j = 0; j < I.n; j++) {
						if (x[j] == 1) {
							printf("\tX_%i = %.2f\n", j, (double)x[j]);
						}
					}
					for (int j = 0; j < I.t; j++) {
						if (y[j] > 0) {
							printf("\tY_%i = %.2f\n", j, (double)y[j]);
						}
					}
				}

				// Compute the binary number represented by the solution.
				int bin = 0;
				int exponent = 0;
				int exponent_sum = 0;
				for (int i = I.n - 1; i > -1; i--) {
					if (x[i] == true) {
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
					for (int i = 0; i < S_greedy.size(); i++) {
						if (S_greedy[i].ID == bin) {
							found = true;
							i = S_greedy.size();
						}
					}
				}

				solution s;

				// Add solution to S
				s.x = x;
				s.y = y;
				s.ID = bin;
				if (found == false) {
					S_greedy.push_back(s);
				}

				// Update 'covered'
				for (int i = covered.size() - 1; i > -1; i--) {
					if (x[covered[i]] == 1) { // If agent 'covered[i]' is selected in the found solution...
						covered.erase(covered.begin() + i);
					}
				}

				// Force that the next solution to be found covers an additional agent
				constraints_added++;
				GRBLinExpr lin;
				for (int i = 0; i < covered.size(); i++) {
					lin += X[covered[i]];
				}
				added_constraints.push_back(model->addConstr(lin >= 1));

			}
			else {
				terminate = true;
			}
		}
		else {
			terminate = true;
		}
	}

	double length = ((double)(clock() - start_time) / CLK_TCK);

	// Remove the constraints that were added in this function:
	for (int i = constraints_added-1; i > -1; i--) {
		model->remove(added_constraints[i]);
	}

	time_greedy_partition = length;
	done_greedy_partition = true;
}

void IPSolver::print_partition(double length) {
	printf("\n Y:\n");
	for (int i = 0; i < I.n; i++) {
		if (Y[i] == 1) {
			printf("\t %i\n", i);
		}
	}
	printf("\n M:\n");
	for (int i = 0; i < I.n; i++) {
		if (M[i] == 1) {
			printf("\t %i\n", i);
		}
	}
	printf("\n N:\n");
	for (int i = 0; i < I.n; i++) {
		if (N[i] == 1) {
			printf("\t %i\n", i);
		}
	}

	printf("\n\n Called solver %i times. (%.3f seconds)\n\t(Number of agents = %i)\n", solver_times, length, I.n);
	printf("\n %i agents in M.\n", M_size);

}

void IPSolver::compare_partition_vs_greedy(bool print) {
	if (done_solve == false) {
		solve(false);
	}
	if (done_partition == false) {
		partition(false);
	}
	if (done_greedy_partition == false) {
		greedy_partition(false);
	}

	// Check which agents are in Y,N,M according to the greedy partition;
	std::vector<int> Y_g (I.n, -1);
	std::vector<int> M_g(I.n, -1);
	std::vector<int> N_g(I.n, -1);

	bool in_Y, in_N;
	for (int i = 0; i < I.n; i++) {
		in_Y = true;
		in_N = true;
		for (int s = 0; s < S_greedy.size(); s++) {
			if (S_greedy[s].x[i] == 1) {
				in_N = false;
			}
			else {
				in_Y = false;
			}
			if ((in_N == false) && (in_Y == false)) {
				// Then the agent must be in M
				M_g[i] = 1;
				Y_g[i] = 0;
				N_g[i] = 0;
			}
		}
		if (in_Y == true) {
			Y_g[i] = 1;
			M_g[i] = 0;
			N_g[i] = 0;
		}
		else if (in_N == true) {
			N_g[i] = 1;
			M_g[i] = 0;
			Y_g[i] = 0;
		}
	}

	// Identify which agents are wrongly classified by the greedy approach
	for (int i = 0; i < I.n; i++) {
		// Only the agents who are truly in M can be wrongly classified:
		if (M[i] == 1) {
			if (M_g[i] == 0) {
				greedy_wrongly_classified.push_back(i);
			}
		}
	}

	if (print) {
		printf("\n %i agents are wrongly classified by the greedy partitioning method.\n", (int) greedy_wrongly_classified.size());
		printf("\t (%i agents in M.)\n", M_size);
		printf("\n The following agents belong to M, but were classified by the greedy method as:\n");
		for (int i = 0; i < greedy_wrongly_classified.size(); i++) {
			int agent = greedy_wrongly_classified[i];
			printf("\t Agent %i: ", agent);
			if (Y_g[agent] == 1) {
				printf("Y\n");
			}
			else if (N_g[agent] == 1) {
				printf("N\n");
			}
		}
	}

}


void IPSolver::find_identical_agents(bool print) {

	// The type of algorithm depends on the type of Y-variables. 
	// For now, we only have an algorithm for when there are no Y-variables.
	// 
	if (I.t == 0) {
		// Use the algorithm ase described in the paper. 
			// Identify all pairs of agents with identical objective function coefficients
			// And for each of them, solve at most two integer programs.

		clock_t start_time = clock();

		// Start by creating a matrix 'equal_obj_coeff'
		std::vector<std::vector<int>> equal_obj_coeff(I.n, std::vector<int>(I.n, 0));

		// Starting from the first agents, fill in all agents with identical objective coefficients.
		for (int i = 0; i < I.n; i++) {
			for (int j = i + 1; j < I.n; j++) {
				if (I.v[i] == I.v[j]) {
					equal_obj_coeff[i][j] = 1;
				}
			}
		}

		// (GUROBI environment already created in constructor)

		int solver_count = 0;
		// Go through all equal pairs
		// We only care about the agents in M!
		for (int i = 0; i < I.n; i++) {
			if (M[i] == 1) {
				for (int j = i + 1; j < I.n; j++) {
					if (M[j] == 1) {
						if (equal_obj_coeff[i][j] == 1) {
							if (identical[i][j] == -1) { // If we didn't fill in the matrix by inference (see last step this for-loop)
								identical[i][j] = 0; // Will be overwritten if needed.

								// Check if there are optimal solutions in which only one of the two agents is selected

								bool check_i_j, check_j_i;
								check_i_j = find_identical_agents_aid(i, j, print);
								solver_count++;

								if (check_i_j == false) {
									check_j_i = find_identical_agents_aid(j, i, print);
									solver_count++;
									if (check_j_i == false) {
										// If both sets of constraints are infeasible, the agents are identical.
										this->identical[i][j] = 1;
									}
								}
							}
						}
					}
				}

				// We can infer some other identical agents: 
					// Imagine we know now that {1,2} is identical and {1,3}, then {2,3} is also identical.
					// Similarly, if we know that {1,2} is identical, and {1,3} is not, then {2,3} is also not identical
				for (int j = i + 1; j < I.n; j++) {
					if (M[j] == 1) {
						for (int k = j + 1; k < I.n; k++) {
							if (M[k] == 1) {
								if (identical[i][j] == 1) {
									if (identical[i][k] == 1) {
										identical[j][k] = 1;
									}
									else {
										identical[j][k] = 0;
									}
								}
							}
						}
					}
				}
			}
		}


		double length = ((double)(clock() - start_time) / CLK_TCK);

		int sum = 0;
		for (int i = 0; i < I.n; i++) {
			for (int j = i + 1; j < I.n; j++) {
				if (identical[i][j] == 1) {
					sum++;
				}
			}
		}

		// Divide the agents in clusters
		identical_cluster(print);

		if (print) {
			printf("\tCalled solver %i times (%.3f seconds).\n\n", solver_count, length);
		}
	}
	else {
		printf("\n We don't have an algorithm (yet) to find identical agents when there are auxiliary Y-variables.\n\n");
	}

	done_find_identical = true;
}

bool IPSolver::find_identical_agents_aid(int i, int j, bool print) {
	// Help function, solves the integer program in which x_i=1, x_j=0 is an optimal solution
	// and returns true if x_j=1,x_i = 0 is also an optimal solution.

	double max_number = 0;
	for (int l = 0; l < I.C.size(); l++) {
		for (int k = 0; k < I.n; k++) {
			if (std::abs(I.A[l][k]) > max_number) {
				max_number = std::abs(I.A[l][k]);
			}
		}
	}

	// Environment used for finding identical agents
	GRBModel* model_identical = new GRBModel(*env_identical);
	GRBVar* X_identical = new GRBVar[I.n];
	for (int k = 0; k < I.n; k++) {
		char name_x[13];
		sprintf_s(name_x, "x_%i", k);
		X_identical[k] = model_identical->addVar(0.0, 1.0, 0.0, GRB_BINARY, name_x);
	}
	GRBVar* S_identical = new GRBVar[I.C.size()];
	int S_counter = 0;
	for (int k = 0; k < I.C.size(); k++) {
		// To reduce the number of constraints, only create s-variable when a_kj > a_ki for constraint k
		if (I.A[k][j] > I.A[k][i]) {
			char name_x[13];
			sprintf_s(name_x, "s_%i", S_counter);
			S_identical[S_counter] = model_identical->addVar(0.0, 1.0, 0.0, GRB_BINARY, name_x);
			S_counter++;
		}
	}

	model_identical->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	try {
		for (int l = 0; l < I.C.size(); l++) {
			GRBLinExpr lin = 0;
			for (int k = 0; k < I.n; k++) {
				if (k != i) {
					if (k != j) {
						lin += I.A[l][k] * X_identical[k];
					}
				}
			}
			model_identical->addConstr(lin + I.A[l][i] <= I.C[l]);
		}

		//model_identical->write("Identical.lp");

		// Objective function original solution
		GRBLinExpr lin = 0;
		for (int k = 0; k < I.n; k++) {
			if (k != i) {
				if (k != j) {
					lin += I.v[k] * X_identical[k];
				}
			}
		}
		model_identical->addConstr(lin + I.v[i] == opt);

		model_identical->addConstr(X_identical[i] == 0);
		model_identical->addConstr(X_identical[j] == 1);

		double epsilon_l, M_l, diff;
		int counter = 0;
		for (int l = 0; l < I.C.size(); l++) {
			if (I.A[l][j] > I.A[l][i]) { // To reduce the number of constraints (see creation s-variables)
				GRBLinExpr lin = 0;
				for (int k = 0; k < I.n; k++) {
					if (k != i) {
						if (k != j) {
							lin += I.A[l][k] * X_identical[k];
						}
					}
				}
				// Compute 'epsilon_l'
				epsilon_l = 2 * max_number; // very large number
				for (int s = 0; s < I.n_var; s++) {
					if (I.A[l][s] != 0) {
						diff = std::abs(I.A[l][s]);
						if (diff < epsilon_l) {
							epsilon_l = diff;
						}
					}
					for (int t = 0; t < I.n; t++) {
						if (I.A[l][s] != I.A[l][t]) {
							diff = std::abs(I.A[l][s] - I.A[l][t]);
							if (diff < epsilon_l) {
								epsilon_l = diff;
							}
						}
					}
				}

				// Compute 'M_l'
				if (I.Y_bool == true) {
					M_l = I.C[l] + epsilon_l;
					for (int s = 0; s < I.n_var; s++) {
						if (s != i) {
							if (I.A[l][s] < 0) {
								M_l += -I.A[l][s];
							}
						}
					}
				}
				else {// If Y-variables can be integer, it is no longer straightforward to find an upper bound...
					// We assume a highest value for Y-variables of 'max_value', but we check later if this is true.
					M_l = I.C[l] + epsilon_l;
					for (int s = 0; s < I.n_var; s++) {
						if (s != i) {
							if (I.A[l][s] < 0) {
								M_l += -max_Y_value * I.A[l][s];
							}
						}
					}
				}

				lin += M_l * S_identical[counter];
				counter++;
				model_identical->addConstr(lin + I.A[l][j] >= I.C[l] + epsilon_l);
			}
		}

		lin = 0;
		for (int l = 0; l < S_counter; l++) {
			lin += 1 * S_identical[l];
		}
		model_identical->addConstr(lin <= S_counter - 1);
	}
	catch (GRBException e) {
		std::cout << "Error number: " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}

	//model_identical->write("Generated Formulations/ModelIdentical.lp");

	model_identical->optimize();

	int status = model_identical->get(GRB_IntAttr_Status);

	bool result;
	if (status == 3) { // Model is infeasible
		result = false;
	}
	else {
		result = true;
		if (false) {
			// Print solution
			std::vector<int> solution;
			for (int k = 0; k < I.m; k++) {
				solution.push_back(X_identical[k].get(GRB_DoubleAttr_X));
			}
			for (int k = 0; k < I.m; k++) {
				if (solution[k] == 1) {
					printf("\tX[%i] = 1\n", k);
				}
			}
		}
	}

	// Delete the model
	delete[] X_identical;
	delete[] S_identical;
	delete model_identical;

	return result;
}

void IPSolver::identical_cluster(bool print) {
	// Create a vector of the agents that are not yet in a cluster
	std::vector<int> not_clustered;
	for (int k = 0; k < I.n; k++) {
		not_clustered.push_back(k);
	}
	cluster_count = 0;

	// First remove all the agents that are not in M (so that are in Y or N)
	for (int k = not_clustered.size() - 1; k > -1; k--) {
		if (M[k] != 1) {
			not_clustered.erase(not_clustered.begin() + k);
		}
	}

	while (not_clustered.size() > 0) {
		int agent = not_clustered[0];
		bool in_cluster = false;

		// Go from end to front
		for (int k = not_clustered.size() - 1; k > -1; k--) {
			if (identical[agent][not_clustered[k]] == 1) {
				// Add agents to cluster
				identical_clustered[agent] = cluster_count;
				identical_clustered[not_clustered[k]] = cluster_count;
				in_cluster = true;

				// Remove agent 'not_clustered[k]' from cluster
				not_clustered.erase(not_clustered.begin() + k);

			}
		}
		if (in_cluster == true) {
			cluster_count++;
		}

		// Remove agent 'not_clustered[0]' from the vector
		not_clustered.erase(not_clustered.begin());
	}

	if (print) {
		printf("\nThere are %i clusters of identical agents in M:\n", cluster_count);
		for (int i = 0; i < cluster_count; i++) {
			printf("\tCLUSTER %i:\n", i + 1);
			for (int k = 0; k < I.n; k++) {
				if (identical_clustered[k] == i) {
					printf("\t  %i\n", k);
				}
			}
		}
	}
}

solution IPSolver::RSD_once(std::vector<int> order, bool print) {
	solution sol;
	RSD_fixed = std::vector<GRBConstr>();
	
	// We only want to permute the agents in 'M', because letting the agents in 'Y' or 'N' 
	// choose among the optimal solutions will not have an impact.
	if (M_size > 0) {
		int counter = 0;
		GRBLinExpr obj = 0.0;
		double delta = 1;

		// We go through 'order' and remove all the agents that are not in M, so that are in Y or N
		for (int i = order.size() - 1; i > -1; i--) {
			if (M[order[i]] == 0) {
				order.erase(order.begin() + i);
			}
		}

		// Store solution of the agents of which we perturbed the objective coefficients
		sol.x = std::vector<double>(I.n, 0);
		sol.y = std::vector<double>(I.t, 0);

		// For the agents in 'Y' and 'N' we already know the solution values.
		for (int i = 0; i < I.n; i++) {
			if (Y[i] == 1) {
				sol.x[i] = (bool) 1;
			}
			else if (N[i] == 1) {
				sol.x[i] = (bool) 0;
			}
		}


		// Taking into account the precision of the solver, we will perturb in different steps
		// For a precision of 1e-6, we can perturb \floor(-log_2(1e-6)) = 19 agents at once
			// This parameter is defined in the function 'InitializeVariables'
		double z;
		for (int j = 0; j < M_size; j=j+number_of_agents_at_once_RSD) {

			obj = 0.0;

			// Build the objective function 
			for (int i = 0; i < I.n; i++) {
				obj += I.v[i] * X[i];
			}

			// Now we go through 'order' from the beginning, and shift the objective function accordingly.
			for (int i = j; i < j+number_of_agents_at_once_RSD; i++) {
				if (i < order.size()) {
					counter++;
					delta = delta * 0.5;
					obj += (delta)*X[order[i]];
				}
				else {
					i = j + number_of_agents_at_once_RSD;
				}
			}

			// Lastly, add the objective coefficients of the Y-variables
			for (int i = 0; i < I.t; i++) {
				obj += I.v[I.n + i] * Y_var[i];
			}

			model->setObjective(obj, GRB_MAXIMIZE);
			//model->write("Generated Formulations/ModelRSD.lp");

			// Optimize
			model->optimize();
			int status = model->get(GRB_IntAttr_Status);
			GRBLinExpr expr;
			if (status != 3) { // If feasible
				for (int i = j; i < j + number_of_agents_at_once_RSD; i++) {
					if (i < order.size()) {
						sol.x[order[i]] = (bool)X[order[i]].get(GRB_DoubleAttr_X);
					}
					else {
						i = j + number_of_agents_at_once_RSD;
					}
				}
				//for (int j = 0; j < I.n; j++) {
					/*if (Y[j] == 1) { // Checking this saves us a bit of time on the expensive .get function
						sol.x[j] = (bool)1;
					}
					else if (N[j] == 1) {
						sol.x[j] = (bool)0;
					}
					else {
						sol.x[j] = (bool)X[j].get(GRB_DoubleAttr_X);
					}*/

					//sol.x[j] = (bool)X[j].get(GRB_DoubleAttr_X);

				//}
				for (int i = 0; i < I.t; i++) {
					sol.y[i] = Y_var[i].get(GRB_DoubleAttr_X);
				}

				z = model->get(GRB_DoubleAttr_ObjVal);

				// Fix the variables by adding constraints:
				for (int i = j; i < j + number_of_agents_at_once_RSD; i++) {
					if (i < order.size()) {
						expr = 0.0;
						expr = X[order[i]];
						//model->write("Generated Formulations/ModelRSD.lp");
						RSD_fixed.push_back(model->addConstr(expr == (int)sol.x[order[i]]));
						//model->write("Generated Formulations/ModelRSD.lp");
					}
					else {
						i = j + number_of_agents_at_once_RSD;
					}
				}
			}
		}


		// Compute the binary number represented by the solution.
		int bin = 0;
		int exponent = 0;
		int exponent_sum = 0;
		for (int i = I.n - 1; i > -1; i--) {
			if (sol.x[i] == true) {
				bin += pow(2, exponent);
				exponent_sum += exponent;
			}
			exponent++;
		}
		if (exponent_sum <= 80) {
			sol.ID = bin;
		}
		else {
			sol.ID = -1;
		}

		// Delete the constraints that were added to fix the variables
		for (int i = M_size - 1; i > -1; i--) {
			//model->write("Generated Formulations/ModelRSD.lp");
			model->remove(RSD_fixed[i]);
			//model->write("Generated Formulations/ModelRSD.lp");

		}
		RSD_fixed.clear();

	}
	else {
		sol.x = std::vector<double>(I.n, 0);
		sol.y = std::vector<double>(I.t, 0);

		// We know for each agent whether they are selected or not
		for (int i = 0; i < I.n; i++) {
			if (Y[i] == 1) {
				sol.x[i] = 1;
				model->addConstr(X[i] == 1);
			}
			else {
				//sol.x[i] = 0; // Already true because of the initialization of sol.x
				model->addConstr(X[i] == 0);
			}
		}

		// Now we solve the model to find the correct values of the Y-variables
		
		// FIND THE CORRECT VALUES OF THE Y-VARIABLES STILL!
		model->optimize();
		int status = model->get(GRB_IntAttr_Status);
		if (status != 3) { // If feasible
			// Store solution
			sol.y = std::vector<double>(I.t, 0);
			for (int j = 0; j < I.t; j++) {
				sol.y[j] = Y_var[j].get(GRB_DoubleAttr_X);
			}
		}

		// Compute the binary number represented by the solution.
		int bin = 0;
		int exponent = 0;
		int exponent_sum = 0;
		for (int i = I.n - 1; i > -1; i--) {
			if (sol.x[i] == true) {
				bin += pow(2, exponent);
				exponent_sum += exponent;
			}
			exponent++;
		}
		if (exponent_sum <= 80) {
			sol.ID = bin;
		}
		else {
			sol.ID = -1;
		}
	}

	return sol;
}

solution IPSolver::RSD_once_no_partition(std::vector<int> order, bool print) {
	// This is the updated function
	solution sol;
	RSD_fixed = std::vector<GRBConstr>();

	// We want to permute all agents, not only the agents in 'M'.

	int counter = 0;
	GRBLinExpr obj = 0.0;
	double delta = 1;

	// Store solution of the agents of which we perturbed the objective coefficients
	sol.x = std::vector<double>(I.n, 0);
	sol.y = std::vector<double>(I.t, 0);

	// Taking into account the precision of the solver, we will perturb in different steps
	// For a precision of 1e-6, we can perturb \floor(-log_2(1e-6)) = 19 agents at once
		// This parameter is defined in the function 'InitializeVariables'
	double z;
	for (int j = 0; j < I.n; j = j + number_of_agents_at_once_RSD) {

		obj = 0.0;

		// Build the objective function 
		for (int i = 0; i < I.n; i++) {
			obj += I.v[i] * X[i];
		}

		// Now we go through 'order' from the beginning, and shift the objective function accordingly.
		for (int i = j; i < j + number_of_agents_at_once_RSD; i++) {
			if (i < order.size()) {
				counter++;
				delta = delta * 0.5;
				obj += (delta)*X[order[i]];
			}
			else {
				i = j + number_of_agents_at_once_RSD;
			}
		}

		// Lastly, add the objective coefficients of the Y-variables
		for (int i = 0; i < I.t; i++) {
			obj += I.v[I.n + i] * Y_var[i];
		}

		model->setObjective(obj, GRB_MAXIMIZE);
		//model->write("Generated Formulations/ModelRSD.lp");

		// Optimize
		model->optimize();
		int status = model->get(GRB_IntAttr_Status);
		GRBLinExpr expr;
		if (status != 3) { // If feasible
			for (int i = j; i < j + number_of_agents_at_once_RSD; i++) {
				if (i < order.size()) {
					sol.x[order[i]] = (bool)X[order[i]].get(GRB_DoubleAttr_X);
				}
				else {
					i = j + number_of_agents_at_once_RSD;
				}
			}
			//for (int j = 0; j < I.n; j++) {
				/*if (Y[j] == 1) { // Checking this saves us a bit of time on the expensive .get function
					sol.x[j] = (bool)1;
				}
				else if (N[j] == 1) {
					sol.x[j] = (bool)0;
				}
				else {
					sol.x[j] = (bool)X[j].get(GRB_DoubleAttr_X);
				}*/

				//sol.x[j] = (bool)X[j].get(GRB_DoubleAttr_X);

			//}
			for (int i = 0; i < I.t; i++) {
				sol.y[i] = Y_var[i].get(GRB_DoubleAttr_X);
			}

			z = model->get(GRB_DoubleAttr_ObjVal);

			// Fix the variables by adding constraints:
			for (int i = j; i < j + number_of_agents_at_once_RSD; i++) {
				if (i < order.size()) {
					expr = 0.0;
					expr = X[order[i]];
					//model->write("Generated Formulations/ModelRSD.lp");
					RSD_fixed.push_back(model->addConstr(expr == (int)sol.x[order[i]]));
					//model->write("Generated Formulations/ModelRSD.lp");
				}
				else {
					i = j + number_of_agents_at_once_RSD;
				}
			}
		}
	}

	// Compute the binary number represented by the solution.
	int bin = 0;
	int exponent = 0;
	int exponent_sum = 0;
	for (int i = I.n - 1; i > -1; i--) {
		if (sol.x[i] == true) {
			bin += pow(2, exponent);
			exponent_sum += exponent;
		}
		exponent++;
	}
	if (exponent_sum <= 80) {
		sol.ID = bin;
	}
	else {
		sol.ID = -1;
	}

	// Delete the constraints that were added to fix the variables
	for (int i = I.n - 1; i > -1; i--) {
		//model->write("Generated Formulations/ModelRSD.lp");
		model->remove(RSD_fixed[i]);
		//model->write("Generated Formulations/ModelRSD.lp");

	}
	RSD_fixed.clear();

	return sol;
}


solution IPSolver::RSD_heuristic_no_partition(std::vector<int> order, bool print) {
	// This is the updated function
	solution sol;
	RSD_fixed = std::vector<GRBConstr>();

	// We want to permute all agents, not only the agents in 'M'.

	int counter = 0;
	GRBLinExpr obj = 0.0;
	double delta = 1;

	// Store solution of the agents of which we perturbed the objective coefficients
	sol.x = std::vector<double>(I.n, 0);
	sol.y = std::vector<double>(I.t, 0);

	// Taking into account the precision of the solver, we will perturb ONCE
	// For a precision of 1e-6, we can perturb \floor(-log_2(1e-6)) = 19 agents at once
		// This parameter is defined in the function 'InitializeVariables'
	double z;
	int j = 0;

	obj = 0.0;

	// Build the objective function 
	for (int i = 0; i < I.n; i++) {
		obj += I.v[i] * X[i];
	}

	// Now we go through 'order' from the beginning, and shift the objective function accordingly.
	for (int i = j; i < j + number_of_agents_at_once_RSD; i++) {
		if (i < order.size()) {
			counter++;
			delta = delta * 0.5;
			obj += (delta)*X[order[i]];
		}
		else {
			i = j + number_of_agents_at_once_RSD;
		}
	}

	// Lastly, add the objective coefficients of the Y-variables
	for (int i = 0; i < I.t; i++) {
		obj += I.v[I.n + i] * Y_var[i];
	}

	model->setObjective(obj, GRB_MAXIMIZE);
	//model->write("Generated Formulations/ModelRSD.lp");

	// Optimize
	model->optimize();
	int status = model->get(GRB_IntAttr_Status);
	GRBLinExpr expr;
	if (status != 3) { // If feasible
		for (int i = 0; i < I.n; i++) {
			sol.x[i] = (bool)X[i].get(GRB_DoubleAttr_X);
		}
		//for (int j = 0; j < I.n; j++) {
			/*if (Y[j] == 1) { // Checking this saves us a bit of time on the expensive .get function
				sol.x[j] = (bool)1;
			}
			else if (N[j] == 1) {
				sol.x[j] = (bool)0;
			}
			else {
				sol.x[j] = (bool)X[j].get(GRB_DoubleAttr_X);
			}*/

			//sol.x[j] = (bool)X[j].get(GRB_DoubleAttr_X);

		//}
		for (int i = 0; i < I.t; i++) {
			sol.y[i] = Y_var[i].get(GRB_DoubleAttr_X);
		}

		z = model->get(GRB_DoubleAttr_ObjVal);
	}

	// Compute the binary number represented by the solution.
	int bin = 0;
	int exponent = 0;
	int exponent_sum = 0;
	for (int i = I.n - 1; i > -1; i--) {
		if (sol.x[i] == true) {
			bin += pow(2, exponent);
			exponent_sum += exponent;
		}
		exponent++;
	}
	if (exponent_sum <= 80) {
		sol.ID = bin;
	}
	else {
		sol.ID = -1;
	}

	return sol;
}


void IPSolver::block_solution(solution sol) {
	if (I.X_bool == true) {
		GRBLinExpr expr = 0.0;
		int counter = 0;
		for (int i = 0; i < I.n; i++) {
			if (sol.x[i] == 1) {
				expr += X[i];
				counter++;
			}
		}
		for (int i = 0; i < I.t; i++) {
			if (sol.y[i] == 1) {
				expr += Y_var[i];
				counter++;
			}
		}
		model->addConstr(expr <= (counter - 1));

		//model->write("Generated Formulations/IPModel.lp");
	}
	else {
		// We are not aware of an elegant method to exclude a single solution when variables can take any integer values...

	}
}

void IPSolver::change_obj_function(std::vector<double> w) {
	GRBLinExpr obj = 0.0;

	for (int i = 0; i < I.n; i++) {
		obj += w[i] * X[i];
	}
	model->setObjective(obj, GRB_MAXIMIZE);
}

std::vector<time_report> IPSolver::compare_time_normal_vs_RSD_variants(int iterations, bool print) {
	std::vector<time_report> TR;
	
	std::vector<int> order;
	for (int i = 0; i < I.n; i++) {
		order.push_back(i); // Just order agents {0,1,2,...,I.n}
	}

	if (done_partition == false) {
		partition(false);
	}

	clock_t start_time = clock();
	solution sol;
	for (int i = 0; i < iterations; i++) {
		sol = RSD_once_no_partition(order, false);
		model->reset();
	}
	double time_RSD_without_partition = ((double)(clock() - start_time) / CLK_TCK);

	start_time = clock();
	for (int i = 0; i < iterations; i++) {
		sol = RSD_once(order, false);
		model->reset();
	}
	double time_RSD = ((double)(clock() - start_time) / CLK_TCK);

	start_time = clock();
	for (int i = 0; i < iterations; i++) {
		solve(false);
		model->reset();
	}
	double time_normal = ((double)(clock() - start_time) / CLK_TCK);

	if (print) {
		printf(" Time normal = %.3f\n", time_normal);
		printf(" Time RSD = %.3f\n", time_RSD);
		printf(" Time RSD (without partition) = %.3f\n", time_RSD_without_partition);
		
	}

	time_report empty;
	empty.label = "Normal";
	empty.t = time_normal / iterations;
	TR.push_back(empty);
	empty.label = "RSD";
	empty.t = time_RSD / iterations;
	TR.push_back(empty);
	empty.label = "RSD_without_partition";
	empty.t = time_RSD_without_partition / iterations;
	TR.push_back(empty);
	return TR;

}

std::vector<time_report> IPSolver::compare_time_normal_vs_RSD_heuristic(int iterations, bool print) {
	std::vector<time_report> TR;

	std::vector<int> order;
	for (int i = 0; i < I.n; i++) {
		order.push_back(i); // Just order agents {0,1,2,...,I.n}
	}

	if (done_partition == false) {
		partition(false);
	}

	clock_t start_time = clock();
	solution sol;
	for (int i = 0; i < iterations; i++) {
		sol = RSD_heuristic_no_partition(order, false);
		model->reset();
	}
	double time_RSD_heuristic = ((double)(clock() - start_time) / CLK_TCK);

	start_time = clock();
	for (int i = 0; i < iterations; i++) {
		sol = RSD_once(order, false);
		model->reset();
	}
	double time_RSD = ((double)(clock() - start_time) / CLK_TCK);

	start_time = clock();
	for (int i = 0; i < iterations; i++) {
		solve(false);
		model->reset();
	}
	double time_normal = ((double)(clock() - start_time) / CLK_TCK);

	if (print) {
		printf(" Time normal = %.3f\n", time_normal);
		printf(" Time RSD = %.3f\n", time_RSD);
		printf(" Time RSD heuristic = %.3f\n", time_RSD_heuristic);

	}

	time_report empty;
	empty.label = "Normal";
	empty.t = time_normal / iterations;
	TR.push_back(empty);
	empty.label = "RSD";
	empty.t = time_RSD / iterations;
	TR.push_back(empty);
	empty.label = "RSD_heuristic";
	empty.t = time_RSD_heuristic / iterations;
	TR.push_back(empty);
	return TR;

}

std::vector<time_report> IPSolver::compare_time_normal_vs_RSD_without_partition(int iterations, bool print) {
	std::vector<time_report> TR;

	std::vector<int> order;
	for (int i = 0; i < I.n; i++) {
		order.push_back(i); // Just order agents {0,1,2,...,I.n}
	}

	clock_t start_time = clock();
	solution sol;
	for (int i = 0; i < iterations; i++) {
		sol = RSD_once_no_partition(order, false);
		model->reset();
	}
	double time_RSD_without_partition = ((double)(clock() - start_time) / CLK_TCK);
start_time = clock();
	for (int i = 0; i < iterations; i++) {
		solve(false);
		model->reset();
	}
	double time_normal = ((double)(clock() - start_time) / CLK_TCK);

	if (print) {
		printf(" Time normal = %.3f\n", time_normal);
		printf(" Time RSD (without partition) = %.3f\n", time_RSD_without_partition);

	}

	time_report empty;
	empty.label = "Normal";
	empty.t = time_normal / iterations;
	TR.push_back(empty);
	empty.label = "RSD_without_partition";
	empty.t = time_RSD_without_partition / iterations;
	TR.push_back(empty);
	return TR;

}

