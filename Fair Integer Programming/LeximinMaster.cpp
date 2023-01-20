#include "LeximinMaster.h"

lottery LeximinMaster::solve(bool print) {
	//model->getEnv().set(GRB_IntParam_OutputFlag, 0);      //comment to see the output of the solver
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	// Add constraint to knapsack to enforce the objective value to be equal to the optimal objective value
	GRBLinExpr expr = 0.0;

	for (int i = 0; i < K->I.n; i++) {
		expr += K->I.v[i] * K->X[i];
	}
	for (int i = 0; i < K->I.t; i++) {
		expr += K->I.v[K->I.n + i] * K->Y_var[i];
	}
	K->model->addConstr(expr == K->opt);

	int iterations = 1;

	// Create a lottery object corresponding to store the solution
	lottery L;

	// We will use two while loops to find the leximin solution.
	// In the inner loop, we will find the weights of the solutions in order to 
		// maximize the selection probability of the worst-off agent in M "that is remaining"
	// In the outer loop, we decide which agents are remaining.
		// As soon as we have found weights that maximize the selection probability of the worst-off agent 
		// we fix the selection probability of that agent, and remove it's corresponding C constraint.
	std::vector<int> M_remaining = K->M;
	int sum = 0;
	for (int i = 0; i < M_remaining.size(); i++) {
		if (M_remaining[i] == 1) {
			sum++;
		}
	}

	clock_t start_time = clock();

	while (sum > 0) {
		// As long as the solution of the pricing problem has a non-negative reduced cost, continue
		double obj_val_pricing = 1;
		double obj_val_master;
		while (obj_val_pricing > 0) {
			model->write("Generated Formulations/LeximinModel.lp");
			if (print) {
				printf("ITERATION %i\n", iterations);
			}
			model->optimize();
			obj_val_master = model->get(GRB_DoubleAttr_ObjVal);
			getDualValues(print);

			// Change the objective of the knapsack problem, which will serve as our pricing problem
			GRBLinExpr obj = 0.0;

			int counter = 0;
			for (int i = 0; i < K->I.n; i++) {
				if (K->M[i] == 1) {
					obj -= dual_C[counter] * K->X[i];
					counter++;
				}
			}
			obj -= dual_C_Sum1[0];
			K->model->setObjective(obj, GRB_MAXIMIZE);
			obj_val_pricing = K->solve(false);

			// Add a new column to the master if the reduced cost is non-negative
			if (obj_val_pricing > 0) {
				addColumn(K->S.back(), print);
				K->block_solution(K->S.back());
			}

			if (print) {
				// If the model was not infeasible
				if (K->model->get(GRB_IntAttr_Status) != 3) {
					printf("\t Reduced cost = %.5f\n\n", obj_val_pricing);
				}
				else {
					printf("\t Pricing problem INFEASIBLE\n\n");
				}
			}
			iterations++;
		}
		// Store the solution in 'L'
		L.p = std::vector<double>(K->I.n, 0);
		L.w = std::vector<double>();
		L.S = std::vector<solution>();
		for (int i = 0; i < w.size(); i++) { // If we use 'K->S.size()' here, then the lastly added solution is also included, but we didn't solve the model with this variable yet
			L.S.push_back(K->S[i]);
			L.w.push_back(w[i].get(GRB_DoubleAttr_X));
		}

		// Now calculate the selection probabilities of the agents
		for (int i = 0; i < w.size(); i++) {
			for (int j = 0; j < K->I.n; j++) {
				L.p[j] += L.w[i] * K->S[i].x[j];
			}
		}

		// Which agents are selected with the maximin probability?
		int counter = -1;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				counter++;
				// To find the index of the corresponding constraint
			}
			double diff = std::abs(L.p[i] - obj_val_master);
			if (diff <= 0.000001) {
				M_remaining[i] = 0;
				// Remove 'MIN_M' from the corresponding constraint
				model->chgCoeff(C[counter], MIN_M, 0.0);

				// Change RHS from the corresponding constraint
				C[counter].set(GRB_DoubleAttr_RHS, obj_val_master);
			}
		}
		model->update();
		model->write("Generated Formulations/LeximinMaster.lp");

		// Count the number of agents for which no selection probability is fixed yet
		sum = 0;
		for (int i = 0; i < M_remaining.size(); i++) {
			if (M_remaining[i] == 1) {
				sum++;
			}
		}
	}

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	return L;
}

void LeximinMaster::addColumn(solution sol, bool print) {
	GRBColumn col;
	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			// It only makes sense to include the agents in 'M'
			if (sol.x[i] == 1) {
				col.addTerm(1, C[counter]);
			}
			counter++;
		}
	}
	col.addTerm(1, C_Sum1);
	columns.push_back(col);

	K->block_solution(sol);

	// Add variable for this
	char name_w[13];
	sprintf_s(name_w, "w_%i", (int) K->S.size() - 1);
	GRBVar empty;
	w.push_back(empty);
	w.back() = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, col, name_w);
	model->write("Generated Formulations/LeximinMaster.lp");
}

void LeximinMaster::getDualValues(bool print) {
	dual_C.clear();
	dual_C_Sum1.clear();

	// Count the number of agents in M
	int size_M = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			size_M++;
		}
	}

	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			dual_C.push_back(C[counter].get(GRB_DoubleAttr_Pi));
			counter++;
		}
	}
	dual_C_Sum1.push_back(C_Sum1.get(GRB_DoubleAttr_Pi));

	counter = 0;
	if (print) {
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				printf("\t\tDual_C_%i = %.2f\n", i, dual_C[counter]);
				counter++;
			}
		}
		printf("\t\tDual_C_Sum1 = %.2f\n", dual_C_Sum1[0]);
	}
}

void LeximinMaster::defineModelConVar(bool print) {
	env = new GRBEnv();
	model = new GRBModel(*env);
	env->set(GRB_StringParam_LogFile, "Generated Log-files/LeximinMaster.log");

	// Define 'MIN_M' variable, which we will maximize
	MIN_M = GRBVar();
	char name[13];
	sprintf_s(name, "MIN_M");
	MIN_M = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, name);

	// First add all the constraints needed
	C = new GRBConstr[K->M.size()];
	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			C[counter] = model->addConstr(-MIN_M >= 0.0);
			counter++;
		}
	}

	C_Sum1 = GRBConstr();
	GRBLinExpr lin = 0.0;
	C_Sum1 = model->addConstr(lin == 1);

	// Now define the columns
	columns = std::vector<GRBColumn>();
	columns.resize(K->S.size());
	for (int s = 0; s < K->S.size(); s++) {
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				// It only makes sense to include the agents in 'M'
				if (K->S[s].x[i] == 1) {
					columns[s].addTerm(1, C[counter]);
				}
				counter++;
			}
		}
		// And for the constraint that forces the sum of the weights to be equal to one
		columns[s].addTerm(1, C_Sum1);

		K->block_solution(K->S[s]);
	}

	model->write("Generated Formulations/LeximinMaster.lp");
	// Define variables
	w = std::vector<GRBVar>(K->S.size());
	for (int i = 0; i < K->S.size(); i++) {
		char name_w[13];
		sprintf_s(name_w, "w_%i", i);
		w[i] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, columns[i], name_w);
	}

	GRBLinExpr obj = MIN_M;
	model->setObjective(obj, GRB_MAXIMIZE);

	dual_C = std::vector<double>();
	dual_C_Sum1 = std::vector<double>();

	if (print) {
		model->write("Generated Formulations/LeximinMaster.lp");
	}
}

LeximinMaster::LeximinMaster(IPSolver* K_in, bool print) {
	K = new IPSolver(K_in, print);

	defineModelConVar(print);
}

LeximinMaster::~LeximinMaster() {
	delete[] C;
	delete K;
	delete model;
	delete env;
}
