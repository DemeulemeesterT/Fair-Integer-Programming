#include "SimplicalDecomposition.h"


lottery SimplicalDecomposition::SD_Nash(bool print) {
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);      //comment to see the output of the solver
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	printf("\n\n\n For Nash, it is very important that every agent in M is at least selected in one of the optimal solutions that is used to initiate the Simplical Decomposition.\n");
	printf(" Therefore, we will initiate the solution set with a greedy covering of the agents.\n\n");
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	int counter;

	// Create a lottery object corresponding to store the solution
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	// Check if greedy partition has been done, and do if needed.
	if (K->done_greedy_partition == false) {
		K->greedy_partition(false);
	}

	// INTRODUCE SPECIFIC VARIABLES AND CONSTRAINTS TO IMPLEMENT NASH
	GRBVar* log_p = new GRBVar[K->M_size];
	// The logarithms of 'p'
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

	// CONSTRAINTS:
	GRBGenConstr* C_log_p = new GRBGenConstr[K->M_size];

	// Let 'log_p[i]' be the log of 'p[i]'
	counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			C_log_p[counter] = model->addGenConstrLog(p[counter], log_p[counter]);
			counter++;
		}
	}
	//model->write("Generated Formulations/SDNashMaster.lp");

	// OBJECTIVE
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

	//model->write("Generated Formulations/SDNashMaster.lp");

	// Add the columns
	columns.resize(K->S_greedy.size());
	for (int s = 0; s < K->S_greedy.size(); s++) {
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				// It only makes sense to include the agents in 'M'
				if (K->S_greedy[s].x[i] == 1) {
					columns[s].addTerm(1, C_p[counter]);
				}
				counter++;
			}
		}
		// And for the constraint that forces the sum of the weights to be equal to one
		columns[s].addTerm(1, C_Sum1);

		K->block_solution(K->S_greedy[s]);

		// Add the solution to the resulting lottery:
		L.S.push_back(K->S_greedy[s]);
	}

	//model->write("Generated Formulations/SDNashMaster.lp");

	for (int i = 0; i < K->S_greedy.size(); i++) {
		char name_w[13];
		sprintf_s(name_w, "w_%i", i);
		w.push_back(model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, columns[i], name_w));
	}

	//model->write("Generated Formulations/SDNashMaster.lp");

	// Add constraint to knapsack to enforce the objective value to be equal to the optimal objective value
	GRBLinExpr expr = 0.0;

	for (int i = 0; i < K->I.n; i++) {
		expr += K->I.v[i] * K->X[i];
	}
	for (int i = 0; i < K->I.t; i++) {
		expr += K->I.v[K->I.n + i] * K->Y_var[i];
	}
	K->model->addConstr(expr == (double) K->opt);

	int iterations = 1;


	// As long as we keep finding better solutions, keep solving the master problem
	double gradient_improvement = 1;
	double obj_val_master = 0;
	double obj_val_pricing = 0;

	clock_t start_time = clock();

	while (gradient_improvement > 0.000001) {
		//model->write("Generated Formulations/SDNashMaster.lp");
		if (print) {
			printf("ITERATION %i\n", iterations);
		}
		model->optimize();
		obj_val_master = model->get(GRB_DoubleAttr_ObjVal);

		// Get the values of the p-variables
		std::vector<double> p_values;
		for (int i = 0; i < dict.size(); i++) {
			p_values.push_back(p[i].get(GRB_DoubleAttr_X));
		}

		// Determine the value of the gradient evaluated at the solution
		std::vector<double> gradient = gradientNash(p_values, print);

		// Change the objective of the knapsack problem, which will serve as our pricing problem
		GRBLinExpr obj = 0.0;
		//GRBQuadExpr obj = 0.0;
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				obj += gradient[counter] * (K->X[i] - p_values[counter]);
				counter++;

				//obj += gradient[counter] * gradient[counter] * K->X[i] * K->X[i] * 0.5;

				// Newton extension, 
				// see p. 229 https://link-springer-com.kuleuven.e-bronnen.be/content/pdf/10.1007%2F978-0-387-74759-0_614.pdf 
				// and  Larsson T, Patriksson M, Rydergren C (1997) Applications
				//of simplicial decomposition with nonlinear column generation to nonlinear network flows.In: Pardalos PM, Hager
				//	WW, Hearn DW(eds) Network Qptimization.Lecture Notes
				//	Economicsand Math Systems.Springer, Berlin, pp 346–373
			}
		}


		K->model->setObjective(obj, GRB_MAXIMIZE);
		//K->model->write("Generated Formulations/IPModel.lp");
		obj_val_pricing = K->solve(false);

		// Compute the gradient_improvement (grad(f(p_values)) * (y - p_values)), where p_values is the solution from the master, and y the solution from the pricing problem
		gradient_improvement = 0;
		for (int i = 0; i < dict.size(); i++) {
			gradient_improvement += gradient[i] * (K->S.back().x[dict[i]] - p_values[i]);
		}

		// Add a new column to the master if the reduced cost is non-negative
		if (gradient_improvement > 0.000001) {
			addColumn(K->S.back(), print);
			L.S.push_back(K->S.back());

			//K->block_solution(K->S.back());
				// This will block entire solution (including agents in Y and y-variables)

			// Block the found solution
			GRBLinExpr lin = 0;
			int counter = 0;
			for (int i = 0; i < K->S.back().x.size(); i++) {
				if (K->M[i] == 1) { // If the agent is in M, that's the only case in which it matters
					if (K->S.back().x[i] == 1) {// If the agent is selected
						lin += K->X[i];
						counter++;
					}
				}
			}
			K->model->addConstr(lin <= (counter - 1));

		}		

		if (print) {
			// If the model was not infeasible
			if (K->model->get(GRB_IntAttr_Status) != 3) {
				printf("\t Gradient improvement = %.5f\n\t", gradient_improvement);
				printf(" Objective function = %.5f\n\n", obj_val_master);

			}
			else {
				printf("\t Pricing problem INFEASIBLE\n\n");
			}
		}
		iterations++;
	}

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	// Store the solution in 'L'
	for (int i = 0; i < w.size(); i++) { // If we use 'K->S.size()' here, then the lastly added solution is also included, but we didn't solve the model with this variable yet
		L.w.push_back(w[i].get(GRB_DoubleAttr_X));
	}

	// Now calculate the selection probabilities of the agents
	for (int i = 0; i < w.size(); i++) {
		for (int j = 0; j < K->I.n; j++) {
			L.p[j] += L.w[i] * L.S[i].x[j];
		}
	}

	delete[] C_log_p;

	if (print) {

	}
	return L;
}

lottery SimplicalDecomposition::Nash_CG(bool print) {
	//model->getEnv().set(GRB_IntParam_OutputFlag, 0);      //comment to see the output of the solver
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	//model->write("Generated Formulations/SDNashMaster.lp");

	if (print) {
		printf("\n\n\n For Nash, it is very important that, for every agent in M,\nthere are at least two optimal solutions in which the utility of the agent differs that are included in the initial subset.\n");
		printf(" Therefore, we will initiate the solution set with the optimal solutions found during the partitioning of the agents.\n\n");

	}
	K->model->getEnv().set(GRB_IntParam_OutputFlag, 0);   //comment to see the output of the solver

	int counter;

	// Create a lottery object corresponding to store the solution
	lottery L;
	L.p = std::vector<double>(K->I.n, 0);
	L.w = std::vector<double>();
	L.S = std::vector<solution>();

	// Check if greedy partition has been done, and do if needed.
	//if (K->done_greedy_partition == false) {
	//	K->greedy_partition(false);
	//}
		// We used this when the X-variables were binary.

	// INTRODUCE SPECIFIC VARIABLES AND CONSTRAINTS TO IMPLEMENT NASH
	GRBVar* log_p = new GRBVar[K->M_size];
	// The logarithms of 'p[i]' - Xmin[i]
	counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			char name_log_p[13];
			sprintf_s(name_log_p, "log_p_%i", i);
			log_p[counter] = model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, name_log_p);
			// Note that each of the variables has objective coefficient 1, as we want to maximize the sum of the log_p variables
			counter++;
		}
	}

	// CONSTRAINTS:
	GRBGenConstr* C_log_p = new GRBGenConstr[K->M_size];

	// Let 'log_p[i]' be the log of 'p[i]'
	counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			C_log_p[counter] = model->addGenConstrLog(p[counter], log_p[counter]);
			counter++;
		}
	}
	//model->write("Generated Formulations/SDNashMaster.lp");

	// OBJECTIVE
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

	//model->write("Generated Formulations/SDNashMaster.lp");

	// Add the columns
	columns.resize(K->S.size());
	for (int s = 0; s < K->S.size(); s++) {
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				// It only makes sense to include the agents in 'M'
				//if (K->S[s].x[i] > K->Xmin[i]) {
					columns[s].addTerm(K->S[s].x[i], C_p[counter]);
				//}
				counter++;
			}
		}
		// And for the constraint that forces the sum of the weights to be equal to one
		columns[s].addTerm(1, C_Sum1);

		//K->block_solution(K->S_greedy[s]);

		// Add the solution to the resulting lottery:
		L.S.push_back(K->S[s]);
	}

	//model->write("Generated Formulations/SDNashMaster.lp");

	for (int i = 0; i < K->S.size(); i++) {
		char name_w[13];
		sprintf_s(name_w, "w_%i", i);
		w.push_back(model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, columns[i], name_w));
	}

	//model->write("Generated Formulations/SDNashMaster.lp");

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


	// As long as we keep finding better solutions, keep solving the master problem
	double gradient_improvement = 1;
	double obj_val_master = 0;
	double obj_val_pricing = 0;

	clock_t start_time = clock();

	bool finished = false;
	while (finished == false) {
		model->write("Generated Formulations/SDNashMaster.lp");
		if (print) {
			printf("ITERATION %i\n", iterations);
		}
		model->optimize();
		int status = model->get(GRB_IntAttr_Status);
		if (status == 3) {
			model->computeIIS();
			model->write("Generated Formulations/Nash_IIS.ilp");
		}
		
		obj_val_master = model->get(GRB_DoubleAttr_ObjVal);

		// Get the values of the p-variables
		std::vector<double> p_values;
		for (int i = 0; i < dict.size(); i++) {
			p_values.push_back(p[i].get(GRB_DoubleAttr_X));
		}

		// Get the values of the w-variables
		std::vector<double> w_values;
		for (int i = 0; i < L.S.size(); i++) {
			w_values.push_back(w[i].get(GRB_DoubleAttr_X));
		}

		double sum_check = 0;
		for (int i = 0; i < w_values.size(); i++) {
			if (w_values[i] > 0) {
				printf("Solution %i\t%.2f\n", i, w_values[i]);
				sum_check += w_values[i];
			}
		}
		printf("Sum = %.4f\n\n", sum_check);

		// Determine the value of the gradient evaluated at the solution
		std::vector<double> gradient = gradientNash(p_values, print);

		// Change the objective of the knapsack problem, which will serve as our pricing problem
		GRBLinExpr obj = 0.0;
		//GRBQuadExpr obj = 0.0;
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				obj += gradient[counter] * K->X[i];
				counter++;

				//obj += gradient[counter] * gradient[counter] * K->X[i] * K->X[i] * 0.5;

				// Newton extension, 
				// see p. 229 https://link-springer-com.kuleuven.e-bronnen.be/content/pdf/10.1007%2F978-0-387-74759-0_614.pdf 
				// and  Larsson T, Patriksson M, Rydergren C (1997) Applications
				//of simplicial decomposition with nonlinear column generation to nonlinear network flows.In: Pardalos PM, Hager
				//	WW, Hearn DW(eds) Network Qptimization.Lecture Notes
				//	Economicsand Math Systems.Springer, Berlin, pp 346–373
			}
		}


		K->model->setObjective(obj, GRB_MAXIMIZE);
		K->model->write("Generated Formulations/IPModel.lp");
		IP_report IP_R = K->solve_return_solution_MENU(print);
			// IP_report contains only the optimal objective function and an optimal solution (see definition IPSolver.h)

		obj_val_pricing = IP_R.opt_obj_value;

		// Store solution
		if (K->model->get(GRB_IntAttr_Status) != 3) { // If pricing problem not infeasible

			// Check for all solutions that are already included in the master problem what their objective value for this gradient would be
			// We only need to check for solutions with a positive weight w
			double best_value = -1e30;
			int j = 0;
			while (j < L.S.size()) {
				if (w_values[j] > 0) {
					double sol_obj = 0.0;
					counter = 0;
					for (int k = 0; k < K->I.n; k++) {
						if (K->M[k] == 1) {
							//if (L.S[j].x[k] > K->Xmin[k]) {
								sol_obj += gradient[counter] * (L.S[j].x[k] - K->Xmin[k]);
								if (print) {
									printf("\t Existing solution obj: %.4f\n", sol_obj);
								}
							//}
							counter++;
						}
					}

					if (sol_obj > best_value) {
						best_value = sol_obj;
					}
					if (best_value > obj_val_pricing - 0.0001) {
						finished = true;
						j = L.S.size(); // Finish the while loop and the search for additional solutions, the optimal solution has been found
						if (print) {
							printf("\n\n The optimal solution has been found.\n");
						}
					}
				}

				j++;
			}

			// Add a new column to the master if the reduced cost is non-negative
			if (finished == false) {
				addColumn(IP_R.s, print);
				L.S.push_back(IP_R.s);
				//addColumn(K->S.back(), print);
				//L.S.push_back(K->S.back());

				//K->block_solution(K->S.back());
					// This will block entire solution (including agents in Y and y-variables)

				// Block the found solution
				/*GRBLinExpr lin = 0;
				int counter = 0;
				for (int i = 0; i < IP_R.s.x.size(); i++) {
					if (K->M[i] == 1) { // If the agent is in M, that's the only case in which it matters
						if (IP_R.s.x[i] == 1) {// If the agent is selected
							lin += K->X[i];
							counter++;
						}
					}
				}
				K->model->addConstr(lin <= (counter - 1));
				//K->model->update();
				*/
			}
		}
		else {
			finished = true;
		}

		if (print) {
			// If the model was not infeasible
			if (K->model->get(GRB_IntAttr_Status) != 3) {
				printf("\t Objective function pricing = %.5f\n\t", obj_val_pricing);
				printf(" Objective function = %.5f\n\n", obj_val_master);

			}
			else {
				printf("\t Pricing problem INFEASIBLE\n\n");
			}
		}
		iterations++;
	}

	L.t = ((double)(clock() - start_time) / CLK_TCK);

	// Store the solution in 'L'
	for (int i = 0; i < w.size(); i++) { // If we use 'K->S.size()' here, then the lastly added solution is also included, but we didn't solve the model with this variable yet
		L.w.push_back(w[i].get(GRB_DoubleAttr_X));
	}

	// Now calculate the selection probabilities of the agents
	for (int i = 0; i < w.size(); i++) {
		for (int j = 0; j < K->I.n; j++) {
			L.p[j] += L.w[i] * L.S[i].x[j];
		}
	}

	delete[] C_log_p;

	if (print) {

	}
	return L;
}

std::vector<double> SimplicalDecomposition::gradientNash(std::vector<double> p_values, bool print) {
	std::vector<double> result;
	// for log
	for (int i = 0; i < dict.size(); i++) {
		result.push_back((double)(1 / p_values[i]));
	}
	return result;
}

void SimplicalDecomposition::addColumn(solution sol, bool print) {
	GRBColumn col;
	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			// It only makes sense to include the agents in 'M'
			if (sol.x[i] > K->Xmin[i]) {
				col.addTerm(sol.x[i], C_p[counter]);
			}
			counter++;
		}
	}
	col.addTerm(1, C_Sum1);
	columns.push_back(col);

	//K->block_solution(sol);

	// Add variable for this
	char name_w[13];
	sprintf_s(name_w, "w_%i", (int)columns.size());
	GRBVar empty;
	w.push_back(empty);
	w.back() = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, col, name_w);
	//model->write("Generated Formulations/SDNashMaster.lp");
}

SimplicalDecomposition::SimplicalDecomposition(IPSolver* K_in, bool print) {
	K = new IPSolver(K_in, print);
	env = new GRBEnv();
	model = new GRBModel(*env);
	env->set(GRB_StringParam_LogFile, "Generated Log-files/SimplicalDecompositionMaster.log");

	// Check if partitioning is done, and do if needed.
	if (K->done_partition == false) {
		K->partition(false);
	}

	// Create the P-variables
	p = new GRBVar[K->M_size];
	dict = std::vector<int>();
	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			dict.push_back(i);
			char name_p[13];
			sprintf_s(name_p, "p_%i", i);
			// This means we know the lowest and highest values of the X-variables in the optimal solutions
			p[counter] = model->addVar(K->Xmin[i], K->Xmax[i], 0.0, GRB_CONTINUOUS, name_p);
			
			counter++;
		}
	}

	// First add all the constraints needed
	C_Sum1 = GRBConstr();
	GRBLinExpr lin = 0.0;
	C_Sum1 = model->addConstr(lin == 1);

	C_p = new GRBConstr[K->M.size()];
	counter = 0;
	for (int i = 0; i < dict.size(); i++) {
		C_p[counter] = model->addConstr(-p[i] - K->Xmin[i] == 0.0);
			// We compare to the 'dystopia point', the minimum value that agent i receives in any of the optimal solutions.
		counter++;
	}

	// Now define the columns
	columns = std::vector<GRBColumn>();

	//model->write("Generated Formulations/SDNashMaster.lp");
	// Define variables
	w = std::vector<GRBVar>();

	//model->write("Generated Formulations/SDNashMaster.lp");

}

SimplicalDecomposition::~SimplicalDecomposition() {
	delete[] C_p;
	delete[] p;
	delete K;
	delete model;
	delete env;
}