#include "LeximinMaster.h"

lottery LeximinMaster::solve(bool print) {
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);      //comment to see the output of the solver
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
	K->model->write("Generated Formulations/IPModel.lp");
	int iterations = 1;

	// Create a lottery object corresponding to store the solution
	lottery L;

	// We will use two while loops to find the leximin_legacy solution.
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
		bool all_solutions_in_S = false;
		// This boolean will become true if all optimal solutions are included in K->S.

		// As long as the solution of the pricing problem has a negative reduced cost, continue
		double obj_val_pricing = -1;
		double obj_val_master;
		while (obj_val_pricing < -0.0001) {
			
			if (print) {
				printf("\nITERATION %i\n\n", iterations);
			}

			GRBLinExpr obj = MIN_M;
			model->setObjective(obj, GRB_MAXIMIZE);

			//model->write("Generated Formulations/LeximinModel.lp");

			model->optimize();
			model->write("Generated Formulations/LeximinMaster.lp");
			obj_val_master = model->get(GRB_DoubleAttr_ObjVal);
			getDualValues(false, print);

			// Change the objective of the knapsack problem, which will serve as our pricing problem
				// To allow for non-binary X-variables, we write the full formulation
				// This means that we divide the X-variables by (Xmax[i] - Xmin[i]) for the agents in that were not yet fixed.
			obj = 0.0;

			int counter = 0;
			for (int i = 0; i < K->I.n; i++) {
				if (K->M[i] == 1) {
					if (M_remaining[i] == 1) {
						obj += dual_C_bound[counter] * K->X[i] / (K->Xmax[i] - K->Xmin[i]);
					}
					else {
						obj += dual_C_bound[counter] * K->X[i];
					}
					counter++;
				}
			}
			obj += dual_C_Sum1[0];
			K->model->setObjective(obj, GRB_MINIMIZE);
			K->model->write("Generated Formulations/IPModel.lp");

			//bool correct_master = check_reduced_cost_S_zero(M_remaining, print);


			IP_report IP_R_pricing = K->solve_return_solution_MENU(false);
			obj_val_pricing = IP_R_pricing.opt_obj_value;
			
			// If the model was not infeasible
			if (K->model->get(GRB_IntAttr_Status) == 3) {
				all_solutions_in_S = true;
				obj_val_pricing = 0; // This will stop the master-pricing iterations
			}
			else { // model was feasible
				// Add a new column to the master if the reduced cost is negative
				if (obj_val_pricing < -0.0001) {
					if (IP_R_pricing.solution_already_in_S == false) {
						addColumn(IP_R_pricing.s, print);
					}
					else {
						// Don't add a new 'w' variable, but block only block the found solution
						K->block_solution(IP_R_pricing.s);
					}
				}
				else {
					// This means that the solution to the primal problem is optimal
					if (IP_R_pricing.solution_already_in_S == false) {
						// Remove the last added solution to K->S, because we will not create a variable for this
						K->S.pop_back();
					}
				}
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

		std::vector<double> w_store(w.size(),0.0);
		for (int i = 0; i < w.size(); i++) {
			w_store[i] = (w[i].get(GRB_DoubleAttr_X));
		}
		

		// Store the solution in 'L'
		L.p = std::vector<double>(K->I.n, 0);
		L.w = std::vector<double>();
		L.S = std::vector<solution>();
		for (int i = 0; i < w.size(); i++) { 
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
		// To determine this, we will again have to solve a series of column generation frameworks
			// This time, one for each agent for whom the selection probabilities are not yet fixed.
		GRBLinExpr obj = Epsilon;
		model->setObjective(obj, GRB_MAXIMIZE);
		int counter = -1;

		// Add constraint to force that 'MIN_M'= found objective value
		C_bound_MIN_M = model->addConstr(MIN_M == obj_val_master);
		model->write("Generated Formulations/LeximinMaster.lp");
		double epsilon_opt = 0.0;

		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				counter++;
				// To find the index of the corresponding constraint
			}
			if (M_remaining[i] == 1) {
				// Change  the corresponding constraint by adding a small epsilon value
				// We will do this by setting the coefficient of the column 'Epsilon' to minus one for that constraint

				double obj_val_pricing_epsilon = -1;
				int iterations_subroutine = 1;
				while (obj_val_pricing_epsilon < -0.0001) {
					if (print) {
						printf("\nITERATION SUBROUTINE %i.%i.%i\n\n", iterations-1, i, iterations_subroutine);
					}

					model->chgCoeff(C_bound[counter], Epsilon, -1.0);
					model->write("Generated Formulations/LeximinModel.lp");
					model->optimize();

					// If optimal solution of Epsilon = 0
						// The leximin selection probability of agent i should be equal to MIN_M
					epsilon_opt = model->get(GRB_DoubleAttr_ObjVal);
					if (epsilon_opt < 0.001) {
						// This means that we are not yet sure whether we found the optimal value of epsilon
						// We will continue looking until we have proof of optimality

							// Do the column generation step to decide whether the found solution is optimal
								// or whether we have to add an additional column and solve again.
							getDualValues(true, print);

							// Change the objective of the knapsack problem, which will serve as our pricing problem
							obj = 0.0;

							// The dual variable of the epsilon-constraint:
							//double test_value = dual_C_bound_MIN_M[0] / (K->Xmax[i] - K->Xmin[i]);
							//obj += dual_C_bound_MIN_M[0] * K->X[i] / (K->Xmax[i] - K->Xmin[i]);

							int counter_inner = 0;
							for (int j = 0; j < K->I.n; j++) {
								if (K->M[j] == 1) {
									if (M_remaining[i] == 1) {
										//test_value = dual_C_bound[counter_inner] / (K->Xmax[j] - K->Xmin[j]);
										obj += dual_C_bound[counter_inner] * K->X[j] / (K->Xmax[j] - K->Xmin[j]);
									}
									else {
										obj += dual_C_bound[counter_inner] * K->X[j];
									}
									counter_inner++;
								}
							}
							obj += dual_C_Sum1[0];
							K->model->setObjective(obj, GRB_MINIMIZE);
							K->model->write("Generated Formulations/IPModel.lp");

							//bool correct_epsilon = check_reduced_cost_S_zero(M_remaining, print);

							IP_report IP_R = K->solve_return_solution_MENU(false);
							obj_val_pricing_epsilon = IP_R.opt_obj_value;
							//obj_val_pricing_epsilon = K->solve(true);

							// If the model was infeasible
							if (K->model->get(GRB_IntAttr_Status) == 3) {
								all_solutions_in_S = true;
								// This means that all optimal solutions are included in the set K->S
								// No more optimal solutions can be found

								// Stop the master-pricing iterations for the epsilon problems
								obj_val_pricing_epsilon = 0;

								// Epsilon*=0 is the optimal solution
									// This means that we have to fix the selection probabilty of agent i
										// to the obtained value of MIN_M

								// We can simply do this by bounding the RHS of the bounding inequality to the optimal objective value of MIN_M...
								C_bound[counter].set(GRB_DoubleAttr_RHS, obj_val_master);
								C_bound[counter].set(GRB_CharAttr_Sense, '=');

								// And by setting the coefficient of MIN_M in that constraint to zero
								model->chgCoeff(C_bound[counter], MIN_M, 0.0);
								//model->write("Generated Formulations/LeximinModel.lp");

								M_remaining[i] = 0;
							}
							
							else {// model was feasible
								// Add a new column to the master if the reduced cost is non-negative
								if (obj_val_pricing_epsilon < -0.0001) {
									if (IP_R.solution_already_in_S == false) {
										addColumn(IP_R.s, print);
									}
									else {
										// Don't add a new 'w' variable, but block only block the found solution
										K->block_solution(IP_R.s);
									}
								}

								else {
									// Epsilon*=0 is the optimal solution
									// This means that we have to fix the selection probabilty of agent i
										// to the obtained value of MIN_M

									// We can simply do this by bounding the RHS of the bounding inequality to the optimal objective value of MIN_M...
									C_bound[counter].set(GRB_DoubleAttr_RHS, obj_val_master);
									C_bound[counter].set(GRB_CharAttr_Sense, '=');

									// And by setting the coefficient of MIN_M in that constraint to zero
									model->chgCoeff(C_bound[counter], MIN_M, 0.0);
									model->write("Generated Formulations/LeximinModel.lp");

									M_remaining[i] = 0;

									if (IP_R.solution_already_in_S == false) {
										// Lastly, we remove the last added solution to K->S, because the K->solve() function will automatically store it
										// while we will not create a variable for this solution
										K->S.pop_back();
									}
								}
							}
						
						if (print) {
							// If the model was not infeasible
							if (K->model->get(GRB_IntAttr_Status) != 3) {
								printf("\t Reduced cost = %.5f\n\n", obj_val_pricing_epsilon);
							}
							else {
								printf("\t Pricing problem INFEASIBLE or ALL optimal solutions are found\n\n");
								//model->write("Generated Formulations/LeximinModel.lp");
								//K->model->computeIIS();
								//K->model->write("Generated Formulations/LeximinModal_IIS.ilp");

							}
						}
					}
					else {
						// Although the objective value of epsilon might not be optimal yet
						// We know that it is at least strictly positive
						// And that agent i is therefore selected with a strictly higher probability than MIN_M

						// Terminate the while-loop
						obj_val_pricing_epsilon = 0;
					}
					iterations_subroutine++;
				}
				// Undo the modifications to the model to solve it again in the next iteration
				model->chgCoeff(C_bound[counter], Epsilon, 0.0);
				//model->write("Generated Formulations/LeximinModel.lp");


			}
		}

		// Remove constraint to force that 'MIN_M' is found objective value.
		model->remove(C_bound_MIN_M);
		//model->write("Generated Formulations/LeximinModel.lp");


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
	//model->write("Generated Formulations/LeximinModel.lp");

	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			// It only makes sense to include the agents in 'M'
			if (sol.x[i] == 1) {
				col.addTerm(sol.x[i]/(K->Xmax[i] - K->Xmin[i]), C_bound[counter]);
				// Formulated like this, the leximin formulation can also be used for
					// non-binary X-variables.
					// For binary X-variables, the coefficient was simply 1.
					// And for agents in M, this will still be the case for binary X-variables,
					// as 1/(1-0) = 1.
			}
			counter++;
		}
	}
	col.addTerm(1, C_Sum1);
	columns.push_back(col);
	//model->write("Generated Formulations/LeximinModel.lp");


	K->block_solution(sol);

	// Block the found solution (only the X-variables)
	GRBLinExpr lin = 0;
	counter = 0;
	for (int i = 0; i < K->S.back().x.size(); i++) {
		if (K->M[i] == 1) { // If the agent is in M, that's the only case in which it matters
			if (K->S.back().x[i] == 1) {// If the agent is selected
				lin += K->X[i];
				counter++;
			}
		}
	}
	K->model->addConstr(lin <= (counter - 1));

	// Add variable for this
	char name_w[13];
	sprintf_s(name_w, "w_%i", (int)K->S.size() - 1);
	GRBVar empty;
	w.push_back(empty);
	w.back() = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, col, name_w);
	model->write("Generated Formulations/LeximinMaster.lp");
}


void LeximinMaster::getDualValues(bool epsilon_constraint, bool print) {
	dual_C_bound.clear();
	dual_C_Sum1.clear();
	dual_C_bound_MIN_M.clear();

	// Count the number of agents in M
	/*int size_M = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			size_M++;
		}
	}
	*/

	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			dual_C_bound.push_back(C_bound[counter].get(GRB_DoubleAttr_Pi));
			counter++;
		}
	}
	dual_C_Sum1.push_back(C_Sum1.get(GRB_DoubleAttr_Pi));

	if (epsilon_constraint == true) {
		dual_C_bound_MIN_M.push_back(C_bound_MIN_M.get(GRB_DoubleAttr_Pi));
	}

	counter = 0;
	if (print) {
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				printf("\t\tDual_C_%i = %.2f\n", i, dual_C_bound[counter]);
				counter++;
			}
		}
		printf("\t\tDual_C_Sum1 = %.2f\n", dual_C_Sum1[0]);

		if (epsilon_constraint == true) {
			printf("\t\tDual_C_bound_MIN_M = %.2f\n", dual_C_bound_MIN_M[0]);
		}
	}
}

void LeximinMaster::defineModelConVar(bool print) {
	env = new GRBEnv();
	model = new GRBModel(*env);
	env->set(GRB_StringParam_LogFile, "Generated Log-files/LeximinMaster.log");

	// Define 'MIN_M' variable, which we will maximize iteratively
	MIN_M = GRBVar();

	char name[13];
	sprintf_s(name, "MIN_M");
	MIN_M = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, name);

	// First add all the constraints needed
	// We will initialize the model by adding this constraint for all agents in M
	C_bound = std::vector<GRBConstr>(K->M_size);
	//C_bound = new GRBConstr[K->M.size()];

	// To allow us to use our model for non-binary X-variables, we will add the following constant term to each of these constraints
		// This term will be zero for binary X-variables.
	double RHS_constant = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {

			RHS_constant += K->Xmin[i] / (K->Xmax[i] - K->Xmin[i]);
		}
	}
	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			sprintf_s(name, "C_bound_%i", i);
			C_bound[counter] = model->addConstr(-MIN_M - RHS_constant >= 0.0, name);
			counter++;
		}
	}


	C_Sum1 = GRBConstr();
	C_bound_MIN_M = GRBConstr();

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
				//if (K->S[s].x[i] == 1) {
				double test_value = K->S[s].x[i] / (K->Xmax[i] - K->Xmin[i]);
					columns[s].addTerm(K->S[s].x[i] / (K->Xmax[i] - K->Xmin[i]), C_bound[counter]);
					// Formulated like this, the leximin formulation can also be used for
						// non-binary X-variables.
						// For binary X-variables, the coefficient was simply 1.
						// And for agents in M, this will still be the case for binary X-variables,
						// as 1/(1-0) = 1.
				//}
				counter++;
			}
		}
		// And for the constraint that forces the sum of the weights to be equal to one
		columns[s].addTerm(1, C_Sum1);

		K->block_solution(K->S[s]);
		K->model->write("Generated Formulations/IPModel.lp");

	}

	//Epsilon_column = GRBColumn();
	// We initialize 'Epsilon_column' with zeroes everywhere
	Epsilon = GRBVar();
	sprintf_s(name, "Epsilon");
	Epsilon = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, name);

	//model->write("Generated Formulations/LeximinMaster.lp");
	// Define variables
	w = std::vector<GRBVar>(K->S.size());
	for (int i = 0; i < K->S.size(); i++) {
		char name_w[13];
		sprintf_s(name_w, "w_%i", i);
		w[i] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, columns[i], name_w);
	}

	GRBLinExpr obj = MIN_M;
	model->setObjective(obj, GRB_MAXIMIZE);

	dual_C_bound = std::vector<double>();
	dual_C_Sum1 = std::vector<double>();

	prob_agents = std::vector<double>();
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			prob_agents.push_back(0);
		}
	}

	//if (print) {
		//model->write("Generated Formulations/LeximinMaster.lp");
	//}
}

bool LeximinMaster::check_reduced_cost_S_zero(std::vector<int> M_remaining, bool print) {
	double red_cost;
	bool all_zero = true;
	for (int t = 0; t < K->S.size(); t++) {
		red_cost = 0;
		int counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				if (M_remaining[i] == 1) {
					red_cost += dual_C_bound[counter] * K->S[t].x[i] / (K->Xmax[i] - K->Xmin[i]);
				}
				else {
					red_cost += dual_C_bound[counter] * K->S[t].x[i];
				}
				counter++;
			}
		}
		red_cost += dual_C_Sum1[0];

		if (abs(red_cost - 0.0) > 0.0000001) {
			all_zero = false;
		}
	}
	return all_zero;
}


LeximinMaster::LeximinMaster(IPSolver* K_in, bool print) {
	K = new IPSolver(K_in, print);
	//K->model->write("Generated Formulations/LeximinMaster.lp");

	defineModelConVar(print);
}

LeximinMaster::~LeximinMaster() {
	delete K;
	delete model;
	delete env;
}

