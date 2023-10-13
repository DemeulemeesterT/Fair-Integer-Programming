#include "IPSolver.h"

void IPSolver::partition_cardinal(bool print) {
	clock_t start_time = clock();

	// Add constraint to enforce that we only find solutions with the optimal objective value
	GRBLinExpr lin = 0;
	for (int j = 0; j < I.n; j++) {
		lin += I.v[j] * X[j];
	}
	for (int j = 0; j < I.t; j++) {
		lin += I.v[I.n + j] * Y_var[j];
	}
	GRBConstr OPT_CONSTRAINT = model->addConstr(lin == (double) opt);
	//model->write("Generated Formulations/IPModel.lp");


	GRBLinExpr obj = 0;

	for (int i = 0; i < I.n; i++) {
		obj = X[i];

		// First, find the minimum value of X[i] in all optimal solutions
		model->setObjective(obj, GRB_MINIMIZE);
		//model->write("Generated Formulations/IPModel.lp");

		solve_partition_cardinal(i, true, print);

		// Next, we find the maximum
		model->setObjective(obj, GRB_MAXIMIZE);
		//model->write("Generated Formulations/IPModel.lp");

		solve_partition_cardinal(i, false, print);
	}

	double length = ((double)(clock() - start_time) / CLK_TCK);

	// Remove the constraint to enforce optimal objective value
	//model->write("Generated Formulations/IPModel.lp");
	model->remove(OPT_CONSTRAINT);
	//model->write("Generated Formulations/IPModel.lp");


	// Determine the number of agents in M
		// This is the number of agents for whom Xmin[i] is not equal to Xmax[i]
	M_size = 0;
	for (int i = 0; i < I.n; i++) {
		if (Xmin[i] != Xmax[i]) {
			M[i] = 1;
			M_size++;
		}
		else {
			M[i] = 0;
		}
	}
	
	if (print) {
		if (M_size > 0) {
			printf("\t\t MIN \t\tMAX\n");
			for (int i = 0; i < I.n; i++) {
				if (Xmin[i] != Xmax[i]) {
					printf("X[%i]\t\t%.2f\t\t%.2f\n", i, Xmin[i], Xmax[i]);
				}
			}
			printf("\n %i agents in M.\n", M_size);
		}
	}

	time_partition = length;

	done_partition = true;
}



double IPSolver::solve_partition_cardinal(int k, bool fill_Xmin, bool print) {
	solver_times++;

	model->optimize();
	int status = model->get(GRB_IntAttr_Status);
	if (status == 3) {
		model->computeIIS();
		model->write("Generated Formulations/UFL_IIS.ilp");
	}

	double z = model->get(GRB_DoubleAttr_ObjVal);

	// Store solution
	std::vector<double> x(I.n, 0);
	if (I.X_integer == true) {
		for (int j = 0; j < I.n; j++) {
			x[j] = (int) X[j].get(GRB_DoubleAttr_X);
		}
	}
	else {
		for (int j = 0; j < I.n; j++) {
			x[j] = X[j].get(GRB_DoubleAttr_X);
		}
	}

	std::vector<double> y(I.t, 0);
	if (I.Y_integer == true) {
		for (int j = 0; j < I.t; j++) {
			y[j] = (int) Y_var[j].get(GRB_DoubleAttr_X);
		}
	}
	else {
		for (int j = 0; j < I.t; j++) {
			y[j] = Y_var[j].get(GRB_DoubleAttr_X);
		}
	}

	// Fill in 'Xmin'/'Xmax'
	if (fill_Xmin == true) {
		Xmin[k] = x[k];
	}
	else {
		Xmax[k] = x[k];
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

	if (solver_times == 1) {
		// Initialize 'opt' if solved for the first time
		opt = z;
	}

	// Calculate objective value
	double obj_val = 0;	// We should calculate it ourselves, because when we use this problem as a pricing problem
							// the obtained objective value is different
	for (int i = 0; i < I.n; i++) {
		obj_val += I.v[i] * x[i];
	}
	for (int i = 0; i < I.t; i++) {
		obj_val += I.v[I.n + i] * y[i];
	}

	if (obj_val >= opt - 0.00001) {
		// Add the solution to S if it is not yet in there
		solution s;

		// Add solution to S
		s.x = x;
		s.y = y;
		s.ID = 0; // We won't use the binary represenation of a solution for cardinal solutions

		if (done_partition == false) {
			check_solution_in_S_cardinal(s, print);
		}
	}
	return z;
}

IP_report IPSolver::solve_return_solution_cardinal(bool print) {
	IP_report IP_R;
	solution s;
	bool found = false;

	std::vector<double> x(I.n, 0);
	std::vector<double> y(I.t, 0);
	
	solver_times++;
	model->write("Generated Formulations/IPModel.lp");
	model->optimize();
	int status = model->get(GRB_IntAttr_Status);
	double z;
	//printf("Status = %i\n", status);
	if (status == 3) {
		model->computeIIS();
		model->write("Generated Formulations/UFL_IIS.ilp");
	}
	if (status != 3) { // If feasible
		z = model->get(GRB_DoubleAttr_ObjVal);
		printf("Optimal solution = %.2f\n", z);
		// Store solution
		if (I.X_integer == true) {
			for (int j = 0; j < I.n; j++) {
				x[j] = (int)X[j].get(GRB_DoubleAttr_X);
			}
		}
		else {
			for (int j = 0; j < I.n; j++) {
				x[j] = X[j].get(GRB_DoubleAttr_X);
			}
		}

		if (I.Y_integer == true) {
			for (int j = 0; j < I.t; j++) {
				y[j] = (int)Y_var[j].get(GRB_DoubleAttr_X);
			}
		}
		else {
			for (int j = 0; j < I.t; j++) {
				y[j] = Y_var[j].get(GRB_DoubleAttr_X);
			}
		}

		if (print) {
			for (int j = 0; j < I.n; j++) {
				if (x[j] > 0) {
					printf("\tX_%i = %.2f\n", j, (double)x[j]);
				}
			}
			for (int j = 0; j < I.t; j++) {
				if (y[j] > 0) {
					printf("\tY_%i = %.2f\n", j, (double)y[j]);
				}
			}
		}

		if (solver_times == 1) {
			// Initialize 'opt' if solved for the first time
			opt = z;
		}

		// Calculate objective value
		double obj_val = 0;	// We should calculate it ourselves, because when we use this problem as a pricing problem
								// the obtained objective value is different
		for (int i = 0; i < I.n; i++) {
			obj_val += I.v[i] * x[i];
		}
		for (int i = 0; i < I.t; i++) {
			obj_val += I.v[I.n + i] * y[i];
		}

		// Add the solution to S if it is not yet in there
		//solution s;

		// Add solution to solution report
		IP_R.s.x = x;
		IP_R.s.y = y;
		IP_R.s.ID = 0; // We won't use the binary represenation of a solution for cardinal solutions

		found = false;
		if (done_partition == false) {
			found = check_solution_in_S_cardinal(s, print);
		}
		if (found == false) {
			S.push_back(IP_R.s);
		}

	}
	else {
		z = 11223344;
		found = true;
	}

	done_solve = true;

	/*
	for (int i = 0; i < I.n; i++) {
		IP_R.s.x[i] = s.x[i];
	}
	for (int i = 0; i < I.t; i++) {
		IP_R.s.y[i] = s.y[i];
	}*/
	//IP_R.s = s;
	IP_R.opt_obj_value = z;
	IP_R.solution_already_in_S = found;
	return IP_R;
}

IP_report IPSolver::solve_return_solution_MENU(bool print) {
	if (this->I.X_bool == true) {
		return solve_return_solution(print);
	}
	else {
		return solve_return_solution_cardinal(print);
	}
}



bool IPSolver::check_solution_in_S_cardinal(solution s, bool print) {
	// This step is only needed during the partitioning phase, and no longer during, e.g., leximin

	// Go through all existing optimal solutions to check whether solution is already in 'S'
	bool found = false;
	for (int i = 0; i < S.size(); i++) {
		bool identical = true;
		for (int j = 0; j < I.n; j++) {
			if (s.x[j] != S[i].x[j]) {
				identical = false;
				j = I.n;
			}
		}

		if (identical == true) {
			found = true;
			i = S.size();
		}
	}

	if (found == false) {
		S.push_back(s);
	}

	return found;
}