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
	GRBConstr OPT_CONSTRAINT = model->addConstr(lin == opt);
	model->write("Generated Formulations/IPModel.lp");


	GRBLinExpr obj = 0;

	for (int i = 0; i < I.n; i++) {
		obj = X[i];

		// First, find the minimum value of X[i] in all optimal solutions
		model->setObjective(obj, GRB_MINIMIZE);
		solve_partition_cardinal(i, true, print);

		// Next, we find the maximum
		model->setObjective(obj, GRB_MAXIMIZE);
		solve_partition_cardinal(i, false, print);
	}

	double length = ((double)(clock() - start_time) / CLK_TCK);

	// Remove the constraint to enforce optimal objective value
	model->write("Generated Formulations/IPModel.lp");
	model->remove(OPT_CONSTRAINT);
	model->write("Generated Formulations/IPModel.lp");


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



double IPSolver::solve_partition_cardinal(int i, bool fill_Xmin, bool print) {
	solver_times++;

	model->optimize();

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
		Xmin[i] = x[i];
	}
	else {
		Xmax[i] = x[i];
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

	if (obj_val == opt) {
		// Add the solution to S if it is not yet in there
		solution s;

		// Add solution to S
		s.x = x;
		s.y = y;
		s.ID = 0; // We won't use the binary represenation of a solution for cardinal solutions
		S.push_back(s);
	}
	return z;
}