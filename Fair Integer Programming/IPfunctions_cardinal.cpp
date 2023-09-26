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

	GRBLinExpr obj = 0;

	for (int i = 0; i < I.n; i++) {
		obj = X[i];

		// First, find the minimum value of X[i] in all optimal solutions
		model->setObjective(obj, GRB_MINIMIZE);
		model->optimize();

		// Store the solution
		Xmin[i] = X[i].get(GRB_DoubleAttr_X);

		// Next, we find the maximum
		model->setObjective(obj, GRB_MAXIMIZE);
		model->optimize();
		Xmax[i] = X[i].get(GRB_DoubleAttr_X);
	}

	double length = ((double)(clock() - start_time) / CLK_TCK);

	// Remove the constraint to enforce optimal objective value
	model->write("Generated Formulations/IPModel.lp");
	model->remove(OPT_CONSTRAINT);
	model->write("Generated Formulations/IPModel.lp");


	// Determine the nunmber of agents in M
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