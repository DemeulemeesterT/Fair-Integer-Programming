#include "SchedulingWeightTard.h"

void SchedulingWeightTard::read_data(std::string name, bool print) {
	if (name == "wt40") {
		n = 40;
	}
	else if (name == "wt50") {
		n = 50;
	}
	else if (name == "wt100") {
		n = 100;
	}
	else {
		printf("\n\nThe requested file cannot be found!\n\n");
		return;
	}

	// From the description of the data:
	// Each file contains the data for 125 instances, listed one after the other.The n processing times are
	// listed first, followed by the n weights, and finally n due dates, for each of the
	// 125 instances in turn.

	SchedWT_inst SWT;
	std::ifstream StFile;
	std::string full_name = "Scheduling instances/" + name + ".txt";
	StFile.open(full_name);

	for (int i = 0; i < 125; i++) {
		SWT.d = std::vector<int>();
		SWT.p = std::vector<int>();
		SWT.w = std::vector<int>();

		int value;
		for (int j = 0; j < n; j++) {
			StFile >> value;
			SWT.p.push_back(value);
		}
		for (int j = 0; j < n; j++) {
			StFile >> value;
			SWT.w.push_back(value);
		}
		for (int j = 0; j < n; j++) {
			StFile >> value;
			SWT.d.push_back(value);
		}

		Sched_I.push_back(SWT);
	}

	if (print) {
		printf("\n The scheduling instances have been successfully read.\n");
	}
}


inst SchedulingWeightTard::generate_instance(int nr, bool export_inst, bool print) {
	inst I;

	I.n = n; // One binary variable for each task (scheduled before due date or not)
	I.t = n; // One integer variable for each task (y_j = time on which task j is scheduled)
	I.n_var = I.n + I.t;
	I.Y_bool = false;
	I.Y_coeff_zero = true;

	I.m = I.n; // One constraint for each task

	// Mixed integer programming formulations for single machine scheduling problems by Keha et al. (2009) 
	// propose an overview of MIP formulations for different scheduling problems.
	// Based on their results, we chose to implement formulation using assignment and positional date variables

	// Objective coefficients
		// Minimizing the weights of the tasks that are planned after the due date would result in objective function
		// min \sum_j (w_j * (1 - x_j)) = min \sum_j (w_j - w_j * x_j)
		// which will find the same optimal objective value as min \sum_j -w_j * x_j
		// which is equivalent to max \sum_j w_j * x_j
			// And we maximize the objective function i IPSolver class
	for (int i = 0; i < I.n; i++) {
		I.v.push_back(Sched_I[nr].w[i]);
	}
	for (int i = 0; i < I.t; i++) {
		I.v.push_back(0);
	}

	// For each task we add the following constraint:
	// y_j <= d_j + (1-x_j) * (\sum_i(p_i) - d_j)
		// If x_j = 1 (task j is scheduled before the due date), this implies that y_j <= d_j
		// If x_j = 0 (task j is scheduled after due date), this implies y_j <= \sum_i p_i (which is automatically satisfied)

	// Capacities of the constraints
	for (int i = 0; i < I.n; i++) {
		I.C.push_back(1);
	}
	for (int i = 0; i < I.n; i++) {
		I.C.push_back(0);
	}


	return I;
}


SchedulingWeightTard::SchedulingWeightTard(std::string name, bool print) {
	read_data(name, print);
}


SchedulingWeightTard::~SchedulingWeightTard() {

}