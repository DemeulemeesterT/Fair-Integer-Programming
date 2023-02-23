#include "SchedulingWeightTard.h"

void SchedulingWeightTard::read_data(std::string name, bool print) {
	if (name == "wt3") {
		n = 3;
	}
	else if (name == "wt4") {
		n = 4;
	}
	else if (name == "wt40") {
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

	data_name = name;

	if (print) {
		printf("\n The scheduling instances have been successfully read.\n");
	}
}


inst SchedulingWeightTard::generate_instance(int nr, bool export_inst, bool print) {
	inst I;

	I.n = n; // One binary variable for each task (scheduled before due date or not)

	int T = 0;
	for (int i = 0; i < I.n; i++) {
		T = T + Sched_I[nr].p[i];
	}
	I.t = T*n; // 'T' binary variables for each task, with T = \sum_j p_j
		// y_jt = 1 if task j starts at time slot t.
	I.n_var = I.n + I.t;
	I.Y_bool = true;
	I.Y_coeff_zero = true;

	I.m = 3 * n + T;

	// Mixed integer programming formulations for single machine scheduling problems by Keha et al. (2009) 
	// propose an overview of MIP formulations for different scheduling problems.
	// Based on their results, we chose to implement formulation using Time index variables [F2]

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

	// Initialize constraint matrix
	I.A = std::vector<std::vector<double>>(I.m, std::vector<double>(I.n_var, 0.0));

	// Fill in
	// First two constraints for each job to ensure that each job is scheduled at exactly once
	// \sum_{t=1}^T y_jt = 1 for all jobs j
		// Because this is an equality constraint, we will add two constraint: one where LHS <= 1 and on where -LHS <= -1
	
	// Formula to access the right y-variable in the constraint matrix
	// y_jt = column 
		// I.n (number of x-variables)
		// + T * j (number of previous y-variables)
		// + t (correct period)

	int counter = I.n; // Will keep track of the correct Y-variable
	int ConNr = 0; // Will keep track of the constraint number
	for (int i = 0; i < I.n; i++) {
		for (int t = 0; t < T; t++) { 
			I.A[ConNr + 2*i][counter + t] = 1;
			I.A[ConNr + 2*i + 1][counter + t] = -1;
		}
		// Fix the capacities of the constraints
		I.C.push_back(1);
		I.C.push_back(-1);

		// Change the y-variables we are considering
		counter = counter + T;
	}
	ConNr = 2 * I.n;

	// Then T constraints to ensure that at most one job can start
	// \sum_j \sum_{s = max(0, t-p_j+1}^t y_js <= 1 for all t = 1...T
	for (int t = 0; t < T; t++) {
		for (int j = 0; j < I.n; j++) {
			int start_s = 0;
			if (t - Sched_I[nr].p[j] + 1 > 0) {
				start_s = t - Sched_I[nr].p[j] + 1;
			}
			for (int s = start_s; s < t + 1; s++) {
				I.A[ConNr][I.n + j * T + s] = 1;
			}
		}
		I.C.push_back(1);
		ConNr++;
	}


	// For each task we add the following constraint:
	// \sum_t(t * y_jt) <= d_j - p_j + (1-x_j) * (\sum_i(p_i) - d_j)
		// If x_j = 1 (task j is scheduled before the due date), this implies that \sum_t(t * y_jt) <= d_j - p_j
		// If x_j = 0 (task j is scheduled after due date), this implies \sum_t(t * y_jt) <= \sum_i p_i - p_j (which is automatically satisfied)
	// This can be rewritten as
	// \sum_t (t * y_jt) + (\sum_i p_i - d_j) * x_j <= \sum_i p_i - p_j

	// CARFEFUL: because the indexing starts from zero in this code, but from 1 in the input files, we use d_j - 1 here!
	for (int i = 0; i < I.n; i++) {
		I.A[ConNr][i] = T - (Sched_I[nr].d[i] - 1);
		for (int t = 0; t < T; t++) {
			I.A[ConNr][I.n + i * T + t] = t;
		}
		ConNr++;
		I.C.push_back(T - Sched_I[nr].p[i]);
	}

	I.data_name = data_name;

	return I;
}


SchedulingWeightTard::SchedulingWeightTard(std::string name, bool print) {
	read_data(name, print);
}


SchedulingWeightTard::~SchedulingWeightTard() {

}