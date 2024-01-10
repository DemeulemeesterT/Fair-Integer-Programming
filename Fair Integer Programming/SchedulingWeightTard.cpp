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
		SWT.r = std::vector<int>(n, 0);

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


inst SchedulingWeightTard::generate_instanceTIF(int nr, bool export_inst, bool print) {
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
	I.X_bool = true;
	I.X_integer = true;

	I.m = 3 * n + T;

	// Mixed integer programming formulations for single machine scheduling problems by Keha et al. (2009) 
	// https://doi.org/10.1016/j.cie.2008.06.008
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

inst SchedulingWeightTard::generate_instanceTIF_WT(int nr, bool release_dates, bool export_inst, bool print) {
	return generate_instanceTIF_WT(Sched_I[nr], release_dates, export_inst, print);
}

inst SchedulingWeightTard::generate_instanceTIF_WT(SchedWT_inst SI, bool release_dates, bool export_inst, bool print) {
	inst I;

	I.n = n; // One variable for each task to represent the utiliy an agent receives from that job:
		// This utility is the inverse of the tardiness.
		// Job before deadline: utility = 0
		// Job after deadline: utility = minus tardiness

	int T = 0;
	for (int i = 0; i < I.n; i++) {
		T = T + SI.p[i];
	}

	if (release_dates == true) {
		// This means we have to take into account release dates
		// Find largest release date
		int r_max = 0;
		for (int i = 0; i < I.n; i++) {
			if (SI.r[i] > r_max) {
				r_max = SI.r[i];
			}
		}

		// Easy estimate of the maximum possible value of T:
			// Plan all jobs from the largest release time
		T = r_max + T;
	}

	I.t = T * n; // 'T' binary variables for each task, with T = \sum_j p_j
		// y_jt = 1 if task j starts at time slot t.
	I.n_var = I.n + I.t;
	I.Y_bool = true;
	I.Y_coeff_zero = true;
	I.X_bool = false;
	I.X_integer = true;

	// Count how many y-variables are before the release dates;
	int y_counter = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < T; j++) {
			if (j < SI.r[i]) {
				y_counter++;
			}
			else {
				j = T;
			}
		}
	}

	I.m = 4 * n + T + y_counter;

	// Mixed integer programming formulations for single machine scheduling problems by Keha et al. (2009) 
	// https://doi.org/10.1016/j.cie.2008.06.008
	// propose an overview of MIP formulations for different scheduling problems.
	// Based on their results, we chose to implement formulation using Time index variables [F2]

	// Objective coefficients
		// Maximizing the weighted utility here is equivalent to minimizing the total weighted tardiness
		// min \sum_j (-w_j * x_j)
		// OR: max \sum_j (w_j * x_j)
		// IPSolver class uses maximization function, so w_j

	for (int i = 0; i < I.n; i++) {
		I.v.push_back(SI.w[i]);
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
			I.A[ConNr + 2 * i][counter + t] = 1;
			I.A[ConNr + 2 * i + 1][counter + t] = -1;
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
			if (t - SI.p[j] + 1 > 0) {
				start_s = t - SI.p[j] + 1;
			}
			for (int s = start_s; s < t + 1; s++) {
				I.A[ConNr][I.n + j * T + s] = 1;
			}
		}
		I.C.push_back(1);
		ConNr++;
	}


	// For each task we add the following constraint:
	// \sum_t(t * y_jt) + p_j <= d_j - x_j
	// \sum_t(t * y_jt) + x_j <= d_j - p_j

	for (int i = 0; i < I.n; i++) {
		I.A[ConNr][i] = 1;
		for (int t = 0; t < T; t++) {
			I.A[ConNr][I.n + i * T + t] = t;
		}
		ConNr++;
		I.C.push_back(SI.d[i] - SI.p[i]);
	}
	
	// BOUNDS
	// For each job set: x_j <= 0
	for (int i = 0; i < I.n; i++) {
		I.A[ConNr][i] = 1;
		I.C.push_back(0);
		ConNr++;
	}

	// And fix that job cannot start before its release date
	// y_jt = <= 0 if t <= r_j - 1
	for (int i = 0; i < I.n; i++) {
		for (int t = 0; t < T; t++) {
			if (t < SI.r[i]) {
				// Before release date, so cannot start
				int index = I.n + T * i + t;
				I.A[ConNr][index] = 1;
				I.C.push_back(0);
				ConNr++;
			}
			else {
				t = T;
			}
		}
	}

	I.data_name = data_name;

	return I;
}


inst SchedulingWeightTard::generate_data_and_instance_TIF_WT(SchedWT_param param, double beta, bool release_dates, bool export_inst, bool print) {
	// Initate random generators
	std::mt19937 generator((unsigned)param.seed);
	//std::uniform_int_distribution<int> distr_w(0, 120);
	std::uniform_int_distribution<int> distr_r(0, (param.n_jobs - 6) * param.common_process_time);
	std::uniform_int_distribution<int> distr_d(param.common_process_time, (param.n_jobs - 5) * param.common_process_time);
	
	// Baptiste et al. (2004)
	double T = 0;
	std::uniform_int_distribution<int> distr_w(1, 1);
	std::uniform_int_distribution<int> distr_p(1, 10);
	


	SchedWT_inst S;
	n = param.n_jobs;
	S.p = std::vector<int>(n, param.common_process_time); // All jobs the same processing time
	S.d = std::vector<int>(n, 0);
	S.w = std::vector<int>(n, 0);
	S.r = std::vector<int>(n, 0);
	data_name = param.name;

	if (release_dates == false) {
		for (int i = 0; i < n; i++) {
			//S.w[i] = (int) distr_w(generator);
			S.w[i] = 1;
			//S.d[i] = distr_d(generator);
			S.r[i] = 0;

			// Baptiste et al. (2004)
			S.p[i] = distr_p(generator);
			
		}
		T = 0;
		for (int i = 0; i < n; i++) {
			T += S.p[i];
		}
		std::uniform_int_distribution<int> distr_d(0, beta * T);
		for (int i = 0; i < n; i++) {
			S.d[i] = S.p[i] + distr_d(generator);
		}

	}

	else {
		// ALL SAME STARTING TIME AT GENERATED RELEASE TIME
		for (int i = 0; i < n; i++) {
			//S.w[i] = (int) distr_w(generator);
			S.w[i] = 1;
			S.r[i] = distr_r(generator);
			S.d[i] = distr_d(generator) + S.r[i];
		}
	}

	Sched_I.push_back(S);

	if (export_inst) {
		export_data(S, print);
	}

	return generate_instanceTIF_WT(S, release_dates, export_inst, print);
}


inst SchedulingWeightTard::generate_instanceLawlerMoore(int nr, bool export_inst, bool print) {
	inst I;
	SchedWT_inst S = order_jobs(nr, print); // Order the jobs in increasing due date

	I.n = n; // One binary variable for each task (scheduled before due date or not)

	I.t = 0; // No auxiliary variables

	I.n_var = I.n + I.t;
	I.Y_bool = true;
	I.Y_coeff_zero = true;

	I.m = n;

	// Objective coefficients
		// Minimizing the weights of the tasks that are planned after the due date would result in objective function
		// min \sum_j (w_j * (1 - x_j)) = min \sum_j (w_j - w_j * x_j)
		// which will find the same optimal objective value as min \sum_j -w_j * x_j
		// which is equivalent to max \sum_j w_j * x_j
			// And we maximize the objective function i IPSolver class
	for (int i = 0; i < I.n; i++) {
		I.v.push_back(S.w[i]);
	}
	for (int i = 0; i < I.t; i++) {
		I.v.push_back(0);
	}

	// Initialize constraint matrix
	I.A = std::vector<std::vector<double>>(I.m, std::vector<double>(I.n_var, 0.0));

	// Constraint
	// for each job k:
	// \sum_j=1^k p_j * x_j <= d_k
		// Where the jobs are ranked from earliest to latest due date
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i + 1; j++) {
			I.A[i][j] = S.p[j];
		}
		I.C.push_back(S.d[i]);
	}

	return I;
}

SchedWT_inst SchedulingWeightTard::order_jobs(int nr, bool print) {
	SchedWT_inst S;
	S.w = std::vector<int>();
	S.p = std::vector<int>();
	S.d = std::vector<int>();

	// Go through the jobs and put them in order of increasing due date
	bool inserted;
	for (int i = 0; i < n; i++) {
		inserted = false;
		if (S.d.size() > 0) {
			for (int j = 0; j < i; j++) {
				if (Sched_I[nr].d[i] < S.d[j]) {
					// Earlier due date, insert before 
					S.d.insert(S.d.begin() + j, Sched_I[nr].d[i]);
					S.p.insert(S.p.begin() + j, Sched_I[nr].p[i]);
					S.w.insert(S.w.begin() + j, Sched_I[nr].w[i]);

					j = i;
					inserted = true;
				}
			}

			// Maybe add it as the last element if it is bigger than all others
			if (inserted == false) {
				S.d.insert(S.d.begin() + i, Sched_I[nr].d[i]);
				S.p.insert(S.p.begin() + i, Sched_I[nr].p[i]);
				S.w.insert(S.w.begin() + i, Sched_I[nr].w[i]);
			}
		}
		else {
			S.d.push_back(Sched_I[nr].d[i]);
			S.w.push_back(Sched_I[nr].w[i]);
			S.p.push_back(Sched_I[nr].p[i]);
		}
	}

	return S;
}


void SchedulingWeightTard::export_data(SchedWT_inst S_in, bool print) {
	std::ofstream S;
 	std::string export_name = "Scheduling instances/Generated/" + data_name + ".dat";
	S.open(export_name);
	S << n << "\n";
	// First, list the processing times
	for (int i = 0; i < n; i++) {
		S << S_in.p[i] << " ";
	}
	S << "\n";

	// Then, list the due dates
	for (int i = 0; i < n; i++) {
		S << S_in.d[i] << " ";
	}
	S << "\n";

	S.close();
}

SchedulingWeightTard::SchedulingWeightTard(std::string name, bool print) {
	read_data(name, print);
}

SchedulingWeightTard::SchedulingWeightTard(bool print) {

}



SchedulingWeightTard::~SchedulingWeightTard() {

}