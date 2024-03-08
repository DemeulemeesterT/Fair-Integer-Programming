#include"Run_configurations.h"

void Method_simulations(int iter, int n_per_method, std::string s, bool print) {
	std::mt19937 generator((unsigned)1);

	// Create random seed vector
	std::vector<int> seeds_vector(iter*n_per_method+1,0);
	std::uniform_int_distribution<int> distr_A(0, 100000);
	for (int i = 0; i < iter*n_per_method + 1; i++) {
		seeds_vector[i] = distr_A(generator);
	}

// Run RSD once
	KidneyExchange* KE = new KidneyExchange("70-instance-15", true);
	int n = KE->n_pairs;

	inst I = KE->generate_instance(true, false);

	IPSolver* K = new IPSolver(I, false);

	K->analyze(false);

	LotteryDesigner* L = new LotteryDesigner(K, false);

	std::vector<lottery> results_RSD = L->compare_methods(s, n_per_method, false, false, seeds_vector[0]);
	std::vector<double> distr_RSD(n, 0.0);
	for (int t = 0; t < n; t++) {
		distr_RSD[t] = results_RSD[0].p[t];
	}

	std::vector<std::vector<double>> max_diff(iter, std::vector<double>(n_per_method, 0.0));
		// One row for each iteration, 
		// One column for each extra sample used in RSD.

	std::vector<lottery> results;
	std::vector<double> distr_sample(n, 0.0);

	// Now do this multiple times
	for (int i = 0; i < iter; i++) {
		for (int j = 0; j < n_per_method; j++) {
			results.push_back(L->compare_methods(s, 1, false, false, seeds_vector[i * n_per_method + j + 1]).back());
			//results.push_back(L->RSD_once(false, seeds_vector[i * n_per_RSD + j + 1]));

			// Distribution after j runs of RSD
			for (int t = 0; t < n; t++) {
				distr_sample[t] = (double)((j * distr_sample[t] +  results.back().p[t])/(j+1));
			}

			// Maximum difference between 'distr_sample' and original RSD
			double max = 0.0;
			double diff = 0.0;
			for (int t = 0; t < n; t++) {
				diff = abs(distr_RSD[t] - distr_sample[t]);
				if (diff > max) {
					max = diff;
				}
			}

			max_diff[i][j] = max;
		}
	}

	// Write away this result
	char filename_front[50];
	char filename_back[50];
	std::string filename_distr = "_";
	std::string filename;
	sprintf_s(filename_front, "Generated Results/Simulations");
	sprintf_s(filename_back, ".csv");
	std::string letter;
	std::ofstream O;

	if (s == "R") {
		filename_distr = "_RSD";
	}
	if (s == "H") {
		filename_distr = "_RSDheuristic";
	}
	if (s == "U") {
		filename_distr = "_UNIFORM";
	}
	if (s == "O") {
		filename_distr = "_RE-INDEX";
	}
	if (s == "P") {
		filename_distr = "PERTURB";
	}
	filename = filename_front + filename_distr + filename_back;

	O.open(filename);
	// First column contains number of RSD samples considered - 1
	// All remaining columns contain the maximum difference for the different iterations
	for (int j = 0; j < n_per_method; j++) {
		O << j << ",";
		for (int i = 0; i < iter; i++) {
			O << max_diff[i][j] << ",";
		}
		O << "\n";
	}
	O.close();
}