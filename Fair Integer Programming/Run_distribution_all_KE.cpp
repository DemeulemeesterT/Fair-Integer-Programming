#include"Run_configurations.h"

std::vector<ReportDistribution*> run_distribution_all_KE(std::string s, double time_limit, int iterations, bool print) {
	// Prepare to go through all instances
	std::vector<std::vector<int>> P;
	P.push_back({ 10,50 }); // The first number is the number of agents, the second the ID of the instance
	P.push_back({ 20,50 });
	P.push_back({ 30,50 });
	P.push_back({ 40,65 });
	P.push_back({ 50,50 });
	P.push_back({ 60,50 });
	P.push_back({ 70,50 });

	// Change this conditional on the name of the distribution you're looking for.
	char filename_front[50];
	char filename_back[50];
	std::string filename_distr = "_";
	std::string filename;
	sprintf_s(filename_front, "Generated Results/ReportKE");
	sprintf_s(filename_back, ".csv");
	std::string letter;
	std::ofstream O;
	std::ofstream Q;

	for (int k = 0; k < s.size(); k++) {
		letter = s[k];
		// Decide the correct name of the file
		if (s[k] == 'C') {
			filename_distr = "_NASH";
		}
		if (s[k] == 'R') {
			filename_distr = "_RSD";
		}
		if (s[k] == 'L') {
			filename_distr = "_LEXIMAX";
		}
		if (s[k] == 'U') {
			filename_distr = "_UNIFORM";
		}
		if (s[k] == 'O') {
			filename_distr = "_RE-INDEX";
		}
		if (s[k] == 'P') {
			filename_distr = "PERTURB";
		}
		if (s[k] == 'T') {
			filename_distr = "RSD-ONCE";
		}
		filename = filename_front + filename_distr + filename_back;

		O.open(filename);
		O << "Name, n, |M|, time solution, time partition, time distribution, generated solutions, used solutions, geometric mean, arithmetic mean, minimum selection probability" << "\n";
		O.close();

		// Now create a separate file that contains the distributions themselves
		std::string filename_long = filename_front + filename_distr + "_distr" + filename_back;
		Q.open(filename_long);
		Q << "Name, n, |M|, distribution" << "\n";
		Q.close();
	}

	std::vector<ReportDistribution*> R;

	char numstr1[21];
	char numstr2[21];
	std::string name;
	std::string name_instance = "-instance-";
	for (int p = 0; p < P.size(); p++) {
		for (int i = 1; i <= P[p][1]; i++) {
			sprintf_s(numstr1, "%i", P[p][0]);
			sprintf_s(numstr2, "%i", i);
			name = numstr1 + name_instance + numstr2;
			KidneyExchange* KE = new KidneyExchange(name, false);
			inst I = KE->generate_instance(false, false);
			IPSolver* K = new IPSolver(I, false);

			// Find the time to solve one instance
			clock_t start_time = clock();
			K->solve(print);
			double time_solve = ((double)(clock() - start_time) / CLK_TCK);

			// Find the time to find partition
			start_time = clock();
			K->partition(print);
			double time_partition = ((double)(clock() - start_time) / CLK_TCK);

			LotteryDesigner* L = new LotteryDesigner(K, print);
			std::vector<lottery> L_vector = L->compare_methods(s, iterations, false);
				// The first element of 'L_vector' will contain the lottery we asked for
				
			////////////////////////////////////////////////////////////
			// Find a way to include time limits into this!
			////////////////////////////////////////////////////////////

			letter;

			for (int k = 0; k < s.size(); k++) {
				letter = s[k];

				// Decide the correct name of the file
				if (s[k] == 'C') {
					filename_distr = "_NASH";
				}
				if (s[k] == 'R') {
					filename_distr = "_RSD";
				}
				if (s[k] == 'L') {
					filename_distr = "_LEXIMAX";
				}
				if (s[k] == 'U') {
					filename_distr = "_UNIFORM";
				}
				if (s[k] == 'O') {
					filename_distr = "_RE-INDEX";
				}
				if (s[k] == 'P') {
					filename_distr = "PERTURB";
				}
				if (s[k] == 'T') {
					filename_distr = "RSD-ONCE";
				}

				filename = filename_front + filename_distr + filename_back;
				std::string filename_long = filename_front + filename_distr + "_distr" + filename_back;

				ReportDistribution* empty = new ReportDistribution(name, K, L_vector[k], print);
				empty->time_single_solution = time_solve;
				empty->time_partition = time_partition;
				empty->time_distribution = L_vector[k].t;

				std::ofstream O;
				O.open(filename, std::fstream::app);
				O << name << ", " << P[p][0] << ", ";
				O << K->M_size << ", ";
				O << empty->time_single_solution << ",";
				O << empty->time_partition << ",";
				O << empty->time_distribution << ",";
				O << empty->K->solver_times << ",";

				// Count number of solutions with strictly positive weight
				int count = 0;
				for (int t = 0; t < L_vector[k].S.size(); t++) {
					if (L_vector[k].w[t] > 0.0001) {
						count++;
					}
				}
				O << count << ",";

				// Compute some metrics
				double Nash_product = 1.0;
				double mean = 0.0;
				double min_selection_prob = 1.0;
				for (int t = 0; t < K->I.n; t++) {
					if (K->M[t] == 1) {
						Nash_product = Nash_product * L_vector[k].p[t];
						mean += L_vector[k].p[t];
						if (L_vector[k].p[t] < min_selection_prob) {
							min_selection_prob = L_vector[k].p[t];
						}
					}
				}
				if (K->M_size > 0) {
					Nash_product = pow(Nash_product, ((double)1 / (double) K->M_size));
					mean = mean / (double)K->M_size;
				}
				else {
					Nash_product = 0.0;
					mean = 0.0;
				}
				O << Nash_product << "," << mean << "," << min_selection_prob << ",";
				O << "\n";

				R.push_back(empty);


				O.close();
				O.clear();

				Q.open(filename_long, std::fstream::app);
				Q << name << ", " << P[p][0] << ", ";
				Q << K->M_size << ", ";
				for (int t = 0; t < K->I.n; t++) {
					Q << L_vector[k].p[t] << ",";
				}
				Q << "\n";

				Q.close();
				Q.clear();

				delete empty;
			}
			delete L;
			delete K;
			delete KE;
			L_vector.clear();
		}
	}

	return R;
}

