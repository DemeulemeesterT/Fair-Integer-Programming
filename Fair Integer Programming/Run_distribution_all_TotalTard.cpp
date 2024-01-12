#include"Run_configurations.h"

std::vector<ReportDistribution*> run_distribution_all_TotalTard(std::string s, int inst_generated, unsigned seed, double time_limit, int iterations, bool print) {
	
	
	// Prepare to go through all instances
	std::vector<int> N_vector;
	std::vector<double> beta_vector;
	N_vector = { 30 };
	beta_vector = { 0.5 };

	// Create vector with seeds, and export it
	std::mt19937 generator((unsigned)seed);
	std::uniform_int_distribution<int> distr(-100000000, 100000000);
	std::vector<std::vector<std::vector<int>>> seed_vector(N_vector.size(), std::vector<std::vector<int>>(beta_vector.size(), std::vector<int>(inst_generated, 0)));

	for (int p = 0; p < N_vector.size(); p++) {
		for (int b = 0; b < beta_vector.size(); b++) {
			for (int i = 0; i < inst_generated; i++) {
				seed_vector[p][b][i] = distr(generator);
			}
		}
	}


	// Change this conditional on the name of the distribution you're looking for.
	char filename_front[50];
	char filename_back[50];
	std::string filename_distr = "_";
	std::string filename;
	sprintf_s(filename_front, "Generated Results/ReportTT");
	sprintf_s(filename_back, ".csv");
	std::string letter;
	std::ofstream O;
	std::ofstream Q;
	std::ofstream XMIN;
	std::ofstream XMAX;

	// Export seeds
	char numstr1[21];
	char numstr2[21];
	char numstr3[21];
	std::string name;
	std::string name_instance = "TT_";
	std::string name_conn = "_";
	std::ofstream SEEDS;
	std::string seeds_name = "_seeds";
	std::string seeds_extension = ".txt";
	std::string filename_seeds = filename_front + seeds_name + seeds_extension;
	SEEDS.open(filename_seeds);
	for (int p = 0; p < N_vector.size(); p++) {
		for (int b = 0; b < beta_vector.size(); b++) {
			for (int i = 0; i < inst_generated; i++) {
				sprintf_s(numstr1, "%i", N_vector[p]);
				sprintf_s(numstr2, "%.2f", beta_vector[b]);
				sprintf_s(numstr3, "%i", i);
				name = name_instance + numstr1 + name_conn + numstr2 + name_conn + numstr3;
				SEEDS << name << " " << seed_vector[p][b][i] << "\n";
			}
		}
	}
	SEEDS.close();

	for (int k = 0; k < s.size(); k++) {
		letter = s[k];
		// Decide the correct name of the file
		if (s[k] == 'C') {
			filename_distr = "_NASH";
		}
		if (s[k] == 'R') {
			filename_distr = "_RSD";
		}
		if (s[k] == 'H') {
			filename_distr = "_RSDheuristic";
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
		O << "Name, n, |M|, time solution, time partition, time distribution, generated solutions, used solutions, geometric mean (rel.), arithmetic mean (rel.), minimum selection probability (rel.), geometric mean (abs.), arithmetic mean (abs.), minimum selection probability (abs.)" << "\n";
		O.close();
		// 'rel.' refers to (distr_i - xmin_i)/(xmax_i - xmin_i)
		// 'abs.' refers to (distr_i - xmin_i)

	// Now create a separate file that contains the distributions themselves
		std::string filename_long = filename_front + filename_distr + "_distr" + filename_back;
		Q.open(filename_long);
		Q << "Name, n, |M|, distribution" << "\n";
		Q.close();

		// Store Xmin and Xmax
		std::string filename_XMIN = filename_front + filename_distr + "_XMIN" + filename_back;
		XMIN.open(filename_XMIN);
		XMIN << "Name, n, |M|, XMIN" << "\n";
		XMIN.close();
		std::string filename_XMAX = filename_front + filename_distr + "_XMAX" + filename_back;
		XMAX.open(filename_XMAX);
		XMAX << "Name, n, |M|, XMAX" << "\n";
		XMAX.close();
	}

	std::vector<ReportDistribution*> R;

	for (int p = 0; p < N_vector.size(); p++) {
		for (int b = 0; b < beta_vector.size(); b++) {
			for (int i = 0; i < inst_generated; i++) {
				sprintf_s(numstr1, "%i", N_vector[p]);
				sprintf_s(numstr2, "%.2f", beta_vector[b]);
				sprintf_s(numstr3, "%i", i);
				name = name_instance + numstr1 + name_conn + numstr2 + name_conn + numstr3;

				SchedWT_param S;
				S.n_jobs = N_vector[p];
				S.seed = seed_vector[p][b][i];
				S.name = name;
				bool release_dates = false;
				double beta = beta_vector[b];

				SchedulingWeightTard* SWT = new SchedulingWeightTard(false);
				inst I = SWT->generate_data_and_instance_TIF_WT(S, beta, release_dates, true, false);

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
					if (s[k] == 'H') {
						filename_distr = "_RSDheuristic";
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
					std::string filename_XMIN = filename_front + filename_distr + "_XMIN" + filename_back;
					std::string filename_XMAX = filename_front + filename_distr + "_XMAX" + filename_back;

					ReportDistribution* empty = new ReportDistribution(name, K, L_vector[k], print);
					empty->time_single_solution = time_solve;
					empty->time_partition = time_partition;
					empty->time_distribution = L_vector[k].t;

					std::ofstream O;
					O.open(filename, std::fstream::app);
					O << name << ", " << N_vector[p] << ", ";
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
					double Nash_product_rel = 1.0;
					double Nash_product_abs = 1.0;
					double mean_rel = 0.0;
					double mean_abs = 0.0;
					double min_selection_prob_rel = 1.0;
					double min_selection_prob_abs = 1.0;

					// For relative
					for (int t = 0; t < K->I.n; t++) {
						if (K->M[t] == 1) {
							Nash_product_rel = Nash_product_rel * (L_vector[k].p[t] - L->K->Xmin[t]) / (L->K->Xmax[t] - L->K->Xmin[t]);
							mean_rel += (L_vector[k].p[t] - L->K->Xmin[t]) / (L->K->Xmax[t] - L->K->Xmin[t]);
							if ((L_vector[k].p[t] - L->K->Xmin[t]) / (L->K->Xmax[t] - L->K->Xmin[t]) < min_selection_prob_rel) {
								min_selection_prob_rel = (L_vector[k].p[t] - L->K->Xmin[t]) / (L->K->Xmax[t] - L->K->Xmin[t]);
							}
						}
					}

					if (K->M_size > 0) {
						Nash_product_rel = pow(Nash_product_rel, ((double)1 / (double)K->M_size));
						mean_rel = mean_rel / (double)K->M_size;
					}
					else {
						Nash_product_rel = 0.0;
						mean_rel = 0.0;
					}

					// For absolute
					for (int t = 0; t < K->I.n; t++) {
						if (K->M[t] == 1) {
							Nash_product_abs = Nash_product_abs * (L_vector[k].p[t] - L->K->Xmin[t]);
							mean_abs += (L_vector[k].p[t] - L->K->Xmin[t]);
							if ((L_vector[k].p[t] - L->K->Xmin[t]) < min_selection_prob_abs) {
								min_selection_prob_abs = (L_vector[k].p[t] - L->K->Xmin[t]);
							}
						}
					}

					if (K->M_size > 0) {
						Nash_product_abs = pow(Nash_product_abs, ((double)1 / (double)K->M_size));
						mean_abs = mean_abs / (double)K->M_size;
					}
					else {
						Nash_product_abs = 0.0;
						mean_abs = 0.0;
					}


					O << Nash_product_rel << "," << mean_rel << "," << min_selection_prob_rel << ",";
					O << Nash_product_abs << "," << mean_abs << "," << min_selection_prob_abs << ",";
					O << "\n";

					R.push_back(empty);


					O.close();
					O.clear();

					Q.open(filename_long, std::fstream::app);
					Q << name << ", " << N_vector[p] << ", ";
					Q << K->M_size << ", ";
					for (int t = 0; t < K->I.n; t++) {
						Q << L_vector[k].p[t] << ",";
					}
					Q << "\n";

					Q.close();
					Q.clear();

					XMIN.open(filename_XMIN, std::fstream::app);
					XMIN << name << ", " << N_vector[p] << ", ";
					XMIN << K->M_size << ", ";
					for (int t = 0; t < K->I.n; t++) {
						XMIN << L->K->Xmin[t] << ",";
					}
					XMIN << "\n";

					XMIN.close();
					XMIN.clear();

					XMAX.open(filename_XMAX, std::fstream::app);
					XMAX << name << ", " << N_vector[p] << ", ";
					XMAX << K->M_size << ", ";
					for (int t = 0; t < K->I.n; t++) {
						XMAX << L->K->Xmax[t] << ",";
					}
					XMAX << "\n";

					XMAX.close();
					XMAX.clear();

					delete empty;
				}
				delete L;
				delete K;
				delete SWT;
				L_vector.clear();
			}
		}
	}

	return R;
}

