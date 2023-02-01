#include"Run_configurations.h"

std::vector<ReportDistribution> run_distribution_all_KE(double time_limit, int iterations, bool print) {
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
	std::ofstream O;
	O.open("Generated Results/ReportDistributionsKE.csv");
	O << "Name, n, |M|, time solution, time partition, time distribution" << "\n";
	O.close();

	std::vector<ReportDistribution> R;

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
			std::vector<lottery> L_vector = L->compare_methods("N", iterations);
				// The first element of 'L_vector' will contain the lottery we asked for
				
			////////////////////////////////////////////////////////////
			// Find a way to include time limits into this!
			////////////////////////////////////////////////////////////

			ReportDistribution empty(name, K, L_vector[0], print);
			empty.time_single_solution = time_solve;
			empty.time_partition = time_partition;
			empty.time_distribution = L_vector[0].t;
			R.push_back(empty);
			
			std::ofstream O;
			O.open("Generated Results/ReportDistributionsKE.csv", std::fstream::app);
			O << name << ", " << P[p][0] << ", ";
			O << K->M_size << ", ";
			O << empty.time_single_solution << ",";
			O << empty.time_partition << ",";
			O << empty.time_distribution << ",";
			O << "\n";

			O.close();
			delete K;
			delete L;
			delete KE;
		}
	}

	return R;
}

