#include"Run_configurations.h"

void compare_time_normal_RSD(int iter, bool print) {
	// Prepare to go through all instances of kidney exchange
	std::vector<std::vector<int>> P;
	P.push_back({ 10,50 }); // The first number is the number of agents, the second the ID of the instance
	P.push_back({ 20,50 });
	P.push_back({ 30,50 });
	P.push_back({ 40,65 });
	P.push_back({ 50,50 });
	P.push_back({ 60,50 });
	P.push_back({ 70,50 });

	std::ofstream O;
	O.open("Generated Results/RSD_vs_Normal_KE.csv");
	O << "Name; Agents; t_Normal; t_RSD_no_partition; t_RSD; t_partition" << "\n";
	O.close();

	char numstr1[21];
	char numstr2[21];
	std::string name;
	std::string name_instance = "-instance-";
	for (int p = 0; p < P.size(); p++) {
		for (int i = 1; i <= P[p][1]; i++) {
			reportGreedy empty;
			sprintf_s(numstr1, "%i", P[p][0]);
			sprintf_s(numstr2, "%i", i);
			name = numstr1 + name_instance + numstr2;
			KidneyExchange* KE = new KidneyExchange(name, false);
			inst I = KE->generate_instance(false, false);
			IPSolver* K = new IPSolver(I, false);
			std::vector<time_report> TR = K->compare_time_normal_vs_RSD_variants(iter, false);

			std::ofstream O;
			O.open("Generated Results/RSD_vs_Normal_KE.csv", std::fstream::app);
			O << name << "; " << K->I.n << "; " << TR[0].t << "; " << TR[2].t << "; " << TR[1].t << "; " << K->time_partition << "\n";
			
			O.close();

			delete K;
			delete KE;
		}
	}

}
