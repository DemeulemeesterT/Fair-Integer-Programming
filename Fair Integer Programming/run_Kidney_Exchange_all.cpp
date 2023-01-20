#include "Run_configurations.h"

std::vector<reportGreedy> run_Kidney_exchange_all(bool print) {

	// Prepare to go through all instances
	std::vector<std::vector<int>> P;
	P.push_back({ 10,50 }); // The first number is the number of agents, the second the ID of the instance
	P.push_back({ 20,50 });
	P.push_back({ 30,50 });
	P.push_back({ 40,65 });
	P.push_back({ 50,50 });
	P.push_back({ 60,50 });
	P.push_back({ 70,50 });

	std::ofstream O;
	O.open("Generated Results/Greedy_KE.csv");
	O << "Name, Wrongly classified, |M|, time partition, time greedy partition, Agents" << "\n";
	O.close();

	std::vector<reportGreedy> RG;
	

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
			K->compare_partition_vs_greedy(false);

			std::ofstream O;
			O.open("Generated Results/Greedy_KE.csv", std::fstream::app);
			O << name << ", " << K->greedy_wrongly_classified.size() << ", ";
			O << K->M_size << ", ";
			O << K->time_partition << ", ";
			O << K->time_greedy_partition;
			for (int j = 0; j < K->greedy_wrongly_classified.size(); j++) {
				O << ", " << K->greedy_wrongly_classified[j];
			}
			O << "\n";

			empty.name = name;
			empty.n_wrongly_class = K->greedy_wrongly_classified.size();
			empty.n_M = K->M_size;
			empty.agents_wrong = K->greedy_wrongly_classified;
			empty.n_S = K->S.size();
			empty.n_S_greedy = K->S_greedy.size();
			empty.time_partition = K->time_partition;
			empty.time_greedy_partition = K->time_greedy_partition;
			RG.push_back(empty);

			O.close();
			delete K;
			delete KE;
		}
	}
	
	return RG;
}