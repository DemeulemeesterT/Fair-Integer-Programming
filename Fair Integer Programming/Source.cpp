#include "Run_configurations.h"

using namespace std;

int main()
{
	//run_distribution_all_KE("H", 3600.0, 1, false);
	//run_distribution_all_SWT("H", 50, 3600.0, 1000, false);
	//compare_time_normal_RSD_SWT(1000, 50, false);
	
	/*
	SchedulingWeightTard* SWT = new SchedulingWeightTard("wt40", false);
	inst I = SWT->generate_instanceLawlerMoore(0, false, true);
	IPSolver* K = new IPSolver(I, true);
	K->solve(true);
	LotteryDesigner* L = new LotteryDesigner(K, true);
	L->compare_methods("RLCU", 1000, true, false, 0);
	delete K;
	delete SWT;
	*/
	

	
	//compare_time_normal_RSD(50, false);

	//std::vector<reportGreedy> RG;
	//RG = run_Kidney_exchange_all(true);

	
	KidneyExchange* KE = new KidneyExchange("10-instance-2", true);
	//KidneyExchange* KE = new KidneyExchange("70-instance-3", true);

	inst I = KE->generate_instance(true, true);
	I.X_bool = false;
	IPSolver* K = new IPSolver(I, true);
	//K->compare_time_normal_vs_RSD_variants(50, true);
	

	//const std::string file_name = "MIPLIB instances/sorrell8.mps";
	//IPSolver* K = new IPSolver(file_name, true);

	//IPSolver* K = new IPSolver(true);
	//K->partition(true);
	K->analyze(true);
	//K->compare_partition_vs_greedy(true);
	//K->compare_time_normal_vs_RSD_without_partition(5, true);

	LotteryDesigner* L = new LotteryDesigner(K, true);
	L->compare_methods("HCR", 10, true, false, 0);
	
	//SimplicalDecomposition* SD = new SimplicalDecomposition(K, true);
	//SD->SD_Nash(true);

	//delete SD;
	delete K;
	//delete KE;
	
}