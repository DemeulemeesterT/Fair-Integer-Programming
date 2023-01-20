#include "Run_configurations.h"

using namespace std;

int main()
{
	
	//compare_time_normal_RSD(50, false);

	//std::vector<reportGreedy> RG;
	//RG = run_Kidney_exchange_all(true);

	KidneyExchange* KE = new KidneyExchange("70-instance-19", true);
	//KidneyExchange* KE = new KidneyExchange("20-instance-1", true);

	inst I = KE->generate_instance(true, true);
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
	L->compare_methods("N", 200, true, 0);
	
	//SimplicalDecomposition* SD = new SimplicalDecomposition(K, true);
	//SD->SD_Nash(true);

	//delete SD;
	delete K;
	//delete KE;


}