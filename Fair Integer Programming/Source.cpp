#include "Run_configurations.h"

using namespace std;

int main()
{
	KidneyExchange* KE = new KidneyExchange("40-instance-48", true);

	inst I = KE->generate_instance(true, false);

	IPSolver* K = new IPSolver(I, true);

	K->analyze(false);

	LotteryDesigner* L = new LotteryDesigner(K, false);

	std::vector<lottery> results = L->compare_methods("ULRC", 1000, true, false, 0);

	delete K;


	// CODE TO RUN INSTANCES IN BATCHES & EXPORT RESULTS:
	//run_distribution_all_TotalTard("RO", 50, 0, 3600.0, 1000, false);
	//run_distribution_all_KE("L", 3600.0, 1, false);
	//run_distribution_all_SWT("H", 50, 3600.0, 1000, false);
	//compare_time_normal_RSD_SWT(1000, 50, false);
	
	
	// CODE TO RUN WEIGHTED TARDINESS EXAMPLE
	//SchedulingWeightTard* SWT = new SchedulingWeightTard("wt40", false);
	//inst I = SWT->generate_instanceTIF_WT(1, false, true);
	/*SchedWT_param S;
	//S.common_process_time = 1;
	S.n_jobs = 15;
	S.seed = -5658298;
	S.name = "test";
	bool release_dates = false;
	double beta = 0.25;

	SchedulingWeightTard* SWT = new SchedulingWeightTard(false);
	inst I = SWT->generate_data_and_instance_TIF_WT(S, beta, release_dates, true, true);
	IPSolver* K = new IPSolver(I, true);
	K->solve(true);
	K->model->write("Generated Formulations/IPModel.lp");
	K->analyze(true);
	LotteryDesigner* L = new LotteryDesigner(K, true);
	L->compare_methods("R", 1000, true, false, 0);
	delete K;
	delete SWT;
	*/
	
}