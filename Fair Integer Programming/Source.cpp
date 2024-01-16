#include "Run_configurations.h"

using namespace std;

int main()
{
	//run_distribution_all_TotalTard("C", 50, 0, 3600.0, 1000, false);
	//run_distribution_all_KE("L", 3600.0, 1, false);
	//run_distribution_all_SWT("H", 50, 3600.0, 1000, false);
	//compare_time_normal_RSD_SWT(1000, 50, false);
	
	
	//SchedulingWeightTard* SWT = new SchedulingWeightTard("wt40", false);
	//inst I = SWT->generate_instanceTIF_WT(1, false, true);
	SchedWT_param S;
	//S.common_process_time = 1;
	S.n_jobs = 10;
	S.seed = 80757630;
	S.name = "test_TT_10_0.25_16";
	bool release_dates = false;
	double beta = 0.25;

	SchedulingWeightTard* SWT = new SchedulingWeightTard(false);
	inst I = SWT->generate_data_and_instance_TIF_WT(S, beta, release_dates, true, true);
	IPSolver* K = new IPSolver(I, true);
	K->solve(true);
	K->model->write("Generated Formulations/IPModel.lp");
	K->analyze(true);
	LotteryDesigner* L = new LotteryDesigner(K, true);
	L->compare_methods("CLR", 100, true, false, 0);
	delete K;
	delete SWT;
	
	
	
	
	//compare_time_normal_RSD(50, false);

	//std::vector<reportGreedy> RG;
	//RG = run_Kidney_exchange_all(true);


	// Uncapacitated Facility Location
	/*Ufl_parameters U;
	U.m = 40;
	U.n = 80;
	U.f_min = 1000;
	U.f_max = 1000;
	U.c_min = 100;
	U.c_max = 1000;
	U.nr_indiff_classes = -1;
	U.filename = "Test";
	U.seed = 125382;

	Ufl myUFL(U, false);*/
	//Ufl myUFL("cap134", true);
	//inst I = myUFL.generate_formulation(false, true);
	
	
	// INSTANCE TO TEST CARDINAL MODEL
	/*inst I;
	I.m = 1;
	I.n = 3;
	I.t = 0;
	I.v = { 4,2,1 };
	I.C = { 8 };
	I.A = { { 4,2,1 } };
	I.n_var = 3;
	I.Y_bool = false;
	I.Y_coeff_zero = true;
	I.Y_integer = false;
	I.X_bool = false;
	I.X_integer = true;
	*/
	
	
	/*KidneyExchange* KE = new KidneyExchange("30-instance-32", true);
	//KidneyExchange* KE = new KidneyExchange("70-instance-3", true);

	inst I = KE->generate_instance(true, false);
	

	IPSolver* K = new IPSolver(I, true);
	//K->model->write("Generated Formulations/IPModel.lp");
	//K->model->set(GRB_IntParam_PoolSearchMode, 2);
	//K->model->set(GRB_IntParam_PoolSolutions, 10);
	//K->model->set(GRB_IntParam_OutputFlag, 1);
	//K->model->optimize();
	
	//K->compare_time_normal_vs_RSD_variants(50, true);
	

	//const std::string file_name = "MIPLIB instances/sorrell8.mps";
	//IPSolver* K = new IPSolver(file_name, true);

	//IPSolver* K = new IPSolver(true);
	//K->partition(true);
	K->analyze(false);


	//K->compare_partition_vs_greedy(true);
	//K->compare_time_normal_vs_RSD_without_partition(5, true);
	LotteryDesigner* L = new LotteryDesigner(K, false);

	std::vector<lottery> results = L->compare_methods("LR", 100, true, false, 0);
	
	//SimplicalDecomposition* SD = new SimplicalDecomposition(K, true);
	//SD->SD_Nash(true);

	//delete SD;
	//delete K;
	//delete KE;
	*/
}