#include "IPSolver.h"
#include "KidneyExchange.h"
#include "SchedulingWeightTard.h"
#include "LotteryDesigner.h"
#include "SimplicalDecomposition.h"
#include "ReportDistribution.h"

// Contains a couple of different runs, like this we don't have to write everything in 'source.cpp'

#pragma once

struct reportGreedy {
	std::string name;
	int n_wrongly_class;
	int n_M;
	std::vector<int> agents_wrong;
	int n_S;
	int n_S_greedy;
	double time_partition;
	double time_greedy_partition;
};

std::vector<reportGreedy> run_Kidney_exchange_all(bool print);

std::vector<ReportDistribution*> run_distribution_all_KE(std::string s, double time_limit, int iterations, bool print);
	// Can be used to find distributions for multiple Kidney Exchange instances sequentially
	// Note 22/02/2023: the time_limit parameter is not implemented yet
		// The choice of the instances to be run can be controlled at the beginning of 'Run_distribution_all_KE.cpp'

std::vector<ReportDistribution*> run_distribution_all_SWT(std::string s, int nr_instances, double time_limit, int iterations, bool print);
	// Can be used to find distributions for multiple Scheduling Weighted Tardiness instances sequentially
	// Note 22/02/2023: the time_limit parameter is not implemented yet
		// The choice of the instances to be run can be controlled at the beginning of 'Run_distribution_all_SWT.cpp'


void compare_time_normal_RSD(int iter, bool print);
	// Compare time solving one instance, solving RSD, and solving RSDheuristic for KIDNEY EXCHANGE

void compare_time_normal_RSD_SWT(int iter, bool print);
	// Compare time solving one instance, solving RSD, and solving RSDheuristic for MINIMIZING THE WEIGHTED NUMBER OF TARDY JOBS




