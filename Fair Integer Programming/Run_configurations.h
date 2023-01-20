#include "IPSolver.h"
#include "KidneyExchange.h"
#include "LotteryDesigner.h"
#include "SimplicalDecomposition.h"

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

void compare_time_normal_RSD(int iter, bool print);
