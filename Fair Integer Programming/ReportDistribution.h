#include"IPSolver.h"

#pragma once
class ReportDistribution
{
public:
	std::string name;
	IPSolver* K;

	// Computation times
	double time_partition;
	double time_distribution;
	double time_single_solution;

	// Solutions and lottery
	lottery L;

	ReportDistribution(IPSolver* K_in, bool print);
	~ReportDistribution();
};

