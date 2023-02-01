#include "ReportDistribution.h"

ReportDistribution::ReportDistribution(std::string name_in, IPSolver* K_in, lottery L_in, bool print) {
	// Initialize the variables
	name = name_in;
	K = new IPSolver(K_in, print);
	
	L = lottery();
	L.p = L_in.p;
	L.w = L_in.w;
	L.S = L_in.S;

	time_partition = -2.0;
	time_distribution = -2.0;
	time_single_solution = -2.0;
}

ReportDistribution::~ReportDistribution() {
	delete K;
}
