#include "ReportDistribution.h"

ReportDistribution::ReportDistribution(IPSolver* K_in, bool print) {
	// Initialize the variables
	name = "";
	K = new IPSolver(K_in, print);
	time_partition = -2.0;
	time_distribution = -2.0;
	time_single_solution = -2.0;

}

ReportDistribution::~ReportDistribution() {
	delete K;
}
