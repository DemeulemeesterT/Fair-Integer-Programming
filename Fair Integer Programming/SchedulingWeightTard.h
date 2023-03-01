#include "IPSolver.h"

#pragma once

struct SchedWT_inst {
	std::vector<int> w;
	// The weights of the jobs in the objective function

	std::vector<int> d;
	// The due dates of the jobs

	std::vector<int> p;
	// The processing times of the jobs

};


class SchedulingWeightTard
{
public:
	// FUNCTIONS
	void read_data(std::string name, bool print);
		// Will fill the 'I'-vector, which contains the instances

	inst generate_instanceTIF(int nr, bool export_inst, bool print);
		// Will create an instance object for the i-th instance
		// The formulation that it will describe is a time-indexed formulation
			// This formulation is rather slow to solve...

	inst generate_instanceLawlerMoore(int nr, bool export_inst, bool print);
		// Will create an instance for the i-th instance
		// The formulation is based on Lawler & Moore (1969) - A Functional Equation and its Application to Resource Allocation and Sequencing Problems

	SchedWT_inst order_jobs(int nr, bool print);
		// Orders the jobs in the i-th instance according to earliest due date


	// GENERAL VARIABLES
	std::vector<SchedWT_inst> Sched_I;
		// Each element of this vector contains the data for an instance

	std::vector<SchedWT_inst> Sched_I_ordered;
		// Orders the jobs according to earliest due date


	int n; 
		// The number of jobs

	std::string data_name;
		// Name of the instance


	// CONSTRUCTORS AND DESTRUCTORS
	SchedulingWeightTard(std::string name, bool print);
		// This constructor will read file 'name.txt'
		// Instances are obtained from http://people.brunel.ac.uk/~mastjjb/jeb/orlib/wtinfo.html 

	~SchedulingWeightTard();

};

