#include "IPSolver.h"

#pragma once

struct SchedWT_inst {
	std::vector<int> w;
	// The weights of the jobs in the objective function

	std::vector<int> d;
	// The due dates of the jobs

	std::vector<int> p;
	// The processing times of the jobs

	std::vector<int> r;
	// Release dates of the jobs

};

struct SchedWT_param {
	// Parameters to generate a WT instance using settings from paper by van den Akker et al. (2010)
	int n_jobs;
	int common_process_time;
	unsigned seed = 0;
	std::string name;
};

class SchedulingWeightTard
{
public:
	// FUNCTIONS
	void read_data(std::string name, bool print);
		// Will fill the 'I'-vector, which contains the instances

// FUNCTIONS FOR MINIMIZING THE NUMBER OF TARDY JOBS //
	inst generate_instanceTIF(int nr, bool export_inst, bool print);
		// Will create an instance object for the i-th instance
		// The formulation that it will describe is a time-indexed formulation
			// This formulation is rather slow to solve...

	inst generate_instanceLawlerMoore(int nr, bool export_inst, bool print);
		// Will create an instance for the i-th instance
		// The formulation is based on Lawler & Moore (1969) - A Functional Equation and its Application to Resource Allocation and Sequencing Problems
	// CAREFUL! THIS FORMULATION USES A DOMINANCE RULE TO FIND AN OPTIMAL SOLUTION
		// This means that this formulation cannot find all of the optimal solutions
		// As a result, our methods will find a less fair distribution than would be possible when
			// optimizing over all of the optimal solutions.


	SchedWT_inst order_jobs(int nr, bool print);
		// Orders the jobs in the i-th instance according to earliest due date


// FUNCTIONS FOR MINIMIZING THE WEIGHTED TARDINESS
	inst generate_instanceTIF_WT(int nr, bool export_inst, bool print);
	// Will create an instance object for the i-th instance
	// The formulation that it will describe is a time-indexed formulation
		// Similar to 'generate_instancesTIF()', but binary x_j variables are replaced by integer variables
	inst generate_instanceTIF_WT(SchedWT_inst SI, bool export_inst, bool print);

	inst generate_data_and_instance_TIF_WT(SchedWT_param param, bool export_inst, bool print);
	// Generate data as in paper van den Akker et al. (2010), Section 9

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

	SchedulingWeightTard(bool print);
		// Can be used to just create object and generate data later.

	~SchedulingWeightTard();

};

