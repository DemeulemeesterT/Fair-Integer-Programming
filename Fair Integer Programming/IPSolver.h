#pragma once

#include <iostream>
#include <string>
#include <random>
#include <time.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "gurobi_c++.h"
#include <cmath>
#include "PartitionCallback.h"

// An IP instance
struct inst {
	int m; // Number of knapsacks/constraints
	int n; // Number of agents/X-variables: binary variables that represent agents who we want to treat fairly in the selection of an optimal solution
	int t; // Number of additional Y-variables: auxiliary binary/integer variables 
	std::vector<double> v; // Objective coefficients, first for the X-variables, then for the Y-variables
	std::vector<double> C; // Capacities of the knapsacks
	std::vector<std::vector<double>> A; // m x (n + t) matrix of constraints (again, first for X-variables, then for Y-variables)
	std::string data_name;
	int n_var; // The total number of variables = n + t
	bool Y_bool; // 'true' if Y-variables are binary, 'false' if Y-variables are integer 
	bool Y_coeff_zero; // 'true' if Y-variables all have coefficient zero in the objective function, 'false' otherwise.
	bool Y_integer; // 'true' if Y-variables are integer;
	bool X_bool; // 'true' if X-variables are binary (dichotomous preferences), 'false' if X-variables are integer or real (cardinal preferences)
	bool X_integer; // 'true' if X-variables are integer, 'false' if X-variables are real values.
};

// One optimal solution to an instance
struct solution {
	std::vector<double> x;
	std::vector<double> y;
	int ID; // This is a unique identifier of the solution. Is equal to -1 if ID is too large to store.
};

// This structure will represent a lottery over the optimal solutions
struct lottery {
	std::vector<solution> S; // The set of optimal solutions
	std::vector<double> w; // The weights of the optimal solutions in the lottery
	std::vector<double> p; // The selection probabilities of the agents
		// NOTE: this code was first written for binary X-variables (dichotomous preferences).
			// That's why we use the term 'probabilities' to refer to the expected value of the decision variable
			// in the final distribution.
			// In general, 'p' refers to the vector containing the expected utilities of the agents
	double t; // The solution time needed to find this lottery.
};

// This structure will return the time, and a label
struct time_report {
	std::string label;
	double t;
};

struct MIPLIB_instance {
	std::string name; // Name of the instance
	std::string solution_info; // =opt= for optimal solution, =best= for best-known solution, =inf= for infeasible, =unkn= for unknown
	double obj_val;
};

struct IP_report {
	double opt_obj_value;
	solution s;
	bool solution_already_in_S;
		// Is 'true' if found solution was already in S.
};

class IPSolver
{
public:
// GUROBI VARIABLES
	GRBEnv* env;
	GRBModel* model;
	GRBVar* X;
	GRBVar* Y_var;
	
	GRBEnv* env_identical;
		// Environment used to find identical agents
		// Created here to avoid message 'Academic license - ...' being printed every time.


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// BINARY X-VARIABLES & GENERAL ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

// GENERAL VARIABLES
	int opt;
	
	int solver_times;
		// Counts the number of times the solver was called

	int number_of_agents_at_once_RSD;
	
	inst I;
		// Sets that contain the agents that are always, sometimes, never selected in a solution
			// -1 means it is unknown
			// 0 means it is not an element
			// 1 means it is an element
	
	std::vector<int> Y;
	std::vector<int> M;
	std::vector<int> N;

	int M_size;
		// The number of agents in M.
	
	std::vector<std::vector<int>> identical;
		// Binary matrix for the identical agents
		// This is an upper triangular matrix, all elements in the lower triangle are zero.
	
	std::vector<int> identical_clustered;
		// Contains cluster id for each agents, -1 if not in a cluster
	
	int cluster_count;
	
	std::vector<solution> S; 
		// The set of all solutions that we have found so far

	std::vector<solution> S_greedy;
		// The set of all solutions found when doing the greedy partition like in
			// Farnadi, G., St-Arnaud, W., Babaki, B., & Carvalho, M. (2021, May). Individual Fairness in Kidney Exchange Programs. 
				// In Proceedings of the AAAI Conference on Artificial Intelligence (Vol. 35, No. 13, pp. 11496-11505).

	std::vector<int> greedy_wrongly_classified;

	double time_partition;
	double time_greedy_partition;

	bool done_partition;
	// 'true' if partition is done
	bool done_greedy_partition;
	bool done_find_identical;
	bool done_solve;

	std::vector<MIPLIB_instance> MIPLIB_summary;
		// Contains information about the MIPLIB instances (name, optimal or best-known solutions...)

// GUROBI constraints
	std::vector<GRBConstr> RSD_fixed;
		// This vector will contain the constraints in which a variable is fixed to the found value in RSD
		


// FUNCTIONS
	double solve(bool print);
		// Solves the IP instance to optimality

	double solve_partition(bool print);
		// Same as 'solve', but used in the partitioning step. 
		// Contains an additional callback that causes the solver the stop once the upper bound 
		// on the objective function is lower than the optimal objective function.

	IP_report solve_return_solution(bool print);
		// Same as 'solve' but returns optimal solution and the optimal objective value

	void analyze(bool print);
		// Calls 'solve', 'partition' and 'find_identical agents'
	
	void partition(bool print);
		// Find sets Y, M and N
			// If the X-variables are not binary, this will find the lowest and the highest value
			// that each of the agents experiences in any of the optimal solutions.

	void greedy_partition(bool print);
		// Find greedy partition based on algorithm by 
			// Farnadi, G., St-Arnaud, W., Babaki, B., & Carvalho, M. (2021, May). Individual Fairness in Kidney Exchange Programs. 
				// In Proceedings of the AAAI Conference on Artificial Intelligence (Vol. 35, No. 13, pp. 11496-11505).
	
	void compare_partition_vs_greedy(bool print);

	void find_identical_agents(bool print);
		// Find identical agents

	solution RSD_once(std::vector<int> order, bool print);
		// Takes as an input an ordering of the students
		// and will perturb the objective function accordingly so the agents get selected in this order (among the optimal solutions)
		// NOTE: 'order' does not have to contain all agents necessarily;
			// the agents that are not in 'order' will not have their objective coefficients permuted.
			// Only the agents in M who appear in order will have their objective coefficients permuted

	solution RSD_once_no_partition(std::vector<int> order, bool print);
		// Same as 'RSD_once', but does not take into account the partition of the agents into Y, N, and M.
		// Objective coefficients of the agents in Y or N can also be permuted in this function.

	solution RSD_heuristic_no_partition(std::vector<int> order, bool print);
	// Same as 'RSD_once_no_partition', but only perturbs the objective function once, for the first 19 agents in the random order.
	
	void block_solution(solution sol);

	std::vector<time_report> compare_time_normal_vs_RSD_variants(int iterations, bool print);
	std::vector<time_report> compare_time_normal_vs_RSD_heuristic(int iterations, bool print);
	std::vector<time_report> compare_time_normal_vs_RSD_without_partition(int iterations, bool print);

	// (DE)CONSTRUCTORS

	// There are four constructors in this class:
		// The first one will prompt to provide the name of the instance you want to work on
			// or will generate one randomly based on your input
		// The second one will create an IPSolver object based on an instance
		// The third one will create an IPSolver object based on another IPSolver object 'K_in' (used to clone a class instance)
		// The fourth one will create an IPSolver object based on an .mps file that contains a MIP formulation
	IPSolver(bool print);
	IPSolver(inst I_in, bool print);
	IPSolver(IPSolver* K_in, bool print);
	IPSolver(const std::string& file_name, bool print);

	~IPSolver();


//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// NON-BINARY X-VARIABLES /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

// VARIABLES
	std::vector<double> Xmin, Xmax;
		// Minimum and maximum value each of the X-variables obtain in any of the optimal solutions
		// These two vectors can contain all information that was stored in N, M, and Y.
			// If i \in M, then Xmin[i] = 0 and Xmax[i] = 1.
			// If i \in N, then Xmin[i] = Xmax[i] = 0.
			// If i \in Y, then Xmin[i] = Xmax[i] = 1.


// FUNCTIONS
	void partition_cardinal(bool print);
		// Solves the 'partition' problem for non-binary variables
		// Will find minimal and maximal value each variable obtains in any of the optimal solutions

	double solve_partition_cardinal(int k, bool fill_Xmin, bool print);
		// Solves a modified version of the problem during the partitioning when X-variables are not binary (cardinal preferences)
		// We will find value of variable X[k], and store it in Xmin/Xmax (depending on 'fill_Xmin')

	void check_solution_in_S_cardinal(solution s, bool print);
		// Will add solution 's' to set S if it is not already there.



private:
	// AUXILIARY FUNCTIONS
	void change_obj_function(std::vector<double> w);
	void print_partition(double length);

	bool find_identical_agents_aid(int i, int j, bool print);
		// Help function, solves the integer program in which x_i=1, x_j=0 is an optimal solution
		// and returns true if x_j=1,x_i = 0 is also an optimal solution.
	
	void identical_cluster(bool print);
		// Help function, will identify clusters of identical agents, based on 'identical'.

	// Auxiliary functions for constructors
	void initializeGurobi(bool print);
	void initializeVariables(inst I_in, bool print);

	int max_Y_value;
		// Auxiliary variable used to determine an appropriate M-value in th constraints to find identical agents when Y can be integer.
	
	int max_Y_increment;
		// If 'M_l' appears to be not high enough in one of the constraints, we will raise the max_Y_value by this much.
		// Both parameters are initialized in 'initializeVariables' function.

	void read_MIPLIB_info(bool print);
	
};

