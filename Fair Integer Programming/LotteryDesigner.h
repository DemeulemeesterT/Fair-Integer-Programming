#include "IPSolver.h"
#include "LeximaxMaster.h"
#include "SimplicalDecomposition.h"

#pragma once
class LotteryDesigner
{
public:
	// Variables
	IPSolver* K;
	IPSolver* K_initial;

	// Main methods
	lottery uniform(int iterations, bool print);
	// 'uniform' selects each agent with a probability equal to the fraction of the optimal solutions in which they are selected
	lottery leximax(bool print);
	// 'leximax' finds a lottery that maximizes the lowest selection probability among the agents in M, and then the second-lowest...
	lottery RSD(int iterations, bool print, unsigned seed = 123456789);
	// 'RSD' randomly orders the agents and let them discard the solutions in which they are not selected one by one in this order.
	lottery RSD_once(bool print, unsigned seed = 123456789);
	// Only runs RSD once.

	lottery rename_variables(int iterations, bool print, unsigned seed = 123456789);
	// Re-index the variables and solve the program again with the variables in this order.
	lottery perturb_objective(int iterations, bool print, unsigned seed = 123456789);
	// Randomly add small perturbations to the objective function
	lottery Nash(int iterations, bool print);
	// 'Nash' finds the point in the convex hull of the optimal solutions with the maximum geometric mean.
		// In other words, it finds the point p = (p_1, ..., p_n) such that p_1 * p_2 * ... * p_n is maximized

	lottery Nash_LP_relax(bool print);
		// Find the solution of maximum Nash welfare over the LP relaxation of the original knapsack problem.
		// We are not sure if the resuling probabilities lie in the convex hull of the optimal solutions...

	std::vector<lottery> compare_methods(std::string s, int iterations, bool print_report, bool interactive_console = false, unsigned seed = 123456789);
	// The 'input' string contains the letters of the methods we want to see compared
	// It will print the results in this order
	// Letters are
		// U = Uniform
		// L = Leximax
		// R = RSD
		// O = Rename Variables ("Order")
		// P = Perturb objective
	//Will output a vector of the found lotteries, in the demanded order
	// 'print_report' is true: a table with the resulting distributions will be printed in the console.

// Constructors and destructors
	LotteryDesigner(IPSolver* K_in, bool print);
	~LotteryDesigner();

private:
	// Auxiliary methods
	void print_lottery(lottery L);
	int count_non_zero_weights(lottery L);
	double lowest_selection_prob_M(lottery L); // Finds the lowest selection probability for the agents in M
	double highest_selection_prob_M(lottery L);
	void reset_K(bool print);
	std::vector<std::vector<int>> order_agents(std::vector<int> v, int iterations, bool print, unsigned seed = 123456789);
	// Returns a vector of different orders of the agents in the vector 'v'.
		// If the number of agents in M is smaller than or equal to 8, it returns all orders
		// Otherwise, it returns a sample of the orders, size 'iterations'
};

