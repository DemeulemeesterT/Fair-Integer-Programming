#include "IPSolver.h"

#pragma once
class SimplicalDecomposition
{
	/*Simplicial Decomposition (SD)
	
	(from Hearn, D. W., Lawphongpanich, S., & Ventura, J. A. (1987). Restricted simplicial decomposition: Computation and extensions. 
	In Computation Mathematical Programming (pp. 99-118). Springer, Berlin, Heidelberg.)

	For a MINIMIZATION problem: 
	Step O: Let x' be a feasible point of S and S' = empty (or some initial subset).
	Step 1 (Subproblem): Let y = arg min{grad(f(x'))y: y \in S}. If grad(f(x'))(y - x') >= O, stop
	and x' is optimal. Otherwise, set S' = S' U {y}, and go to Step 2.
	Step 2 (Master problem): Let x'=argmin{f(x): x \in Conv(S')}, and purge S' of
	all extreme points with zero weight in the expression of x' as a convex combination
	of elements of S'. Return to Step 1.
		Purgin the unused extremal points is optional.
	*/
public:
// VARIABLES
	IPSolver* K;
	// The solutions in 'K->S' will be used as the variables in the master problem

	GRBEnv* env;
	GRBModel* model;

	std::vector<GRBVar> w;
	// w[s] = Probability of being selected for solution s
	GRBVar* p;
	// p[i] = Probability of being selected for agent 'dict[p]'
	// P[i] variables are only created for agents i \in M (who are selected in some but not in all of the optimal solutions)
		// CAREFUL! For non-binary X-variables, the variable p[i] will be equal to the expected
		// value of the X[i]-variable in the resulting distribution
		// MINUS the dystopia point of agent i, namely the minimal outcome she experiences in any of the optimal solutions.

	std::vector<int> dict;
	// dict[i] contains which agent is represented by variable p[i]

	GRBConstr C_Sum1;
	// Sum of the weights equals one
	
	GRBConstr* C_p;
	// Constraints linking the P-variables (selection probabilities of the agents), with the W-variables (selection probabilities of the solutions)

	std::vector<GRBColumn> columns;
	// Contains the constraint coefficients by column (solution)

// FUNCTIONS
	lottery SD_Nash(bool print);
		// Returns the lottery that maximizes the Nash Social Welfare solution using the simplical decomposition technique

	lottery Nash_CG(bool print);
		// Adopts the column generation procedure similar to Flanigan et al. (2021)
		// Two differences with SD_Nash:
			// Objective function of the pricing problem is simply the sum of the selected dual variables
			// Stopping criterion is different: check if there's a solution that's already included
				// as a variable in the restricted master problem with a higher objective value. If yes, STOP. (p. 19 of supplementary material Flanigan et al. (2021)

	std::vector<double> gradientNash(std::vector<double> p_values, bool print);
		// Returns the value of the gradient of the Nash Social Welfare objective function, evaluated at point 'p_values'

	void addColumn(solution sol, bool print);

// CONSTRUCTORS AND DECONSTRUCTORS
	SimplicalDecomposition(IPSolver* K_in, bool print);
		// Importantly, we will only add the columns (i.e. a subset of the optimal solutions) to the model
		// when we know which objective function we will optimize.
		// The reason for this is that some objective function require the initial set of columns to satisfy some criteria in order to function properly.
		// For the Nash Social Welfare objective function, for example, which is the product of the selection probabilities of the agents in M,
			// all agents in M should be selected in at least one of the optimal solutions.

	~SimplicalDecomposition();
};

