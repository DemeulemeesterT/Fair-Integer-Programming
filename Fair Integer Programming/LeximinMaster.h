#include "IPSolver.h"

#pragma once
class LeximinMaster
{
public:
	IPSolver* K;
	// The solutions in 'K->S' will be used as the variables in the master problem

	GRBEnv* env;
	GRBModel* model;

	std::vector<GRBVar> w;
	// w[s] = Probability of being selected for solution s
	GRBVar MIN_M;
	// This variable contains the lowest selection probability for the agents in M
	
	// THE FOLLOWING CONSTRAINTS contain one element for each agent in M
		// e.g., C_bound[3] represents the C_bound constraint for the third agent in M.
	std::vector<GRBConstr> C_bound;
	// Constraints that determine the constraints to link w with MIN_M.

	GRBConstr C_Sum1;
	// Sum of the weights equals one

	GRBConstr C_bound_MIN_M;
	// Constraint will be used to decide which agents is assigned with obtained selection probabilities

	std::vector<GRBColumn> columns;
	// Contains the constraint coefficients by column (solution)

	GRBColumn Epsilon_column;
	// This variable is used to determine whether an agent should be selected with a certain probability
	GRBVar Epsilon;
	// Variable associated with 'Epsilon_column'


	std::vector<double> dual_C_bound;

	std::vector<double> dual_C_Sum1;

	std::vector<double> dual_C_bound_MIN_M;

	std::vector<double> prob_agents;
	// Contains the probabilities with which the agents in M are selected in the leximin distribution

	lottery solve(bool print);
	void getDualValues(bool epsilon_constraint, bool print);
		// 'epsilon_constraint' is true when we are calling the dual variables when the epsilon variable is present in the model
	void addColumn(solution sol, bool print);


	// Auxiliary functions
	void defineModelConVar(bool print);

	LeximinMaster(IPSolver* K_in, bool print);
	~LeximinMaster();
};

