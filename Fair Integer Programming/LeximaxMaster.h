#include "IPSolver.h"

#pragma once
class LeximaxMaster
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
	
	GRBConstr* C;
	// Constraints that determine the constraints to link w with MIN_M.
	GRBConstr C_Sum1;
	// Sum of the weights equals one

	std::vector<GRBColumn> columns;
	// Contains the constraint coefficients by column (solution)

	GRBColumn Epsilon_column;
	// This variable is used to determine whether an agent should be selected with a certain probability
	GRBVar Epsilon;
	// Variable associated with 'Epsilon_column'


	std::vector<double> dual_C;
	std::vector<double> dual_C_Sum1;

	std::vector<double> prob_agents;
	// Contains the probabilities with which the agents in M are selected in the leximax distribution

	lottery solve(bool print);
	void getDualValues(bool print);
	void addColumn(solution sol, bool print);


	// Auxiliary functions
	void defineModelConVar(bool print);

	LeximaxMaster(IPSolver* K_in, bool print);
	~LeximaxMaster();
};

