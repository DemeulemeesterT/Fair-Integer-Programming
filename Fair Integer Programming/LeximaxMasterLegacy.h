#include "IPSolver.h"

#pragma once
class LeximaxMasterLegacy
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

	std::vector<double> dual_C;
	std::vector<double> dual_C_Sum1;

	lottery solve(bool print);
	void getDualValues(bool print);
	void addColumn(solution sol, bool print);


	// Auxiliary functions
	void defineModelConVar(bool print);

	LeximaxMasterLegacy(IPSolver* K_in, bool print);
	~LeximaxMasterLegacy();
};

