The data files are formatted as follows:

FIRST LINE
- Number of constraints
- Number of X-variables (binary variables that each represent an agent who we want to treat fairly in the selection of the optimal solution)
- Number of Y-variables (binary/integer auxiliary variables)
- B or I: by writing I, you declare the Y-variables to be integer. By writing B, you declare the Y-variables to be binary

SECOND LINE
Coefficients in the objective function (first the coefficients of the X-variables, then the coefficients of the Y-variables)

THIRD LINE
Right-hand sides of the constraints, where we assume that the weighted sum of the decision variables is LESSER THAN OR EQUAL to the right-hand side.

REMAINING LINES
The coefficient matrix of the constraints, where we again list the coefficients of the X-variables before those of the Y-variables.