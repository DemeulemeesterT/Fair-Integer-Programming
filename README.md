# Fair integer programming under dichotomous preferences

This code can be used to compute distributions over the optimal solutions of an integer programming formulation when each of the agents involved can be represented by a single variable (which can be binary, integer, or real).

The distributions that can be computed without changes to the code are:
* Uniform
* Leximin
* Random Serial Dictatorship (RSD)
* Nash Maximum Welfare

Moreover, it can be used to compute the Perturb and the Re-index heuristics. We refer the reader to the corresponding paper for more details about these distribution rules, and the underlying algorithms.

## Initialization
The code uses the Gurobi optimization software (version 10.0.0), and calls it from C++ by using Visual Studio 2019. Gurobi offers a free [academic license][2]. More information about how to configure a Gurobi C++ project with Microsoft Visual Studio can be found [here][1].

[1]: https://support.gurobi.com/hc/en-us/articles/360013194392-How-do-I-configure-a-new-Gurobi-C-project-with-Microsoft-Visual-Studio-2017-
[2]: https://www.gurobi.com/academia/academic-program-and-licenses/

## Use
To use the code, the integer programming formulation must be transformed into an object of the structure `inst`. An instance of an `inst` object contains the following elements
* `m`: number of constraints in the model (integer)
* `n`: number of agents, which we refer to as $X$-variables (integer)
* `t`: number of auxiliary variables, which we refer to as $Y$-variables (integer)
* `n_var`: total number of variables (`n_var`$=n+t$)
* `v`: vector of objective coefficients, which first contains the coefficients of the $X$-variables, and then of the $Y$-variables
* `C`: the right-hand sides of the constraints
* `A`: the constraint matrix, which is of dimension $m\cdot(n+t)$
* `Y_bool`: boolean which is `true` if the $Y$-variables are binary
* `Y_coeff_zero`: boolean which is `true` if the $Y$-variables are all equal to zero
* `Y_integer`: boolean which is `true` if the $Y$-variables are all integers
* `X_bool`: boolean which is `true` if the $X$-variables are all binary
* `X_integer`: boolean which is `true` if the $X$-variables are all integer
* `data_name`: name of the formulation (optional)




Most parameters for the code can be modified from the `Source.cpp` file. In that file, you can also declare for which type of symmetric $\alpha$-hedonic game and for which parameter values you want to generate an instance, if it exists. 
* `q_vector` contains the values of $q$-size core stability for which you want to generate an instance.
* `c_vector` contains the values of $c$-improvement core stability for which you want to generate an instance (denoted by $k$ in the paper).
* `n_vector` contains the sizes of the blocking coalition (denoted by $m$ in the paper).

You can either use the function `alpha_hedonic_game()` and specify the $\alpha$-vector,or use the already existing functions for symmetric fraction hedonic games, symmetric modified fractional hedonic games, and symmetric additively separable hedonic games already exists.

In the code, the default setting is that $u_i(\mathcal{C}) = 1$ for all agents $i$ in the blocking coalition, and that the utilities of the agents, which are denoted by the `X`-variables in the code, are real numbers greater than zero. However, both assumptions can be changed in the `alpha_hedonic_game.cpp` file (or the correct file for alternative subclasses of symmetric hedonic games).

* To change $u_i(\mathcal{C})$, comment out the `model.addConstr(V[i] == 1);` line. To restrict the solution space, we recommend fixing `V[i]`, which represents $u_i(\mathcal{C})$, for at least one of the agents.
* To restrict the utilities to be binary, for example, replace the line `X[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name_x);` by `X[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_BINARY, name_x);`. One could also impose alternative lower or upper bounds on the allowed utilities by modifying the first two values in the declaration of the `X`-variables.
