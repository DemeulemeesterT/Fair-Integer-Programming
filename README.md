# Fair integer programming

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

## Minimal working example for kidney exchange instances
This is a minimal working example, which can simply be pasted into the `source.cpp` file:  

```
using namespace std;

int main()
{
	KidneyExchange* KE = new KidneyExchange("40-instance-48", true);

	inst I = KE->generate_instance(true, false);
	
	IPSolver* K = new IPSolver(I, true);
	
	K->analyze(false);
	
	LotteryDesigner* L = new LotteryDesigner(K, false);

	std::vector<lottery> results = L->compare_methods("ULRC", 1000, true, false, 0);
	
	delete K;	
}
```

## Use
### Creating an `inst` object for your formulation
To use the code, the integer programming formulation must be transformed into an object of the structure `inst`. An instance of an `inst` object contains the following elements
* `m`: number of constraints in the model (integer)
* `n`: number of agents, which we refer to as $X$-variables (integer)
* `t`: number of auxiliary variables, which we refer to as $Y$-variables (integer)
* `n_var`: total number of variables (`n_var`=`n`+`t`)
* `v`: vector of objective coefficients, which first contains the coefficients of the $X$-variables, and then of the $Y$-variables
* `C`: the right-hand sides of the constraints
* `A`: the constraint matrix, which is of dimension `m`(`n`+`t`)
* `Y_bool`: boolean which is `true` if the $Y$-variables are binary
* `Y_coeff_zero`: boolean which is `true` if the $Y$-variables are all equal to zero
* `Y_integer`: boolean which is `true` if the $Y$-variables are all integers
* `X_bool`: boolean which is `true` if the $X$-variables are all binary
* `X_integer`: boolean which is `true` if the $X$-variables are all integer
* `data_name`: name of the formulation (optional)

There are several ways to create an object of the `inst` structure, but this process is problem-specific and formulation-specific. The reason for this is that it is not clear from a given formulation which variables should be considered as $X$-variables (because they represent the agents in a certain problem), and which variables should be considered as auxiliary variables. For the cycle formulation of the kidney exchange problem, a code has been written to read data from a file (`KidneyExchange::read_data()`) and create the corresponding `inst` object (`KidneyExchange::generate_instance()`). This code can be used to guide you into writing your own code for a specific problem and a specific formulation.

### Creating an `IPSolver*` object from your `inst` object
Once you have your formulation as an `inst` object `I`, you should create an object of the class `IPSolver*` (this is a dynamic object) through the following command:  

```
IPSolver* K = new IPSolver(I, print);
```  

(The boolean `print` refer to whether or not you want to print output of the function on screen for the remainder of the code.)

### Analyzing your formulation
Now that we have the `IPSolver*` object, we can start analyzing our formulation.  
* For binary $X$-variables (dichotomous preferences), we can find out using Proposition 1 from the paper which variables are selected in all, none, or in some but not in all of the optimal solutions.
* For non-binary $X$-variables (cardinal preferences), we can find out the largest and the smallest value of each of the $X$-variables in any of the optimal solutions.

For both analyses, we can simply call the following function of our `IPSolver*` object `K`:  

```
K->analyze(print);
```

### Finding the distributions over the optimal solutions
Lastly, we want to find distributions over the optimal solutions. This can be done easily by the following two commands:

```
LotteryDesigner* L = new LotteryDesigner(K, print);
L->compare_methods(string, 10, interactive_display, print, 0);
```

* `string` is a string that contains information about which distributions you want to compute. The most relevant distributions are referred to by the following letters:
     * `U`: Uniform
     * `L`: Leximin
     * `C`: Nash Maximum Welfare
     * `R`: Random Serial Dictatorship (RSD)  
  A string such as `LCR` will compute the leximin, Nash, and RSD distributions, for example.
* `interactive_display` is a boolean variable. If entered as `true`, you can see the menu of possible distributions to choose from.

The resulting output could look this for an instance with binary $X$-variables, where, for example, the uniform, leximin, Nash, and RSD distributions have been computed.  

![Example Output Distributions](https://github.com/DemeulemeesterT/Fair-Integer-Programming/assets/59369043/954b0065-070a-4991-a18c-f5221bfc5675)


### Extracting selection probabilities of optimal solutions
To view the selection probabilities of the optimal solutions for the computed distributions, you can query information from the created vector of `lottery` objects that is returned by the function `compare_methods` in the object of the class `Lotterydesigner*`.




