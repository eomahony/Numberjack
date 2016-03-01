/* cost_functions.mzn
-------------------------
Provides a default decomposition of cost functions
for weighted csp solvers. These map tuples in 
extensional form to cost values that are then added to an overall violation.

Example:
x, y | cost
0  0 | 4
0  1 | 3
1  0 | 2 
1  1 | 4

for a binary cost function. The assignment (x -> 1, y -> 0) would yield minimal overall
violation (2).

Toulbar2 has dedicated algorithms for such cost functions and can be accessed via Numberjack.
For all other solvers, this decomposition provides an encoding using table constraints
to represent arbitrary (finite-domain, up to ternary) functions

AS, ISSE Augsburg, part of MiniBrass, 2016 
*/


% costs should provide the cost value for each domain value of x (in increasing order) 
% e.g. x \in 0..2 ->, costs = [4, 1, 3] would be interpreted as f(x -> 0) = 4, f(x -> 1) = 1, f(x -> 2) = 3
predicate cost_function_unary(var int: x, array[int] of int: costs, var int: costVariable );

% costs should be provided for each element in the cartesian product of two domains, in increasing lexicographic order
% e.g. 
% x, y | cost
% 0  0 | 4
% 0  1 | 3
% 1  0 | 2 
% 1  1 | 4
% would be represented by the cost vector [4, 3, 2, 4]
predicate cost_function_binary(var int: x, var int: y, array[int] of int: costs, var int: costVariable );

% costs should be provided for each element in the cartesian product of three domains, in increasing lexicographic order
% e.g. 
% x, y, z | cost
% 0  0  0 | 5
% 0  0  1 | 2
% 0  1  0 | 1
% 0  1  1 | 3
% 1  0  0 | 4
% 1  0  1 | 6 
% 1  1  0 | 2 
% 1  1  1 | 3
% would be represented by the cost vector [5, 2, 1, 3, 4, 6, 2, 3], a minimal-cost assignment would be 0  1, 0 
predicate cost_function_ternary(var int: x, var int: y, var int: z, array[int] of int: costs, var int: costVariable );

% The nary variant requires default_costs for all tuples and allows to list explicit exceptions, exceptions[i] should correspond to costs[i]
% e.g. 
% x, y | cost
% 0  0 | 4
% 0  1 | 3
% 1  0 | 2 
% 1  1 | 4
% could be represented by vars = [x,y], default_costs = 4, exceptions = [| 0,1 | 1, 0 |], costs = [3, 2] 
predicate cost_function_nary(array[int] of var int: vars, int: default_costs, array[int,int] of int: exceptions, array[int] of int: costs, var int: costVariable);
