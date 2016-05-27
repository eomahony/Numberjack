from __future__ import division, print_function
from Numberjack import *


# This example demonstrates finding all solutions to a given model.
# solver.startNewSearch() should be called to set up internal data structures in
# the solver first, then call solver.getNextSolution() repeatedly until it
# returns a value other than SAT.


def solve(param):
    num_solutions = 0
    N = param["N"]
    decsionvars = VarArray(N)  # Array of N Boolean variables
    model = Model(Sum(decsionvars) == (N // 2))

    solver = model.load(param["solver"])
    solver.startNewSearch()  # Must be called to set up internal data structures
    while solver.getNextSolution() == SAT:
        values = [x.get_value() for x in decsionvars]
        print("Solution:", values)
        num_solutions += 1

    print("Found a total of %d solutions." % num_solutions)


if __name__ == '__main__':
    default = {'solver': 'MiniSat', 'N': 4}
    param = input(default)
    solve(param)
