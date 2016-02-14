from __future__ import division, print_function
from Numberjack import *


# This example demonstrates finding maximally different solutions. Each
# iteration forbids any previous solution and minimises similarity to the
# previous solution. The example model simply asks for N Boolean variables, with
# (N/2) of them set to 1, but any other satisfaction model could be used also.


def solve(param):
    N = param["N"]
    previoussolutions = []

    while True:
        decsionvars = VarArray(N)  # Array of N Boolean variables
        model = Model(Sum(decsionvars) == (N // 2))

        # Forbid previous solutions
        for values in previoussolutions:
            model += Disjunction([x != v for x, v in zip(decsionvars, values)])

        # Minimise the similarity to the last solution
        if previoussolutions:
            model += Minimise(Sum([x == v for x, v in zip(decsionvars, previoussolutions[-1])]))

        solver = model.load(param["solver"])
        solver.solve()

        if solver.is_sat():
            values = [x.get_value() for x in decsionvars]
            print("Solution:", values)
            previoussolutions.append(values)

        elif solver.is_unsat():
            print("No more solutions.")
            break

    print("Found a total of %d solutions." % len(previoussolutions))


if __name__ == '__main__':
    default = {'solver': 'Mistral', 'N': 6}
    param = input(default)
    solve(param)
