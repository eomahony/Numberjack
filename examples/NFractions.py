# N Fractions - non-generalized version
#
# You must assign the digits 1 to 9 to the variables
# A, B, C, D, E, F, G, H, I, such that each digit may
# occur only once, and the following equation is satisfied:
#
#  A      D      G
# ---- + ---- + ---- = 1
#  BC     EF     HI
#
# CSPlib Problem 041 - http://www.csplib.org/Problems/prob041/

from __future__ import print_function
from Numberjack import *


def modelNFractions():
    a, b, c, d, e, f, g, h, i = VarArray(9, 1, 9)

    model = Model(
        # Rearrange the equation to avoid the usage of floats.
        a * (10 * e + f) * (10 * h + i) +
        d * (10 * b + c) * (10 * h + i) +
        g * (10 * b + c) * (10 * e + f) ==
        (10 * b + c) * (10 * e + f) * (10 * h + i),

        # Each digit must occurr once
        AllDiff([a, b, c, d, e, f, g, h, i]))

    return (a, b, c, d, e, f, g, h, i, model)


def solve(param):
    (a, b, c, d, e, f, g, h, i, model) = modelNFractions()
    solver = model.load(param["solver"])
    solver.solve()

    if solver.is_sat():
        print("(" + str(a) + " / " + str(b) + str(c) + ") + " + \
              "(" + str(d) + " / " + str(e) + str(f) + ") + " + \
              "(" + str(g) + " / " + str(h) + str(i) + ") = 1")
    elif solver.is_unsat():
        print("Unsatisfiable")
    else:
        print("Timed out")


if __name__ == '__main__':
    default = {'solver': 'Mistral', 'verbose': 1}
    param = input(default)
    solve(param)

