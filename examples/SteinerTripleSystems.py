# Steiner Triple Systems

# The ternary Steiner problem of order n consists of finding a set of n(n-1)/6
# triples of distinct integer elements in {1,...,n} such that any two triples
# have at most one common element. It is a hypergraph problem coming from
# combinatorial mathematics where n modulo 6 has to be equal to 1 or 3.

# CSPLib Problem 044 - http://www.csplib.org/Problems/prob044/

from __future__ import print_function, division
from Numberjack import *


def get_model(N):
    assert N >= 3, "Error: N must be at least 3."
    assert N % 6 in [1, 3], "Error: N % 6 must be 1 or 3."

    N_ROWS = N * (N - 1) // 6
    matrix = Matrix(N_ROWS, N, 0, 1)

    model = Model()
    for row in matrix.row:
        model += Sum(row) == 3

    for i in range(N_ROWS-1):
        for j in range(i + 1, N_ROWS):
            model += Sum([matrix[i][k] * matrix[j][k] for k in range(N)]) <= 1

        # Symmetry breaking
        model += LeqLex(matrix.row[i], matrix.row[i+1])

    return matrix, model


def solve(param):
    matrix, model = get_model(param['N'])

    solver = model.load(param['solver'])
    solver.setHeuristic(param["var"], param["val"])
    solver.setVerbosity(param['verbose'])
    solver.solveAndRestart()

    if solver.is_sat():
        for row in matrix:
            triple = [i + 1 for i, x in enumerate(row) if x.get_value() == 1]
            print(triple)
    elif solver.is_unsat():
        print("Unsatisfiable")
    else:
        print("Unknown")
    print(solver.getNodes())


if __name__ == '__main__':
    default = {'N': 7, 'solver': 'Mistral', 'verbose': 0, "var": "Lex", "val": "Lex"}
    param = input(default)
    solve(param)
