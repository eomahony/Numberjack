# Killer Sudoku
#
# Killer Sudoku is a Sudoku puzzle with a twist.
#
# The standard Sudoku puzzle is a 9x9 matrix to be filled out
# in such a way that each row, column and each of the 3x3
# sub-matrices would contain the numbers 1 to 9.
#
# In Killer Sudoku, there is one extra constraint. There are
# so-called "cages" which are much like the 3x3 sub-matrices.
# They may contain at most one of each number (1 to 9) and the
# numbers' sum must be equal to the cage's number (see CSPLib link).
#
# CSPlib Problem 057 - http://www.csplib.org/Problems/prob057/

from __future__ import print_function
from Numberjack import *


def get_model(param):
    N = param['N']
    cages = parsecages(examplecages)

    grid = Matrix(N*N, N*N, 1, N*N)

    model = Model(
        [AllDiff(row) for row in grid.row],
        [AllDiff(col) for col in grid.col],

        # Contents of sub-matrices must be distinct.
        [AllDiff(grid[x:x + N, y:y + N]) for x in range(0, N*N, N)
            for y in range(0, N * N, N)],
    )

    # Contents of cages must be distinct.
    for cagetotal, cellids in cages:
        cagecells = [grid.flat[x] for x in cellids]
        model += AllDiff(cagecells)
        model += Sum(cagecells) == cagetotal

    return grid, model


def solve(param):
    grid, model = get_model(param)

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.solve()

    if solver.is_sat():
        print(str(grid))
    elif solver.is_unsat():
        print("Unsatisfiable")
    else:
        print("Unknown")


def parsecages(textblob):
    cages = []
    for line in textblob.split("\n"):
        if not line.strip():
            continue
        bits = list(map(int, line.split()))
        cages.append((bits[0], bits[1:]))
    return cages


examplecages = """
3 0 1
15 2 3 4
22 5 13 14 22
4 6 15
16 7 16
15 8 17 26 35
25 9 10 18 19
17 11 12
9 20 21 30
8 23 32 41
20 24 25 33
6 27 36
14 28 29
17 31 40 49
17 34 42 43
13 37 38 46
20 39 48 57
12 44 53
27 45 54 63 72
6 47 55 56
20 50 59 60
6 51 52
10 58 66 67 75
14 61 62 70 71
8 64 73
16 65 74
15 68 69
13 76 77 78
17 79 80
"""


if __name__ == '__main__':
    default = {'N': 3, 'solver': 'Mistral', 'verbose': 0}
    param = input(default)
    solve(param)

