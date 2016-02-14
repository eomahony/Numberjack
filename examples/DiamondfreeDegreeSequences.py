# Diamond-free Degree Sequences
#
# Fill in a binary matrix of size n * n in such a way that
# - For every grouping of four rows, the sum of their non-symmetrical
#   values is less than or equal to 4,
# - No rows contain just zeroes,
# - Every row has a sum modulo 3,
# - The sum of the matrix is modulo 12.
# - No row R contains a 1 in its Rth column.
#
# Note on first constraint in model:
# A group of four nodes can have at most four edges between them.
# Since the matrix for this model will represent the adjacency
# matrix of the graph, we need to take into consideration the fact
# that the matrix will be symmetrical.
#
# CSPLib Problem 050 - http://www.csplib.org/Problems/prob050/
from __future__ import print_function
from Numberjack import *
from itertools import combinations


def get_model(N):
    # By definition a and b will have the same cardinality:
    matrix = Matrix(N, N)

    model = Model(
        # No rows contain just zeroes.
        [Sum(row) > 0 for row in matrix.row],

        # Every row has a sum modulo 3.
        [Sum(row) % 3 == 0 for row in matrix.row],

        # The sum of the matrix is modulo 12.
        Sum(matrix.flat) % 12 == 0,

        # No row R contains a 1 in its Rth column.
        [matrix[row][row] == 0 for row in range(N)])

    # Every grouping of 4 rows can have at most a sum of 4 between them.
    for a, b, c, d in combinations(range(N), 4):
        model += Sum([matrix[a][b], matrix[a][c], matrix[a][d],
                      matrix[b][c], matrix[b][d], matrix[c][d]]) <= 4

    # Undirected graph
    for i in range(N):
        model += matrix[i][i] == 0  # No looping arcs

        for j in range(N):
            model += matrix[i][j] == matrix[j][i]

    # Symmetry breaking
    for i in range(N-1):
        model += LeqLex(matrix.row[i], matrix.row[i+1])

    return matrix, model


def solve(param):
    N = param['N']

    matrix, model = get_model(N)

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.solve()

    if solver.is_sat():
        print(str(matrix) + '\n')

        print("Degree sequence:",)
        for row in matrix:
            print(sum([x.get_value() for x in row]))
    elif solver.is_unsat():
        print("Unsatisfiable")
    else:
        print("Timed out")


if __name__ == '__main__':
    default = {'N': 10, 'solver': 'MiniSat', 'verbose': 1}
    param = input(default)
    solve(param)

