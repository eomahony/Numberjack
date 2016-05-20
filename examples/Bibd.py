from Numberjack import *

# Balanced Incomplete Block Design (BIBD) --- CSPLib prob028

# A BIBD is defined as an arrangement of v distinct objects into b blocks such
# that each block contains exactly k distinct objects, each object occurs in
# exactly r different blocks, and every two distinct objects occur together in
# exactly lambda blocks. Another way of defining a BIBD is in terms of its
# incidence matrix, which is a v by b binary matrix with exactly r ones per row,
# k ones per column, and with a scalar product of lambda 'l' between any pair of
# distinct rows.


def get_model(v, b, r, k, l):
    matrix = Matrix(v, b)
    model = Model(
        [Sum(row) == r for row in matrix.row],  # every row adds up to r
        [Sum(col) == k for col in matrix.col],  # every column adds up to k

        # the scalar product of every pair of columns adds up to l
        [Sum([(row[col_i] * row[col_j]) for row in matrix.row]) == l
            for col_i in range(v) for col_j in range(col_i)],
    )
    return matrix, model


def solve(param):
    matrix, model = get_model(param['v'], param['b'], param['r'], param['k'], param['l'])

    if param['solver'] == 'Mistral':
        model += [matrix.row[i] <= matrix.row[i+1] for i in range(param['v']-1)]
        model += [matrix.col[i] <= matrix.col[i+1] for i in range(param['b']-1)]

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])
    solver.solve()

    out = ''
    if solver.is_sat():
        out += str(matrix)
    elif solver.is_unsat():
        out += "UNSAT"
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out


default = {'solver': 'Mistral', 'v': 7, 'b': 7, 'r': 3, 'k': 3, 'l': 1, 'verbose': 0, 'tcutoff': 10}


if __name__ == '__main__':
    param = input(default)
    print(solve(param))
