from Numberjack import *


# Magic Square --- CSPLib prob019

# A Magic Square of order N is an N x N matrix of values from 1 to N^2, where
# each run, column, and diagonal sum to the same value. This value can be
# calculated as N * (N^2 + 1) / 2.


def get_model(N):
    sum_val = N * (N * N + 1) // 2  # This is what all the columns, rows and diagonals must add up tp

    square = Matrix(N, N, 1, N * N)
    model = Model(
        AllDiff(square),

        [Sum(row) == sum_val for row in square.row],
        [Sum(col) == sum_val for col in square.col],

        Sum([square[a, a] for a in range(N)]) == sum_val,
        Sum([square[a, N - a - 1] for a in range(N)]) == sum_val
    )
    return square, model


def solve(param):
    square, model = get_model(param['N'])
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setHeuristic(param['var'], param['val'], param['rand'])
    solver.setTimeLimit(param['cutoff'])

    if param['restart'] == 'yes':
        solver.solveAndRestart()
    else:
        solver.solve()

    out = ''
    if solver.is_sat():
        out = str(square)
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out


default = {'solver': 'Mistral', 'N': 4, 'var': 'MinDomain',
           'val': 'RandomMinMax', 'restart': 'yes', 'rand': 2, 'verbose': 0, 'cutoff': 10}


if __name__ == '__main__':
    param = input(default)
    print(solve(param))
