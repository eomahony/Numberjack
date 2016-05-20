from Numberjack import *


# Costas Array

# The costas array problem is to place N points on an N * N board such that each
# row and column contains only one point, and the pairwise distances between
# points is also distinct. i.e. such that each row of the triangular distance
# matrix constains no repeat distances.


def get_model(N):
    # Create the variables
    sequence = VarArray(N, 1, N)

    # State the model
    model = Model(
        AllDiff(sequence),
        [AllDiff([sequence[j] - sequence[j + i + 1] for j in range(N - i - 1)]) for i in range(N - 2)]
    )
    return sequence, model


def printCostasTriangle(sequence):
    N = len(sequence)
    out = ''.join([str(int(var.get_value())).rjust(3) for var in sequence])+'\n'
    for i in range(N-1):
        out += ''.join([str(int(sequence[j].get_value())-int(sequence[j+i+1].get_value())).rjust(3)
                        for j in range(N - i - 1)])+'\n'
    return out


def solve(param):
    sequence, model = get_model(param['N'])

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])

    solver.solveAndRestart()

    out = ''
    if solver.is_sat():
        out += printCostasTriangle(sequence)
    out += '\nNodes: ' + str(solver.getNodes())
    return out


default = {'solver': 'Mistral', 'N': 6, 'verbose': 1, 'tcutoff': 3}


if __name__ == '__main__':
    param = input(default)
    print(solve(param))
