from Numberjack import *


# Golomb Ruler --- CSPLib prob006

# A Golomb ruler may be defined as a set of m marks/integers 0 = a_1 < ... < a_m
# such that the pairwise differences between marks are distinct. The objective
# is to find optimal (minimum length) rulers.


def get_model(param):
    m = param['marks']
    n = 2 ** (m - 1)

    marks = VarArray(m, n, 'm')
    distance = [marks[i] - marks[j] for i in range(1, m) for j in range(i)]

    model = Model(
        Minimise(marks[-1]),  # objective function

        [marks[i-1] < marks[i] for i in range(1, m)],
        AllDiff(distance),
        marks[0] == 0,  # symmetry breaking

        [distance[i * (i - 1) // 2 + j] >= ruler[i - j] for i in range(1, m) for j in range(0, i - 1) if (i - j < m)]
    )

    return marks, model


def solve(param):
    marks, model = get_model(param)

    solver = model.load(param['solver'], marks)
    solver.setHeuristic(param['heuristic'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])

    solver.solve()

    out = ''
    if solver.is_sat():
        out = str(marks)
    out += ('\nNodes: ' + str(solver.getNodes()))
    out += '\nTime: ' + str(solver.getTime())
    return out


ruler = (0, 1, 3, 6, 11, 17, 25, 34, 44, 55, 72, 85, 106, 127)
default = {'solver': 'Mistral', 'marks': 6, 'heuristic': 'Impact', 'verbose': 0, 'tcutoff': 60}


if __name__ == '__main__':
    param = input(default)
    print(solve(param))
