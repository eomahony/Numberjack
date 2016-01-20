# All-Interval Series
#
# Arrange the numbers 1 to n in a vector s, such that the
# interval vector v=(|s2-s1|,|s3-s2|,...,|sn-sn-1|) is a
# permutation of {1,2,...,n-1}.
#
# CSPLib Problem 007 - http://www.csplib.org/Problems/prob007/
from __future__ import print_function

from Numberjack import *


def get_model(N):
    series = VarArray(N, 0, N - 1)

    model = Model(
        AllDiff(series),
        AllDiff([Abs(series[i] - series[i - 1]) for i in range(1, N)]))

    return series, model


def solve(param):
    N = param['N']

    series, model = get_model(N)

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.solve()

    if solver.is_sat():
        print("Solution:", str(series))
        print("Absolute differences:")
        print([abs(series[i].get_value() - series[i - 1].get_value())
               for i in range(1, N)])
    elif solver.is_unsat():
        print("Unsatisfiable")
    else:
        print("Timed out")


if __name__ == '__main__':
    default = {'N': 12, 'solver': 'MiniSat', 'verbose': 1}
    param = input(default)
    solve(param)

