# Magic Sequences - Magic Square variant
#
# A magic sequence of length n is a sequence of integers x[0]
# to x[n-1], such that for all i in 0 to n-1, the number i occurs
# exactly x[i] times in the sequence.
#
# For instance, 6,2,1,0,0,0,1,0,0,0 is a magic sequence since 0
# occurs 6 times in it, 1 occurs twice, 2 occurs once, ...
#
# CSPlib Problem 019 - http://www.csplib.org/Problems/prob019/

from __future__ import print_function
from Numberjack import *
import sys
if sys.version_info[0] > 2:
    xrange = range

def get_model(N):
    seq = VarArray(N, N)

    model = Model()
    for i in xrange(N):
        model += seq[i] == Sum([seq[j] == i for j in xrange(N)])

    return seq, model


def solve(param):
    N = param['N']

    seq, model = get_model(N)

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.solve()

    if solver.is_sat():
        print(str(seq))
    elif solver.is_unsat():
        print("Unsatisfiable")
    else:
        print("Timed out")


if __name__ == '__main__':
    default = {'N': 10, 'solver': 'Mistral', 'verbose': 0}
    param = input(default)
    solve(param)

