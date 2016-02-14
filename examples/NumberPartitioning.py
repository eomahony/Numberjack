# Number Partitioning
#
# In this problem, the numbers 1 to N must be arranged in two
# different groups A and B in such a way that:
# 1. A and B have the same cardinality.
# 2. The sum of numbers in A equals the sum of numbers in B.
# 3. The sum of the squares of numbers in A equals the sum of
# the squares of numbers in B.
#
# CSPLib Problem 049 - http://www.csplib.org/Problems/prob049/
#
# Note: There is no solution for N < 8! Also, there is no solution
# if N is not a multiple of 4! (Source: CSPLib - see above link.)

from __future__ import print_function, division
from Numberjack import *


def get_model(N):
    assert N % 2 == 0, "N should be even"
    # A and B will have the same cardinality:
    A = VarArray(N // 2, 1, N)
    B = VarArray(N // 2, 1, N)
    sumtotal = N*(N+1)//4

    model = Model(
        # Each of the numbers 1 to N must be present exactly once.
        AllDiff([x for x in A+B]),

        # The sum of numbers in A equals the sum of numbers in B.
        Sum(A) == Sum(B),

        Sum(A) == sumtotal,
        Sum(B) == sumtotal,

        # The sum of the squares of numbers in A equals the sum of
        # the squares of numbers in B.
        Sum([x*x for x in A]) == Sum([y*y for y in B])
    )

    # Symmetry breaking
    model += A[0] == 1
    model += A[0] < B[0]
    for i in range(N // 2 - 1):
        model += A[i] < A[i + 1]
        model += B[i] < B[i + 1]

    return A, B, model


def solve(param):
    N = param['N']

    A, B, model = get_model(N)

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.solveAndRestart()

    if solver.is_sat():
        a = [x.get_value() for x in A]
        b = [x.get_value() for x in B]
        print("A: " + str(A), sum(a))
        print("B: " + str(B), sum(b))
    elif solver.is_unsat():
        print("Unsatisfiable")
    else:
        print("Timed out")
    print("%d nodes" % solver.getNodes())


if __name__ == '__main__':
    default = {'N': 8, 'solver': 'Mistral', 'verbose': 0}
    param = input(default)
    solve(param)

