# 015: Schur"s Lemma
#
# Given n balls, labelled 1 to n, put them into 3 boxes such
# that for any triple of balls (x,y,z) with x+y=z, not all are
# in the same box.
#
# CSPLib Problem 015 - http://www.csplib.org/Problems/prob015/

from __future__ import print_function
from Numberjack import *


def get_model(N, C):
    balls = VarArray(N, C)

    model = Model()

    for i in range(1, N + 1):
        for j in range(1, N - i + 1):
                model += Disjunction([
                    balls[i - 1] != balls[j - 1],
                    balls[i - 1] != balls[i + j - 1],
                    balls[j - 1] != balls[i + j - 1]])

    return balls, model


def solve(param):
    balls, model = get_model(param["N"], param["C"])

    solver = model.load(param["solver"])
    solver.setVerbosity(param["verbose"])
    solver.solve()

    if solver.is_sat():
        print(str(balls))
    elif solver.is_unsat():
        print("Unsatisifiable")
    else:
        print("Timed out")


if __name__ == "__main__":
    default = {"N": 12, "C": 3, "solver": "MiniSat", "verbose": 1}
    param = input(default)
    solve(param)

