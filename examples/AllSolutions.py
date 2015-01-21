from Numberjack import *


def solve():
    xs = VarArray(5, 1, 5)
    model = Model(AllDiff(xs))

    solver = model.load("Mistral")
    solver.startNewSearch()
    while solver.getNextSolution() == SAT:
        print xs
        print "Nodes:", solver.getNodes()
        print "Backtracks:", solver.getBacktracks(), "\n"


if __name__ == '__main__':
    solve()
