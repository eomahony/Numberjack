from Numberjack import *
import SCIP
import MipWrapper

vars = VarArray(4, 1,4)

N = 4
sequence = vars

diffs1 = [vars[i] - vars[i+1] for i in range(3)]
diffs2 = [vars[0] - vars[2], vars[1] - vars[3]]

model = Model(
        AllDiff(vars),
        AllDiff(diffs1),
        AllDiff(diffs2)
        )

print model

solver = SCIP.Solver(model)
print solver.solve()

print vars
print [d.get_value() for d in diffs1]
print [d.get_value() for d in diffs2]

#MipWrapper.Solver(model).solve()