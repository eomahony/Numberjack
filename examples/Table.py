from Numberjack import *


def solve(param):
    x1, x2, x3 = VarArray(3, 1, 3)
    model = Model()

    model += Table([x1, x2], [[1, 1], [2, 2], [3, 3]], type='conflict')
    model += Table([x2, x3], [[1, 2], [2, 3]], type='support')

    print model

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])

    solver.solve()
    if solver.is_sat():
        print x1.get_value(), x2.get_value(), x3.get_value()


default = {'N': 3, 'solver': 'Mistral', 'verbose': 0, 'tcutoff': 30}


if __name__ == '__main__':
    param = input(default)
    solve(param)
