from Numberjack import *


def get_model():
    # Create the model
    model = Model()

    #Variable definitions, we need a variable for every letter
    s, m = VarArray(2, 1, 9)  # These can't be zero as they are the start of a word
    e, n, d, o, r, y = VarArray(6, 0, 9)  # These can

    model.add(      s*1000 + e*100 + n*10 + d +
                    m*1000 + o*100 + r*10 + e ==
          m*10000 + o*1000 + n*100 + e*10 + y)

    # Post the all different constraint on all the variables
    model.add(AllDiff((s, e, n, d, m, o, r, y)))

    return s, e, n, d, m, o, r, y, model


def solve(param):
    s, e, n, d, m, o, r, y, model = get_model()

    # Load up model into solver
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])

    # Now Solve
    solver.solve()
    if solver.is_sat():
        print("    %d %d %d %d" % (s.get_value(), e.get_value(), n.get_value(), d.get_value()))
        print("+   %d %d %d %d" % (m.get_value(), o.get_value(), r.get_value(), e.get_value()))
        print("= %d %d %d %d %d" % (m.get_value(), o.get_value(), n.get_value(), e.get_value(), y.get_value()))


if __name__ == '__main__':
    default = {'solver': 'Mistral', 'verbose': 0}
    param = input(default)
    solve(param)
