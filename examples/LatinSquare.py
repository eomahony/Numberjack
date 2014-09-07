from Numberjack import *


def get_model(pls_filename):
    """
        Returns a QCP model for the partial Latin Square contained in
        'pls_filename'.
    """

    with open(pls_filename, "rt") as f:
        header = f.readline().split()
        n = int(header[-1].strip())

        model = Model()
        matrix = []
        columns = [[] for i in xrange(n)]

        for i, row in enumerate(f):
            bits = map(int, row.split())
            row = []
            for j, v in enumerate(bits):
                var = None
                if v >= 0:
                    var = Variable([v], "x%d.%d" % (i, j))
                else:
                    var = Variable(n, "x%d.%d" % (i, j))
                row.append(var)
                columns[j].append(var)
            model += AllDiff(row)
            matrix.append(row)

        for col in columns:
            model += AllDiff(col)
    return matrix, model


def setencodings(model, encoding, constraint):
    """
        Iterates throught all constrainsts in 'model' to sets all occurences of
        the class 'constraint' to use the encoding object 'encoding'.
    """
    for c in model.constraints:
        if isinstance(c, constraint):
            c.encoding = encoding


def solve(param):
    matrix, model = get_model(param["pls"])
    defaultencoding = None

    # Uncommend these two lines to use a custom SAT encoding.
    # setencodings(model, NJEncodings[param['alldiffencoding']], AllDiff)
    # defaultencoding = NJEncodings[param["encoding"]]
    solver = model.load(param["solver"], encoding=defaultencoding)

    # Uncomment one of these lines if you're using a SAT/MIP/MIP solver to
    # save the model in the respective formats.
    # solver.output_cnf("%s.cnf" % param["pls"])
    # solver.output_lp("%s.lp" % param["pls"])
    # solver.output_mps("%s.mps" % param["pls"])
    solver.setTimeLimit(param["tcutoff"])
    solver.setVerbosity(param["verbose"])
    solver.solve()

    f = sys.stdout
    if solver.is_sat():
        print >> f, "s SATISFIABLE"
        for row in matrix:
            print >> f, "v", " ".join("%2d" % v.get_value() for v in row)
    elif solver.is_unsat():
        print >> f, "s UNSATISFIABLE"
    else:
        print >> f, "s UNKNOWN"
    print >> f, "c Nodes %d" % solver.getNodes()
    print >> f, 'c Failures', solver.getFailures()
    print >> f, 'c NumVariables', solver.getNumVariables()
    print >> f, 'c NumConstraints', solver.getNumConstraints()
    print >> f, "c EncodeTime %.5f" % solver.load_time
    print >> f, "c SolveTime %.5f" % solver.getTime()


if __name__ == "__main__":
    default = {
        "solver": "MiniSat",
        "pls": "data/qwh-5-5.pls",
        "encoding": "directorder",
        "alldiffencoding": "direct",
        "tcutoff": 3600, 'verbose': 1,
    }
    param = input(default)
    solve(param)
