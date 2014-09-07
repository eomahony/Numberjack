from Numberjack import *

# An example of mixing the SAT encoding of a model, each variable is encoded in
# a different manner and the AllDiff is encoded in multiple ways.
# See NJEncodings for a list of predefined EncodingConfiguration.


def solve(param):
    model = Model()
    v1, v2, v3 = VarArray(3, 3)

    # Pure direct/sparse
    v1.encoding = EncodingConfiguration(direct=True, order=False)

    # Regular
    v2.encoding = EncodingConfiguration(direct=True, order=True)

    # Full regular
    v3.encoding = EncodingConfiguration(direct=False, order=True)

    # Encode the supports of a disequality
    neq = v1 != v2
    neq.encoding = EncodingConfiguration(conflict=False, support=True, order=False)
    model += neq

    # Take the default for this
    model += v3 > v2

    # Encode an AllDiff using the pairwise decomposition into disequalities,
    # and a ladder encoding of an at-most-one over the values, and add
    # pigeon hole constraints
    variables = VarArray(3, 3)
    cons1 = AllDiff(variables)
    cons1.encoding = EncodingConfiguration(
        alldiff_encoding=AllDiffEncoding.PairwiseDecomp |
        AllDiffEncoding.LadderAMO | AllDiffEncoding.PigeonHole)
    model += cons1

    solver = model.load(param["solver"])
    solver.output_cnf("mixedencoding.cnf")
    solver.setTimeLimit(param["tcutoff"])
    solver.setVerbosity(param["verbose"])
    solver.solve()

    if solver.is_sat():
        print v1.get_value(), v2.get_value(), v3.get_value()
        print [v.get_value() for v in variables]
    elif solver.is_unsat():
        print "UNSATISFIABLE"
    else:
        print "UNKNOWN"


if __name__ == "__main__":
    default = {
        "solver": "MiniSat", "tcutoff": 3600, 'verbose': 1,
    }
    param = input(default)
    solve(param)
