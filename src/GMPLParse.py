from Numberjack import *
import OsiGlpk as Osi, os, sys

vars = {}

def getNJExp(mexp, ident):
    exp_type = mexp.get_type()
    if exp_type == 'var':
        return getNJVar(mexp, ident)
    else:
        return getNJPred(mexp)

def getNJPred(mexp):
    exp_type = mexp.get_type()

    children = []
    for i in range(mexp.get_arity()):
        children.append(getNJExp(mexp.get_child(i), 0))

    if exp_type == "le":
        return children[0] <= mexp.get_parameter(0)
    elif exp_type == "ge":
        return children[0] >= mexp.get_parameter(0)
    elif exp_type == "sum":
        weights = []
        for i in range(mexp.get_arity()):
            weights.append(mexp.get_parameter(i))
        return Sum(children, weights)
    elif exp_type == "minimise":
        return Minimise(children[0])
    else:
        print "Error: Failed to parse expression:", type
        exit(1)

def parseModel(gmplfile, datafile=None):
    vars = {}
    parser = Osi.Solver(Model())
    parser.load_gmpl(gmplfile, datafile)
    prse = parser.solver

    model = Model()

    n = prse.num_expression()
    for i in range(n):
        expr = getNJExp(prse.get_expression(i), 0)
        model.add(expr)
    return model

def getNJVar(mexp, ident):
    ident = mexp.getVariableId();
    if not vars.has_key(ident):
        if mexp.is_continuous():
            vars[ident] = Variable(mexp.get_min(),
                                   mexp.get_max())
        else:
            vars[ident] = Variable(int(mexp.get_min()),
                                   int(mexp.get_max()))
    return vars[ident]

if __name__ == "__main__":
    from OsiClp import Solver
    gmplfile = None
    datafile = None
    if len(sys.argv) >= 2:
        gmplfile = sys.argv[1]
        if len(sys.argv) == 3:
            datafile = sys.argv[-1]
        model = parseModel(gmplfile, datafile)

        print model
        s = Solver(model)
        print s.solve()
        print model.variables
        if len(model.get_exprs()) == 0:
            print "ERROR!!!!!!!!!!!!!!!!"
