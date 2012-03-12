from Numberjack import *
import os
import sys
import OsiGlpk as Osi

vars = {}


def getNJExp(mexp, ident):
    exp_type = mexp.get_type()
    if exp_type == 'var':
        return getNJVar(mexp, ident)
    else:
        return getNJPred(mexp)

def ifwhole(val):
    return val if val % 1 != 0 else int(val)

def getNJPred(mexp):
    exp_type = mexp.get_type()

    children = []
    for i in range(mexp.get_arity()):
        children.append(getNJExp(mexp.get_child(i), 0))

    if exp_type == "le":
        return children[0] <= ifwhole(mexp.get_parameter(0))
    elif exp_type == "ge":
        return children[0] >= ifwhole(mexp.get_parameter(0))
    elif exp_type == "sum":
        weights = []
        for i in range(mexp.get_arity()):
            weights.append(ifwhole(mexp.get_parameter(i)))
        return Sum(children, weights)
    elif exp_type == "minimise":
        return Minimise(children[0])
    else:
        print "Error: Failed to parse expression:", type
        exit(1)


def parseModel(filename, filetype='mod', datafile=None):
    print (filename, filetype, datafile)
    vars = {}
    parser = Osi.Solver(Model())
    if filetype == 'mod':
        parser.load_gmpl(filename, datafile)
    elif filetype == 'mps':
        parser.load_mps(filename, 'mps')
    elif filetype == 'lp':
        parser.load_lp(filename, 0)
    else:
        print "Unknown filetype, exiting"
        return None
    prse = parser.solver

    model = Model()

    n = prse.num_expression()
    for i in range(n):
        expr = getNJExp(prse.get_expression(i), 0)
        model.add(expr)
    return model


def getNJVar(mexp, ident):
    ident = mexp.getVariableId()
    if ident not in vars:
        if mexp.is_continuous():
            vars[ident] = Variable(mexp.get_min(),
                                   mexp.get_max(),
                                   mexp.name())
        else:
            vars[ident] = Variable(int(mexp.get_min()),
                                   int(mexp.get_max()),
                                   mexp.name())
    return vars[ident]


if __name__ == "__main__":
    from Mistral import Solver
    if len(sys.argv) >= 2:
        datafile = None
        filename = sys.argv[1]
        parts = filename.split('.')
        extensions = parts[-2:len(parts)]
        # if file is compressed ignore
        if extensions[1] == 'gz':
            filetype = extensions[0]
        else:
            filetype = extensions[1]

        if len(sys.argv) == 3:
            datafile = sys.argv[-1]
        model = parseModel(filename, filetype, datafile)
        if len(model.get_exprs()) == 0:
            print "ERROR!!!!!!!!!!!!!!!!"
            exit()

        #print model
        s = Solver(model)
        print s.solve()
        print {v.name(): v.get_value() for v in model.variables}
