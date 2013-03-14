from Numberjack import *


def getNJExp(variables, mexp, ident):
    t = mexp.get_type()
    if t == 'var':
        return getNJVar(variables, mexp, ident)
    else:
        return getNJPred(variables, mexp)


def getNJPred(variables, mexp):
    t = mexp.get_type()

    children = []
    for i in range(mexp.get_arity()):
        children.append(getNJExp(variables, mexp.get_child(i), 0))

    if t == "eq":
        return children[0] == children[1]
    elif t == "sum":
        return Sum(children)
    elif t == "prec":  # I think this it an le
        return children[0] <= children[1]
    elif t == "mul":
        # Check this for SCIP
        return children[0] * children[1]
    elif t == "ne":
        return children[0] != children[1]
    elif t == "or":
        return children[0] | children[1]
    elif t == "abs":
        return Max([children[0], -children[0]])
    elif t == "max":
        return Max(children)
    elif t == "min":
        return Min(children)
    elif t == "neg":
        return -children[0]
    else:
        print "Error: Failed to parse expression:", t
        exit(1)


def parseModel(xmlfile):
    import Mistral

    variables = {}
    parser = Mistral.Solver(Model())
    parser.load_xml(xmlfile, 0)
    prse = parser.solver
    model = Model()

    n = prse.num_expression()
    for i in range(n):
        expr = getNJExp(variables, prse.get_expression(i), 0)
        model.add(expr)

    return model


def getNJVar(variables, mexp, ident):
    variables[mexp.getVariableId()] = Variable(mexp.get_min(), mexp.get_max())
    return variables[mexp.getVariableId()]


if __name__ == "__main__":
    import sys

    xmlfile = sys.argv[-1]
    if xmlfile[-3:] != 'xml':
        print "argument's probably wrong, mate"
    else:
        #print "Parseing Model", xmlfile
        model = parseModel(xmlfile)
        print model
        if len(model.get_exprs()) == 0:
            print "ERROR!!!!!!!!!!!!!!!!"
