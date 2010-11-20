from Numberjack import *
import Mistral, os, sys

vars = {}

def getNJExp(mexp, ident):
    type = mexp.get_type()
    if type == 'var':
        return getNJVar(mexp, ident)
    else:
        return getNJPred(mexp)
            
def getNJPred(mexp):
    type = mexp.get_type()
    
    children = []
    for i in range(mexp.get_arity()):
        children.append(getNJExp(mexp.get_child(i), 0))
        
    if type == "eq":
        return children[0] == children[1]
    elif type == "sum":
        return Sum(children)
    elif type == "prec": # I think this it an le
        return children[0] <= children[1]
    elif type == "mul":
        # Check this for SCIP
        return children[0] * children[1]
    elif type == "ne":
        return children[0] != children[1]
    elif type == "or":
        return children[0] | children[1]
    elif type == "abs":
        return Max([children[0], -children[0]])
    elif type == "max":
        return Max(children)
    elif type == "min":
        return Min(children)
    elif type == "neg":
        return -children[0]
    else:
        print "Error: Failed to parse expression:", type
        exit(1)
        
def parseModel(xmlfile):
    vars = {}
    parser = Mistral.Solver(Model())
    parser.load_xml(xmlfile, 0)
    prse = parser.solver

    model = Model()

    n = prse.num_expression()
    for i in range(n):
        expr = getNJExp(prse.get_expression(i), 0)
        model.add(expr)
        
    return model
    
def getNJVar(mexp, ident):
    vars[mexp.getVariableId()] = Variable(mexp.get_min(),
                                          mexp.get_max())
    return vars[mexp.getVariableId()]

if __name__ == "__main__":
    xmlfile = sys.argv[-1]
    if xmlfile[-3:] != 'xml':
        print "argument's probably wrong, mate"
    else:
        #print "Parseing Model", xmlfile
        model = parseModel(xmlfile)
        print model
        if len(model.get_exprs()) == 0:
            print "ERROR!!!!!!!!!!!!!!!!"
