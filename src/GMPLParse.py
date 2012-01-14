from Numberjack import *
import OsiGlpk as Osi, os, sys
import OsiGlpk

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
        
def parseModel(gmplfile):
    vars = {}
    parser = Osi.Solver(Model())
    parser.load_gmpl(gmplfile)
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
        vars[ident] = Variable(mexp.get_min(),
                               mexp.get_max())
    return vars[ident]

if __name__ == "__main__":
    gmplfile = sys.argv[-1]
    if gmplfile[-3:] != 'mod':
        print "argument's probably wrong, mate"
    else:
        #print "Parseing Model", gmplfile
        model = parseModel(gmplfile)
        print model
        s = Osi.Solver(model)
        print s.solve()
        print model.variables
        
        if len(model.get_exprs()) == 0:
            print "ERROR!!!!!!!!!!!!!!!!"
