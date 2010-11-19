from Numberjack import *
import Mistral, MiniSat

def get_model():
    # Create the model
    model = Model()

    #Variable definitions, we need a variable for every letter
    s,m = (Variable(1,9) for val in range(2)) # These can't be zero as 
                                              # they are the start of a word
    e,n,d,o,r,y = (Variable(0,9) for val in range(6)) # These can

    model.add(      s*1000 + e*100 + n*10 + d +
                    m*1000 + o*100 + r*10 + e ==
          m*10000 + o*1000 + n*100 + e*10 + y)

    # Post the all different constraint on all the variables
    model.add(AllDiff((s,e,n,d,m,o,r,y)))

    return s,e,n,d,m,o,r,y,model


def solve(param):
    s,e,n,d,m,o,r,y,model = get_model()

    # Load up model into solver
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])

    # Now Solve
    out = ''
    if solver.solve():
        out += "     " + str(s) + ' ' + str(e) + ' ' + str(n) + ' ' + str(d) + ' ' +  '\n'
        out += " +  " + ' ' + str(m) + ' ' + str(o) + ' ' + str(r) + ' ' + str(e) + ' ' +  '\n'
        out +=  "= " + ' ' + str(m) + ' ' + str(o) + ' ' + str(n) + ' ' + str(e) + ' ' + str(y) + ' ' +  '\n'
    return out


solvers = ['Mistral']
default = {'solver':'Mistral', 'verbose':1, 'tcutoff':3}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)

