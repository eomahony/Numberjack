from Numberjack import *

def get_model(N):
    # Create the variables
    sequence = VarArray(N,1,N)
    
    # State the model
    model = Model(
        AllDiff(sequence),
        [AllDiff([sequence[j] - sequence[j+i+1] for j in range(N-i-1)]) for i in range(N-2)] 
        )
    
    return sequence,model
    
def printCostasTriangle(sequence):
    N = len(sequence)
    print ''.join([str(int(var.get_value())).rjust(3) for var in sequence])
    for i in range(N-1):
        print ''.join([str(int(sequence[j].get_value())-int(sequence[j+i+1].get_value())).rjust(3) 
                       for j in range(N-i-1)])

def solve(param):
    sequence,model = get_model(param['N'])
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbosity'])
    solver.setTimeLimit(3)
    solver.solve()
    
    if solver.is_sat():
        printCostasTriangle(sequence)
        print '\n\t', solver.getNodes(), 'nds  ', solver.getTime(), 's'    


solvers = ['Mistral', 'MiniSat', 'SCIP', 'Walksat']
default = {'solver':'Mistral', 'N':6, 'verbosity':1}

if __name__ == '__main__':
    param = input(default) 
    solve(param)


