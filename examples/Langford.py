from Numberjack import *

def get_model(M,N):
    X = [Variable(1, N*M-(i+2)*(M-1)) for i in range(N)]
    
    model = Model(
        AllDiff( [X[i]+((i+2)*j) for j in range(M) for i in range(N)] ),
        X[0] > X[1] # Break symmetry
    )
    return (X,model)

def solve(param):
    (X,model) = get_model(param['M'],param['N'])
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])

    out = '' 
    if param['solver'] == 'Mistral':
        solver.startNewSearch();
        langford_number = 0
        while solver.getNextSolution() == SAT:
            out += (printLangford(param['M'],X)+'\n')
            langford_number += 1
        out += ('L('+str(param['M'])+','+str(param['N'])+') = '+str(langford_number)+'\n')
    else:
        if solver.solve(): out += (printLangford(param['M'],X)+'\n')
        else: out += 'No solution\n'
    out += ('Nodes: ' + str(solver.getNodes()))
    return out    


def printLangford(M,X):
    N = len(X)
    sequence = [0]*(M*N)
    for i in range(N):
        for j in range(M):
            sequence[X[i].get_value()+j*(i+2)-1] = (i+1)
    return str(sequence)


solvers = ['Mistral', 'MiniSat', 'SCIP', 'Walksat']
default = {'solver':'Mistral', 'N':7, 'M':2, 'verbose':1, 'tcutoff':3}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)

