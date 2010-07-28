from Numberjack import *

def model_langford_sequence(M,N):
    X = [Variable(1, N*M-(i+2)*(M-1)) for i in range(N)]
    
    model = Model(
        AllDiff( [X[i]+((i+2)*j) for j in range(M) for i in range(N)] ),
        X[0] > X[1] # Break symmetry
    )
    return (X,model)

def compute_langford_number(param):
    (X,model) = model_langford_sequence(param['M'],param['N'])
    solver = model.load(param['solver'])

    if param['solver'] == 'Mistral':
        solver.startNewSearch();
        langford_number = 0
        while solver.getNextSolution() == SAT:
            printLangford(param['M'],X)
            langford_number += 1
        print 'L('+str(param['M'])+','+str(param['N'])+') =', langford_number
    else:
        if solver.solve(): printLangford(param['M'],X)
        else: print 'No solution'

    print 'Nodes:', solver.getNodes(), ' Time:', solver.getTime()

def printLangford(M,X):
    N = len(X)
    sequence = [0]*(M*N)
    for i in range(N):
        for j in range(M):
            sequence[X[i].get_value()+j*(i+2)-1] = (i+1)
    print sequence

compute_langford_number(input({'solver':'Mistral', 'N':10, 'M':3}))

