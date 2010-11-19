from Numberjack import *
import MipWrapper

def get_model(N):
    # Create the variables
    sequence = VarArray(N,1,N)
    
    seqs = [[sequence[j] - sequence[j+i+1] for j in range(N-i-1)] for i in range(N-2)]
    
    # State the model
    model = Model(
        AllDiff(sequence),
        #[AllDiff([sequence[j] - sequence[j+i+1] for j in range(N-i-1)]) for i in range(N-2)] 
        [ AllDiff(seq) for seq in seqs ]
        )
    
    return sequence,model, seqs
    
def printCostasTriangle(sequence):
    N = len(sequence)
    out = ''.join([str(int(var.get_value())).rjust(3) for var in sequence])+'\n'
    for i in range(N-1):
        out += ''.join([str(int(sequence[j].get_value())-int(sequence[j+i+1].get_value())).rjust(3)
                        for j in range(N-i-1)])+'\n'
    return out

def solve(param):
    sequence,model, seqs = get_model(param['N'])
        
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])
    
    res = solver.solve()  

    out = ''
    if solver.is_sat():
        out += printCostasTriangle(sequence)
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out



solvers = ['Mistral', 'MiniSat', 'SCIP', 'Walksat']
default = {'solver':'Mistral', 'N':6, 'verbose':1, 'tcutoff':3}


if __name__ == '__main__':
    param = input(default) 
    print solve(param)


