from Numberjack import *

def get_model(N):
    queens = [Variable(N) for i in range(N)]
    model  = Model( 
        AllDiff( queens ),
        AllDiff( [queens[i] + i for i in range(N)] ),
        AllDiff( [queens[i] - i for i in range(N)] ) 
        )
    return (queens,model)

def solve(param):
    (queens,model) = get_model(param['N'])
    solver = model.load(param['solver'])
    solver.setHeuristic(param['heuristic'], param['value'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])
    solver.solve()

    out = ''
    if solver.is_sat() and param['print'] == 'yes':
        out += print_chessboard(queens)
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out    

def print_chessboard(queens):
    out = '+---'*len(queens)+'+\n'
    for queen in queens:
        out += ('|   '*queen.get_value()+'| Q |'+'   |'*(len(queens)-1-queen.get_value())+'\n'+'+---'*len(queens)+'+\n')
    return out

solvers = ['Mistral', 'MiniSat', 'SCIP', 'Walksat']
default = {'solver':'Mistral', 'N':6, 'heuristic':'MinDomainMinVal',
           'print':'yes', 'value':'Lex', 'verbose':1, 'tcutoff':3}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)
