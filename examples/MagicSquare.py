from Numberjack import *

def get_model(N):
    sum_val = N*(N*N+1)/2 # This is what all the columns, rows and diagonals must add up tp

    square = Matrix(N,N,1,N*N)
    model = Model(
            AllDiff( square.flat ),
        
            [Sum(row) == sum_val for row in square.row],
            [Sum(col) == sum_val for col in square.col],
            
            Sum([square[a,a] for a in range(N)]) == sum_val,
            Sum([square[a,N-a-1] for a in range(N)]) == sum_val
            )
    return (square,model)
    
def solve(param):
    (square,model) = get_model(param['N'])
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setHeuristic(param['var'], param['val'], param['rand'])
    solver.setTimeLimit(param['cutoff'])

    if param['restart'] == 'yes':
        solver.solveAndRestart();
    else:
        solver.solve()

    out = ''
    if solver.is_sat():
        out = str(square)
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out    

solvers = ['Mistral', 'MiniSat', 'SCIP', 'Walksat']
default = {'solver':'Mistral', 'N':4, 'var':'MinDomain',
           'val':'RandomMinMax', 'restart':'yes', 'rand':2, 'verbose':1, 'cutoff':10}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)
