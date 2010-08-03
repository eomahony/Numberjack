from Numberjack import *

def model_queens(N):
    queens = [Variable(N) for i in range(N)]
    model  = Model( 
        AllDiff( queens ),
        AllDiff( [queens[i] + i for i in range(N)] ),
        AllDiff( [queens[i] - i for i in range(N)] ) 
        )
    return (queens,model)

def solve_queens(param):
    (queens,model) = model_queens(param['N'])
    solver = model.load(param['solver'])
    solver.setHeuristic(param['heuristic'], param['value'])
    solver.solve()
    if solver.is_sat() and param['print'] == 'yes':
        print_chessboard(queens)
    print 'Nodes:', solver.getNodes(), ' Time:', solver.getTime()

def print_chessboard(queens):
    separator = '+---'*len(queens)+'+'
    for queen in queens:
        print separator
        print '|   '*queen.get_value()+'| Q |'+'   |'*(len(queens)-1-queen.get_value())
    print separator

solve_queens(input({'solver':'Mistral', 'N':10, 'heuristic':'MinDomainMinVal',
                    'print':'yes', 'value':'Lex'}))

