from Numberjack import *

# N-Queens

# The N-Queens problem is the probelm of placing N queens on an N x N chess
# board such that no two queens are attacking each other. A queen is attacking
# another if it they are on the same row, same column, or same diagonal.


def get_model(N):
    queens = VarArray(N, N)
    model = Model(
        AllDiff(queens),
        AllDiff([queens[i] + i for i in range(N)]),
        AllDiff([queens[i] - i for i in range(N)])
    )
    return queens, model


def solve(param):
    queens, model = get_model(param['N'])
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
    out = '+---' * len(queens) + '+\n'
    for queen in queens:
        out += ('|   '*queen.get_value()+'| Q |'+'   |'*(len(queens)-1-queen.get_value())+'\n'+'+---'*len(queens)+'+\n')
    return out


default = {'solver': 'Mistral', 'N': 6, 'heuristic': 'MinDomainMinVal',
           'print': 'yes', 'value': 'Lex', 'verbose': 0, 'tcutoff': 30}


if __name__ == '__main__':
    param = input(default)
    print(solve(param))
