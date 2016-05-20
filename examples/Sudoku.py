from Numberjack import *


def get_model(N, clues):
    grid = Matrix(N*N, N*N, 1, N*N)

    sudoku = Model([AllDiff(row) for row in grid.row],
                   [AllDiff(col) for col in grid.col],
                   [AllDiff(grid[x:x + N, y:y + N]) for x in range(0, N*N, N) for y in range(0, N * N, N)],
                   [(x == int(v)) for x, v in zip(grid.flat, "".join(open(clues)).split()) if v != '*']
                   )
    return grid, sudoku


def solve(param):
    N = param['N']
    clues = param['file']

    grid, sudoku = get_model(N, clues)

    solver = sudoku.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])

    solver.solve()
    out = ""
    if solver.is_sat():
        out = str(grid)

    return out


default = {'N': 3, 'solver': 'Mistral', 'file': 'data/sdk.txt', 'verbose': 0, 'tcutoff': 30}


if __name__ == '__main__':
    param = input(default)
    print(solve(param))
