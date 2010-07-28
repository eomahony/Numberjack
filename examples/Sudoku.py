from Numberjack import *

param = input({'N':3,'solver':'Mistral','file':'data/sdk.txt'})
N = param['N']
clues = param['file']

grid = Matrix(N*N,N*N,1,N*N,'cell_')

sudoku = Model( [AllDiff(row) for row in grid.row],
                [AllDiff(col) for col in grid.col],
                [AllDiff(grid[x:x+N, y:y+N].flat) for x in range(0,N*N,N) for y in range(0,N*N,N)],
                [(x == int(v)) for (x,v) in zip(grid.flat, "".join(open(clues)).split() ) if v != '*']
               )

solver = sudoku.load(param['solver'])
solver.solve()
print solver.is_sat()
print grid

