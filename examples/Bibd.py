from Numberjack import *

def model_bibd(v, b, r, k, l):
    matrix = Matrix(v,b)
    model = Model( 
        [Sum(row) == k for row in matrix.row], # every row adds up to k        
        [Sum(col) == r for col in matrix.col], # every column adds up to r
        
        # the scalar product of every pair of columns adds up to l
        [Sum([(row[col_i] & row[col_j]) for row in matrix.row]) == l  
         for col_i in range(v) for col_j in range(col_i)]
        )

    return (matrix,model)

def solve_bibd(param):
    (matrix,model) = model_bibd(*param('v','b','r','k','l')) 

    if param['solver'] == 'Mistral':
        model += [matrix.row[i] <= matrix.row[i+1] for i in range(param['v']-1)]
        model += [matrix.col[i] <= matrix.col[i+1] for i in range(param['b']-1)]

    solver = model.load(param['solver'])
    if solver.solve():
        print matrix
    print 'Nodes:', solver.getNodes(), ' Time:', solver.getTime()

solve_bibd(input({'solver':'Mistral','v':7,'b':7,'r':3,'k':3,'l':1}))

