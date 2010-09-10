#! /usr/bin/env python

from Numberjack import *

def model_magic_square(N):
    sum_val = N*(N*N+1)/2 # The magic number
    
    square = Matrix(N,N,1,N*N)
    water = Matrix(N,N,1,N*N) # water[a][b] stands for the water level in square a,b

    # Sum of water depth
    objective = (Sum(water.flat) - (N*N*(N*N+1)/2))

    model = Model(
        AllDiff( square.flat ),
        
        [Sum(row) == sum_val for row in square.row],
        [Sum(col) == sum_val for col in square.col],
        
        Sum([square[a][a] for a in range(N)]) == sum_val,
        Sum([square[a][N-a-1] for a in range(N)]) == sum_val,

        # symmetry breaking
        square[0][0] < square[0][N-1],
        square[0][0] < square[N-1][N-1],
        square[0][N-1] < square[N-1][0],

        # first, the rim
        water.row[0] == square.row[0],
        water.row[N-1] == square.row[N-1],
        water.col[0] == square.col[0],
        water.col[N-1] == square.col[N-1],
        # then the inner cells (max between their own height and of the water level around)
        [water[a][b] == Max((square[a][b], 
                            Min((water[a-1][b], water[a][b-1], 
                                 water[a+1][b], water[a][b+1])))
                            ) for a in range(1,N-1) for b in range(1,N-1)],

        # add the objective
        Maximise(objective)
        )

    return (square, water, objective, model)



    
def basic_solve(param):

    N = param['N']
    square,water,objective,model = model_magic_square(N)

    #solver = Mistral.Solver(model)
    solver = model.load(param['solver'])
    solver.setRandomSeed(11041979)

    solver.setVerbosity(param['verbose'])
    solver.setHeuristic(param['var'], param['val'], abs(param['rand']))
    solver.setNodeLimit(param['cutoff'])
    solver.setTimeLimit(param['tcutoff'])

    solver.solveAndRestart(GEOMETRIC, param['base'], param['factor'], param['decay']);
        
    #print 'Objective:', objective.get_value(), 'Nodes:', solver.getNodes(), ' Time:', solver.getTime()
    return (Solution(square),Solution(water),solver)
    

def print_water_level(water_level):
    import pylab

    N = len(water_level)
    cmap = pylab.cm.get_cmap('jet', N*N) 
    im1 = pylab.imshow(water_level, cmap=cmap, interpolation='nearest')
    pylab.colorbar()
    pylab.show()





def print_water(squares, waters):
    M = len(squares)
    N = len(squares[0])

    import pylab

    water_quantity = [[[waters[k][i][j] - squares[k][i][j] for j in range(N)] for i in range(N)] for k in range(M)]
    water_level = waters

    labels = ['Solid']
    labels.extend(['Water ('+str(i+1)+')' for i in range(M)])
    
    fracs = [(N*N)*(N*N+1)/2]
    fracs.extend([sum([sum(row) for row in waterq]) for waterq in water_quantity])

    explode=[0]
    explode.extend([0.1 for i in range(M)])


    #cmap = pylab.cm.get_cmap('RdYlBu',N*N) 
    cmap = pylab.cm.get_cmap('Blues',N*N) 


    for m in range(M):
        pylab.subplot(M+1,2,2*m+1)
        im1 = pylab.imshow(water_level[m], cmap=cmap, interpolation='nearest', origin='lower')
        if m == 0: pylab.title('Level', bbox={'facecolor':'0.9', 'pad':5})

        pylab.subplot(M+1,2,2*m+2)
        im2 = pylab.imshow(water_quantity[m], cmap=cmap, interpolation='nearest', origin='lower')
        if m == 0: pylab.title('Depth', bbox={'facecolor':'0.9', 'pad':5})

    pylab.subplot(M+1,2,2*M+1)
    pylab.pie(fracs, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)

    pylab.subplot(M+1,2,2*M+2)
    pylab.colorbar()

    pylab.axis('off')
    pylab.show()



def get_wall(square):
    N = len(square.row) ## order of the magic square
    No3 = int(float(N)/3) ## floor of N/3

    wall = []
    wall.extend(square[0][No3:N-No3]) ## middle third of the first row
    wall.extend(square[N-1][No3:N-No3]) ## middle third of the last row
    wall.extend(square.col[0][No3:N-No3]) ## middle third of the first column
    wall.extend(square.col[N-1][No3:N-No3]) ## middle third of the last column
    wall.extend([square[1+i][N-No3+i] for i in range(No3-1)]) ## small diagonal
    wall.extend([square[N-No3+i][N-2-i] for i in range(No3-1)]) ## small diagonal
    wall.extend([square[N-2-i][No3-1-i] for i in range(No3-1)]) ## small diagonal
    wall.extend([square[No3-1-i][1+i] for i in range(No3-1)]) ## small diagonal

    return wall


def solve_wall(param):

    N = param['N']
    square,water,objective,model = model_magic_square(N)
    wall = get_wall(square)

    solver = model.load(param['solver'])
    #solver = Mistral.Solver(model)
    solver.setRandomSeed(11041979)

    solver.setVerbosity(param['verbose'])
    solver.setHeuristic(param['var'], param['val'], abs(param['rand']))
    solver.setTimeLimit(param['tcutoff'])

    square_sol = None
    water_sol = None

    water_level = 0
    margin = 0
    while margin <= 1:
        #print '\n\nsolve with margin', margin
        solver.setNodeLimit(param['cutoff']/2)

        ## save 
        solver.save()

        ## post the wall
        for bric in wall: solver.post(bric > N*N-(len(wall)+margin))

        ## implied constraint
        for cell in (set(square.flat) - set(wall)): solver.post(cell <= (N*N-(len(wall)-margin)))

        ## solve
        solver.solveAndRestart(GEOMETRIC, param['base'], param['factor'], param['decay']);
        print solver.is_sat()
        if solver.is_sat() and objective.get_value() > water_level:
            water_level = objective.get_value()
            square_sol = Solution(square)
            water_sol = Solution(water)

        ## relax the wall
        solver.reset()        
        solver.undo()
        margin += 1

    #print 'Objective:', (water_level), 'Nodes:', solver.getNodes(), ' Time:', solver.getTime()
    return (square_sol,water_sol,solver)




def solve(param):
    square,water,solver = None,None,None
    if param['algo'] == 'basic':
        square,water,solver = basic_solve(param)
        if param['print'] == 'yes': print_water([square],[water])
    elif param['algo'] == 'wall':
        square,water,solver = solve_wall(param)
        if param['print'] == 'yes': print_water([square],[water])
    else:
        square1,water1,solver = basic_solve(param)
        square2,water2,solver = solve_wall(param)
        if param['print'] == 'yes': print_water([square1,square2],[water1,water2])
    out = ''
    if solver.is_sat():
        out = str(square)+'\n\n'+str(water)+'\n'
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out   


solvers = ['Mistral']
default = {'N':4, 'var':'DomainOverWDegree', 'proba':0.8,
           'val':'RandomSplit', 'restart':'yes', 'rand':5, 
           'verbose':1, 'cutoff':30000, 'factor':1.2, 'base':64, 
           'decay':0.0,'algo':'basic','solver':'Mistral', 
           'print':'no','tcutoff':3}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)
