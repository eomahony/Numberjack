from Numberjack import *

def costas_array(N):
    # Create the variables
    sequence = VarArray(N,1,N)
    
    # State the model
    model = Model(
        AllDiff(sequence),
        [AllDiff([sequence[j] - sequence[j+i+1] for j in range(N-i-1)]) for i in range(N-2)] 
        )
    
    return sequence,model
    
def printCostasTriangle(sequence):
    N = len(sequence)
    print ''.join([str(int(var.get_value())).rjust(3) for var in sequence])
    for i in range(N-1):
        print ''.join([str(int(sequence[j].get_value())-int(sequence[j+i+1].get_value())).rjust(3) 
                       for j in range(N-i-1)])



param = input({'solver':'Mistral', 'N':6, 'verbosity':1})

sequence,model = costas_array(param['N'])
solver = model.load(param['solver'])
solver.setVerbosity(param['verbosity'])
solver.solve()
    
printCostasTriangle(sequence)
print '\n\t', solver.getNodes(), 'nds  ', solver.getTime(), 's'

