from Numberjack import *

def get_model(data):
    orders = Matrix(data.get("Orders"), data.get("Orders"))
    colours = Matrix(data.get("Colours"), data.get("Orders"))
    slabs = VarArray(data.get("Orders"), data.get("SlabSizes"))
    
    model = Model(
        # Objective
        Minimise(Sum(slabs)),
        # Maintain slab capacity
        [Sum(row, data.get("OrderWeights")) <= slab_var for (row, slab_var) in zip(orders, slabs)],
        # Make sure all orders are filles
        [Sum(col) == 1 for col in orders.col],
        # Make sure that at most p colours are assigned to each slab
        [Sum(row) <= data.get("p") for row in colours],
        # Channel between orders and colours
        [[order_var <= row_colour[colour] for (order_var, colour) in zip(col_order, data.get("OrderColours"))]
            for (col_order, row_colour) in zip(orders, colours.col)],
    )

    print model

    return (orders, colours, slabs, model)

def solve(param):
    data = SteelMillData()

    (orders, colours, slabs, model) = get_model(data)
    if param['solver'] == 'Mistral':
        # Break symmetries
        model.add([slabs[i] <= slabs[i+1] for i in range(0, data.get("Orders")-1)])

    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])
    solver.solve()
    print "Orders:\n", orders
    print "Colours:\n", colours
    print "Slabs:\n", slabs
    print 'Nodes:', solver.getNodes(), ' Time:', solver.getTime()

class SteelMillData:
    
    def __init__(self):
        self.SlabSizes = [0,1,3,4]
        self.Orders = 9
        self.Colours = 6
        self.OrderWeights = [2,3,1,1,1,1,1,2,1]
        self.OrderColours = [0,1,1,3,4,4,4,5,5]
        self.p = 2
    
    def get(self, name):
        if hasattr(self, name):
            return getattr(self, name)
        print "No Such Data!"
        return None


solvers = ['Mistral', 'MiniSat', 'SCIP', 'Walksat']
default = {'solver':'Mistral', 'verbose':1, 'tcutoff':3}


if __name__ == '__main__':
    param = input(default) 
    solve(param)

