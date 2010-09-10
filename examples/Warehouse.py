#! /usr/bin/env python

from Numberjack import *


def solve(param):
    data = WareHouseParser(param['data'])
    lib = __import__(param['solver'])
    cutoff = param['cutoff']

    WareHouseOpen = VarArray(data.get("NumberOfWarehouses"))
    
    ShopSupplied = Matrix(data.get("NumberOfShops"),
                          data.get("NumberOfWarehouses"))
    
    # Cost of running warehouses
    warehouseCost = Sum(WareHouseOpen, data.get("WareHouseCosts"))
    
    # Cost of shops using warehouses
    transpCost = Sum([ Sum(varRow, costRow) for (varRow, costRow) in zip(ShopSupplied, data.get("SupplyCost"))])
    
    obj = warehouseCost + transpCost
    
    model = Model(
        # Objective function
        Minimise(obj), 
        # Channel from store opening to store supply matrix
        [[var <= store for var in col] for (col, store) in zip(ShopSupplied.col, WareHouseOpen)],
        # Make sure every shop if supplied by one store
        [Sum(row) == 1 for row in ShopSupplied.row],
        # Make sure that each store does not exceed it's supply capacity
        [Sum(col) <= cap for (col, cap) in zip(ShopSupplied.col, data.get("Capacity"))]
    )

    solver = lib.Solver(model)
    solver.setNodeLimit(cutoff)
    solver.setHeuristic('Impact')
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])
    solver.solve()

    return "\nFINAL COST: " + str(obj.get_value())
    
class WareHouseParser:
    
    def __init__(self, file):
        lines = open(file).readlines() 
        self.NumberOfWarehouses = int(lines[0][4:-2])
        self.NumberOfShops = int(lines[1][4:-2])
        self.FixedCost = int(lines[2][6:-2])
        self.WareHouseCosts = [self.FixedCost for val in range(self.NumberOfWarehouses)]
        self.SupplyCost = eval("["+
                               " ".join((map((lambda a: a[:-1]+","),
                                                lines[4:-1]))) +
                               "]")
        # There was a fixed capacity for all the warehouse problems
        self.Capacity = [4 for val in range(self.NumberOfWarehouses)]
        
    def get(self, name):
        if hasattr(self, name):
            return getattr(self, name)
        print name + " \t No Such Data!"
        return None
    

solvers = ['Mistral', 'SCIP']
default = {'solver':'Mistral', 'data':'data/cap44.dat.txt', 'cutoff':5000, 'verbose':1, 'tcutoff':3}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)
