#! /usr/bin/env python

# Warehouse Location Problem
#
# Given a set of existing shops and a candidate set of warehouses to be opened,
# the warehouse location problem consists of choosing which warehouses to be
# opened and consequently the respective shops which each one will supply. There
# is a cost associated with opening each warehouse, as well as a supply cost for
# each warehouse-shop supply pair. The objective is to minimise the total cost
# of warehouse operations and supply costs.
#
# CSPLib Problem 034 - http://www.csplib.org/Problems/prob034/

from __future__ import print_function
from Numberjack import *
import re


def solve(param):
    data = WareHouseParser(param['data'])

    WareHouseOpen = VarArray(data.NumberOfWarehouses)
    ShopSupplied = Matrix(data.NumberOfShops, data.NumberOfWarehouses)

    # Cost of running warehouses
    warehouseCost = Sum(WareHouseOpen, data.WareHouseCosts)

    # Cost of shops using warehouses
    transpCost = Sum([Sum(varRow, costRow) for (varRow, costRow) in zip(ShopSupplied, data.SupplyCost)])

    obj = warehouseCost + transpCost

    model = Model(
        # Objective function
        Minimise(obj),
        # Channel from store opening to store supply matrix
        [[var <= store for var in col] for (col, store) in zip(ShopSupplied.col, WareHouseOpen)],
        # Make sure every shop if supplied by one store
        [Sum(row) == 1 for row in ShopSupplied.row],
        # Make sure that each store does not exceed it's supply capacity
        [Sum(col) <= cap for (col, cap) in zip(ShopSupplied.col, data.Capacity)]
    )

    solver = model.load(param['solver'])
    # solver.setNodeLimit(param['cutoff'])
    solver.setHeuristic('DomainOverWDegree', 'Guided')
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])
    solver.solveAndRestart()

    if solver.is_sat():
        if solver.is_opt():
            print("Optimal")

        print("Total cost: ", str(obj.get_value()))
        print("Nodes:", solver.getNodes())
        print("SolveTime:", solver.getTime())

    elif solver.is_unsat():
        print("Unsatisfiable")
    else:
        print("Unknown")


class WareHouseParser(object):
    "Parses and stores the data for a warehouse location problem instance."

    def __init__(self, filename):
        with open(filename, "rt") as f:
            alltext = f.read()
            matchnbw = re.search(r"NbW=(?P<NumberOfWarehouses>\d+);", alltext)
            matchnbs = re.search(r"NbS=(?P<NumberOfShops>\d+);", alltext)
            matchfixed = re.search(r"fixed=(?P<FixedCost>\d+);", alltext)
            self.NumberOfWarehouses = int(matchnbw.groupdict()["NumberOfWarehouses"])
            self.NumberOfShops = int(matchnbs.groupdict()["NumberOfShops"])
            self.FixedCost = int(matchfixed.groupdict()["FixedCost"])

            self.SupplyCost = []
            matchsupply = re.search(r"SupplyCost=\[(?P<supplystr>.*)\];", alltext, re.MULTILINE | re.DOTALL)
            supplylines = matchsupply.groupdict()["supplystr"].strip().split("\n")
            for supplyline in supplylines:
                costs = map(int, re.findall(r"\d+", supplyline))
                self.SupplyCost.append(costs)

        self.WareHouseCosts = [self.FixedCost for val in range(self.NumberOfWarehouses)]
        # # There was a fixed capacity for all the warehouse problems
        self.Capacity = [4 for val in range(self.NumberOfWarehouses)]


default = {'solver': 'Mistral', 'data': 'data/cap44.dat', 'cutoff': 50000, 'verbose': 1, 'tcutoff': 30}

if __name__ == '__main__':
    param = input(default)
    solve(param)
