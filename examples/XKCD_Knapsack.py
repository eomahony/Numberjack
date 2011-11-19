#!/usr/bin/env python

"""
Knapsack problem from http://xkcd.com/287/
"""

from Numberjack import *


def get_model(param):
    num_prices = 6
    price = [215, 275, 335, 355, 420, 580]
    apitizers = ["Mixed Fruit", "French Fries", "Side Salad", "Hot Wings", "Mozzarella Sticks", "Sample Plate"]
    total = 1505
    max_qty = 100

    quantities = VarArray(num_prices, 0, max_qty, 'quantities')
    model = Model(
        Sum([quantities[i] * price[i] for i in xrange(num_prices)]) == total
    )

    return num_prices, price, quantities, apitizers, model


def solution_string(num_prices, price, quantities, apitizers):
    return "\n".join("%s x %s ($%.2lf)" % (quantities[i], apitizers[i], price[i] / 100.0) for i in xrange(num_prices))


def solve(param):
    num_prices, price, quantities, apitizers, model = get_model(param)
    solver = model.load(param['solver'], quantities)

    if param['solver'] == 'Mistral':
        solver.startNewSearch()
        while solver.getNextSolution() == SAT:
            print "\nSOLUTION:\n", solution_string(num_prices, price, quantities, apitizers)
    else:
        if solver.solve():
            print solution_string(num_prices, price, quantities, apitizers)
        else:
            print "No solution"


solvers = ['Mistral', 'SCIP']
default = {'solver': 'Mistral', 'verbose': 1, 'tcutoff': 3}

if __name__ == '__main__':
    param = input(default)
    solve(param)
