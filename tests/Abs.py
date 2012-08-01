#!/usr/bin/env python

from Numberjack import *

def test_abs():
    x = Variable(-10, 10)
    model = Model(Abs(x) < 3)
    
    solver = model.load("Mistral")
    solver.startNewSearch();
    
    return " ".join(str(x) for s in solver.solutions() if solver.is_sat())
        

if __name__ == "__main__":
    print "Solutions:", test_abs()