#! /usr/bin/env python

try:
    import sys

    print sys.argv[1]

    solver = open(sys.argv[1], 'a')
    
    solver_code = '''
    
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False):
        Numberjack.NBJ_STD_Solver.__init__(self, "%s", model, X, FD)
    
                ''' % (sys.argv[1][:-3].split("/")[1])

    #for line in footer:
    #    solver.write(line)
    solver.write(solver_code)

except:
    print 'abort'
    pass
