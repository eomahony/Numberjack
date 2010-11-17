#! /usr/bin/env python

try:
    import sys

    solver = open(sys.argv[1], 'a')

    sol_name = sys.argv[1][:-3].split("/")[-1]
    wrap_name = sol_name
    if len(sys.argv) > 2: wrap_name = sys.argv[2]

    solver_code = '''
    
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1):
        Numberjack.NBJ_STD_Solver.__init__(self, "%s", "%s", model, X, FD,
            clause_limit)
    
                ''' % (sol_name,wrap_name)

    #for line in footer:
    #    solver.write(line)
    solver.write(solver_code)
    solver.close()


except:
    print 'abort'
    pass
