#! /usr/bin/env python

import os
import threading
import subprocess
from sets import Set


solvers = Set([])
for line in subprocess.Popen( ['find', 'solvers', '-name', '*.slv'], stdout=subprocess.PIPE ).stdout:
    solvers.add( line[8:-5] )

oldsolvers = Set([])
for line in subprocess.Popen( ['find', 'results', '-name', '*.res'], stdout=subprocess.PIPE ).stdout:
    if line[8:-5] not in solvers:
        oldsolvers.add( line[8:-5] )

sat = Set([])
unsat = {}.fromkeys( oldsolvers )
path = {}.fromkeys( solvers )

problem = []

for solver in solvers:
    path[solver] = {}
    for line in open( 'results/'+solver+'.pth' ):
        path[solver][line.split()[0]] = line.split()[1]

for solver in oldsolvers:
    unsat[solver] = Set([])
    for line in open( 'results/'+solver+'.res' ):        
        if line.split()[1] == 'SATISFIABLE':
            sat.add( line.split()[0] )
        elif line.split()[1] == 'UNSATISFIABLE':
            unsat[solver].add( line.split()[0] )

def checkSolution( solution, bench ):
    process = subprocess.Popen( ['tools/CSP-XML-parser/tools/CSP-verifier', bench], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    if process.communicate( solution )[0] == 'OK\n':
        return True
    return False

def checkOldProofs( bench ):
    print '\t +CHECK', bench.ljust(49),
    problem = False
    for solver in oldsolvers:
        if bench in unsat[solver]:
            print ' /!\\ WARNING, the old proof by', solver, 'was wrong /!\\',
            problem.append( 'the old proof by '+solver+' was wrong' )
            problem = True
    if problem == False:
        print 'OK'
    else:
        print ''


for solver in solvers:
    print solver
    solutions = {}
    for line in open( 'results/'+solver+'.sol' ):
        solutions[line.split()[0]] = 's SATISFIABLE\nv '+line[line.find(' ')+1:]
    for line in open( 'results/'+solver+'.res' ):
        bench = line.split()[0]
        outcome = line.split()[1]
        if outcome == 'SATISFIABLE':
            print '\t+CHECK', bench.ljust(50),
            if checkSolution( solutions[bench], path[solver][bench] ):
                print 'OK'
                checkOldProofs( bench )
            else:
                problem.append( solver+' gives a wrong certificate for '+bench )
                print ' /!\\ WARNING,', solver, 'gives a wrong certificate for', bench, ' /!\\'
        else:
            print '\t-CHECK', bench.ljust(50),
            if bench in sat:
                problem.append( solver+' gives a wrong proof for '+bench )
                print ' /!\\ WARNING,', solver, 'gives a wrong proof for', bench, ' /!\\'
            else:
                print 'OK'

print problem







