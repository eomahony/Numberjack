#! /usr/bin/env python

import sys
import subprocess
from sets import Set

all_paths = sys.path
all_paths.append('.')
#install_paths = Set([])


files = Set(['Mistral', 'mistral_lib', 'SCiP', 'SolverExample', 'Numberjack', 'MiniSat', 'Gecode', 'Casper', 'ThreadSolver', 'BenchMarks', 'Checker', 'Constraints', 'Decomp', '__init__', 'Modelling', 'setup_so', 'SolverInterface', 'SolverTests', 'Solving', 'Variables', 'WrapperTests', 'XCSPOut', 'XCSPOutTest'])

for path in all_paths:    
    try:
        #print 'ls', path
        ls = subprocess.Popen(['ls', '-1', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ls.wait()
        for pylib in ls.stdout:
            for mod in files:
                if pylib.find(mod) != -1:
                    #print pylib, pylib[:1], pylib[-4:]
                    if ((pylib[:1] == '_' and pylib[-4:] == '.so\n') or pylib.find('.py') != -1 or pylib.find('egg') != -1):
                        print 'rm', '-rf', path+'/'+pylib[:-1]
                        rm = subprocess.Popen(['rm', '-rf', path+'/'+pylib[:-1]])
                        rm.wait()
    except:
        continue

print 'rm', '-rf', './setup.py'
rm = subprocess.Popen(['rm', '-rf', './setup.py'])
rm.wait()
