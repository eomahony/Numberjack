from __future__ import print_function
from Numberjack import *
import pickle
import subprocess


ls = subprocess.Popen(['ls', '-1'], stdout=subprocess.PIPE)

param = input({'update':'no', 'solver':'Mistral'})

outcome = {}

reference = {}

try:
    reference = pickle.load(open('.unit_test_reference.pkl'))
except:
    pass

for example in ls.stdout:
    if example[-4:-1] == '.py' and example[:-4] != 'unit_tests':

        lib_str = example[:-4]
        outcome[lib_str] = {}

        print('\n**************************************')
        print('******', lib_str)
        print('**************************************')
        lib = __import__(lib_str)
        lib_param = lib.default
        lib_param['verbose'] = 0
        for solver in lib.solvers:
            print('run', lib_str, 'with', solver)
            lib_param['solver'] = solver
            outcome[lib_str][solver] = lib.solve(lib_param)

for lib in list(reference.keys()):
    print('checking', lib+':')
    for solver in list(reference[lib].keys()):
        print('\t ->', solver, end=' ')
        if lib not in outcome:
            print('Warning: no record for', lib)
        else:
            if solver not in outcome[lib]:
                print('Warning: no record for', solver, 'on', lib)
            else:
                if outcome[lib][solver] != reference[lib][solver]:
                    print('Warning: discrepancy on ', lib, 'with', solver+':')
                    print('was:')
                    print(reference[lib][solver])
                    print('now:')
                    print(outcome[lib][solver])
                    if param['update'] == 'yes':
                        print('... saving')
                        reference[lib][solver] = outcome[lib][solver]
                else:
                    print('ok')

for lib in list(outcome.keys()):
    for solver in list(outcome[lib].keys()):
        if lib not in reference:
            print(lib, 'is a new problem... saving')
            reference[lib] = outcome[lib]
        else:
            if solver not in reference[lib]:
                print(solver, 'did not previously run on', lib+'... saving')
                reference[lib][solver] = outcome[lib][solver]


pickle.dump(reference, open('.unit_test_reference.pkl', 'w'))






