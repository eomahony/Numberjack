from Numberjack import *
import cPickle
import subprocess


ls = subprocess.Popen(['ls', '-1'], stdout=subprocess.PIPE)

param = input({'update':'no', 'solver':'Mistral'})

outcome = {}

reference = {}

try:
    reference = cPickle.load(open('.unit_test_reference.pkl'))
except:
    pass

for example in ls.stdout:
    if example[-4:-1] == '.py' and example[:-4] != 'unit_tests':

        lib_str = example[:-4]
        outcome[lib_str] = {}

        print '\n**************************************'
        print '******', lib_str
        print '**************************************'
        lib = __import__(lib_str)
        lib_param = lib.default
        lib_param['verbose'] = 0
        for solver in lib.solvers:
            print 'run', lib_str, 'with', solver
            lib_param['solver'] = solver
            outcome[lib_str][solver] = lib.solve(lib_param)

for lib in reference.keys():
    print 'checking', lib+':'
    for solver in reference[lib].keys():
        print '\t ->', solver,
        if not outcome.has_key(lib):
            print 'Warning: no record for', lib
        else:
            if not outcome[lib].has_key(solver):
                print 'Warning: no record for', solver, 'on', lib
            else:
                if outcome[lib][solver] != reference[lib][solver]:
                    print 'Warning: discrepancy on ', lib, 'with', solver+':'
                    print 'was:'
                    print reference[lib][solver]
                    print 'now:'
                    print outcome[lib][solver]
                    if param['update'] == 'yes': 
                        print '... saving'
                        reference[lib][solver] = outcome[lib][solver]
                else:
                    print 'ok'

for lib in outcome.keys():
    for solver in outcome[lib].keys():
        if not reference.has_key(lib):
            print lib, 'is a new problem... saving'
            reference[lib] = outcome[lib]
        else:
            if not reference[lib].has_key(solver):
                print solver, 'did not previously run on', lib+'... saving'
                reference[lib][solver] = outcome[lib][solver]


cPickle.dump(reference, open('.unit_test_reference.pkl', 'w'))






