from Numberjack import *

import subprocess
ls = subprocess.Popen(['ls', '-1'], stdout=subprocess.PIPE)

for example in ls.stdout:
    if example[-4:-1] == '.py' and example[:-4] != 'unit_tests':

        lib_str = example[:-4]

        print '**************************************'
        print '******', lib_str
        print '**************************************'
        lib = __import__(lib_str)
        param = lib.default
        for solver in lib.solvers:
            print 'run', lib_str, 'with', solver
            param['solver'] = solver
            lib.solve(param)
            print ''





