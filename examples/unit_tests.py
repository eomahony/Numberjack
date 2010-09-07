from Numberjack import *

import subprocess
ls = subprocess.Popen(['ls', '-1'], stdout=subprocess.PIPE)

for example in ls.stdout:
    if example[-4:-1] == '.py' and example[:-4] != 'unit_tests':
        lib = __import__(example[:-4])
        param = lib.default
        for solver in lib.solvers:
            print 'run', example[:-4], 'with', solver
            param['solver'] = solver
            lib.solve(param)
            print ''





