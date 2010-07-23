#! /usr/bin/env python

import subprocess

modules = []
list = subprocess.Popen(['ls', '-1'], stdout=subprocess.PIPE)
for line in list.stdout:
    if line[-4:] == '.py\n':
        modules.append(line[:-4])

from distutils.core import setup

setup(name='Numberjack'+modules[0],
      version='1.0',
      py_modules=modules,
      )


import sys
all_paths = sys.path
for path in all_paths:    
    try:
        #print 'ls', path
        ls = subprocess.Popen(['ls', '-1', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ls.wait()
        for pylib in ls.stdout:
            pylib = pylib[:-1]
            #print '\t', pylib,
            for module in modules:
                #print '=', module, '?',
                if pylib[:-3] == module:                    
                    #print 'Moving', '_'+module+'.so to', path+'/'
                    cp = subprocess.Popen(['cp', '_'+module+'.so', path+'/'], stderr=subprocess.PIPE)
                    cp.wait()
                    break
                #print ''
    except:
        continue


