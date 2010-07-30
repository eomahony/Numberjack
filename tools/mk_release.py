#! /usr/bin/env python

import subprocess
import os


version_number = '0'
release_number = '1'
svn_revision = ''


subprocess.Popen(['git', 'pull'], stdout=subprocess.PIPE).wait()
subprocess.Popen(['git', 'commit', '-a', '-m', "archive"]).wait()
#subprocess.Popen(['svn', 'update'], stdout=subprocess.PIPE).wait()
#subprocess.Popen(['svn', 'commit', '-m', "archive"]).wait()
#update = subprocess.Popen(['svn', 'update'], stdout=subprocess.PIPE)
#update.wait()

"""
for line in update.stdout:
    print line,
    splitted = line.split()
    if len(splitted) == 3:
        if splitted[0] == 'At' and splitted[1] == 'revision':
            svn_revision = splitted[2][:-1]
"""

subprocess.Popen(['date', '"+%y-%m-%d"'], stdout=subprocess.PIPE).wait()

suffix = version_number+'.'+release_number+'.'+svn_revision

print suffix

exit(1)
            
subprocess.Popen( ['svn', 'export', 
                   'http://4c110.ucc.ie:8000/svn/Numberjack', 
                   'releases/Numberjack-'+suffix] ).wait()

#svn export http://4c110.ucc.ie:8000/svn/Roadef2010

os.chdir( './releases/' )

the_solvers = ['mistral', 'gecode', 'casper', 'minisat', 'scip']


# remove useless things from mistral
subprocess.Popen( ['rm', '-r', 
                   'Numberjack-'+suffix+'/solvers/mistral/data']).wait()
subprocess.Popen( ['rm', '-r', 
                   'Numberjack-'+suffix+'/solvers/mistral/tools']).wait()

# make Numberjack archive
subprocess.Popen( ['tar', '-cvzf', 
                   'Numberjack_complete-'+suffix+'.tgz', 
                   'Numberjack-'+suffix] ).wait()

for sol in the_solvers:
    # make Mistral archive
    subprocess.Popen( ['mv', 'Numberjack-'+suffix+'/solvers/'+sol+'/', './'] )
    subprocess.Popen( ['tar', '-cvzf', 
                       'nbj_'+sol+'-1.'+str(int(svn_revision))+'.tgz', 
                       sol] ).wait()    
    # remove the exported files
    subprocess.Popen( ['rm', '-r', sol] )


# make archive without back-end solvers
subprocess.Popen( ['tar', '-cvzf', 
                   'Numberjack-'+suffix+'.tgz', 
                   'Numberjack-'+suffix] ).wait()

# remove the exported files
subprocess.Popen( ['rm', '-r', 'Numberjack-'+suffix] )

