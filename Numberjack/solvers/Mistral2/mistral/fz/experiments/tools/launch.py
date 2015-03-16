#! /usr/bin/env python

import os
import signal
import threading
import subprocess
import sys
import timing
from threading import Thread

N = 1;

input = 'dir'
source = '/storagespace/store1/XML/Normalised/csp/'
timelimit = 60
key = None
solverdir = './'
bSeed = 11041979
eSeed = 11041980
step = 1

lock = threading.Lock()

def getparameters( paramlist ):
    if len(paramlist) >= 1:
        global N
        global key
        global input
        global source
        global timelimit
        global solverdir
        global bSeed
        global eSeed
        global step
        offset=2
        if paramlist[0][:2] == '-j':
            N = int(paramlist[0][2:])
            if N > 8:
                N = 8
            offset=1
        elif paramlist[0] == '-key':
            if len(paramlist) >= 2:
                key = paramlist[1]
        elif paramlist[0] == '-solverdir':
            if len(paramlist) >= 2:
                solverdir = paramlist[1]
        elif paramlist[0] == '-seed':
            stuple = []
            for value in paramlist[1:]:
                if value[0] == '-': 
                    break
                else: 
                    stuple.append( value )
            if len(stuple) > 0:
                bSeed = int(stuple[0])
            if len(stuple) > 1:
                eSeed = int(stuple[1])
            if len(stuple) > 2:
                step = int(stuple[2])
            offset = 1+len(stuple)
        elif paramlist[0] == '-dir':
            input = 'dir'
            if len(paramlist) >= 2:
                source = paramlist[1]
        elif paramlist[0] == '-file':
            input = 'list'
            if len(paramlist) >= 2:
                source = paramlist[1]
        elif paramlist[0] == '-timelimit':
            if len(paramlist) >= 2:
                timelimit = int(paramlist[1])
        else: print 'unknown keyword: \"'+paramlist[0]+'\". The next param (\"'+paramlist[1]+'\") is discarded '
        getparameters( paramlist[offset:] )

getparameters( sys.argv[1:] )



dateproc = subprocess.Popen( ['date', '+%y%m%d%H%M%S'], stdout=subprocess.PIPE )
date = dateproc.stdout.readline()[:-1]
dateproc.stdout.close()
hostproc = subprocess.Popen( ['hostname'], stdout=subprocess.PIPE )
host = hostproc.stdout.readline()[:-1]
hostproc.stdout.close()
key = 'exp_'+host+'_'+date+'_'+key

sys.stderr.write( 'Used key: '+key+'\n' )

key_file = open( 'results/'+key+'.key', 'w' )

def kill( process ):
    if process.poll() == None:
        os.kill( process.pid, signal.SIGKILL )

# get solver names
solvers = []
process = subprocess.Popen( ['find', solverdir, '-name', '*.slv'], stdout=subprocess.PIPE )
for line in process.stdout:
    solvers.append( line[line.rfind('/')+1:-5] )
process.stdout.close()

key_file.write( str(len(solvers)) )
for solver in solvers:
    key_file.write( ' '+solver )
key_file.write( ' '+str(timelimit) )
key_file.write( ' '+str(N) )
key_file.flush()


# get command lines for each solver
commandlines = {}.fromkeys( solvers )
for solver in solvers:
    commandlines[solver] = []
    command = open(solverdir+"/"+solver+".slv").readline().split()
    for i in range( len(command) ):
        j = command[i].find('HOME')
        if j >= 0:
            command[i] = (command[i][:j]+'src/'+command[i][(j+4):])
    commandlines[solver].extend( command )

# get bench's paths
benchpaths = []
if input == 'list': # from list
    for line in open( source ):
        benchpaths.append( line[:-1] )
else: # form directory
    benchproc = subprocess.Popen( ['find', source, '-name', '*.*'], stdout=subprocess.PIPE )
    for line in benchproc.stdout:
        benchpaths.append( line[:-1] )
    benchproc.stdout.close()


def limit( process, cutoff ):
    t = threading.Timer( cutoff, kill, [process] )
    t.start()
    return t

#print 'THE BENCHS:'
#print benchpaths

#print 'THE SOLVERS:'
#print solvers

jobs = []
for bench in benchpaths:
    #if len(bench) > 1:
    for seed in range(bSeed, eSeed, step):
        for solver in solvers:
            jobs.append( [solver, bench, str(seed)] )
GLOBAL_COUNTER = 0;

print 'THE JOBS:'
print jobs


class ExperimentRunner(Thread):
    solvers = None
    id = None
    solution_file = None
    result_file = None
    error_file = None
    path_file = None
    stat_file = None
    log_file = None
    chrono = None
    def __init__ (self,id,slvs):
        Thread.__init__(self)
        self.id = id
        self.solution_file = {}.fromkeys( solvers )
        self.result_file   = {}.fromkeys( solvers )
        self.error_file    = {}.fromkeys( solvers )
        self.path_file     = {}.fromkeys( solvers )
        self.stat_file     = {}.fromkeys( solvers )
        self.log_file      = {}.fromkeys( solvers )
        
        self.chrono = timing.PTimer()
        self.solvers = slvs

        for solver in self.solvers:
            sys.stderr.write('       Core'+str(self.id)+' open files for '+solver+'\n') 
            self.solution_file[solver] = open( 'results/'+key+'_'+solver+'.sol'+str(self.id), 'w' ) 
            self.result_file[solver] = open( 'results/'+key+'_'+solver+'.res'+str(self.id), 'w' ) 
            self.error_file[solver] = open( 'results/'+key+'_'+solver+'.err'+str(self.id), 'w' ) 
            self.path_file[solver] = open( 'results/'+key+'_'+solver+'.pth'+str(self.id), 'w' ) 
            self.stat_file[solver] = open( 'results/'+key+'_'+solver+'.sta'+str(self.id), 'w' ) 
            self.log_file[solver] = open( 'results/'+key+'_'+solver+'.log'+str(self.id), 'w' ) 


    def run(self):
        global lock
        while 1:

    ##############################################################################
    ########### HERE WE NEED TO LOCK ALL THESE FILES IF WE MULTITHREAD ###########
    ##############################################################################
            '''
            lock.acquire()
            try:
            '''
            global GLOBAL_COUNTER
            self.idx = GLOBAL_COUNTER
            GLOBAL_COUNTER += 1
            
            '''
            finally:
                lock.release()
            '''
    ###########################################################
    ########### HERE WE CAN UNLOCK THE RESULT FILES ###########
    ###########################################################
                
            if self.idx >= len(jobs): break

            solver = jobs[self.idx][0]
            bench = jobs[self.idx][1]
            seed = jobs[self.idx][2]

            print str(int(float(self.idx)/float(len(jobs))*100.0)).rjust(3)+'% - Core'+str(self.id)+' run '+solver+' on '+bench[bench.rfind('/')+1:]+' with seed ' + seed

            solution = 'None';
            result = '---------';
            stat = '\n';

            #print commandlines[solver]

            command = [commandlines[solver][0]];
            command.append( '--seed' )
            command.append( str(seed) )
            for word in commandlines[solver][1:]:
                if word != 'BENCHNAME':
                    command.append( word )
                else:
                    command.append( bench )

            #print command

#################################################
            print '     ',
            for word in command:
                print word,
            print '\n'
            self.chrono.start()
            process = subprocess.Popen( command, bufsize=32768, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
            clock = limit( process, timelimit )
            process.wait()
            self.chrono.finish()
            output, error = process.communicate() #process.stdout
#################################################

            clock.cancel()
            seedbench = bench[(bench.rfind('/')+1):]
            dotidx = bench.rfind('.')
            if seed != '0':
                if dotidx >= 0:
                    seedbench = seedbench[:-4]+('^'+str(seed));
                else:
                    seedbench = seedbench+('^'+str(seed));


            #self.solution_file[solver].flush()
            #self.result_file[solver].flush()
            #self.error_file[solver].flush()
            #self.stat_file[solver].flush()
            #self.path_file[solver].flush()
            #self.solution_file[solver].write( seedbench+' ' )
            #self.result_file[solver].write( seedbench+' ' )
            #self.error_file[solver].write( seedbench+' ' )
            #self.stat_file[solver].write( seedbench+' ' )
            #self.path_file[solver].write( seedbench+' '+bench+'\n' )

            solution_out = seedbench+' '
            result_out = seedbench+' '
            error_out = seedbench+' '
            stat_out = seedbench+' '
            path_out = seedbench+' '+bench+'\n'


            line = ''
            for char in output:
                if char == '\n':
                    self.log_file[solver].write( line+'\n' )
                    data_line = line.split()
                    if len(data_line) > 0:
                        if data_line[0] == 'v':
                            solution = line[(line.find('v')+1):]
                        elif data_line[0] == 's':
                            result = data_line[1];
                        elif data_line[0] == 'd':
                            stat = line[(line.find('d')+1):]+' '+stat
                    line = ''
                else:
                    line += char

            #self.solution_file[solver].write( solution+'\n' )
            solution_out += (solution+'\n')
            #self.result_file[solver].write( result+'\n' )    
            result_out += (result+'\n')
            stat_out += ('TIME '+str(float(int(self.chrono.seconds() * 1000))/1000)+stat)
            #self.stat_file[solver].write( stat )

            line = ''
            for char in error:
                if char == '\n':
                    #self.error_file[solver].write( line, )
                    error_out += line
                    line = ''
                else:
                    line += char
            #self.error_file[solver].write( '\n' )
            error_out += '\n'

            process.stdout.close()
            process.stderr.close()

            #lock.acquire()
            self.solution_file[solver].flush()
            self.result_file[solver].flush()
            self.error_file[solver].flush()
            self.stat_file[solver].flush()
            self.path_file[solver].flush()
            self.log_file[solver].flush()
            self.solution_file[solver].write( solution_out )
            self.result_file[solver].write( result_out )
            self.error_file[solver].write( error_out )
            self.stat_file[solver].write( stat_out )
            self.path_file[solver].write( path_out )
            #lock.release()

        for solver in self.solvers:
            sys.stderr.write('       Core'+str(self.id)+' close files for '+solver+'\n') 
            self.solution_file[solver].close();
            self.result_file[solver].close();
            self.error_file[solver].close();
            self.path_file[solver].close();
            self.stat_file[solver].close();
            self.log_file[solver].close();


ER = []
for i in range(N):
    ER.append( ExperimentRunner(i, solvers) )
    ER[i].start()






