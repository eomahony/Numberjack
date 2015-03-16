#! /usr/bin/env python

import time
import sys
import subprocess
import threading
from threading import Thread
from sets import Set
import math


GLOBAL_COUNTER = 0
GLOBAL_END = 0
LAST_THREAD = 1
KEY2 = None

lock = threading.Lock()


def pretty( val, prec, width, ul=False, bf=False ):

    if prec == 0:
        string_rep = str(int(val))
    else:
        fvalue = float(val)
    #    ivalue = round(fvalue)
    #    if ivalue == fvalue:
    #        string_rep = str(int(ivalue))
    #    else:
        string_rep = str(fvalue)
        dotidx = string_rep.find('.')
        if dotidx >= 0:
            #if fvalue > 50.0:
            #    string_rep = string_rep[:dotidx]
            #else:
            #    string_rep = string_rep[:dotidx+3]
            while len(string_rep) <= (dotidx+prec):
               string_rep += '0'; 
            string_rep = string_rep[:dotidx+prec+1]        
        else:
            string_rep += '.'
            for i in range(prec):
                string_rep += '0'

    if bf:
        string_rep = '\\textbf{'+string_rep+'}'
    if ul:
        string_rep = '\underline{'+string_rep+'}'
    return string_rep.rjust(width)



def checkSolution( solution, bench ):
    process = subprocess.Popen( ['tools/CSP-XML-parser/tools/CSP-verifier', bench], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    answer = process.communicate( solution )[0]
    if answer == 'OK\n':
        return True
    return False


class solutionChecker(Thread):    
    def __init__ (self, jobs, verified, problem, id, parser):
        Thread.__init__(self)
        self.id = id
        self.jobs = jobs
        self.verified = verified
        self.problem = problem
        self.parser = parser

    def run(self):
        global GLOBAL_END
        global LAST_THREAD
        global lock
        try:
            while 1:
                ##############################################################################
                ########### HERE WE NEED TO LOCK ALL THESE FILES IF WE MULTITHREAD ###########
                ##############################################################################
                lock.acquire()
                try:
                    global GLOBAL_COUNTER
                    idx = GLOBAL_COUNTER
                    GLOBAL_COUNTER += 1
                finally:
                    lock.release()
                ###########################################################
                ########### HERE WE CAN UNLOCK THE RESULT FILES ###########
                ###########################################################

                    if idx >= GLOBAL_END:                    
                        break

                    solver = self.jobs[idx][0]
                    solution = self.jobs[idx][1]
                    path = self.jobs[idx][2]
                    bench = path[path.rfind('/')+1:]

                    print str(int(float(idx)/float(GLOBAL_END)*100.0)).rjust(3)+'% Core'+str(self.id), 'check', solver+'\'s solution for', bench.ljust(50)
                    sys.stdout.flush()
                    correct = False
                    process = subprocess.Popen( ['tools/CSP-XML-parser/tools/CSP-verifier', path], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                    answer = process.communicate( solution )[0]
                    if answer == 'OK\n':
                        correct = True

                    #lock.acquire()
                    #try:
                    if correct == True:
                        self.verified[solver].add( bench )
                    else:
                        print ' /!\\ WARNING,', solver, 'gives a wrong certificate for', bench, ' /!\\'
                        sys.stdout.flush()
                        self.problem[solver].append( bench )        
                    #finally:
                    #    lock.release()
        finally:
            lock.acquire()
            try:
                LAST_THREAD -= 1
                print 'core'+str(self.id), 'stop', str(LAST_THREAD)
                if LAST_THREAD == 0:
                    print 'core'+str(self.id), ' is the last thread'
                    print 'core'+str(self.id), ' call the parser'
                    self.parser.printWrongCertificate();
            finally:
                lock.release()





class solutionChecker2(Thread):    
    def __init__ (self, jobs, id):
        Thread.__init__(self)
        self.id = id
        self.jobs = jobs
        self.verified = {}
        self.problem = {}

    def run(self):
        global GLOBAL_END
        global lock
        global KEY2
        while 1:
            ##############################################################################
            ########### HERE WE NEED TO LOCK ALL THESE FILES IF WE MULTITHREAD ###########
            ##############################################################################
            lock.acquire()

            global GLOBAL_COUNTER
            idx = GLOBAL_COUNTER
            GLOBAL_COUNTER += 1

            lock.release()
            ###########################################################
            ########### HERE WE CAN UNLOCK THE RESULT FILES ###########
            ###########################################################

            if idx >= GLOBAL_END:                    
                break

            solver = self.jobs[idx][0]
            solution = self.jobs[idx][1]
            path = self.jobs[idx][2]
            bench = path[path.rfind('/')+1:]


            print str(int(float(idx)/float(GLOBAL_END)*100.0)).rjust(3)+'% Core'+str(self.id), 'check', solver+'\'s solution for', bench.ljust(50)
            sys.stdout.flush()
            correct = False
            process = subprocess.Popen( ['tools/CSP-XML-parser/tools/CSP-verifier', path], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            answer = process.communicate( solution )[0]
            print '     Core'+str(self.id), 'has finished with', solver+'\'s solution for', bench.ljust(50)
            sys.stdout.flush()
            if answer == 'OK\n':
                if self.verified.has_key( solver ) == False:
                    self.verified[solver] = Set([])
                self.verified[solver].add( bench )
            else:
                print ' /!\\ WARNING,', solver, 'gives a wrong certificate for', bench, ' /!\\'
                sys.stdout.flush()
                if self.problem.has_key( solver ) == False:
                    self.problem[solver] = []
                self.problem[solver].append( bench )        
        print self.problem
        print self.verified




class Parser:
    '''
    # parameters 
    printlatex = False
    printsolved = False
    printeasy = 0
    printhard = 0
    inputSolvers = []
    list = False
    sample = []
    summary = True
    verif = False
    inputStats = []
    inputStatsSet = Set([])
    keys = []
    timelimit = 0
    showerror = Set([])
    reference = None
    delta = Set([])
    dtype = Set(['res', 'sta'])
    nThreads = 1

    # raw data merged result files
    raw_data = {}.fromkeys( ['err', 'log', 'pth', 'res', 'sol', 'sta'] )

    # utils
    setOfSolvers = Set([])
    solvers = []
    result = None
    stat = None
    benchmark = None # the set of instances for each solver
    solved_by_one = None # the set of instances solved by at least one solver
    solved_bench = None # the set of instances solved by each solver
    error_bench = None # the set of instances prompting an error, for each solver
    stat_keyword = None # the set of stats in output for each solver
    statistics = None # the set of stats that need to be displayed for each solver
    common_bench = None # the set of instances on which every solver ran
    ordered_bench = None # the list of common instances, in lexico order
    common_solved = None # the set of instances solved by every solver
    solver_length = 10;
    bench_length = 10;

    # verification
    jobs = None
    wrong_certificate = None
    wrong_proof = None
    verified = None
    '''

    def __init__ (self,paramlist):
        self.onlyaverage = False;
        self.latex_args = []
        self.timelimit = 0
        self.printlatex = False
        self.printsolved = False
        self.printeasy = 0
        self.printhard = 0
        self.list = False
        self.summary = True
        self.verif = False
        self.reference = None
        self.nThreads = 1
        self.inputGBench = []
        self.inputBBench = []
        self.inputSolvers = []
        self.inputStats = []
        self.inputStatsSet = Set([])
        self.keys = []
        self.showerror = Set([])
        self.delta = Set([])
        self.dtype = Set(['res', 'sta'])
        self.raw_data = {}.fromkeys( ['err', 'log', 'pth', 'res', 'sol', 'sta'] )
        self.setOfSolvers = Set([])
        self.solvers = []
        self.sample = {}
        self.sample_stats = []
        self.emph = {}
        self.emph_stats = []
        self.underlineoptimal = False
        self.result = None
        self.stat = None
        self.benchmark = None # the set of instances for each solver
        self.solved_by_one = None # the set of instances solved by at least one solver
        self.solved_bench = None # the set of instances solved by each solver
        self.error_bench = None # the set of instances prompting an error, for each solver
        self.stat_keyword = None # the set of stats in output for each solver
        self.statistics = None # the set of stats that need to be displayed for each solver
        self.common_bench = None # the set of instances on which every solver ran
        self.ordered_bench = None # the list of common instances, in lexico order
        self.common_solved = None # the set of instances solved by every solver

        self.precision = {}
        self.replace = {}

        self.solver_length = 10
        self.bench_length = 10

        # verification
        self.jobs = None
        self.wrong_certificate = None
        self.wrong_proof = None
        self.verified = None

        self.getparameters( paramlist )


    def getparameters( self, paramlist ):
        if len(paramlist) >= 1:
            global LAST_THREAD
            global GLOBAL_END
            global KEY2
            increment = 2
            if paramlist[0] == '-key':
                for key in paramlist[1:]:
                    if key[0] == '-': 
                        break
                    else: 
                        self.keys.append( key )
                KEY2 = self.keys[0]
                increment = 1+len(self.keys)
            elif paramlist[0] == '-solvers':
                for solver in paramlist[1:]:
                    if solver[0] == '-': 
                        break
                    else: 
                        self.inputSolvers.append( solver )
                increment = 1+len(self.inputSolvers)
            elif paramlist[0] == '-benchmarks':
                for bench in paramlist[1:]:
                    if bench[0] == '-': 
                        break
                    else: 
                        if bench[0] != '^':
                            self.inputGBench.append( bench )
                        else:
                            self.inputBBench.append( bench[1:] )
                increment = 1+len(self.inputGBench)+len(self.inputBBench)
            elif paramlist[0][:2] == '-j':
                self.nThreads = int(paramlist[0][2:])
                if self.nThreads > 8:
                    self.nThreads = 8
                LAST_THREAD = self.nThreads
                increment=1
            elif paramlist[0] == '-list':
                self.list = True
                increment = 1
            elif paramlist[0] == '-sample':
                stat = 'TIME'
                increment = 1
                for word in paramlist[1:]:
                    if word[0] == '-': 
                        break
                    elif word == 'avg' or word == 'min' or word == 'max' or word == 'dev' or word == 'med' or word == 'sig': 
                        self.sample[stat].append(word)
                        increment += 1
                    else:
                        increment += 1
                        stat = word
                        if not self.sample.has_key(stat):
                            self.sample[stat] = []
                            self.sample_stats.append( stat )
                            self.inputStats.append( stat )
                            self.inputStatsSet.add( stat )
            elif paramlist[0] == '-emph':
                stat = 'TIME'
                increment = 1
                for word in paramlist[1:]:
                    if word[0] == '-': 
                        break
                    elif word == 'max' or word == 'min':
                        self.emph[stat] = word
                        increment += 1
                    else:
                        increment += 1
                        stat = word
                        if not self.emph.has_key(stat):
                            self.emph[stat] = None
                            self.emph_stats.append( stat )
                            self.inputStats.append( stat )
                            self.inputStatsSet.add( stat )
            elif paramlist[0] == '-latex':                
                self.printlatex = True
                increment = 1
                cur_args = []
                while len(paramlist[increment:]):
                    word = paramlist[increment]
                    if word == '+':
                        self.latex_args.append(cur_args)
                        cur_args = []
                    else:
                        cur_args.append(word)
                    increment += 1
                self.latex_args.append(cur_args)
                return

            elif paramlist[0] == '-solved':
                self.printsolved = True
                increment = 1
            elif paramlist[0][:5] == '-easy':
                self.printeasy = float(paramlist[0][5:])
                increment = 1
            elif paramlist[0][:5] == '-uopt':
                self.underlineoptimal = True
                increment = 1
            elif paramlist[0][:8] == '-average':
                self.onlyaverage = True
                increment = 1
            elif paramlist[0][:5] == '-hard':
                self.printhard = float(paramlist[0][5:])
                increment = 1
            elif paramlist[0][:6] == '-verif':
                if len(paramlist[0]) > 6:
                    GLOBAL_END = int(paramlist[0][6:])
                increment = 1
                self.verif = True
            elif paramlist[0] == '-delta':
                increment = 1
                for solver in paramlist[1:]:
                    if solver[0] == '-': 
                        break
                    elif self.reference == None:
                        self.reference = solver
                        increment = 2;
                    else:
                        self.delta.add( solver )
                increment += len(self.delta)   
            elif paramlist[0] == '-showerror':
                for solver in paramlist[1:]:
                    if solver[0] == '-': 
                        break
                    else:
                        self.showerror.add( solver )
                increment += 1+len(self.showerror)  
            elif paramlist[0] == '-nosummary':
                self.summary = False
                increment = 1
            elif paramlist[0] == '-stats':
                increment = 1
                for st in paramlist[1:]:
                    if st[0] == '-': 
                        break
                    else: 
                        self.inputStats.append( st )
                        self.inputStatsSet.add( st )
                    increment = 1+increment
            elif paramlist[0] == '-precision':
                increment = 1
                for i in range(0, len(paramlist[1:]), 2):
                    if paramlist[1+i][0] == '-' or i>=len(paramlist[1:]): 
                        break
                    else: 
                        self.precision[paramlist[1+i]] = int(paramlist[2+i])
                    increment = 2+increment
            elif paramlist[0] == '-replace':
                increment = 1
                for i in range(0, len(paramlist[1:]), 3):
                    if paramlist[1+i][0] == '-' or i>=len(paramlist[1:]): 
                        break
                    else: 
                        #print paramlist[1+i], '<-', paramlist[2+i], ',' , paramlist[3+i]
                        self.replace[paramlist[1+i]] = (int(paramlist[2+i]), paramlist[3+i])
                    increment = 3+increment
            else: print 'unknown keyword: \"'+paramlist[0]+'\". The next param (\"'+paramlist[1]+'\") is discarded '
            self.getparameters( paramlist[increment:] )


    def significance( self, benches, solver1, solver2, statistic ):
        if solver1 == solver2:
            return 0
        else:
            rscript = open('/storagespace/store1/CPAIExpe/tmpRscript.R', 'w')
            rscript.write('t = t.test( x = c(' )
            seq = '';
            for bench in benches:
                seq += (self.stat[solver1][bench][statistic]+',')
            rscript.write( seq[:-1]+'), y = c(' )
            seq = '';
            for bench in benches:
                seq += (self.stat[solver2][bench][statistic]+',')
            rscript.write( seq[:-1]+'), alternative = c("two.sided"), paired = TRUE )\n' )
            rscript.write( '1-t$p.value\n' )
            rscript.close()
            R = subprocess.Popen( ['/home/ehebrard/util/R-2.9.0/bin/Rscript', 'tmpRscript.R'], stdout=subprocess.PIPE, stderr=subprocess.PIPE )
            x = float( R.stdout.readline().split()[1] )
        #x = 1
            return x


    def average( self, benches, solver, statistic ):
        avg = 0.0
        for bench in benches:
            avg += float(self.stat[solver][bench][statistic])
        avg /= len(benches)
        return avg

    def minimum( self, benches, solver, statistic ):
        tmin = sys.maxint
        for bench in benches:
            x = float(self.stat[solver][bench][statistic])
            if x < tmin:
                tmin = x
        return tmin

    def maximum( self, benches, solver, statistic ):
        tmax = -sys.maxint
        for bench in benches:
            x = float(self.stat[solver][bench][statistic])
            if x > tmax:
                tmax = x
        return tmax

    def median( self, benches, solver, statistic ):
        seq = []
        for bench in benches:
            seq.append( float(self.stat[solver][bench][statistic]) )
        seq.sort()
        return seq[len(seq)/2]

    def deviation( self, benches, solver, statistic ):
        avg = 0.0
        for bench in benches:
            avg += float(self.stat[solver][bench][statistic]) 
        avg /= len(benches)
        dev = 0.0
        for bench in benches:
            x = abs( float(self.stat[solver][bench][statistic]) - avg )
            dev += (x*x)
        return math.sqrt(dev)



    def _average( self, stats ):
        avg = 0.0
        for stat in stats:
            avg += stat
        avg /= len(stats)
        return avg

    def _minimum( self, stats ):
        tmin = sys.maxint
        for stat in stats:
            x = stat
            if x < tmin:
                tmin = x
        return tmin

    def _maximum( self, stats ):
        tmax = -sys.maxint
        for stat in stats:
            x = stat
            if x > tmax:
                tmax = x
        return tmax

    def _median( self, stats ):
        seq = []
        for stat in stats:
            seq.append( stat )
        seq.sort()
        return seq[len(seq)/2]

    def _deviation( self, stats ):
        avg = 0.0
        for stat in stats:
            avg += stat 
        avg /= len(stats)
        dev = 0.0
        for stat in stats:
            x = abs( stat - avg )
            dev += (x*x)
        return math.sqrt(dev)

    def initialize(self):                
        N = {}.fromkeys( self.keys )
        solverByHost = {}.fromkeys( self.keys )
        hostnames = {}.fromkeys( self.keys )

        ############################################################
        # We try to make something out of the mess of hostnames/keys/solver-names/threads
        ############################################################
        for key in self.keys: # We get all result-file matching the keys
            # then we find out what nodes were used and store this info
            hostnames[key] = Set([])
            for line in subprocess.Popen( ['find', 'results', '-name', 'exp_*_'+key+'.key'], stdout=subprocess.PIPE ).stdout:
                hostnames[key].add( line[:-5] )

            # N will store the number of threads used for each run
            N[key] = {}.fromkeys( hostnames[key] )
            solverByHost[key] = {}.fromkeys( hostnames[key] ) 

            for host in hostnames[key]:
                # the rest of the data we need is in the key-files
                data = open( host+'.key' ).readline().split()
                nSolvers = int(data[0])
                N[key][host] = int(data[nSolvers+2])
                if self.timelimit == 0:
                    self.timelimit = int(data[nSolvers+1])
                else:
                    if self.timelimit != int(data[nSolvers+1]):
                        print 'WARNING!!', 'inconsistent timelimit on different nodes'
                slv = Set(data[1:nSolvers+1])
                solverByHost[key][host] = slv
                self.setOfSolvers |= slv

        ############################################################
        # We init the list of solvers, intersection between the user
        # defined ones and all solvers encountered in the key-files
        ############################################################
        if len(self.inputSolvers) == 0:
            self.solvers.extend( self.setOfSolvers )
        else:
            for solver in self.inputSolvers:
                if solver in self.setOfSolvers:
                    self.solvers.append( solver )
        self.setOfSolvers &= Set(self.solvers)
        for key in self.keys:
            for host in hostnames[key]:
                solverByHost[key][host] &= self.setOfSolvers

        self.benchmark = {}.fromkeys( self.solvers )
        for solver in self.solvers:
            if self.solver_length < len(solver):
                self.solver_length = len(solver)
            self.benchmark[solver] = Set([])

        if self.verif == True:
            self.dtype.add('sol')
        if (self.verif == True) | (self.printsolved == True) | ((self.printeasy + self.printhard) > 0) | (len(self.showerror) > 0) | (len(self.delta) > 0):
            self.dtype.add('pth')

        ############################################################
        # We collect the data into the struct 'raw_data'
        ############################################################        
        bench = None
        for dat in self.dtype:
            repeat = {}.fromkeys( self.solvers )
            for solver in self.solvers:
                repeat[solver] = {}
            self.raw_data[dat] = {}.fromkeys( self.solvers )
            for solver in self.solvers:
                self.raw_data[dat][solver] = {}
            for key in self.keys:
                for host in hostnames[key]:
                    for solver in solverByHost[key][host]:
                        for id in range(N[key][host]):
                            for line in open( host+'_'+solver+'.'+dat+str(id) ):
                                data_line = line.split()
                                bench = data_line[0]
                                if repeat[solver].has_key( bench ):
                                    repeat[solver][bench] += 1
                                    bench = (bench[:-4]+'_'+str(repeat[solver][bench])+bench[-4:])
                                else:
                                    repeat[solver][bench] = 0
                                self.benchmark[solver].add( bench )
                                if dat == 'sol':
                                    self.raw_data[dat][solver][bench] = line[len(data_line[0]):]
                                else:
                                    self.raw_data[dat][solver][bench] = data_line[1:]

        self.result = {}.fromkeys( self.solvers )
        self.stat = {}.fromkeys( self.solvers )
        self.solved_bench = {}.fromkeys( self.solvers )
        self.error_bench = {}.fromkeys( self.solvers )
        self.stat_keyword = {}.fromkeys( self.solvers )
        self.wrong_certificate = {}.fromkeys( self.solvers )
        self.wrong_proof = {}.fromkeys( self.solvers )

        for solver in self.solvers:
            self.result[solver] = {}.fromkeys( self.benchmark[solver] )
            self.stat[solver] = {}.fromkeys( self.benchmark[solver] )
            self.solved_bench[solver] = Set([])
            self.error_bench[solver] = Set([])
            self.stat_keyword[solver] = Set(['PROOF'])
            self.wrong_certificate[solver] = []
            self.wrong_proof[solver] = []
            outcome = ''
            for bench in self.benchmark[solver]:
                self.stat[solver][bench] = {}
                sta_line = self.raw_data['sta'][solver][bench]
                for i in range( 0, len(sta_line)-1, 2 ):
                    self.stat[solver][bench][sta_line[i]] = sta_line[i+1]
                    self.stat_keyword[solver].add( sta_line[i] )
                    outcome = self.raw_data['res'][solver][bench][0]
                self.result[solver][bench] = outcome
                self.stat[solver][bench]['PROOF'] = 0
                if( outcome[-11:] == 'SATISFIABLE' or 
                    outcome[-11:] == 'UNSATISFIABLE' or 
                    outcome[-11:] == 'OPTIMAL' 
                    ):
                    self.stat[solver][bench]['PROOF'] = 1
                    self.solved_bench[solver].add( bench )
                elif outcome[-11:] == 'SUBOPTIMAL':
                    self.solved_bench[solver].add( bench )
                elif outcome != 'UNKNOWN':
                    if len(sta_line) >= 2:
                        rtime = float(sta_line[1])
                        if rtime < (self.timelimit/2):
                            self.result[solver][bench] = 'error--------'
                            self.error_bench[solver].add( bench )


        """
        print 'raw_data:\n', self.raw_data
        
        print 'result:\n', self.result

        print 'solved:\n', self.solved_bench

        print 'error:\n', self.error_bench

        sys.exit(1)
        """

        ############################################################
        # Compute the set of benches on which the comparison will be made
        ############################################################
        
        if len(self.inputGBench)+len(self.inputBBench) == 0:
            self.common_bench = self.benchmark[self.solvers[0]].copy();
            #print self.common_bench
            for solver in self.solvers[1:]:
                #print self.benchmark[solver]
                self.common_bench &= self.benchmark[solver]
        else:
            self.common_bench = Set([])
            tmp_common_bench = self.benchmark[self.solvers[0]].copy();
            for solver in self.solvers[1:]:
                tmp_common_bench &= self.benchmark[solver]

            for bench in tmp_common_bench:
                fine = (len(self.inputGBench) == 0)
                if not fine:
                    for tag in self.inputGBench:
                        if bench.find(tag) >= 0:
                            fine = True
                            break
                if fine and len(self.inputBBench) > 0:
                    for tag in self.inputBBench:
                        if bench.find(tag) >= 0:
                            fine = False
                            break
                if fine:
                    self.common_bench.add(bench)


        self.ordered_bench = []
        for bench in self.common_bench:
            benchname = bench
            begin = benchname.rfind( 'normalized-' );
            if begin >= 0:
                benchname = benchname[begin+11:]
            dot = benchname.rfind( '.' )
            if dot >= 0:
                benchname = benchname[:dot]
            circ = benchname.rfind( '^' )
            if circ >= 0:
                benchname = benchname[:circ]
            if self.bench_length < len(benchname):
                self.bench_length = len(benchname)
            self.ordered_bench.append( bench )
        self.ordered_bench.sort();
        self.solved_bench[self.solvers[0]] &= self.common_bench
        self.error_bench[self.solvers[0]] &= self.common_bench
        self.common_solved = self.solved_bench[self.solvers[0]].copy();
        self.solved_by_one = self.solved_bench[self.solvers[0]].copy();
        for solver in self.solvers[1:]:
            self.solved_by_one |= self.solved_bench[solver]
            self.solved_bench[solver] &= self.common_bench
            self.error_bench[solver] &= self.common_bench
            self.common_solved &= self.solved_bench[solver]    


        ############################################################
        # Compute the set of statistics that will be used
        ############################################################
        self.statistics = {}.fromkeys( self.solvers )
        all_stats = Set([])
        for solver in self.solvers:
            self.statistics[solver] = (self.inputStatsSet & self.stat_keyword[solver])
            all_stats = (all_stats | self.stat_keyword[solver])
        for st in all_stats:
            if not self.precision.has_key(st):
                self.precision[st] = 2
        if not self.precision.has_key('PRD'):
            self.precision['PRD'] = 4
        if not self.precision.has_key('sig'):
            self.precision['sig'] = 3
    

    def printList(self):
        if self.printlatex == False:
            print '\t+------+'
            print '\t| LIST |'
            print '\t+------+'
        separator = '=='
        for i in range(self.bench_length):
            separator += '='
        for solver in self.solvers:            
            for i in range( (9+13*len(self.statistics[solver])) ):
                separator += '='

        print ' '.rjust(self.bench_length), 
        if self.printlatex == False:
            print '|',
        else:
            print '\\multirow{2}{*}{Instance}',
        for solver in self.solvers:
            justif = 4+(13*len(self.statistics[solver]))        
            if self.printlatex:
                solvername = solver.replace('_', '\_')
                print ' & \\multicolumn{', len(self.statistics[solver])+1 ,'}{c|}{', solvername, '}', 
            else:
                print ' ', solver.rjust(justif), '|',
              
        if self.printlatex:
            print '\\\\ \n',
        else:  
            print '\n', separator

        if self.printlatex == False:
            print 'instance '.rjust(self.bench_length),
        for solver in self.solvers:
            if self.printlatex:
                print '&',
            else:
                print '|', 
            for st in self.statistics[solver]:
                if self.printlatex:
                    print st, '&',
                else:
                    print st.rjust(12),
            print 'RES'.rjust(6), 
        if self.printlatex:
            print '',
        else:
            print '|' 
        if self.printlatex:
            print '\\\\ \n \\hline \n',
        else: 
            print separator
        for bench in self.ordered_bench:
            benchname = bench
            begin = benchname.rfind( 'normalized-' );
            if begin >= 0:
                benchname = benchname[begin+11:]
            dot = benchname.rfind( '.' )
            if dot >= 0:
                benchname = benchname[:dot]
            circ = benchname.rfind( '^' )
            if circ >= 0:
                benchname = benchname[:circ]
            if self.printlatex:
                benchname = benchname.replace('_', '\_')
                print benchname.ljust(self.bench_length+1), '&',
            else:
                print benchname.ljust(self.bench_length+1)+'|',
            for solver in self.solvers:
                if bench in self.solved_bench[solver]:
                    for st in self.statistics[solver]:
                        if self.stat[solver][bench].has_key( st ):
                            print pretty(self.stat[solver][bench][st], self.precision[st], 12),
                        else:
                            print ' '.rjust(12),
                        if self.printlatex:
                            print '&',
                else:
                    for st in self.statistics[solver]:
                        print ' '.rjust(12),
                if self.result[solver][bench] == 'SATISFIABLE' : print 'SAT'.rjust(6),
                elif self.result[solver][bench] == 'OPTIMAL' : print 'OPT'.rjust(6),
                elif self.result[solver][bench] == 'SUBOPTIMAL' : print 'subopt'.rjust(6),
                elif self.result[solver][bench] == 'UNSATISFIABLE' : print 'unsat'.rjust(6),
                elif self.result[solver][bench] == 'error--------' : print '<error>'.rjust(6),
                else : print '<> '.rjust(6),
                if self.printlatex == False:
                    print '|', 
            if self.printlatex:
                print '\\\\ \n',
            else:
                print ''
        if self.printlatex == False:
            print separator
            print '\n'


    def printSample(self):
        if self.printlatex == False and self.onlyaverage == False: 
            print '\t+--------+'
            print '\t| Sample |'
            print '\t+--------+'

        new_ub = []
        new_proof = []


        sample_length = 0
        for stat in self.sample_stats:
            sample_length += len(self.sample[stat])

        crossaverage = {}
        for solver in self.solvers:
            for stat in self.sample_stats:
                crossaverage[solver+stat+'PRD'] = []
                for method in self.sample[stat]:
                    crossaverage[solver+stat+method] = []
                    

        if self.printlatex == False:
            print ' '.rjust(self.bench_length), 
            separator = '=='
            for i in range(self.bench_length):
                separator += '='
            for solver in self.solvers:            
                for i in range( (2+13*sample_length) ):
                    separator += '='
            print '|',
        else:
            print '\\multirow{2}{*}{Instance}',

        for solver in self.solvers:
            justif = (13*sample_length)-3
            if self.printlatex:
                solvername = solver.replace('_', '\_')
                print '& \\multicolumn{', (sample_length) ,'}{c|}{', solvername, '}', 
            else:
                print ' ', solver.rjust(justif), '|',

        if self.printlatex:
            print '~\\\\ \n',
        else:  
            print '\n', separator

        if self.printlatex == False:
            print 'instance '.rjust(self.bench_length),

        for solver in self.solvers:
            if not self.printlatex:
                print '|', 
            for st in self.sample_stats:
                if self.printlatex:
                    print '& \\multicolumn{', len(self.sample[st]), '}{c|}{', st.title(), '}',
                else:
                    for method in self.sample[st]:
                        print (st[:5]+' '+method[:3]).rjust(12),

        if self.printlatex:
            print '~\\\\ \n \hline \n',
            for solver in self.solvers:
                for st in self.sample_stats:
                    for method in self.sample[st]:
                        print '&', method[:3],
            print '~\\\\ \n',
        else:
            print '|' 
            print separator

        bag = {}
        keys = []
        for bench in self.ordered_bench:
            if bench in self.common_solved:
                bench_key = bench[:bench.rfind('^')];
                if bag.has_key( bench_key ) == False:
                    bag[bench_key] = Set([])
                    keys.append(bench_key)
                bag[bench_key].add( bench )
        
        for bench_key in keys:
            benchname = bench_key
            begin = benchname.rfind( 'normalized-' );
            if begin >= 0:
                benchname = benchname[begin+11:]

            if self.onlyaverage == False:
                if self.printlatex:
                    bname = benchname.replace('_', '\_')
                    print bname.ljust(self.bench_length),
                else:
                    print benchname.ljust(self.bench_length),
                if not self.printlatex:
                    print '|',

            values = {}.fromkeys(self.sample_stats)
            for st in self.sample_stats:
                values[st] = {}.fromkeys(self.sample[st])
                for method in self.sample[st]:
                    values[st][method] = {}.fromkeys(self.solvers)
            if self.underlineoptimal:
                values['OPTIMAL'] = {}
                values['OPTIMAL']['max'] = {}.fromkeys(self.solvers)
                

            for solver in self.solvers:
                scope = (bag[bench_key] & self.solved_bench[solver])
                if self.underlineoptimal:
                    values['OPTIMAL']['max'][solver] = self.maximum(scope, solver, 'OPTIMAL')
                for st in self.sample_stats:
                    for method in self.sample[st]:
                        x = 0.0
                        if method == 'avg':
                            x = self.average(scope, solver, st)
                        elif method == 'min':
                            x = self.minimum(scope, solver, st)
                        elif method == 'max':
                            x = self.maximum(scope, solver, st)
                        elif method == 'med':
                            x = self.median(scope, solver, st)
                        elif method == 'dev':
                            x = self.deviation(scope, solver, st)
                        elif method[:3] == 'sig':
                            if method[4:] != solver:
                                x = self.significance(scope, solver, method[4:], st)
                            else:
                                x = 0
                        values[st][method][solver] = x

            for st in self.emph_stats:
                for method in self.sample[st]:
                    if self.emph[st] == 'min':
                        values[st][method]['best'] = min((values[st][method][solver], solver) for solver in self.solvers)[1]
                    elif self.emph[st] == 'max':
                        values[st][method]['best'] = max((values[st][method][solver], solver) for solver in self.solvers)[1]


            ubstat = 'UPPERBOUND'
            if 'OBJECTIVE' in self.emph_stats:
                ubstat = 'OBJECTIVE'


            #print benchname, ubstat, values['OPTIMAL']['max']['best'][:3], values[ubstat]['min']['best'][:3], 

            if 'OPTIMAL' in self.emph_stats:
                if values['OPTIMAL']['max']['best'][:3] == 'now' or values['OPTIMAL']['max']['best'][:3] == 'osp':
                    isin = True
                    #print 'optimal!', 
                    for solver in self.solvers:
                        if solver[:3] != 'now' and solver[:3] != 'osp' and values['OPTIMAL']['max'][solver] == values['OPTIMAL']['max'][values['OPTIMAL']['max']['best']]:
                            #print '...also', solver,
                            isin = False
                    if isin:
                        new_proof.append(benchname)
            #print ''


            if ubstat in self.emph_stats:
                if values[ubstat]['min']['best'][:3] == 'now' or values[ubstat]['min']['best'][:3] == 'osp':
                    isin = True
                    #print 'best!', 
                    for solver in self.solvers:
                        if solver[:3] != 'now' and solver[:3] != 'osp' and values[ubstat]['min'][solver] == values[ubstat]['min'][values[ubstat]['min']['best']]:
                            #print '...also', solver,
                            isin = False
                    if isin:
                        new_ub.append(benchname)
            #print ''

                

            best_objective = min(values[st][self.sample[st][0]][solver] for solver in self.solvers)
            for solver in self.solvers:
                scope = (bag[bench_key] & self.solved_bench[solver])
                for st in self.sample_stats:
                    x = 0
                    if best_objective != 0:
                        x = float(abs(values[st][self.sample[st][0]][solver]-best_objective))/float(best_objective)
                    crossaverage[solver+st+'PRD'].append( x )
                    for method in self.sample[st]:
                        crossaverage[solver+st+method].append( values[st][method][solver] )


            if self.onlyaverage == False:
                for solver in self.solvers:
                    scope = (bag[bench_key] & self.solved_bench[solver])
                    #for st in self.statistics[solver]:
                    for st in self.sample_stats:
                        for method in self.sample[st]:
                            x = values[st][method][solver]
                            underlined = False
                            bold = False
                            if self.printlatex:
                                if self.underlineoptimal == True and (st == 'OBJECTIVE' or st == 'UPPERBOUND') and method == 'min' and values['OPTIMAL']['max'][solver] == 1:
                                    underlined = True
                                if st in self.emph:
                                    if method in self.emph[st]:
                                        if solver == values[st][method]['best'] or values[st][method][values[st][method]['best']] == x:
                                            bold = True
                                print '~&~',
                            if self.printlatex and self.replace.has_key(st):
                                if x>=self.replace[st][0]:
                                    print self.replace[st][1],
                                else:
                                    print pretty(x, self.precision[st], 12, underlined, bold),
                            else:
                                print pretty(x, self.precision[st], 12, underlined, bold),
                    if not self.printlatex:
                        print '|',
            if self.onlyaverage == False: 
                if self.printlatex:
                    print '~\\\\'
                else:
                    print ''

            """
            if self.onlyaverage == False:
                if self.printlatex:
                    bname = benchname.replace('_', '\_')
                    print bname.ljust(self.bench_length),
                else:
                    print benchname.ljust(self.bench_length),
                if not self.printlatex:
                    print '|',
            if self.onlyaverage == False:
                for solver in self.solvers:
                    scope = (bag[bench_key] & self.solved_bench[solver])
                    for st in self.sample_stats:
                        for method in self.sample[st]:

                            if self.printlatex:
                                print '~&~',
                            print pretty(x, self.precision+2, 12),
                            if self.printlatex: #and best:
                                print '}',
                    if not self.printlatex:
                        print '|',
            if self.onlyaverage == False: 
                if self.printlatex:
                    print '\\\\'
                else:
                    print ''
            """


        if not self.printlatex:
            if self.onlyaverage == False: 
                print separator
            print 'average: '.ljust(self.bench_length+1)+'|',
            for solver in self.solvers:
                #for st in self.statistics[solver]:
                for st in self.sample_stats:
                    for method in self.sample[st]:
                        print pretty(sum(crossaverage[solver+st+method])/len(crossaverage[solver+st+method]), self.precision[st], 12),
                print '|',
            print '\n',
            """
            print 'PRD: '.ljust(self.bench_length+1)+'|',
            for solver in self.solvers:
                for st in self.sample_stats:
                    if len(crossaverage[solver+st+'PRD']) > 0:
                        print pretty(100*sum(crossaverage[solver+st+'PRD'])/len(crossaverage[solver+st+'PRD']), self.precision['PRD'], len(self.sample[st])*12),
                    else:
                        print '- '.rjust(len(self.sample[st])*12+1),
                print '|',
            print '\n'+separator                        
            print '\n'
            """
            print separator+'\n'
        else:
            print '\hline \n average',
            for solver in self.solvers:
                for st in self.sample_stats:
                    for method in self.sample[st]:
                        print '~&~', pretty(sum(crossaverage[solver+st+method])/len(crossaverage[solver+st+method]), self.precision[st], 12),
            print '~\\\\'

            print '\hline \n PRD',
            for solver in self.solvers:
                for st in self.sample_stats:
                    for mt in self.sample[st]:
                        print '~&~',
                if len(crossaverage[solver+st+'PRD']) > 0:
                    print pretty(100*sum(crossaverage[solver+st+'PRD'])/len(crossaverage[solver+st+'PRD']), self.precision['PRD'], 12),
            print '~\\\\'


            if 'OPTIMAL' in self.emph_stats and 'UPPERBOUND' in self.emph_stats:
                seq1 = []
                seq2 = []
                for solver in self.solvers:
                    #print crossaverage[solver+'UPPERBOUND'+'PRD']

                    seq1.append(sum(crossaverage[solver+'OPTIMALmax'])/len(crossaverage[solver+'OPTIMALmax']))
                    seq2.append(sum(crossaverage[solver+'UPPERBOUND'+'PRD'])/len(crossaverage[solver+'UPPERBOUND'+'PRD']))
                best1 = max(seq1)
                best2 = min(seq2)

                print  '\n',
                for i in range(len(seq1)):
                    bf = False
                    if seq1[i] == best1:
                        bf = True
                    print '~&~', pretty(seq1[i], self.precision['OPTIMAL'], 12, False, bf),

                    bf = False
                    if seq2[i] == best2:
                        bf = True
                    print '~&~', pretty(seq2[i], self.precision['PRD'], 12, False, bf),
                print '~\\\\'


            print '\nnew proofs:', len(new_proof)
            for bench in new_proof:
                print bench

            print '\nnew bests:', len(new_ub)
            for bench in new_ub:
                print bench




    def printTotal(self, thresh):
        if self.printlatex == False and self.onlyaverage == False: 
            print '\t+-------+'
            print '\t| Total |'
            print '\t+-------+'

        sample_length = 0
        for stat in self.sample_stats:
            sample_length += len(self.sample[stat])

        values = {}
        for stat in self.sample_stats:
            values[stat] = []


        bag = {}
        keys = []
        for bench in self.ordered_bench:
            if bench in self.common_solved:
                bench_key = bench[:bench.rfind('^')];
                if bag.has_key( bench_key ) == False:
                    bag[bench_key] = Set([])
                    keys.append(bench_key)
                bag[bench_key].add( bench )
        
        for bench_key in keys:
            for st in self.sample_stats:
                for solver in self.solvers:
                    scope = (bag[bench_key] & self.solved_bench[solver])
                    for bench in scope:
                        if float(self.stat[solver][bench]['TIME']) < thresh:
                            values[st].append( float(self.stat[solver][bench][st]) )

        
        separator = '========';
        for st in self.sample_stats:
            separator += '|='
            for method in self.sample[st]:
                separator += '============='

        print ' total ',
       
        for st in self.sample_stats:
            print '|', 
            for method in self.sample[st]:
                print (st[:5]+' '+method[:3]).rjust(12),
        print '\n', separator


        print str(len(values[self.sample_stats[0]])).rjust(5), ' ',
        for st in self.sample_stats:
            print '|', 
            for method in self.sample[st]:
                x = 0
                if method == 'avg':
                    x = self._average(values[st])
                elif method == 'min':
                    x = self._minimum(values[st])
                elif method == 'max':
                    x = self._maximum(values[st])
                elif method == 'med':
                    x = self._median(values[st])
                elif method == 'dev':
                    x = self._deviation(values[st])
                print pretty(x, self.precision[st], 12),
        print '\n', separator        


    def printLatexAverages( self, benches, solver ):
        usd = False
        for st in self.inputStats:
            if st in self.statistics[solver]:
                if usd == True:
                    print '~&~', pretty(self.average( benches, solver, st ), self.precision[st], 12),#+'\\%',
                    usd = False
                else:
                    print '~&~', pretty(self.average( benches, solver, st ), self.precision[st], 12),
                    usd = True
            else:
                print '~&~', ' '.rjust(12),


    def printAverages( self, benches, solver ):
        for st in self.inputStats:
            if st in self.statistics[solver]:
                print '|', pretty(self.average( benches, solver, st ), self.precision[st], 12),
            else:
                print '|', ' '.rjust(12),

    def getAverages( self, benches, solver ):
        line = ''
        for st in self.inputStats:
            if st in self.statistics[solver]:
                line += ('| '+pretty(self.average( benches, solver, st ), self.precision[st], 12)+' ')
            else:
                line += ('| '+' '.rjust(12)+' ')
        return line

    def printSummary(self):

        if self.printlatex == False:
            ordered_solvers = []
            for solver in self.solvers:
                ordered_solvers.append( '='+(str(len(self.common_bench) - len(self.solved_bench[solver]))).rjust(6)+'|'+solver )
            ordered_solvers.sort();
            for i in range(len(ordered_solvers)):
                ordered_solvers[i] = ordered_solvers[i][ordered_solvers[i].find('|')+1:]

            print '\t+-------------------------------+'
            print '\t| SUMMARY SELF -', str(len(self.common_bench)).rjust(4), 'instances |'
            print '\t+-------------------------------+'

            if len(self.common_bench) > 0 :

                print ' '.rjust(self.solver_length), '|', '#solved'.rjust(12)+'     |', '#errors'.rjust(7),
                for st in self.inputStats:
                    print '|', st.rjust(12),
                print '|'
                separator = '==============================='#=============================='
                for i in range(self.solver_length):
                    separator += '='
                for st in self.inputStats:
                    separator += '==============='
                print separator
                for solver in ordered_solvers: #self.solvers:
                    print solver.rjust(self.solver_length),
                    print '|', str(len(self.solved_bench[solver])).rjust(4), '/',  str(len(self.common_bench)).ljust(4), (str((len(self.solved_bench[solver]) * 100) / len(self.common_bench))+'%').rjust(4), '|'+str(len(self.error_bench[solver])).rjust(8), 
                    self.printAverages( self.solved_bench[solver], solver )
                    print '|'
                print separator

            print '\n'
            print '\t+---------------------------------+'
            print '\t| SUMMARY COMMON -', str(len(self.common_solved)).rjust(4), 'instances |'
            print '\t+---------------------------------+'


            sigmethod = None
            for method in self.sample:
                if method[:3] == 'sig':
                    sigmethod = method[4:]
                    break

            if len(self.common_solved) > 0 :

                for i in range(len(ordered_solvers)):
                    the_solver = ordered_solvers[i]
                    ordered_solvers[i] = (self.getAverages(self.common_solved, ordered_solvers[i])+'$'+ordered_solvers[i])
                    if sigmethod != None:
                        sig = pretty(self.significance(self.common_solved, the_solver, sigmethod, self.inputStats[0]), self.precision['sig'], 5)
                        ordered_solvers[i] += ('@'+sig)
                        

                print ' '.rjust(self.solver_length), 
                for st in self.inputStats:
                    if sigmethod == None:
                        print '|', st.rjust(12),
                    else:
                        print '|', st.rjust(18),
                print ' |'
                separator = '==='
                for i in range(self.solver_length):
                    separator += '='
                for st in self.inputStats:
                    separator += '==============='
                    if sigmethod != None:
                        separator += '======'
                print separator


                ordered_solvers.sort()
                for solver in ordered_solvers: #self.solvers:
                    sig = ''
                    sol = ''
                    sta = ''
                    if sigmethod:
                        sig = solver[solver.find('@')+1:]
                        sol = solver[solver.find('$')+1:solver.find('@')].rjust(self.solver_length)
                    else:
                        sol = solver[solver.find('$')+1:].rjust(self.solver_length)
                    sta = solver[:solver.find('$')]
                    print sol, sta, 
                    if sigmethod != None:
                        print sig,
                    print '|'

                    #print solver[:solver.find('$')],
                    #print solver.rjust(30),
                    #self.printAverages( self.common_solved, solver )
                    #print '|'
                print separator

            else:

                print 'No data'
            print '\n'

        else:
            
            if len(self.common_solved) > 0 :

                for solver in self.solvers:
                    self.printLatexAverages(self.common_solved, solver)
                print '\\\\\n'


    def printEasyInstances(self):
        for bench in self.common_solved:
  #          print bench,
            max = 0;
            for solver in self.solvers:
                t = float(self.stat[solver][bench]['TIME'])
 #               print t,
                if t > max:
                    max = t
#            print '=>', max
            if max <= self.printeasy:
                print self.raw_data['pth'][self.solvers[0]][bench][0]

    def printHardInstances(self):
        hard = self.common_bench.copy()
        for bench in self.common_solved:
  #          print bench,
            max = 0;
            for solver in self.solvers:
                t = float(self.stat[solver][bench]['TIME'])
 #               print t,
                if t > max:
                    max = t
#            print '=>', max
            if max <= self.printhard:
                hard.remove(bench)
        for bench in hard:
            print self.raw_data['pth'][self.solvers[0]][bench][0]

    def printSolvedInstances(self):
        for bench in self.solved_by_one:
            print self.raw_data['pth'][self.solvers[0]][bench][0]

    def printErrors(self):
        for solver in self.showerror:
            #subprocess.Popen( ['rm', '-rf', 'results/'+solver+'.errorfiles/'] )
            #subprocess.Popen( ['mkdir', 'results/'+solver+'.errorfiles/'] )
            for bench in self.error_bench[solver]:
                #subprocess.Popen( ['cp', self.raw_data['pth'][solver][bench][0], 'results/'+solver+'.errorfiles/'] )
                #for k in self.raw_data['pth'].keys():
                #    print k
                #for k in self.raw_data['pth'][solver].keys():
                #    print k
                print self.raw_data['pth'][solver][bench][0]
                

    def printDelta(self):
        for solver in self.delta:
            subprocess.Popen( ['rm', '-rf', 'results/'+solver+'.deltafiles/'] )
            subprocess.Popen( ['mkdir', 'results/'+solver+'.deltafiles/'] )
            for bench in (self.benchmark[reference] & self.benchmark[solver]):
                if (result[reference][bench][-11:] == "SATISFIABLE") & (result[solver][bench][-11:] != "SATISFIABLE"):
                    subprocess.Popen( ['cp', self.raw_data['pth'][solver][bench][0], 'results/'+solver+'.deltafiles/'] )
                    print self.raw_data['pth'][solver][bench][0]



    def verifyCertificate(self):
        for solver in self.solvers:
            print solver
            for bench in self.solved_bench[solver]:
                if self.result[solver][bench] == 'SATISFIABLE':
                    solution = 's SATISFIABLE\nv '+self.raw_data['sol'][solver][bench]
                    print '\t+CHECK', bench.ljust(50),
                    if checkSolution( solution, self.raw_data['pth'][solver][bench][0] ):
                        print 'OK'
                    else:
                        print ' /!\\ WARNING,', solver, 'gives a wrong certificate for', bench, ' /!\\'
                        self.wrong_certificate[solver].append(bench)

    def verifyProof(self):
        for solver in self.solvers:
            print solver
            for bench in self.solved_bench[solver]:
                if self.result[solver][bench] == 'UNSATISFIABLE':
                    print '\t+CHECK', bench.ljust(50),
                    solvedByOther = False
                    for other in self.solvers:
                        if (self.result[other][bench] == 'SATISFIABLE') & (bench not in self.wrong_certificate[other]):
                            print ' /!\\ WARNING,', solver, 'gives a wrong proof for', bench, ' /!\\'
                            self.wrong_proof[solver].append(bench)
                        else:
                            print 'OK'


    def verifyCertificate(self):
        self.jobs = []
        self.verified = {}.fromkeys( self.solvers )
        for solver in self.solvers:
            self.verified[solver] = Set([])
            process = subprocess.Popen( ['ls', 'results/verified_'+self.keys[0]+'_'+solver+'.ver'], stderr=subprocess.PIPE)
            errput = process.stderr.readline()[:17]
            if errput != 'ls: cannot access':
                for bench in open('results/verified_'+self.keys[0]+'_'+solver+'.ver'):
                    self.verified[solver].add( bench[:-1] )
            unverified = (self.solved_bench[solver] - self.verified[solver])
            for bench in unverified:
                if self.result[solver][bench] == 'SATISFIABLE':
                    self.jobs.append( [solver, 's SATISFIABLE\nv '+self.raw_data['sol'][solver][bench], self.raw_data['pth'][solver][bench][0]] )
        checkers = []
        global GLOBAL_END
        print GLOBAL_END
        if (GLOBAL_END == 0) | (GLOBAL_END > len(self.jobs)):
            GLOBAL_END = len(self.jobs)
        print GLOBAL_END
        for i in range(self.nThreads):
            checkers.append( solutionChecker(self.jobs, self.verified, self.wrong_certificate, i, self) )
            #checkers.append( solutionChecker(self.jobs, i) )
        for checker in checkers:
            checker.start()
        #while 1:
        #    time.sleep ( 10 )
        #    for checker in checkers:
        #        print checker.isAlive()
        


    def printWrongCertificate(self):
        for solver in self.solvers:            
            if solver not in self.wrong_certificate.keys():
                print solver, 'gives no wrong certificate :-)'
            elif len(self.wrong_certificate[solver]) == 0:
                print solver, 'gives no wrong certificate :-)'
        for solver in self.wrong_certificate.keys():
            for bench in self.wrong_certificate[solver]:
                print solver, 'gives a wrong certificate for', bench[bench.rfind('/')+1:]
        for solver in self.verified.keys():
            verif_file = open( 'results/verified_'+self.keys[0]+'_'+solver+'.ver', 'w' )
            for bench in self.verified[solver]:
                verif_file.write( bench+'\n' )


    def run(self):
        if len(self.sample_stats) > 0:
            self.printSample()
        if self.list == True:
            self.printList()
        if self.summary == True:
            self.printSummary()
        if len(self.showerror) > 0:
            self.printErrors()
        if len(self.delta) > 0:    
            self.printDelta()
        if self.verif:
            self.verifyCertificate()
        if self.printsolved:
            self.printSolvedInstances()
        if self.printeasy > 0:
            self.printEasyInstances()
        if self.printhard > 0:
            self.printHardInstances()


        #self.printTotal(1200);


        #for solver in self.wrong_certificate.keys():
        #    for bench in self.wrong_certificate[solver]:
        #        print solver, 'gives a wrong certificate for', bench
        #for solver in self.wrong_proof.keys():
        #    for bench in self.wrong_proof[solver]:
        #        print solver, 'gives a wrong proof for', bench


parser = Parser(sys.argv[1:])
#parser.getparameters( sys.argv[1:] )
if len(parser.latex_args) == 0:    
    parser.initialize()
    parser.run()
else:

    print "\\documentclass{llncs}"
    print "\\usepackage{amssymb}"
    print "\\usepackage{latexsym}"
    print "\\usepackage{epsfig}"
    print "\\usepackage[dvips]{color}"
    print "\\usepackage{multirow}"
    print "\\begin{document}"
    print "\\newlength{\\halftextwidth}"
    print "\\setlength{\\halftextwidth}{0.47\\textwidth}"
    print "\\def\\halffigsize{2.75in}"
    print "\\def\\negvspace{0in}"
    print "\\def\\posvspace{0em}"
    print "\\def\\msetl{\\mbox{$\\{ \\! \\{$}}"
    print "\\def\\msetr{\\mbox{$\\} \\! \\}$}}"
    print "\\addtolength{\\textwidth}{0.5in}"
    print "\\addtolength{\\topmargin}{-0.5in}"
    print "\\addtolength{\\textheight}{0.5in}"

    for args in parser.latex_args:

        new_parser = Parser(args)
        new_parser.printlatex = True
        #new_parser.keys = []
        #new_parser.keys.append(key)
        new_parser.initialize()

        print "\\begin{table*}[t]"
        #print "\\begin{scriptsize}"
        print "\\caption{",
        for stat in new_parser.inputStats:
            print stat.title(),
        print ": ",
        for solver in new_parser.solvers:
            print solver.replace('_','-'), 
        print ".}"
        print "\\begin{center}"
        print "\\begin{tabular}{|l||",

        if len(new_parser.sample_stats) > 0:
            for i in range(len(new_parser.solvers)):
                for st in new_parser.sample_stats:
                    for method in new_parser.sample[st]:
                        print 'r',
                    print '|',
                print '|',
        else:
            for i in range(len(new_parser.solvers)):
                for j in range(len(new_parser.inputStats)+1):
                    print 'r',
                print '|',
        print "}"
        print "\\hline \\hline"

        '''
        the_keys = parser.keys
        last_key = the_keys[0]
        for key in the_keys:

            new_parser = Parser(sys.argv[1:])
            new_parser.keys = []
            new_parser.keys.append(key)
            new_parser.initialize()
            new_parser.run()
        '''

        new_parser.run()

        print "\\hline\\hline"
        print "\\end{tabular}"
        print "\\end{center}"
        #print "\\end{scriptsize}"
        print "\\end{table*}"


    print "\\end{document}"


'''
    print "\\multirow{2}{*}{Benchmark} & Training Set ",
    for solver in parser.solvers:
        print " & \\multicolumn{2}{c|}{", solver.replace('_','-'), "}",
    print "\\\\"
    print "& prop. ",
    for solver in parser.solvers:
        for stat in parser.inputStats:
            print "&", stat[:5],
    print "\\\\"
    print "\\hline \\hline"
'''

'''
        if last_key[:last_key.rfind('.')] != key[:key.rfind('.')]:
            print "\\hline"
            last_key = key
        print key[3:key.rfind('.')-1], '&', key[key.rfind('.')-1:],
'''
