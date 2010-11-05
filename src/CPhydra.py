#! /usr/bin/env python

'''

CPhydra module in python

-eomahony

'''
import sys
from Numberjack import *
import Mistral
import cPickle
import os
import networkx

from random import *

epsilon = 0.0001

class Case:
    def __init__(self, name='no_name', f={}, b={}):
        self.name = name
        self.features = f
        self.beliefs = b

    def compare(self, case):
        common_keys = (frozenset(self.features.keys()) & 
                       frozenset(case.features.keys()))

        distance = 0
        for k in common_keys:
            d = (self.features[k] - case.features[k])
            distance += (d*d)

        return 1+distance        

    def expected_solve_time(self,solver,default=1800.0):
        if self.beliefs.has_key(solver): return self.beliefs[solver]
        else: return default

    def __str__(self):
        sep = '+'+'='*58+'+'
        return_str = sep+'\n|'+('case name : ').rjust(29)+self.name.ljust(29)+'|\n|'
        sorted_keys = self.features.keys()
        sorted_keys.sort()
        for k in sorted_keys:
            return_str += ((k+' : ').rjust(29)+str(self.features[k]).ljust(29)+'|\n|')
        return_str += (' beliefs:').ljust(58)+'|\n|'
        sorted_keys = self.beliefs.keys()
        sorted_keys.sort()
        for k in sorted_keys:
            return_str += ((k+' : ').rjust(29)+str(self.beliefs[k]).ljust(29)+'|\n|')
        return return_str[:-1]+sep+'\n'

class Casebase(list):
    '''
    This will store the casebase
    '''
    def __init__(self, cases=[]):
        list.__init__(self, cases)

    def __str__(self):
        return '\n'.join(map(str, self))

    def get_probabilities(self, case, solver, k=10):
        global epsilon
        W = 0
        data = {}
        cases = self.get_similar(case, k)
        for d,c in cases:
            if d>1+epsilon:
                t = c.expected_solve_time(solver)
                if data.has_key(t): data[t] += 1/d 
                else: data[t] = 1/d
                W += 1/d
        runtimes = data.keys()
        runtimes.sort()
        proba = [(0,0)]
        p = 0
        for r in runtimes: 
            p += data[r]
            proba.append((r, p/W))
        return proba

    def get_similar(self, case, k=10):
        comparisons = [(case.compare(c), c) for c in self]
        comparisons.sort()
        #return [case for (_, case) in comparisons[:k]]
        return comparisons[:k]
    
    def save(self,file):
        cPickle.dump(self, open(file, 'w'))

    def add(self,case):
        '''
        This method adds a case to the base,
        or updates the current case if a similar exists 
        '''
        (distance,most_similar) = min([(case.compare(c), c) for c in self])
        if distance <= epsilon: # we assume that it is the same instance
            for k in case.keys(): most_similar[k] = case[k]
        else: self.append(case)

def load_casebase(file): 
    return Casebase(cPickle.load(open(file)))
    
def load_case(file):
    parser = Mistral.Solver(Model())
    #print file
    parser.load_xml(file)
    return Case(file[file.rfind('/')+1:], 
                parser.get_static_features(),
                parser.get_dynamic_features())


def load_graph(file):
    G = networkx.Graph()

    parser = Mistral.Solver(Model())
    parser.load_xml(file)
    parser.extract_graph()
    for x in range(parser.numNodes()):
        for y in parser.get_neighbors(x):
            G.add_edge(x,y)

    return G
    '''
    pos=networkx.graphviz_layout(G,prog="neato")
    networkx.draw(G,
                  pos,
                  node_size=40,
                  node_color=1,
                  vmin=0.0,
                  vmax=1.0,
                  with_labels=False
                  )
    '''
    
def classify(case):
    if case.features['max_arity'] <= 2:
        if case.features['percent_ext'] > 0:
            return "binaryExt"
        else:
            return "binaryInt"
    elif case.features['percent_global'] > 0:
        return "global"
    elif case.features['percent_ext'] > 0:
        return "nAryExt"
    else:
        return "nAryInt"
    
def schedule(cases, solvers, total_time):
    def weight(case):
        d = case[0]
        if d >= 1000: return 1
        elif d <= 0.01: return 100000
        else: return int(1000/d)  

    model = Model()

    allocated_times = VarArray(len(solvers), 0, total_time)
    prob_solved = VarArray(len(cases))

    model.add(Sum(allocated_times) == total_time)

    for prob, case in zip(prob_solved, cases):
        expected_times = [case[1].expected_solve_time(solver) for solver in solvers]
        model.add(prob == (Sum([stime >= int(etime+1) for stime, etime in 
                                zip(allocated_times, expected_times)]) >= 1))

    model.add(Maximise(Sum(prob_solved, [weight(case) for case in cases])))

    solver = Mistral.Solver(model, allocated_times)
    solver.setHeuristic('MinDomain','RandomSplit',3)

    solver.solveAndRestart()

    return dict(zip(solvers, map(value, allocated_times)))
    
def split_casebase_in_half(casebase, seed_val=0):
    cat_count = {"global":[], "binaryExt":[],
                "binaryInt":[], "nAryInt":[],
                "nAryExt":[]}
    
    for case in casebase:
        cat_count[classify(case)].append(case)
    
    half1 = []
    half2 = []
    
    seed(seed_val)
    for list in cat_count.values():
        shuffle(list)
        half1.extend(list[:len(list)/2])
        half2.extend(list[len(list)/2:])
    return (half1, half2)
    
def test_grid_schedule(casebase, solvers, total_time, machines):
    import math
    half1, half2 = split_casebase_in_half(casebase)
    
    print len(half1)
    print total_time
    for case in half1:
        beliefs = case.beliefs
        min_time = min(beliefs.items(), key=lambda x: x[1])
        print int(math.ceil(min_time[1]))
        
    #sch = schedule_grid_computation(half1, solvers, total_time, machines, Casebase(half2))
    
def scheudle_grid(case_dir, solvers, total_time, machines, casebase):
    # 1. Parse features vectors for every case
    cases = []
    for file in os.listdir(case_dir):
        parser = Mistral.Solver()
        parser.load_xml(file)
        cases.append(Case(file[file.rfind('/')+1:], 
                parser.get_static_features(),
                parser.get_dynamic_features()))
    return schedule_grid_computation(cases, solvers, total_time, machines, casebase)
    
def schedule_grid_computation(cases, solvers, total_time, machines, casebase):
   
    total_time *= machines
    
    # 2. Now retrieve similar cases for each of the cases
    nn_cases = [(case, casebase.get_similar(case)) for case in cases]
    
    # 3. Scheduling time
    model = Model()
    
    # Set up the SolTask variables
    tasks = []
    for solver in solvers:
        for case in nn_cases:
            #print map(classify, (map(lambda x: x[1], case[1])))
            tasks.append(SolTask(solver, case[0], case[1], total_time))
            
    # Add constraint on end times
    model.add([task.end_time <= total_time for task in tasks])
            
    # Add the cumulative (maybe)
    #for t1, t2 in pair_of(tasks):
    #    model.add(t1.start_time + t1.duration < t2.start_time or
    #              t2.start_time + t2.duration < t1.start_time)
    model.add(Sum([task.duration for task in tasks]) <= total_time)

    # Add the objective function
    model.add(Maximise(Sum([task.power for task in tasks])))
    
    #print model
    solver = Mistral.Solver(model)
    solver.setVerbosity(1)
    solver.setTimeLimit(2*60)
    solver.setHeuristic('MinDomain','RandomSplit',3)
    solver.solveAndRestart()
    
    print "".join(map(str, tasks))
        
class SolTask:
    '''
    Class to handle scheduling across multiple machines
    '''
    def __init__(self, solver, problem, scases, max_end_time):
        self.solver = solver
        self.problem = problem
        self.start_time = Variable(0, max_end_time)
        self.duration = Variable(0, max_end_time)
        self.end_time = self.start_time + self.duration
        self.power = Sum([self.duration >=
                          int(case[1].expected_solve_time(solver))+1
                          # Hacked to an int to keep mistral happy
                          for case in scases])
        
    def __str__(self):
        return (self.solver + " on " + self.problem.name + " from " +
                str(self.start_time) + " to " + str(self.end_time) +
                " ( " + str(self.duration) + " )\n")


CBR = load_casebase('./tools/casebase.pkl')
#print CBR[:-1]

all_easy = 0
easy_nothard = 0
easy_hard = 0
all_middle = 0
noteasy_hard = 0
all_hard = 0

for case in CBR[:-1]:
    print case.beliefs['mistral'], case.beliefs['choco'], case.beliefs['abscon']
    mistral = 2
    choco = 2
    abscon = 2
    if case.beliefs['mistral'] <= 10: mistral = 0
    elif case.beliefs['mistral'] <= 1799: mistral = 1
    if case.beliefs['choco'] <= 10: choco = 0
    elif case.beliefs['choco'] <= 1799: choco = 1
    if case.beliefs['abscon'] <= 10: abscon = 0
    elif case.beliefs['abscon'] <= 1799: abscon = 1
    
    if max([mistral,choco,abscon]) == 0: all_easy += 1
    elif max([mistral,choco,abscon]) == 1:
        if min([mistral,choco,abscon]) == 0: easy_nothard +=1
        else: all_middle += 1
    elif min([mistral,choco,abscon]) == 0: easy_hard +=1
    elif min([mistral,choco,abscon]) == 1: noteasy_hard +=1
    else: all_hard += 1

print 1, all_easy
print 2, easy_nothard
print 3, all_middle
print 4, easy_hard
print 5, noteasy_hard
print 6, all_hard
    
"""
if __name__ == "__main__":

    params = input({'-extract-graph':'no',
                    '-allocate-time':'no', 
                    '-update':'no', 
                    '-expected-time':'no', 
                    '-expected-distribution':'no', 
                    'solvers':['abscon','choco','mistral'], 
                    'file':'no', 'dir':'.', 'features':[], 
                    'values':[], 
                    'name':'target', 'time':1800,
                    'cputimes':[1800.0, 1800.0, 1800.0],
                    'prob-dir': 'none',
                    'machines': '1',
                    '-grid-allocate':'no',
                    'get-features':'no'})

    cbr_file = params['dir']+'/tools/casebase.pkl'
    CBR = load_casebase(cbr_file)

    target = None
    target_file = None
    if params['file'][0] != '/':
        target_file = params['dir']+'/'+params['file']
    else:
        target_file = params['file']

    solvers = []
    if type(params['solvers']) == str:
        solvers.append(params['solvers'])
    else:
        solvers.extend(params['solvers'])

    cputimes = []
    if type(params['cputimes']) == float:
        cputimes.append(params['cputimes'])
    else:
        cputimes.extend(params['cputimes'])

    if params['file'] != 'no':
        target = load_case(target_file)
    elif params['prob-dir'] != 'none':
        prob_dir = params['prob-dir']
    elif len(params['features']) > 0:
        print "Feats:", len(params['features'])
        print "Values:",len(params['values'])
        if len(params['features']) != len(params['values']):
            print 'warning: features names and values do not align'
        if len(solvers) != len(cputimes):
            print 'warning: solvers names and cputimes do not align'
        target = Case(params['name'], 
                      dict(zip(params['features'],map(float,params['values']))),
                      dict(zip(solvers, map(float, cputimes))))
    else:
        print 'need an xml file or a set of features/values - exiting'
        exit(1)

    if params['-allocate-time'] == 'yes':
        similar_cases = CBR.get_similar(target)
        print schedule(similar_cases, solvers, params['time'])
    elif params['-update'] == 'yes':
        CBR.add(target)
        CBR.save(cbr_file)
    elif params['-expected-time'] == 'yes':
        '''
        for solver in solvers:
            distribution = CBR.get_probabilities(target, solver, 30)
            (ot, op) = distribution[0] 
            for t,p in distribution[1:]:
                print '[',ot,'..',t,'[: ', op
                ot = t
                op = p
            print '[',ot,'.. inf [: ', op

        '''
        similar_cases = CBR.get_similar(target)
        print dict(zip(solvers,[(sum([c[1].expected_solve_time(s) 
                                      for c in similar_cases])/len(similar_cases)) 
                                for s in solvers]))
        
    elif params['get-features'] == 'yes':
        print " ".join(map(str, target.features.values()))

    elif params['-expected-distribution'] == 'yes':
        
        for solver in solvers:
            distribution = CBR.get_probabilities(target, solver, 30)
            (ot, op) = distribution[0] 
            for t,p in distribution[1:]:
                print '[',ot,'..',t,'[: ', op
                ot = t
                op = p
            print '[',ot,'.. inf [: ', op
    
    elif params['-grid-allocate'] == 'yes':
        # Do clever stuff
        #sch = scheudle_grid(prob_dir, solvers, int(params['time']),
        #                    int(params['machines']), CBR)
        test_grid_schedule(CBR, solvers, int(params['time']), int(params['machines']))
    elif params['-extract-graph'] == 'yes':
        G = load_graph(target_file)
        print G.edges()
    else:
        print 'use one of --alocate-time, --update or --expected-time\n(or -help)'
"""

