#! /usr/bin/env python

import sys
from Numberjack import *


###############################################
#######   Class JSP: problem instance   #######
###############################################
class JSP:
    def __init__(self,data_file):

        stream = open(data_file)
        (n, m) = stream.readline().split()[:2]
        self.nJobs = int(n)
        self.nMachines = int(m)

        stream.readline()
        self.job = []
        self.machine = [[] for i in range(self.nMachines)]
        self.m = [[None]*self.nMachines for i in range(self.nJobs)]

        for i in range(self.nJobs):
            self.job.append([int(elt) for elt in (stream.readline()[:-1]).split()])

        stream.readline()
        for i in range(self.nJobs):
            machines = (stream.readline()[:-1]).split()
            for j in range(len(machines)):
                self.machine[int(machines[j])-1].append((i,j))
                self.m[i][j] = (int(machines[j])-1)

    def __str__(self):
        return( '\n'+str(self.job)+'\n\n'+str(self.machine)+'\n\n'+str(self.m)+'\n' )

    def lower_bound(self):
        longest_job = max([sum(job) for job in self.job])
        longest_machine = max([sum([self.job[i][j] for (i,j) in mac]) for mac in self.machine])
        return max([longest_job, longest_machine])

    def upper_bound(self):
        M_job = [0]*self.nJobs
        M_machine = [0]*self.nMachines

        for i in range(self.nMachines):
            for j in range(self.nJobs):
                start_time = max(M_job[j], M_machine[self.m[j][i]])
                M_job[j] = start_time+self.job[j][i]
                M_machine[self.m[j][i]] = start_time+self.job[j][i]
                
        return max(max(M_job), max(M_machine))


###############################################
###### Class JSP_Model: constraint model ######
###############################################
class JSP_Model(Model):

    def __init__(self,jsp):
        Model.__init__(self)

        C_max = jsp.upper_bound()

        Tasks = Matrix([[Task(C_max, p) for p in job] for job in jsp.job])
        Machines = [UnaryResource([Tasks[m] for m in machine]) for machine in jsp.machine]

        for task in Tasks.row: 
            self += [task[i] < task[i+1] for i in range(jsp.nMachines-1)]
        self += Machines

        self.sequence = sum(Machines, [])
        self.tasks = Tasks.flat
        self.Jobs = Tasks
    

###############################################
##############  function solve   ##############
###############################################
def dichotomic_step(model, solver, C_max, best_solution, verb, cutoff):
    solver.save()
    for task in model.tasks:
        solver.post(task < C_max)

    if best_solution != None: solver.guide(best_solution)
    solver.setNodeLimit(cutoff)
    solver.solveAndRestart(GEOMETRIC, 256, 1.3)

    outcome = (None, None, C_max)
    if solver.is_sat():
        outcome = (True, solver.get_solution(), max([task.get_min() + task.duration for task in model.tasks]))
        if verb>0: print '   SAT :-)', 
    elif solver.is_unsat():
        outcome = (False, None, C_max)
        if verb>0: print ' UNSAT :-)', 
    else:
        if verb>0: print ' ABORT :-(', 

    if verb>0: print str(solver.getTime()).rjust(8), 's', str(solver.getNodes()).rjust(10), 'nds'
 
    solver.reset()

    solver.undo()

    return outcome


###############################################
########  function dichotomic search   ########
###############################################
def dichotomic_search(model, solver, max_infeasible, min_feasible, verb, cutoff):
    if verb>0: print 'start dichotmic search', cutoff

    lb = max_infeasible
    ub = min_feasible
    best_solution = None

    while lb+1 < ub:
        C_max = int((lb + ub) / 2)
        
        if verb>0: print 'c   current bounds:', ('['+str(lb+1)+'..'+str(ub)+']').rjust(16), ' solve', str(C_max).ljust(6),

        (feasible, solution, C_max) = dichotomic_step(model, solver, C_max, best_solution, verb, cutoff)

        if feasible:
            ub = C_max
            best_solution = solution
        else:
            lb = C_max
            if feasible is not None: 
                max_infeasible = C_max

    min_feasible = ub
    return (max_infeasible, min_feasible, best_solution)


###############################################
########   function branch and bound   ########
###############################################
def branch_and_bound(model, lib, max_infeasible, min_feasible, verb, best=None):
    C_max = Variable(max_infeasible+1, min_feasible)
    for task in model.tasks: 
        model.add(task < C_max)
    model.add( Minimise(C_max) )

    if lib == 'Mistral': solver = model.load( lib, model.sequence )
    else: solver = model.load( lib )

    solver.setVerbosity(verb-1)
    solver.setHeuristic('Scheduling')
    if best is not None: solver.guide(best)

    
    outcome = (max_infeasible+1, None)
    print solver.solve() #solver.solveAndRestart():

    if solver.is_sat():
        best = solver.get_solution()

    if solver.is_opt():
        if verb>0: print 'c   Found optimal solution:', C_max.solution()
        outcome = (C_max.get_value(), C_max.get_value(), best)
    else:
        if verb>0: print 'c   Best C_max:', C_max.solution()
        outcome = (max_infeasible+1, C_max.get_value(), best)
    return outcome


###############################################
###########  main solver function   ###########
###############################################
def solve(param):
    jsp = JSP(param['data'])
    lib = param['solver']
    verb = param['verbose']

    model = JSP_Model(jsp)
    
    if lib == 'Mistral': solver = model.load(lib, model.sequence)
    else: solver = model.load(lib)

    solver.setHeuristic('Scheduling', 'Promise')
    solver.setVerbosity(param['verbose']-1)

    (lb, ub) = (jsp.lower_bound()-1, jsp.upper_bound())

    (lb, ub, best) = dichotomic_search(model, solver, lb, ub, verb, param['tcutoff'])
    if verb>0: print 'start branch & bound in ['+str(lb)+'..'+str(ub)+']'
    if lb+1 < ub: (lb, ub, best) = branch_and_bound(model, lib, lb, ub, verb, best)

    ## finalize the solution (tasks)
    solver.reset()
    if lib == 'Mistral': 
        for disjunct in model.sequence: solver.post(disjunct == best[disjunct])
        solver.propagate()
        for task in model.tasks:
            solver.post(task == task.get_min())
            solver.propagate()
    best = solver.get_solution()


    schedule = [[-1]*ub for job in jsp.job]
    index = 0
    for machine in jsp.machine:
        index += 1
        for m in machine:
            start = model.Jobs[m].get_value()
            for i in range(model.Jobs[m].duration):
                schedule[m[0]][start+i] = index

    out = ''
    if solver.is_sat():
        out = str(schedule)
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out    


    if param['print'] == 'yes':
        ###############################################
        ############# Output (Matplotlib) #############
        ###############################################
        print '\n display schedule'

        width = 60
        print_schedule = []
        for row in schedule:
            print_schedule.extend([row]*width)

        import pylab
        pylab.yticks( pylab.arange(int(width/2), width*(len(jsp.job)+1), width), ['job'+str(i+1) for i in range(len(jsp.job))] )
        cmap = pylab.cm.get_cmap('jet', len(jsp.machine)+1)
        cmap.set_under(color='w')
        im1 = pylab.imshow(print_schedule, cmap=cmap, interpolation='nearest', vmin=0)
        #pylab.colorbar()
        pylab.show()


solvers = ['Mistral', 'MiniSat']
default = {'solver':'Mistral', 'data':'data/tiny_jsp.txt', 'print':'no', 'verbose':1, 'tcutoff':3}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)




