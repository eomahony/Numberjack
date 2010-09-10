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
        ## matrix of job's durations
        self.job = []
        ## indices of the jobs organized by machines
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



def get_model(jsp):
    ###############################################
    ##############      Model        ##############
    ###############################################
    lb = jsp.lower_bound()
    ub = jsp.upper_bound()

    C_max = Variable(lb, ub, 'C_max')
    Jobs = Matrix([[Task(ub, p) for p in job] for job in jsp.job])

    model = Model( 
        [UnaryResource([Jobs[m] for m in machine]) for machine in jsp.machine],
        [[job[i] < job[i+1] for i in range(jsp.nMachines-1)] for job in Jobs],

        [job[-1] < C_max for job in Jobs],
        Minimise( C_max )
        )

    return C_max,Jobs,model


def solve(param):
    ###############################################
    ##############     Solving       ##############
    ###############################################
    jsp = JSP(param['data'])

    C_max,Jobs,model = get_model(jsp)
    solver = model.load(param['solver'])

    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])

    if sys.argv[-1] == 'scheduling':
        solver.setHeuristic('Scheduling', 'Promise', 2)
        solver.solveAndRestart(GEOMETRIC, 256, 1.3)
    else:
        solver.solve()

    schedule = [[-1]*C_max.get_value() for job in jsp.job]
    index = 0
    for machine in jsp.machine:
        index += 1
        for m in machine:
            for i in range(Jobs[m].duration):
                start = Jobs[m].get_value()
                schedule[m[0]][start+i] = index

    ###############################################
    ############# Output (Matplotlib) #############
    ###############################################
    if param['print'] == 'yes':
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

    out = ''
    if solver.is_sat():
        out = str(schedule)
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out    


###############################################
##############      Input        ##############
###############################################
solvers = ['Mistral', 'MiniSat']
default = {'solver':'Mistral', 'data':'data/tiny_jsp.txt', 'print':'no', 'verbose':1, 'tcutoff':3}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)

