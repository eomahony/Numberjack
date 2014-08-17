from Numberjack.ExternalSolver import ExternalCNFSolver
from Numberjack import NBJ_STD_Solver
import re


class Riss3gSolver(ExternalCNFSolver):

    def __init__(self):
        super(Riss3gSolver, self).__init__()
        self.solverexec = "riss3g.sh"

        self.info_regexps = {  # See doc on ExternalSolver.info_regexps
            'nodes': (re.compile(r'^decisions[ ]+:[ ]+(?P<nodes>\d+)[ ]'), int),
            'failures': (re.compile(r'^conflicts[ ]+:[ ]+(?P<failures>\d+)[ ]'), int),
            'propags': (re.compile(r'^propagations[ ]+:[ ]+(?P<propags>\d+)[ ]'), int),
            'time': (re.compile(r'^CPU time[ ]+:[ ]+(?P<time>\d+\.\d+)[ ]'), float),
        }

    def build_solver_cmd(self):
        return "%(solverexec)s %(filename)s" % vars(self)


class Solver(NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        NBJ_STD_Solver.__init__(self, "Riss3g", "SatWrapper", model, None, FD, clause_limit, encoding)
        self.solver_id = model.getSolverId()
        self.solver.set_model(model, self.solver_id, self.Library, solver=self)
