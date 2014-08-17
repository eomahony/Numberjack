from Numberjack.ExternalSolver import ExternalCNFSolver
from Numberjack import NBJ_STD_Solver
import re


class ClaspCNFSolver(ExternalCNFSolver):

    def __init__(self):
        super(ClaspCNFSolver, self).__init__()
        self.solverexec = "clasp"

        # self.info_regexps = {  # See doc on ExternalSolver.info_regexps
        #     'nodes': (re.compile(r'^decisions[ ]+:[ ]+(?P<nodes>\d+)[ ]'), int),
        #     'failures': (re.compile(r'^conflicts[ ]+:[ ]+(?P<failures>\d+)[ ]'), int),
        # }

    def build_solver_cmd(self):
        # The Verbosity that we pass down to the solver should be at least 1 so
        # that we can parse information like number of nodes, failures, etc.
        return "%(solverexec)s %(filename)s" % vars(self)


class Solver(NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        NBJ_STD_Solver.__init__(self, "ClaspCNF", "SatWrapper", model, None, FD, clause_limit, encoding)
        self.solver_id = model.getSolverId()
        self.solver.set_model(model, self.solver_id, self.Library, solver=self)
