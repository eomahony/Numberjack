from Numberjack.ExternalSolver import ExternalCNFSolver
from Numberjack import NBJ_STD_Solver


class BreakIDGlucoseSolver(ExternalCNFSolver):

    def __init__(self):
        super(BreakIDGlucoseSolver, self).__init__()
        self.solverexec = "runBreakIDGlucose.sh"

    def build_solver_cmd(self):
        return "%(solverexec)s %(filename)s %(tempdir)s" % vars(self)


class Solver(NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        NBJ_STD_Solver.__init__(self, "BreakIDGlucose", "SatWrapper", model, None, FD, clause_limit, encoding)
        self.solver_id = model.getSolverId()
        self.solver.set_model(model, self.solver_id, self.Library, solver=self)
