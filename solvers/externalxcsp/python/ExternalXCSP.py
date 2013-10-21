from ExternalSolver import ExternalSolver
import Numberjack
import sys


class ExternalXCSPIntVariable(object):

    def __init__(self, nj_var):
        self.nj_var = nj_var
        self.value = None

    def get_value(self):
        print "ExternalXCSPIntVariable.get_value", str(self.value)
        return self.value


class ExternalXCSPSolver(ExternalSolver):

    def __init__(self):
        super(ExternalXCSPSolver, self).__init__()

    def set_model(self, model, solver_id, solver_name, solver):
        self.model = model
        self.output_model()

        for nj_var, i in sorted(self.out_object.njvar_mapping.iteritems(), key=lambda (k, v): v):
            my_var = ExternalXCSPIntVariable(nj_var)
            nj_var.setVar(solver_id, solver_name, my_var, new_solver=solver)
            nj_var.solver = solver
            self.variables.append(my_var)

    def output_model(self):
        from XCSPOut import XCSPOutput
        self.out_object = XCSPOutput(self.model)
        self.out_object.output(self.filename)

    def parse_output(self, output):
        """
            Parses the solver output, which should conform to the output of the
            CSP Solver competitions. http://cpai.ucc.ie/09/
        """
        values = []
        for line in output.split("\n"):
            line = line.strip()
            first_two = line[:2]
            if len(line) == 0:
                continue

            if first_two == "s ":
                print line
                if "UNSATISFIABLE" in line or "UNSAT" in line:
                    self.sat = Numberjack.UNSAT
                elif "SATISFIABLE" in line or "SAT" in line:
                    self.sat = Numberjack.SAT
            elif first_two == "v ":
                print "Valued:", line
                values.extend(map(int, line[2:].split()))
            elif first_two == "d " or first_two == "c ":
                self.parse_solver_info_line(line[2:])

        if self.sat == Numberjack.SAT:
            for i, variable in enumerate(self.variables):
                print "Setting", variable, "equal to:", values[i]
                variable.value = values[i]

    def parse_solver_info_line(self, line):
        bits = [a.strip() for a in line.split() if len(a.strip()) > 0]  # Strip additional whitespace
        for i in xrange(0, len(bits), 2):
            try:
                name, value = bits[i], bits[i + 1]
                print name, value
                if name == "NODES" or name == "NDS":
                    self.nodes = int(value)
                elif name == "BACKTRACKS":
                    self.backtracks = int(value)
                elif name == "CHECKS":
                    self.checks = int(value)
                elif name == "FAILURES":
                    self.failures = int(value)
            except Exception as e:
                print >> sys.stderr, str(e)
                pass


class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        # if X:
        #     import sys
        #     print >> sys.stderr, "Warning explicitly specifying the decision variables for an external CSP solver is currently not supported."

        # We pass an empty model to NBJ_STD_Solver to prevent it from loading
        # and trying to decompse each expression.
        Numberjack.NBJ_STD_Solver.__init__(self, "ExternalXCSP", "ExternalXCSP", None, None, FD, clause_limit, encoding)
        self.solver_id = model.getSolverId()
        self.model = model
        self.solver.set_model(model, self.solver_id, self.Library, solver=self)
