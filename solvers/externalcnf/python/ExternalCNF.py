from ExternalSolver import ExternalSolver
from SatWrapper import SatWrapperSolver
import Numberjack
import datetime


class ExternalCNFSolver(ExternalSolver, SatWrapperSolver):

    def __init__(self):
        # We need to call init on both parent classes, super won't do this
        # if they both inherit from base class object.
        ExternalSolver.__init__(self)
        SatWrapperSolver.__init__(self)

    def set_model(self, model, solver_id, solver_name, solver):
        self.model = model
        print "Outputting to:", self.filename
        SatWrapperSolver.output_cnf(self, self.filename)

    def parse_output(self, output):
        """
            Parses the solver output, which should conform to the output of the
            SAT Solver competitions.
            http://www.satcompetition.org/2009/format-solvers2009.html
        """
        from SatWrapper import SatWrapperIntArray
        print "c Parse output"
        values = SatWrapperIntArray()
        for line in output.split("\n"):
            line = line.strip()
            first_two = line[:2]
            if len(line) == 0:
                continue

            if first_two == "s ":
                if "UNSATISFIABLE" in line or "UNSAT" in line:
                    self.sat = Numberjack.UNSAT
                elif "SATISFIABLE" in line or "SAT" in line:
                    self.sat = Numberjack.SAT
            elif first_two == "v ":
                t = datetime.datetime.now()
                for v in map(int, line[2:].split()):
                    if v != 0:
                        values.add(v)
                print "time to add this value line: %.4f" % (datetime.datetime.now() - t).total_seconds()

        #     elif first_two == "d " or first_two == "c ":
        #         self.parse_solver_info_line(line[2:])

        t = datetime.datetime.now()
        self.store_solution(values)
        print "Time to store solution: %.4f" % (datetime.datetime.now() - t).total_seconds()
