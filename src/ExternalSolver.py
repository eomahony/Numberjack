import Numberjack
import subprocess as sp
import threading
import signal
import os


class Command(object):
    """
        Class for running a subprocess in isolation, optionally with time and
        memory limits.
    """

    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None
        self.exitcode = None
        self.timed_out = False
        self.timing = {}

    def parse_timing(self):
        import re
        time_re = re.compile(r"^(?P<type>(real|user|sys))\t(?P<minutes>[\d]+)m(?P<seconds>[\d]+.[\d]+)s$")
        for line in self.stderr.split("\n"):
            m = time_re.search(line)
            if m:
                d = m.groupdict()
                seconds = float(d['minutes']) * 60.0 + float(d['seconds'])
                self.timing[d['type']] = seconds

    def run(self, timelimit=0, mem_limit=None):
        def target():
            cmd = "time %s" % self.cmd
            if mem_limit:
                cmd = "ulimit -v %d; %s" % (mem_limit, cmd)  # Assumes bash
            self.process = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, preexec_fn=os.setpgrp)

            self.stdout, self.stderr = self.process.communicate()
            self.exitcode = self.process.returncode
            self.parse_timing()

        thread = threading.Thread(target=target)
        thread.start()
        if timelimit > 0:
            thread.join(float(timelimit))
        else:
            thread.join()

        if thread.is_alive():
            self.timed_out = True
            # Send the TERM signal to the process group
            os.killpg(self.process.pid, signal.SIGTERM)
            thread.join(3.0)
            if thread.is_alive():
                # Thread is still alive after TERM signal, send KILL signal.
                os.killpg(self.process.pid, signal.SIGKILL)
                self.process.kill()
                thread.join()


def print_commented(blob):
    for line in blob.split("\n"):
        if not line.startswith("c "):
            print "c",
        print line


class ExternalSolver(object):

    def __init__(self):
        self.model = None
        self.variables = []
        self.filename = self.generate_filename()
        self.timelimit = 0
        self.mem_limit = None
        self.out_object = None
        self.sat = Numberjack.UNKNOWN
        self.opt = False
        self.backtracks = -1
        self.nodes = 0
        self.failures = -1
        self.checks = -1
        self.propags = -1
        self.time = -1
        self.verbosity = 0  # Currently not passed down to external solvers
        self.var_heuristic = None  # Currently not passed down to external solvers
        self.val_heuristic = None  # Currently not passed down to external solvers
        self.randomization = None  # Currently not passed down to external solvers

    def generate_filename(self):
        import tempfile
        tf = tempfile.NamedTemporaryFile(delete=False)
        return tf.name

    def set_model(self, model, solver_id, solver_name, solver):
        pass
        # self.model = model
        # self.output_model()

        # for nj_var, i in sorted(self.out_object.njvar_mapping.iteritems(), key=lambda (k, v): v):
        #     my_var = ExternalSolverIntVariable(nj_var)
        #     nj_var.setVar(solver_id, solver_name, my_var, new_solver=solver)
        #     nj_var.solver = solver
        #     self.variables.append(my_var)

    def build_solver_cmd(self):
        pass

    def solve(self, *args, **kwargs):
        cmd = self.build_solver_cmd()
        print "Will run cmd:", cmd
        c = Command(cmd)
        c.run(timelimit=self.timelimit)
        print "Job finished"
        print_commented(c.stdout)
        print_commented(c.stderr)
        self.parse_output(c.stdout)
        return self.is_sat()

    def solveAndRestart(self, *args, **kwargs):
        return self.solve(*args, **kwargs)

    def output_model(self):
        from XCSPOut import XCSPOutput
        self.out_object = XCSPOutput(self.model)
        self.out_object.output(self.filename)

    def parse_output(self, output):
        pass
        # values = []
        # for line in output.split("\n"):
        #     line = line.strip()
        #     first_two = line[:2]
        #     if len(line) == 0:
        #         continue

        #     if first_two == "s ":
        #         print line
        #         if "UNSATISFIABLE" in line or "UNSAT" in line:
        #             self.sat = Numberjack.UNSAT
        #         elif "SATISFIABLE" in line or "SAT" in line:
        #             self.sat = Numberjack.SAT
        #     elif first_two == "v ":
        #         print "Valued:", line
        #         values.extend(map(int, line[2:].split()))
        #     elif first_two == "d " or first_two == "c ":
        #         self.parse_solver_info_line(line[2:])

        # if self.sat == Numberjack.SAT:
        #     for i, variable in enumerate(self.variables):
        #         print "Setting", variable, "equal to:", values[i]
        #         variable.value = values[i]

    # def parse_solver_info_line(self, line):
    #     bits = [a.strip() for a in line.split() if len(a.strip()) > 0]  # Strip additional whitespace
    #     for i in xrange(0, len(bits), 2):
    #         try:
    #             name, value = bits[i], bits[i + 1]
    #             print name, value
    #             if name == "NODES" or name == "NDS":
    #                 self.nodes = int(value)
    #             elif name == "BACKTRACKS":
    #                 self.backtracks = int(value)
    #             elif name == "CHECKS":
    #                 self.checks = int(value)
    #             elif name == "FAILURES":
    #                 self.failures = int(value)
    #         except Exception as e:
    #             print >> sys.stderr, str(e)
    #             pass

    def is_sat(self):
        return self.sat == Numberjack.SAT

    ## Returns True iff the solver proved unsatisfiability
    def is_unsat(self):
        return self.sat == Numberjack.UNSAT

    def is_opt(self):
        return self.opt

    def setVerbosity(self, verbosity):
        self.verbosity = verbosity

    def setHeuristic(self, var_h, val_h, randomization):
        self.var_heuristic = var_h
        self.val_heuristic = val_h
        self.randomization = randomization

    def setTimeLimit(self, timelimit):
        self.timelimit = timelimit

    def getBacktracks(self):
        return self.backtracks

    def getNodes(self):
        return self.nodes

    def getFailures(self):
        return self.failures

    def getChecks(self):
        return self.checks

    def getPropags(self):
        return self.propags

    def getTime(self):
        return self.time

    def printPython(self):
        return repr(self)


# class Solver(Numberjack.NBJ_STD_Solver):
#     def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
#         # if X:
#         #     import sys
#         #     print >> sys.stderr, "Warning explicitly specifying the decision variables for an external CSP solver is currently not supported."

#         # We pass an empty model to NBJ_STD_Solver to prevent it from loading
#         # and trying to decompse each expression.
#         Numberjack.NBJ_STD_Solver.__init__(self, "ExternalSolver", "ExternalSolver", None, None, FD, clause_limit, encoding)
#         self.solver_id = model.getSolverId()
#         self.model = model
#         self.solver.set_model(model, self.solver_id, self.Library, solver=self)
