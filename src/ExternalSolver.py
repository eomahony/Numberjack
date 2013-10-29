from SatWrapper import SatWrapperSolver
import Numberjack
import subprocess as sp
import threading
import datetime
import signal
import atexit
import sys
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
        self.stdout = ""
        self.stderr = ""

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

        def tidy_up(*args, **kwargs):
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

        # Set handlers for term and interupt signals
        signal.signal(signal.SIGTERM, tidy_up)
        signal.signal(signal.SIGINT, tidy_up)

        thread.start()
        if timelimit > 0:
            thread.join(float(timelimit))
        else:
            thread.join()
        tidy_up()


def print_commented(blob, comment_prefix="c"):
    for line in blob.split("\n"):
        if not line.startswith("%s " % comment_prefix):
            print comment_prefix,
        print line


## Base class for using an external solver binary that doesn't have a native interface.
#
#  Provides common functions that are used during 
#  See also ExternalCNF and ExternalXCSP which subclass this.
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
        self.threads = 1

        # info_regexps should be a dictionary containing the attributes to parse
        # from the solver's output. The dictionary key is the attribute name,
        # the value for that key should be a tuple of two items: a regular
        # expression and a function which takes a string as parameter and
        # returns the data you would like to store. The regular expression
        # should have a symbolic group name matching the key of the dictionary
        # entry and will be stored as a property on self, i.e. it is equivalent
        # to self.key = cast_func(group_match_str)
        self.info_regexps = {}

        # Currently not passed down to external solvers since we require a
        # specific level of output to parse their output. Instead, we control
        # the level this output we print based on this value.
        self.verbosity = 0

        # None of these are currently being passed down to external solvers:
        self.var_heuristic = None
        self.val_heuristic = None
        self.randomization = None

        # Register the clean_up function to be run when python is exiting, there
        # may be a better time for this...
        atexit.register(ExternalSolver.clean_up, self)

    def generate_filename(self):
        import tempfile
        tf = tempfile.NamedTemporaryFile(delete=False)
        return tf.name

    def clean_up(self):
        if self.filename:
            try:
                os.remove(self.filename)
            except Exception:
                pass  # shhh

    def set_model(self, model, solver_id, solver_name, solver):
        pass

    def build_solver_cmd(self):
        pass

    def solve(self, *args, **kwargs):
        cmd = self.build_solver_cmd()
        print "c Running:", cmd
        c = Command(cmd)
        c.run(timelimit=self.timelimit)
        print "c External solver finished."
        print_commented(c.stdout)
        print_commented(c.stderr)
        self.parse_output(c.stdout)
        return self.is_sat()

    def solveAndRestart(self, *args, **kwargs):
        return self.solve(*args, **kwargs)

    def parse_output(self, output):
        pass

    def parse_solver_info_line(self, line):
        for key, (regexp, cast_func) in self.info_regexps.iteritems():
            match = regexp.match(line)
            if match:
                setattr(self, key, cast_func(match.groupdict()[key]))

    def is_sat(self):
        return self.sat == Numberjack.SAT

    # Returns True iff the solver proved unsatisfiability, not if no solution was found.
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

    def setThreadCount(self, num_threads):
        if num_threads >= 1:
            self.threads = num_threads

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


# ---------- External CNF ----------


class ExternalCNFSolver(ExternalSolver, SatWrapperSolver):

    def __init__(self):
        # We need to call init on both parent classes, super won't do this
        # if they both inherit from base class object.
        ExternalSolver.__init__(self)
        SatWrapperSolver.__init__(self)

    def set_model(self, model, solver_id, solver_name, solver):
        self.model = model
        # print "c minimise_obj:", str(self.minimise_obj), self.minimise_obj.get_min(), self.minimise_obj.get_max()
        # print "c maximise_obj:", str(self.maximise_obj)
        print "c Outputting to:", self.filename
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
                # t = datetime.datetime.now()
                for v in map(int, line[2:].split()):
                    if v != 0:
                        values.add(v)
                # print "time to add this value line: %.4f" % (datetime.datetime.now() - t).total_seconds()

        #     elif first_two == "d " or first_two == "c ":
        #         self.parse_solver_info_line(line[2:])

        t = datetime.datetime.now()
        self.store_solution(values)
        print "c Time to store solution: %.4f" % (datetime.datetime.now() - t).total_seconds()


# ---------- External XCSP ----------


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


# class Solver(Numberjack.NBJ_STD_Solver):
#     def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
#         # if X:
#         #     import sys
#         #     print >> sys.stderr, "Warning explicitly specifying the decision variables for an external CSP solver is currently not supported."

#         # We pass an empty model to NBJ_STD_Solver to prevent it from loading
#         # and trying to decompse each expression.
#         Numberjack.NBJ_STD_Solver.__init__(self, "ExternalXCSP", "ExternalXCSP", None, None, FD, clause_limit, encoding)
#         self.solver_id = model.getSolverId()
#         self.model = model
#         self.solver.set_model(model, self.solver_id, self.Library, solver=self)
