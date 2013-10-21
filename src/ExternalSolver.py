import Numberjack
import subprocess as sp
import threading
import signal
import atexit
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
        self.verbosity = 0  # Currently not passed down to external solvers
        self.var_heuristic = None  # Currently not passed down to external solvers
        self.val_heuristic = None  # Currently not passed down to external solvers
        self.randomization = None  # Currently not passed down to external solvers

        # info_regexps should be a dictionary containing the attributes to parse
        # from the solver's output. The dictionary key is the attribute name,
        # the value for that key should be a tuple of two items: a regular
        # expression and a function which takes a string as parameter and
        # returns the data you would like to store. The regular expression
        # should have a symbolic group name matching the key of the dictionary
        # entry and will be stored as a property on self, i.e. it is equivalent
        # to self.key = cast_func(group_match_str)
        self.info_regexps = {}

        # Register the clean_up function to be run when python is exiting, there may be a better time for this...
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

