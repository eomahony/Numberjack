import subprocess as sp
import signal
import threading
import os


SIGTERM_TIMEOUT = 1.0


class Command(object):
    def __init__(self, cmd, memlimit=None):
        self.cmd = cmd
        self.memlimit = memlimit
        self.process = None
        self.stdout = None
        self.stderr = None
        self.exitcode = None
        self.timed_out = False

    def run(self, timeout=None):
        def target():
            cmd = self.cmd
            if self.memlimit:
                cmd = "ulimit -v %d; %s" % (self.memlimit, cmd)
            self.process = sp.Popen(cmd,
                                    stdout=sp.PIPE, stderr=sp.PIPE,
                                    shell=True, preexec_fn=os.setpgrp)
            self.stdout, self.stderr = self.process.communicate()
            self.exitcode = self.process.returncode

        thread = threading.Thread(target=target)
        thread.start()
        thread.join(float(timeout))
        if thread.is_alive():
            self.timed_out = True

            # Send the TERM signal to all the process groups
            os.killpg(self.process.pid, signal.SIGTERM)
            thread.join(SIGTERM_TIMEOUT)
            if thread.is_alive():
                # Send the KILL signal if the process hasn't exited by now.
                os.killpg(self.process.pid, signal.SIGKILL)
                self.process.kill()
                thread.join()
