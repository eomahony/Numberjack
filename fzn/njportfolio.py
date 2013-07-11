import subprocess as sp
import datetime
import signal
import threading
import os
import sys


result_poll_timeout = 0.2


def run_cmd(process_name, pid_queue, result_queue, cmd, memlimit):
    if memlimit:
        cmd = "ulimit -v %d; %s" % (memlimit, cmd)
    process = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE,
                       shell=True, preexec_fn=os.setpgrp)
    # Tell the parent the process pid, so it can be killed later
    pid_queue.put(process.pid, True)
    stdout, stderr = process.communicate()
    exitcode = process.returncode
    if exitcode == 0:
        result_queue.put([True, exitcode, process_name, stdout, stderr], True)
    else:
        result_queue.put([False, exitcode, process_name, stdout, stderr], True)


def njportfolio(njfilename, cores, timeout, memlimit):
    from multiprocessing import Queue, cpu_count
    from Queue import Empty

    start_time = datetime.datetime.now()

    if cores <= 0 or cores > cpu_count():
        cores = cpu_count()

    solvers = ['Mistral', 'Gurobi', 'Toulbar2', 'SCIP', 'MiniSat']
    result_queue = Queue()
    pid_queue = Queue()
    threads = []
    for i, s in enumerate(solvers):
        if i >= cores:
            break
        cmd = "python %s -solver %s -tcutoff %d -threads 1" % (njfilename, s, int(timeout))
        process_name = s
        thread = threading.Thread(target=run_cmd, args=(process_name, pid_queue, result_queue, cmd, int(memlimit / cores)))
        threads.append(thread)
        thread.start()
        print "% Launching:", cmd

    num_finished = 0
    finished_names = []
    while True:
        try:
            success, exitcode, process_name, stdout, stderr = result_queue.get(True, result_poll_timeout)
            num_finished += 1
            finished_names.append(process_name)
            if success:
                t = (datetime.datetime.now() - start_time).total_seconds()
                print "%% Solver %s finished after %.1f" % (process_name, t)
                print stdout
                break
            else:
                print "%% Failed: %s exitcode: %d" % (process_name, exitcode)
                print_commented(stdout)
                print_commented(stderr)

            if num_finished == len(threads):
                break
        except Empty:
            # Nothing posted to the result_queue yet.
            pass

    # Tidy up, kill all processes that did not finish.
    num_pids_seen = 0
    while not pid_queue.empty() and num_pids_seen < len(threads):
        try:
            pid = pid_queue.get()
            num_pids_seen += 1
            os.killpg(pid, signal.SIGKILL)
        except Empty:
            pass
        except OSError:
            pass  # Process already finished.
    print "%% Total time in njportfolio: %.1f" % (datetime.datetime.now() - start_time).total_seconds()


def print_commented(text):
    for line in text.split("\n"):
        print "%", line


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print >> sys.stderr, "Usage: python %s njfilename cores timeout memlimit" % sys.argv[0]
    njportfolio(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))
