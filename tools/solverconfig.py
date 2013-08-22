#!/usr/bin/env python

import os
import sys


"""This script will attempt to find the home directory for a named solver."""


def get_solver_home(solvername):
    from distutils.spawn import find_executable

    if solvername == "CPLEX":
        # Try for environmental variable first
        env_path = os.getenv('CPLEXDIR')
        if env_path and len(env_path.strip()) > 0:
            return env_path

        # Try to find the cplex binary in the PATH
        ex_path = find_executable('cplex')
        if ex_path:
            ex_path = os.path.realpath(ex_path)  # Expand symbolic links if any
            ex_dir = os.path.dirname(ex_path)  # Path to the bin directory
            return os.path.abspath(os.path.join(ex_dir, os.pardir, os.pardir))
    else:
        print >> sys.stderr, "Error unknown solver name", solvername

    return ""


def get_cflags(solvername):
    solverhome = get_solver_home(solvername)
    if solvername == "CPLEX":
        concertdir = get_concert_dir(solverhome)
        cplexincdir = os.path.join(solverhome, "include")
        concertincdir = os.path.join(concertdir, "include")

        cplexstd = "-O -fPIC -fexceptions -DNDEBUG -DIL_STD -I%s -I%s" % (cplexincdir, concertincdir)
        if sys.platform == "darwin":
            return cplexstd
        elif sys.platform.startswith('linux'):
            return "-fno-strict-aliasing %s" % (cplexstd)
        else:
            print >> sys.stderr, "Error unknown platform", sys.platform
            sys.exit(1)
    else:
        print >> sys.stderr, "Error unknown solvername", solvername
        sys.exit(1)


def get_lflags(solvername):
    solverhome = get_solver_home(solvername)
    if solvername == "CPLEX":
        concertdir = get_concert_dir(solverhome)
        cplexlibdir = get_cplex_lib_dir(solverhome)
        concertlibdir = get_cplex_lib_dir(concertdir)

        cplexstd = "-L%s -L%s -lconcert -lilocplex -lcplex -lm -lpthread" % (cplexlibdir, concertlibdir)
        if sys.platform == "darwin":
            return "%s -framework CoreFoundation -framework IOKit" % cplexstd
        elif sys.platform.startswith('linux'):
            return cplexstd
        else:
            print >> sys.stderr, "Error unknown platform", sys.platform
            sys.exit(1)
    else:
        print >> sys.stderr, "Error unknown solvername", solvername
        sys.exit(1)


def get_concert_dir(cplexdir):
    return os.path.abspath(os.path.join(cplexdir, os.pardir, "concert"))


def get_cplex_lib_dir(cplexdir):
    """
    CPLEX lib structure: cplexdir/lib/ARCHITECTURE/BUILDTYPE/
    This function will walk the cplexdir/lib folder to find the
    first folder with a library in it.
    """
    libfolder = os.path.join(cplexdir, "lib")
    for dirpath, dirnames, filenames in os.walk(libfolder):
        for f in filenames:
            if f.endswith(".a"):
                return dirpath
    print >> sys.stderr, "Error could not find the lib folder for %s" % cplexdir
    sys.exit(1)


def exit_usage():
    print >> sys.stderr, "Usage: %s CPLEX (home|cflags|lflags)" % sys.argv[0]
    sys.exit(1)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        exit_usage()
    else:
        solvername = sys.argv[1]
        query = sys.argv[2]

        if query == "home":
            print get_solver_home(solvername)
        elif query == "cflags":
            print get_cflags(solvername)
        elif query == "lflags":
            print get_lflags(solvername)
        else:
            print >> sys.stderr, "Error unknown query argument:", query
            exit_usage()
