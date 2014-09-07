# from Numberjack.ExternalSolver import which
import os


def available_solvers():
    excludedsolvers = ['__init__', 'MipWrapper', 'SatWrapper']
    solvers = set()
    thisdir = os.path.dirname(os.path.abspath(__file__))
    for basename in os.listdir(thisdir):
        if not basename.endswith(".py") and not basename.endswith(".pyc"):
            continue
        modulename = basename.replace(".pyc", "").replace(".py", "")

        if modulename in solvers:
            continue  # Skip duplicates

        # If the module has an is_available attritube, call it, but if not
        # assume it is available. Typically the external solvers have this func
        try:
            njsolvers = "Numberjack.solvers"
            solvermod = "%s.%s" % (njsolvers, modulename)
            s = __import__(solvermod, fromlist=[njsolvers])
            solver = getattr(s, "%sSolver" % modulename, None)
            if solver is not None:
                testsolver = solver()
                if hasattr(testsolver, "is_available"):
                    if not testsolver.is_available():
                        continue
                solvers.add(modulename)
        except Exception as e:
            pass

    return list(sorted(solvers))
