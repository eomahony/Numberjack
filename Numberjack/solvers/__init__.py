import os


def available_solvers():
    excludedsolvers = ['__init__', 'MipWrapper', 'SatWrapper']
    solvers = set()
    thisdir = os.path.dirname(os.path.abspath(__file__))
    for basename in os.listdir(thisdir):
        if not basename.endswith(".py") and not basename.endswith(".pyc"):
            continue
        modulename = basename.replace(".pyc", "").replace(".py", "")
        if modulename not in excludedsolvers:
            solvers.add(modulename)

    return list(sorted(solvers))
