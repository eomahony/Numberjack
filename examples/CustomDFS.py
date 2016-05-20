from Numberjack import *


def lexvar(variables):
    "Returns the first unassigned variable in 'variables'."
    for x in variables:
        if x.get_size() > 1:
            return x
    return None


def binarydfs(solver, variables, variableselection=lexvar):
    """
        Performs binary branching depth first search.
        'variables' should be a list or a VarArray of the decision variables.
    """

    proven_infeasibility = False
    prevvar, prevval = None, None
    depth = 0

    def log(s):
        print("  " * depth, s)

    while not proven_infeasibility:
        if solver.propagate():  # Propagate to a fixed point. Returns False if there was a conflict
            # left branch
            log("Variable domains at depth %d:" % depth)
            for x in variables:
                log(x)

            x = variableselection(variables)
            if x is None:  # Complete, no more variables to assign
                return solver.get_solution()
            v = x.get_min()

            prevvar, prevval = x.name(), v  # Save for debugging output
            log("Decision: %s == %d" % (prevvar, prevval))

            solver.save()
            solver.post(x == v)
            depth += 1

        else:  # right branch
            proven_infeasibility = not solver.undo()
            if not proven_infeasibility:
                solver.deduce()  # Invert the previous decision, i.e. x != v
                log("Deduce: %s != %d" % (prevvar, prevval))
                depth -= 1

    return None


def main():
    v1, v2 = VarArray(2, 10)
    m = Model(v1 != v2)
    # m = Model(v1 + v2 > 4, v1 > v2)  # Force backtracking with lex var ordering

    # s = m.load("Mistral")
    s = m.load("MiniSat")
    s.startNewSearch()
    binarydfs(s, [v1, v2])
    if s.is_sat():
        print("\nSolution:")
        print("v1 = %d, v2 = %d" % (v1.get_value(), v2.get_value()))


if __name__ == "__main__":
    main()
