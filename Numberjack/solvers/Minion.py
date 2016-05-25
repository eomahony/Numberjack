# This file represents the interface from Numberjack to the Minion constraint
# solver http://minion.sourceforge.net
# If you have minion installed on your machine and the 'minion' binary is in
# your path, then you can use Minion with Numberjack.
#
# Some Numberjack constraints still to be implemented: Mod, Abs, Table There are
# more constraints supported by Minion that are not in Numberjack yet, It would
# be goot to support these in Numberjack but decompositions should also be added
# so that they can be used with other solvers. In particular reifyimply would be
# nice.

from __future__ import print_function, division

from Numberjack.ExternalSolver import ExternalSolver
from Numberjack import NBJ_STD_Solver, SAT, UNSAT
import re


class Minion_Expression(object):

    def __init__(self, *args):
        self.nbj_ident = None
        self.varname = None
        self.solver = None
        self.lb = self.ub = None

    def add(self, solver, toplevel):
        pass

    def has_been_added(self):
        return self.solver is not None

    def get_min(self):
        return self.lb

    def get_max(self):
        return self.ub

    def get_size(self):
        return self.ub - self.lb + 1


class Minion_IntVar(Minion_Expression):

    def __init__(self, *args):
        super(Minion_IntVar, self).__init__()
        self.domain = self.value = None

        if len(args) == 0:  # Boolean
            self.lb, self.ub = 0, 1
        elif len(args) == 3:
            self.lb, self.ub, self.nbj_ident = args
        elif len(args) == 2:
            if hasattr(args[0], "__iter__"):
                self.domain, self.nbj_ident = args
                self.lb, self.ub = min(self.domain), max(self.domain)
            else:
                self.lb, self.ub = args
                self.nbj_ident = -1
        else:
            raise Exception("Unknown constructor for %s" % str(type(self)))

    def add(self, solver, toplevel):
        if self.lb == self.ub:
            if self.lb is None or self.ub is None:
                raise Exception(
                    "Error Minion_IntVar.add was called on a variable without "
                    " a lower bound or upper bound. nbj_ident:%d" %
                    self.nbj_ident)
            self.value = self.lb
            return str(self.lb)

        if not self.has_been_added():
            self.varname = "x%d" % (solver.variablecount)  # Minion variable name
            solver.variablecount += 1
            self.solver = solver

            if self.lb == 0 and self.ub == 1:
                solver.create_variable(self, self.varname, "BOOL %s" % self.varname)
            elif self.domain and len(self.domain) != (self.ub - self.lb) + 1:
                dom_str = csvstr(list(map(str, self.domain)))
                solver.create_variable(self, self.varname, "SPARSEBOUND %s {%s}" % (self.varname, dom_str))
            else:
                solver.create_variable(self, self.varname, "DISCRETE %s {%d..%d}" % (self.varname, self.lb, self.ub))
        return self

    def get_value(self):
        return self.value


class MinionIntArray(list):
    add = list.append
    size = list.__len__


class MinionExpArray(list):
    add = list.append
    size = list.__len__


class Minion_binop(Minion_Expression):

    def __init__(self, arg1, arg2):
        super(Minion_binop, self).__init__()
        self.vars = [arg1, arg2]
        self.lb, self.ub = 0, 1

    def add(self, solver, toplevel):
        if not self.has_been_added():
            self.solver = solver
            self.vars[0] = self.vars[0].add(solver, False)
            if isinstance(self.vars[1], Minion_Expression):
                self.vars[1] = self.vars[1].add(solver, False)
            if not toplevel:
                self.auxvar = Minion_IntVar().add(solver, False)
                self.varname = varname(self.auxvar)
        return self


class Minion_ne(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_ne, self).add(solver, toplevel)
            consstr = "diseq(%s, %s)" % (varname(self.vars[0]), varname(self.vars[1]))
            if toplevel:
                solver.print_constraint(consstr)
            else:
                solver.print_constraint("reify(%s, %s)" % (consstr, varname(self)))
        return self


class Minion_eq(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_eq, self).add(solver, toplevel)
            consstr = "eq(%s, %s)" % (varname(self.vars[0]), varname(self.vars[1]))
            if toplevel:
                solver.print_constraint(consstr)
            else:
                solver.print_constraint("reify(%s, %s)" % (consstr, varname(self)))
        return self


class Minion_lt(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_lt, self).add(solver, toplevel)
            consstr = "ineq(%s, %s, -1)" % (varname(self.vars[0]), varname(self.vars[1]))
            if toplevel:
                solver.print_constraint(consstr)
            else:
                solver.print_constraint("reify(%s, %s)" % (consstr, varname(self)))
        return self


class Minion_le(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_le, self).add(solver, toplevel)
            consstr = "ineq(%s, %s, 0)" % (varname(self.vars[0]), varname(self.vars[1]))
            if toplevel:
                solver.print_constraint(consstr)
            else:
                solver.print_constraint("reify(%s, %s)" % (consstr, varname(self)))
        return self


class Minion_gt(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_gt, self).add(solver, toplevel)
            consstr = "ineq(%s, %s, -1)" % (varname(self.vars[1]), varname(self.vars[0]))
            if toplevel:
                solver.print_constraint(consstr)
            else:
                solver.print_constraint("reify(%s, %s)" % (consstr, varname(self)))
        return self


class Minion_ge(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_ge, self).add(solver, toplevel)
            consstr = "ineq(%s, %s, 0)" % (varname(self.vars[1]), varname(self.vars[0]))
            if toplevel:
                solver.print_constraint(consstr)
            else:
                solver.print_constraint("reify(%s, %s)" % (consstr, varname(self)))
        return self


class Minion_or(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_or, self).add(solver, toplevel)

            # watched-or requires the constraints to be specified in brackets,
            # possible to change to this later? for now decompose to
            # Sum(reified variables) >= 1
            s = Minion_Sum(self.vars)
            return Minion_ge(s, 1).add(solver, toplevel)

        return self


class Minion_and(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            if toplevel:
                self.vars[0].add(solver, toplevel)
                self.vars[1].add(solver, toplevel)
            else:
                super(Minion_and, self).add(solver, toplevel)

                # watched-and requires the constraints to be specified in
                # brackets, possible to change to this later? for now decompose
                # to Sum(reified variables) >= 1
                s = Minion_Sum(self.vars)
                return Minion_eq(s, len(self.vars)).add(solver, toplevel)

        return self


class Minion_mul(Minion_binop):

    def __init__(self, *args):
        super(Minion_mul, self).__init__(*args)
        var1, var2 = self.vars
        l1, u1 = var1.get_min(), var1.get_max()

        if isinstance(var2, Minion_Expression):
            l2, u2 = var2.get_min(), var2.get_max()
        else:
            if not isinstance(var2, int):
                raise Exception("Multiplication must be either by an expression or int, got '%s'." % str(type(var2)))
            l2 = u2 = var2

        self.lb = min(l1 * l2, l1 * u2, u1 * l2, u1 * u2)
        self.ub = max(l1 * l2, l1 * u2, u1 * l2, u1 * u2)

    def add(self, solver, toplevel):
        assert not toplevel, "Constraint is only valid as a sub-expression."
        if not self.has_been_added():
            super(Minion_mul, self).add(solver, toplevel)
            self.auxvar = Minion_IntVar(self.lb, self.ub).add(solver, False)
            self.varname = varname(self.auxvar)
            solver.print_constraint("product(%s, %s, %s)" % (varname(self.vars[0]), varname(self.vars[1]), varname(self)))
        return self


class Minion_div(Minion_binop):

    def __init__(self, *args):
        super(Minion_div, self).__init__(*args)
        var1, var2 = self.vars
        l1, u1 = var1.get_min(), var1.get_max()

        if isinstance(var2, Minion_Expression):
            l2, u2 = var2.get_min(), var2.get_max()
        else:
            if not isinstance(var2, int):
                raise Exception("Division must be either by an expression or int, got '%s'." % str(type(var2)))
            l2 = u2 = var2

        self.lb = min(l1 / l2, l1 / u2, u1 / l2, u1 / u2)
        self.ub = max(l1 / l2, l1 / u2, u1 / l2, u1 / u2)

    def add(self, solver, toplevel):
        assert not toplevel, "Constraint is only valid as a sub-expression."
        if not self.has_been_added():
            super(Minion_div, self).add(solver, toplevel)
            self.auxvar = Minion_IntVar(self.lb, self.ub).add(solver, False)
            self.varname = varname(self.auxvar)
            solver.print_constraint("div(%s, %s, %s)" % (varname(self.vars[0]), varname(self.vars[1]), varname(self)))
        return self


class Minion_AllDiff(Minion_Expression):

    def __init__(self, *args):
        super(Minion_AllDiff, self).__init__()
        if len(args) == 1:
            self.vars = args[0]
        elif len(args) == 2:
            self.vars = [args[0], args[1]]

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_AllDiff, self).add(solver, toplevel)
            for i in range(len(self.vars)):
                if isinstance(self.vars[i], Minion_Expression):
                    self.vars[i] = self.vars[i].add(solver, False)

            if len(self.vars) == 1:  # Just return the variable
                return self.vars[0]
            elif len(self.vars) == 2:
                # Replace a binary alldiff with a disequality
                ne = Minion_ne(*self.vars)
                return ne.add(solver, toplevel)
            else:
                assert toplevel, "Constraint not implemented as a sub-expression/reified yet."

                solver.print_constraint("gacalldiff([%s])" % (csvstr(list(map(varname, self.vars)))))
        return self


class Minion_LeqLex(Minion_Expression):

    def __init__(self, children):
        self.vars = children
        super(Minion_LeqLex, self).__init__()

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_LeqLex, self).add(solver, toplevel)
            for i in range(len(self.vars)):
                self.vars[i] = self.vars[i].add(solver, False)

            assert toplevel, "Constraint not implemented as a sub-expression/reified yet."
            names = [varname(x) for x in self.vars]
            solver.print_constraint(
                "lexleq([%s], [%s])" %
                (csvstr(names[:len(names)/2]), csvstr(names[len(names)/2:])))
        return self


class Minion_LessLex(Minion_Expression):

    def __init__(self, children):
        self.vars = children
        super(Minion_LessLex, self).__init__()

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_LessLex, self).add(solver, toplevel)
            for i in range(len(self.vars)):
                self.vars[i] = self.vars[i].add(solver, False)

            assert toplevel, "Constraint not implemented as a sub-expression/reified yet."
            names = [varname(x) for x in self.vars]
            solver.print_constraint(
                "lexless([%s], [%s])" %
                (csvstr(names[:len(names)/2]), csvstr(names[len(names)/2:])))
        return self


class Minion_Gcc(Minion_Expression):

    def __init__(self, children, vals, lb_card, ub_card):
        self.vars = children
        self.vals = vals
        self.lb_card = lb_card
        self.ub_card = ub_card
        super(Minion_Gcc, self).__init__()

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_Gcc, self).add(solver, toplevel)
            for i in range(len(self.vars)):
                self.vars[i] = self.vars[i].add(solver, False)

            assert toplevel, "Constraint not implemented as a sub-expression/reified yet."
            names = [varname(x) for x in self.vars]
            value_str = csvstr(self.vals)
            auxvariables = [Minion_IntVar(l, u) for l, u in zip(self.lb_card, self.ub_card)]
            for i in range(len(auxvariables)):
                auxvariables[i] = auxvariables[i].add(solver, False)
            vec_str = csvstr(list(map(varname, auxvariables)))
            solver.print_constraint("gccweak([%s], [%s], [%s])" %
                                    (csvstr(names), value_str, vec_str))
        return self


class Minion_Sum(Minion_Expression):

    def __init__(self, *args):
        super(Minion_Sum, self).__init__()
        self.offset = 0
        self.weights = None
        self.auxvar = None

        if len(args) >= 1 and len(args) <= 3:
            if hasattr(args[0], '__iter__'):
                self.vars = args[0]
            else:
                self.vars = [args[0]]

            if len(args) >= 2:
                self.weights = args[1]

            if len(args) == 3:
                self.offset = args[2]

        elif len(args) == 4:
            self.vars = [args[0], args[1]]
            self.weights = args[2]
            self.offset = args[3]

        else:
            raise Exception("Invalid constructor to Minion_Sum args: %s" % str(args))

        if self.weights:
            self.lb = sum(w * x.get_min() if w >= 0 else w * x.get_max() for w, x in zip(self.weights, self.vars)) + self.offset
            self.ub = sum(w * x.get_max() if w >= 0 else w * x.get_min() for w, x in zip(self.weights, self.vars)) + self.offset
        else:
            self.lb = sum(x.get_min() for x in self.vars) + self.offset
            self.ub = sum(x.get_max() for x in self.vars) + self.offset

    def add(self, solver, toplevel):
        assert not toplevel, "Constraint is only valid as a sub-expression."
        if not self.has_been_added():
            super(Minion_Sum, self).add(solver, toplevel)
            for i in range(len(self.vars)):
                self.vars[i] = self.vars[i].add(solver, False)

            if len(self.vars) == 1 and self.offset == 0 and \
                    (self.weights is None or self.weights[0] == 1):
                return self.vars[0]

            names = [varname(x) for x in self.vars]
            self.auxvar = Minion_IntVar(self.lb, self.ub).add(solver, False)
            self.varname = varname(self.auxvar)

            if self.offset != 0:
                # To achieve (a+b+c + offset = x) with Minion, we must use
                # auxiliary expressions: (a+b+c = y) and (x = y + offset).
                # Auxiliary sum is not needed if there is only one variable 'a'
                if len(self.vars) > 1:
                    args = [self.vars]
                    if self.weights:
                        args.append(self.weights)
                    auxsum = Minion_Sum(*args).add(solver, False)
                    self.vars = [auxsum]
                    self.weights = None
                    names = [varname(auxsum)]

                solver.print_constraint("ineq(%s, %s, %d)" % (varname(self), names[0], self.offset))
                solver.print_constraint("ineq(%s, %s, %d)" % (names[0], varname(self), -self.offset))

            else:
                varvecstr = csvstr(names)
                if self.weights and any(x != 1 for x in self.weights):  # Weighted
                    constantvecstr = csvstr(self.weights)
                    solver.print_constraint("weightedsumgeq([%s], [%s], %s)" % (constantvecstr, varvecstr, varname(self)))
                    solver.print_constraint("weightedsumleq([%s], [%s], %s)" % (constantvecstr, varvecstr, varname(self)))
                else:
                    solver.print_constraint("sumgeq([%s], %s)" % (varvecstr, varname(self)))
                    solver.print_constraint("sumleq([%s], %s)" % (varvecstr, varname(self)))
            return self.auxvar

        return self

    def get_value(self):
        return self.auxvar.get_value()


class Minion_Max(Minion_Expression):

    def __init__(self, *args):
        super(Minion_Max, self).__init__()
        if len(args) == 1:
            self.vars = args[0]
        elif len(args) == 2:
            self.vars = [args[0], args[1]]
        else:
            raise Exception("Invalid constructor to Minion_Max args: %s" % str(args))

        self.lb = max(x.get_min() for x in self.vars)
        self.ub = max(x.get_max() for x in self.vars)

    def add(self, solver, toplevel):
        assert not toplevel, "Constraint is only valid as a sub-expression."
        if not self.has_been_added():
            super(Minion_Max, self).add(solver, toplevel)
            for i in range(len(self.vars)):
                self.vars[i] = self.vars[i].add(solver, False)

            if len(self.vars) == 1:
                return self.vars[0]

            names = [varname(x) for x in self.vars]
            auxvar = Minion_IntVar(self.lb, self.ub).add(solver, False)
            solver.print_constraint("max([%s], %s)" % (csvstr(names), varname(auxvar)))
            return auxvar

        return self


class Minion_Min(Minion_Expression):

    def __init__(self, *args):
        super(Minion_Min, self).__init__()
        if len(args) == 1:
            self.vars = args[0]
        elif len(args) == 2:
            self.vars = [args[0], args[1]]
        else:
            raise Exception("Invalid constructor to Minion_Min args: %s" % str(args))

        self.lb = min(x.get_min() for x in self.vars)
        self.ub = min(x.get_max() for x in self.vars)

    def add(self, solver, toplevel):
        assert not toplevel, "Constraint is only valid as a sub-expression."
        if not self.has_been_added():
            super(Minion_Min, self).add(solver, toplevel)
            for i in range(len(self.vars)):
                self.vars[i] = self.vars[i].add(solver, False)

            if len(self.vars) == 1:
                return self.vars[0]

            names = [varname(x) for x in self.vars]
            auxvar = Minion_IntVar(self.lb, self.ub).add(solver, False)
            solver.print_constraint("min([%s], %s)" % (csvstr(names), varname(auxvar)))
            return auxvar

        return self


class Minion_Element(Minion_Expression):

    def __init__(self, *args):
        super(Minion_Element, self).__init__()
        self.indexvar = None

        if len(args) == 1:
            self.vars = args[0][:-1]
            self.indexvar = args[0][-1]
        else:
            raise Exception("Invalid constructor to Minion_Element args: %s" % str(args))

        # Compute the lower and upper bound for the auxiliary variable
        lowind = max(0, self.indexvar.get_min())
        highind = min(len(self.vars), self.indexvar.get_max())
        self.lb = min(self.vars[i].get_min() for i in range(lowind, highind))
        self.ub = min(self.vars[i].get_max() for i in range(lowind, highind))

    def add(self, solver, toplevel):
        assert not toplevel, "Constraint is only valid as a sub-expression."
        if not self.has_been_added():
            super(Minion_Element, self).add(solver, toplevel)
            for i in range(len(self.vars)):
                self.vars[i] = self.vars[i].add(solver, False)
            self.indexvar = self.indexvar.add(solver, False)

            names = [varname(x) for x in self.vars]
            auxvar = Minion_IntVar(self.lb, self.ub).add(solver, False)
            solver.print_constraint("watchelement([%s], %s, %s)" % (csvstr(names), varname(self.indexvar), varname(auxvar)))

            return auxvar

        return self


class Minion_Minimise(Minion_Expression):

    def __init__(self, var):
        super(Minion_Minimise, self).__init__()
        self.var = var

    def add(self, solver, toplevel):
        assert toplevel, "Constraint is not valid as a sub-expression."
        if not self.has_been_added():
            super(Minion_Minimise, self).add(solver, toplevel)
            self.var = self.var.add(solver, False)
            solver.print_search("MINIMISING %s" % varname(self.var))
        return self


class Minion_Maximise(Minion_Expression):

    def __init__(self, var):
        super(Minion_Maximise, self).__init__()
        self.var = var

    def add(self, solver, toplevel):
        assert toplevel, "Constraint is not valid as a sub-expression."
        if not self.has_been_added():
            super(Minion_Maximise, self).add(solver, toplevel)
            self.var = self.var.add(solver, False)
            solver.print_search("MAXIMISING %s" % varname(self.var))
        return self


class MinionSolver(ExternalSolver):

    HEADER, VARIABLES, CONSTRAINTS, SEARCH = 0, 1, 2, 3

    def __init__(self):
        super(MinionSolver, self).__init__()
        self.solverexec = "minion"  # executable name to check availability
        self.name_var_map = {}  # Maps an output variable name back to the Variable object
        self.last_section = MinionSolver.HEADER
        self.variablecount = self.constraintcount = 0
        self.cmdlineoptdict = {}

        self.info_regexps = {  # See doc on ExternalSolver.info_regexps
            'nodes': (re.compile(r'^Nodes:[ ]+(?P<nodes>\d+)$'), int),
            'time': (re.compile(r'^Solve Time:[ ]+(?P<time>\d+\.\d+)$'), float),
        }
        self.f = open(self.filename, "wt")
        print("MINION 3", file=self.f)

    def build_solver_cmd(self):
        # The Verbosity that we pass down to the solver should be at least 1 so
        # that we can parse information like number of nodes, failures, etc.
        return "%(solverexec)s %(filename)s " % vars(self) + self.build_cmdlineoptions()

    def build_cmdlineoptions(self):
        s = ""
        for k, v in self.cmdlineoptdict.items():
            s += " " + str(k)
            if v is not None:
                s += " " + str(v)
        return s

    def add(self, expr):
        expr.add(self, True)

    def initialise(self, searchvars=None):
        if searchvars:
            self.print_search("VARORDER [%s]" % csvstr(list(map(varname, searchvars))))

    def solve(self, *args, **kwargs):
        print("**EOF**", file=self.f)
        self.f.close()

        if self.verbosity >= 2:
            self.setOption("-verbose", None)

        if self.verbosity >= 3:
            with open(self.filename, "rt") as f:
                for line in f:
                    print(line, end=' ')

        return super(MinionSolver, self).solve(*args, **kwargs)

    def create_variable(self, localvarobj, name, s):
        self.name_var_map[name] = localvarobj
        self.print_variable(s)

    def print_variable(self, s):
        if self.last_section != MinionSolver.VARIABLES:
            self.last_section = MinionSolver.VARIABLES
            print("**VARIABLES**", file=self.f)  # FIXME switching back and forth
        print(s, file=self.f)

    def print_constraint(self, s):
        if self.last_section != MinionSolver.CONSTRAINTS:
            self.last_section = MinionSolver.CONSTRAINTS
            print("**CONSTRAINTS**", file=self.f)
        self.constraintcount += 1
        print(s, file=self.f)

    def print_search(self, s):
        if self.last_section != MinionSolver.SEARCH:
            self.last_section = MinionSolver.SEARCH
            print("**SEARCH**", file=self.f)
        print(s, file=self.f)

    def getNumVariables(self):
        return self.variablecount

    def getNumConstraints(self):
        return self.constraintcount

    def setRandomSeed(self, seed):
        self.setOption("-randomseed", seed)

    def setNodeLimit(self, limit):
        self.setOption("-nodelimit", limit)

    def setTimeLimit(self, limit):
        self.setOption("-cpulimit", limit)

    def setOption(self, name, value=None):
        self.cmdlineoptdict[name] = value

    def parse_output(self, output):
        minionvarid = 0
        # Assumes variables are printed in the order x1, x2, ...
        for line in output.split("\n"):
            line = line.strip()
            # print repr(line)
            if line.startswith("Sol: "):
                name = "x%d" % (minionvarid)  # Minion variable name
                minionvarid += 1
                val = int(line.split(" ")[-1])
                # print "setting %s to %d" % (name, val)
                assert name in self.name_var_map, "Unknown variable in solver's output %s" % varname
                self.name_var_map[name].value = val

            elif line.startswith("Solutions Found: "):
                sols = int(line.split(" ")[-1])
                if sols > 0:
                    self.sat = SAT
            elif line.startswith("Problem solvable?: "):
                if line.endswith(" no"):
                    self.sat = UNSAT

            elif line.startswith("Solution Number: "):
                minionvarid = 0  # Reset variable counter for the next solution

            else:
                self.parse_solver_info_line(line)


def csvstr(l):
    return ",".join(str(v) for v in l)


def varname(x):
    if isinstance(x, Minion_Expression):
        assert x.varname is not None, "Error varname not set %s %s" % (str(x), str(x.nbj_ident))
        return x.varname
    return str(x)


class Solver(NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        NBJ_STD_Solver.__init__(self, "Minion", "Minion", model, X, FD, clause_limit, encoding)
