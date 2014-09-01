from Numberjack.ExternalSolver import ExternalSolver
from Numberjack import NBJ_STD_Solver, Variable, SAT, UNSAT
import sys
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
        # print "IntVar", len(args), repr(args)
        self.domain = self.value = None

        if len(args) == 0:  # Boolean
            self.lb, self.ub = 0, 1
        elif len(args) == 3:
            self.lb, self.ub, self.nbj_ident = args
        elif len(args) == 2:
            self.domain, self.nbj_ident = args
        else:
            raise Exception("Unknown constructor for %s" % str(type(self)))

    def add(self, solver, toplevel):
        if self.lb == self.ub:
            return str(lb)

        if not self.has_been_added():
            self.varname = "x%d" % (solver.variable_id)  # Minion variable name
            solver.variable_id += 1
            self.solver = solver

            if self.lb == 0 and self.ub == 1:
                solver.create_variable(self, self.varname, "BOOL %s" % self.varname)
            elif self.domain and len(self.domain) != (self.ub - self.lb) + 1:
                dom_str = csvstr(map(str, self.domain))
                solver.create_variable(self, self.varname, "SPARSEBOUND %s {%s}" % (self.varname, dom_str))
            else:
                solver.create_variable(self, self.varname, "DISCRETE %s {%d..%d}" % (self.varname, self.lb, self.ub))
        return self

    def get_value(self):
        return self.value


# class MinionIntArray(list):
#     pass


class MinionExpArray(list):
    add = list.append
    # def add(self, x):
    #     self.append(x)


class Minion_binop(Minion_Expression):

    def __init__(self, arg1, arg2):
    # def __init__(self, minionconname, arg1, arg2):
    #     self.minionconname = minionconname
        self.vars = [arg1, arg2]
        super(Minion_binop, self).__init__()

    def add(self, solver, toplevel):
        if not self.has_been_added():
            self.solver = solver
            # assert toplevel, "Constraint %s not supported as a sub-expression." % self.minionconname
            self.vars[0] = self.vars[0].add(solver, False)
            if isinstance(self.vars[1], Minion_Expression):
                self.vars[1] = self.vars[1].add(solver, False)

            # solver.print_constraint("%s(%s, %s)" % (self.minionconname, varname(self.vars[0]), varname(self.vars[1])))
        return self


class Minion_ne(Minion_binop):

    # def __init__(self, *args):
    #     super(Minion_ne, self).__init__("diseq", *args)

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_ne, self).add(solver, toplevel)
            solver.print_constraint("diseq(%s, %s)" % (varname(self.vars[0]), varname(self.vars[1])))
        return self


class Minion_eq(Minion_binop):

    # def __init__(self, *args):
    #     super(Minion_eq, self).__init__("eq", *args)

    def add(self, solver, toplevel):
        if not self.has_been_added():
            assert toplevel, "Constraint implemented as a sub-expression/refieid yet."
            super(Minion_eq, self).add(solver, toplevel)
            solver.print_constraint("eq(%s, %s)" % (varname(self.vars[0]), varname(self.vars[1])))
        return self


class Minion_lt(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            assert toplevel, "Constraint implemented as a sub-expression/refieid yet."
            super(Minion_eq, self).add(solver, toplevel)
            solver.print_constraint("ineq(%s, %s, -1)" % (varname(self.vars[0]), varname(self.vars[1])))
        return self


class Minion_le(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            assert toplevel, "Constraint implemented as a sub-expression/refieid yet."
            super(Minion_le, self).add(solver, toplevel)
            solver.print_constraint("ineq(%s, %s, 0)" % (varname(self.vars[0]), varname(self.vars[1])))
        return self


class Minion_gt(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            assert toplevel, "Constraint implemented as a sub-expression/refieid yet."
            super(Minion_gt, self).add(solver, toplevel)
            solver.print_constraint("ineq(%s, %s, -1)" % (varname(self.vars[1]), varname(self.vars[0])))
        return self


class Minion_ge(Minion_binop):

    def add(self, solver, toplevel):
        if not self.has_been_added():
            assert toplevel, "Constraint implemented as a sub-expression/refieid yet."
            super(Minion_ge, self).add(solver, toplevel)
            solver.print_constraint("ineq(%s, %s, 0)" % (varname(self.vars[1]), varname(self.vars[0])))
        return self


class Minion_AllDiff(Minion_Expression):

    def __init__(self, children):
        self.vars = children
        super(Minion_AllDiff, self).__init__()

    def add(self, solver, toplevel):
        if not self.has_been_added():
            super(Minion_AllDiff, self).add(solver, toplevel)
            for i in xrange(len(self.vars)):
                self.vars[i] = self.vars[i].add(solver, False)

            if len(self.vars) == 1:
                print "Singular variable alldiff, why??"
                return self.vars[0]
            elif len(self.vars) == 2:
                ne = Minion_ne(*self.vars)
                return ne.add(solver, toplevel)
            else:
                assert toplevel, "Constraint implemented as a sub-expression/refieid yet."

                solver.print_constraint("gacalldiff([%s])" % (csvstr(map(varname, self.vars))))
        return self


# class ExternalIntVariable(object):

#     def __init__(self, nj_var):
#         self.nj_var = nj_var
#         self.value = None

#     def get_value(self):
#         # print "Get value", self.nj_var.name()
#         return self.value

#     def get_min(self):
#         return self.nj_var.lb
#     #     return self.nj_var.get_min(solver=None)

#     def get_max(self):
#         return self.nj_var.ub
#     #     return self.nj_var.get_max(solver=None)

#     def get_size(self):
#         return self.get_max() - self.get_min() + 1


class MinionSolver(ExternalSolver):

    HEADER, VARIABLES, CONSTRAINTS = 0, 1, 2

    def __init__(self):
        super(MinionSolver, self).__init__()
        self.solverexec = "minion"
        self.name_var_map = {}  # Maps an output variable name back to the Variable object
        self.expr_name_map = {}  # Maps an Expression object to the output variable name
        self.last_section = MinionSolver.HEADER
        self.variable_id = 0

        self.info_regexps = {  # See doc on ExternalSolver.info_regexps
            'nodes': (re.compile(r'^Nodes:[ ]+(?P<nodes>\d+)$'), int),
            'time': (re.compile(r'^Solve Time:[ ]+(?P<time>\d+\.\d+)$'), float),
            # 'failures': (re.compile(r'^conflicts[ ]+:[ ]+(?P<failures>\d+)[ ]'), int),
        }
        self.f = open(self.filename, "wt")
        print >> self.f, "MINION 3"

    def build_solver_cmd(self):
        # The Verbosity that we pass down to the solver should be at least 1 so
        # that we can parse information like number of nodes, failures, etc.
        return "%(solverexec)s %(filename)s" % vars(self)

    def add(self, expr):
        # print "Solver add"
        # print type(expr), str(expr)
        expr.add(self, True)

    def initialise(self, searchvars=None):
        print "initialise", str(searchvars)
        pass  # FIXME, could add the search vars

    def solve(self, *args, **kwargs):
        print >> self.f, "**EOF**"
        print "calling solve"
        self.f.close()

        # DEBUG
        with open(self.filename, "rt") as f:
            for line in f:
                print line,

        super(MinionSolver, self).solve(*args, **kwargs)

    # def set_model(self, model, solver_id, solver_name, solver):
    #     self.solver_id = solver_id
    #     self.solver_name = solver_name
    #     self.model = model
    #     self.model.close(solver)
    #     self.output_model()
    #
    # def output_model(self):
    #     print "c Outputting to:", self.filename
    #     with open(self.filename, "wt") as f:
    #         print >> f, "MINION 3"
    #         self.output_variables(f)
    #         self.output_constraints(f)
    #         print >> f, "**EOF**"

    #     # DEBUG
    #     with open(self.filename, "rt") as f:
    #         for line in f:
    #             print line,

    def create_variable(self, localvarobj, name, s):
        self.name_var_map[name] = localvarobj
        self.print_variable(s)

    def print_variable(self, s):
        if self.last_section != MinionSolver.VARIABLES:
            self.last_section = MinionSolver.VARIABLES
            print >> self.f, "**VARIABLES**"  # FIXME switching back and forth
        print >> self.f, s

    def print_constraint(self, s):
        if self.last_section != MinionSolver.CONSTRAINTS:
            self.last_section = MinionSolver.CONSTRAINTS
            print >> self.f, "**CONSTRAINTS**"
        print >> self.f, s

    # def create_variable(self, f, lb, ub, domain=None, v=None):
    #     name = "x%d" % (self.variable_id)  # Minion variable name
    #     self.variable_id += 1
    #     if v:
    #         # Create a wrapper variable that numberjack will call get_value on
    #         # which will be associated with the Minion variable
    #         # my_var = ExternalIntVariable()
    #         my_var = ExternalIntVariable(v)
    #         v.setVar(self.solver_id, self.solver_name, my_var, new_solver=self)
    #         v.solver = self

    #         self.name_var_map[name] = my_var
    #         self.expr_name_map[v] = name

    #     if lb == ub:
    #         return str(lb)
    #     elif lb == 0 and ub == 1:
    #         self.print_variable(f, "BOOL %s" % name)
    #     elif domain and len(domain) != (ub - lb) + 1:
    #         dom_str = ",".join(str(x) for x in domain)
    #         self.print_variable(f, "SPARSEBOUND %s {%s}" % (name, dom_str))
    #     else:
    #         self.print_variable(f, "DISCRETE %s {%d..%d}" % (name, lb, ub))
    #     return name

    # def create_aux_variable(self, f, expr):
    #     print "# creating aux variable for", str(expr)
    #     return self.create_variable(f, expr.lb, expr.ub, v=expr)

    # def output_variable(self, f, v):
    #     return self.create_variable(f, v.lb, v.ub, domain=v.domain_, v=v)

    # def output_variables(self, f):
    #     for v in self.model.variables:
    #         self.output_variable(f, v)

    # def output_constraints(self, f):
    #     for c in self.model.get_exprs():
    #         self.output_expr(f, c, toplevel=True)

    # def output_expr(self, f, e, toplevel=True):
    #     """"
    #         Outputs the expression 'e' to minion format. If toplevel is False
    #         then will reify the expression and return the name of the
    #         auxiliary Boolean variable that it was reified to.
    #     """
    #     def getchildname(x):
    #         # print "# getchildname", type(x), str(x)
    #         if type(x) in [int]:
    #             return str(x)
    #         elif isinstance(x, Variable):
    #             # Return the minion name for this variable
    #             return self.expr_name_map[x]
    #         else:
    #             return self.output_expr(f, x, toplevel=False)
    #             # print >> sys.stderr, "Error need to get the reified version of this constraint."
    #             # sys.exit(1)

    #     op = e.get_operator()
    #     print "# op:", op, toplevel, "e:", e, "children:", e.children
    #     names = [getchildname(child) for child in e.children]
    #     print "#", names

    #     # -------------------- Top level constraints --------------------
    #     if toplevel:

    #         if op == "ne":
    #             self.print_constraint(f, "diseq(%s, %s)" % tuple(names))

    #         elif op == "eq":
    #             self.print_constraint(f, "eq(%s, %s)" % tuple(names))

    #         elif op == "lt":
    #             self.print_constraint(f, "ineq(%s, %s, -1)" % tuple(names))

    #         elif op == "le":
    #             self.print_constraint(f, "ineq(%s, %s, 0)" % tuple(names))

    #         elif op == "gt":
    #             self.print_constraint(f, "ineq(%s, %s, -1)" % (names[1], names[0]))

    #         elif op == "ge":
    #             self.print_constraint(f, "ineq(%s, %s, 0)" % (names[1], names[0]))

    #         elif op == "AllDiff":
    #             self.print_constraint(f, "gacalldiff([%s])" % (csvstr(names)))

    #         elif op == "LeqLex":
    #             self.print_constraint(f, "lexleq([%s], [%s])" % (csvstr(names[:len(names)/2]), csvstr(names[len(names)/2:])))

    #         elif op == "LessLex":
    #             self.print_constraint(f, "lexless([%s], [%s])" % (csvstr(names[:len(names)/2]), csvstr(names[len(names)/2:])))

    #         elif op == "Gcc":
    #             print "# Gcc Parameters:", e.parameters
    #             value_str = csvstr(e.parameters[0])
    #             vec_str = csvstr(self.create_variable(f, l, u) for l, u in zip(e.parameters[1], e.parameters[2]))
    #             self.print_constraint(f, "gccweak([%s], [%s], [%s])" % (csvstr(names), value_str, vec_str))

    #         else:
    #             print >> sys.stderr, "# UNKNOWN top level constraint", op, c
    #             sys.exit(1)

    #     # -------------------- Sub-expressions --------------------
    #     else:
    #         aux_name = self.create_aux_variable(f, e)
    #         if op == "Abs":
    #             self.print_constraint(f, "abs(%s, %s)" % (aux_name, names[0]))

    #         elif op == "div":
    #             self.print_constraint(f, "div(%s, %s, %s)" % (names[0], names[1], aux_name))

    #         elif op == "Element":
    #             self.print_constraint(f, "element([%s], %s, %s)" % (csvstr(names[:-1]), names[-1], aux_name))

    #         elif op == "Sum":
    #             flat_coefs, offset = e.parameters

    #             if offset != 0:
    #                 assert(len(names) == 1, "asdf")  # FIXME
    #                 self.print_constraint(f, "ineq(%s, %s, %d)" % (aux_name, names[0], offset))
    #                 self.print_constraint(f, "ineq(%s, %s, %d)" % (names[0], aux_name, -offset))
    #                 print >> sys.stderr, "Error: translation of Sum with offset not implemented yet."
    #                 print >> sys.stderr, offset, flat_coefs, names
    #                 sys.exit(1)

    #             constantvecstr = csvstr(flat_coefs)
    #             varvecstr = csvstr(names)
    #             if any(lambda x: x != 1 for x in flat_coefs):  # Weighted
    #                 self.print_constraint(f, "weightedsumgeq([%s], [%s], %s)" % (constantvecstr, varvecstr, aux_name))
    #                 self.print_constraint(f, "weightedsumleq([%s], [%s], %s)" % (constantvecstr, varvecstr, aux_name))
    #             else:  # Unweighted FIXME test
    #                 self.print_constraint(f, "sumgeq([%s], %s)" % (varvecstr, aux_name))
    #                 self.print_constraint(f, "sumleq([%s], %s)" % (varvecstr, aux_name))

    #         else:
    #             print >> sys.stderr, "# UNKNOWN sub-expression", op, e
    #             sys.exit(1)
    #         return aux_name

    #     return None

    def parse_output(self, output):
        print "c Parse output"
        minionvarid = 0
        # Assumes variables are printed in
        for line in output.split("\n"):
            line = line.strip()
            # print repr(line)
            if line.startswith("Sol: "):
                name = "x%d" % (minionvarid)  # Minion variable name
                minionvarid += 1
                val = int(line.split(" ")[-1])
                # print "setting %s to %d" % (name, val)

                self.name_var_map[name].value = val
                # value = self.name_var_map[name]
                # var.lb = var.ub = val

            elif line.startswith("Solutions Found: "):
                sols = int(line.split(" ")[-1])
                if sols == 0:
                    self.sat = UNSAT
                elif sols > 0:
                    self.sat = SAT

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
        # self.solver_id = model.getSolverId()
        # self.solver.set_model(model, self.solver_id, self.Library, solver=self)
