from Numberjack.ExternalSolver import ExternalSolver
from Numberjack import NBJ_STD_Solver, Variable, SAT, UNSAT
import sys
import re


class ExternalIntVariable(object):

    def __init__(self, nj_var):
        self.nj_var = nj_var
        self.value = None

    def get_value(self):
        # print "Get value", self.nj_var.name()
        return self.value

    def get_min(self):
        return self.nj_var.lb
    #     return self.nj_var.get_min(solver=None)

    def get_max(self):
        return self.nj_var.ub
    #     return self.nj_var.get_max(solver=None)

    def get_size(self):
        return self.get_max() - self.get_min() + 1


class Minion_IntVar(object):
    pass


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

    def build_solver_cmd(self):
        # The Verbosity that we pass down to the solver should be at least 1 so
        # that we can parse information like number of nodes, failures, etc.
        return "%(solverexec)s %(filename)s" % vars(self)

    def set_model(self, model, solver_id, solver_name, solver):
        self.solver_id = solver_id
        self.solver_name = solver_name
        self.model = model
        self.model.close(solver)
        self.output_model()

    def output_model(self):
        print "c Outputting to:", self.filename
        with open(self.filename, "wt") as f:
            print >> f, "MINION 3"
            self.output_variables(f)
            self.output_constraints(f)
            print >> f, "**EOF**"

        # DEBUG
        with open(self.filename, "rt") as f:
            for line in f:
                print line,

    def print_variable(self, f,  s):
        if self.last_section != MinionSolver.VARIABLES:
            self.last_section = MinionSolver.VARIABLES
            print >> f, "**VARIABLES**"  # FIXME switching back and forth
        print >> f, s

    def print_constraint(self, f, s):
        if self.last_section != MinionSolver.CONSTRAINTS:
            self.last_section = MinionSolver.CONSTRAINTS
            print >> f, "**CONSTRAINTS**"
        print >> f, s

    def create_variable(self, f, lb, ub, domain=None, v=None):
        name = "x%d" % (self.variable_id)  # Minion variable name
        self.variable_id += 1
        if v:
            # Create a wrapper variable that numberjack will call get_value on
            # which will be associated with the Minion variable
            # my_var = ExternalIntVariable()
            my_var = ExternalIntVariable(v)
            v.setVar(self.solver_id, self.solver_name, my_var, new_solver=self)
            v.solver = self

            self.name_var_map[name] = my_var
            self.expr_name_map[v] = name

        if lb == ub:
            return str(lb)
        elif lb == 0 and ub == 1:
            self.print_variable(f, "BOOL %s" % name)
        elif domain and len(domain) != (ub - lb) + 1:
            dom_str = ",".join(str(x) for x in domain)
            self.print_variable(f, "SPARSEBOUND %s {%s}" % (name, dom_str))
        else:
            self.print_variable(f, "DISCRETE %s {%d..%d}" % (name, lb, ub))
        return name

    def create_aux_variable(self, f, expr):
        print "# creating aux variable for", str(expr)
        return self.create_variable(f, expr.lb, expr.ub, v=expr)

    def output_variable(self, f, v):
        return self.create_variable(f, v.lb, v.ub, domain=v.domain_, v=v)

    def output_variables(self, f):
        for v in self.model.variables:
            self.output_variable(f, v)

    def output_constraints(self, f):
        for c in self.model.get_exprs():
            self.output_expr(f, c, toplevel=True)

    def output_expr(self, f, e, toplevel=True):
        """"
            Outputs the expression 'e' to minion format. If toplevel is False
            then will reify the expression and return the name of the
            auxiliary Boolean variable that it was reified to.
        """
        def getchildname(x):
            # print "# getchildname", type(x), str(x)
            if type(x) in [int]:
                return str(x)
            elif isinstance(x, Variable):
                # Return the minion name for this variable
                return self.expr_name_map[x]
            else:
                return self.output_expr(f, x, toplevel=False)
                # print >> sys.stderr, "Error need to get the reified version of this constraint."
                # sys.exit(1)

        op = e.get_operator()
        print "# op:", op, toplevel, "e:", e, "children:", e.children
        names = [getchildname(child) for child in e.children]
        print "#", names

        # -------------------- Top level constraints --------------------
        if toplevel:

            if op == "ne":
                self.print_constraint(f, "diseq(%s, %s)" % tuple(names))

            elif op == "eq":
                self.print_constraint(f, "eq(%s, %s)" % tuple(names))

            elif op == "lt":
                self.print_constraint(f, "ineq(%s, %s, -1)" % tuple(names))

            elif op == "le":
                self.print_constraint(f, "ineq(%s, %s, 0)" % tuple(names))

            elif op == "gt":
                self.print_constraint(f, "ineq(%s, %s, -1)" % (names[1], names[0]))

            elif op == "ge":
                self.print_constraint(f, "ineq(%s, %s, 0)" % (names[1], names[0]))

            elif op == "AllDiff":
                self.print_constraint(f, "gacalldiff([%s])" % (csvstr(names)))

            elif op == "LeqLex":
                self.print_constraint(f, "lexleq([%s], [%s])" % (csvstr(names[:len(names)/2]), csvstr(names[len(names)/2:])))

            elif op == "LessLex":
                self.print_constraint(f, "lexless([%s], [%s])" % (csvstr(names[:len(names)/2]), csvstr(names[len(names)/2:])))

            elif op == "Gcc":
                print "# Gcc Parameters:", e.parameters
                value_str = csvstr(e.parameters[0])
                vec_str = csvstr(self.create_variable(f, l, u) for l, u in zip(e.parameters[1], e.parameters[2]))
                self.print_constraint(f, "gccweak([%s], [%s], [%s])" % (csvstr(names), value_str, vec_str))

            else:
                print >> sys.stderr, "# UNKNOWN top level constraint", op, c
                sys.exit(1)

        # -------------------- Sub-expressions --------------------
        else:
            aux_name = self.create_aux_variable(f, e)
            if op == "Abs":
                self.print_constraint(f, "abs(%s, %s)" % (aux_name, names[0]))

            elif op == "div":
                self.print_constraint(f, "div(%s, %s, %s)" % (names[0], names[1], aux_name))

            elif op == "Element":
                self.print_constraint(f, "element([%s], %s, %s)" % (csvstr(names[:-1]), names[-1], aux_name))

            elif op == "Sum":
                flat_coefs, offset = e.parameters

                if offset != 0:
                    assert(len(names) == 1, "asdf")  # FIXME
                    self.print_constraint(f, "ineq(%s, %s, %d)" % (aux_name, names[0], offset))
                    self.print_constraint(f, "ineq(%s, %s, %d)" % (names[0], aux_name, -offset))
                    print >> sys.stderr, "Error: translation of Sum with offset not implemented yet."
                    print >> sys.stderr, offset, flat_coefs, names
                    sys.exit(1)

                constantvecstr = csvstr(flat_coefs)
                varvecstr = csvstr(names)
                if any(lambda x: x != 1 for x in flat_coefs):  # Weighted
                    self.print_constraint(f, "weightedsumgeq([%s], [%s], %s)" % (constantvecstr, varvecstr, aux_name))
                    self.print_constraint(f, "weightedsumleq([%s], [%s], %s)" % (constantvecstr, varvecstr, aux_name))
                else:  # Unweighted FIXME test
                    self.print_constraint(f, "sumgeq([%s], %s)" % (varvecstr, aux_name))
                    self.print_constraint(f, "sumleq([%s], %s)" % (varvecstr, aux_name))

            else:
                print >> sys.stderr, "# UNKNOWN sub-expression", op, e
                sys.exit(1)
            return aux_name

        return None

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


class Solver(NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        NBJ_STD_Solver.__init__(self, "Minion", "Minion", None, None, FD, clause_limit, encoding)
        self.solver_id = model.getSolverId()
        self.solver.set_model(model, self.solver_id, self.Library, solver=self)
