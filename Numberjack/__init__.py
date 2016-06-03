# Copyright 2009 - 2016 Insight Centre for Data Analytics, UCC
from __future__ import print_function, division


UNSAT, SAT, UNKNOWN, LIMITOUT = 0, 1, 2, 4
LUBY, GEOMETRIC = 0, 1
MAXCOST = 100000000

from .solvers import available_solvers
import weakref
import datetime
import types
import sys
#SDG: extend recursive limit for predicate decomposition
sys.setrecursionlimit(10000)
#SDG: needed by the default eval method in BinPredicate
import operator

# Python3 compatibility fixes
if sys.version_info[0] > 2:
    # Python3 does not have a separate `long` type; instead, it only has `int`.
    # https://www.python.org/dev/peps/pep-0237/
    long = int

val_heuristics = ['Lex', 'AntiLex', 'Random', 'RandomMinMax', 'DomainSplit', 'RandomSplit', 'Promise', 'Impact', 'No', 'Guided']
var_heuristics = ['No', 'MinDomain', 'Lex', 'AntiLex', 'MaxDegree', 'MinDomainMinVal', 'Random', 'MinDomainMaxDegree', 'DomainOverDegree', 'DomainOverWDegree', 'DomainOverWLDegree', 'Neighbour', 'Impact', 'ImpactOverDegree', 'ImpactOverWDegree', 'ImpactOverWLDegree', 'Scheduling']


def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str) and not issubclass(type(el), Expression):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def numeric(x):
    tx = type(x)
    return tx is int or tx is float

def lt_with_none(x, y):
    """
    Emulate the behavior of < in python2, which allows comparison against None with the following semantics:
        None < None: False
        None < int: True
        int < None: False
    """
    if y is None:
        return False
    elif x is None:
        return True
    else:
        return x < y

def gt_with_none(x, y):
    """
    Emulate the behavior of > in python2, which allows comparison against None with the following semantics:
        None > None: False
        None > int: False
        int > None: True
    """
    if x is None:
        return False
    elif y is None:
        return True
    else:
        return x > y

# Numberjack exceptions:


class ConstraintNotSupportedError(Exception):
    """
    Raised if the solver being loaded does not support the constraint, and no
    decomposition is available for the constraint. For example in the case of
    loading a divison expression with a Mixed Integer Programming solver.
    """

    def __init__(self, value, solver=None):
        self.value = value
        self.solver = solver

    def __str__(self):
        return "ERROR: Constraint %s not supported by solver %s and no decomposition is available." % (self.value, self.solver)


class UnsupportedSolverFunction(Exception):
    """
    Raised if a solver does not support a particular API call.
    """

    def __init__(self, solver, func_name, msg=""):
        self.solver = solver
        self.__name__ = func_name
        self.msg = msg

    def __str__(self):
        return "ERROR: The solver %s does not support the function '%s'. %s" % (self.solver, self.__name__, self.msg)


class InvalidEncodingException(Exception):
    """
    Raised if an invalid encoding was specified, for example if no domain
    encoding is turned on.
    """

    def __init__(self, msg=""):
        self.msg = msg

    def __str__(self):
        return "ERROR: Invalid encoding configuration. %s" % self.msg


class InvalidConstraintSpecification(Exception):
    """
    Raised in the case of the invalid use of a constraint.
    """

    def __init__(self, msg=""):
        self.msg = msg

    def __str__(self):
        return "ERROR: Invalid constraint specification. %s" % self.msg


class ModelSizeError(Exception):
    """
    Raised if the size of a model has grown excessively large when decomposing
    some constraints for a solver.
    """

    def __init__(self, value, solver=None):
        self.value = value
        self.solver = solver

    def __str__(self):
        return "ERROR: Model decomposition size too big %s for solver %s." % (self.value, self.solver)


# Numberjack domain and expressions:


class Domain(list):
    def __init__(self, arg1, arg2=None):
        """
        \internal
        This class is used to wrap the domain of variables
        in order to print them and/or iterate over values

        Initialised from a list of values, or a lower and an upper bound
        """
        if arg2 is None:
            list.__init__(self, arg1)
            self.sort()
            self.is_bound = False
        else:
            list.__init__(self, [arg1, arg2])
            self.is_bound = True
        self.current = -1

    def next(self):
        """
        \internal
        Wrapper for __next__ for python2 compatibility
        """
        return self.__next__()

    def __next__(self):
        """
        \internal
        Returns the next value when iterating
        """
        self.current += 1
        if self.is_bound:
            if self[0] + self.current > self[-1]:
                raise StopIteration
            else:
                return self[0] + self.current
        else:
            if self.current >= list.__len__(self):
                raise StopIteration
            else:
                return list.__getitem__(self, self.current)

    def __str__(self):
        """
        \internal
        """
        if self.is_bound:
            lb = self[0]
            ub = self[-1]
            if lb + 1 == ub and type(lb) is int:
                return '{' + str(lb) + ',' + str(ub) + '}'
            else:
                return '{' + str(lb) + '..' + str(ub) + '}'

        def extend(idx):
            x = self[idx]
            y = x
            idx += 1
            while idx < len(self):
                if type(self[idx]) is int and self[idx] == y + 1:
                    y = self[idx]
                else:
                    break
                idx += 1
            return (x, y, idx)

        ret_str = '{'
        idx = 0
        while idx < len(self):
            if idx > 0:
                ret_str += ','
            (x, y, idx) = extend(idx)
            ret_str += str(x)
            if type(x) is int and x + 1 < y:
                ret_str += ('..' + str(y))
            elif x != y:
                ret_str += (',' + str(y))

        return ret_str + '}'


class Expression(object):
    """
    Base class from which all expressions and variables inherit.

    :param str operator: the name of this expression operator or variable name.
    """

    def __init__(self, operator):
        #self.mod = None
        self.ident = -1
        self.operator = operator

        # This is the stuff for maintaining multiple representations of the
        # model among different solvers
        self.var_list = []
        self.encoding = None
        self.solver = None

    def __iter__(self):
        return self.get_domain()

    def get_solver(self):
        """
        Returns the solver with which this expression was last loaded.

        :return: The last loaded solver, or `None` if it has not been loaded
            anywhere.
        :rtype: `NBJ_STD_Solver`
        """
        if getattr(self, 'solver', False):
            return self.solver
        else:
            return None

    def initial(self):
        """
        Returns a string representing the initial domain of the expression. For
        example:

        .. code-block:: python

            var1 = Variable(0, 10)
            print var1.initial()
            >>> x0 in {0..10}

        :return: A String representation of original expression definition
        :rtype: str
        """
        output = self.name()
        if self.domain_ is None:
            output += ' in ' + str(Domain(self.lb, self.ub))
        else:
            output += ' in ' + str(Domain(self.domain_))
        return output

    def domain(self, solver=None):
        """
        Returns a string representing the current domain of the expression.

        :param `NBJ_STD_Solver` solver: If specified, the solver for
            which this expression has been loaded. If not specified, the solver
            that has most recenlty loaded the expression will be used.
        """
        output = self.name() + ' in ' + str(self.get_domain(solver=solver))
        return output

    ## Returns a string containing the value of the expression
    # @param solver Solver from which expression solution will be sourced
    # @return String representation of expression solution
    #
    #    solution(self, solver=None) :- Returns a string representing the current solution
    #    of the expression in the solver specified. If no solver is specified then
    #    the returned string represents the solution to the expression in the solver
    #    that has most recently loaded and solved the expression.
    #
    def solution(self, solver=None):
        """
        .. deprecated:: 1.1
           Instead you should use :func:`get_value`, this function is equivalent
           to calling :func:`str` on that.

        Returns a string containing the solution value of the expression. For a
        native representation of the solution value, use :func:`get_value`

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            solution will be sourced, if `None` then the most recently loaded
            solver is used.
        :return: String representation of the expressions solution.
        :rtype: str
        """
        return str(self.get_value(solver))

    def name(self):
        """
        The name of the variable or the operator if this is an expression.
        """
        return self.operator

    def __str__(self):
        if self.is_built() and self.solver.is_sat():
            return self.solution()
        else:
            return self.domain()

    def is_str(self):
        lb = getattr(self, 'lb', None)
        if lb is not None:
            return not numeric(lb)
        return False

    def getVar(self, solver_id):
        # \internal
        return self.var_list[solver_id - 1]

    def setVar(self, solver_id, solver_name, variable, new_solver=None):
        # \internal
        if (solver_id - 1) < len(self.var_list):
            self.var_list[solver_id - 1] = variable
        else:
            self.var_list.append(variable)

    def has_children(self):
        # \internal
        return hasattr(self, 'children')

    def has_parameters(self):
        # \internal
        return hasattr(self, 'parameters')

    def is_built(self, solver=None):
        if solver is None:
            return len(self.var_list) > 0
        else:
            return solver.solver_id - 1 < len(self.var_list)

    def is_var(self):
        return not issubclass(type(self), Predicate)

    def close(self):
        # \internal
        if self.has_children():
            for child in self.children:
                tc = type(child)
                if tc not in [int, int, float, str, bool]:
                    child.close()

    def get_domain(self, solver=None):
        """
        Creates a new :class:`Domain` instance representing the current domain
        of the expression.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            domain will be sourced, if `None` then the most recently loaded
            solver is used.
        :return: The current domain of the expression.
        :rtype: Domain
        """
        if self.is_built(solver) and not isinstance(self.lb, float):
            if solver is None:
                solver = self.solver
            lb, ub = self.get_min(solver), self.get_max(solver)
            if self.get_size(solver) == (ub - lb + 1):
                dom = list(range(lb, ub + 1))
            else:
                # we should make that more efficient by using the underlying
                # solvers to iterate
                dom = [lb]
                while True:
                    v = solver.next(self, dom[-1])
                    if v <= dom[-1]:
                        break
                    else:
                        dom.append(v)
            return Domain(dom)
        elif self.domain_ is not None:
            return Domain(self.domain_)
        else:
            return Domain(self.lb, self.ub)

    def get_value(self, solver=None):
        """
        The current value of the expression. Should be used to retrieve the
        assigned value of a variable and the value of expressions like the
        objective function.

        In the case of variables, it may not be passed in to the solver if it is
        not involved in a non-trivial constraint. For example, `x <= 1`, `x`
        will not get added (by this constraint alone) if it has a upper bound
        which is less or equal to 1. The variable's lower bound will be returned
        as the value in this case.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            value will be sourced, if `None` then the most recently loaded
            solver is used.
        :return: The current value of the expression.
        :rtype: The same as the original domain, either `int`, `float`, or
            `str`.
        """
        has_value = False
        if self.is_built(solver):
            if self.solver.is_sat():
                has_value = True
        value = None
        # In the case of a variable not being created in the interface, return
        # lb as per the doc above.
        if len(self.var_list) == 0 or \
                (solver and ((solver.solver_id - 1) < len(self.var_list) or
                 (solver.solver_id - 1) == 0)):
            has_value = False
            value = self.lb

        if has_value:
            if solver is not None:
                var = self.var_list[solver.solver_id - 1]
            else:
                var = self.var_list[-1]

            if self.is_str():
                value = self.model.strings[var.get_value()]
            else:
                value = var.get_value()

            if isinstance(self, Variable):
                value = type(self.lb)(value)

        return value

    def get_size(self, solver=None):
        """
        The current size of the expression's domain.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            domain size will be sourced, if `None` then the most recently loaded
            solver is used.
        :return: The size of the expression's domain.
        :rtype: `int`
        """
        if solver is not None:
            if self.is_built(solver):
                return self.var_list[solver.solver_id - 1].get_size()
            else:
                return self.ub - self.lb + 1
        elif self.is_built():
            return self.var_list[-1].get_size()
        else:
            return self.ub - self.lb + 1

    def get_min(self, solver=None):
        """
        Current lower bound of the expression.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            lower bound will be sourced, if `None` then the most recently loaded
            solver is used.
        :return: The current lower bound of the expression.
        :rtype: The same as the original domain, either `int`, `float`, or
            `str`.
        """
        the_min = self.lb
        if solver is not None:
            if self.is_built(solver):
                the_min = self.var_list[solver.solver_id - 1].get_min()
        elif self.is_built():
            the_min = self.var_list[-1].get_min()
        if self.is_str():
            return self.model.strings[the_min]
        return the_min

    def get_max(self, solver=None):
        """
        Current upper bound of variable.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            upper bound will be sourced, if `None` then the most recently loaded
            solver is used.
        :return: The current upper bound of the variable.
        :rtype: The same as the original domain, either `int`, `float`, or
            `str`.
        """
        the_max = self.ub
        if solver is not None:
            if self.is_built(solver):
                the_max = self.var_list[solver.solver_id - 1].get_max()
        elif self.is_built():
            the_max = self.var_list[-1].get_max()
        if self.is_str():
            return self.model.strings[the_max]
        return the_max

    # SDG: methods to access initial lb and ub and domain of the expression or
    # one of its children
    def get_ub(self, child=None):
        if (child == None):
            return self.ub
        else:
            return self.children[child].ub if issubclass(type(self.children[child]), Expression) else self.children[child]

    def get_lb(self, child=None):
        if (child == None):
            return self.lb
        else:
            return self.children[child].lb if issubclass(type(self.children[child]), Expression) else self.children[child]

    # not safe! FIXME
    def get_domain_tuple(self):
        if self.is_str():
            tmp_domain = sorted([self.model.string_map[value] for value in self.domain_])
            return (tmp_domain[0], tmp_domain[len(tmp_domain) - 1], tmp_domain)
        else:
            return (self.lb, self.ub, self.domain_)

    def get_children(self):
        # \internal
        if self.has_children():
            return self.children
        else:
            return None

    def get_operator(self):
        # \internal
        return self.operator

    def __and__(self, pred):
        return And([self, pred])

    def __rand__(self, pred):
        return And([self, pred])

    def __or__(self, pred):
        return Or([self, pred])

    def __ror__(self, pred):
        return Or([self, pred])

    def __add__(self, pred):
        return Sum([self, pred], [1, 1])

    def __radd__(self, pred):
        return Sum([pred, self], [1, 1])

    def __sub__(self, pred):
        var = Sum([self, pred], [1, -1])
        var.name = '(' + str(self) + '-' + str(pred) + ')'
        return var

    def __rsub__(self, pred):
        return Sum([pred, self], [1, -1])

    def __div__(self, pred):
        return Div([self, pred])

    def __rdiv__(self, pred):
        return Div([pred, self])

    def __mul__(self, pred):
        return Mul([self, pred])

    def __rmul__(self, pred):
        return Mul([self, pred])

    def __mod__(self, pred):
        return Mod([self, pred])

    def __rmod__(self, pred):
        return Mod([pred, self])

    def __eq__(self, pred):
        return Eq([self, pred])

    def __ne__(self, pred):
        return Ne([self, pred])

    def __lt__(self, pred):
        return Lt([self, pred])

    def __gt__(self, pred):
        return Gt([self, pred])

    def __le__(self, pred):
        return Le([self, pred])

    def __ge__(self, pred):
        return Ge([self, pred])

    def __neg__(self):
        return Neg([self])

    def __contains__(self, v):
        if self.is_built():
            return self.var_list[-1].contain(v)
        elif self.domain_ is None:
            return self.lb <= v <= self.ub
        else:
            return v in self.domain_

    def __hash__(self):
        return id(self)


class Model(object):
    """
    Model object which stores the variables and constraints. The constraints
    declarations are trees, whose internal nodes are predicates or constraints
    and leaves are variables.

    Model can be initialized with any number of arguments. Each argument will be
    treated as an :class:`Expression` or a list of :class:`Expression` to be
    added into the model. If no argument is given the Model will be initially
    empty. An :class:`Expression` can be subsequently added using the method
    :func:`add` or the operator '+='.
    """

    def __init__(self, *expr):
        ## \internal - List of expressions (predicate trees) that where added to
        #  the model
        self.__expressions = []

        #: A :class:`VarArray`: containing the leaves of the predicate trees.
        self.variables = VarArray([])

        #: A :class:`VarArray`: containing the roots of the predicate trees.
        self.constraints = VarArray([])

        #SDG: initial bounds for optimization problems
        ## Initial upper bound for a minimization problem
        self.upper_bound = MAXCOST
        ## Initial lower bound for a maximization problem
        self.lower_bound = -MAXCOST

        ## \internal - Before giving an expression to a solver, or before
        #  printing it, it needs to
        self.closed = 0

        ## \internal - Every new solver get assigned an unique id
        self.current_id = 0

        ## \internal - Initialise from an expression?
        if len(expr) > 0:
            self.add_prime(expr)

    def getSolverId(self):
        # \internal - generates an ident for each new solver
        self.current_id += 1
        return self.current_id

    def add(self, *expr):
        """Add an expresion, or a list/tuple/dict of expressions to the model.

        :param expr: Any number of (or nested lists of) Expression instances.
        """
        self.add_prime(expr)

    def add_prime(self, expr):
        ## \internal - Used to distinguish between a single Expression and a
        #  list of Expressions
        if issubclass(type(expr), list):
            for exp in expr:
                self.add_prime(exp)
        elif issubclass(type(expr), tuple):
            for exp in expr:
                self.add_prime(exp)
        elif issubclass(type(expr), dict):
            for key in expr:
                self.add_prime(exp[key])
        else:
            self.__expressions.append(expr)
        self.close_exp()

    def __iadd__(self, *expr):
        """Can be used to add an expression or a collection of expressions to
        the model like: `model += expression`

        :param expr: Any number of (or nested lists of) Expression instances.
        """
        self.add_prime(expr)
        return self

    def add_expression(self, exp, level):
        ## \internal - add the Expression tree to the model and assign identifiers to the nodes
        # this expression is new, choose an identifiant for it
        te = type(exp)
        if te not in [int, int, float, str, bool]:
            ## THIS IS BUGGY, WE CANNOT ADD THE SAME VARIABLE TO SEVERAL MODELS
            if exp.ident == -1:
            #if exp.mod != self:
                #exp.mod = self
                if exp.get_children() is None:
                    if exp.is_var():
                        exp.ident = len(self.variables)
                        self.variables.append(exp)
                else:  # it is a constraint
                    exp.ident = -2 - len(self.constraints)
                    self.constraints.append(exp)
                    for child in exp.get_children():
                        self.add_expression(child, level + 1)

    def close_exp(self):
        ## \internal - close() is used to fire up preprocessing requiring
        #  knowledge about the whole model
        for i in range(self.closed, len(self.__expressions)):
            self.add_expression(self.__expressions[i], 0)
        self.closed = len(self.__expressions)

    def close(self, solver=None):
        ## \internal - close() is used to fire up preprocessing requiring
        #  knowledge about the whole model and the solver used.

        #SDG: check if it is an optimization problem
        if not any([issubclass(type(expr), Minimise) or issubclass(type(expr), Maximise) or issubclass(type(expr), CostFunction) for expr in self.__expressions]):
            self.upper_bound = 1
        if solver is not None and solver.Library is 'Toulbar2' and self.upper_bound is not None:
            solver.setOption('updateUb', str(self.upper_bound))

        if self.closed == len(self.__expressions):
            tmp_strings = []
            for var in self.variables:
                if var.is_str():
                    var.model = self
                    for value in var.domain_:
                        tmp_strings.append(value)

            self.strings = sorted(set(tmp_strings))
            self.string_map = {}.fromkeys(self.strings)

            for k in range(len(self.strings)):
                self.string_map[self.strings[k]] = k

        self.closed += 1

        if self.closed == len(self.__expressions) + 1:
            for expr in self.get_exprs():
                expr.close()

        if getattr(solver, "Library", None) == 'Toulbar2':
            #print self   #SDG: VERY USEFUL FOR DEBUGGING
            occur = {}
            def rec_occur(expr):
                if not(issubclass(type(expr), Expression)):
                    return
                if expr.is_var():
                    occur[expr] = occur.get(expr, 0) + 1
                else:
                    for j,subexpr in enumerate(expr.children):
                        rec_occur(subexpr)

            for expr in self.__expressions:
                rec_occur(expr)

            def rec_functional(objexpr, objvar, j, minimization):
                if not(issubclass(type(objvar), Expression)):
                    return
                if objvar.is_var():
                    if objvar.lb==objvar.ub or occur[objvar] != 2:  # avoid replacing non intermediate variables
                        return
                    # replace [Predicate(obj,..),Eq(Sum([obj]+vars, [+-1]+coefs),expr)] by [Predicate(Sum(vars+[expr],[-+]coefs+[-1]),..),Eq(0,0)]
                    # and [Predicate(obj,..),Eq(expr,Sum([obj]+vars, [+-1]+coefs))] by [Predicate(Sum(vars+[expr],[-+]coefs+[-1]),..),Eq(0,0)]
                    objconstr = [expr for expr in self.__expressions if issubclass(type(expr), Eq) and ((issubclass(type(expr.children[0]), Sum) and any([expr.children[0].children[u] is objvar for u in range(len(expr.children[0].children))])) or (issubclass(type(expr.children[1]), Sum) and any([expr.children[1].children[u] is objvar for u in range(len(expr.children[1].children))])))]
                    if (len(objconstr)==1):
                        if issubclass(type(objconstr[0].children[0]), Sum):
                            mysum = 0
                        else:
                            mysum = 1
                        pos = filter(lambda u: objconstr[0].children[mysum].children[u] is objvar, range(len(objconstr[0].children[mysum].children)))[0]
                        coefobj = objconstr[0].children[mysum].parameters[0][pos]
                        if (coefobj != -1 and coefobj != 1):
                            return
                        del objconstr[0].children[mysum].children[pos]
                        del objconstr[0].children[mysum].parameters[0][pos]
                        coefeq = objconstr[0].children[1-mysum]
                        if not issubclass(type(coefeq), int) or coefeq != 0:
                            objconstr[0].children[mysum].children.append(Variable(coefeq,coefeq,str(coefeq)) if issubclass(type(coefeq), int) else coefeq)
                            objconstr[0].children[mysum].parameters[0].append(-1)
                        if (coefobj==1):
                            for u in range(len(objconstr[0].children[mysum].children)):
                                objconstr[0].children[mysum].parameters[0][u] = -objconstr[0].children[mysum].parameters[0][u]
                        objexpr.children[j] = objconstr[0].children[mysum]
                        #print "REPLACE",objvar,"by",objexpr.children[j],"in",objexpr
                        objconstr[0].children[0] = Variable(0,0,'0')
                        objconstr[0].children[1] = 0
                        occur[objvar] -= 2
                        rec_functional(objexpr, objexpr.children[j], j, minimization)
                        if issubclass(type(objexpr.children[j]), Expression): objexpr.children[j].close()
                    else:
                        # replace [Predicate(obj,..),Eq(obj,expr)] by [Predicate(expr,..),Eq(0,0)]
                        # and [Predicate(obj,..),Eq(expr,obj)] by [Predicate(expr,..),Eq(0,0)]
                        objconstr = [expr for expr in self.__expressions if issubclass(type(expr), Eq) and ((expr.children[0] is objvar) or (expr.children[1] is objvar))]
                        if (len(objconstr)==1):
                            if (objconstr[0].children[0] is objvar):
                                objexpr.children[j] = objconstr[0].children[1]
                            else:
                                objexpr.children[j] = objconstr[0].children[0]
                            #print "REPLACE",objvar,"by",objexpr.children[j],"in",objexpr
                            objconstr[0].children[0] = Variable(0,0,'0')
                            objconstr[0].children[1] = 0
                            occur[objvar] -= 2
                            rec_functional(objexpr, objexpr.children[j], j, minimization)
                            if issubclass(type(objexpr.children[j]), Expression): objexpr.children[j].close()
                        elif minimization is not None:
                            # ONLY in the objective function: replace [Predicate(obj,..),Table([obj]+vars,tuples,'support')] by [Predicate(Function(vars,dict),..), Table([],[],'support')]
                            objconstr = [expr for expr in self.__expressions if issubclass(type(expr), Table) and any([expr.children[u] is objvar for u in range(len(expr.children))]) and (expr.parameters[1] == 'support')]
                            if (len(objconstr)==1):
                                pos = filter(lambda u: objconstr[0].children[u] is objvar, range(len(objconstr[0].children)))[0]
                                dictionary = {}
                                for t in objconstr[0].parameters[0]:
                                    mytuple = t[:pos]+t[pos+1:]
                                    if minimization:
                                        dictionary[mytuple] = min(dictionary.get(mytuple, MAXCOST), t[pos])
                                    else:
                                        dictionary[mytuple] = max(dictionary.get(mytuple, -MAXCOST), t[pos])
                                objexpr.children[j] = Function(objconstr[0].children[:pos] + objconstr[0].children[pos+1:], dictionary, MAXCOST if minimization else -MAXCOST)
                                objconstr[0].children = []
                                objconstr[0].parameters = [[],0]
                                occur[objvar] -= 2
                                rec_functional(objexpr, objexpr.children[j], j, minimization)
                                if issubclass(type(objexpr.children[j]), Expression): objexpr.children[j].close()
                else:
                    for j,var in enumerate(objvar.children):
                        rec_functional(objvar, var, j, minimization)

            objexpr = [expr for expr in self.__expressions if issubclass(type(expr), Minimise) or issubclass(type(expr), Maximise)]
            if (len(objexpr)==1 and issubclass(type(objexpr[0].children[0]), Expression) and objexpr[0].children[0].is_var()):
                objvar = objexpr[0].children[0]

                # replace Eq('objective',obj) or Eq('obj',obj) by Eq(0,0)  #SDG: VERY HUGLY!!! (avoid creating an objective variable just for Minizinc output purposes)
                objconstr = [expr for expr in self.__expressions if issubclass(type(expr), Eq) and expr.children[0].is_var() and (expr.children[0].name()=='objective' or expr.children[0].name()=='obj') and (expr.children[1] is objvar)]
                if (len(objconstr)==1):
                    objconstr[0].children[0] = Variable(0,0,'0')
                    objconstr[0].children[1] = 0
                    occur[objvar] -= 1

                # remove first intermediate variables in the objective function
                rec_functional(objexpr[0], objvar, 0, issubclass(type(objexpr[0]), Minimise))


            # remove intermediate variables in the constraints if possible
            for expr in self.__expressions:
                if issubclass(type(expr), Predicate) and not(issubclass(type(expr), (Minimise, Maximise))):
                    for (j,var) in enumerate(expr.children):
                        if not( issubclass(type(expr), Eq) and (issubclass(type(var), Sum) or (issubclass(type(var), Expression) and var.is_var())) ):
                            rec_functional(expr, var, j, None)

            # remove dummy equations or dummy Table
            pos = 0
            while (pos < len(self.__expressions)):
                expr = self.__expressions[pos]
                if (issubclass(type(expr), Eq) and issubclass(type(expr.children[0]), Variable) and issubclass(type(expr.children[1]), int) and expr.children[0].get_lb()==expr.children[0].get_ub()==expr.children[1]) or (issubclass(type(expr), Table) and len(expr.children)==0):
                    del self.__expressions[pos]
                else:
                    pos += 1

            #print self   #SDG: VERY USEFUL FOR DEBUGGING

    def __str__(self):
        ## \internal - print
        mod = 'assign:\n  '
        for var in self.variables:
            mod += var.domain() + '\n  '
        mod += '\nsubject to:\n  '
        for con in self.__expressions:
            mod += con.__str__() + '\n  '
        return mod

    def get_exprs(self):
        ## \internal - return the list of Expressions
        return self.__expressions

    def load(self, solvername, X=None, encoding=None):
        """
        The solver is passed as a string, the corresponding module is imported,
        a Solver object created, initialised, and returned.

        :param str solvername: the name of the solver being loaded. Should be
               one of the modules in :mod:`Numberjack.solvers`.
        :param str X: the decision variables.
        :type X: :class:`list` or :class:`VarArray`
        :param EncodingConfiguration encoding: An :class:`EncodingConfiguration`
            instance defining the default encoding for expressions.
        :type encoding: :class:`EncodingConfiguration`
        :raises ImportError: if the named solver could not be loaded.
        :returns: an instance of a :class:`NBJ_STD_Solver` subclass.
        """
        try:
            solverspkg = "Numberjack.solvers"
            solverstring = "%s.%s" % (solverspkg, solvername)
            lib = __import__(solverstring, fromlist=[solverspkg])
            solver = lib.Solver(self, X, encoding=encoding)
        except ImportError:
            raise ImportError(
                "ERROR: Failed during import, wrong module name? (%s)" %
                solvername)
        return solver

    def solve_with(self, library, encoding=None):
        """
        .. deprecated:: 1.1
           Instead you should use :func:`load` first and call solve on that
           solver object instead.

        The solver is passed as a string, the corresponding module is
        imported, a Solver object created, initialised and called.
        A Solution object (dictionary: var -> val) is returned, if the
        Model is unsatisfiable, or the solver fails to solve it, the
        Solution will be empty (None)
        """
        solver = self.load(library, encoding)
        solver.solve()
        return solver.get_solution()

    def set_upper_bound(self, ub):
        """
        For weighted CSPs, sets the initial upper bound.

        :param int ub: The initial upper bound.
        """
        self.upper_bound = str(ub)

    def set_lower_bound(self, lb):
        """
        For weighted CSPs, sets the initial lower bound.

        :param int ub: The initial lower bound.
        """
        self.lower_bound = lb


class Variable(Expression):
    """
    Creates a new variable. The following tables shows example calls to the
    constructor which results in different kinds of variables.

    .. code-block:: python

        Variable()  # Binary variable
        Variable(N)  # Variable in the domain of {0, N-1}
        Variable('x')  # Binary variable called 'x'
        Variable(N, 'x')  # Variable in the domain of {0, N-1} called 'x'
        Variable(l,u)  # Variable in the domain of {l, u}
        Variable(l,u, 'x')  # Variable in the domain of {l, u} called 'x'
        Variable(list)  # Variable with domain specified as a list
        Variable(list, 'x')  # Variable with domain specified as a list called x

    To create a continuous variable just use float values instead of integer
    values when specifying the domain of the variable.

    .. note::
        Typically, just the Mixed Integer Programming solvers support continuous
        valued variables, other solvers use finite domains.

    The variable's domain may also be specified as a list of strings, and the
    interface between the solver's representation will be handled by Numberjack.
    For example in a map-colouring problem, we may have something like the
    following:

    .. code-block:: python

        v = Variable(['red', 'green', 'blue'])
        # ... load and solve the model
        v.get_value()  # Returns 'red'
    """

    def __init__(self, argopt1=None, argopt2=None, argopt3=None):
        domain = None
        lb = 0
        ub = 1
        name = 'x'

        if argopt3 is not None:
            lb = argopt1
            ub = argopt2
            name = argopt3
        elif argopt2 is not None:
            if type(argopt2) is str:
                if numeric(argopt1):
                    ub = argopt1 - 1
                    lb = type(ub)(lb)  # Ensure lb has the same datatype as ub
                else:
                    domain = sorted(argopt1)
                    lb = domain[0]
                    ub = domain[-1]
                name = argopt2
            else:
                lb = argopt1
                ub = argopt2
        elif argopt1 is not None:
            if type(argopt1) is str:
                name = argopt1
            elif numeric(argopt1):
                ub = argopt1 - 1
                lb = type(ub)(lb)  # Ensure lb has the same datatype as ub
            else:
                domain = sorted(argopt1)
                lb = domain[0]
                ub = domain[-1]

        tlb = type(lb)
        tub = type(ub)
        if tlb not in [int, int, float, str]:
            raise TypeError("Warning lower bound of %s is not an int or a float or a string" % name)
        elif tub not in [int, int, float, str]:
            raise TypeError("Warning upper bound of %s is not an int or a float or a string" % name)
        elif type(name) is not str:
            raise TypeError("Warning name variable is not a string")
        elif lb > ub:
            raise ValueError("Warning lower bound (%r) of %s greater than upper bound (%r)" % (lb, name, ub))

        Expression.__init__(self, name)
        self.domain_ = domain
        self.lb = lb
        self.ub = ub


class VarArray(list):
    """
    A VarArray is a list of :class:`Expression` objects. Various methods are
    overloaded to allow easy declaration, formatted printing, and syntactic
    sugars for modelling. The following tables shows example calls to the
    constructor which results in different kinds of variable arrays.

    .. code-block:: python

        VarArray(l)  # creates an array from a list l
        VarArray(n)  # creates an array of n Boolean variables
        VarArray(n, 'x')  # creates an array of n Boolean variables with
                          # names 'x0..xn-1'
        VarArray(n, m)  # creates an array of n variables with domains [0..m-1]
        VarArray(n, m, 'x')  # creates an array of n variables with domains
                             # [0..m-1] and names 'x0..xn-1'
        VarArray(n, d)  # creates an array of n variables with domains specified
                        # in the list 'd'
        VarArray(n, d, 'x')  # creates an array of n variables with domains
                             # specified in the list 'd' and names 'x0..xn-1'
        VarArray(n, l, u, 'x')  # creates an array of n variables with domains
                                # [l..u] and names 'x0..xn-1'
        VarArray(n, l, u)  # creates an array of n variables with domains [l..u]

    VarArray's allow you to state :class:`Element` and lexicographic ordering
    constraints over a sequence of variables using, respectively the operator
    '[]' and '<', '>', '<=', '>='. For instance, given two VarArray X and Y, and
    an Expression x:

    .. code-block:: python

        X[x]  # returns an Element expression, that is, a variable equal to the
              # xth element of the array X
        X <= Y  # returns a LeqLex constraint between X and Y
    """

    def __init__(self, n, optarg1=None, optarg2=None, optarg3=None):
        domain = None
        if hasattr(n, '__iter__'):
            list.__init__(self, n)
            return
        else:
            lb = 0
            ub = 1
            name = 'x'
            if optarg1 is not None:
                if type(optarg1) is str:
                    name = optarg1
                elif type(optarg2) is int or type(optarg2) is float:
                    lb = optarg1
                    ub = optarg2
                    if optarg3 is not None:
                        name = optarg3
                else:
                    if issubclass(type(optarg1), list) or issubclass(type(optarg1), list):
                        domain = optarg1
                        domain.sort()
                        lb = domain[0]
                        ub = domain[-1]
                    else:
                        ub = optarg1 - 1
                    if optarg2 is not None:
                        name = optarg2
        names = name
        if type(name) is str:
            names = [name + str(i) for i in range(int(n))]
        if domain is None:
            self.__init__([Variable(lb, ub, names[i]) for i in range(n)])
        else:
            self.__init__([Variable(domain, names[i]) for i in range(n)])

    def initial(self):
        """
        Returns a string representing the initial definition of the content of
        the array.

        :rtype: str
        """
        return "[" + ", ".join([var.initial() for var in self]) + "]"

    def domain(self, solver=None):
        """
        Returns a string representing the current state of the content of the
        array.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            state will be sourced, if `None` then the most recently loaded
            solver is used.
        :rtype: str
        """
        return "[" + ", ".join([var.domain(solver) for var in self]) + "]"

    def name(self):
        """
        Returns a string containing a brief view of the content of the array.

        :rtype: str
        """
        return "[" + ", ".join([var.name() for var in self]) + "]"

    def solution(self, solver=None):
        """
        .. deprecated:: 1.1
           Instead you should use :func:`Expression.get_value` on each item and
           call :func:`str` on that.

        Returns a string containing the valuation of the content of the array.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            state will be sourced, if `None` then the most recently loaded
            solver is used.
        :rtype: str
        """
        return "[" + ", ".join([var.solution(solver) for var in self]) + "]"

    def __str__(self):
        return "[" + ", ".join([var.__str__() for var in self]) + "]"

    def __getitem__(self, expr):
        if type(expr) is int:
            return list.__getitem__(self, expr)
        else:
            return Element(self, expr)

    def __getslice__(self, i, j):
        return VarArray(list.__getslice__(self, i, j))

    def __lt__(self, other):
        """
        Syntactic sugar for the lexicographic order constraint :class:`LessLex`
        so it can be specified on two VarArray like so: `X < Y`

        :param VarArray other: Another VarArray of the same length.
        :rtype: LessLex
        """
        return LessLex(self, other)

    def __le__(self, other):
        """
        Syntactic sugar for the lexicographic order constraint :class:`LeqLex`
        so it can be specified on two VarArray like so: `X <= Y`

        :param VarArray other: Another VarArray of the same length.
        :rtype: LeqLex
        """
        return LeqLex(self, other)

    def __gt__(self, other):
        """
        Syntactic sugar for the lexicographic order constraint :class:`LessLex`
        so it can be specified on two VarArray like so: `X > Y`

        :param VarArray other: Another VarArray of the same length.
        :rtype: LessLex
        """
        return LessLex(other, self)

    def __ge__(self, other):
        """
        Syntactic sugar for the lexicographic order constraint :class:`LeqLex`
        so it can be specified on two VarArray like so: `X >= Y`

        :param VarArray other: Another VarArray of the same length.
        :rtype: LeqLex
        """
        return LeqLex(other, self)

    def __eq__(self, other):
        """
        Syntactic sugar for the equality constraint `X == Y`.

        :param VarArray other: Another VarArray of the same length.
        :rtype: A list of equality (:class:`Eq`) expressions.
        """
        return [Eq((x, y)) for x, y in zip(self, other)]


class Matrix(list):
    """
    A Matrix is a two-dimensional list of variables or
    :class:`.Expression` objects. Various methods are overloaded to
    allow easy declaration, formatted printing, and syntactic sugars for
    modelling. The following tables shows example calls to the constructor
    which results in different kinds of matrices.

    .. code-block:: python

        Matrix(l)  # creates a Matrix from a list l
        Matrix(n, m)  # creates a n x m Matrix of Boolean variables
        Matrix(n, m, 'x')  # creates a n x m Matrix of Boolean variables with
                           # names 'x0.0..xn-1.m-1'
        Matrix(n, m, u)  # creates a n x m Matrix of variables with domains
                         # [0..u-1]
        Matrix(n, m, u, 'x')  # creates a n x m Matrix of variables with
                              # domains [0..u-1] and names 'x0.0..xn-1.m-1'
        Matrix(n, m, l, u)  # creates a n x m Matrix of variables with domains
                            # [l..u]
        Matrix(n, m, l, u, 'x')  # creates a n x m Matrix of variables with
                                 # domains [l..u] and names 'x0.0..xn-1.m-1'

    Matrices feature specific handlers to access (subsets of) rows and columns.
    The fields `row`, `col`, and `flat` respectively refer to the list of rows,
    columns and cells in the matrix. For instance:

    .. code-block:: python

        m = Matrix(5,4,1,3,'cell_')
        print m
        >>> [[cell_0.0, cell_0.1, cell_0.2, cell_0.3],
        >>>  [cell_1.0, cell_1.1, cell_1.2, cell_1.3],
        >>>  [cell_2.0, cell_2.1, cell_2.2, cell_2.3],
        >>>  [cell_3.0, cell_3.1, cell_3.2, cell_3.3],
        >>>  [cell_4.0, cell_4.1, cell_4.2, cell_4.3]]
        print m.row
        >>> [[cell_0.0, cell_0.1, cell_0.2, cell_0.3],
        >>>  [cell_1.0, cell_1.1, cell_1.2, cell_1.3],
        >>>  [cell_2.0, cell_2.1, cell_2.2, cell_2.3],
        >>>  [cell_3.0, cell_3.1, cell_3.2, cell_3.3],
        >>>  [cell_4.0, cell_4.1, cell_4.2, cell_4.3]]
        print m.col
        >>> [[cell_0.0, cell_1.0, cell_2.0, cell_3.0, cell_4.0],
        >>>  [cell_0.1, cell_1.1, cell_2.1, cell_3.1, cell_4.1],
        >>>  [cell_0.2, cell_1.2, cell_2.2, cell_3.2, cell_4.2],
        >>>  [cell_0.3, cell_1.3, cell_2.3, cell_3.3, cell_4.3]]
        print m.flat
        >>> [cell_0.0, cell_0.1, cell_0.2, cell_0.3, cell_1.0, cell_1.1, ...]

    Matrices support Element constraints on row, column or flatten views.
    """

    def __init__(self, optarg1=None, optarg2=None, optarg3=None, optarg4=None, optarg5=None):
        n = 1
        m = 1
        lb = 0
        ub = 1
        name = 'x'

        self.row = None   # accessor to the list of rows
        self.col = None   # accessor to the list of columns
        self.flat = None  # accessor to the list of cells

        if optarg2 == None:
            if optarg1 != None:
                # BH: This could create rows with varying numbers of columns if given a list with different values.
                #     Should this be allowed? If so, then we need to verify any assumptions being made in this code.
                list.__init__(self, [VarArray(row, "%s%d." % (name, i)) for i, row in enumerate(optarg1)])
            else:
                list.__init__(self)
                return
        else:
            n = optarg1
            m = optarg2
            if optarg5 is not None:
                lb = optarg3
                ub = optarg4
                name = optarg5
            elif optarg4 is not None:
                if type(optarg4) is str:
                    name = optarg4
                    ub = optarg3 - 1
                else:
                    ub = optarg4
                    lb = optarg3
            elif optarg3 is not None:
                if type(optarg3) is str:
                    name = optarg3
                else:
                    ub = optarg3 - 1
            list.__init__(self, [VarArray(m, lb, ub, name + str(j) + '.') for j in range(int(n))])
        self.row = self
        self.col = Matrix()
        for column in zip(*self):
            self.col.append(VarArray(column))
        self.col.col = self
        self.col.row = self.col
        self.flat = VarArray([var for row in self for var in row])
        self.col.flat = self.flat

    def initial(self):
        """
        Returns a string representing the initial definition of the content of
        the matrix.

        :rtype: str
        """
        return "[" + ",\n ".join([row.initial() for row in self]) + "]"

    def domain(self, solver=None):
        """
        Returns a string representing the current state of the content of the
        matrix.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            state will be sourced, if `None` then the most recently loaded
            solver is used.
        :rtype: str
        """
        return "[" + ",\n ".join([row.domain(solver) for row in self]) + "]"

    def name(self):
        """
        Returns a string containing a brief view of the content of the matrix.

        :rtype: str
        """
        return "[" + ",\n ".join([row.name() for row in self]) + "]"

    def solution(self, solver=None):
        """
        .. deprecated:: 1.1
           Instead you should use :func:`Expression.get_value` on each item and
           call :func:`str` on that.

        Returns a string containing the valuation of the content of the array.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            state will be sourced, if `None` then the most recently loaded
            solver is used.
        :rtype: str
        """
        return "[" + ",\n ".join([row.solution(solver) for row in self]) + "]"

    def __str__(self):
        return "[" + ",\n ".join([row.__str__() for row in self]) + "]"

    def __getitem__(self, i):
        if type(i) is int:
            return list.__getitem__(self, i)
        elif type(i) is tuple:
            if type(i[0]) is int:
                return list.__getitem__(self.row, i[0]).__getitem__(i[1])
            elif type(i[1]) is int:
                return list.__getitem__(self.col, i[1]).__getitem__(i[0])
            elif type(i[0]) is slice:
                aux = Matrix(list.__getitem__(self, i[0])).col
                aux = Matrix(list.__getitem__(aux, i[1])).col
                return aux
            else:
                return Element(self.flat, (i[0] * len(self.col)) + i[1])
        elif type(i) is slice:
            return list.__getitem__(self.row, i)
        else:
            return MatrixWrapper(i, self)

    def __getslice__(self, i, j):
        return Matrix(list.__getslice__(self, i, j))


class MatrixWrapper(list):

    def __init__(self, var, matrix):
        self.var = var
        self.matrix = matrix

    def __getitem__(self, item):
        return self.matrix[self.var, item]

    def __str__(self):
        return str(self.var) + " th index of " + str(self.matrix)


## Class that all constraints inherit from
#    All constraints in Numberjack extend the Predicate class. It provides
#    accessors to get information about the predicate trees and the variables
#    the constraints constrain.
#
#    A given predicate can have a different meaning when posted at the top-level
#    or nested in a predicate tree. For instance:
#
# \code
#    from Numberjack import *
#    x = Variable(1,5)
#    y = Variable(1,4)
#    x_lt_y = (x<y)
#    m1 = Model( x_lt_y )
#    print m1
#    >>> assign:
#    >>>   x0 in {1..5}
#    >>>   x1 in {1..4}
#    >>>
#    >>> subject to:
#    >>>   (x0 < x1)
#
#    x = Variable(1,5)
#    y = Variable(1,4)
#    x_gt_y = (x>y)
#    m2 = Model( x_lt_y | x_gt_y )
#    print m2
#    >>> assign:
#    >>>   x0 in {1..5}
#    >>>   x1 in {1..4}
#    >>>
#    >>> subject to:
#    >>>   ((x0 < x1) or (x0 > x1))
# \endcode
#
#    - In the first Model (m1), the object x_lt_y is understood as a precedence Constraint
#      between the Variables x and y
#
#    - In the second Model (m2), the same object x_lt_y is understood as a Boolean variable,
#      whose truth value corresponds to the relation (x<y) and that can be constrained,
#      here with an "Or" constraint.
class Predicate(Expression):
    """
    All constraints in Numberjack extend the this class. It provides accessors
    to get information about the predicate trees and the variables the
    constraints constrain.
    """  # FIXME add doc about creating subclasses for custom constraints

    def __init__(self, vars, op):
        Expression.__init__(self, op)
        self.set_children(vars)

    def set_children(self, children):
        #self.children = children

        ## List of children of the predicate
        #self.children = [child for child in children]
        self.children = flatten(children)

    def initial(self):
        """
        Returns a string representing the initial definition of the content of
        the predicate and all its children. For example:

        .. code-block:: python

            var1 = Variable(0, 10)
            constraint = var1 < 10
            print constraint.initial()
            >>> (x0 < 10)

        :rtype: str
        """
        save_str = Expression.__str__
        Expression.__str__ = Expression.initial
        output = self.__str__()
        Expression.__str__ = save_str
        return output

    def domain(self, solver=None):
        """
        Returns a string representing the current domain of the predicate.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            state will be sourced, if `None` then the most recently loaded
            solver is used.
        :rtype: str
        """
        save_str = Expression.__str__
        Expression.__str__ = lambda x: x.domain(solver)
        output = self.__str__()
        Expression.__str__ = save_str
        return output

    def solution(self, solver=None):
        """
        .. deprecated:: 1.1
           Instead you should use :func:`Expression.get_value` and call
           :func:`str` on that.

        Returns a string containing the valuation of the predicate.

        :param `NBJ_STD_Solver` solver: If specified, the solver from which the
            state will be sourced, if `None` then the most recently loaded
            solver is used.
        :rtype: str
        """
        save_str = Expression.__str__
        Expression.__str__ = lambda x: x.solution(solver)
        output = self.__str__()
        Expression.__str__ = save_str
        return output

    def name(self):
        """
        Returns a string that represents the name of the Predicate and the name
        of all its children.

        :rtype: str
        """
        return self.__str__()

    def __str__(self):
        save_str = Expression.__str__
        Expression.__str__ = Expression.name
        output = self.operator + "(" + ", ".join(map(str, self.children)) + ")"
        Expression.__str__ = save_str
        return output


class BinPredicate(Predicate):
    """
    All binary predicates such as And, LessThan and GreaterThan extend this
    class. They are separated from the base Predicate class to facilitate
    ease of print representations of the predicates.
    """

    def __init__(self, vars, op):
        Predicate.__init__(self, vars, op)

    def get_symbol(self):
        return 'x'

    def __str__(self):
        save_str = Expression.__str__
        Expression.__str__ = Expression.name
        output = '(' + str(self.children[0]) + ' ' + self.get_symbol() + ' ' + str(self.children[1]) + ')'
        Expression.__str__ = save_str
        return output

    #SDG: generic eval method using Predicate operator name
    def eval(self, x,y):
        try:
            return int(getattr(operator, self.operator)(x,y))
        except AttributeError:
            return int(eval(str(x) + ' ' + self.get_symbol() + ' ' + str(y)))


class And(BinPredicate):
    """
    Logical 'and' expression. Can be used at both the top-level to specify an
    'and' constraint, or reified to equal the truth value of the relation.

    .. code-block:: python

        var1 = Variable()  # Binary variable
        var2 = Variable()  # Binary variable

        model.add(var1 & var2)  # Post var1 And var2 constraint

        var1 = Variable()  # Binary variable
        var2 = Variable()  # Binary variable
        var3 = Variable()  # Binary variable

        model.add( var3 == (var1 & var2) )  # Used as an expression
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "and")
        self.lb = min(self.get_lb(0), self.get_lb(1))
        self.ub = min(self.get_ub(0), self.get_ub(1))

    def get_symbol(self):
        return '&'


class Or(BinPredicate):
    """
    Logical 'or' expression. Can be used at both the top-level to specify an
    'or' constraint, or reified to equal the truth value of the relation.

    .. code-block:: python

        var1 = Variable()  # Binary variable
        var2 = Variable()  # Binary variable

        model.add(var1 | var2)  # Post var1 Or var2 constraint

        var1 = Variable()  # Binary variable
        var2 = Variable()  # Binary variable
        var3 = Variable()  # Binary variable

        model.add( var3 == (var1 | var2) )  # Used as an expression
    """
    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "or")
        self.lb = max(self.get_lb(0), self.get_lb(1))    #SDG: initialize lb,ub
        self.ub = max(self.get_ub(0), self.get_ub(1))

    def get_symbol(self):
        return 'or'


class Div(BinPredicate):
    """
    Division expression to equal the integral division of the two operands.
    Cannot be used as a top-level constraint.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)

        divexp1 = var2 / var1
        divexp2 = var2 / 10

    .. warning::

        Cannot be used with all solvers and
        :class:`.ConstraintNotSupportedError` will be raised when loading the
        model if this is the case.
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "div")
        #SDG: initialize lb,ub
        if (self.get_lb(1) < 0 and self.get_ub(1) > 0):
            self.lb = min(self.get_lb(0), -1 * self.get_ub(0))   #SDG: Warning! It assumes var2 can be equal to -1 or 1
            self.ub = max(self.get_ub(0), -1 * self.get_lb(0))
        elif (self.get_ub(1) < 0):
            self.lb = min(self.get_lb(0) // self.get_ub(1), self.get_ub(0) // self.get_ub(1))
            self.ub = max(self.get_lb(0) // self.get_ub(1), self.get_ub(0) // self.get_ub(1))
        elif (self.get_lb(1) > 0):
            self.lb = min(self.get_lb(0) // self.get_lb(1), self.get_ub(0) // self.get_lb(1))
            self.ub = max(self.get_lb(0) // self.get_lb(1), self.get_ub(0) // self.get_lb(1))
        else:
            self.lb = None
            self.ub = None

    def get_symbol(self):
        return '/'


class Mul(BinPredicate):
    """
    Multiplication expression to equal the multiplication of the two operands.
    Cannot be used as a top-level constraint.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)

        mulxp1 = var2 * var1
        mulexp2 = var2 * 10

    .. warning::

        Cannot be used with all solvers and
        :class:`.ConstraintNotSupportedError` will be raised when loading the
        model if this is the case.
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "mul")
        #SDG: initialize lb,ub
        self.lb = min(self.get_lb(0) * self.get_lb(1), self.get_lb(0) * self.get_ub(1), self.get_ub(0) * self.get_lb(1), self.get_ub(0) * self.get_ub(1))
        self.ub = max(self.get_lb(0) * self.get_lb(1), self.get_lb(0) * self.get_ub(1), self.get_ub(0) * self.get_lb(1), self.get_ub(0) * self.get_ub(1))

    def get_symbol(self):
        return '*'


class Mod(BinPredicate):
    """
    Modulus expression to equal the modulo two expressions or an expression and
    a constraint. Cannot be used as a top-level constraint.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)

        modexp1 = var2 % var1
        modexp2 = var2 % 10

    .. note::

        For MIP and SAT, the constraint is encoded such that the remainder takes
        the sign of the numerator, as per the C standard. This differs from
        Python where the remainder takes the sign of the denominator.

    .. warning::

        Cannot be used with all solvers and
        :class:`.ConstraintNotSupportedError` will be raised when loading the
        model if this is the case.
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "mod")
        #SDG: initialize lb,ub
        if (self.get_ub(1) > 0):
            self.lb = 0
            self.ub = self.get_ub(1) - 1
        else:
            self.lb = None   #SDG: Warning! lb and ub undefined if var2 can be negative
            self.ub = None

    def get_symbol(self):
        return '%'


class Eq(BinPredicate):
    """
    Equality expression between two expressions, or an expression and a
    constant. It can be used as either a top-level hard constraint or reified
    as a sub-expression.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)
        var3 = Variable()

        model.add(var1 == var2)
        model.add(var1 == 5)

        model.add( var3 == (var1 == var2) )
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "eq")
        #SDG: initialize lb,ub
        self.lb = int((self.get_lb(0) == self.get_ub(0)) and (self.get_lb(1) == self.get_ub(1)) and (self.get_lb(0) == self.get_lb(1)))
        self.ub = int(not(lt_with_none(self.get_ub(0), self.get_lb(1)) or gt_with_none(self.get_lb(0), self.get_ub(1))))

    def get_symbol(self):
        return '=='


class Ne(BinPredicate):
    """
    Disequality expression between two expressions, or an expression and a
    constant. It can be used as either a top-level hard constraint or reified as
    a sub-expression.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)
        var3 = Variable()

        model.add(var1 != var2)
        model.add(var1 != 5)

        model.add( var3 != (var1 != var2) )
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "ne")
        #SDG: initialize lb,ub
        self.lb = int((self.get_ub(0) < self.get_lb(1)) or (self.get_lb(0) > self.get_ub(1)))
        self.ub = int(not((self.get_lb(0) == self.get_ub(0)) and (self.get_lb(1) == self.get_ub(1)) and (self.get_lb(0) == self.get_lb(1))))

    def get_symbol(self):
        return '!='  #SDG: operator '=/=' does not belong to python language


class Lt(BinPredicate):
    """
    Less-than expression between two expressions, or an expression and a
    constant. It can be used as either a top-level hard constraint or reified as
    a sub-expression.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)
        var3 = Variable()

        model.add(var1 < var2)
        model.add(var1 < 5)

        model.add( var3 < (var1 < var2) )
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "lt")
        #SDG: initialize lb,ub
        self.lb = int(self.get_ub(0) < self.get_lb(1))
        self.ub = int(not(self.get_lb(0) >= self.get_ub(1)))

    def get_symbol(self):
        return '<'


class Gt(BinPredicate):
    """
    Greater-than expression between two expressions, or an expression and a
    constant. It can be used as either a top-level hard constraint or reified as
    a sub-expression.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)
        var3 = Variable()

        model.add(var1 > var2)
        model.add(var1 > 5)

        model.add( var3 > (var1 > var2) )
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "gt")
        #SDG: initialize lb,ub
        self.lb = int(self.get_lb(0) > self.get_ub(1))
        self.ub = int(not(self.get_ub(0) <= self.get_lb(1)))

    def get_symbol(self):
        return '>'


class Le(BinPredicate):
    """
    Less than or equal expression between two expressions, or an expression and
    a constant. It can be used as either a top-level hard constraint or reified
    as a sub-expression.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)
        var3 = Variable()

        model.add(var1 <= var2)
        model.add(var1 <= 5)

        model.add( var3 <= (var1 <= var2) )
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "le")
        #SDG: initialize lb,ub
        self.lb = int(self.get_ub(0) <= self.get_lb(1))
        self.ub = int(not(self.get_lb(0) > self.get_ub(1)))

    def get_symbol(self):
        return '<='


class Ge(BinPredicate):
    """
    Greater than or equal expression between two expressions, or an expression
    and a constant. It can be used as either a top-level hard constraint or
    reified as a sub-expression.

    .. code-block:: python

        var1 = Variable(0, 10)
        var2 = Variable(0, 100)
        var3 = Variable()

        model.add(var1 >= var2)
        model.add(var1 >= 5)

        model.add( var3 >= (var1 >= var2) )
    """

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "ge")
        #SDG: initialize lb,ub
        self.lb = int(self.get_lb(0) >= self.get_ub(1))
        self.ub = int(not(self.get_ub(0) < self.get_lb(1)))

    def get_symbol(self):
        return '>='


class Neg(Predicate):
    """
    Negate expression, used to negate another expression. It is equivalent to
    multiplying by -1.

    .. code-block:: python

        var = Variable(1, 10)
        model.add(-var < 3)
    """

    def __init__(self, vars):
        Predicate.__init__(self, vars, "neg")
        self.lb = -(self.get_ub(0))      #SDG: initialize lb/ub
        self.ub = -(self.get_lb(0))

    def get_min(self, solver=None):
        return -1 * self.children[0].get_max(solver)

    def get_max(self, solver=None):
        return -1 * self.children[0].get_min(solver)

    def __str__(self):
        return '-' + str(self.children[0])

    def decompose(self):
        return [self.children[0] * -1]


class Abs(Predicate):
    """
    Absolute expression, represents the absolute value of an expression or
    variable.

    .. code-block:: python

        var = Variable(-5, 5)
        model.add(Abs(var) < 3)
    """

    def __init__(self, vars):
        Predicate.__init__(self, [vars], "Abs")
        if (self.get_lb(0) < 0 and self.get_ub(0) > 0):      #SDG: initialize lb/ub
            self.lb = 0
            self.ub = max(abs(self.get_lb(0)),abs(self.get_ub(0)))
        else:
            self.lb = min(abs(self.get_lb(0)),abs(self.get_ub(0)))
            self.ub = max(abs(self.get_lb(0)),abs(self.get_ub(0)))

    def __str__(self):
        return "Abs(" + str(self.children[0]) + ")"

    def decompose(self):
        return [Max([self.children[0], Neg([self.children[0]])])]


class Table(Predicate):
    """
    Table constraint explicityly specifying the list of allowed or forbidden
    tuples. Must be used as a top-level constraint, it cannot be used as a
    sub-expression.

    :param list vars: the variables to be constrained by the constraint.
    :param tuples: list of tuples used for the table constraint.
    :param type: type of table constraint, either support of conflict
    """

    def __init__(self, vars, tuples=[], type='support'):
        Predicate.__init__(self, vars, "Table")
        self.parameters = [[tuple for tuple in tuples], type]
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    def addSupport(self, tuple):
        """
        Adds in a support tuple to the list of tuples, if the Table constraint is
        of type support. If the table constraint is of type conflict then if the
        provided tuple is a conflict tuple in the tuple list, it is removed.

        :param tuple: tuple to be added.
        """
        if self.parameters[1] == 'support':
            self.parameters[0].append(tuple)
        else:
            self.parameters[0].remove(tuple)

    def addConflict(self, tuple):
        """
        Adds in a conflict tuple to the list of tuples, if the Table constraint
        is of type conflict. If the table constraint is of type support then if
        the provided tuple is a support tuple in the tuple list, it is removed.

        :param tuple: tuple to be added.
        """
        if self.parameters[1] != 'support':
            self.parameters[0].append(tuple)
        else:
            self.parameters[0].remove(tuple)

    def printTable(self):
        """
        Prints the table of tuples to standard output.
        """
        print(self.parameters[1])
        for var in self.children:
            print(var, end=' ')
        print('\n (' + self.parameters[1] + ')', self.parameters[0])

    def __str__(self):      #SDG: pretty print of Table Predicate
        return self.operator + "([" + ",".join([str(var) for var in self.children]) + "]," + str(self.parameters[0]) + ",'" + self.parameters[1] + "')"


class Sum(Predicate):
    """
    Sum expression with linear coefficients. Numberjack will detect inline sum
    expressions to extract the sum expression. For example, the following three
    statements are equivalent however the last one requires the least amount of
    overhead by Numberjack.

    .. code-block:: python

        2*a + b + 0.5*c + 3*d == e
        Sum([2*a, b, 0.5*c, 3*d]) == e
        Sum([a,b,c,d], [2, 1, 0.5, 3]) == e

    .. note::

        Cannot be used as a top-level constraint, but can be used as the objective funcion.

    :param vars: the variables to be summed.
    :param coefs: list of coefficients, which is [1,1,..,1] by default.
    """

    def __init__(self, vars, coefs=None, offset=0):
        Predicate.__init__(self, vars, "Sum")

        if coefs is None:
            coefs = [1 for var in self.children]

        self.parameters = [coefs, offset]
        #SDG: initial bounds
        self.lb = sum(c*self.get_lb(i) if (c >= 0) else c*self.get_ub(i) for i,c in enumerate(coefs)) + offset
        self.ub = sum(c*self.get_ub(i) if (c >= 0) else c*self.get_lb(i) for i,c in enumerate(coefs)) + offset

    def close(self):
        # This handles the scalar constraint, i.e. with weights

        Predicate.close(self)

        def extract_sum(var, coef):

            """
                Function that extracts A + B + C into SUM(A,B,C) also works with coefs
                It's kinda clever
            """

            if hasattr(var, 'operator'):
                if var.operator == "Sum":
                    res = []
                    for (s_var, s_coef) in zip(var.get_children(), var.parameters[0]):
                        res.extend(extract_sum(s_var, s_coef * coef))
                    if var.parameters[1] != 0:
                        res.append((var.parameters[1], 1))
                    return res
                elif var.operator == "mul":
                    if (var.get_children()[0].operator == "var" and
                        (type(var.get_children()[1]) == int or
                         type(var.get_children()[1]) == float)):
                        return [(var.get_children()[0], var.get_children()[1] * coef)]
                    elif (type(var.get_children()[1]) == int or
                         type(var.get_children()[1]) == float):

                        return [(new_var, new_coef * var.get_children()[1] * coef)
                            for (new_var, new_coef) in extract_sum(var.get_children()[0], 1)]

                    else:
                        return [(var, 1 * coef)]  # It is quadratic?
                else:
                    return [(var, 1 * coef)]
            else:
                return [(var, 1 * coef)]

        # This is where it should go looking for +s
        set_vars = set([])
        list_vars = []
        map_coefs = {}
        offset = self.parameters[1]

        for (var, coef) in zip(self.children, self.parameters[0]):
            list = extract_sum(var, coef)
            for (nVar, nCoef) in list:
                if type(nVar) == int or type(nVar) == float:
                    offset += (nVar * nCoef)
                else:
                    if nVar in set_vars:
                        map_coefs[nVar] += nCoef
                    else:
                        set_vars.add(nVar)
                        list_vars.append(nVar)
                        map_coefs[nVar] = nCoef

        flat_vars = []
        flat_coefs = []
        for nVar in list_vars:
            if map_coefs[nVar] != 0:
                flat_vars.append(nVar)
                flat_coefs.append(map_coefs[nVar])

        self.set_children(flat_vars)
        self.parameters = [flat_coefs, offset]

    def __str__(self):
        #print len(self.children)
        op = '('
        if len(self.parameters[0]):
            if self.parameters[0][0] != 1:
                op += (str(self.parameters[0][0]) + '*')
            op += (self.children[0].__str__())
        for i in range(1, len(self.children)):
            if self.parameters[0][i] == 1:
                op += (' + ' + self.children[i].__str__())
            elif self.parameters[0][i] == -1:
                op += (' - ' + self.children[i].__str__())
            elif self.parameters[0][i] > 0:
                op += (' + ' + str(self.parameters[0][i]) + '*' + self.children[i].__str__())
            elif self.parameters[0][i] < 0:
                op += (' - ' + str(-self.parameters[0][i]) + '*' + self.children[i].__str__())
        if self.parameters[1] > 0:
            op += (' + ' + str(self.parameters[1]))
        elif self.parameters[1] < 0:
            op += (' - ' + str(-self.parameters[1]))
        return op + ')'

    def decompose(self):
        def addition(X):
            if len(X) == 1:
                return X[0]
            else:
                return Add([X[0], addition(X[1:])])   #SDG: use specific Add BinPredicate instead of Sum
        return [addition([(child if coef is 1 else (child * Variable(coef,coef,str(coef)))) for child, coef in zip(self.children, self.parameters[0])] + [Variable(e,e,str(e)) for e in self.parameters[1:] if e is not 0])]


class OrderedSum(Predicate):
    """
    Conjunction of a chain of precedence with a sum expression (without linear coefficients)

    The following:
        OrderedSum([a,b,c,d], l, u)

    is logically equivalent to:

        Sum([a,b,c,d]) >= l
        Sum([a,b,c,d]) <= u
        a >= b
        b >= c
        c >= d

    :param vars: the variables to be summed/sequenced.
    :param l: lower bound of the sum
    :param u: upper bound of the sum
    """

    def __init__(self, vars, l, u):
        Predicate.__init__(self, vars, "OrderedSum")

        self.parameters = [l, u]

    def close(self):
        Predicate.close(self)

    def __str__(self):
        op = str(self.parameters[0]) + ' <= ('+(self.children[0].__str__())
        for i in range(1, len(self.children)):
            op += (' + ' + self.children[i].__str__())
        return op + ') <= ' + str(self.parameters[1])


class Product(Predicate):
    """
    Syntactic sugar for a product expression on a list of :class:`.Expression`.
    `Product([x, y, z])` is equivalent to writing `(x * y * z)` but is more
    flexible in that it excepts a variable number of expressions.

    This expression is not supported by any solver at the moment, so it will be
    decomposed to chained product expressions.

    :param vars: the variables or expressions which will be multiplied by each other.
        This should be a :class:`.VarArray` or `list` with at least one item.

    .. note::

        Cannot be used as a top-level constraint, only as a sub-expression.

    .. versionadded:: 1.2.0
    """

    def __init__(self, vars):
        Predicate.__init__(self, vars, "Product")
        self.lb = None
        self.ub = None

    def __str__(self):
        return "(%s)" % " * ".join(map(str, self.children))

    def decompose(self):
        if len(self.children) == 0:
            return []
        elif len(self.children) == 1:
            return [self.children[0]]
        elif len(self.children) == 2:
            return [self.children[0] * self.children[1]]
        else:
            from operator import mul
            ret = [reduce(mul, self.children[1:], self.children[0])]
            return ret


class AllDiff(Predicate):
    """
    All-different constraint on a list of :class:`.Expression`, enforces that
    each takes a different value.

    :param vars: the variables or expressions which must take different values.
        This should be a :class:`.VarArray` or `list` with at least two items.

    .. note::

        Can only be used as a top-level constraint, not reified.
    """

    def __init__(self, vars, type=None):
        Predicate.__init__(self, vars, "AllDiff")
        if len(vars) < 2:
            raise InvalidConstraintSpecification("AllDiff requires a list of at least 2 expressions.")
        if type != None:
            self.parameters = [type]
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    #def __str__(self):
    #    return " AllDiff(" + " ".join(map(str, self.children)) + " ) "


class AllDiffExcept0(Predicate):
    """
    All-different except zero constraint on a list of :class:`.Expression`,
    enforces that each takes a different value, except those which take the
    value 0.

    :param vs: the variables or expressions which must take different values.
        This should be a :class:`.VarArray` or `list` with at least two items.

    .. note::

        Can only be used as a top-level constraint, not reified.
    """

    def __init__(self, vs):
        Predicate.__init__(self, vs, "AllDiffExcept0")
        if len(vs) < 2:
            raise InvalidConstraintSpecification("AllDiff requires a list of at least 2 expressions.")

    def decompose(self):
        from itertools import combinations
        return [Disjunction([x == 0, y == 0, x != y]) for x, y in combinations(self.children, 2)]


class Gcc(Predicate):
    """
    The Global Cardinality Constraint limits the number of times that certain
    values can be used within a set of variables. For example, we might want the
    value 1 to occur at least once and at most twice, the value 2 to occur
    exactly twice, the value 3 at most four times, and so on.

    :param vars: the variables which are being constrained. This should be a
        :class:`.VarArray` or `list` with at least two items.
    :param dict cards: A dictionary mapping each constrained value to a two item
        tuple for the lower and upper bounds on the number of occurrences of
        that value.

    .. code-block:: python

        X = VarArray(5, 1, 4)
        cards = {1: (1, 2), 2: (2, 2), 3: (0, 3), 4: (1, 2)}
        model = Model(Gcc(X,cards))

    .. note::

        Can only be used as a top-level constraint, not reified.
    """

    def __init__(self, vars, cards):
        Predicate.__init__(self, vars, "Gcc")
        values = list(cards.keys())
        values.sort()
        lb = []
        ub = []
        for val in values:
            lb.append(cards[val][0])
            ub.append(cards[val][1])
        self.parameters = [values, lb, ub]
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    def __str__(self):
        save_str = Expression.__str__
        Expression.__str__ = Expression.name
        output = " Gcc(" + " ".join(map(str, self.children)) + " | "
        for v, l, u in zip(*(self.parameters)):
            output += str(v) + ' in [' + str(l) + ',' + str(u) + '] '
        Expression.__str__ = save_str
        return output + ')'

    def decompose(self):
        X = self.children
        decomp = []
        for val, l, u in zip(self.parameters[0], self.parameters[1], self.parameters[2]):
            card = Variable(l, u)
            decomp.append((card == Cardinality(X, val)))
        return decomp


class Max(Predicate):
    """
    The maximum value of a set of Variables.

    :param vars: the variables or expressions. This should be a
        :class:`.VarArray` or `list` with at least two items.

    .. note::

        Cannot be used as a top-level constraint, only as a sub-expression.
    """

    def __init__(self, vars):
        Predicate.__init__(self, vars, "Max")
        #SDG:  initial bounds
        self.lb = max(self.get_lb(i) for i in range(len(vars)))
        self.ub = max(self.get_ub(i) for i in range(len(vars)))

    def get_min(self, solver=None):
        return max([x.get_min(solver) if type(x) not in [int, long, float, str, bool] else x for x in self.children])

    def get_max(self, solver=None):
        return max([x.get_max(solver) if type(x) not in [int, long, float, str, bool] else x for x in self.children])

    def get_value(self, solver=None):
        return max([x.get_value(solver) if type(x) not in [int, long, float, str, bool] else x for x in self.children])

    def decompose(self):
        X = self.children
        M = Variable(self.get_min(), self.get_max(), 'Max')
        decomp = [M]
        decomp.extend([M >= x for x in X])
        decomp.append(Disjunction([M <= x for x in X]))
        return decomp

    #def __str__(self):
    #    return " MAX ( " + " ".join(map(str, self.children)) + " ) "


class Min(Predicate):
    """
    The minimum value of a set of Variables.

    :param vars: the variables or expressions. This should be a
        :class:`.VarArray` or `list` with at least two items.

    .. note::

        Cannot be used as a top-level constraint, only as a sub-expression.
    """
    def __init__(self, vars):
        Predicate.__init__(self, vars, "Min")
        #SDG:  initial bounds
        self.lb = min(self.get_lb(i) for i in range(len(vars)))
        self.ub = min(self.get_ub(i) for i in range(len(vars)))

    def get_min(self, solver=None):
        return min([x.get_min(solver) if type(x) not in [int, long, float, str, bool] else x for x in self.children])

    def get_max(self, solver=None):
        return min([x.get_max(solver) if type(x) not in [int, long, float, str, bool] else x for x in self.children])

    def get_value(self, solver=None):
        return min([x.get_value(solver) if type(x) not in [int, long, float, str, bool] else x for x in self.children])

    def decompose(self):
        X = self.children
        M = Variable(self.get_min(), self.get_max(), 'Min')
        decomp = [M]
        decomp.extend([M <= x for x in X])
        decomp.append(Disjunction([M >= x for x in X]))
        return decomp

    #def __str__(self):
    #    return " MIN ( " + " ".join(map(str, self.children)) + " ) "


class Element(Predicate):
    """
    Given an integer Variable \e index and a VarArray \e vars, Element is the
    Predicate holding the value of the variable at index `index` of the array
    `vars`.

    :param vars: the variables or expressions. This should be a
        :class:`.VarArray` or `list` with at least two items.

    .. code-block:: python

        vars = VarArray(5, 1, 4)
        index = Variable(0, 4)
        elt1 = Element(vars, index)
        elt2 = vars[index]

    .. note::

        Cannot be used as a top-level constraint, only as a sub-expression.
    """

    def __init__(self, vars, index):
        children = list(vars)
        children.append(index)
        Predicate.__init__(self, children, "Element")
        #SDG:  initial bounds
        self.lb = min(self.get_lb(i) for i in range(len(vars)))
        self.ub = max(self.get_ub(i) for i in range(len(vars)))


## Boolean Clause
class Clause(Predicate):
    def __init__(self, *vars):
        Predicate.__init__(self, [], "Clause")
        polarity = []
        self.children = []
        for literal in vars:
            if literal.operator == 'neg':
                polarity.append(0)
                self.children.append(literal.children[0])
            else:
                polarity.append(1)
                self.children.append(literal)
        self.parameters = [polarity]
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    def add(self, literal):
        if literal.operator == 'neg':
            self.parameters[0].append(0)
            self.children.append(literal.children[0])
        else:
            self.parameters[0].append(1)
            self.children.append(literal)

    def __str__(self):
        ret_str = "Clause("
        for i in range(len(self.children)):
            ret_str += ' '
            if self.parameters[0][i] == 1:
                ret_str += str(self.children[i])
            else:
                ret_str += ('~' + str(self.children[i]))
        return ret_str + ' )'


class LessLex(Predicate):
    """
    Less-than lexicographic ordering constraint between two lists of
    expressions.

    :param vars_1: the fist list of variables or expressions.
    :param vars_1: the second list of variables or expressions.

    .. note::

        Can only be used as a top-level constraint, not reified.
    """

    def __init__(self, vars_1, vars_2):
        children = list(vars_1)
        children.extend(vars_2)
        Predicate.__init__(self, children, "LessLex")
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    def __str__(self):
        length = len(self.children) // 2

        toprint = '[' + str(self.children[0])
        for i in range(1, length):
            toprint += (', ' + str(self.children[i]))
        toprint += '] < [' + str(self.children[length])
        for i in range(length + 1, 2 * length):
            toprint += (', ' + str(self.children[i]))
        return toprint + ']'


class LeqLex(Predicate):
    """
    Less-than-or-equal lexicographic ordering constraint between two lists of
    expressions.

    :param vars_1: the fist list of variables or expressions.
    :param vars_1: the second list of variables or expressions.

    .. note::

        Can only be used as a top-level constraint, not reified.
    """

    def __init__(self, vars_1, vars_2):
        children = list(vars_1)
        children.extend(vars_2)
        Predicate.__init__(self, children, "LeqLex")
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    def __str__(self):
        length = len(self.children) // 2

        toprint = '[' + str(self.children[0])
        for i in range(1, length):
            toprint += (', ' + str(self.children[i]))
        toprint += '] <= [' + str(self.children[length])
        for i in range(length + 1, 2 * length):
            toprint += (', ' + str(self.children[i]))
        return toprint + ']'


class Maximise(Predicate):
    """
    Maximisation objective function, sets the goal of search to be the
    maximisation of its arguments.

    :param vars: The :class:`.Variable` or :class:`.Expression` to be maximized.
    """

    def __init__(self, vars):
        Predicate.__init__(self, [vars], "Maximise")
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    #def __str__(self):
    #    return " Maximise ( " + " ".join(map(str, self.children)) + " ) "


def Maximize(var):
    """
    Alias for American spelling of :class:`.Maximise`.
    """
    return Maximise(var)


class Minimise(Predicate):
    """
    Minimisation objective function, sets the goal of search to be the
    minimisation of its arguments.

    :param vars: The :class:`.Variable` or :class:`.Expression` to be minimized.
    """

    def __init__(self, vars):
        Predicate.__init__(self, [vars], "Minimise")
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    #def __str__(self):
    #    return " Minimise ( " + " ".join(map(str, self.children)) + " ) "


def Minimize(var):
    """
    Alias for American spelling of :class:`.Minimise`.
    """
    return Minimise(var)


class Disjunction(Predicate):
    """
    Disjunction specifying that at least one of the sub-expressions should be
    true. This can be a top-level constraint or reified as a sub-expression.

    :param vars: the variables or expressions. This should be a
        :class:`.VarArray` or `list` with at least two items.
    """

    def __init__(self, vars):
        Predicate.__init__(self, vars, "OR")
        self.lb = max(self.get_lb(i) for i in range(len(vars)))    #SDG: initialize lb,ub
        self.ub = max(self.get_ub(i) for i in range(len(vars)))

    def decompose(self):
        return [Sum(self.children) > 0]


class Conjunction(Predicate):
    """
    Conjunction specifying that all the sub-expressions should be true. This
    should only be used as a reified sub-expression, otherwise, there is an
    implicity conjunction across all top-level constraints anyway.

    :param vars: the variables or expressions. This should be a
        :class:`.VarArray` or `list` with at least two items.
    """

    def __init__(self, vars):
        Predicate.__init__(self, vars, "AND")
        self.lb = min(self.get_lb(i) for i in range(len(vars)))    #SDG: initialize lb,ub
        self.ub = min(self.get_ub(i) for i in range(len(vars)))

    def decompose(self):
        return [Sum(self.children) == len(self.children)]


# BH (2014/10/15): disabled as it appears to be buggy and not used anywhere.
# class Convex(Predicate):
#     def __init__(self, vars):
#         Predicate.__init__(self, [var for var in vars], "Convex")
#         self.lb = None  #SDG: initial lb/ub undefined
#         self.ub = None

#     def __str__(self):
#         return "[" + " ".join(map(str, self.children)) + "] is row-convex"

#     def decompose(self):
#         ### BUGGY!!

#         print "Decomposing Row convexity constraint", self
#         X = self.children
#         n = len(X)
#         first = Variable(n)
#         last = Variable(n)
#         decomposition = [X[i] <= (first <= i) for i in range(n)]
#         decomposition.extend([X[i] <= (last >= i) for i in range(n)])
#         decomposition.extend([((first <= i) & (i <= last)) <= X[i] for i in range(n)])
#         decomposition.append(first <= last)

#         print VarArray(decomposition)
#         return decomposition


class Cardinality(Predicate):
    """
    Counts the number of expressions which have been assigned a specific value.

    :param vars: the variables or expressions. This should be a
        :class:`.VarArray` or `list` with at least two items.
    :param int value: the value for which the cardinality of is being counted.

    .. note::

        Cannot be used as a top-level constraint, only as a sub-expression.
    """
    def __init__(self, vars, value):
        Predicate.__init__(self, [var for var in vars], "Card")
        self.parameters = [value]
        #SDG: initial lb/ub
        self.lb = sum(((value == x.get_lb()) and (value == x.get_ub())) for x in vars)
        self.ub = sum(((value >= x.get_lb()) and (value <= x.get_ub())) for x in vars)

    def __str__(self):
        return "card of " + str(self.parameters[0]) + " in [" + " ".join(map(str, self.children)) + "]"

    def decompose(self):
        X = self.children
        val = self.parameters[0]
        return [Sum([x == val for x in X])]


# BH (2014/10/15): disabled as it does not appear to be used anywhere
# class Cmp2(BinPredicate):
#     ## expression equal to -1 if x_0 < x_1, 1 if x_1 < x_0 and 0 otherwise
#     def __init__(self, vars):
#         BinPredicate.__init__(self, vars, "Cmp")

#     def get_symbol(self):
#         return 'cmp'

#     def decompose(self):
#         X = self.children
#         return [(X[1] < X[0]) - (X[0] < X[1])]


## @defgroup sched_group Scheduling constructs
# The scheduling constructs and constraints
#  @{
#

## A special class for simple representations of scheduling tasks.
#
#    The Task class allows for simplified modeling of tasks in scheduling
#    applications. It encapsulates the earliest start time, latest end time
#    (makespan), and duration.
#
#    There are various ways of declaring a Task:
#
#    - M = Task() creates a Task with an earliest start time of 0, latest end time of 1, and duration 1
#    - M = Task(ub) creates a Task with an earliest start time of 0, latest end time of 'ub', and duration 1
#    - M = Task(ub, dur) creates a Task with an earliest start time of 0, latest end time of 'ub', and duration 'dur'
#    - M = Task(lb, ub, dur) creates a Task with an earliest start time of 'lb', latest end time of 'ub', and duration 'dur'
#
#    When the model is solved, @get_value() returns the start
#    time of the task.
#
class Task(Expression):
    """
    The Task class allows for simplified modeling of tasks in scheduling
    applications. It encapsulates the earliest start time, latest end time
    (makespan), and duration. The following tables shows example calls to the constructor
    which results in different kinds of matrices.

    .. code-block:: python

        Task()  # creates a Task with an earliest start time of 0, latest end
                # time of 1, and duration 1
        Task(ub)  # creates a Task with an earliest start time of 0, latest end
                  # time of 'ub', and duration 1
        Task(ub, dur)  # creates a Task with an earliest start time of 0,
                       # latest end time of 'ub', and duration 'dur'
        Task(lb, ub, dur)  # creates a Task with an earliest start time of 'lb',
                           # latest end time of 'ub', and duration 'dur'

    When the model is solved, :func:`Numberjack.Expression.get_value` returns
    the start time of the task.
    """

    def __init__(self, arg1=None, arg2=None, arg3=None):

        lb = 0
        ub = 1
        self.duration = 1
        if arg1 != None:
            if arg2 != None:
                if arg3 != None:
                    lb = arg1
                    ub = arg2 - arg3
                    self.duration = arg3
                else:  # only 2 args, read as <makespan,duration>
                    ub = arg1 - arg2
                    self.duration = arg2
            else:  # only 1 arg, read as the makespan
                ub = arg2 - 1
        Expression.__init__(self, "t")
        self.lb = lb
        self.ub = ub
        self.domain_ = None

    def __str__(self):
        if self.is_built() and self.solver.is_sat():
            return str(self.get_value())
        else:
            ident = str(self.ident)
            if self.ident == -1:
                ident = ''
            return 't' + str(ident) + ': [' + str(self.get_min()) + ':' + str(self.duration) + ':' + str(self.get_max() + self.duration) + ']'

    def __lt__(self, pred):
        """
        Creates a precedence expression on the task. Can be used to specify
        precedence between two tasks or that a task has finished before a
        certain time.

        :param pred: if `pred` is an `int` then constrains that the task is
            finished by `pred`, if `pred` is another :class:`.Task` instance,
            then creates a top-level :class:`.Precedence` constraint.
        """
        if type(pred) is int:
            return Le([self, pred - self.duration])
        return Precedence(self, pred, self.duration)

    def __gt__(self, pred):
        if type(pred) is int:
            return Gt([self, pred])
        return pred.__lt__(self)

    def requires(self, resource):
        resource.add(self)

    def reset(self, makespan):
        self.ub = makespan - self.duration


class Precedence(Predicate):
    """
    Precedence constraint enforces a certain gap between the start time of two
    tasks. Equivalent to `task_i + duration <= task_j`.

    :param Task task_i: the first task.
    :param Task task_j: the second task.
    :param int dur: the gap to maintain between the two tasks.
    """
    def __init__(self, task_i, task_j, dur):
        Predicate.__init__(self, [task_i, task_j], "Precedence")
        self.parameters = [dur]
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    def decompose(self):
        return [((self.children[0] + self.parameters[0]) <= self.children[1])]

    def __str__(self):
        return str(self.children[0]) + ' + ' + str(self.parameters[0]) + ' <= ' + str(self.children[1])


class NoOverlap(Predicate):
    """
    Binary disjunctive constraint enforces that two tasks do not overlap.
    Equivalent to `(task_i + duration_i <= task_j) | (task_j + duration_j <=
    task_i)` .

    :param Task task_i: the first task.
    :param Task task_j: the second task.
    :param int dur_i: the duration of the first task, if `None` then
        `task_i.duration` will be used.
    :param int dur_j: the duration of the second task, if `None` then
        `task_j.duration` will be used.
    """

    def __init__(self, task_i, task_j, dur_i=None, dur_j=None):
        if dur_i == None:
            dur_i = task_i.duration
        if dur_j == None:
            dur_j = task_j.duration
        Predicate.__init__(self, [task_i, task_j], "NoOverlap")
        self.parameters = [dur_i, dur_j]
        self.lb = None  #SDG: initial lb/ub undefined
        self.ub = None

    def decompose(self):
        return [(((self.children[0] + self.parameters[0]) <= (self.children[1])) |
                 ((self.children[1] + self.parameters[1]) <= (self.children[0])))]

    def __str__(self):
        return str(self.children[0]) + ' + ' + str(self.parameters[0]) + ' <= ' + str(self.children[1]) + ' OR ' + str(self.children[1]) + ' + ' + str(self.parameters[1]) + ' <= ' + str(self.children[0])


class UnaryResource(list):
    """
    Unary resource constraint ensures that only one of the specified list of
    tasks are running at each time point. An optional distance between tasks can
    be specified also.

    :param arg: a list of :class:`.Task` instances.
    :param int distance: optional distance between tasks.
    """

    def __init__(self, arg=[], distance=0):
        list.__init__(self, [NoOverlap(arg[i], arg[j], arg[i].duration + distance, arg[j].duration + distance) for i in range(1, len(arg)) for j in range(i)])
        self.tasks = [task for task in arg]
        self.distance = distance

    def add(self, new_task):
        """
        Add an additional task to the existing list of tasks obeying this unary
        resource.

        :param Task new_task: the additional task to include.
        """
        for task in self.tasks:
            self.append(NoOverlap(new_task, task, new_task.duration + self.distance, task.duration + self.distance))
        self.tasks.append(new_task)

    def __str__(self):
        return "[" + " ".join(map(str, self.tasks)) + "] share a unary resource"  # "+" ".join(map(str, self))


# BH (2014/10/15): disabled as it does not appear to be used anywhere and
# appears to duplicate UnaryResource
# class UnaryResourceB(Predicate):
#     def __init__(self, tasks, distance=0):
#         Predicate.__init__(self, [task for task in tasks], "UnaryResource")
#         self.distance = distance

#     def add(self, new_task):
#         self.children.append(new_task)

#     def __str__(self):
#         return "[" + " ".join(map(str, self.children)) + "] share a unary resource"

#     def decompose(self):
#         return [NoOverlap(task1, task2)  # , task1.duration+self.distance, task2.duration+self.distance)
#                 for task1, task2 in pair_of(self.children)]


class ParamList(dict):
    def __init__(self, X):
        dict.__init__(self, X)

    def __call__(self, *args):
        return [self.__getitem__(arg) for arg in args]


def input(default):
    """
    Read command line arguments from the user. This is useful to establish some
    default parameters of your model and solving process, and to subsequently
    allow these to be easily changed by specifying a command line argument.

    The `default` argument allows you to specify the list of allowed options, as
    well as their default values. Any option given on the command line that is
    not present in this list will raise an error. Values specified on the
    command line will be coerced into the same data type as is given for that
    option in `default`.

    :param dict default: a dictionary of the valid options and their default
        values.
    :return: a dictionary of the options from `default`, possibly with updated
        values from the command line.
    :rtype: dict

    For example, if launching I can change the default parameters of a model
    like so:

    .. code-block:: bash

        python numberjackfile.py -solver MiniSat -N 10

    .. code-block:: python

        # numberjackfile.py
        default = {'N': 5, 'solver': 'Mistral', 'tcutoff': 30}
        param = input(default)
        # param will be a dict {'N': 10, 'solver': 'MiniSat', 'tcutoff': 30}

    .. deprecated:: 1.1

        This function will be renamed or replaced in 2.0 to avoid the naming
        clash with the builtin input function when imported with `*`.
    """
    import sys

    param_list = ParamList(default)
    option = None
    params = []

    commandline = [arg for arg in sys.argv[1:]]
    commandline.append('-end_argument')

    for arg in commandline:
        #print arg

        if arg[0] == '-' and arg != '-1.0':  # new argument
            # first take previous param into account
            #print 'end of option:', params
            if option != None and option != '1.0':
                if len(params) == 0:
                    param_list[option] = 'yes'
                elif len(params) == 1:
                    if type(param_list[option]) == int:
                        #print 'int'
                        param_list[option] = int(params[0])
                    elif type(param_list[option]) == float:
                        #print 'float'
                        param_list[option] = float(params[0])
                    else:
                        #print 'string'
                        param_list[option] = params[0]
                else:
                    if len(param_list[option]) > 0:
                        if type(param_list[option][0]) == int:
                        #print 'int'
                            np = [int(p) for p in params]
                            params = np
                        elif type(param_list[option][0]) == float:
                        #print 'float'
                            np = [float(p) for p in params]
                            params = np
                    param_list[option] = params

            option = arg[1:]
            #print 'new option', option

            if option != 'end_argument':
                if option not in param_list:

                    #print 'unknwn option'
                    if option == 'h' or option == '-h' or option == 'help' or option == '-help':
                        #print 'help'
                        the_keys = list(param_list.keys())
                        the_keys.sort()
                        for key in the_keys:
                            print(('-' + key).ljust(20) + ":", end=' ')
                            if type(param_list[key]) == int:
                                print('int'.ljust(14), end=' ')
                            if type(param_list[key]) == float:
                                print('float'.ljust(14), end=' ')
                            elif type(param_list[key]) == str:
                                print('string'.ljust(14), end=' ')
                            elif hasattr(param_list[key], '__iter__'):
                                print('list', end=' ')
                                if len(param_list[key]) > 0:
                                    if type(param_list[key][0]) == int:
                                        print('of int   ', end=' ')
                                    elif type(param_list[key][0]) == float:
                                        print('of float ', end=' ')
                                    else:
                                        print('of string', end=' ')
                                else:
                                    print('of string', end=' ')
                            print(' (default=' + str(param_list[key]) + ')')
                        exit(1)
                    else:
                        print('Warning: wrong parameter name, ignoring arguments following', option)
                else:
                    #print 'fine option init param list'
                    params = []
        else:
            params.append(arg)
            #print 'input param:', params

    return param_list


def pair_of(l):
    #return [pair for k in range(1, len(l)) for pair in zip(l, l[k:])]
    pairs = []
    for j in range(1,len(l)):
        for i in range(j):
            pairs.append((l[i],l[j]))
    return pairs



def value(x):
    return x.get_value()


def total_seconds(td):
    # Python 2.6 doesn't have timedelta.total_seconds
    f = getattr(td, 'total_seconds', None)
    if f:
        return f()
    else:
        return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 1e6) / 1e6


def load_in_decompositions():
    from . import Decomp

    # First add the constraints
    for attr_name in dir(Decomp):
        attr = getattr(Decomp, attr_name)
        if hasattr(attr, "__name__") and attr.__name__.find("decompose") is -1:
            if not hasattr(sys.modules[__name__], attr.__name__):
                setattr(sys.modules[__name__], attr.__name__, attr)

    # Now load up the decompositions, introspection is great huh :)
    for attr_name in dir(sys.modules[__name__]):
        if hasattr(Decomp, "decompose_" + attr_name):
            setattr(getattr(sys.modules[__name__], attr_name),
                    "decompose", getattr(Decomp, "decompose_" + attr_name))


load_in_decompositions()


class Solution(list):
    """
    Class which will extract a list of the solution values for a given list of
    variables. Solution values are order by the order in which the variables are
    given in `vars`. This functionality is equivalent to calling
    :func:`Expression.get_value` on each variable, but it is just wrapped here.

    :param vars: a Matrix, VarArray, or list of variables to extract the
        solution for.
    """
    def __init__(self, vars):
        list.__init__(self)
        if issubclass(type(vars), Matrix):
            self.variables = vars.flat
            for row in vars.row:
                self.append([x.get_value() for x in row])
        else:
            self.variables = vars
            self.extend([x.get_value() for x in vars])
        self.dico = {}.fromkeys(self.variables)
        for x in self.variables:
            self.dico[x] = x.get_value()

    def __getitem__(self, i):
        if type(i) is int:
            return list.__getitem__(self, i)
        else:
            return self.dico[i]

    def __contains__(self, x):
        return x in self.dico

    def __str__(self):
        if len(self) == 0:
            return '[]'
        elif type(self[0]) is list:
            return '[' + ',\n '.join([str(row) for row in self]) + ']'
        else:
            return list.__str__(self)


class Nogood(object):
    def __init__(self, clause, solver):
        self.clause = clause
        self.solver = solver

    def __str__(self):
        return self.solver.print_clause(self.clause)


keep_alive = []
class NBJ_STD_Solver(object):
    """
        Generic solver class which will be subclassed by the solver interfaces.
        Provides common functionality which will be used in each of the solver
        interfaces.

        .. note::

            The user should not need to instantiate this class, instead
            :func:`Numberjack.Model.load` should be used to return an instance
            of a subclass. However the following methods can be used on that
            instance.
    """

    def __init__(self, Library, Wrapper, model=None, X=None, FD=False,
                 clause_limit=-1, encoding=None):
        self.decomposition_store = []
        self.enc_config_cache = {}
        self.free_memory = None
        self.verbosity = 0

        # Time to load the model into the solver. Can be used to track the time
        # to encode time for linearization and SAT. Only if the model has been
        # passed to the constructor of the solver though or model.load() is
        # used.
        self.load_time = None

        # self.solver = getattr(sys.modules[Library], Library + "Solver", None)()
        solverpkg = "Numberjack.solvers"
        libstr = "%s.%s" % (solverpkg, Library)
        wrapstr = "%s.%s" % (solverpkg, Wrapper)
        self.solver = getattr(sys.modules[libstr], Library + "Solver", None)()

        if hasattr(self.solver, "setClauseLimit"):
            self.solver.setClauseLimit(clause_limit)

        self.LibraryPacakge = libstr
        self.WrapperPackage = wrapstr
        self.Library = Library
        self.Wrapper = Wrapper

        self.ExpArray = getattr(sys.modules[wrapstr],
                                Wrapper + "ExpArray", None)

        self.IntArray = getattr(sys.modules[wrapstr],
                                Wrapper + "IntArray", None)

        self.DoubleArray = getattr(sys.modules[wrapstr],
                                   Wrapper + "DoubleArray", None)

        self.IntVar = getattr(sys.modules[wrapstr],
                              Wrapper + "_IntVar", None)

        self.FloatVar = getattr(sys.modules[wrapstr],
                                Wrapper + "_FloatVar", None)

        self.EncodingConfiguration = getattr(sys.modules[wrapstr],
                                             "EncodingConfiguration", None)

        if "_" + Library in sys.modules:
            self.free_memory = getattr(sys.modules["_" + Library],
                                       "delete_" + Library + "Solver", None)

        self.variables = None

        if model is not None:
            var_array = None
            self.solver_id = model.getSolverId()
            self.model = weakref.proxy(model)
            self.model.close(self)   #SDG: needs to know for which solver the model is built
            if self.EncodingConfiguration:
                if not encoding:
                    encoding = EncodingConfiguration()
                self.solver.encoding = self.getEncodingConfiguration(encoding)

            # Load each expression and its variables into the solver. Record
            # time to do so.
            loadstart = datetime.datetime.now()
            for expr in self.model.get_exprs():
                self.solver.add(self.load_expr(expr))
            self.load_time = total_seconds(datetime.datetime.now() - loadstart)

            if X is not None:
                self.variables = VarArray(flatten(X))
                var_array = self.ExpArray()
                for x in self.variables:
                    if len(x.var_list) > self.solver_id - 1:
                        var_array.add(x.var_list[self.solver_id - 1])
                if FD:
                    self.solver.forceFiniteDomain(var_array)
                if var_array.size() > 0:
                    self.solver.initialise(var_array)
            else:
                self.variables = weakref.proxy(self.model.variables)
                self.solver.initialise()

    def add_to_store(self, x):
        if issubclass(type(x), Predicate):
            self.decomposition_store.append(x)
            for child in x.children:
                self.add_to_store(child)

    #def getIntVar(self, lb, ub, ident):
    def getIntVar(self, arg1, arg2, argopt1=None):
        var = None
        try:
            if argopt1 is None:
                var = self.IntVar(arg1, arg2)
            else:
                var = self.IntVar(arg1, arg2, argopt1)
        except:
            raise Exception("ERROR while creating variable")
        return var

    def getFloatVar(self, lb, ub, ident):
        var = None
        try:
            var = self.FloatVar(lb, ub, ident)
        except:
            raise ValueError("ERROR: Solver does not support real variables")
        return var

    def getEncodingConfiguration(self, enc_config):
        # Returns the C++ EncodingConfiguration equivalent of enc_config from
        # the wrapper.
        if not self.EncodingConfiguration:
            raise UnsupportedSolverFunction(self.Library, "EncodingConfiguration", "This solver does not support custom encoding settings.")

        if enc_config not in self.enc_config_cache:
            try:
                self.enc_config_cache[enc_config] = self.EncodingConfiguration(
                    enc_config.direct, enc_config.order,
                    enc_config.conflict, enc_config.support,
                    enc_config.amo_encoding, enc_config.alldiff_encoding)
            except Exception as e:
                raise e
        return self.enc_config_cache[enc_config]

    def load_expr(self, expr):
        # print('load', expr , expr.get_lb() if issubclass(type(expr),Expression) else None, expr.get_ub() if issubclass(type(expr),Expression) else None)   #SDG: VERY USEFUL FOR DEBUGGING
        if type(expr) is str:
            return self.model.string_map[expr]
        if type(expr) is bool:
            return int(expr)
        if type(expr) in [int, int, float]:
            # It is a constant, handle appropriatly
            return expr
        #if not expr.has_children():
        if expr.is_var():
            # It is a leaf
            if expr.is_var():  # Just to be sure
                alreadyBuild = False
                if expr.is_built(self):
                    #if expr.get_solver() == self: # Checks if I have already assigned
                    alreadyBuild = True

                if alreadyBuild:
                    return expr.var_list[self.solver_id - 1]
                else:
                    # It is probably a variable
                    (lb, ub, domain) = expr.get_domain_tuple()

                    var = None
                    if domain is None or (ub - lb + 1) == len(domain):
                        if type(lb) is int:
                            var = self.getIntVar(lb, ub, expr.ident)
                        else:
                            var = self.getFloatVar(lb, ub, expr.ident)
                    else:
                        w_array = self.IntArray()
                        for val in domain:
                            w_array.add(val)
                        var = self.getIntVar(w_array, expr.ident)

                    if expr.encoding:
                        var.encoding = self.getEncodingConfiguration(expr.encoding)
                    expr.setVar(self.solver_id, self.Library, var, self)
                    expr.solver = self
                    return var
            else:
                raise Exception("Problem, no such type exists in converting models")
        else:
#            factory = getattr(sys.modules[self.Library],
#                              "%s_%s" % (self.Library, expr.get_operator()), None)
            factory = getattr(sys.modules[self.WrapperPackage],
                              "%s_%s" % (self.Wrapper, expr.get_operator()), None)

            if factory is not None:
                arguments = None
                if len(expr.get_children()) <= 2:
                    arguments = [self.load_expr(child) for child in expr.get_children()]
                else:
                    var_array = self.ExpArray()
                    for child in expr.get_children():
                        myarg = self.load_expr(child)
                        var_array.add(myarg)
                        keep_alive.append(myarg)
                    arguments = [var_array]
                if expr.has_parameters():  # != None: # assumes an array of integers
                    for param in expr.parameters:
                        if hasattr(param, '__iter__') and not isinstance(param, str): # in python3, strings *do* have an __iter__ method
                            w_array = None
                            if any((type(w) == float for w in param)):
                                w_array = self.DoubleArray()
                            else:
                                w_array = self.IntArray()

                            if len(param) and hasattr(param[0], '__iter__'): # expr.get_operator() is "Table": #SdG: also need for PostNary
                                for w in param:
                                    for v in w:
                                        w_array.add(v)
                                arguments.append(w_array)
                            else:
                                for w in param:
                                    w_array.add(w)
                                arguments.append(w_array)
                        else:
                            arguments.append(param)
                try:
                    var = factory(*arguments)
                except NotImplementedError as e:
                    print("Error the solver does not support this expression:", str(expr), file=sys.stderr)
                    print("Type:", type(expr), file=sys.stderr)
                    print("Children:", str(expr.children), map(type, expr.children), file=sys.stderr)
                    print("Params:", str(getattr(expr, 'parameters', None)), file=sys.stderr)
                    raise e

                if expr.encoding:
                    var.encoding = self.getEncodingConfiguration(expr.encoding)
                expr.setVar(self.solver_id, self.Library, var, self)
                expr.solver = self
                return var
            else:
                return self.decompose_expression(expr)

    def decompose_expression(self, expr):
        if hasattr(expr, "decompose"):
            expr_list = expr.decompose()
            #print expr_list   #SDG: VERY USEFUL FOR DEBUGGING
            obj_exp = []
            #SDG: all decomposed expressions except the first are assumed to be constraints (ie, at the top-level)
            for exp in expr_list[1:]:
                exp.encoding = expr.encoding
                obj = self.load_expr(exp)
                obj_exp.append(obj)
                self.solver.add(obj)

            expr.solver = self
            #expr_list[0].close()
            for exp in expr_list:
                exp.close()
                self.add_to_store(exp)
            decomp = self.load_expr(expr_list[0])
            return decomp
        else:
            raise ConstraintNotSupportedError(expr.get_operator(), self.Library)

    def solve(self):
        """
        Calls solve on the underlying solver.

        Captures :exc:`KeyboardInterrupt` or :exc:`SystemExit` signals and
        returns None in this case.

        :return: `True` if the solver found a satisfiable solution, `False`
            otherwise.
        """
        try:
            if self.solver.solve() == SAT:
                return True
            return False
        except (KeyboardInterrupt, SystemExit):
            print('Program Interrupted')
            return

    def solveAndRestart(self, policy=GEOMETRIC, base=64, factor=1.3, decay=0.0, reinit=-1):
        """
        Calls solve with restarts on the underlying solver.

        .. deprecated:: 1.1
           Restarting is the default in most underlying solvers, currently only
           Mistral version 1.55 will not perform restarts unless started with
           this method. Instead you should use :func:`solve` which will use
           restarts by default in other solvers.
        """
        if reinit == -1:
            if self.solver.solveAndRestart(policy, base, factor, decay) == SAT:
                return True
        else:
            if self.solver.solveAndRestart(policy, base, factor, decay, reinit) == SAT:
                return True
        return False

    def startNewSearch(self):
        "Initialise structures for a depth first search."
        self.solver.startNewSearch()

    def getNextSolution(self):
        "Search for the next solution"
        return self.solver.getNextSolution()

    def solutions(self):
        "A generator which will yield `True` until no other solution exists."
        while self.getNextSolution():
            yield True  # Could return something more useful??

    def next(self, exp, v):
        # \internal used to get the next value in the domain of 'exp' after v
        #return self.solver.next(exp.var_list[self.solver_id-1], v)
        return exp.var_list[self.solver_id - 1].next(v)

    ##@name Search programming methods
    # @{

    def propagate(self):
        """
        Tell the solver to reach a fixed point.

        :return: `False` if an inconsistency is found, `True` otherwise.
        """
        return self.solver.propagate()

    def save(self):
        "Tell the solver to save the current state."
        self.solver.save()

    def undo(self, bj=1):
        """
        Tell the solver to restore its state to the one last enqueued by
        'save()'
        """
        return self.solver.undo(bj)

    def post(self, exp):
        """
        Tell the solver to add a constraint in the current state.

        :param exp: should be a unary constraint, for example `x == 5`.
        """
        self.solver.post(exp.operator, exp.children[0].var_list[self.solver_id - 1], exp.children[1])

    def deduce(self, exp=None):
        """
        Tell the solver to post the negation of the last decision in the current
        states.
        """
        if exp is None:
            self.solver.deduce()
        else:
            self.solver.post(exp.operator, exp.children[0].var_list[self.solver_id - 1], exp.children[1])

    def deduce_print(self, lvl):
        x = self.variables[self.solver.get_decision_id()]
        print(lvl * ' ', x.domain(self))
        self.solver.deduce()

    ## Follow a 'left' branch in a binary serach tree (eq. to save() + post())
    def branch_left(self, exp):
        self.solver.save()
        self.solver.post(exp.operator, exp.children[0].var_list[self.solver_id - 1], exp.children[1])

    ## Follow a 'right' branch in a binary serach tree (eq. to undo() + deduce())
    def branch_right(self):
        return self.solver.branch_right()

    def reset(self, full=False):
        """
        Resets the data structure of the solver to the initial state.

        :param bool full: whether the top-level deduction should be undone too.
        """
        self.solver.reset(full)

    ## @}

    def analyze_conflict(self):
        btlevel = self.solver.analyze_conflict()
        if btlevel < 0:
            btlevel = None
        return (self.solver.get_learnt_clause(), btlevel)

    def print_all_clauses(self):
        for i in range(self.solver.nbClauses()):
            self.solver.get_clause(i)
            self.print_clause()

    def get_last_nogood(self):
        #self.solver.get_last_nogood(i)
        self.print_clause()

    def print_clause(self):
        def get_literal(i):
            var = self.solver.get_nogood_var(i)
            val = self.solver.get_nogood_val(i)
            type = self.solver.get_nogood_type(i)
            sign = self.solver.get_nogood_sign(i)
            lit = str(self.variables[var])

            if type == 0:
                if sign == 0:
                    lit += ' == '
                else:
                    lit += ' != '
            else:
                if sign == 0:
                    lit += ' <= '
                else:
                    lit += ' > '
            lit += str(val)
            return lit

        print('(', end=' ')
        for i in range(self.solver.get_nogood_size()):
            if i > 0:
                print('or', end=' ')
            print(get_literal(i), end=' ')
        print(')')

    def sacPreprocess(self, type):
        self.solver.sacPreprocess(type)
        return (not (self.solver.is_unsat()))

    ##@name Parameter tuning methods
    # @{

    def setHeuristic(self, arg1, arg2='No', arg3=0):
        """
        Sets the variable and value heuristics.

        .. note::
            Currently only supports setting the heuristics for Mistral solvers
            but a generic method for all solvers is in the works.
        """
        var_name = arg1
        val_name = arg2
        randomization = arg3
        if type(val_name) is int:
            randomization = val_name
            val_name = 'No'
        if var_name not in var_heuristics:
            print('c Warning: "' + var_name + '" unknown, use MinDomain instead')
            print('c legal variable orderings: ', var_heuristics)
        if val_name not in val_heuristics:
            print('c Warning: "' + val_name + '" unknown, use Lex instead')
            print('c legal value orderings: ', val_heuristics)
        self.solver.setHeuristic(str(var_name), str(val_name), randomization)

    def setFailureLimit(self, cutoff):
        """
        Sets a limit on the number of failures encountered before aborting
        search.
        """
        self.solver.setFailureLimit(cutoff)

    def setTimeLimit(self, cutoff):
        """
        Sets a limit on the CPU time before aborting search.
        """
        self.solver.setTimeLimit(cutoff)

    def setNodeLimit(self, cutoff):
        """
        Sets a limit on the number of nodes explored before aborting search.
        """
        self.solver.setNodeLimit(cutoff)

    def setVerbosity(self, degree):
        """
        Sets the verbosity level of the solver.
        """
        self.verbosity = max(0, degree)
        self.solver.setVerbosity(degree)

    def setThreadCount(self, num_threads):
        """
        Sets the number of threads a solver should use.
        """
        f = getattr(self.solver, 'setThreadCount', None)
        if f:
            f(num_threads)
        else:
            if self.verbosity > 0:
                print("Warning: this solver does not support " \
                    "the ability to specify a thread count.", file=sys.stderr)

    def setOptimalityGap(self, gap):
        """
        Sets the target relative optimality gap tolerance.
        """
        if hasattr(self.solver, 'setOptimalityGap'):
            return self.solver.setOptimalityGap(gap)
        else:
            raise UnsupportedSolverFunction(
                self.Library, "setOptimalityGap", "This solver does not support"
                " setting the optimility gap.")
        return None

    def setRandomSeed(self, seed):
        """
        Sets the initial random seed.
        """
        self.solver.setRandomSeed(seed)

    def setWorkMem(self, mb):
        """
        Set the limit of working memory, only used for CPLEX.
        """
        if hasattr(self.solver, 'setWorkMem'):
            return self.solver.setWorkMem(mb)
        else:
            raise UnsupportedSolverFunction(
                self.Library, "setWorkMem", "This solver does not support"
                " setting the work memory.")
        return None

    def getWorkMem(self):
        """
        Get the limit of working memory, only used for CPLEX.
        """
        if hasattr(self.solver, 'getWorkMem'):
            return self.solver.getWorkMem()
        else:
            raise UnsupportedSolverFunction(
                self.Library, "getWorkMem", "This solver does not support"
                " getting the work memory.")
        return None

    def setOption(self,func,param=None):
        """
        Sets an option in Toulbar2 whose name is passed as the first parameter,
        and value as a second one.
        """
        try:
            function = getattr(self.solver,func)
        except AttributeError:
            print("Warning: "+func+" option does not exist in this solver!")
        else:
            if param is None:
                function()
            else:
                function(param)

    ## @}

    def setRandomized(self, degree):
        self.solver.setRandomized(degree)

    def addNogood(self, vars, vals):
        var_array = self.ExpArray()
        for var in vars:
            var_array.add(var.var_list[self.solver_id - 1])
        val_array = self.IntArray()
        for val in vals:
            val_array.add(val)
        self.solver.addNogood(var_array, val_array)

    def setAntiLex(self, vars):
        var_array = self.ExpArray()
        for var in vars:
            var_array.add(var.var_list[self.solver_id - 1])
        self.solver.setAntiLex(var_array)

    def guide(self, vars, vals=None, probs=[]):
        var_array = self.ExpArray()
        val_array = self.IntArray()
        pro_array = self.DoubleArray()

        if vals is not None:
            for var in vars:
                var_array.add(var.var_list[self.solver_id - 1])
            for val in vals:
                val_array.add(val)
        else:
            for var in vars.variables:
                var_array.add(var.var_list[self.solver_id - 1])
                if var in vars:
                    val_array.add(vars[var])

        for pro in probs:
            pro_array.add(pro)
        self.solver.guide(var_array, val_array, pro_array)

    #def backtrackTo(self, level):
    #    self.solver.backtrackTo(level)
    #def assign(self,x,v):
    #    self.solver.assign(x.var,v)
    #def upOneLevel(self):
    #    self.solver.upOneLevel()

    def setLowerBounds(self, vars, lb):
        var_array = self.ExpArray()
        lob_array = self.IntArray()
        for (x, l) in zip(vars, lb):
            var_array.add(x.var_list[self.solver_id - 1])
            lob_array.add(l)
        self.solver.setLowerBounds(var_array, lob_array)

    def setUpperBounds(self, vars, ub):
        var_array = self.ExpArray()
        upb_array = self.IntArray()
        for (x, u) in zip(vars, ub):
            var_array.add(x.var_list[self.solver_id - 1])
            upb_array.add(u)
        self.solver.setUpperBounds(var_array, upb_array)

    def setRestartNogood(self):
        self.solver.setRestartNogood()

    ##@name Accessors
    # @{

    def get_solution(self):
        """
        Extract a :class:`.Solution` object from the solver representing the
        list of assigned values.

        :rtype: :class:`.Solution`
        """
        self.solver.store_solution()
        solution = Solution(self.variables)
        return solution

    def is_opt(self):
        """
        Returns `True` if the solver found a solution and proved its
        optimality, `False` otherwise.
        """
        return self.solver.is_opt()

    def is_sat(self):
        "Returns `True` if the solver found a solution, `False` otherwise."
        return self.solver.is_sat()

    def is_unsat(self):
        """
        Returns `True` if the solver proved unsatisfiability, `False` otherwise.
        """
        return self.solver.is_unsat()

    def getOptimum(self):
        """
        Returns the current best solution cost from Toulbar2.

        :raises UnsupportedSolverFunction: if called on a solver other than
            Toulbar2.
        """
        if hasattr(self.solver, 'getOptimum'):
            return self.solver.getOptimum()
        else:
            raise UnsupportedSolverFunction(
                self.Library, "getOptimum", "This solver does not support "
                "getOptimum, this is a Toulbar2 function only.")

    def getOptimalityGap(self):
        """
        Returns the optimality gap from the solver. Valid for MIP solvers only.

        :raises UnsupportedSolverFunction: if called on a non MIP solver.
        """
        if hasattr(self.solver, 'getOptimalityGap'):
            return self.solver.getOptimalityGap()
        else:
            raise UnsupportedSolverFunction(
                self.Library, "getOptimalityGap", "This solver does not "
                "support getting the optimility gap.")

    def getBacktracks(self):
        "Returns the number of backtracks performed during the last search."
        return self.solver.getBacktracks()

    def getNodes(self):
        "Returns the number of nodes explored during the last search."
        return self.solver.getNodes()

    def getFailures(self):
        "Returns the number of failures encountered during the last search."
        return self.solver.getFailures()

    def getPropags(self):
        """
        Returns the number of constraint propagations performed during the last
        search.
        """
        return self.solver.getPropags()

    def getTime(self):
        "Returns the CPU time required for the last search."
        return self.solver.getTime()

    ## @}

    def getChecks(self):
        "Returns the number of constraint checks. for the last search."
        return self.solver.getChecks()

    def printStatistics(self):
        """
        Asks the solver to print some basic statistics about its last search.

        .. deprecated:: 1.1
        """
        print('')
        self.solver.printStatistics()
        print('')

    def getNumVariables(self):
        """
        Get the number of variables that have been created in the underlying
        solver. This figure can be different to the number of variables that you
        created in your model. For SAT and MIP solvers, this figure will be the
        number of Boolean variables which had to be created during the encoding
        step, including any auxiliary variables.
        """
        return self.solver.getNumVariables()

    def getNumConstraints(self):
        """
        Get the number of constraints that have been created in the underlying
        solver. This figure can be different to the number of constraints that
        you created in your model. For SAT solvers it will be the number of CNF
        clauses created by the encoding, for MIP solvers it will be the number
        of linear expressions created.
        """
        return self.solver.getNumConstraints()

    def load_xml(self, file, type=4):
        """
        This function only allows you to load an XCSP instance into the Mistral
        solver. You should use the :mod:`Numberjack.XCSP` module to build a
        generic Numberjack model from an XCSP instance, which can be loaded with
        other underlying solvers.

        .. deprecated:: 1.1
            Use :mod:`Numberjack.XCSP` instead.
        """
        self.solver.load_xml(file, type)

    def load_mps(self, filename, extension):
        """
        Asks the underlying MIP solver to load an MPS file.

        :param filename: the path to the file.
        :param extension: the file's extension.
        :raises UnsupportedSolverFunction: if called on a non MIP solver.
        """
        if not hasattr(self.solver, 'load_mps'):
            raise UnsupportedSolverFunction(
                str(type(self)), "load_mps", "Please load the model using a "
                "MIP solver to use this functionality.")
        self.solver.load_mps(filename, extension)

    def load_gmpl(self, filename, data=None):
        """
        Asks the underlying MIP solver to load a GMPL file, possibly with a
        separate data file.

        :param filename: the path to the file.
        :param data: optional path to a data file.
        :raises UnsupportedSolverFunction: if called on a non MIP solver.
        """
        if not hasattr(self.solver, 'load_gmpl'):
            raise UnsupportedSolverFunction(
                str(type(self)), "load_gmpl", "Please load the model using a "
                "MIP solver to use this functionality.")
        if data == None:
            self.solver.load_gmpl(filename)
        else:
            self.solver.load_gmpl(filename, data)

    def load_lp(self, filename, epsilon):
        """
        Asks the underlying MIP solver to load an LP file, possibly with a
        separate data file.

        :param filename: the path to the file.
        :param epsilon: epsilon
        :raises UnsupportedSolverFunction: if called on a non MIP solver.
        """
        if not hasattr(self.solver, 'load_lp'):
            raise UnsupportedSolverFunction(
                str(type(self)), "load_lp", "Please load the model using a "
                "MIP solver to use this functionality.")
        self.solver.load_lp(filename, epsilon)

    def shuffle_cnf(self, *args, **kwargs):
        """
        Shuffle the internal CNF representation before writing it to a file.
        This renames the variables, shuffles their order in each clause, and
        shuffles the ordering of the clauses. This currently has no affect on
        the built-in MiniSat or WalkSat solvers since clauses are added directly
        via their API when they are generated but can be used with any of the
        other external file based SAT solvers. This should be called before
        :meth:`.output_cnf`.

        :param int seed: The seed for the random number generator.
        :raises UnsupportedSolverFunction: if called on a non SAT-based solver.
        """
        if hasattr(self.solver, 'shuffle_cnf'):
            self.solver.shuffle_cnf(*args, **kwargs)
        else:
            raise UnsupportedSolverFunction(
                str(type(self)), "shuffle_cnf", "Please load the model using "
                "a SAT solver to use this functionality.")

    def output_cnf(self, filename):
        """
        Output the CNF representation of a model to a file. The model must have
        been loaded with a SAT-based solver.

        :param str filename: The filename of where to output the CNF file.
        :raises UnsupportedSolverFunction: if called on a non SAT-based solver.
        """
        from Numberjack.solvers.SatWrapper import SatWrapperSolver as sws
        if not issubclass(type(self.solver), sws):
            raise UnsupportedSolverFunction(
                str(type(self)), "output_cnf", "Please load the model using a "
                "SAT solver to use this functionality.")
        self.solver.output_cnf(filename)

    def output_lp(self, filename):
        """
        Output the LP representation of a model to a file. The model must have
        been loaded with a MIP-based solver.

        :param str filename: The filename of where to output the LP file.
        :raises UnsupportedSolverFunction: if called on a non MIP-based solver.
        """
        if hasattr(self.solver, 'output_lp'):
            if not filename.endswith(".lp"):
                filname = "%s.lp" % filename
            self.solver.output_lp(filename)
        else:
            raise UnsupportedSolverFunction(
                str(type(self)), "output_lp", "This solver does not support "
                "outputing LP files.")

    def output_mps(self, filename):
        """
        Output the MPS representation of a model to a file. The model must have
        been loaded with a MIP-based solver.

        :param str filename: The filename of where to output the MPS file.
        :raises UnsupportedSolverFunction: if called on a non MIP-based solver.
        """
        if hasattr(self.solver, 'output_mps'):
            if not filename.endswith(".mps"):
                filname = "%s.mps" % filename
            self.solver.output_mps(filename)
        else:
            raise UnsupportedSolverFunction(
                str(type(self)), "output_mps", "This solver does not support "
                "outputing MPS files.")

    def num_vars(self):
        """
        .. deprecated:: 1.1
            Use :meth:`.getNumVariables` intead.
        """
        if hasattr(self.solver, 'num_vars'):
            return self.solver.num_vars()

        raise UnsupportedSolverFunction(
            str(type(self)), "num_vars", "This functionality has been "
            "deprecated, use getNumVariables.")

    def extract_graph(self):
        self.solver.extract_graph()

    def numNodes(self):
        """
        .. deprecated:: 1.1
            Use :meth:`.getNodes` intead.
        """
        return self.solver.numNodes()

    def degree(self, var):
        return self.solver.get_degree(var)

    def get_neighbors(self, x):
        neighbors = []
        for y in range(self.solver.degree(x)):
            neighbors.append(self.solver.get_neighbor(x, y))
        return neighbors

    def get_static_features(self):
        feats = {}
        for i in range(12):
            feats[self.solver.get_feature_name(i)] = self.solver.get_feature(i)
        for i in range(16, 36):
            feats[self.solver.get_feature_name(i)] = self.solver.get_feature(i)
        return feats

    def get_dynamic_features(self):
        feats = {}
        for i in range(12, 16):
            feats[self.solver.get_feature_name(i)] = self.solver.get_feature(i)
        return feats

    def get_features(self):
        """
        Compute and return the 36 CPHydra features. The feature set includes 32
        static features regarding the contraints, domains, etc, and 4 dynamic
        features about weights, nodes, propagations computed after a 2 second
        run of Mistral. Note that these can only be compute when the XCSP
        instance has been loaded directly by Mistral like in the following
        example:

        .. code-block:: python

            model = Model()
            solver = model.load('Mistral')
            solver.load_xml(xml_filename)
            features = solver.get_features()

        .. deprecated:: 1.1

            This will be replaced with more extensive and flexible functionality
            in future releases.
        """
        feats = {}
        for i in range(36):
            feats[self.solver.get_feature_name(i)] = self.solver.get_feature(i)
        return feats

    def __str__(self):
        #setActive(self)
        #return ''
        self.solver.printPython()
        return ' '

    def delete(self):
        if self.free_memory:
            self.free_memory(self.solver)


def enum(*sequential):
    enums = dict(list(zip(sequential, (2 ** i for i in range(len(sequential))))))
    return type('Enum', (), enums)


# This enum ordering must be the same as that specified in the enums
# EncodingConfiguration::AMOEncoding and AllDiffEncoding in SatWrapper.hpp
AMOEncoding = enum('Pairwise', 'Ladder')
AllDiffEncoding = enum('PairwiseDecomp', 'LadderAMO', 'PigeonHole')


class EncodingConfiguration(object):
    """
    Specifies a configuration for the SAT encoding that expressions will take
    when translated to conjunctive normal form for the SAT solver. Support for
    configuring the MIP encoding could be added later.

    At least one of ``direct`` or ``order`` should be ``True``. If both are set
    to ``True``, then the domain will be encoded using both representations and
    chanelled between each other. The default is to have both enabled,
    corresponding to the so called regular encoding.

    At least one of ``conflict`` or ``support`` should be ``True``, these
    specify which form of a constraint is encoded into CNF.

    :param bool direct: Whether the direct (or sparse) encoding of domains
        should be generated.
    :param bool order: Whether the order encoding of the domains should be
        generated.
    :param bool conflict: Whether the conflict clauses of a constraint should be
        generated.
    :param bool support: Whether the support clauses of a constraint should be
        generated.
    :param int amo_encoding: The at-most-one encoding to be used. An enum is
        defined with the supported possibilities: ``AMOEncoding.Pairwise``, and
        ``AMOEncoding.Ladder``.
    :param int alldiff_encoding: The encoding used when encoding all-different
        constraints. The possibilites are defined in the ``AllDiffEncoding``
        enum and can be binary or'd with each other to be passed as a single
        int, like so: ``AllDiffEncoding.PairwiseDecomp |
        AllDiffEncoding.LadderAMO | AllDiffEncoding.PigeonHole``.
    """

    def __init__(self, direct=True, order=True, conflict=True, support=False,
                 amo_encoding=AMOEncoding.Pairwise,
                 alldiff_encoding=AllDiffEncoding.PairwiseDecomp):
        # Domain encodings
        self.direct = direct
        self.order = order

        # Constraint encoding
        self.conflict = conflict
        self.support = support

        # At Most One encoding.
        self.amo_encoding = amo_encoding

        # All Different encoding.
        self.alldiff_encoding = alldiff_encoding

        # Check validity of the encoding config
        if not self.direct and not self.order:
            raise InvalidEncodingException(
                "Domains must be encoded using at least one encoding: "
                "direct|order.")

        if not self.conflict and not self.support:
            raise InvalidEncodingException(
                "Constraints must be encoded using at least one encoding: "
                "conflict|support.")

        if not self.amo_encoding & AMOEncoding.Pairwise and \
           not self.amo_encoding & AMOEncoding.Ladder:
            raise InvalidEncodingException(
                "Invalid at-most-one encoding specified: %s" %
                (str(self.amo_encoding)))

        # if self.amo_encoding & AMOEncoding.Pairwise and not self.direct:
        #     raise InvalidEncodingException("Domains must be encoded using the direct encoding if using the pairwise AMO encoding.")

        # if self.alldiff_encoding & AllDiffEncoding.PairwiseDecomp and not self.direct:
        #     raise InvalidEncodingException("The direct encoding must be enabled if the pairwise decomposition all different is used.")

    # Make EncodingConfiguration hashable so that it can be used as a dictionary
    # key for the cache of encoding configs during translation to SAT.
    def __hash__(self):
        return hash((self.direct, self.order, self.conflict, self.support, self.amo_encoding, self.alldiff_encoding))

    def __eq__(self, other):
        return (self.direct == other.direct) and \
               (self.order == other.order) and \
               (self.conflict == other.conflict) and \
               (self.support == other.support) and \
               (self.amo_encoding == other.amo_encoding) and \
               (self.alldiff_encoding == other.alldiff_encoding)

    def __str__(self):
        return "EncodingConfig<direct:%r, order:%r, conflict:%r, support:%r, amo:%r, alldiff:%r>" % (
            self.direct, self.order, self.conflict, self.support, self.amo_encoding, self.alldiff_encoding)


NJEncodings = {
    # "directsupport": EncodingConfiguration(direct=True, order=False, conflict=True, support=True, amo_encoding=AMOEncoding.Pairwise),
    "direct": EncodingConfiguration(direct=True, order=False, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.PairwiseDecomp),
    "support": EncodingConfiguration(direct=True, order=False, conflict=False, support=True, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.PairwiseDecomp),
    "order": EncodingConfiguration(direct=False, order=True, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.PairwiseDecomp),
    "directorder": EncodingConfiguration(direct=True, order=True, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.PairwiseDecomp),
    "pairwiseconflictandsupport": EncodingConfiguration(direct=True, order=False, conflict=True, support=True, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.PairwiseDecomp),
    "ladder_direct": EncodingConfiguration(direct=True, order=False, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.LadderAMO),
    "ladder_support": EncodingConfiguration(direct=True, order=False, conflict=False, support=True, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.LadderAMO),
    "ladder_directorder": EncodingConfiguration(direct=True, order=True, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.LadderAMO),
    "pairwise_and_ladder_direct": EncodingConfiguration(direct=True, order=False, conflict=True, support=False, alldiff_encoding=AllDiffEncoding.PairwiseDecomp | AllDiffEncoding.LadderAMO),
    "pairwise_ladder_pigeon": EncodingConfiguration(direct=True, order=True, conflict=True, support=False, alldiff_encoding=AllDiffEncoding.PairwiseDecomp | AllDiffEncoding.LadderAMO | AllDiffEncoding.PigeonHole),
    "pairwisesupport_ladder_pigeon": EncodingConfiguration(direct=True, order=True, conflict=False, support=True, alldiff_encoding=AllDiffEncoding.PairwiseDecomp | AllDiffEncoding.LadderAMO | AllDiffEncoding.PigeonHole),
    "pairwise_order_pigeon": EncodingConfiguration(direct=False, order=True, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.PairwiseDecomp | AllDiffEncoding.PigeonHole),
    "pairwise_directorder_pigeon": EncodingConfiguration(direct=True, order=True, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.PairwiseDecomp | AllDiffEncoding.PigeonHole),
    "ladder_directorder_pigeon": EncodingConfiguration(direct=True, order=True, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise, alldiff_encoding=AllDiffEncoding.LadderAMO | AllDiffEncoding.PigeonHole),
}
