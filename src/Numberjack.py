##@mainpage Numberjack  @authors Eoin O'Mahony, Emmanuel Hebrard & Barry O'Sullivan 
#
# \section intro_sec What is numberjack?
#
# Numberjack is a modelling package written in Python for constraint
# programming. Python benefits from a large and active programming
# community, Numberjack is therefore a perfect tool to embed CP
# technology into larger applications. It is designed to support a
# number of underlying C/C++ solvers as egg files, that is, seamlessly
# and efficiently. Currently, there are four available back-ends: a MIP
# solver (SCIP), a SAT solver (MiniSat) and two CP solvers (Mistral and
# Gecode).
#
# \section license_sec License:
#
#  Numberjack is a constraint satisfaction and optimisation library
#  Copyright (C) 2009 Cork Constraint Computation Center, UCC
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#  The authors can be contacted electronically at 
#  numberjack.support@gmail.com


UNSAT     =  0;
SAT       =  1;
UNKNOWN   =  2;
LIMITOUT  =  3;
LUBY      =  0; 
GEOMETRIC =  1;

#from sets import Set
import types, sys

val_heuristics = ['Lex', 'AntiLex', 'Random', 'RandomMinMax', 'DomainSplit', 'RandomSplit', 'Promise', 'Impact', 'No']
var_heuristics = ['No', 'MinDomain', 'Lex', 'AntiLex', 'MaxDegree', 'MinDomainMinVal', 'Random', 'MinDomainMaxDegree', 'DomainOverDegree', 'DomainOverWDegree', 'DomainOverWLDegree', 'Neighbour', 'Impact', 'ImpactOverDegree', 'ImpactOverWDegree', 'ImpactOverWLDegree', 'Scheduling']




def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring) and not issubclass(type(el), Expression):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result



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
            list.__init__(self, [arg1,arg2])
            self.is_bound = True
        self.current = -1

    def next(self):
        """
        \internal
        Returns the next value when iterating
        """
        self.current += 1
        if self.is_bound:
            if self[0]+self.current > self[-1]: raise StopIteration
            else: return self[0]+self.current
        else:
            if self.current >= list.__len__(self): raise StopIteration
            else: return list.__getitem__(self, self.current)

    def __str__(self):
        """
        \internal
        """
        if self.is_bound:
            lb = self[0]
            ub = self[-1]
            if lb+1 == ub and type(lb) is int:
                return '{'+str(lb)+','+str(ub)+'}'
            else:
                return '{'+str(lb)+'..'+str(ub)+'}'

        def extend(idx):
            x = self[idx]
            y = x
            idx += 1
            while idx < len(self):
                if type(self[idx]) is int and self[idx] == y+1: y=self[idx]
                else: break
                idx += 1
            return (x,y,idx)

        ret_str = '{'
        idx = 0
        while idx < len(self):
            if idx > 0: ret_str += ','
            (x,y,idx) = extend(idx)
            ret_str += str(x)
            if type(x) is int and x+1 < y: ret_str += ('..'+str(y))
            elif x != y: ret_str += (','+str(y))

        return ret_str+'}'


## Base Expression class from which everything inherits
#
#    All variables and constraints are expressions. Variables are
#    just expressions where predicates extend the expression class to add
#    more functionality.
#   
class Expression(object):
    
    def __init__(self, operator):
        self.ident = -1
        self.operator = operator
        
        # This is the stuff for maintaining multiple representations of the
        # model among different solvers
        self.var_list = []

    def __iter__(self):
        return self.get_domain()

    def get_solver(self):
        if getattr(self, 'solver', False):
            return self.solver
        else:
            return None


    ## Returns a string representing the initial domain of the expression
    # @return String representation of original expression definition
    #
    #    This method returns a string representing the original definition
    #    of the expression.
    #    
    #    For example:
    # \code
    #    var1 = Variable(0, 10)
    #    print var1.initial()
    #    >>> x0 in {0..10}
    # \endcode
    def initial(self):
        output = self.name()
        if self.domain_ is None: output += ' in ' + str(Domain(self.lb, self.ub))
        else: output += ' in ' + str(Domain(self.domain_))
        return output
    
    ## Returns a string representing the current domain of the expression
    # @param solver Solver from which current expression domain will be sourced
    # @return String representing of current expression domain
    #
    #    domain(self, solver=None) :- Returns a string representing the current domain
    #    of the expression in the solver specified. If no solver is specified then
    #    the returned string represents the domain of the expression in the solver
    #    that has most recently loaded the expression.
    #
    def domain(self,solver=None):
        output = self.name() + ' in ' + str(self.get_domain())
        return output
    
    ## Returns a string containing the value of the expression
    # @param solver Solver from which expression solution will be sourced
    # @return String represent-ion of expression solution
    #
    #    solution(self, solver=None) :- Returns a string representing the current solution
    #    of the expression in the solver specified. If no solver is specified then
    #    the returned string represents the solution to the expression in the solver
    #    that has most recently loaded and solved the expression.
    #
    def solution(self, solver=None):
        return str(self.get_value(solver))

    def name(self):
        output = self.operator
        if (output == 'x' or output == 't') and self.ident >= 0:
            output += str(self.ident)
        return output 

    def __str__(self):
        if self.is_built() and self.solver.is_sat():
            return self.solution()
        else: return self.domain()

    def is_str(self):
        return hasattr(self, 'model')
 
    def getVar(self, solver_id):
        return self.var_list[solver_id-1]
    
    def setVar(self, solver_id, solver_name, variable, new_solver=None):
        if (solver_id-1) < len(self.var_list):
            self.var_list[solver_id-1] = variable
        else:
            self.var_list.append(variable)

    def has_children(self):
        return hasattr(self, 'children')

    def has_parameters(self):
        return hasattr(self, 'parameters')

    def is_built(self, solver=None):
        if solver is None:
            return len(self.var_list) > 0
        else: 
            return solver.solver_id - 1 < len(self.var_list)

    def is_var(self):
        return not issubclass(type(self), Predicate)
        
    def close(self):
        if self.has_children():
            for child in self.children:
                tc = type(child)
                if tc is not int and tc is not str and tc is not bool:
                    child.close()

    def get_domain(self, solver=None):
        if self.is_built(solver):
            if solver is None: solver = self.solver
            if self.get_size(solver) == (self.get_max(solver) - self.get_min(solver) + 1):
                dom = range(self.get_min(solver), self.get_max(solver)+1)
            else:
                # we should make that more efficient by using the underlying solvers to iterate
                #dom = [v for v in range(self.get_min(solver), self.get_max(solver)+1) if v in self] 
                dom = [self.get_min(solver)]
                while True: 
                    v = solver.next(self, dom[-1])
                    if v == dom[-1]: break
                    else: dom.append(v)
            return Domain(dom)
        elif self.domain_ is not None: 
            return Domain(self.domain_)
        else: 
            return Domain(self.lb, self.ub) #range(self.lb, self.ub+1))

    ## Current value of the expression
    # @param solver solver reference for solver from which current value will be sourced
    # @return current value of expression
    #
    #    get_value(self, solver=None) :- Returns the current value of the expression
    #    in the specified solver. Provided the solver has found a solution. If no solution
    #    has been found None is returned.
    #    
    #    If not solver is specified then the solver that has loaded the expression
    #    most recently is used. 
    #
    def get_value(self, solver=None):
        has_value = False        
        if self.is_built(solver):
            if self.solver.is_sat():
                has_value = True
        value = None
        
        if solver is not None:
            var = self.var_list[solver.solver_id-1]
        else: var = self.var_list[-1]    
            
        if has_value:
            if self.is_str():
                value = self.model.strings[var.get_value()]
            else:
                value = var.get_value()
        return value


    def get_size(self, solver=None):
        if solver is not None: 
            if self.is_built(solver):
                return self.var_list[solver.solver_id-1].get_size()
            else:
                return self.ub - self.lb + 1
        elif self.is_built():
            return self.var_list[-1].get_size()
        else:
            return self.ub - self.lb + 1

    ## Current lower bound of variable
    # @param solver solver reference for solver from which current bound will be sourced
    # @return current lower bound of expression
    #
    #    get_min(self, solver=None) :- Returns the current lower bound of the expression
    #    in the specified solver.
    #    
    #    If not solver is specified then the solver that has loaded the expression
    #    most recently is used. 
    #
    def get_min(self, solver=None):
        the_min = self.lb
        if solver is not None: 
            if self.is_built(solver):
                the_min = self.var_list[solver.solver_id-1].get_min()
        elif self.is_built():
            the_min = self.var_list[-1].get_min()
        if self.is_str():
            return self.model.strings[the_min]
        return the_min

    ## Current upper bound of variable
    # @param solver solver reference for solver from which current bound will be sourced
    # @return current upper bound of expression
    #
    #    get_max(self, solver=None) :- Returns the current upper bound of the expression
    #    in the specified solver.
    #    
    #    If not solver is specified then the solver that has loaded the expression
    #    most recently is used. 
    #
    def get_max(self, solver=None):
        the_max = self.ub
        if solver is not None: 
            if self.is_built(solver):
                the_max = self.var_list[solver.solver_id-1].get_max()
        elif self.is_built():
            the_max = self.var_list[-1].get_max()
        if self.is_str():
            return self.model.strings[the_max]
        return the_max

    # not safe!
    def get_ub(self):
        #print "DEPRECATED: Expression.get_ub()"
        return self.ub
    
    # not safe!
    def get_lb(self):
        #print "DEPRECATED: Expression.get_lb()"
        return self.lb

    # not safe!
    def get_domain_tuple(self):
        if self.is_str():
            tmp_domain = sorted([self.model.string_map[value] for value in self.domain_])
            return (tmp_domain[0], tmp_domain[len(tmp_domain)-1], tmp_domain)
        else:
            return (self.lb, self.ub, self.domain_)

    def get_children(self):
        if self.has_children():
            return self.children
        else:
            return None
        
    def get_operator(self):
        return self.operator

    def __and__(self, pred):
        return And([self,pred])

    def __or__(self, pred):
        return Or([self,pred])

    def __add__(self, pred):
        return Sum([self, pred], [1, 1])

    def __radd__(self, pred):
        return Sum([pred,self], [1, 1])
        
    def __sub__(self, pred):
        var = Sum([self, pred], [1, -1])
        var.name = '('+str(self)+'-'+str(pred)+')'
        return var

    def __rsub__(self, pred):
        return Sum([pred, self], [1, -1])

    def __div__(self, pred):
        return Div([self,pred])

    def __rdiv__(self, pred):
        return Div([pred,self])
        
    def __mul__(self, pred):
        return Mul([self,pred])

    def __rmul__(self, pred):
        return Mul([self,pred])

    def __mod__(self, pred):
        return Mod([self,pred])

    def __rmod__(self, pred):
        return Mod([pred,self])
        
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
        if self.is_built(): return self.var_list[-1].contain(v)
        elif self.domain_ is None: return self.lb <= v <= self.ub
        else: return v in self.domain_

            

## @defgroup mod_group Modelling constructs
# The modelling constructs
#  @{
#

## Model object
#
#    Stores the variables and constraints.
#    The constraints declarations are trees, whose internal nodes 
#    are predicates or constraints and leaves are variables.
#
#    - self.variables contains the leaves of these trees.
#
#    - self.constraints contains the internal nodes of these trees.
#
#    The constructor of Model can be called with any number of arguments.
#    Each argument will be treated as an Expression or a list of 
#    Expressions to be added into the model. If no argument is given
#    the Model will be initially empty.
#    If X and Y are lists of Expressions, then the following declarations are all valid:
#
# \code
#    m = Model()
#    m = Model( X[0]<Y[0] )
#    m = Model( X[0]<Y[0], Y[0]!=Z[0] )
#    m = Model( [x<y for (x,y) in zip(X,Y)] )
# \endcode
#
#    Expressions and lists of Expressions can be subsequently added 
#    using the method self.add() or the operator '+='.
#
class Model(object):


    def __init__(self, *expr):
        ## \internal - List of expressions (predicate trees) that where added to the model
        self.__expressions = []

        ## Leaves of the predicate trees 
        self.variables = VarArray([])

        ## Roots of the predicate trees
        self.constraints = VarArray([])

        ## \internal - Before giving an expression to a solver, or before printing it, it needs to
        self.closed = 0

        ## \internal - Every new solver get assigned an unique id
        self.current_id = 0

        ## \internal - Initialise from an expression?
        if len(expr) > 0:
            self.add_prime(expr)
    


    def getSolverId(self):
        # \internal \internal - generates an ident for each new solver
        self.current_id += 1
        return self.current_id
        
    ## Add an expresion, or a list/tuple/dict of expressions to the model
    # @param *expr Any number of (nested lists of) Expressions
    def add(self, *expr):
        self.add_prime(expr)
       
    def add_prime(self, expr):
        ## \internal - Used to distinguish between a single Expression and a list of Expressions
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
        ## \internal - '+=' operator
        self.add_prime(expr)
        return self

    def add_expression(self, exp, level):
        ## \internal - add the Expression tree to the model and assign identifiers to the nodes
        # this expression is new, choose an identifiant for it
        te = type(exp)
        if te is not int and te is not str and te is not bool:
            if exp.ident == -1:            
                if exp.get_children() is None:
                    if exp.is_var():
                        exp.ident = len(self.variables)
                        self.variables.append(exp)
                else: # it is a constraint
                    exp.ident = -2-len(self.constraints)
                    self.constraints.append(exp)
                    for child in exp.get_children():
                        self.add_expression(child,level+1)

    def close_exp(self):
        ## \internal - close() is used to fire up preprocessing requiring knowledge about the whole model
        for i in range(self.closed, len(self.__expressions)):
            self.add_expression(self.__expressions[i], 0)
        self.closed = len(self.__expressions)
        
    def close(self):
        ## \internal - close() is used to fire up preprocessing requiring knowledge about the whole model
        if self.closed == len(self.__expressions):
            tmp_strings = []
            for var in self.variables:
                if var.get_lb() is None:
                    var.model = self
                    for value in var.domain_:
                        tmp_strings.append(value)

            self.strings = sorted(set(tmp_strings))
            self.string_map = {}.fromkeys(self.strings)
                
            for k in range(len(self.strings)):
                self.string_map[self.strings[k]] = k

        self.closed += 1

        if self.closed == len(self.__expressions)+1:
            for expr in self.get_exprs():
                expr.close()
        
    def __str__(self):
        ## \internal - print
        mod = 'assign:\n  '
        for var in self.variables:
            mod += var.domain()+'\n  '
        mod += '\nsubject to:\n  '
        for con in self.__expressions:
            mod += con.__str__()+'\n  '
        return mod

    def get_exprs(self):
        ## \internal - return the list of Expressions
        return self.__expressions
    
    def is_native_model(self):
        ## \internal - distinguish with NumberjackSolver models (deprecated?)
        return False

    ## Return a Solver for that Model
    # @param library A string standing for a back end solver ('Mistral', 'MiniSat', 'SCIP')
    def load(self,library,X=None):
        """
        The solver is passed as a string, the corresponding module is
        imported, a Solver object created, initialised and returned
        """
        try:
            lib =  __import__(library)
            solver = lib.Solver(self,X)
        except ImportError:
            print "ERROR: Failed during import, wrong module name? ("+library+")"
            sys.exit(1)
        return solver

    ## Solve and returns a solution 
    # @param library A string standing for a back end solver ('Mistral', 'MiniSat', 'SCIP')
    # @return A Solution object
    def solve_with(self,library):
        """
        The solver is passed as a string, the corresponding module is
        imported, a Solver object created, initialised and called.
        A Solution object (dictionary: var -> val) is returned, if the
        Model is unsatisfiable, or the solver fails to solve it, the
        Solution will be empty (None)
        """
        solver = self.load(library)
        solver.solve()
        return solver.get_solution()


## Variable creation function
# Creates a variable, see below for usage. To create a continuous variable
# for solvers capable of handling continuous variables use float values
# instead of integer values when specifying the domain of the variable.

class Variable(Expression):

    def __init__(self, argopt1=None, argopt2=None, argopt3=None):
        '''
         Variable() :- Binary variable
         Variable(N) :- Variable in the domain of {0, N-1}
         Variable('x') :- Binary variable called 'x'
         Variable(N, 'x') :- Variable in the domain of {0, N-1} called 'x'
         Variable(l,u) :- Variable in the domain of {l, u}
         Variable(l,u, 'x') :- Variable in the domain of {l, u} called 'x'
         Variable(list) :- Variable with domain specified as a list
         Variable(list, 'x') :- Variable with domain specified as a list called 'x'
        '''

        def numeric(x): 
            tx = type(x)
            return tx is int or tx is float

        domain=None
        lb=0
        ub=1
        name = 'x'

        if argopt3 is not None:
            lb = argopt1
            ub = argopt2
            name  = argopt3
        elif argopt2 is not None:
            if type(argopt2) is str:
                if numeric(argopt1):
                    ub = argopt1-1
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
                ub = argopt1-1
            else: 
                domain = sorted(argopt1)
                lb = domain[0]
                ub = domain[-1]            

        if type(lb) is not int and type(lb) is not float and type(lb) is not str:
            print "Warning lower bound of %s is not an int or a float or a string" % name
            exit(1)
        elif type(ub) is not int and type(ub) is not float and type(lb) is not str:
            print "Warning upper bound of %s is not an int or a float or a string" % name
            exit(1)
        elif type(name) is not str:
            print "Warning name variable is not a string"
            exit(1)
        elif lb > ub:
            print "Warning lower bound of %s greater than upper bound " % name
            exit(1)

        Expression.__init__(self, name)
        self.domain_ = domain
        self.lb = lb
        self.ub = ub


'''
def Variable(argopt1=None, argopt2=None, argopt3=None):
    
     Variable() :- Binary variable
     Variable(N) :- Variable in the domain of {0, N-1}
     Variable('x') :- Binary variable called 'x'
     Variable(N, 'x') :- Variable in the domain of {0, N-1} called 'x'
     Variable(l,u) :- Variable in the domain of {l, u}
     Variable(l,u, 'x') :- Variable in the domain of {l, u} called 'x'
     Variable(list) :- Variable with domain specified as a list
     Variable(list, 'x') :- Variable with domain specified as a list called 'x'
    

    def numeric(x): 
        tx = type(x)
        return tx is int or tx is float

    domain=None
    lb=0
    ub=1
    name = 'x'

    if argopt3 is not None:
        lb = argopt1
        ub = argopt2
        name  = argopt3
    elif argopt2 is not None:
        if type(argopt2) is str:
            if numeric(argopt1):
                ub = argopt1-1
            else: domain = sorted(argopt1)
            name = argopt2
        else:
            lb = argopt1
            ub = argopt2
    elif argopt1 is not None: 
        if type(argopt1) is str:
            name = argopt1
        elif numeric(argopt1): 
            ub = argopt1-1
        else: domain = sorted(argopt1)

    exp = Expression(name)

    exp.domain_ = domain
    exp.lb = lb
    exp.ub = ub

    return exp
'''


## Array of Expression
#
#    A VarArray is a list of Expression objects. 
#    Various methods are overloaded to allow easy declaration, formatted printing 
#    and modelling syntactic sugars
#    
#    There are various ways of declaring a VarArray:
#
#    - X = VarArray(l) creates an array from a list l
#    - X = VarArray(n) creates an array of n Boolean variables
#    - X = VarArray(n, 'x') creates an array of n Boolean variables with names 'x0..xn-1'
#    - X = VarArray(n, m, 'x') creates an array of n variables with domains [0..m-1] and names 'x0..xn-1'
#    - X = VarArray(n, m) creates an array of n variables with domains [0..m-1] 
#    - X = VarArray(n, d, 'x') creates an array of n variables with domains d and names 'x0..xn-1'
#    - X = VarArray(n, d) creates an array of n variables with domains d 
#    - X = VarArray(n, l, u, 'x') creates an array of n variables with domains [l..u] and names 'x0..xn-1'
#    - X = VarArray(n, l, u) creates an array of n variables with domains [l..u] 
#
#    VarArray's allow to state Element and Lex Ordering constraint over a sequence
#    of variables using, respectively the operator '[]' and '<', '>', '<=', '>='.
#    For instance, given two VarArray X and Y, and an Expression x:
#    
#    X[x] returns an Element Expression, that is, a variable equal to the xth element of the array X
#    X <= Y returns a LeqLex Constraint between X and Y
#
class VarArray(list):

    def __init__(self, n, optarg1=None, optarg2=None, optarg3=None):
        domain = None
        if hasattr(n,'__iter__'):
            list.__init__(self, n)
            return
        else:
            lb = 0
            ub = 1
            name = 'x'
            if optarg1 is not None:
                #if type(optarg1) is str:
                if hasattr(optarg1, '__iter__'):
                    name = optarg1
                elif type(optarg2) is int or type(optarg2) is float:
                    lb = optarg1
                    ub = optarg2
                    if optarg3 is not None: name = optarg3
                else:
                    if issubclass(type(optarg1), list) or issubclass(type(optarg1), list):
                        domain = optarg1
                        domain.sort()
                        lb = domain[0]
                        ub = domain[-1]
                    else: ub = optarg1-1
                    if optarg2 is not None: name = optarg2
        names = name
        if type(name) is str:
            names = [name+str(i) for i in range(n)]
        if domain is None:
            self.__init__([Variable(lb, ub, names[i]) for i in range(n)])
        else:
            self.__init__([Variable(domain, names[i]) for i in range(n)])

    ## Returns a string representing the initial definition of the content of the arrray 
    #@return string
    def initial(self):
        return "[" + ", ".join([var.initial() for var in self]) + "]"
    
    ## Returns a string representing the current state of the content of the array
    #@param solver The linked Solver 
    #@return string
    def domain(self, solver=None):
        return "[" + ", ".join([var.domain(solver) for var in self]) + "]"
    
    ## Returns a string containing a brief view of the content of the array
    #@return string
    def name(self):
        return "[" + ", ".join([var.name() for var in self]) + "]"

    ## Returns a string containing the valuation of the content of the array
    #@param solver The linked Solver 
    #@return string
    def solution(self, solver=None):
        return "[" + ", ".join([var.solution(solver) for var in self]) + "]"

    def __str__(self):
        return "[" + ", ".join([var.__str__() for var in self]) + "]"

    def __getitem__(self, expr):
        if type(expr) is int:
            return list.__getitem__(self, expr)
        else:
            return Element(self, expr)

    def __getslice__(self, i, j):
        return VarArray(list.__getslice__(self,i,j))
    
    ## Syntactic sugar for the Lex Order Constraint X < Y
    #@param other Another VarArray of same length
    #@return LessLex Expression
    def __lt__(self, other):
        return LessLex(self, other)

    ## Syntactic sugar for the Lex Order Constraint X <= Y
    #@param other Another VarArray of same length
    #@return LessLex Expression
    def __le__(self, other):
        return LeqLex(self, other)

    ## Syntactic sugar for the Lex Order Constraint X > Y
    #@param other Another VarArray of same length
    #@return LessLex Expression
    def __gt__(self, other):
        return LessLex(other, self)

    ## Syntactic sugar for the Lex Order Constraint X >= Y
    #@param other Another VarArray of same length
    #@return LessLex Expression
    def __ge__(self, other):
        return LeqLex(other, self)

    ## Syntactic sugar for the Equality Constraint X == Y
    #@param other Another VarArray of same length
    #@return a list of Equality Expressions
    def __eq__(self, other):
        return [Eq((x,y)) for (x,y) in zip(self,other)]
        

## Matrix of Expression
#
#    A Matrix is a two-dimensional list of Expression objects. 
#    Various methods are overloaded to allow easy declaration, formatted printing 
#    and modelling syntactic sugars
#    
#    There are various ways of declaring a Matrix:
#
#    - \code M = Matrix(l) \endcode creates a Matrix from a list l
#    - M = Matrix(n, m) creates a n x m Matrix of Boolean variables
#    - M = Matrix(n, m, 'x') creates a n x m Matrix of Boolean variables with names 'x0_0..xn-1_m-1'
#    - M = Matrix(n, m, u) creates a n x m Matrix of variables with domains [0..u-1] 
#    - M = Matrix(n, m, u, 'x') creates a n x m Matrix of variables with domains [0..u-1] and names 'x0_0..xn-1_m-1'
#    - M = Matrix(n, m, l, u) creates a n x m Matrix of variables with domains [l..u]
#    - M = Matrix(n, m, l, u, 'x') creates a n x m Matrix of variables with domains [l..u] and names 'x0_0..xn-1_m-1'
#
#    Matrices feature specific handlers to access (subsets of) rows and columns. 
#    The fields 'row', 'col' and 'flat' respectively refer to the list of rows, 
#    columns and cell in the matrix. For instance:
#
# \code
#    m = Matrix(5,4,1,3,'cell_')
#    print m
#    >>> [[cell_0.0, cell_0.1, cell_0.2, cell_0.3],
#    >>>  [cell_1.0, cell_1.1, cell_1.2, cell_1.3],
#    >>>  [cell_2.0, cell_2.1, cell_2.2, cell_2.3],
#    >>>  [cell_3.0, cell_3.1, cell_3.2, cell_3.3],
#    >>>  [cell_4.0, cell_4.1, cell_4.2, cell_4.3]]
#    print m.row
#    >>> [[cell_0.0, cell_0.1, cell_0.2, cell_0.3],
#    >>>  [cell_1.0, cell_1.1, cell_1.2, cell_1.3],
#    >>>  [cell_2.0, cell_2.1, cell_2.2, cell_2.3],
#    >>>  [cell_3.0, cell_3.1, cell_3.2, cell_3.3],
#    >>>  [cell_4.0, cell_4.1, cell_4.2, cell_4.3]]
#    print m.col
#    >>> [[cell_0.0, cell_1.0, cell_2.0, cell_3.0, cell_4.0],
#    >>>  [cell_0.1, cell_1.1, cell_2.1, cell_3.1, cell_4.1],
#    >>>  [cell_0.2, cell_1.2, cell_2.2, cell_3.2, cell_4.2],
#    >>>  [cell_0.3, cell_1.3, cell_2.3, cell_3.3, cell_4.3]]
#    print m.flat
#    >>> [cell_0.0, cell_0.1, cell_0.2, cell_0.3, cell_1.0, cell_1.1, cell_1.2, cell_1.3, cell_2.0, cell_2.1, cell_2.2, cell_2.3, cell_3.0, cell_3.1, cell_3.2, cell_3.3, cell_4.0, cell_4.1, cell_4.2, cell_4.3]
# \endcode
#
#    Matrices support Element constraints on row, column or flatten views.
class Matrix(list):

    ## accessor to the list of rows
    row = None

    ## accessor to the list of columns
    col = None

    ## accessor to the list of cells
    flat = None


    def __init__(self, optarg1=None, optarg2=None, optarg3=None, optarg4=None, optarg5=None):
        n = 1
        m = 1
        lb = 0
        ub = 1
        name = 'x'

        if optarg2 == None:
            if optarg1 != None:
                list.__init__(self, [VarArray(row) for row in optarg1])
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
                    ub = optarg3-1
                else:
                    ub = optarg4
                    lb = optarg3
            elif optarg3 is not None:
                if type(optarg3) is str: name = optarg3
                else: ub = optarg3-1
            list.__init__(self, [VarArray(m, lb, ub, name+str(j)+'.') for j in range(n)])
        self.row = self
        self.col = Matrix()
        for column in zip(*self):
            self.col.append(VarArray(column))
        self.col.col = self
        self.col.row = self.col
        self.flat = VarArray([var for row in self for var in row])
        self.col.flat = self.flat
        
    ## Returns a string representing the initial definition of the content of the arrray 
    #@return string
    def initial(self):
        return "[" + ",\n ".join([row.initial() for row in self]) + "]"
        
    ## Returns a string representing the current state of the content of the matrix
    #@param solver The linked Solver 
    #@return string
    def domain(self, solver=None):
        return "[" + ",\n ".join([row.domain(solver) for row in self]) + "]"

    ## Returns a string containing a brief view of the content of the matrix
    #@return string
    def name(self):
        return "[" + ",\n ".join([row.name() for row in self]) + "]"

    ## Returns a string containing the valuation of the content of the matrix
    #@param solver The linked Solver 
    #@return string
    def solution(self, solver=None):
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
                aux = Matrix(list.__getitem__(self,i[0])).col
                aux = Matrix(list.__getitem__(aux,i[1])).col
                return aux
            else: 
                return Element(self.flat, (i[0] * len(self.col)) + i[1])
        else:
            return MatrixWrapper(i, self)


    def __getslice__(self, i, j):
        return Matrix(list.__getslice__(self,i,j))

     
class MatrixWrapper(list):
    
    def __init__(self, var, matrix):
        self.var = var
        self.matrix = matrix
        
    def __getitem__(self, item):
        return self.matrix[self.var, item]
        
    def __str__(self):
        return str(self.var) + " th index of " + str(self.matrix)

##  @}    




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

    def __init__(self, vars, op):
        Expression.__init__(self, op)
        self.set_children(vars)

    def set_children(self, children):
        #self.children = children
        
        ## List of children of the predicate 
        #self.children = [child for child in children]
        self.children = flatten(children)


    ## Returns a string representing the initial definition of the Predicate
    # @return String representation of original predicate definition
    #
    #    Returns a string that represents the initial definition of the Predicate
    #    and all its children.
    #    
    # \code    
    #    var1 = Variable(0, 10)
    #    constraint = var1 < 10
    #    model = Model(constraint)
    #    solver = Solver(model)
    #    solver.solve()
    #    .
    #    .
    #    print constraint.initial()
    #    >>>(x0 in {0..10} < 10)
    # \endcode
    #
    def initial(self):
        save_str = Expression.__str__
        Expression.__str__ = Expression.initial
        output = self.__str__()
        Expression.__str__ = save_str
        return output
    
    ## Returns a string representing the current representation of the predicate
    # @param solver solver from which current domains of leaves will be sourced
    # @return String representation of current predicate state
    #
    #    Returns a string that represents the current value of the predicate and its'
    #    children in the specified solver. If no solver is specified then the solver
    #    that has loaded the expression the solver most recently is used.
    #    
    # \code
    #    var1 = Variable(0, 10)
    #    constraint = var1 < 10
    #    model = Model(constraint)
    #    solver = Solver(model)
    #    solver.solve()
    #    .
    #    .
    #    print constraint.domain()
    #    >>> (x0 in {0..9} < 10)
    # \endcode
    #
    def domain(self, solver=None):
        save_str = Expression.__str__
        Expression.__str__ = lambda x: x.domain(solver)
        output = self.__str__()
        Expression.__str__ = save_str
        return output
    
    ## Returns a string representing the current solution of the predicate
    # @param solver solver from which current domains of leaves will be sourced
    # @return String representation of the current predicate solution
    #
    #    Returns a string that represents the current solution of the predicate and its'
    #    children in the specified solver. If no solver is specified then the solver
    #    that has loaded the expression the solver most recently is used.
    #
    # \code
    #    var1 = Variable(0, 10)
    #    constraint = var1 < 10
    #    model = Model(constraint)
    #    solver = Solver(model)
    #    solver.solve()
    #    .
    #    .
    #    print constraint.solution()
    #    >>> (0 < 10)
    # \endcode
    #
    def solution(self, solver=None):
        save_str = Expression.__str__
        Expression.__str__ = lambda x: x.solution(solver)
        output = self.__str__()
        Expression.__str__ = save_str
        return output
    
    ## Returns a string representing the name of the Predicate
    # @return String representation of predicate definition
    #
    #    Returns a string that represents the name of the Predicate
    #    and the name of all its children.
    #    
    # \code    
    #    var1 = Variable(0, 10)
    #    constraint = var1 < 10
    #    model = Model(constraint)
    #    solver = Solver(model)
    #    solver.solve()
    #    .
    #    .
    #    print constraint.name()
    #    >>>(x0 < 10)
    # \endcode
    #
    def name(self):
        return self.__str__()

    def __str__(self):
        save_str = Expression.__str__
        Expression.__str__ = Expression.name
        output = self.operator + "(" + ", ".join(map(str, self.children)) + ")"
        Expression.__str__ = save_str
        return output

## Class that all binary predicates inherit from
#
#    All binary predicates such as And, LessThan and GreaterThan extend this
#    class. They are separated from the base Predicate class to facilitate
#    ease of print representations of the predicates.
#
class BinPredicate(Predicate):
    
    def __init__(self, vars, op):
        Predicate.__init__(self, vars, op)

    def get_symbol(self): return 'x'

    def __str__(self):
        save_str = Expression.__str__
        Expression.__str__ = Expression.name
        output = '('+str(self.children[0])+' '+self.get_symbol()+' '+str(self.children[1])+')'
        Expression.__str__ = save_str
        return output

## And expression
#
# \note
#   - Top-level: 'And' Constraint
#   - Nested: Equal to the truth value of the 'And' relation
#
# \code    
#    var1 = Variable() :- Binary variable
#    var2 = Variable() :- Binary variable
#    
#    model.add(var1 & var2) :- Post var1 And var2 constraint
#    
#    var1 = Variable() :- Binary variable
#    var2 = Variable() :- Binary variable
#    var3 = Variable() :- Binary variable
#    
#    model.add( var3 == (var1 & var2) ) :- Used as an expression
# \endcode
#
class And(BinPredicate):
    
    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "and")

    def get_symbol(self): return '&'

## @defgroup cons_group Constraint list
# The list of constraints defined in Numberjack
#  @{
#

## Or expression
#
# \note
#   - Top-level: 'Or' Constraint
#   - Nested: Equal to the truth value of the 'Or' relation
#    
# \code
#    var1 = Variable() :- Binary variable
#    var2 = Variable() :- Binary variable
#    
#    model.add(var1 | var2) :- Post var1 Or var2 constraint
#    
#    var1 = Variable() :- Binary variable
#    var2 = Variable() :- Binary variable
#    var3 = Variable() :- Binary variable
#    
#    model.add( var3 == (var1 | var2) ) :- Used as an expression
# \endcode
#
class Or(BinPredicate):
    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "or")
    
    def get_symbol(self): return 'or'

## Div expression
#
# \note
#   - Top-level: Can not be used as top-level Constraint
#   - Nested: Equal to the integral division of the operands
#
# \warning Can not be used with all solvers
#
#    Div expression can be used to divide two expressions or an expression and
#    a constraint.
#    
# \code
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    
#    divexp1 = var2 / var1
#    divexp2 = var2 / 10
# \endcode
#
class Div(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "div")
    
    def get_symbol(self): return '/'

## Mul expression
#
# \note
#   - Top-level: Can not be used as top-level Constraint
#   - Nested: Equal to the multiplication of the operands
# 
# \warning Can not be used with all solvers
#
#    Mul expression can be used to multiply two expressions or an expression and
#    a constraint.
#
# \code    
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    
#    mulxp1 = var2 * var1
#    mulexp2 = var2 * 10
# \endcode
#
class Mul(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "mul")
    
    def get_symbol(self): return '*'

## Mod expression
# \note
#   - Top-level: Can not be used as top-level Constraint
#   - Nested: Equal to the modulo of the operands
#
# \warning Can not be used with all solvers
#
#    Mod expression can be used to modulo two expressions or an expression and
#    a constraint.
#
# \code    
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    
#    modxp1 = var2 % var1
#    modexp2 = var2 % 10
# \endcode
#
class Mod(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "mod")
    
    def get_symbol(self): return '%'

## Equal expression
#
# \note
#   - Top-level: Equality Constraint
#   - Nested: Equal to the truth value of the Equality relation between the arguments
#
#    Eq expression can be used to create an equality expression between two
#    expressions or an expression and a constraint. It can be used as either
#    a hard constraint or as an expression.
#
# \code    
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    var3 = Variable()
#    
#    model.add(var1 == var2)
#    model.add(var1 == 5)
#
#    model.add( var3 == (var1 == var2) )    
# \endcode
#
class Eq(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "eq")
    
    def get_symbol(self): return '=='

## Not Equal expression
#
# \note
#   - Top-level: Disequality Constraint
#   - Nested: Equal to the truth value of the Disequality relation between the arguments
#
#    Neq expression can be used to create an inequality expression between two
#    expressions or an expression and a constraint. It can be used as either
#    a hard constraint or as an expression.
#
# \code    
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    var3 = Variable()
#    
#    model.add(var1 != var2)
#    model.add(var1 != 5)
#
#    model.add( var3 != (var1 != var2) )    
# \endcode
#
class Ne(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "ne")
    
    def get_symbol(self): return '=/='

## Less than expression
#
# \note
#   - Top-level: Less than Constraint
#   - Nested: Equal to the truth value of the less than relation between the arguments
#
#    Lt expression can be used to create an less than expression between two
#    expressions or an expression and a constraint. It can be used as either
#    a hard constraint or as an expression.
#
# \code    
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    var3 = Variable()
#    
#    model.add(var1 < var2)
#    model.add(var1 < 5)
#
#    model.add( var3 < (var1 < var2) )
# \endcode
#
class Lt(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "lt")
    
    def get_symbol(self): return '<'

## Greater than expression
#
# \note
#   - Top-level: Greater than Constraint
#   - Nested: Equal to the truth value of the Greater than relation between the arguments
#
#    Gt expression can be used to create an greater than expression between two
#    expressions or an expression and a constraint. It can be used as either
#    a hard constraint or as an expression.
#
# \code    
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    var3 = Variable()
#    
#    model.add(var1 > var2)
#    model.add(var1 > 5)
#
#    model.add( var3 > (var1 > var2) )    
# \endcode
#
class Gt(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "gt")
    
    def get_symbol(self): return '>'

## Less than or equal expression
#
# \note
#   - Top-level: Less or equal than Constraint
#   - Nested: Equal to the truth value of the less than or equal relation between the arguments
#
#    Le expression can be used to create an less than or equal expression between two
#    expressions or an expression and a constraint. It can be used as either
#    a hard constraint or as an expression.
#
# \code    
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    var3 = Variable()
#    
#    model.add(var1 <= var2)
#    model.add(var1 <= 5)
#
#    model.add( var3 <= (var1 <= var2) )    
# \endcode
#
class Le(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "le")
    
    def get_symbol(self): return '<='

## Greater than or equal expression
#
# \note
#   - Top-level: Greater than or equal Constraint
#   - Nested: Equal to the truth value of the greater than or equal relation between the arguments
#
#    Ge expression can be used to create an greater than or equal expression between two
#    expressions or an expression and a constraint. It can be used as either
#    a hard constraint or as an expression.
#
# \code    
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#    var3 = Variable()
#    
#    model.add(var1 >= var2)
#    model.add(var1 >= 5)
#
#    model.add( var3 >= (var1 >= var2) )    
# \endcode
#
class Ge(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "ge")
    
    def get_symbol(self): return '>='
    
## Negate expression
#
# \note
#   - Top-level: Cannot be used as a top level predicate
#   - Nested: Expression representing the negation of its' argument
#
#    Neg expression, it is used to negate another expression. It is equivalent
#    to multiplying by -1.
#
# \code    
#    var = Variable(1, 10)
#    
#    model.add(-var < 3)
# \endcode
#
class Neg(Predicate):

    def __init__(self, vars):
        Predicate.__init__(self, vars, "neg")

    def __str__(self):
        return '-'+str(self.children[0])
        
## Table constraint
#
# \note
#   - Top-level: Table constrint
#   - Nested: Cannot be used as a nested predicate
#
#    Table constraint used to handle list of allowed or forbidden tuples.
#
class Table(Predicate):
    
    ## Table constraint constructor
    # @param vars variables to be constrained by the constraint
    # @param tuples tuples used for the table constraint
    # @param type type of table constraint, e.g. support of conflict
    def __init__(self, vars, tuples=[], type='support'):
        Predicate.__init__(self, vars, "Table")
        self.parameters = [[tuple for tuple in tuples], type]

    ## Adds another support tuple to the list of tuples
    # @param tuple tuple to be added
    #
    #    Adds in a support tuple to the list of tuples if the Table constraint
    #    is of type support. If the table constraint is of type conflict then
    #    in the provided tuple is a conflict tuple in the tuple list it is removed.
    #
    def addSupport(self, tuple):

        if self.parameters[1] == 'support':
            self.parameters[0].append(tuple)
        else:
            self.parameters[0].remove(tuple)

    ## Adds another conflict tuple to the list of tuples
    # @param tuple tuple to be added
    #
    #    Adds in a conflict tuple to the list of tuples if the Table constraint
    #    is of type conflict. If the table constraint is of type support then
    #    in the provided tuple is a support tuple in the tuple list it is removed.
    #
    def addConflict(self, tuple):

        if self.parameters[1] != 'support':
            self.parameters[0].append(tuple)
        else:
            self.parameters[0].remove(tuple)
        
    ## Prints the table of tuples    
    def printTable(self):
        print self.parameters[1]
        for var in self.children:
            var.print_me(),
        print '\n', self.parameters[0]





## Sum constraint
#
# \note
#   - Top-level: Cannot be used as a top level constraint
#   - Nested: Expression representing the sum of the predicates linear coefficients
#
#    Sum Constraint with linear coefficients.
#
class Sum(Predicate):
    
    ## Table constraint constructor
    # @param vars variables to be summed
    # @param coefs list of coefficients ([1,1,..1] by default)
    def __init__(self, vars, coefs = None):
        Predicate.__init__(self, vars, "Sum")
        
        if coefs is None:
            coefs = [1 for var in vars]
    
        self.parameters = [coefs, 0]

    def close(self):

        """
        This handles the scalar constraint, i.e. with weights
        """
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
                        res.extend(extract_sum(s_var, s_coef*coef))
                    if var.parameters[1] != 0:
                        res.append((var.parameters[1], 1))
                    return res
                elif var.operator == "mul":
                    if (var.get_children()[0].operator == "var" and
                        (type(var.get_children()[1]) == types.IntType or
                         type(var.get_children()[1]) == types.FloatType)):
                        return [(var.get_children()[0], var.get_children()[1] * coef)]
                    elif (type(var.get_children()[1]) == types.IntType or
                         type(var.get_children()[1]) == types.FloatType):

                        return [(new_var, new_coef*var.get_children()[1]*coef)
                            for (new_var, new_coef) in extract_sum(var.get_children()[0], 1)]

                    else:
                        return [(var, 1*coef)] # It is quadratic?
                else:
                    return [(var, 1*coef)]
            else:
                return [(var, 1*coef)]

            
        # This is where it should go looking for +s
        set_vars = set([])        
        list_vars = []
        map_coefs = {}
        offset = self.parameters[1]

        for (var, coef) in zip(self.children, self.parameters[0]):
            list = extract_sum(var, coef) 
            for (nVar, nCoef) in list:
                if type(nVar) == types.IntType or type(nVar) == types.FloatType:
                    offset += (nVar*nCoef)
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
                op += (str(self.parameters[0][0])+'*')
            op += (self.children[0].__str__())
        for i in range(1,len(self.children)):
            if self.parameters[0][i] == 1:
                op += (' + '+self.children[i].__str__())
            elif self.parameters[0][i] == -1:
                op += (' - '+self.children[i].__str__())
            elif self.parameters[0][i] > 0:
                op += (' + '+str(self.parameters[0][i])+'*'+self.children[i].__str__())
            elif self.parameters[0][i] < 0:
                op += (' - '+str(-self.parameters[0][i])+'*'+self.children[i].__str__())
        if self.parameters[1] > 0:
            op += (' + '+str(self.parameters[1]))
        elif self.parameters[1] < 0:
            op += (' - '+str(-self.parameters[1]))
        return op+')'


    def decompose(self):
        def addition(X):
            if len(X) == 1: return X[0]
            else: return X[0] + self.addition(X[1:])
        return [addition([child * coef for child, coef in zip(self.children, self.parameters)])]

## All-Different Constraint
#
# \note
#   - Top-level: All Different Constraint
#   - Nested: Cannot be used as a nested predicate
#
#    All-Different Constraint on a list of Expressions
#
class AllDiff(Predicate):

    def __init__(self, vars, type=None):
        Predicate.__init__(self, vars, "AllDiff")
        if type != None:
            self.parameters = [type]
    
    #def __str__(self):
    #    return " AllDiff(" + " ".join(map(str, self.children)) + " ) "

## Global Cardinality Constraint
#
# \note
#   - Top-level: Cardinality Constraint Constraint
#   - Nested: Cannot be used as a nested predicate
#
#    Global Cardinality Constraint on a list of Expressions
#
class Gcc(Predicate):

    ## Global Cardinality constraint constructor
    # @param vars variables 
    # @param cards dictionary: value -> pairs of integers (l,u)
    #
    # The cardinalities connot be Variables. 
    # The parameter cards is a dictionary mapping a set of values
    # to lower and upper bounds. For instance:
    #
    # \code
    #X = VarArray(5,1,4)
    #cards = {1:(1,2), 2:(2,2), 3:(0,3), 4:(1,2)}
    #model = Model( Gcc(X,cards) )
    #print model
    #>>> assign:
    #>>>   x0 in {1..4}
    #>>>   x1 in {1..4}
    #>>>   x2 in {1..4}
    #>>>   x3 in {1..4}
    #>>>   x4 in {1..4}
    #>>>   
    #>>> subject to:
    #>>>   Gcc(x0 x1 x2 x3 x4 | 1 in [1,2] 2 in [2,2] 3 in [0,3] 4 in [1,2] )
    # \endcode
    #
    def __init__(self, vars, cards):
        Predicate.__init__(self, vars, "Gcc")
        values = cards.keys()
        values.sort()
        lb = []
        ub = []
        for val in values:
            lb.append(cards[val][0])
            ub.append(cards[val][1])
        self.parameters = [values, lb, ub]
        
    def __str__(self):
        save_str = Expression.__str__
        Expression.__str__ = Expression.name
        output = " Gcc(" + " ".join(map(str, self.children)) + " | " 
        for (v,l,u) in zip(*(self.parameters)):
            output += str(v)+' in ['+str(l)+','+str(u)+'] '
        Expression.__str__ = save_str
        return output+')'

    def decompose(self):
        X = self.children
        decomp = []
        for val, l, u in zip(self.parameters[0], self.parameters[1], self.parameters[2]):
            card = Variable(l,u)
            decomp.append((card == Cardinality(X,val)))
        return decomp

## Max Constraint
#
# \note
#   - Top-level: Cannot be used as a top level constraint
#   - Nested: Expression representing the maximum of its arguments
#
# Predicate holding the maximum value of a set of Variables
#
class Max(Predicate):

    def __init__(self, vars):
        Predicate.__init__(self, vars, "Max")

    def get_min(self, solver=None):
        return max([x.get_min(solver) for x in self.children])

    def get_max(self, solver=None):
        return max([x.get_max() for x in self.children])

    def decompose(self):
        X = self.children
        M = Variable(self.get_min(), self.get_max(), 'Max')
        decomp = [M]
        decomp.extend([M >= x for x in X])
        decomp.append(Disjunction([M <= x for x in X]))
        return decomp

    #def __str__(self):
    #    return " MAX ( " + " ".join(map(str, self.children)) + " ) "


## Min Constraint
#
# \note
#   - Top-level: Cannot be used as a top level constraint
#   - Nested: Expression representing the minimum of its arguments
#
# Predicate holding the maximum value of a set of Variables
#
class Min(Predicate):
    def __init__(self, vars):
        Predicate.__init__(self, vars, "Min")

    def get_min(self, solver=None):
        return min([x.get_min(solver) for x in self.children])

    def get_max(self, solver=None):
        return  min([x.get_max(solver) for x in self.children])

    def decompose(self):
        X = self.children
        M = Variable(self.get_min(), self.get_max(), 'Min')
        decomp = [M]
        decomp.extend([M <= x for x in X])
        decomp.append(Disjunction([M >= x for x in X]))
        return decomp

    #def __str__(self):
    #    return " MIN ( " + " ".join(map(str, self.children)) + " ) "

## Element Constraint
#
# \note
#   - Top-level: Cannot be used as a top level constraint
#   - Nested: Element constraint
#
# Given an integer Variable \e index and a VarArray \e vars,
# \e Element is the Predicate holding the value of the variable at index \e index of the array \e vars.
#
# \code
#vars = VarArray(5,1,4,'var')
#index = Variable(0,4)
#elt1 = Element(vars, index)
#elt2 = vars[index]
#print elt1
#>>> Element(var0, var1, var2, var3, var4, x)
#print elt2
#>>> Element(var0, var1, var2, var3, var4, x)
# \endcode 
#
class Element(Predicate):
    def __init__(self, vars, index):
        children = list(vars)
        children.append(index)
        Predicate.__init__(self, children, "Element")
        
        
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

    def add(self, literal):
        if literal.operator == 'neg':
            self.parameters[0].append(0)
            self.children.append(literal.children[0])
        else:
            self.parameters[0].append(1)
            self.children.append(literal)

    def __str__(self):
        ret_str = "Clause(";
        for i in range(len(self.children)):
            ret_str += ' '
            if self.parameters[0][i] == 1:
                ret_str += str(self.children[i])
            else:
                ret_str += ('~'+str(self.children[i]))
        return ret_str+' )';
    
## Less Lexicographic Ordering Constraint
class LessLex(Predicate):
    def __init__(self, vars_1, vars_2):
        children = list(vars_1)
        children.extend(vars_2)
        Predicate.__init__(self, children, "LessLex")
        
    def __str__(self):
        length = len(self.children)/2

        toprint = '['+str(self.children[0])
        for i in range(1,length):
            toprint += (', '+str(self.children[i]))
        toprint += '] < ['+str(self.children[length])
        for i in range(length+1,2*length):
            toprint += (', '+str(self.children[i]))
        return toprint+']'

## Leq Lexicographic Ordering Constraint
class LeqLex(Predicate):
    def __init__(self, vars_1, vars_2):
        children = list(vars_1)
        children.extend(vars_2)
        Predicate.__init__(self, children, "LeqLex")
        
    def __str__(self):
        length = len(self.children)/2

        toprint = '['+str(self.children[0])
        for i in range(1,length):
            toprint += (', '+str(self.children[i]))
        toprint += '] <= ['+str(self.children[length])
        for i in range(length+1,2*length):
            toprint += (', '+str(self.children[i]))
        return toprint+']'

## Maximisation objective function
#
# Sets the goal of search to be the maximisation of its' argument
#
class Maximise(Predicate):
    def __init__(self, vars):
        Predicate.__init__(self, [vars], "Maximise")
        
    #def __str__(self):
    #    return " Maximise ( " + " ".join(map(str, self.children)) + " ) "

## Alias for American spelling
def Maximize(var):
    return Maximise(var)

## Minimisation objective function
#
# Sets the goal of search to be the minimisation of its' argument
#
class Minimise(Predicate):
    def __init__(self, vars):
        Predicate.__init__(self, [vars], "Minimise")

    #def __str__(self):
    #    return " Minimise ( " + " ".join(map(str, self.children)) + " ) "

## Alias for American spelling
def Minimize(var):
    return Minimise(var)
## @}


class Disjunction(Predicate):

    def __init__(self, vars):
        Predicate.__init__(self, vars, "OR")

    def decompose(self):
        return [Sum(self.children) > 0]


class Convex(Predicate):
    def __init__(self, vars):
        Predicate.__init__(self, [var for var in vars], "Convex")
        
    def __str__(self):
        return "[" + " ".join(map(str, self.children)) + "] is row-convex"
        
    def decompose(self):
        ### BUGGY!!

        print "Decomposing Row convexity constraint", self
        X = self.children
        n = len(X)
        first = Variable(n)
        last = Variable(n)
        decomposition = [X[i] <= (first <= i) for i in range(n)]
        decomposition.extend([X[i] <= (last >= i) for i in range(n)])
        decomposition.extend([((first <= i) & (i <= last)) <= X[i] for i in range(n)])
        decomposition.append(first <= last)

        print VarArray(decomposition)
        return decomposition

class Cardinality(Predicate):
    def __init__(self, vars, value):
        Predicate.__init__(self, [var for var in vars], "Card")
        self.parameters = [value]

    def __str__(self):
        return "card of "+str(self.parameters[0])+" in [" + " ".join(map(str, self.children)) + "]"
        
    def decompose(self):
        X = self.children
        val = self.parameters[0]
        return [Sum([x == val for x in X])]


class Cmp2(BinPredicate):
    ## expression equal to -1 if x_0 < x_1, 1 if x_1 < x_0 and 0 otherwise
    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "Cmp")

    def get_symbol(self): return 'cmp'
        
    def decompose(self):
        X = self.children
        return [(X[1]<X[0])-(X[0]<X[1])]



## @defgroup sched_group Scheduling constructs
# The scheduling constructs and constraints
#  @{
#

'''
class Job(VarArray):
    def __init__(self, L):
        VarArray.__init__(self, n)
'''   


## Task specific Variable
class Task(Expression):
# arg1: lower bound
# arg2: upper bound
# arg3: duration

# arg1: upper bound
# arg2: duration

    def __init__(self, arg1=None, arg2=None, arg3=None):
        
        lb = 0
        ub = 1
        self.duration = 1
        if arg1 != None:
            if arg2 != None:
                if arg3 != None:
                    lb = arg1
                    ub = arg2-arg3
                    self.duration = arg3
                else: # only 2 args, read as <makespan,duration>
                    ub = arg1-arg2
                    self.duration = arg2
            else: # only 1 arg, read as the makespan
                ub = arg2-1
        Expression.__init__(self, "t")
        self.lb = lb
        self.ub = ub
        self.domain_ = None

    def __str__(self):
        if self.is_built() and self.solver.is_sat():
            return str(self.get_value())
        else:
            ident = str(self.ident)
            if self.ident == -1: ident = ''
            return 't'+str(ident)+': ['+str(self.get_min())+':'+str(self.duration)+':'+str(self.get_max()+self.duration)+']'
        
    def __lt__(self, pred):
        if type(pred) is int:
            return Le([self, pred-self.duration])
        return Precedence(self, pred, self.duration)
    
    def __gt__(self, pred):
        if type(pred) is int:
            return Gt([self, pred])
        return pred.__lt__(self)

    def requires(self, resource):
        resource.add(self)

    def reset(self, makespan):
        self.ub = makespan-self.duration

## Precedence Constraint
class Precedence(Predicate):
    def __init__(self, task_i, task_j, dur):
        Predicate.__init__(self, [task_i, task_j], "Precedence")
        self.parameters = [dur]

    def decompose(self):
        return [ ((self.children[0] + self.parameters[0]) <= self.children[1]) ]

    def __str__(self):
        return str(self.children[0])+' + '+str(self.parameters[0])+' <= '+str(self.children[1])


## Binary disjunctive Constraint
class NoOverlap(Predicate):

    def __init__(self, task_i, task_j, dur_i=None, dur_j=None):
        if dur_i == None: dur_i = task_i.duration
        if dur_j == None: dur_j = task_j.duration
        Predicate.__init__(self, [task_i, task_j], "NoOverlap")
        self.parameters = [dur_i,dur_j]
        
    def decompose(self):
        return [ ( ((self.children[0] + self.parameters[0]) <= (self.children[1])) |
                   ((self.children[1] + self.parameters[1]) <= (self.children[0])) )]

    def __str__(self):
        return str(self.children[0])+' + '+str(self.parameters[0])+' <= '+str(self.children[1])+' OR '+str(self.children[1])+' + '+str(self.parameters[1])+' <= '+str(self.children[0])

'''
## Not sure what that is?
class Interval(list):
    
    def __init__(self, lhs, rhs):
        list.__init__(self)
        self.lhs = lhs
        self.rhs = rhs    
        
    def __contains__(self, y):
        if type(y) is not types.IntType and type(y) is not types.FloatType:
            constraint = (y >= self.lhs) & (y <= self.rhs)
            self.append(constraint)
            return constraint
        else:
            constraint = (self.lhs <= y) & (self.rhs >= y)
            self.append(constraint)
            return constraint
        
    def isin(self, y):
        return self.__contains__(y)
'''

## Unary resource constraint
class UnaryResource(list):
    
    def __init__(self, arg=[], distance=0):
        list.__init__(self, [NoOverlap(arg[i], arg[j], arg[i].duration+distance, arg[j].duration+distance) for i in range(1,len(arg)) for j in range(i)])
        self.tasks = [task for task in arg]
        self.distance = distance

    def add(self, new_task):
        for task in self.tasks:
            self.append( NoOverlap(new_task, task, new_task.duration+self.distance, task.duration+self.distance) )            
        self.tasks.append(new_task)
        
    def __str__(self):
        return "[" + " ".join(map(str, self.tasks)) + "] share a unary resource" #: "+" ".join(map(str, self))
        
class UnaryResourceB(Predicate):
    def __init__(self, tasks, distance=0):
        Predicate.__init__(self, [task for task in tasks], "UnaryResource")
        self.distance = distance
        
    def add(self, new_task):
        self.children.append(new_task)
        
    def __str__(self):
        return "[" + " ".join(map(str, self.children)) + "] share a unary resource"
        
    def decompose(self):
        return [NoOverlap(task1, task2) #, task1.duration+self.distance, task2.duration+self.distance)
                for task1, task2 in pair_of(self.children)]
        
## @}

class ParamList(dict):
    def __init__(self, X):
        dict.__init__(self,X)

    def __call__(self, *args):
        return [self.__getitem__(arg) for arg in args]


def input(default):
    import sys

    param_list = ParamList(default)
    option = None
    params = []

    commandline = [arg for arg in sys.argv[1:]]
    commandline.append('-end_argument')

    for arg in commandline:
        #print arg

        if arg[0] == '-' and arg != '-1.0': # new argument
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
                    param_list[option] = params;

            option = arg[1:]
            #print 'new option', option

            if option != 'end_argument':
                if not param_list.has_key(option):
                    
                    #print 'unknwn option'
                    if option == 'h' or option == '-h' or option == 'help' or option == '-help':
                        #print 'help'
                        the_keys = param_list.keys()
                        the_keys.sort()
                        for key in the_keys:
                            print ('-'+key).ljust(20)+":",
                            if type(param_list[key]) == int:
                                print 'int'.ljust(14),
                            if type(param_list[key]) == float:
                                print 'float'.ljust(14),
                            elif type(param_list[key]) == str:
                                print 'string'.ljust(14),
                            elif hasattr(param_list[key],'__iter__'):
                                print 'list',
                                if len(param_list[key]) > 0:
                                    if type(param_list[key][0]) == int:
                                        print 'of int   ',
                                    elif type(param_list[key][0]) == float:
                                        print 'of float ',
                                    else:
                                        print 'of string',
                                else:
                                    print 'of string',
                            print ' (default='+str(param_list[key])+')'
                        exit(1)
                    else:
                        print 'Warning: wrong parameter name, ignoring arguments following', option
                else:
                    #print 'fine option init param list'
                    params = []
        else:
            params.append(arg)
            #print 'input param:', params

    return param_list
                

def pair_of(l):
    return [pair for k in range(1,len(l)) for pair in zip(l, l[k:])]
    
def value(x): return x.get_value()

def load_in_decompositions():
    import Decomp
    
    # First add the constraints
    for attr_name in dir(Decomp):
        attr = getattr(Decomp, attr_name)
        if hasattr(attr, "__name__") and attr.__name__.find("decompose") is -1:
            if not hasattr(sys.modules[__name__], attr.__name__):
                setattr(sys.modules[__name__], attr.__name__, attr)
                
    # Now load up the decompositions, introspection is great huh :)
    for attr_name in dir(sys.modules[__name__]):
        if hasattr(Decomp, "decompose_"+attr_name):
            setattr(getattr(sys.modules[__name__], attr_name),
                    "decompose", getattr(Decomp, "decompose_"+attr_name) )

load_in_decompositions()


class Solution(list):
    def __init__(self,vars):
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
            return list.__getitem__(self,i)
        else:
            return self.dico[i]

    def __contains__(self,x):
        return self.dico.has_key(x)
    
    def __str__(self):
        if len(self) == 0: return '[]'
        elif type(self[0]) is list:
            return '['+',\n '.join([str(row) for row in self])+']'
        else: return list.__str__(self)
    #    output = 
    #    #output = str(self.variables[0])+' = '+str(self[self.variables[0]])
    #    #for var in self.variables[1:]:
    #    #    output += ', '+str(var)+' = '+str(self[var])
    #    #return output

class Nogood(object):
    def __init__(self, clause, solver):
        self.clause = clause
        self.solver = solver

    def __str__(self):
        return self.solver.print_clause(self.clause)


## @defgroup sol_group Solving methods
# The Solvers
#  @{
#

## Generic Solver Class
class NBJ_STD_Solver(object):
    def __init__(self, Library, Wrapper, model=None, X=None, FD=False,
                 clause_limit=-1):
        
        

        self.decomposition_store = []

        self.solver = getattr(sys.modules[Library],
                              Library+"Solver", None)()
        
        if hasattr(self.solver, "setClauseLimit"):
            self.solver.setClauseLimit(clause_limit);
        
        self.Library = Library
        
        self.Wrapper = Wrapper

        """
        self.ExpArray = getattr(sys.modules[Library],
                                Library+"ExpArray", None)
        
        self.IntArray = getattr(sys.modules[Library],
                                Library+"IntArray", None)
        
        self.DoubleArray = getattr(sys.modules[Library],
                                   Library+"DoubleArray", None)
        
        self.IntVar = getattr(sys.modules[Library],
                              Library+"_IntVar", None)
        
        self.FloatVar = getattr(sys.modules[Library],
                                Library+"_FloatVar", None)
        """
           
        self.ExpArray = getattr(sys.modules[Wrapper],
                                Wrapper+"ExpArray", None)
        
        self.IntArray = getattr(sys.modules[Wrapper],
                                Wrapper+"IntArray", None)
        
        self.DoubleArray = getattr(sys.modules[Wrapper],
                                   Wrapper+"DoubleArray", None)
        
        self.IntVar = getattr(sys.modules[Wrapper],
                              Wrapper+"_IntVar", None)
        
        self.FloatVar = getattr(sys.modules[Wrapper],
                                Wrapper+"_FloatVar", None)

                     
        self.free_memory = getattr(sys.modules["_"+Library],
                                   "delete_"+Library+"Solver", None)
        
        self.variables = None
        
        self.solver_id = model.getSolverId()

        if model != None: 
            var_array = None
            self.model = model
            self.model.close()
            for expr in self.model.get_exprs():
                self.solver.add(self.load_expr(expr))

            if X != None:

                self.variables = VarArray(flatten(X))

                var_array = self.ExpArray()
                for x in self.variables:
                    if len(x.var_list) > self.solver_id-1:
                        var_array.add(x.var_list[self.solver_id-1])
                if FD:
                    self.solver.forceFiniteDomain(var_array)
                if var_array.size() > 0:
                    self.solver.initialise(var_array)
            else:
                self.variables = self.model.variables
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
                var = self.IntVar(arg1,arg2)
            else:
                var = self.IntVar(arg1,arg2,argopt1)
        except:
            print "ERROR while creating variable"
            sys.exit(1)
        return var

    def getFloatVar(self, lb, ub, ident):
        var = None
        try:
            var = self.FloatVar(lb,ub,ident)
        except:
            print "ERROR: Solver does not support real variables"
            sys.exit(1)
        return var

    def load_expr(self, expr):
        #print 'load', expr
        if type(expr) is str:
            return self.model.string_map[expr]
        if type(expr) is bool:
            return int(expr)
        if type(expr) is int or type(expr) is float:
            # It is a constant, handle appropriatly
            return expr
        #if not expr.has_children():
        if expr.is_var():
            # It is a leaf
            if expr.is_var(): # Just to be sure
                alreadyBuild = False
                if expr.is_built(self):
                    #if expr.get_solver() == self: # Checks if I have already assigned
                    alreadyBuild = True

                if alreadyBuild:
                    return expr.var_list[self.solver_id-1]
                else:
                    # It is probably a variable
                    (lb,ub,domain) = expr.get_domain_tuple()

                    var = None
                    if domain is None:
                        if type(lb) is int:
                            var = self.getIntVar(lb, ub, expr.ident)
                        else:
                            var = self.getFloatVar(lb, ub, expr.ident)
                    elif ub-lb+1 == len(domain):
                        if type(lb) is int:
                            var = self.getIntVar(lb, ub, expr.ident)
                        else:
                            var = self.getFloatVar(lb, ub, expr.ident)
                    else:
                        w_array = self.IntArray()
                        for val in domain:
                            w_array.add(val)
                        var = self.getIntVar(w_array, expr.ident)

                    expr.setVar(self.solver_id, self.Library, var, self)
                    expr.solver = self
                    return var
            else:
                raise Error("Problem, no such type exists in converting models")
        else:
#            factory = getattr(sys.modules[self.Library],
#                              "%s_%s" % (self.Library, expr.get_operator()), None)
            factory = getattr(sys.modules[self.Wrapper],
                              "%s_%s" % (self.Wrapper, expr.get_operator()), None)

            if factory is not None:
            
                arguments = None
                if len(expr.get_children()) <= 2:
                    arguments = [self.load_expr(child) for child in expr.get_children()]
                else:
                    var_array = self.ExpArray()
                    for child in expr.get_children():
                        var_array.add(self.load_expr(child))                
                    arguments = [var_array]
                if expr.has_parameters(): # != None: # assumes an array of integers
                    for param in expr.parameters:
                        if hasattr(param, '__iter__'):
                            if expr.get_operator() is "Table":
                                w_array = self.IntArray()
                                for w in param:
                                    for v in w:
                                        w_array.add(v)
                                arguments.append(w_array)    
                            else:
                                w_array = self.IntArray()
                                for w in param:
                                    w_array.add(w)
                                arguments.append(w_array)
                        else:
                            arguments.append(param)
    
                var = factory(*arguments)
                expr.setVar(self.solver_id, self.Library, var, self)
                expr.solver = self
                return var
            else:
                return self.decompose_expression(expr)

    def decompose_expression(self, expr):
        if hasattr(expr, "decompose"):
            expr_list = expr.decompose()
            obj_exp = []
            for exp in expr_list[1:]: 
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
            print "ERROR: Solver does not support \"%s\" and no decomposition is available" % (expr.get_operator())
            exit(1)
            raise ConstraintNotSupportedError(expr.get_operator(), self.Library)

        
    ##@name Solving methods
    # @{

    ## Solves the extracted Model
    def solve(self):
        try:
            if self.solver.solve() == SAT:
                return True
            return False
        except (KeyboardInterrupt, SystemExit):
            print 'Program Interrupted'
            return

    ## Solves using restarts
    def solveAndRestart(self, policy=GEOMETRIC, base=64, factor=1.3, decay=0.0, reinit=-1):
        if reinit == -1:
            if self.solver.solveAndRestart(policy,base,factor,decay) == SAT:
                return True
        else:
            if self.solver.solveAndRestart(policy,base,factor,decay,reinit) == SAT:
                return True
        return False

    ## Initialise structures for a depth first search
    def startNewSearch(self):
        self.solver.startNewSearch()

    ## Search for the next solution
    def getNextSolution(self):
        return self.solver.getNextSolution()

    ## @}

    def next(self, exp, v):
        #return self.solver.next(exp.var_list[self.solver_id-1], v)
        return exp.var_list[self.solver_id-1].next(v)

    ##@name Search programming methods
    # @{

    ## Tell the solver to reach a fixed point
    #@return False if an inconsistency is found, True otherwise
    def propagate(self):
        return self.solver.propagate()

    ## Tell the solver to save the current state. 
    def save(self):
        self.solver.save()

    ## Tell the solver to restore its state to the one last enqueued by 'save()'
    def undo(self,bj=1):
        return self.solver.undo(bj)

    ## Tell the solver to add a constraint in the current state (it should be unary for now)
    def post(self, exp):
        self.solver.post(exp.operator, exp.children[0].var_list[self.solver_id-1], exp.children[1])

    ## Tell the solver to post the negation of the last decision in the current states
    def deduce(self, exp=None):
        if exp is None: self.solver.deduce()
        else: self.solver.post(exp.operator, exp.children[0].var_list[self.solver_id-1], exp.children[1])
    def deduce_print(self,lvl):
        x = self.variables[self.solver.get_decision_id()]
        print lvl*' ', x.domain(self)
        self.solver.deduce()


    ## Follow a 'left' branch in a binary serach tree (eq. to save() + post())
    def branch_left(self,exp):
        self.solver.save()
        self.solver.post(exp.operator, exp.children[0].var_list[self.solver_id-1], exp.children[1])

    ## Follow a 'right' branch in a binary serach tree (eq. to undo() + deduce())
    def branch_right(self):
        return self.solver.branch_right()

    ## Resets the data structure of the solver to the initial state 
    # @param full Boolean stating if the top-level deduction should be undone too
    def reset(self, full=False):
        self.solver.reset(full)

    ## @}

    def analyze_conflict(self):
        btlevel = self.solver.analyze_conflict()
        if btlevel < 0: btlevel = None
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
            
        print '(',
        for i in range(self.solver.get_nogood_size()):
            if i > 0: print 'or',
            print get_literal(i),
        print ')'




    def sacPreprocess(self,type):
        self.solver.sacPreprocess(type)
        return (not (self.solver.is_unsat()))

    ##@name Parameter tuning methods
    # @{

    ## Sets the variable and value heuristics
    def setHeuristic(self, arg1, arg2='No', arg3=0):
        var_name = arg1
        val_name = arg2
        randomization = arg3
        if type(val_name) is int:
            randomization = val_name
            val_name = 'No'
        if var_name not in var_heuristics:
            print 'c Warning: "'+var_name+'" unknown, use MinDomain instead' 
            print 'c legal variable orderings: ', var_heuristics
        if val_name not in val_heuristics:
            print 'c Warning: "'+val_name+'" unknown, use Lex instead' 
            print 'c legal value orderings: ', val_heuristics
        self.solver.setHeuristic(str(var_name), str(val_name), randomization)
   
    ## Sets a limit on the number of failures encountered before aborting search
    def setFailureLimit(self,cutoff):  
        self.solver.setFailureLimit(cutoff)

    ## Sets a limit on the CPU time before aborting search
    def setTimeLimit(self,cutoff):
        self.solver.setTimeLimit(cutoff)

    ## Sets a limit on the number of nodes explored before aborting search
    def setNodeLimit(self,cutoff):
        self.solver.setNodeLimit(cutoff)

    ## Sets the verbosity level of the solver
    def setVerbosity(self, degree):
        self.solver.setVerbosity( degree )

    ## Sets the initial random seed
    def setRandomSeed(self, seed):
        self.solver.setRandomSeed(seed);

    ## @}
    
    def setRandomized(self, degree):
        self.solver.setRandomized(degree)

    def addNogood(self, vars, vals):
        var_array = self.ExpArray()
        for var in vars:
            var_array.add(var.var_list[self.solver_id-1])
        val_array = self.IntArray()
        for val in vals:
            val_array.add(val)
        self.solver.addNogood(var_array, val_array)
        
    def setAntiLex(self, vars):
        var_array = self.ExpArray()
        for var in vars:
            var_array.add(var.var_list[self.solver_id-1])
        self.solver.setAntiLex(var_array)

    def guide(self, vars, vals=None, probs=[]):
        
        var_array = self.ExpArray()
        val_array = self.IntArray()
        pro_array = self.DoubleArray()

        if vals is not None:
            for var in vars:
                var_array.add(var.var_list[self.solver_id-1])
            for val in vals:
                val_array.add(val)
        else:
            for var in vars.variables:
                var_array.add(var.var_list[self.solver_id-1])
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
        tmp_array = []
        var_array = self.ExpArray()
        lob_array = self.IntArray()
        for (x,l) in zip(vars,lb):
            var_array.add(x.var_list[self.solver_id-1])
            lob_array.add(l)
        self.solver.setLowerBounds(var_array, lob_array)
    def setUpperBounds(self, vars, ub):
        tmp_array = []
        var_array = self.ExpArray()
        upb_array = self.IntArray()
        for (x,u) in zip(vars,ub):
            var_array.add(x.var_list[self.solver_id-1])
            upb_array.add(u)
        self.solver.setUpperBounds(var_array, upb_array)
    def setRestartNogood(self):
        self.solver.setRestartNogood()


    ##@name Accessors
    # @{

    ## Extract a Solution object from the solver
    def get_solution(self):
        self.solver.store_solution()
        solution = Solution(self.variables)
        return solution

    ## Returns True iff the solver found a solution and proved its optimality
    def is_opt(self):
        return self.solver.is_opt()

    ## Returns True iff the solver found a solution
    def is_sat(self):
        return self.solver.is_sat()

    ## Returns True iff the solver proved unsatisfiability
    def is_unsat(self):
        return self.solver.is_unsat()

    ## Returns the number of backtracks performed during the last search
    def getBacktracks(self):
        return self.solver.getBacktracks()

    ## Returns the number of nodes explored during the last search
    def getNodes(self):
        return self.solver.getNodes()
    
    ## Returns the number of failures encountered during the last search
    def getFailures(self):
        return self.solver.getFailures()

    ## Returns the number of constraint propagations performed during the last search
    def getPropags(self):
        return self.solver.getPropags()

    ## Returns the CPU time required for the last search
    def getTime(self):
        return self.solver.getTime()

    ## @}

    def getChecks(self):
        return self.solver.getChecks()

    def printStatistics(self):
        print ''
        self.solver.printStatistics()
        print ''

    def getNumVariables(self):
        return self.solver.getNumVariables()

    def getNumConstraints(self):
        return self.solver.getNumConstraints()

    def load_xml(self,file,type=4):
        self.solver.load_xml(file,type)
    
    def num_vars(self):
        return self.solver.num_vars()

    def extract_graph(self):
        self.solver.extract_graph()

    def numNodes(self):
        return self.solver.numNodes()

    def degree(self, var):
        return self.solver.get_degree(var)

    def get_neighbors(self, x):
        neighbors = []
        for y in range(self.solver.degree(x)):
            neighbors.append(self.solver.get_neighbor(x,y))
        return neighbors
            
    def get_static_features(self):
        feats = {}
        for i in range(12):
            feats[self.solver.get_feature_name(i)] = self.solver.get_feature(i)
        for i in range(16,36):
            feats[self.solver.get_feature_name(i)] = self.solver.get_feature(i)
        return feats
    
    def get_dynamic_features(self):
        feats = {}
        for i in range(12,16):
            feats[self.solver.get_feature_name(i)] = self.solver.get_feature(i)
        return feats

    def __str__(self):
        #setActive(self)
        #return ''
        self.solver.printPython()
        return ' '

    def delete(self):
        self.free_memory( self.solver )

##  @}


"""

Numberjack exceptions:

"""
import exceptions

class ConstraintNotSupportedError(exceptions.Exception):
    def __init__(self, value, solver):
        self.value = value
        self.solver = solver
        
    def __str__(self):
        return "ERROR: Constraint %s not supported by solver %s" % (self.value, self.solver)

