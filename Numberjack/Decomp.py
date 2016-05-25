from __future__ import division, print_function

MODELMAXSIZE = 1000000

from Numberjack import *

'''
Custom constraints section
'''

class MyConstraint(Expression):

    def __init__(self, vars):
        Expression.__init__(self, "MyConstraint")

        self.set_children(vars)

        print("This is the MyConstraint method")

    def decompose(self):
        '''
        Decompose must return either a list containing a list of expressions
        '''

        constraint_list = []
        variable = Variable(0, 1)

        return (variable, constraint_list)

class MyAllDiff(Expression):

    def __init__(self, vars):
        Expression.__init__(self, "MyAllDiff")
        self.set_children(vars)

    def decompose(self):
        return [var1 != var2 for var1, var2 in pair_of(self.children)]

class MyAddTwo(Expression):

    def __init__(self, vars):
        Expression.__init__(self, "MyAddTwo")
        self.set_children(vars)

    def decompose(self):
        return [self.children[0] + self.children[1]]

class MyAddThree(Expression):
    def __init__(self, vars):
        Expression.__init__(self, "MyAddThree")
        self.set_children(vars)

class MySum(Expression):

    def __init__(self, vars):
        Expression.__init__(self, "MySum")
        self.set_children(vars)

    def addition(self,X):
        if len(X) == 1: return X[0]
        else: return X[0] + self.addition(X[1:])

    def decompose(self):
        return [self.addition(self.children)]


# SDG: new BinPredicate for decomposing Sum into additions
## Add expression
#
# \note
#   - Top-level: Can not be used as top-level Constraint
#   - Nested: Equal to the sum of the operands
#
#    Add expression can be used when decomposing a Sum.
#
# \code
#    var1 = Variable(0, 10)
#    var2 = Variable(0, 100)
#
#    addexp = Add(var1,var2)
# \endcode
#
class Add(BinPredicate):

    def __init__(self, vars):
        BinPredicate.__init__(self, vars, "add")
        self.lb = self.get_lb(0) + self.get_lb(1)
        self.ub = self.get_ub(0) + self.get_ub(1)

    def get_symbol(self):
        return '+'

# SDG: new Predicate for expressing arbitrary functions
## Function
#
# \note
#   - Top-level: Cannot be used as a top level constraint
#   - Nested: Expression representing an arbitrary function
#
#    Function predicate expressed by a dictionary.
#
class Function(Predicate):

    ## Function predicate constructor
    # @param vars variables involved by the predicate
    # @param dictionary a dictionary of tuples (keys) with their associated results (set of "(value_var1,value_var2,...):result")
    # @param defval default result value (0 by default)
    def __init__(self, vars, dictionary={}, defval=0):
        Predicate.__init__(self, vars, "Function")
        d = []
        self.lb = defval
        self.ub = defval
        for t,res in dictionary.items():
            d.append((t,res))
            if (res < self.lb):
                self.lb = res
            if (res > self.ub):
                self.ub = res
        self.parameters = [dict(d), defval]

    def decompose(self):
        #SDG: Warning! It assumes the default value is never used (dict entries = cartesian product of self.children)
        res = Variable(min(e for e in self.parameters[0].values()), max(e for e in self.parameters[0].values()))
        tuples = []
        for key,val in self.parameters[0].items():
            tuples.append(tuple([val] + list(key)))
        return [res, Table([res] + self.children, tuples, 'support')]

# SDG: new Predicate for expressing cost functions
## Base class of cost functions
#
class CostFunction(Predicate):
    def __init__(self, vars, operator):
        Predicate.__init__(self, vars, operator)
        self.lb = None
        self.ub = None

    def __str__(self):
        return super(CostFunction, self).__str__() + str(self.parameters)


## PostNullary Constraint
#
# \note
#   - Top-level: PostNullary Constraint
#   - Nested: Cannot be used as a nested predicate
#
class PostNullary(CostFunction):
    def __init__(self, cost):
        CostFunction.__init__(self, [], "PostNullary")
        self.parameters = [cost]

    def decompose(self):
        return [Minimize(cost)]


## PostUnary Constraint
#
# \note
#   - Top-level: PostUnary Constraint
#   - Nested: Cannot be used as a nested predicate
#
class PostUnary(CostFunction):
    def __init__(self, var, costs):
        CostFunction.__init__(self, [var], "PostUnary")
        self.parameters = [costs]

    def decompose(self):
        obj = Variable(0,max(self.parameters[0]))
        return [Table([obj] + self.get_children(), [(w,i) for i,w in enumerate(self.parameters[0])], 'support'), Minimize(obj)]

## PostBinary Constraint
#
# \note
#   - Top-level: PostBinary Constraint
#   - Nested: Cannot be used as a nested predicate
#
# \code
#   var1 = Variable(1,3)
#   var2 = Variable(1,3)
#   post = PostBinary(var1,var2,[5,3,2,2,2,2,2,1,2])
#   model = Model(post)
#   print model
#>>> assign:
#>>>   x0 in {1..3}
#>>>   x1 in {1..3}
#>>>
#>>> subject to:
#>>>   PostBinary(x0, x1)
# \endcode
#
class PostBinary(CostFunction):

    ## PostBinary constraint constructor
    # @param var1 first variable
    # @param var2 second variable
    # @param costs table of costs
    #
    # \note
    #   - Top-level: PostBinary Constraint
    #   - Nested: Cannot be used as a nested predicate
    #
    def __init__(self, var1, var2, costs):
        CostFunction.__init__(self, [var1, var2], "PostBinary")
        self.parameters = [costs]

    def decompose(self):
        obj = Variable(0,max(self.parameters[0]))
        var2size = self.get_children()[1].get_size()
        return [Table([obj] + self.get_children(), [(w, i // var2size , i % var2size) for i,w in enumerate(self.parameters[0])], 'support'), Minimize(obj)]

## PostTernary Constraint
#
# \note
#   - Top-level: PostTernary Constraint
#   - Nested: Cannot be used as a nested predicate
#
# \code
#   var1 = Variable(1,2)
#   var2 = Variable(1,2)
#   var3 = Variable(1,2)
#   post = PostTernary(var1,var2,var3,[5,3,2,2,2,2,2,1])
#   model = Model(post)
#   print model
#>>> assign:
#>>>   x0 in {1..2}
#>>>   x1 in {1..2}
#>>>   x2 in {1..2}
#>>>
#>>> subject to:
#>>>   PostTernary(x0, x1, x2)
# \endcode
#
class PostTernary(CostFunction):

    ## PostTernary constraint constructor
    # @param var1 first variable
    # @param var2 second variable
    # @param var3 third variable
    # @param costs table of costs
    #
    # \note
    #   - Top-level: PostTernary Constraint
    #   - Nested: Cannot be used as a nested predicate
    #
    def __init__(self, var1, var2, var3, costs):
        CostFunction.__init__(self, [var1, var2, var3], "PostTernary")
        self.parameters = [costs]

    def decompose(self):
        obj = Variable(0,max(self.parameters[0]))
        var2size = self.get_children()[1].get_size()
        var3size = self.get_children()[2].get_size()
        return [Table([obj] + self.get_children(), [(w, i // var2size // var3size , i // var3size, i % var3size) for i,w in enumerate(self.parameters[0])], 'support'), Minimize(obj)]

## PostNary Constraint
# @param vars list of variables
# @param arity number of variables in the list
# @param default_cost cost used by default if an assignment
#
# \note
#   - Top-level: PostNary Constraint
#   - Nested: Cannot be used as a nested predicate
#
# \code
#   var1 = Variable(2)
#   var2 = Variable(2)
#   var3 = Variable(2)
#   var4 = Variable(2)
#   post = PostNary([var1,var2,var3,var4],4,0)
#   for i in range(2):
#     for j in range(2):
#       for k in range(2):
#         for l in range(2):
#           post.add([i,j,k,l], reduce(op.xor, [i,j,k,l])) # simulates a soft xor constraint on four variables
#   model = Model(post)
#   print model
#>>> assign:
#>>>   x0 in {0,1}
#>>>   x1 in {0,1}
#>>>   x2 in {0,1}
#>>>   x3 in {0,1}
#>>>   
#>>> subject to:
#>>>   PostNary(x0, x1, x2, x3)[4, 0, [[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 1], [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 1, 0], [0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, 0], [1, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 1], [1, 1, 1, 0], [1, 1, 1, 1]], [0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0]]
# \endcode
class PostNary(CostFunction):
    def __init__(self, vars, arity, default_cost):
        CostFunction.__init__(self, vars, "PostNary")
        self.parameters = [arity, default_cost,[],[]]

    def add(self, tupleIndex, cost):
            self.parameters[2].append(tupleIndex)
            self.parameters[3].append(cost)

## PostWSum Constraint
#
# @param vars variables array
# @param arity constraint arity
# @param semantics semantic constraint
# @param baseCost baseCost constraint
# @param comparator comparator constraint
# @param rightRes right result
#
# \note
#   - Top-level: PostWSum Constraint
#   - Nested: Cannot be used as a nested predicate
#
# \code
# vars = VarArray(5,1,4)
# post = PostWSum(vars,5,'hard','1000','==',5)
# model = Model(post)
# \endcode
class PostWSum(CostFunction):
    def __init__(self, vars, arity, semantics, baseCost, comparator, rightRes):
        CostFunction.__init__(self, vars, "PostWSum")
        self.parameters = [arity, semantics, baseCost, comparator, rightRes]

## PostWVarSum Constraint
#
# \note
#   - Top-level: PostWSum Constraint
#   - Nested: Cannot be used as a nested predicate
#
class PostWVarSum(CostFunction):
    def __init__(self, vars, arity, semantics, baseCost, comparator, rightVar):
        CostFunction.__init__(self, vars, "PostWVarSum")
        self.parameters = [arity, semantics, baseCost, comparator]
        self.children.append(rightVar)

## PostWAmong Constraint
#
# \note
#   - Top-level: PostWAmong Constraint
#   - Nested: Cannot be used as a nested predicate
#
class PostWAmong(CostFunction):
    def __init__(self, vars, arity, semantics, baseCost, rightVar = None):
        CostFunction.__init__(self, [vars], "PostWAmong")
        self.parameters = [arity, semantics, baseCost, []]
        if rightVar != None :
            self.children.append(rightVar)
        else :
            self.rightvar = None

    def addValues(self, values):
        [self.parameters[3].append(value) for value in values]

    def setBounds (self, lb, ub):
        if self.rightvar == None :
            self.parameters.append(lb)
            self.parameters.append(ub)
        else :
            print("You're using WVarAmong constraint, bounds are set by the Variable passed in the last parameter of the constraint constructor!")

## Regular Constraint
#
# \note
#   - Top-level: Regular Constraint
#   - Nested: Cannot be used as a nested predicate
#
class Regular(CostFunction):
    def __init__(self, vars, arity, nbStates, type = None, measureCost = None):
        CostFunction.__init__(self, vars, "Regular")
        self.parameters = [arity, nbStates,[],[],[]]
        if type != None:
            self.parameters.append(type)
            self.parameters.append(measureCost)
        else:
            self.parameters.append([])
            self.parameters.append([])
            self.parameters.append([])

    def initialStates(self, state, cost=None):
        self.parameters[2].append(state)
        if cost != None:
            self.parameters[5].append(cost)

    def acceptingStates(self, state, cost=None):
        self.parameters[3].append(state)
        if cost != None:
            self.parameters[6].append(cost)

    def transitions(self, start, symbol, end, cost=None):
        self.parameters[4].append([start, symbol, end])
        if cost != None:
            self.parameters[7].append(cost)

## Same Constraint
#
# \note
#   - Top-level: Same Constraint
#   - Nested: Cannot be used as a nested predicate
#
#    Same Constraint on two lists of Expressions (with same size)
#
class Same(CostFunction):

    def __init__(self, varsL, varsR, type=None, semantics=None, baseCost=None):
        CostFunction.__init__(self, [varsL, varsR], "Same")
        if type != None:
            self.parameters = [type, semantics]
            if baseCost != None:
                self.parameters.append(baseCost)

## PostWSameGcc Constraint
#
# \note
#   - Top-level: PostWSameGcc Constraint
#   - Nested: Cannot be used as a nested predicate
#
class PostWSameGcc(CostFunction):
    def __init__(self, varsL, varsR, cards, type, semantics, baseCost):
        CostFunction.__init__(self, [varsL, varsR], "PostWSameGcc")
        values = list(cards.keys())
        values.sort()
        lb = []
        ub = []
        for val in values:
            lb.append(cards[val][0])
            ub.append(cards[val][1])
        self.parameters = [values, lb, ub, type, semantics, baseCost]

## PostWOverlap Constraint
#
# \note
#   - Top-level: PostWOverlap Constraint
#   - Nested: Cannot be used as a nested predicate
#
class PostWOverlap(CostFunction):
    def __init__(self, vars, arity, semantics, baseCost, comparator, rightRes):
        CostFunction.__init__(self, vars, "PostWOverlap")
        self.parameters = [arity, semantics, baseCost, comparator, rightRes]


'''
Custom decomposition section
'''
def decompose_AllDiff(self):
    return [var1 != var2 for var1, var2 in pair_of(self.children)]

#SDG: decompostion of Element constraint using disjunctions (not Berge acyclic!)
def decompose_Element(self):
    u = set()
    for e in self.children[:-1]:
        u = u | set([e] if type(e) is int else list(range(e.lb, e.ub + 1)))
    res = Variable(list(u))
    # ret = ([res, (self.children[-1] >= 0), (self.children[-1] < len(self.children)-1)] + [((res == e) | (self.children[-1] != i)) for i, e in enumerate(self.children[:-1])])
    # print "Decomposed:"
    # for x in map(str, ret):
    #     print x
    # return ret
    return ([res, (self.children[-1] >= 0), (self.children[-1] < len(self.children)-1)] + [((res == e) | (self.children[-1] != i)) for i, e in enumerate(self.children[:-1])])

#SDG: decompose LeqLex and LessLex
def decompose_LessLex(self):
    length = len(self.children) // 2
    def lexico(vars1,vars2,i):
        if (i == len(vars1)-1):
            return (vars1[i] < vars2[i])
        else:
            return (vars1[i] < vars2[i] | And([(vars1[i] == vars2[i]), lexico(vars1,vars2,i+1)]))
    return [lexico(self.children[:length], self.children[length:], 0)]

def decompose_LeqLex(self):
    length = len(self.children) // 2
    def lexico(vars1,vars2,i):
        if (i == len(vars1)-1):
            return (vars1[i] <= vars2[i])
        else:
            return (vars1[i] <= vars2[i] | And([(vars1[i] == vars2[i]), lexico(vars1,vars2,i+1)]))
    return [lexico(self.children[:length], self.children[length:], 0)]

#SDG: util functions for automatic decomposition
def cartesian_product(expr):
    def prod( iterable ):
        p= 1
        for n in iterable:
            p *= n
        return p

    res = prod(e.get_size() for e in get_scope(expr))
    if (res > MODELMAXSIZE):
        raise ModelSizeError(res)
    return res

def get_arity(expr):
    return len(get_scope(expr))

def get_scope(expr):
    if type(expr) in [int, int, float, bool, str]:
        return []
    elif expr.is_var():
        if expr.get_lb() == expr.get_ub():
            return []
        else:
            return [expr]
    else:
        res = []
        for e in expr.children:
            res.extend(get_scope(e))
        return list(set(res))

def evaluate(expr, assignment):
    if type(expr) in [int, int, float, bool, str]:
        return expr
    elif expr.is_var():
        if (expr.get_lb()==expr.get_ub()):
            return expr.get_lb()
        else:
            return assignment[expr]
    elif issubclass(type(expr), BinPredicate):
        return expr.eval(evaluate(expr.children[0],assignment),evaluate(expr.children[1],assignment))
    elif issubclass(type(expr), Abs):
        return abs(evaluate(expr.children[0],assignment))
    elif issubclass(type(expr), Neg):
        return -(evaluate(expr.children[0],assignment))
    elif issubclass(type(expr), Max):
        return max([evaluate(e,assignment) for e in expr.children])
    elif issubclass(type(expr), Min):
        return min([evaluate(e,assignment) for e in expr.children])
    elif issubclass(type(expr), Sum):
        return sum([evaluate(e,assignment) * expr.parameters[0][i] for (i,e) in enumerate(expr.children)])  + expr.parameters[1] # AS: these are additional constants that are moved from children to parameters in model.close()
    elif issubclass(type(expr), Function):
        return expr.parameters[0].get(tuple([evaluate(e,assignment) for e in expr.children]), expr.parameters[1])
    elif issubclass(type(expr), Disjunction):
        return (sum([evaluate(e,assignment) for e in expr.children])>0)
    elif issubclass(type(expr), Conjunction):
        return (sum([evaluate(e,assignment) for e in expr.children])==len(expr.children))
    elif issubclass(type(expr), Element):
        return evaluate(expr.children[evaluate(expr.children[-1],assignment)],assignment)
    else:
        raise ConstraintNotSupportedError(type(expr))

#SDG: automatic decompostion of any BinPredicate into a Table constraint with support tuples
def decompose_BinPredicate(self):
    arity = get_arity(self)
    scope = get_scope(self)
    res = Variable(self.get_lb(),self.get_ub())
    if arity is 0:
        val = evaluate(self, dict([]))
        res2 = Variable(val,val,str(val))
        return [res2]
    elif arity is 1:
        return [res, Table([res, scope[0]], [(evaluate(self, dict([(scope[0],val0)])), val0) for val0 in (scope[0].domain_ if scope[0].domain_ is not None else range(scope[0].get_lb(), scope[0].get_ub()+1))])]
    elif arity is 2:
        product = cartesian_product(self)  # check if decomposition is feasible in size
        return [res, Table([res, scope[0], scope[1]], [(evaluate(self, dict([(scope[0],val0), (scope[1],val1)])), val0, val1) for val0 in (scope[0].domain_ if scope[0].domain_ is not None else range(scope[0].get_lb(), scope[0].get_ub()+1)) for val1 in (scope[1].domain_ if scope[1].domain_ is not None else range(scope[1].get_lb(), scope[1].get_ub()+1))])]
    elif arity is 3:
        product = cartesian_product(self)  # check if decomposition is feasible in size
        return [res, Table([res, scope[0], scope[1], scope[2]], [(evaluate(self, dict([(scope[0],val0), (scope[1],val1), (scope[2],val2)])), val0, val1, val2) for val0 in (scope[0].domain_ if scope[0].domain_ is not None else range(scope[0].get_lb(), scope[0].get_ub()+1)) for val1 in (scope[1].domain_ if scope[1].domain_ is not None else range(scope[1].get_lb(), scope[1].get_ub()+1)) for val2 in (scope[2].domain_ if scope[2].domain_ is not None else range(scope[2].get_lb(), scope[2].get_ub()+1))])]
#    elif arity is 4:
#        product = cartesian_product(self)  # check if decomposition is feasible in size
#        return [res, Table([res, scope[0], scope[1], scope[2], scope[3]], [(evaluate(self, dict([(scope[0],val0), (scope[1],val1), (scope[2],val2), (scope[3],val3)])), val0, val1, val2, val3) for val0 in (scope[0].domain_ if scope[0].domain_ is not None else xrange(scope[0].get_lb(), scope[0].get_ub()+1)) for val1 in (scope[1].domain_ if scope[1].domain_ is not None else xrange(scope[1].get_lb(), scope[1].get_ub()+1)) for val2 in (scope[2].domain_ if scope[2].domain_ is not None else xrange(scope[2].get_lb(), scope[2].get_ub()+1)) for val3 in (scope[3].domain_ if scope[3].domain_ is not None else xrange(scope[3].get_lb(), scope[3].get_ub()+1))])]
    else:
        size = (self.get_ub(0)+1-self.get_lb(0)) * (self.get_ub(1)+1-self.get_lb(1))
        if (size > MODELMAXSIZE):
            raise ModelSizeError(size)
        tuples = [(self.eval(val0, val1), val0, val1) for val0 in (self.children[0].domain_ if issubclass(type(self.children[0]), Expression) and self.children[0].is_var() and self.children[0].domain_ is not None else range(self.get_lb(0), self.get_ub(0)+1)) for val1 in (self.children[1].domain_ if issubclass(type(self.children[1]), Expression) and self.children[1].is_var() and self.children[1].domain_ is not None else range(self.get_lb(1), self.get_ub(1)+1))]
        return [res, Table([res, self.children[0] if issubclass(type(self.children[0]), Expression) else Variable(self.children[0],self.children[0],str(self.children[0])), self.children[1] if issubclass(type(self.children[1]), Expression) else Variable(self.children[1],self.children[1],str(self.children[1]))], tuples)]

#SDG: decompose maximise into minimise
def decompose_Maximise(self):
    return [Minimise(Neg([self.children[0]]))]

#SDG: decompose minimise into cost functions
def decompose_Minimise(self):
    if type(self.children[0]) is int:
        return [PostNullary(self.children[0])]
    elif self.children[0].is_var():
        return [PostUnary(self.children[0], [c for c in range(self.children[0].get_lb(), self.children[0].get_ub()+1)])]
    elif issubclass(type(self.children[0]), Sum):
        res = []
        for i in range(len(self.children[0].children)):
            res.extend(decompose_Minimise(Minimise(Mul([self.children[0].children[i], self.children[0].parameters[0][i]]))))
        return res + [decompose_Minimise(Minimise(e)) for e in self.children[0].parameters[1:] if e is not 0]
    else:
        product = cartesian_product(self.children[0])  # check if decomposition is feasible in size
        arity = get_arity(self.children[0])
        scope = get_scope(self.children[0])
        #print "decompose ", self.children[0], " size: ", product, " arity: ", arity, " scope: ", scope
        if arity is 0:
            return [PostNullary(evaluate(self.children[0], dict([])))]
        elif arity is 1:
            costs = [evaluate(self.children[0], dict([(scope[0],val0)])) for val0 in (scope[0].domain_ if scope[0].domain_ is not None else range(scope[0].get_lb(), scope[0].get_ub()+1))]
            return [PostUnary(scope[0], costs)]
        elif arity is 2:
            costs = [evaluate(self.children[0], dict([(scope[0],val0), (scope[1],val1)])) for val0 in (scope[0].domain_ if scope[0].domain_ is not None else range(scope[0].get_lb(), scope[0].get_ub()+1)) for val1 in (scope[1].domain_ if scope[1].domain_ is not None else range(scope[1].get_lb(), scope[1].get_ub()+1))]
            return [PostBinary(scope[0], scope[1], costs)]
        else:
            obj = Variable(self.children[0].get_lb(), self.children[0].get_ub())
            return decompose_Minimise(Minimise(obj)) + [obj == self.children[0]]
