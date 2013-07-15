#import Numberjack

MODELMAXSIZE = 1000000

from Numberjack import *

'''
Custom constraints section
'''

class MyConstraint(Expression):
    
    def __init__(self, vars):
        Expression.__init__(self, "MyConstraint")
        
        self.set_children(vars)
    
        print "This is the MyConstraint method"
    
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
        lb0 = self.children[0] if type(self.children[0]) is int else self.children[0].lb
        ub0 = self.children[0] if type(self.children[0]) is int else self.children[0].ub
        lb1 = self.children[1] if type(self.children[1]) is int else self.children[1].lb
        ub1 = self.children[1] if type(self.children[1]) is int else self.children[1].ub
        self.lb = lb0 + lb1
        self.ub = ub0 + ub1

    def get_symbol(self):
        return '+'


# SDG: new Predicate for expressing cost functions
## Base class of cost functions
#
class CostFunction(Predicate):
    def __init__(self, vars, operator):
        Predicate.__init__(self, vars, operator)
        self.lb = None
        self.ub = None


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
# \warning Can only be used with Toulbar2 solver
#
#

class PostBinary(CostFunction):
    
    ## PostBinary constraint constructor
    # @param var1 first variable
    # @param var2 second variable 
    # @param costs table of costs
    #
	# \note
	#   - Top-level: PostTernary Constraint
	#   - Nested: Cannot be used as a nested predicate
	#
    # \code
    #var = Variable(1,3)
	#var1 = Variable(1,3)
	#post = PostBinary(var,var1,[5,3,2,2,2,2,2,1,2])
	#model = Model(post)
    #print model
    #>>> assign:
    #>>>   x0 in {1..3}
    #>>>   x1 in {1..3}
    #>>>   
    #>>> subject to:
    #>>>   PostBinary(x0, x1)
    # \endcode
    #
    
    def __init__(self, var1, var2, costs):
        CostFunction.__init__(self, [var1, var2], "PostBinary")
        self.parameters = [costs]

    def decompose(self):
        obj = Variable(0,max(self.parameters[0]))
        var2size = self.get_children()[1].get_size()
        return [Table([obj] + self.get_children(), [(w, i / var2size , i % var2size) for i,w in enumerate(self.parameters[0])], 'support'), Minimize(obj)]

## PostTernary Constraint
# @param var first variable
# @param var1 second variable 
# @param var2 third variable 
# @param costs table of costs
#
# \note
#   - Top-level: PostTernary Constraint
#   - Nested: Cannot be used as a nested predicate
#
# \code
# var = Variable(2)
# var1 = Variable(2)
# var2 = Variable(2)
# post = PostTernary(var,var1,var2,[5,3,2,2,2,2,2,1])
# model = Model(post)
# \endcode
      
class PostTernary(CostFunction):
    def __init__(self, var1, var2, var3, costs):
        CostFunction.__init__(self, [var1, var2, var3], "PostTernary")
        self.parameters = [costs]

    def decompose(self):
        obj = Variable(0,max(self.parameters[0]))
        var2size = self.get_children()[1].get_size()
        var3size = self.get_children()[2].get_size()
        return [Table([obj] + self.get_children(), [(w, i / var2size / var3size , i / var3size, i % var3size) for i,w in enumerate(self.parameters[0])], 'support'), Minimize(obj)]

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
            print "You're using WVarAmong constraint, bounds are set by the Variable passed in the last parameter of the constraint constructor!"

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
        
## PostNary Constraint
#
# \note
#   - Top-level: PostNary Constraint
#   - Nested: Cannot be used as a nested predicate
#      
class PostNary(CostFunction):
    def __init__(self, vars, arity, default_cost):
        CostFunction.__init__(self, vars, "PostNary")
        self.parameters = [arity, default_cost,[],[]]
        
    def add(self, tupleIndex, cost):
            self.parameters[2].append(tupleIndex)
            self.parameters[3].append(cost)

## PostWSameGcc Constraint
#
# \note
#   - Top-level: PostWSameGcc Constraint
#   - Nested: Cannot be used as a nested predicate
#
class PostWSameGcc(CostFunction):
    def __init__(self, varsL, varsR, cards, type, semantics, baseCost):
        CostFunction.__init__(self, [varsL, varsR], "PostWSameGcc")
        values = cards.keys()
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
        u = u | set([e] if type(e) is int else range(e.lb, e.ub + 1))
    res = Variable(list(u))
    return ([res, (self.children[-1] >= 0), (self.children[-1] < len(self.children)-1)] + [((res == (Variable(e,e,str(e)) if type(e) is int else e)) | (self.children[-1] != i)) for i, e in enumerate(self.children[:-1])])

#SDG: decompose LeqLex and LessLex
def decompose_LessLex(self):
    length = len(self.children) / 2
    def lexico(vars1,vars2,i):
        if (i == len(vars1)-1):
            return (vars1[i] < vars2[i])
        else:
            return (vars1[i] < vars2[i] | And([(vars1[i] == vars2[i]), lexico(vars1,vars2,i+1)]))
    return [lexico(self.children[:length], self.children[length:], 0)]

def decompose_LeqLex(self):
    length = len(self.children) / 2
    def lexico(vars1,vars2,i):
        if (i == len(vars1)-1):
            return (vars1[i] <= vars2[i])
        else:
            return (vars1[i] <= vars2[i] | And([(vars1[i] == vars2[i]), lexico(vars1,vars2,i+1)]))
    return [lexico(self.children[:length], self.children[length:], 0)]

#SDG: automatic decompostion of any BinPredicate into a Table constraint with support tuples
def decompose_BinPredicate(self):
    tuples = [(self.eval(val1, val2), val1, val2) for val1 in range(self.get_lb(0), self.get_ub(0)+1) for val2 in range(self.get_lb(1), self.get_ub(1)+1)]
#    res = Variable(min(t[0] for t in tuples), max(t[0] for t in tuples))
    res = Variable(self.get_lb(),self.get_ub())
    return ([res, Table([res, self.children[0] if issubclass(type(self.children[0]), Expression) else Variable(self.children[0],self.children[0],str(self.children[0])), self.children[1] if issubclass(type(self.children[1]), Expression) else Variable(self.children[1],self.children[1],str(self.children[1]))], tuples)])

#SDG: decompose maximise into minimise
def decompose_Maximise(self):
    return [Minimise(Neg([self.children[0]]))]

#SDG: decompose minimise into cost functions
def decompose_Minimise(self):
    def cartesian_product(expr):
        def prod( iterable ):
            p= 1
            for n in iterable:
                p *= n
            return p

        if type(expr) in [int, long, float, bool, str]:
            return 1
        elif expr.is_var():
            return expr.get_size()
        else:
            res = prod(cartesian_product(e) for e in expr.children)
            if (res > MODELMAXSIZE):
                raise ModelSizeError("Memory limit exceeded: cannot decompose this expression!")
            return res
        
    def get_arity(expr):
        if type(expr) in [int, long, float, bool, str]:
            return 0
        elif expr.is_var():
            if expr.get_lb() == expr.get_ub():
                return 0
            else:
                return 1
        else:
            return sum(get_arity(e) for e in expr.children)
        
    def get_scope(expr):
        if type(expr) in [int, long, float, bool, str]:
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
            return res
        
    def evaluate(expr, assignment):
        if type(expr) in [int, long, float, bool, str]:
            return expr
        elif expr.is_var():
            return assignment[expr]
        elif issubclass(type(expr), BinPredicate):
            return expr.eval(evaluate(expr.children[0],assignment),evaluate(expr.children[1],assignment))
        else:
            raise ModelSizeError("Unknown predicate: cannot decompose this expression!")
        
    if type(self.children[0]) is int:
        return [PostNullary(self)]
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
#        print "decompose ", self.children[0], " size: ", product, " arity: ", arity, " scope: ", scope
        if arity is 0:
            return [PostNullary(self.children[0])]
        elif arity is 1:
            costs = [evaluate(self.children[0], dict([(scope[0],val1)])) for val1 in range(scope[0].get_lb(), scope[0].get_ub()+1)]
            return [PostUnary(scope[0], costs)]
        else:
            raise ModelSizeError("Arity too large: cannot decompose this expression!")
