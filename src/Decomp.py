#import Numberjack

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


'''
Custom decomposition section
'''
def decompose_AllDiff(self):
    return [var1 != var2 for var1, var2 in pair_of(self.children)]
    
    
    
