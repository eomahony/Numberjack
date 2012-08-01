'''

This module will contain test classes for Numberjack's linear programming
components. Namely the <, <=, >= == and +,-,*,/ operations

'''
from Numberjack import *
import unittest

class LinearBinTest:
    def __init__(self):
        self.functs = {"__eq__":[0],
                       "__lt__":[-1],
                       "__le__":[0, -1],
                       "__gt__":[1],
                       "__ge__":[0, 1]}
    
    def set_up(self):
        self.model = Model()
        #print self.model
        self.variables = [Variable(l, u) for l, u in self.domains]
        self.model.add(getattr(getattr(self.variables[0], self.art_op)(self.variables[1]),
                               self.comp_op)(self.rhs))
        
    def test_assign(self):
        #print self.model
        return getattr(getattr(self.variables[0].get_value(),
                               self.art_op)(self.variables[1].get_value()),
                       "__cmp__")(self.rhs) in self.functs[self.comp_op]

class Test(LinearBinTest):
    def __init__(self, name, comp, art, rhs = 20, domains = [(0, 10), (0, 10)] ):
        LinearBinTest.__init__(self)
        self.test_name = name
        self.comp_op = comp
        self.art_op = art
        self.domains = domains
        self.rhs = rhs
        self.set_up()

# Test for success
TestSuccess = [
    ["Test Eq 1", "__eq__", "__add__"],
    ["Test Eq 2", "__eq__", "__sub__", 5],
    ["Test Le 1", "__lt__", "__add__"],
    ["Test Le 2", "__lt__", "__sub__"],
    ["Test Leq 1", "__le__", "__add__"],
    ["Test Leq 2", "__le__", "__sub__"],
    ["Test Ge 1", "__gt__", "__add__", 19],
    ["Test Ge 2", "__gt__", "__sub__", 4],
    ["Test Geq 1", "__ge__", "__add__"],
    ["Test Geq 2", "__ge__", "__sub__", 5],
]

# Test for failure
TestFailures = [
    ["Test Eq 1", "__eq__", "__add__", 21],
    ["Test Eq 2", "__eq__", "__sub__", 11],
    ["Test Le 1", "__lt__", "__add__", -1],
    ["Test Le 2", "__lt__", "__sub__", -11],
    ["Test Leq 1", "__le__", "__add__", -1],
    ["Test Leq 2", "__le__", "__sub__", -11],
    ["Test Ge 1", "__gt__", "__add__", 20],
    ["Test Ge 2", "__gt__", "__sub__", 11],
    ["Test Geq 1", "__ge__", "__add__", 21],
    ["Test Geq 2", "__ge__", "__sub__", 11],
]

# Test Boundary
TestBoundary = [
    ["Test Eq 1", "__eq__", "__add__", -20, [(-10, 0),(-10, 0)]],
    ["Test Eq 2", "__eq__", "__sub__", 0, [(-10, 0),(-10, 0)]],
    ["Test Le 1", "__lt__", "__add__", 0, [(-10, 0),(-10, 0)]],
    ["Test Le 2", "__lt__", "__sub__", 0, [(-10, 0),(-10, 0)]],
    ["Test Leq 1", "__le__", "__add__", 0, [(-10, 0),(-10, 0)]],
    ["Test Leq 2", "__le__", "__sub__", 0, [(-10, 0),(-10, 0)]],
    ["Test Ge 1", "__gt__", "__add__", -1, [(-10, 0),(-10, 0)]],
    ["Test Ge 2", "__gt__", "__sub__", -1, [(-10, 0),(-10, 0)]],
    ["Test Geq 1", "__ge__", "__add__", -1, [(-10, 0),(-10, 0)]],
    ["Test Geq 2", "__ge__", "__sub__", -1, [(-10, 0),(-10, 0)]],
]

class LinearTest(unittest.TestCase):
    
    solver = None
            
    def testSuccess1(self):
        te = Test(*TestSuccess[0])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testSuccess2(self):
        te = Test(*TestSuccess[1])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testSuccess3(self):
        te = Test(*TestSuccess[2])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testSuccess4(self):
        te = Test(*TestSuccess[3])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testSuccess5(self):
        te = Test(*TestSuccess[4])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testSuccess6(self):
        te = Test(*TestSuccess[5])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testSuccess7(self):
        te = Test(*TestSuccess[6])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testSuccess8(self):
        te = Test(*TestSuccess[7])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())

    def testSuccess9(self):
        te = Test(*TestSuccess[8])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testSuccess10(self):
        te = Test(*TestSuccess[9])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
            
    def testFailure1(self):
        te = Test(*TestFailures[0])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
        
    def testFailure2(self):
        te = Test(*TestFailures[1])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
        
    def testFailure3(self):
        te = Test(*TestFailures[2])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
        
    def testFailure4(self):
        te = Test(*TestFailures[3])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
        
    def testFailure5(self):
        te = Test(*TestFailures[4])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
        
    def testFailure6(self):
        te = Test(*TestFailures[5])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
        
    def testFailure7(self):
        te = Test(*TestFailures[6])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
    
    def testFailure8(self):
        te = Test(*TestFailures[7])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
        
    def testFailure9(self):
        te = Test(*TestFailures[8])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
        
    def testFailure10(self):
        te = Test(*TestFailures[9])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertFalse(solver.solve())
            
    def testBoundary1(self):
        te = Test(*TestBoundary[0])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testBoundary2(self):
        te = Test(*TestBoundary[1])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testBoundary3(self):
        te = Test(*TestBoundary[2])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testBoundary4(self):
        te = Test(*TestBoundary[3])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testBoundary5(self):
        te = Test(*TestBoundary[4])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testBoundary6(self):
        te = Test(*TestBoundary[5])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testBoundary7(self):
        te = Test(*TestBoundary[6])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
    
    def testBoundary8(self):
        te = Test(*TestBoundary[7])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testBoundary9(self):
        te = Test(*TestBoundary[8])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
    def testBoundary10(self):
        te = Test(*TestBoundary[9])
        solver = LinearTest.solver(te.model)
        solver.setVerbosity(0)
        self.assertTrue(solver.solve())
        self.assertTrue(te.test_assign())
        
def run():
    unittest.main()
            
if __name__ == "__main__":
    #for solver in ['SCIP', 'Mistral', 'OsiClp', 'OsiGlpk', 'OsiDylp', 'OsiVol', 'OsiSym']: 
    for solver in ['Mistral']:
        print solver
        solver_module = __import__(solver)
        LinearTest.solver = solver_module.Solver
        unittest.main(exit=False)
