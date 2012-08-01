'''

This module will contain test classes for Numberjack's logical
components. Namely the & and | operations

'''
from Numberjack import *
import unittest

class LogicalTest(unittest.TestCase):
    solver = None
    
    def testAndSucc(self):
        m = Model()
        var1, var2 = [Variable(), Variable()]
        m.add(var1 & var2)
        s = LogicalTest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var1.get_value() == 1 and var2.get_value() == 1)
    
    def testAndFail(self):
        m = Model()
        var1, var2 = [Variable(), Variable(0, 0)]
        m.add(var1 & var2)
        s = LogicalTest.solver(m)
        self.assertFalse(s.solve())
    
    def testOrSucc(self):
        m = Model()
        var1, var2 = [Variable(), Variable()]
        m.add(var1 | var2)
        s = LogicalTest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var1.get_value() == 1 or var2.get_value() == 1)
    
    def testOrFail(self):
        m = Model()
        var1, var2 = [Variable(0, 0), Variable(0, 0)]
        m.add(var1 | var2)
        s = LogicalTest.solver(m)
        self.assertFalse(s.solve())
    
    def testAndReif1(self):
        m = Model()
        var1, var2 = [Variable(), Variable()]
        m.add(var1 & var2 == 0)
        s = LogicalTest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var1.get_value() == 0 or var2.get_value() == 0)
        
    def testAndReif2(self):
        m = Model()
        var1, var2 = [Variable(), Variable()]
        m.add(var1 & var2 == 1)
        s = LogicalTest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var1.get_value() == 1 and var2.get_value() == 1)
    
    def testOrReif1(self):
        m = Model()
        var1, var2 = [Variable(), Variable()]
        m.add(var1 | var2 == 0)
        s = LogicalTest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var1.get_value() == 0 and var2.get_value() == 0)
        
    def testOrReif2(self):
        m = Model()
        var1, var2 = [Variable(), Variable()]
        m.add(var1 | var2 == 1)
        s = LogicalTest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var1.get_value() == 1 or var2.get_value() == 1)

if __name__ == "__main__":
    import Mistral
    LogicalTest.solver = Mistral.Solver
    unittest.main()