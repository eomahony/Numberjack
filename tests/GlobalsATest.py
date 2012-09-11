'''

This module should handle testing of global constraints that recieve a large
amound of use. Specifically the AllDiff

'''
from Numberjack import *
import unittest
import sys
import random

class GlobalsATest(unittest.TestCase):
    solver = None
    
    '''
    AllDiff Tests
    '''
    def testAllDiffBinarySucc(self):
        m = Model()
        var1, var2 = [Variable([0]), Variable([1])]
        m.add(AllDiff([var1,var2]))
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue([var1.get_value(),var2.get_value()].count(1) == 1)
    
    def testAllDiffBinaryFail(self):
        m = Model()
        var1, var2 = [Variable([0]), Variable([0])]
        m.add(AllDiff([var1,var2]))
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
    
    def testAllDiffBooleanSucc(self):
        m = Model()
        var1, var2 = [Variable(), Variable()]
        m.add(AllDiff([var1,var2]))
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue([var1.get_value(),var2.get_value()].count(1) == 1)

    #def testAllDiffSingleSucc(self):
    #    m = Model()
    #    var = Variable(0,10)
    #    m.add(AllDiff([var]))
    #    s = GlobalsATest.solver(m)
    #    self.assertTrue(s.solve())
    #    self.assertTrue(var.get_value() == 0)
    #
    #def testAllDiffSingleBooleanSucc(self):
    #    m = Model()
    #    var = Variable()
    #    m.add(AllDiff([var]))
    #    s = GlobalsATest.solver(m)
    #    self.assertTrue(s.solve())
    #    self.assertTrue(var.get_value() == 0)

    def testAllDiffMulipleBooleanFail(self):
        m = Model()
        varray = [Variable() for i in range(3)]
        m.add(AllDiff(varray))
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
    
    def testAllDiffLargeVarraySucc(self):
        m = Model()
        varray = [Variable(0,301) for i in range(300)]
        m.add(AllDiff(varray))
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(set([i.get_value() for i in varray]) == set(range(0,300)))
    
    def testAllDiffLargeVarraySucc(self):
        m = Model()
        varray = [Variable(0,301) for i in range(400)]
        m.add(AllDiff(varray))
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
    
    '''
    Sum Tests
    '''
    def testSumSingleBooleanSucc(self):
        m = Model()
        var = Variable()
        m.add(Sum([var])==1)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var.get_value() == 1)
        
    def testSumSingleBooleanFail(self):
        m = Model()
        var = Variable()
        m.add(Sum([var])==2)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
        
    def testSumMultiBooleanSucc(self):
        m = Model()
        varray = [Variable() for i in range(5)]
        m.add(Sum(varray)==5)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(sum([v.get_value() for v in varray]) == 5)
    
    def testSumMultiBooleanFail(self):
        m = Model()
        varray = [Variable() for i in range(5)]
        m.add(Sum(varray)==30)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
    
    def testSumSingleVarSucc(self):
        m = Model()
        var = Variable(0,5)
        m.add(Sum([var])==1)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var.get_value() == 1)
        
    def testSumSingleVarFail(self):
        m = Model()
        var = Variable(0,5)
        m.add(Sum([var])==20)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
    
    def testSumMultiVarSucc(self):
        m = Model()
        varray = [Variable(0,5) for i in range(10)]
        m.add(Sum(varray)==30)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue( [v.get_value() for v in varray].count(5)==6 )
        
    def testSumMultiVarFail(self):
        m = Model()
        varray = [Variable(0,5) for i in range(10)]
        m.add(Sum(varray)==100)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
        
    def testSumLargeBooleanArraySucc(self):
        m = Model()
        varray = [Variable() for i in range(5000)]
        m.add(Sum(varray) == 4500)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(sum([v.get_value() for v in varray]) == 4500)
        
    def testSumLargeBolleanArrayFail(self):
        m = Model()
        varray = [Variable() for i in range(5000)]
        m.add(Sum(varray) == 5500)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
        
    def testSumLargeVarArraySucc(self):
        m = Model()
        varray = [Variable(0,20) for i in range(5000)]
        m.add(Sum(varray) == 32000)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(sum([v.get_value() for v in varray]) == 32000)
        
    def testSumLargeVarArrayFail(self):
        m = Model()
        varray = [Variable(0,20) for i in range(5000)]
        m.add(Sum(varray) == 320000)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
        
    def testSumHugeDomainArraySucc(self):
        m = Model()
        varray = [Variable(0,1000000) for i in range(300)]
        m.add(Sum(varray) == 123456789)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(sum([v.get_value() for v in varray]) == 123456789)
        
    def testWeightedSumSingleSucc(self):
        m = Model()
        var = Variable(0,10)
        m.add(Sum([var],[7])==21)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var.get_value() == 3)

    def testWeightedSumSingleFail(self):
        m = Model()
        var = Variable(0,10)
        m.add(Sum([var],[7])==72)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
    
    def testWeightedSumSingleBooleanSucc(self):
        m = Model()
        var = Variable()
        m.add(Sum([var],[4]) == 4)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        self.assertTrue(var.get_value() == 1)

    def testWeightedSumSingleBooleanFail(self):
        m = Model()
        var = Variable()
        m.add(Sum([var],[4]) == 3)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
        
    def testWeightedSumMultiVarSucc(self):
        m = Model()
        varray = [ Variable(0,20) for i in range(5) ]
        weights = [2,6,2,7,8]
        m.add(Sum(varray,weights) == 30)
        s = GlobalsATest.solver(m)
        self.assertTrue(s.solve())
        count = 0
        for i in range(len(varray)):
            count += varray[i].get_value()*weights[i]
        self.assertTrue(count == 30)

    def testWeightedSumMultiVarFail(self):
        m = Model()
        varray = [ Variable(0,20) for i in range(30) ]
        weights = [2,6,2,7,5,6,2,6,7,9] + [1 for i in range(20)]
        m.add(Sum(varray,weights) == 16000)
        s = GlobalsATest.solver(m)
        self.assertFalse(s.solve())
        
    #def testElementSingleBooleanSucc(self):
    #def testElementMultiBooleanSucc(self):
    #def testElementLargeBooleanSucc(self):
    #def testElementSingleVarSucc(self):
    #def testElementMultiVarSucc(self):
    #def testElementLargeVarSucc(self):
    
    #def testGccSingleBooleanSucc(self):
    #def testGccMultiBooleanFail(self):
    #def testGccLargeBooleanFail(self):
    #def testGccSingleVarSucc(self):
    #def testGccMultiVarSucc(self):
    #def testGccLargeVarSucc(self):
    
if __name__ == "__main__":
    import Mistral
    GlobalsATest.solver = Mistral.Solver
    unittest.main()
