'''
Created on Feb 11, 2016

@author: Alexander Schiendorfer, ISSE Augsburg
'''
import unittest

from Numberjack import *;
from Numberjack.Decomp import *;


class TestDecompEvaluate(unittest.TestCase):
    

    def setUp(self):
        self.x,self.y,self.z = VarArray(3,0,2)
        x, y, z = self.x, self.y, self.z

        self.model = Model()
        self.assignment = dict([(x,0), (y,1), (z,2)])


    def tearDown(self):
        pass

    def testLiteral(self):
        self.model.close()
        res = evaluate(5, self.assignment)
        self.assertEqual(res, 5, "Should be literal 5")
        res = evaluate(5.16, self.assignment)
        self.assertEqual(res, 5.16, "Should be literal 5.16")
        res = evaluate(True, self.assignment)
        self.assertTrue(res, "Should be literal True")
        res = evaluate(False, self.assignment)
        self.assertFalse(res, "Should be literal False")
        res = evaluate("giuseppe", self.assignment)
        self.assertEqual(res, "giuseppe", "Should be literal giuseppe")
        
        
    def testSum(self):
        self.expr2 = self.x + 1
        self.expr3 = self.x + self.y + self.z - 1
        self.model.add(self.expr2)
        self.model.add(self.expr3)
        self.model.close()
        
        res = evaluate(self.expr2, self.assignment)
        # x is assigned to 0, so (x+1) should be 1
        self.assertEqual(res, 1, "Should be 1")
        
        # x -> 0, y -> 1, z -> 2, thus x + y + z - 1 should be 2
        res = evaluate(self.expr3, self.assignment)
        self.assertEqual(res, 2, "Should be 1")
     
    def testMaxMin(self):
        self.expr2 = Min([self.x, 1]) # should be 0
        self.expr3 = Max([self.x, self.y, self.z + 1]) # should be 3
        self.model.add(self.expr2)
        self.model.add(self.expr3)
        self.model.close()
        
        res = evaluate(self.expr2, self.assignment)
        self.assertEqual(res, 0, "Should be 0")
        
        res = evaluate(self.expr3, self.assignment)
        self.assertEqual(res, 3, "Should be 1")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSum']
    unittest.main()