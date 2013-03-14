'''

This module is responsible for testing the SAT encodings of Numberjack.

'''
from Numberjack import *
import unittest


class SATEncodingTest(unittest.TestCase):
    solver = None
    # encoding = EncodingConfiguration(direct=True, order=False, conflict=True, support=True, amo_encoding=AMOEncoding.Pairwise)   # Direct & Support Encoding
    # encoding = EncodingConfiguration(direct=True, order=False, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise)  # Direct Encoding
    # encoding = EncodingConfiguration(direct=True, order=False, conflict=False, support=True, amo_encoding=AMOEncoding.Pairwise)  # Support Encoding
    # encoding = EncodingConfiguration(direct=False, order=True, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise)  # Order Encoding
    encoding = EncodingConfiguration(direct=True, order=True, conflict=True, support=False, amo_encoding=AMOEncoding.Pairwise)  # Direct & Order Encoding

    # Tests checking that the overloaded Variable constructor conforms
    # to its definition in the documentation.

    def testEq1(self):
        v1 = Variable(1, 5)
        m = Model(v1 == 3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 3)

    def testEq2(self):
        v1 = Variable(1, 5)
        v2 = Variable(3, 8)
        m = Model(v1 == v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v2.get_value())

    def testEq3(self):
        v1 = Variable(1, 4)
        v2 = Variable(4, 6)
        m = Model(v1 == v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v2.get_value())

    def testEqGap(self):
        v1 = Variable([1, 3, 5, 7])
        v2 = Variable(4, 8)
        m = Model(v1 == v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v2.get_value())

    # def testConstantDomain(self):
    #     v1 = Variable(1, 4)
    #     v2 = Variable([2])
    #     m = Model(v1 == v2)
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.solve()
    #     self.assertTrue(s.is_sat())
    #     self.assertEqual(v1.get_value(), v2.get_value())

    def testOffsetDomain(self):
        v1 = Variable(1, 4)
        v2 = Variable(1, 5)
        m = Model(v1 + 3 == v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() + 3, v2.get_value())

    def testNe1(self):
        v1 = Variable(1, 2)
        m = Model(v1 != 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertNotEqual(v1.get_value(), 1)

    def testNe2(self):
        v1 = Variable(1, 2)
        v2 = Variable(1, 2)
        m = Model(v1 != v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertNotEqual(v1.get_value(), v2.get_value())

    def testGt(self):
        v1 = Variable(1, 3)
        v2 = Variable(1, 3)
        m = Model(v1 > v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreater(v1.get_value(), v2.get_value())

    def testGe(self):
        v1 = Variable(1, 3)
        v2 = Variable(1, 3)
        m = Model(v1 >= v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreaterEqual(v1.get_value(), v2.get_value())

    def testLt(self):
        v1 = Variable(1, 3)
        v2 = Variable(1, 3)
        m = Model(v1 < v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), v2.get_value())

    def testLe(self):
        v1 = Variable(1, 3)
        v2 = Variable(1, 3)
        m = Model(v1 <= v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLessEqual(v1.get_value(), v2.get_value())

    def testLtGap(self):
        v1 = Variable([1, 3, 5])
        v2 = Variable([3])
        m = Model(v1 < v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), v2.get_value())

    def testAddEqConstant(self):
        v1 = Variable(5)
        v2 = Variable(5)
        m = Model(v1 + v2 == 6)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() + v2.get_value(), 6)

    def testAddEq(self):
        v1 = Variable(5)
        v2 = Variable(5)
        v3 = Variable(5)
        m = Model(v1 + v2 == v3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() + v2.get_value(), v3.get_value())

    def testSubEq(self):
        v1 = Variable(5)
        v2 = Variable(5)
        v3 = Variable(5)
        m = Model(v1 - v2 == v3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() - v2.get_value(), v3.get_value())

    def testAddMultiple(self):
        v1 = Variable(5)
        v2 = Variable(5)
        v3 = Variable(5)
        v4 = Variable(5)
        v5 = Variable(5)
        m = Model(v1 + v2 + v3 + v4 == v5)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() + v2.get_value() + v3.get_value() + v4.get_value(), v5.get_value())

    def testSumEq(self):
        v1 = Variable(5)
        v2 = Variable(5)
        v3 = Variable(5)
        m = Model(Sum([v1, v2]) == v3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() + v2.get_value(), v3.get_value())

    def testSumArray(self):
        vs = VarArray(4, 1, 5)
        v_sum = Variable(5)
        m = Model(Sum(vs) == v_sum)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(sum(v.get_value() for v in vs), v_sum.get_value())

    def testMulConstant(self):
        v1 = Variable(5)
        v2 = Variable(5)
        m = Model(v1 * 3 == v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() * 3, v2.get_value())

    # FIXME Encoding multiplication between expressions is not supported yet.
    # def testMulVars(self):
    #     v1 = Variable(5)
    #     v2 = Variable(5)
    #     m = Model(v1 * v2 == 3)
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.solve()
    #     self.assertTrue(s.is_sat())
    #     self.assertEqual(v1.get_value() * v2.get_value(), 3)

    # def testMulVars(self):
    #     v1 = Variable(5)
    #     v2 = Variable(5)
    #     v3 = Variable(5)
    #     m = Model(v1 * v2 == v3)
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.solve()
    #     self.assertTrue(s.is_sat())
    #     self.assertEqual(v1.get_value() * v2.get_value(), v3.get_value())

    # def testDivConstant(self):
    #     v1 = Variable(5)
    #     v2 = Variable(5)
    #     m = Model(v1 / 3 == v2)
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.solve()
    #     self.assertTrue(s.is_sat())
    #     self.assertEqual(v1.get_value() / 3, v2.get_value())

    def testAllDiff(self):
        vs = VarArray(5, 1, 5)
        m = Model(AllDiff(vs))
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        values = [v.get_value() for v in vs]
        self.assertItemsEqual(range(1, 6), values)

    def testMaximise(self):
        v1 = Variable(5)
        v2 = Variable(5)
        m = Model(v1 < v2, Maximise(v1))
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 4)

    def testMinimise(self):
        v1 = Variable(5)
        v2 = Variable(5)
        m = Model(v1 < v2, Minimise(v2))
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 0)
        self.assertEqual(v2.get_value(), 1)

    def testFormula(self):
        a, b, c = VarArray(3)
        m = Model(((a == True) | (b == False)),
                  ((b == True) | (a == False)),
                  (c == True))
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertTrue(c.get_value())
        self.assertEqual(a.get_value(), b.get_value())

    def testFormula2(self):
        a, b, c = VarArray(3)
        m = Model(((a == True) | (b == False)) &
                  ((b == True) | (a == False)) &
                  (c == True))
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertTrue(c.get_value())
        self.assertEqual(a.get_value(), b.get_value())

    # FIXME Encoding of Table constraint is not supported yet.
    # def testTable(self):
    #     v1, v2 = VarArray(2, 1, 3)
    #     t = Table([v1, v2], [[1, 1], [2, 2], [3, 3]])
    #     m = Model(t)
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.solve()
    #     self.assertTrue(s.is_sat())
    #     self.assertNotEqual(v1.get_value(), v2.get_value())

    def testMax(self):
        v1 = Variable(5)
        v2 = Variable(5)
        m = Model(Max([v1, v2]) < 2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), 2)
        self.assertLess(v2.get_value(), 2)

    def testMaxNegVars(self):
        v1 = Variable(-5, -1)
        v2 = Variable(-5, -1)
        m = Model(Max([v1, v2]) < -2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), -2)
        self.assertLess(v2.get_value(), -2)

    # def testLargeDomain(self):
    #     v1 = Variable(10000)
    #     v2 = Variable(10000)
    #     m = Model(v1 == v2)
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.output_cnf("large_domain.cnf")
    #     s.solve()
    #     self.assertTrue(s.is_sat())
    #     import sys
    #     print >> sys.stderr, v1.get_value(), v2.get_value()

    # def testAllDiffLargeDomain(self):
    #     vs = VarArray(100, 1, 10000)
    #     m = Model(AllDiff(vs))
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.output_cnf("large_alldiff.cnf")
    #     s.solve()
    #     self.assertTrue(s.is_sat())

    def testAbsEq(self):
        v1 = Variable(-4, -1)
        m = Model(Abs(v1) == 2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), -2)

    def testAbsGt(self):
        v1 = Variable(-4, -1)
        m = Model(Abs(v1) > 2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), -2)

    def testAbsLt(self):
        v1 = Variable(-4, -1)
        m = Model(Abs(v1) < 2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreater(v1.get_value(), -2)
