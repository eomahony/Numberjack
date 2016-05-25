'''

This module is responsible for testing the SAT encodings of Numberjack.

'''
from Numberjack import *
from copy import copy
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

    def testSmallDomainEq(self):
        v1 = Variable(2)
        m = Model(v1 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 1)

    def testSmallDomainLt(self):
        v1 = Variable(2)
        m = Model(v1 < 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), 1)

    def testSmallDomainLe(self):
        v1 = Variable(2)
        m = Model(v1 <= 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLessEqual(v1.get_value(), 1)

    def testSmallDomainGt(self):
        v1 = Variable(2)
        m = Model(v1 > 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 1)

    def testSmallDomainGe(self):
        v1 = Variable(2)
        m = Model(v1 >= 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 1)

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

    def testBooleanAddEqConstant(self):
        v1 = Variable()
        v2 = Variable()
        m = Model(v1 + v2 == 2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() + v2.get_value(), 2)

    def testBooleanAddEq(self):
        v1 = Variable()
        v2 = Variable()
        v3 = Variable(2)
        m = Model(v1 + v2 == v3, v3 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v3.get_value(), 1)
        self.assertEqual(v1.get_value() + v2.get_value(), v3.get_value())

    def testBooleanAddEq2(self):
        v1 = Variable()
        v2 = Variable()
        v3 = Variable(3)
        m = Model(v1 + v2 == v3, v3 == 2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() + v2.get_value(), v3.get_value())

    def testNegBooleanAddEq(self):
        v1 = Variable()
        v2 = Variable()
        m = Model(v1 + v2 >= 1, -1 * v1 + -1 * v2 >= -1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreaterEqual(v1.get_value() + v2.get_value(), 1)
        self.assertGreaterEqual(-1 * v1.get_value() + -1 * v2.get_value(), -1)

    def testAddGapsSat(self):
        v1 = Variable([0, 4, 8])
        v2 = Variable()
        v3 = Variable(10)
        m = Model(v1 + v2 == v3, v3 == 5)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() + v2.get_value(), v3.get_value())
        self.assertEqual(v3.get_value(), 5)

    def testAddGapsUnsat(self):
        v1 = Variable([0, 4, 8])
        v2 = Variable()
        v3 = Variable(10)
        m = Model(v1 + v2 == v3, v3 == 3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_unsat())

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

    def testSumSingleCoef(self):
        v1 = Variable(5)
        m = Model(Sum([v1], [-2]) == -4)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() * -2, -4)

    def testSumCoefEq(self):
        v1 = Variable(5)
        v2 = Variable(5)
        v3 = Variable(5)
        m = Model(Sum([v1, v2], [2, 2]) == v3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() * 2 + v2.get_value() * 2, v3.get_value())

    def testSumCoefNeg(self):
        v1 = Variable(5)
        v2 = Variable(5)
        m = Model(Sum([v1, v2], [-1, -1]) == -3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() * -1 + v2.get_value() * -1, -3)

    def testMulConstant(self):
        v1 = Variable(5)
        v2 = Variable(5)
        m = Model(v1 * 3 == v2)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() * 3, v2.get_value())

    def testMulVars(self):
        v1 = Variable(5)
        v2 = Variable(5)
        m = Model(v1 * v2 == 3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() * v2.get_value(), 3)

    def testMulPredicate(self):
        v1 = Variable(5)
        v2 = Variable(5)
        v3 = Variable(5)
        m = Model(v1 * v2 == v3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() * v2.get_value(), v3.get_value())

    def testMulGapsSat(self):
        v1 = Variable([0, 4, 8])
        v2 = Variable()
        v3 = Variable(10)
        m = Model(v1 * v2 == v3, v3 == 4)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value() * v2.get_value(), v3.get_value())
        self.assertEqual(v3.get_value(), 4)

    def testMulGapsUnsat(self):
        v1 = Variable([0, 4, 8])
        v2 = Variable()
        v3 = Variable(10)
        m = Model(v1 * v2 == v3, v3 == 3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_unsat())

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
        self.assertEqual(set(list(range(1, 6))), set(values))

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

    def testTableConflict(self):
        v1, v2 = VarArray(2, 1, 3)
        t = Table([v1, v2], [[1, 1], [2, 2], [3, 3]], type="conflict")
        m = Model(t)
        print(m)
        e = copy(SATEncodingTest.encoding)
        e.direct = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertNotEqual(v1.get_value(), v2.get_value())

    def testTableSupport(self):
        v1, v2 = VarArray(2, 1, 3)
        t = Table([v1, v2], [[1, 1], [2, 2], [3, 3]], type="support")
        m = Model(t)
        e = copy(SATEncodingTest.encoding)
        e.direct = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v2.get_value())

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
    #     v1 = Variable(1000)
    #     v2 = Variable(1000)
    #     m = Model(v1 == v2)
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.solve()
    #     self.assertTrue(s.is_sat())

    # def testAllDiffLargeDomain(self):
    #     vs = VarArray(100, 1, 1000)
    #     m = Model(AllDiff(vs))
    #     s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
    #     s.solve()
    #     self.assertTrue(s.is_sat())

    # ---------------- Absolute ----------------

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

    # ---------------- Reification of inequalities ----------------
    # Note that some of these tests require that the order encoding of the domain
    # is enabled. For now, this will just take a copy and set the bit to true,
    # but should find a better solution for parameterized test-cases in future.

    def testReifEqTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 == 3), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 1)

    def testReifEqFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 == 3), v2 == 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertNotEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 0)

    def testReifNeTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 != 3), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertNotEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 1)

    def testReifNeFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 != 3), v2 == 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 0)

    def testReifLeTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 <= 3), v2 == 1)
        e = copy(SATEncodingTest.encoding)
        e.order = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLessEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 1)

    def testReifLeFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 <= 3), v2 == 0)
        e = copy(SATEncodingTest.encoding)
        e.order = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreater(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 0)

    def testReifLtTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 < 3), v2 == 1)
        e = copy(SATEncodingTest.encoding)
        e.order = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 1)

    def testReifLtFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 < 3), v2 == 0)
        e = copy(SATEncodingTest.encoding)
        e.order = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreaterEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 0)

    def testReifGeTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 >= 3), v2 == 1)
        e = copy(SATEncodingTest.encoding)
        e.order = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreaterEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 1)

    def testReifGeFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 >= 3), v2 == 0)
        e = copy(SATEncodingTest.encoding)
        e.order = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 0)

    def testReifGtTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 > 3), v2 == 1)
        e = copy(SATEncodingTest.encoding)
        e.order = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreater(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 1)

    def testReifGtFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        m = Model(v2 == (v1 > 3), v2 == 0)
        e = copy(SATEncodingTest.encoding)
        e.order = True
        s = SATEncodingTest.solver(m, encoding=e)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLessEqual(v1.get_value(), 3)
        self.assertEqual(v2.get_value(), 0)

    def testReifAndTrue(self):
        v1 = Variable(10)
        v2 = Variable(10)
        m = Model(
            v1 == v2,
            ((v1 + 1) == v2) | ((v1 == 9) & (v2 == 9)),
        )
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v2.get_value())

    def testReifOrTrue(self):
        v1 = Variable(10)
        v2 = Variable(10)
        m = Model(
            v1 == v2,
            ((v1 + 1) == v2) | ((v1 == 9) | (v2 == 9)),
        )
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertTrue((v1.get_value() + 1 == v2.get_value()) or
                        (v1.get_value() == 9 or v2.get_value() == 9))

    def testOr(self):
        v1 = Variable(10)
        v2 = Variable(10)
        m = Model((v1 == 2) | (v2 == 2))
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertTrue(v1.get_value() == 2 or v2.get_value() == 2)

    def testOr2(self):
        v1 = Variable(10)
        v2 = Variable(10)
        m = Model((v1 == 2) | False)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertTrue(v1.get_value() == 2 or v2.get_value() == 2)

    # ---------------- Reification of inequalities between two expressions ----------------

    def testReif2ExprEqTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 == v3 + 2), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 1)

    def testReif2ExprEqFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 == v3 + 2), v2 == 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertNotEqual(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 0)

    def testReif2ExprNeTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 != v3 + 2), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertNotEqual(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 1)

    def testReif2ExprNeFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 != v3 + 2), v2 == 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 0)

    def testReif2ExprLeTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 <= v3 + 2), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLessEqual(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 1)

    def testReif2ExprLeFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 <= v3 + 2), v2 == 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreater(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 0)

    def testReif2ExprLtTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 < v3 + 2), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 1)

    def testReif2ExprLtFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 < v3 + 2), v2 == 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreaterEqual(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 0)

    def testReif2ExprGeTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 >= v3 + 2), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreaterEqual(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 1)

    def testReif2ExprGeFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 >= v3 + 2), v2 == 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLess(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 0)

    def testReif2ExprGtTrue(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 > v3 + 2), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertGreater(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 1)

    def testReif2ExprGtFalse(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 > v3 + 2), v2 == 0)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertLessEqual(v1.get_value(), v3.get_value() + 2)
        self.assertEqual(v2.get_value(), 0)

    def testReif2ExprMul(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(1, 5)
        v4 = Variable(1, 5)
        m = Model(v2 == (v1 == v3 * v4), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v3.get_value() * v4.get_value())
        self.assertEqual(v2.get_value(), 1)

    def testReif2ExprFactorDomain(self):
        v1 = Variable(5)
        v2 = Variable()
        v3 = Variable(5)
        m = Model(v2 == (v1 == v3 * 2), v2 == 1)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), v3.get_value() * 2)
        self.assertEqual(v2.get_value(), 1)

    # ---------------- Modulus ----------------

    def testModConstant(self):
        v1 = Variable(10)
        v2 = Variable(5)
        m = Model(v2 == (v1 % 3), v1 == 4)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), 4)
        self.assertEqual(v2.get_value(), 1)

    def testModNegConstant(self):
        v1 = Variable(-5, 5)
        v2 = Variable(-5, 5)
        m = Model(v2 == (v1 % 3), v1 == -5)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), -5)
        self.assertEqual(v2.get_value(), -2)

    def testModTwoVars(self):
        v1 = Variable(-10, 10)
        v2 = Variable(-5, 5)
        v3 = Variable(-5, 5)
        m = Model(v3 == (v1 % v2), v1 == -5, v2 == 3)
        s = SATEncodingTest.solver(m, encoding=SATEncodingTest.encoding)
        s.solve()
        self.assertTrue(s.is_sat())
        self.assertEqual(v1.get_value(), -5)
        self.assertEqual(v2.get_value(), 3)
        self.assertEqual(v3.get_value(), -2)
