'''

This module is responsible for testing core/internal elements of Numberjack.

'''
from Numberjack import *
import unittest


class CoreTest(unittest.TestCase):
    solver = None

    # Tests checking that the overloaded Variable constructor conforms
    # to its definition in the documentation.

    def testVariableBinary(self):
        # Variable() :- Binary variable
        v = Variable()
        self.assertEqual(v.lb, 0)
        self.assertEqual(v.ub, 1)

    def testVariableRange(self):
        # Variable(N) :- Variable in the domain of {0, N-1}
        v = Variable(10)
        self.assertEqual(v.lb, 0)
        self.assertEqual(v.ub, 9)

    def testVariableFloatRange(self):
        v = Variable(10.0)
        self.assertEqual(type(v.lb), float)
        self.assertEqual(type(v.ub), float)
        self.assertAlmostEqual(v.lb, 0.0)
        self.assertAlmostEqual(v.ub, 9.0)

    def testVariableName(self):
        # Variable('x') :- Binary variable called 'x'
        v = Variable('y')
        self.assertEqual(v.name(), 'y')

    def testVariableUbName(self):
        # Variable(N, 'x') :- Variable in the domain of {0, N-1} called 'x'
        v = Variable(10, 'y')
        self.assertEqual(v.name(), 'y')
        self.assertEqual(v.lb, 0)
        self.assertEqual(v.ub, 9)

    def testVariableLowerUpper(self):
        # Variable(l,u) :- Variable in the domain of {l, u}
        v = Variable(0, 1)
        self.assertEqual(v.lb, 0)
        self.assertEqual(v.ub, 1)

        v = Variable(10, 25)
        self.assertEqual(v.lb, 10)
        self.assertEqual(v.ub, 25)

        v = Variable(15, 15)
        self.assertEqual(v.lb, 15)
        self.assertEqual(v.ub, 15)

    def testVariableLowerUpperFloat(self):
        v = Variable(0.0, 2.0)
        self.assertAlmostEqual(v.lb, 0.0)
        self.assertAlmostEqual(v.ub, 2.0)

    def testVariableRangeName(self):
        # Variable(l,u, 'x') :- Variable in the domain of {l, u} called 'x'
        v = Variable(15, 20, 'y')
        self.assertEqual(v.lb, 15)
        self.assertEqual(v.ub, 20)
        self.assertEqual(v.name(), 'y')

    def testVariableList(self):
        # Variable(list) :- Variable with domain specified as a list
        list_domain = [2, 8, 6, 4]
        v = Variable(list_domain)
        self.assertEqual(v.lb, 2)
        self.assertEqual(v.ub, 8)
        self.assertEqual(set(list_domain), set(v.get_domain()))

    def testVariableListName(self):
        # Variable(list, 'x') :- Variable with domain specified as a list called 'x'
        list_domain = [2, 8, 6, 4]
        v = Variable(list_domain, 'y')
        self.assertEqual(v.name(), 'y')
        self.assertEqual(v.lb, 2)
        self.assertEqual(v.ub, 8)
        self.assertEqual(set(list_domain), set(v.get_domain()))

    def testDomainListIterator(self):
        "Test the iterator code on Domain returns the correct values when given the domain as a list."
        list_domain = list(range(10))
        x = Variable(list_domain)
        list_from_iterator = [v for v in x]
        self.assertEqual(set(list_domain), set(list_from_iterator))

    def testDomainBoundedIterator(self):
        "Test the iterator code on Domain returns the correct values when given the domain as bounds."
        x = Variable(5, 10)
        list_from_iterator = [v for v in x]
        self.assertEqual(set(list(range(5, 11))), set(list_from_iterator))

    def testBuiltExpressionDomainRange(self):
        v = Variable(10)
        model = Model(v)
        solver = CoreTest.solver(model)
        self.assertEqual(set(list(range(10))), set(v.get_domain()))

    def testBuiltExpressionDomainWithGaps(self):
        list_domain = list(range(0, 10, 2))
        v = Variable(list_domain)
        model = Model(v)
        solver = CoreTest.solver(model)
        self.assertEqual(set(list_domain), set(v.get_domain()))

    def testExpressionRangeSize(self):
        v = Variable(10)
        self.assertEqual(v.get_size(), 10)

    def testBuiltExpressionRangeSize(self):
        v = Variable(10)
        model = Model(v)
        solver = CoreTest.solver(model)
        self.assertEqual(v.get_size(), 10)

    def testModelLoad(self):
        "Tests that we can load a solver from a model by name. We assume that Mistral will be available at a minimum."
        import Numberjack.solvers.Mistral as Mistral
        m = Model()
        solver = m.load('Mistral')
        self.assertIsInstance(solver, Mistral.Solver)

    def testModelLoadNonExistantSolver(self):
        m = Model()
        self.assertRaises(ImportError, m.load, 'solverdoesnotexist')

    # Tests for Matrix constructor

    def testMatrixList(self):
        # M = Matrix(l) creates a Matrix from a list l
        l = [1, 2, 3]
        m = Matrix(l)
        self.assertEqual(len(m), len(l))
        for i in range(len(l)):
            self.assertEqual(len(m[i]), l[i])
            for j in range(l[i]):
                v = m[i][j]
                self.assertIsInstance(v, Variable)
                self.assertEqual(v.lb, 0)
                self.assertEqual(v.ub, 1)
                self.assertEqual(v.name(), "x%d.%d" % (i, j))

    def testMatrixBoolean(self):
        # - M = Matrix(n, m) creates a n x m Matrix of Boolean variables
        n, m = 3, 3
        self.verifyMatrixOfVariables(Matrix(n, m), n, m, 0, 1)

    def testMatrixBooleanName(self):
        # - M = Matrix(n, m, 'x') creates a n x m Matrix of Boolean variables with names 'x0.0..xn-1.m-1'
        n, m, name = 3, 3, 'y'
        self.verifyMatrixOfVariables(Matrix(n, m, name), n, m, 0, 1, name_prefix=name)

    def testMatrixUb(self):
        # - M = Matrix(n, m, u) creates a n x m Matrix of variables with domains [0..u-1]
        n, m, ub = 3, 3, 10
        self.verifyMatrixOfVariables(Matrix(n, m, ub), n, m, 0, ub - 1)

    def testMatrixUbName(self):
        # - M = Matrix(n, m, u, 'x') creates a n x m Matrix of variables with domains [0..u-1] and names 'x0.0..xn-1.m-1'
        n, m, ub, name = 3, 3, 10, 'y'
        self.verifyMatrixOfVariables(Matrix(n, m, ub, name), n, m, 0, ub - 1, name_prefix=name)

    def testMatrixLbUb(self):
        # - M = Matrix(n, m, l, u) creates a n x m Matrix of variables with domains [l..u]
        n, m, lb, ub = 3, 3, 5, 10
        self.verifyMatrixOfVariables(Matrix(n, m, lb, ub), n, m, lb, ub)

    def testMatrixLbUbName(self):
        # - M = Matrix(n, m, l, u, 'x') creates a n x m Matrix of variables with domains [l..u] and names 'x0.0..xn-1.m-1'
        n, m, lb, ub, name = 3, 3, 5, 10, 'y'
        self.verifyMatrixOfVariables(Matrix(n, m, lb, ub, name), n, m, lb, ub, name_prefix=name)

    def verifyMatrixOfVariables(self, matrix, n, m, lb, ub, name_prefix='x'):
        self.assertEqual(len(matrix), n)
        for i in range(n):
            self.assertEqual(len(matrix[i]), m)
            for j in range(m):
                v = matrix[i][j]
                self.assertIsInstance(v, Variable)
                self.assertEqual(v.lb, lb)
                self.assertEqual(v.ub, ub)
                self.assertEqual(v.name(), "%s%d.%d" % (name_prefix, i, j))

    def testMatrixSliceRows(self):
        n, m = 5, 5
        matrix = Matrix(n, m)
        m_slice = matrix[1:3]
        print("slice:", repr(m_slice))
        self.assertEqual(len(m_slice), 2)
        for row in m_slice:
            self.assertEqual(len(row), m)

    def testMatrixSliceRowsCols(self):
        n, m = 5, 5
        matrix = Matrix(n, m)
        m_slice = matrix[1:3, 1:4]
        self.assertEqual(len(m_slice), 2)
        for row in m_slice:
            self.assertEqual(len(row), 3)

    def testMatrixRowColIndex(self):
        n, m = 5, 5
        matrix = Matrix(n, m)
        v = matrix[2, 2]
        self.assertEqual(v.name(), "x2.2")
