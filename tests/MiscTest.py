'''

This module is responsible for Testing miscellaneous models in Numberjack

'''
from Numberjack import *
import unittest


class MiscTest(unittest.TestCase):

    solver = None

    def testQueens(self):
        queens = VarArray(6, 6)
        model = Model(
            AllDiff(queens),
            AllDiff([queens[i] + i for i in range(6)]),
            AllDiff([queens[i] - i for i in range(6)])
        )
        solver = MiscTest.solver(model)
        self.assertTrue(solver.solve())

    def testMagicSquare(self):
        n = 5

        square = Matrix(n, n, 1, n * n)

        sum_val = n * (n * n + 1) // 2  # This is what all the columns, rows and diagonals must add up tp

        model = Model(
            AllDiff(square.flat),

            [Sum(row) == sum_val for row in square.row],
            [Sum(col) == sum_val for col in square.col],

            Sum([square[a][a] for a in range(n)]) == sum_val,
            Sum([square[a][n - a - 1] for a in range(n)]) == sum_val
        )

        solver = MiscTest.solver(model)
        self.assertTrue(solver.solve())
        for row in square:
            self.assertTrue(sum([v.get_value() for v in row]) == sum_val)
        for col in square.col:
            self.assertTrue(sum([v.get_value() for v in col]) == sum_val)
        self.assertTrue(sum([v.get_value() for v in [square[a][a] for a in range(n)]]) == sum_val)
        self.assertTrue(sum([v.get_value() for v in [square[a][n - a - 1] for a in range(n)]]) == sum_val)

    def testSendMoreMoney(self):
        # Create the model
        model = Model()

        #Variable definitions, we need a variable for every letter
        s, m = (Variable(1, 9) for val in range(2))  # These can't be zero as they are the start of a word
        e, n, d, o, r, y = (Variable(0, 9) for val in range(6))  # These can

        model.add(      s*1000 + e*100 + n*10 + d +
                        m*1000 + o*100 + r*10 + e ==
              m*10000 + o*1000 + n*100 + e*10 + y)

        model.add(AllDiff((s, e, n, d, m, o, r, y)))  # Post the all different constraint on all the variables

        solver = MiscTest.solver(model)
        self.assertTrue(solver.solve())

        self.assertTrue(sum([v.get_value() * x for v, x in [(s, 1000), (e, 100), (n, 10), (d, 1)]]) +
                        sum([v.get_value() * x for v, x in [(m, 1000), (o, 100), (r, 10), (e, 1)]]) ==
                        sum([v.get_value() * x for v, x in [(m, 10000), (o, 1000), (n, 100), (e, 10), (y, 1)]]))

    def testBIBD(self):
        v, b, r, k, l = [7, 7, 3, 3, 1]
        # Create the variable matrix and the model
        matrix = Matrix(b, v)
        model = Model(
            [Sum(row) == k for row in matrix.row],  # every row adds up to k
            [Sum(col) == r for col in matrix.col],  # every column adds up to r

            # the scalar product of every pair of columns adds up to l
            [Sum([(row[col_i] & row[col_j]) for row in matrix.row]) == l
             for col_i in range(v) for col_j in range(col_i)]
            )
        solver = MiscTest.solver(model)
        self.assertTrue(solver.solve())

    def testCostas(self):
        N = 6
        sequence = VarArray(N, 1, N)
        model = Model(
            AllDiff(sequence),
            [AllDiff([sequence[j] - sequence[j + i + 1] for j in range(N - i - 1)]) for i in range(N - 2)]
        )

        solver = MiscTest.solver(model)
        self.assertTrue(solver.solve())

    def TestLangford(self):
        N, M = [10, 3]
        X = [Variable(1, N * M - (i + 2) * (M - 1)) for i in range(N)]

        model = Model(
            AllDiff([X[i] + ((i + 2) * j) for j in range(M) for i in range(N)]),
            X[0] > X[1]  # Break symmetry
        )

        solver = MiscTest.solver(model)
        self.assertTrue(solver.solve())

    def testGolomb(self):
        nbMarks = 4

        rulerSize = 2 ** (nbMarks - 1)
        marks = VarArray(nbMarks, rulerSize)
        model = Model(
            Minimise(marks[nbMarks - 1]),  # objective function

            [marks[i - 1] < marks[i] for i in range(1, nbMarks)],
            AllDiff([marks[i] - marks[j] for i in range(1, nbMarks) for j in range(i)]),
            marks[0] == 0  # symmetry breaking
        )

        solver = MiscTest.solver(model)
        self.assertTrue(solver.solve())

    #def testOpticalModelling(self):
    #
    #    alarm_matrix = [
    #        [1, 2, 3,                   10],
    #        [                  7          ],
    #        [               6, 7,         ],
    #        [            5, 6, 7,         ],
    #        [   2, 3, 4,          8,    10],
    #        [      3, 4,          8, 9, 10],
    #    ]
    #
    #    Monitors, Nodes, alarm_matrix = [10, 6, alarm_matrix]
    #
    #    class VectorPairAllDiff(Expression):
    #        def __init__(self, row1, row2):
    #            Expression.__init__(self, "VectorAllDiff")
    #            self.children = row1
    #            self.rows = [row1, row2]
    #
    #        def decompose(self):
    #            return [ Sum([var1 != var2 for var1, var2 in zip(self.rows[0], self.rows[1])]) > 0 ]
    #
    #    monitors_on = VarArray(Monitors)
    #    being_monitored = Matrix(Nodes, Monitors)
    #
    #    model = Model()
    #
    #    model.add(Minimise(Sum(monitors_on)))
    #    # Link the monotirs on to being alarmed
    #    model.add( [ monitor == ( Sum(col) >= 1 ) for col, monitor in zip(being_monitored.col, monitors_on) ])
    #    # Link alarm_matrix to variable matrix
    #    for monitored_row, possible_monitor_row in zip(being_monitored, alarm_matrix):
    #        model.add([monitored_row[idx - 1] == 0 for idx in
    #                   [x for x in range(Monitors) if x not in possible_monitor_row]])
    #    # Monotor everything
    #    model.add( [ Sum(row) > 0 for row in being_monitored] )
    #    # Make sure to have different alarms
    #    model.add([VectorPairAllDiff(x1, x2) for x1, x2 in pair_of(being_monitored)])
    #
    #    solver = MiscTest.solver(model)
    #    self.assertTrue(solver.solve())

    def testSensorNetwork(self):

        field_cover = [
            [1, 2, 3         ],
            [1, 2, 3, 4      ],
            [1,    3, 4, 5, 6],
            [1,    3,       6],
            [   2, 3, 4, 5, 6]
        ]
        fields, covers, no_sensors, field_cover = [5, 3, 6, field_cover]

        sensors = Matrix(covers, no_sensors)
        field_covered = Matrix(covers, fields)

        model = Model()

        # Use each sensor only once
        model.add([Sum(col) == 1 for col in sensors.col])

        # Make sure all the fields are covered
        for sens, field in zip(sensors, field_covered):  # For every cover
            model.add([Sum([sens[sens_idx - 1] for sens_idx in field_cover[field_idx]]) >= 1 for field_idx in range(fields)])

        solver = MiscTest.solver(model)
        self.assertTrue(solver.solve())
