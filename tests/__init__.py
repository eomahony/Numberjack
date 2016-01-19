import Numberjack.solvers.Mistral as TestSolver
# import Numberjack.solvers.CPLEX as TestSolver
import Numberjack.solvers.MiniSat as MiniSat

from .LinearTest import LinearTest
from .LogicalTest import LogicalTest
from .GlobalsATest import GlobalsATest
from .MiscTest import MiscTest
from .CoreTest import CoreTest
from .SATEncodingTest import SATEncodingTest

LinearTest.solver = TestSolver.Solver
LogicalTest.solver = TestSolver.Solver
MiscTest.solver = TestSolver.Solver
GlobalsATest.solver = TestSolver.Solver
CoreTest.solver = TestSolver.Solver
SATEncodingTest.solver = MiniSat.Solver
