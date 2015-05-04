import Numberjack.solvers.Mistral as Mistral
import Numberjack.solvers.MiniSat as MiniSat

from .LinearTest import LinearTest
from .LogicalTest import LogicalTest
from .GlobalsATest import GlobalsATest
from .MiscTest import MiscTest
from .CoreTest import CoreTest
from .SATEncodingTest import SATEncodingTest

LinearTest.solver = Mistral.Solver
LogicalTest.solver = Mistral.Solver
MiscTest.solver = Mistral.Solver
GlobalsATest.solver = Mistral.Solver
CoreTest.solver = Mistral.Solver
SATEncodingTest.solver = MiniSat.Solver
