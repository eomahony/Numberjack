import Mistral
from LinearTest import *
from LogicalTest import *
from GlobalsATest import *
from MiscTest import *
from CoreTest import *

LinearTest.solver = Mistral.Solver
LogicalTest.solver = Mistral.Solver
MiscTest.solver = Mistral.Solver
GlobalsATest.solver = Mistral.Solver
CoreTest.solver = Mistral.Solver
