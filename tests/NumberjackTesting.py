
import unittest
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

if __name__ == "__main__":
    
    # First test Mistral
    
    suite1 = unittest.TestLoader().loadTestsFromTestCase(LinearTest)
    suite1.addTests(unittest.TestLoader().loadTestsFromTestCase(LogicalTest))
    suite1.addTests(unittest.TestLoader().loadTestsFromTestCase(MiscTest))
    suite1.addTests(unittest.TestLoader().loadTestsFromTestCase(GlobalsATest))
    suite1.addTests(unittest.TestLoader().loadTestsFromTestCase(CoreTest))
    
    unittest.TextTestRunner(verbosity=2).run(suite1)
    
    # Then test SCIP
    #import SCIP
    #LinearTest.solver = SCIP.Solver
    #LogicalTest.solver = SCIP.Solver
    #MiscTest.solver = SCIP.Solver
    #GlobalsATest.solver = SCIP.Solver    

    #suite1 = unittest.TestLoader().loadTestsFromTestCase(LinearTest)
    #suite1.addTests(unittest.TestLoader().loadTestsFromTestCase(LogicalTest))
    #suite1.addTests(unittest.TestLoader().loadTestsFromTestCase(MiscTest))
    #suite1.addTests(unittest.TestLoader().loadTestsFromTestCase(GlobalsATest))
    #
    #unittest.TextTestRunner(verbosity=2).run(suite1)
    
