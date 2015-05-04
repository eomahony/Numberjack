'''
  Numberjack is a constraint satisfaction and optimisation library
  Copyright (C) 2009 Cork Constraint Computation Center, UCC
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.
  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The authors can be contacted electronically at 
  numberjack.support@gmail.com
'''

import unittest
import sys

from Numberjack import *
from Mistral import Solver

class Test(unittest.TestCase):
    
    def testVariable_backtrack(self):
        '''
        This is a quick test of the ability of the variable to backtrack properly
        '''
        # Removal test
        var = Variable(list(range(0,10)))
        assert(var.can_be_instantiated_to(5))
        var.backtrack_stamp()
        var.remove_value(5)
        assert(not var.can_be_instantiated_to(5))
        var.backtrack()
        assert(var.can_be_instantiated_to(5))
        
        # Assignment test
        var = Variable(list(range(0,10)))
        assert(not var.get_is_instantiated())
        var.backtrack_stamp()
        var.set_value(5)
        assert(not var.can_be_instantiated_to(4))
        assert(var.get_is_instantiated())
        assert(var.get_value() == 5)
        var.backtrack()
        assert(not var.get_is_instantiated())
        assert(var.can_be_instantiated_to(4))
        
    def testEqual(self):
        var1 = Variable([0])
        var2 = Variable(list(range(0,1)))
        
        model = NativeModel()
        
        model.add_variable((var1, var2))
        
        solver = Solver(model)
        assert(solver.solve())


    def testStupid_not_eq(self):
        var1 = Variable(list(range(0,3)))
        var2 = Variable(list(range(0,3)))
        
        model = NativeModel()
        model.add_variable((var1, var2))
        model.add_constraint(NotEqual((var1, var2)))
            
        solver = Solver(model)
        
        assert(solver.solve())    
        assert(var1.get_value() != var2.get_value())
        
    def testPropogationHandler(self):
        var1 = Variable(list(range(0,3)))
        var2 = Variable(list(range(0,3)))
        
        neq = NotEqual((var1, var2)).repr_con
        
        neq.set_up_for_solver()
        
        ph = Solving.PropogationHandler()
        
        ph.add_to_stack(neq)
        
        var1.set_propogation_handler(ph)
        var2.set_propogation_handler(ph)
        
        assert(len(ph._PropogationHandler__prop_stack) == 1)
        
        var1.remove_value(1)
        var2.remove_value(0)
        
        assert(len(ph._PropogationHandler__prop_stack) == 1)
        
                
    def testNot_eq_test2(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        
        model = NativeModel()
        model.add_variable((var1, var2, var3))
                           
        model.add_constraint(NotEqual((var1, var2)))
        model.add_constraint(NotEqual((var1, var3)))
        model.add_constraint(NotEqual((var2, var3)))   
        
        bc = BasicConstraint([var1])
        
        def bc_prop(obj):
            vars = obj.get_variables()
            if vars[0].get_is_instantiated() and vars[0].get_value() == 1:
                obj.fail()
            
        bc.rt_propogate = bc_prop
        
        model.add_constraint(bc)
        
        solver = Solver(model)
        
        assert(solver.solve()) 
        
    def testGEQ(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        model = NativeModel()
        model.add_variable((var1, var2, var3))
        
        model.add_constraint(Geq((var1, var2)))
        model.add_constraint(Geq((var2, var3)))
        model.add_constraint(Geq((var1, var3)))
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d %d " % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testLEQ(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        model = NativeModel()
        model.add_variable((var1, var2, var3))
        
        model.add_constraint(Leq((var1, var2)))
        model.add_constraint(Leq((var2, var3)))
        model.add_constraint(Leq((var1, var3)))
        
        model.add_constraint(Equal((var1, var2)))
        model.add_constraint(NotEqual((var2, var3)))
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d %d " % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testLT(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        model = NativeModel()
        model.add_variable((var1, var2, var3))
        
        model.add_constraint(Lt((var1, var2)))
        model.add_constraint(Lt((var2, var3)))
        model.add_constraint(Lt((var1, var3)))
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d %d " % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testGT(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        model = NativeModel()
        model.add_variable((var1, var2, var3))
        
        model.add_constraint(Gt((var1, var2)))
        model.add_constraint(Gt((var2, var3)))
        model.add_constraint(Gt((var1, var3)))
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d %d " % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testPlus(self):
        var1, var2 = (Variable(list(range(0,4))), Variable(list(range(0,4))))
      
        model = NativeModel()
        model.add_variable((var1, var2))
        
        model.add_constraint(Equal((var1, Plus(var2, 3))))
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d " % (var1.get_value(), var2.get_value())
        
    def testMinus(self):
        var1, var2 = (Variable(list(range(0,4))), Variable(list(range(0,4))))
        
        model = NativeModel()
        model.add_variable((var1, var2))
        
        model.add_constraint(Equal((var1, Minus(var2, 3))))
        
        solver = Solver(model)
        assert(solver.solve())
        
    def testTimes(self):
        var1, var2 = (Variable(list(range(1,4))), Variable(list(range(0,4))))
        
        model = NativeModel()
        
        model.add_constraint(var1 == var2 * 2)
        
        solver = Solver(model)
        assert(solver.solve())
        
    def testVariableBounds(self):
        var = Variable(list(range(0,10)))
        
        assert(var.get_lower() == 0)
        assert(var.get_upper() == 9)
        
        var.backtrack_stamp()
        
        var.set_lower(4)
        assert(var.get_lower() == 4)
        assert(var.get_upper() == 9)
        
        var.set_upper(6)
        assert(var.get_lower() == 4)
        assert(var.get_upper() == 6)
        
        var.backtrack()
        
        assert(var.get_lower() == 0)
        assert(var.get_upper() == 9)
        
    def testSum(self):
        var1 = Variable(list(range(0,3)))
        var2 = Variable(list(range(0,3)))
        var3 = Variable(list(range(3,5)))
        
        model = NativeModel()
        
        
        svar = Sum((var1, var2))
        
        model.add_variable((var1, var2, var3, svar))
        
        model.add_constraint(Equal((svar, var3)))
        
        solver = Solver(model)
        
        assert(solver.solve())
        
        #svar.print_domain()
        
        #var1.print_domain()
        #var2.print_domain()
        #var3.print_domain()
        
        #print "%d %d %d" % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testWeightedSum(self):
        var1 = Variable(list(range(2,3)))
        var2 = Variable(list(range(2,3)))
        var3 = Variable(list(range(3,10)))
        
        model = NativeModel()
        
        
        svar = Sum(((var1, var2), (2,1)))
        
        model.add_variable((var1, var2, var3, svar))
        
        model.add_constraint(Equal((svar, var3)))
        
        solver = Solver(model)
        
        #svar.print_domain()
        
        assert( solver.solve() )
        
        #svar.print_domain()
        
        #var1.print_domain()
        #var2.print_domain()
        
        #var3.print_domain()
        
        #print "%d %d %d" % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testSumAgain(self):
        var1 = Variable(list(range(1,5)))
        var2 = Variable(list(range(1,5)))
        var3 = Variable(list(range(1,10)))
        var4 = Variable(list(range(1,10)))
        
        sum1 = Sum(((var1, var2),(2,2)))
        sum2 = Sum((var3, var4))
        
        model = NativeModel()
        
        model.add_variable((var1, var2, var3, var4, sum1, sum2))
        
        model.add_constraint(Equal((sum1, sum2)))
        
        solver = Solver(model)
        
        assert( solver.solve() )
        
    def testSumAgainLeq(self):
        var1 = Variable(list(range(1,5)))
        var2 = Variable(list(range(1,5)))
        var3 = Variable(list(range(1,10)))
        var4 = Variable(list(range(1,10)))
        
        sum1 = Sum(((var1, var2),(2,2)))
        sum2 = Sum((var3, var4))
        
        model = NativeModel()
        
        model.add_variable((var1, var2, var3, var4, sum1, sum2))
        
        model.add_constraint(Gt((sum1, sum2)))
        
        solver = Solver(model)
        
        assert( solver.solve() )
        
    def testSumConstant(self):
         obj1 = Variable(list(range(0,6)))
         obj2 = Variable(list(range(0,8)))
         obj3 = Variable(list(range(0,11)))
    
         capacity = Variable(list(range(34,35)))
         volumes = [7, 5, 3]
    
         model = NativeModel()
    
         vol_sum = Sum(((obj1, obj2, obj3), volumes))
    
         model.add_variable((obj1, obj2, obj3, vol_sum, capacity))
    
         model.add_constraint(Leq((vol_sum, capacity)))
         
         solver = Solver(model)
        
         '''
         obj1.print_domain()
         obj2.print_domain()
         obj3.print_domain()
         vol_sum.print_domain()
         capacity.print_domain()
         '''
         
         assert( solver.solve() )
         
         '''
         obj1.print_domain()
         obj2.print_domain()
         obj3.print_domain()
         vol_sum.print_domain()
         capacity.print_domain()
         '''
        
    #def testBinSearch(self):
      #  var1 = Variable(range(0,3))
        #assert( var1.find([1,3,4], 2)[0] == 0)
        
    def testSumProp(self):
        var1, var2, var3, var4 = (Variable(list(range(1,5))) for x in range(0,4))
        
        var3 = Variable([1,10])
        var4 = Variable([2,10])
        
        sum1 = Sum((var1, var2))
        sum2 = Sum(((var3, var4), (1,2)))
        
        model = NativeModel()
        
        model.add_variable((var1, var2, var3, var4, sum1, sum2))
        model.add_constraint(Equal((sum1, sum2)))
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print [x.get_value() for x in (var1, var2, var3, var4)]
        
        assert(var1.get_value() + var2.get_value() == var3.get_value() + var4.get_value()*2)
        
    def testVarExtractionNeq(self):
        var1, var2 = (Variable(list(range(0,3))) for x in range(0,2))
        model = NativeModel()
        model.add_constraint(NotEqual((var1, var2)))
        solver = Solver(model)
        assert(solver.solve())
        #print [v.get_value() for v in (var1, var2)]
        
    def testVarExtractionWSum(self):
        var1 = Variable(list(range(2,3)))
        var2 = Variable(list(range(2,3)))
        var3 = Variable(list(range(3,10)))
        model = NativeModel()  
        
        model.add_constraint(Equal((Sum(((var1, var2), (2,1))), var3)))
        
        solver = Solver(model)
        
        #svar.print_domain()
        
        assert( solver.solve() )
        
    '''
    Operator overloading tests
    '''
        
    def testOperatorOVerloading(self):
        var1 = Variable(list(range(0,4)))
        var2 = Variable(list(range(0,4)))
        
        model = NativeModel()
        model.add_constraint(var1 != var2)
                             
        solver = Solver(model)
        assert(solver.solve())
        
    def testEqualOver(self):
        var1 = Variable([0])
        var2 = Variable(list(range(0,1)))
        
        model = NativeModel()
        
        model.add_constraint(var1 == var2)
        
        solver = Solver(model)
        assert(solver.solve())


    def testStupid_not_eqOver(self):
        var1 = Variable(list(range(0,3)))
        var2 = Variable(list(range(0,3)))
        
        model = NativeModel()
        model.add_constraint(var1 != var2)
            
        solver = Solver(model)
        
        assert(solver.solve())    
        assert(var1.get_value() != var2.get_value())
        
        
    def testGEQOver(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        model = NativeModel()
        
        model.add_constraint(var1 >= var2)
        model.add_constraint(var2 >= var3)
        model.add_constraint(var1 >= var3)
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d %d " % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testLEQOver(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        model = NativeModel()
        
        model.add_constraint(var1 <= var2)
        model.add_constraint(var2 <= var3)
        model.add_constraint(var1 <= var3)
        
        model.add_constraint(var1 == var2)
        model.add_constraint(var2 != var3)
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d %d " % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testLTOver(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        model = NativeModel()
        
        model.add_constraint(var1 < var2)
        model.add_constraint(var2 < var3)
        model.add_constraint(var1 < var3)
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d %d " % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testGTOver(self):
        var1, var2, var3 = (Variable(list(range(0,3))) for x in range(0,3))
        model = NativeModel()
        
        model.add_constraint(var1 < var2)
        model.add_constraint(var2 < var3)
        model.add_constraint(var1 < var3)
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d %d " % (var1.get_value(), var2.get_value(), var3.get_value())
        
    def testPlusOver(self):
        var1, var2 = (Variable(list(range(0,4))), Variable(list(range(0,4))))
      
        model = NativeModel()
        
        model.add_constraint(var1 == var2 + 3)
        
        solver = Solver(model)
        assert(solver.solve())
        
        #print "%d %d " % (var1.get_value(), var2.get_value())
        
    def testMinusOver(self):
        var1, var2 = (Variable(list(range(0,4))), Variable(list(range(0,4))))
        
        model = NativeModel()
        
        model.add_constraint(var1 == var2 - 3)
        
        solver = Solver(model)
        assert(solver.solve())
        
    def testTimesOver(self):
        var1, var2 = (Variable(list(range(1,4))), Variable(list(range(0,4))))
        
        model = NativeModel()
        
        model.add_constraint(var1 == var2 * 2)
        
        solver = Solver(model)
        assert(solver.solve())
        
    def testMinusOverNativeModelOp(self):
        var1, var2 = (Variable(list(range(0,4))), Variable(list(range(0,4))))
        
        model = NativeModel()
        
        model << (var1 == var2 - 3)
        
        solver = Solver(model)
        assert(solver.solve())
        
    def testTimesOverNativeModelOp(self):
        var1, var2 = (Variable(list(range(1,4))), Variable(list(range(0,4))))
        
        model = NativeModel()
        
        model << (var1 == var2 * 2)
        
        solver = Solver(model)
        assert(solver.solve())
        
    def testMaximization(self):
        model = NativeModel()
        
        var1, var2 = (Variable(0,4) for x in range(0,2))
        var3 = Variable(0,10)
        
        model << (var1 + var2 == var3)
        
        solver = Solver(model)
        
        assert( solver.maximise(var3) )
        
        #print [v.get_value() for v in (var1, var2, var3)]
        
    def testNegDomains(self):
        var1 = Variable(-3, -1)
        var2 = Variable(-2, 0)
        
        model = NativeModel()

        model << (var1*2 <= var2)
        
        solver = Solver(model)
        
        assert(solver.solve())
        
    '''
    Tests of stuff implemented from minizinc stuff
    '''
    
    def testAtLeastSAT(self):
        model = NativeModel()
        var1, var2 = (Variable(0,5), Variable(0,5))
        
        model << (AtLeast(([var1, var2], 2, 2)))
        
        solver = Solver(model)
        
        assert(solver.solve())
        assert(var1.get_value() is 2)
        assert(var2.get_value() is 2)
        
    def testAtLeastUNSAT_1(self):
        model = NativeModel()
        var1, var2 = (Variable(0,5), Variable(0,5))
        
        model << (AtLeast(([var1, var2], 6, 2)))
        
        solver = Solver(model)
        
        assert(not solver.solve())
        
    def testAtLeastUNSAT_2(self):
        model = NativeModel()
        var1, var2 = (Variable(0,5), Variable(0,5))
        
        model << (AtLeast(([var1, var2], 0, 3)))
        
        solver = Solver(model)
        
        assert(not solver.solve())
        
    def testAtMost(self):
        model = NativeModel()
        var1, var2 = (Variable(2,2), Variable(0,5))
        
        model << (AtMost(([var1, var2], 2, 1)))
        
        solver = Solver(model)
        
        assert(solver.solve())
        
    def testAtMostUNSAT_1(self):
        model = NativeModel()
        var1, var2 = (Variable(2,2), Variable(0,5))
        
        model << (AtMost(([var1, var2], 2, 0)))
        
        solver = Solver(model)
        
        assert(not solver.solve())
        
    def testAtMostUNSAT_2(self):
        model = NativeModel()
        var1, var2 = (Variable(2,2), Variable(0,1))
        
        model << (AtMost(([var1, var2], 2, 2)))
        
        solver = Solver(model)
        
        assert(solver.solve())
        
    def testExactlySAT(self):
        model = NativeModel()
        var1, var2 = (Variable(0,5), Variable(0,5))
        
        model << (Exactly(((var1, var2), 2, 2)))
        
        solver = Solver(model)
        
        assert(solver.solve())
        assert(var1.get_value() is 2)
        assert(var2.get_value() is 2)
        
    def testExactlySAT_2(self):
        model = NativeModel()
        var1, var2 = (Variable(0,5), Variable(0,5))
        
        model << (Exactly(((var1, var2), 2, 1)))
        
        solver = Solver(model)
        
        assert(solver.solve())
        assert((var1.get_value() is 2 or var2.get_value() is 2))
        
    def testExactlyUNSAT_1(self):
        model = NativeModel()
        var1, var2 = (Variable(0,5), Variable(0,1))
        
        model << (Exactly(((var1, var2), 2, 2)))
        
        solver = Solver(model)
        
        assert(not solver.solve())
        
    def testExactlyUNSAT_2(self):
        model = NativeModel()
        var1, var2 = (Variable(0,5), Variable(0,5))
        
        model << (Exactly(((var1, var2), 2, 3)))
        
        solver = Solver(model)
        
        assert(not solver.solve())
        
    def taestElementSAT_1(self):
        model = NativeModel()
        
        i = Variable(0,5)
        k = Variable(1,5)
        coef = [0,1,2,3,4,5]
        # Make sure that coef[k] = i
        
        model << (Element((i, k, coef)))
        
        solver = Solver(model)
        
        assert(solver.solve())
        assert(i.get_value() is coef[k.get_value()])
        
    def testMinusVars(self):
        var1 = Variable(0, 3)
        var2 = Variable(0, 1)
        var3 = Variable(0, 0)
        
        model = NativeModel()
        
        minVar = var1 - var2

        model << (minVar == var3)
        
        solver = Solver(model)
        
        
        assert( solver.solve() )
        
    def testNAnd_1(self):
        var1 = Variable(0, 1)
        var2 = Variable(1, 2)
        var3 = Variable(0, 1)
        
        model = NativeModel()
        
        model << (Nand((var1 != var2, var2 != var3)))
        
        solver = Solver(model)
        
        assert ( solver.solve() )
        
        #print [var.get_value() for var in (var1, var2, var3)]
        
    def testNand_2(self):
        
        var1 = Variable(0, 0)
        var2 = Variable(0, 0)
        var3 = Variable(0, 0)
        
        model = NativeModel()
        
        model << (Nand((var1 == var2, var2 == var3)))
        
        solver = Solver(model)
        
        solved = solver.solve()
        
        assert ( not solved )
        
        
        #print [var.get_value() for var in (var1, var2, var3)]
        
    def testNand_3(self):
        var1 = Variable(1, 1)
        var2 = Variable(0, 2)
        var3 = Variable(1, 1)
        
        model = NativeModel()
        
        model << (Nand((var1 == var2, var1 == var3)))
        
        solver = Solver(model)
        
        assert ( solver.solve() )
        
        #print [var.get_value() for var in (var1, var2, var3)]
        
    def testNot(self):
        var1 = Variable(0, 1)
        var2 = Variable(1, 2)
  
        model = NativeModel()
        
        model << (Not(var1 != var2))
        
        solver = Solver(model)
        
        assert(solver.solve())
        
        #print [var.get_value() for var in (var1, var2)]
        
    def testNotUnsat(self):
        var1 = Variable(0, 1)
        var2 = Variable(2, 3)
  
        model = NativeModel()
        
        model << (Not(var1 != var2))
        
        solver = Solver(model)
        
        assert(not solver.solve())
        
        #print [var.get_value() for var in (var1, var2)]
        
    def testOr_1(self):
        var1 = Variable(0, 1)
        dur1 = 1
        
        var2 = Variable(0, 3)
        dur2 = 2
        
        model = NativeModel()
        
        model.add( ( var1 + dur1 < var2 ) | ( var2 + dur2 < var1 ) )
        
        solver = Solver(model)
        
        assert(solver.solve())
        
        
        
        assert ( (var1.get_value() + dur1 < var2.get_value ) or
            ( var2.get_value() + dur2 < var1.get_value  ) )
        
    def testOr_2(self):
        var1 = Variable(0, 2)
        dur1 = 1
        
        var2 = Variable(0, 4)
        dur2 = 2
        
        var3 = Variable(0, 6)
        dur3 = 3
        
        model = NativeModel()
        
        #model.add( NoOverlap( (var1, var2), (dur1, dur2) ) )
        #model.add( NoOverlap( (var1, var3), (dur1, dur3) ) )
        #model.add( NoOverlap( (var2, var3), (dur2, dur3) ) )
        
        ur = UnaryResource ( [ NoOverlap( (var1, var2), (dur1, dur2) ),
                                     NoOverlap( (var1, var3), (dur1, dur3) ),
                                     NoOverlap( (var2, var3), (dur2, dur3) ) ] )
        
        model.add( ur )
        
        solver = Solver(model)
        
        assert(solver.solve())
        
        assert ( (var1.get_value() + dur1 < var2.get_value ) or
            ( var2.get_value() + dur2 < var1.get_value  ) )
        
    def testAnd_ScalarProd(self):
        vars1 = [Variable(0,1) for x in range(0,5)]
        vars2 = [Variable(0,1) for x in range(0,5)]
        vars3 = [Variable(0,1) for x in range(0,5)]
        vars4 = [Variable(0,1) for x in range(0,5)]
        
        model = NativeModel()
        
        model.add(Sum([ vars1[i] & vars2[i] for i in range(0,5)]) ==
                  Sum([ vars3[i] & vars4[i] for i in range(0,5)]))
        
        model.add( Sum([ vars3[i] & vars4[i] for i in range(0,5)]) == 2 )
        
        solver = Solver(model)
        
        assert( solver.solve() )
        
        '''
        print ""
        print [v.get_value() for v in vars1]
        print [v.get_value() for v in vars2]
        print [v.get_value() for v in vars3]
        print [v.get_value() for v in vars4]
        '''
        
    def testSmallBibd(self):
        v, b, r, k, l = [7,7,3,3,1]
        # This test uses the modelling layer for ease of use
        import Numberjack
        model = Numberjack.Model() # Create the model
        matrix = Numberjack.Matrix(v,b)

        # every row adds up to k
        model.add( [Numberjack.Sum(row) == k for row in matrix.row] )

        # every column adds up to r
        model.add( [Numberjack.Sum(col) == r for col in matrix.col] )

        # the scalar product of every pair of columns adds up to l
    
        model.add( [ Numberjack.Sum([ (row[col_i] & row[col_j]) for row in matrix.row]) == l  
                    for col_i in range(v) for col_j in range(col_i) ] )
        
        solver = Solver(model, [var for col in matrix.col for var in col])
        
        assert(solver.solve())
        
        
    def testAnd_1(self):
        var1 = Variable(0, 1)
        var2 = Variable(1, 2)
        var3 = Variable(0, 1)
        
        model = NativeModel()
        
        model << (And((var1 != var2, var2 != var3)))
        
        solver = Solver(model)
        
        assert ( solver.solve() )
        
        # print [var.get_value() for var in (var1, var2, var3)]
        
    def testAnd_2(self):
        var1 = Variable(0, 1)
        var2 = Variable(0, 1)
        var3 = Variable(1, 1)
        
        model = NativeModel()
        
        model << (And((var1 == var2, var2 == var3)))
        
        solver = Solver(model)
        
        assert ( solver.solve() )
        
        # print [var.get_value() for var in (var1, var2, var3)]
        
    def testAnd_3(self):
        var1 = Variable(1, 1)
        var2 = Variable(2, 2)
        var3 = Variable(1, 1)
        
        model = NativeModel()
        
        model << (And((var1 == var2, var1 == var3)))
        
        solver = Solver(model)
        
        assert ( not solver.solve() )
        
        # print [var.get_value() for var in (var1, var2, var3)]
        
    def testOr(self):
        var1 = Variable(0, 1)
        var2 = Variable(1, 1)
        var3 = Variable(1, 1)
        
        model = NativeModel()
        
        model << (Or((var1 != var2, var2 != var3)))
        
        solver = Solver(model)
        
        assert ( solver.solve() )
        
        # print [var.get_value() for var in (var1, var2, var3)]
        
    def testOr2(self):
        var1 = Variable(0, 1)
        var2 = Variable(1, 2)
        var3 = Variable(0, 1)
        
        model = NativeModel()
        
        model << (Or((var1 == var2, var2 != var3)))
        
        solver = Solver(model)
        
        assert ( solver.solve() )
        
        # print [var.get_value() for var in (var1, var2, var3)]
        
    def testOr3(self):
        var1 = Variable(0, 0)
        var2 = Variable(0, 0)
        var3 = Variable(0, 0)
        
        model = NativeModel()
        
        model << (And((var1 != var2, var2 != var3)))
        
        solver = Solver(model)
        
        assert ( not solver.solve() )
        
        # print [var.get_value() for var in (var1, var2, var3)]
        

    def testFeasTuple_1(self):
        var1 = Variable(0, 0)
        var2 = Variable(0, 0)
        
        tups = ((0,0),)
        
        model = NativeModel()
        
        model.add(Table((var1, var2), tups))
        
        solver = Solver(model)
        
        assert(solver.solve())
        
    def testFeasTuple_2(self):
        var1 = Variable(1, 1)
        var2 = Variable(0, 1)
        
        tups = ((1,0), (0,1))
        
        model = NativeModel()
        
        model.add(Table((var1, var2), tups))
        
        solver = Solver(model)
        
        assert(solver.solve())
        
        #print [var.get_value() for var in (var1, var2)]
        
    def testFeasTuple_3(self):
        var1 = Variable(0, 5)
        var2 = Variable(0, 5)
        
        tups = ((0,6),)
        
        model = NativeModel()
        
        model.add(Table((var1, var2), tups))
        
        solver = Solver(model)
        
        assert(not solver.solve())
        
    def testTable_4(self):
        tuples = [(1,1,2),
          (1,1,3),
          (1,2,1),
          (1,2,2),
          (2,1,1),
          (2,1,2),
          (2,2,1),
          (2,3,1)]
        
        X = [Variable(1,3) for i in range(7)]
        
        model = NativeModel()

        model.add( Table( [X[0], X[1], X[2]], tuples ) )
        model.add( Table( [X[1], X[2], X[3]], tuples ) )
        model.add( Table( [X[2], X[3], X[4]], tuples ) )  
        model.add( Table( [X[3], X[4], X[5]], tuples ) )   
        model.add( Table( [X[4], X[5], X[6]], tuples ) )

        solver = Solver(model)
        
        assert( solver.solve() )

        #print [var.get_value() for var in X] 
        
    def testSingleProp_1(self):
        var1 = Variable(0, 0)
        var2 = Variable(0, 5)
        
        model = NativeModel()
        
        model.add(var2 != var1)
         
        solver = Solver(model)

        solver.propogate()
        
        assert( not var2.can_be_instantiated_to(0) )
        
    def testSingleProp(self):
        var1 = Variable(0, 3)
        var2 = Variable(0, 5)
        
        model = NativeModel()
        
        model.add(var1 > var2)
         
        solver = Solver(model)

        solver.propogate()
        
        #var1.print_domain()
        #var2.print_domain()
        
    def testDecVars(self):
        var1 = Variable(0, 3)
        var2 = Variable(0, 5)
        var3 = Variable(1,2)
        
        model = NativeModel()
        
        model.add(var1 == var2)
        model.add(var1 != var3)
         
        solver = Solver(model, [var1, var2])
        
        assert(solver.solve())
        assert(not var3.get_is_instantiated())
        
    def testHeurSel(self):
        
        var1 = Variable(0, 3)
        var2 = Variable(0, 5)
        var3 = Variable(1,2)
        
        model = NativeModel()
        
        model.add(var1 == var2)
        model.add(var1 != var3)
         
        solver = Solver(model)
        
        solver.setHeuristic("RandomVar", "RandomVal")
        
        assert(solver.solve())
        
    def testElement_1(self):
        var1 = Variable(0, 3)
        var2 = Variable(0, 3)
        var3 = Variable(3, 3)
    
        varI = Variable(0, 5)
        varK = Variable(2,2)
        
        model = NativeModel()
        
        model.add(Element([var1, var2, var3, varK]) == varI)
        
        solver = Solver(model)
        
        assert (solver.solve() )
        
        assert([var1, var2, var3][varK.get_value()] == varI.get_value())

if __name__ == "__main__":
    unittest.main()
    
