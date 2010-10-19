from Numberjack import *
import MipWrapper
import SCIP

# Test <= reif

var1 = Variable(0, 5)
var2 = Variable(6, 10)
var3 = Variable()
model = Model(var3 == (var1 <= var2))
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var3.get_value() == 1)
print "<= reif test1 passed"

var1 = Variable(0, 5)
var2 = Variable(6, 10)
var3 = Variable()
model = Model(var3 == (var1 <= var2),
              var3 < 1)
solver = SCIP.Solver(model)
assert( not solver.solve() )
print "<= reif test2 passed"

var1 = Variable(0, 5)
var2 = Variable(4, 10)
var3 = Variable()
model = Model(var3 == (var1 <= var2),
              var3 < 1)
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var3.get_value() == 0 )
print "<= reif test3 passed"

# Test >= reif

var1 = Variable(0, 5)
var2 = Variable(6, 10)
var3 = Variable()
model = Model(var3 == (var2 >= var1))
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var3.get_value() == 1)
print ">= reif test1 passed"

var1 = Variable(0, 5)
var2 = Variable(6, 10)
var3 = Variable()
model = Model(var3 == (var2 >= var1),
              var3 < 1)
solver = SCIP.Solver(model)
assert( not solver.solve() )
print ">= reif test2 passed"

var1 = Variable(0, 5)
var2 = Variable(4, 10)
var3 = Variable()
model = Model(var3 == (var2 >= var1),
              var3 < 1)
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var3.get_value() == 0 )
print ">= reif test3 passed"

# Test < reif

var1 = Variable(0, 5)
var2 = Variable(6, 10)
var3 = Variable()
model = Model(var3 == (var1 < var2))
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var3.get_value() == 1)
print "< reif test1 passed"

var1 = Variable(0, 5)
var2 = Variable(6, 10)
var3 = Variable()
model = Model(var3 == (var1 < var2),
              var3 < 1)
solver = SCIP.Solver(model)
assert( not solver.solve() )
print "< reif test2 passed"

var1 = Variable(0, 5)
var2 = Variable(4, 10)
var3 = Variable()
model = Model(var3 == (var1 < var2),
              var3 < 1)
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var3.get_value() == 0 )
print "< reif test3 passed"

# Test > reif

var1 = Variable(0, 5)
var2 = Variable(6, 10)
var3 = Variable()
model = Model(var3 == (var2 > var1))
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var3.get_value() == 1)
print "> reif test1 passed"

var1 = Variable(0, 5)
var2 = Variable(6, 10)
var3 = Variable()
model = Model(var3 == (var2 > var1),
              var3 < 1)
solver = SCIP.Solver(model)
assert( not solver.solve() )
print "> reif test2 passed"

var1 = Variable(0, 5)
var2 = Variable(4, 10)
var3 = Variable()
model = Model(var3 == (var2 > var1),
              var3 < 1)
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var3.get_value() == 0 )
print "> reif test3 passed"

# Testing binary constraints with variables and constants

var1 = Variable(0, 10)
model = Model( var1 < 5 )
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var1.get_vallue() < 5 )
print "< constant passed"

var1 = Variable(0, 10)
model = Model( var1 <= 5 )
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var1.get_vallue() <= 5 )
print "<= constant passed"

var1 = Variable(0, 10)
model = Model( var1 > 5 )
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var1.get_vallue() > 5 )
print "> constant passed"

var1 = Variable(0, 10)
model = Model( var1 >= 5 )
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var1.get_vallue() >= 5 )
print ">= constant passed"

var1 = Variable(0, 10)
model = Model( var1 == 5 )
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var1.get_vallue() == 5 )
print "== constant passed"

var1 = Variable(0, 10)
model = Model( var1 != 5 )
solver = SCIP.Solver(model)
assert( solver.solve() )
assert( var1.get_vallue() != 5 )
print "<= constant passed"






