#! /usr/bin/env python
import datetime
from Numberjack import *
from Numberjack.Decomp import PostBinary, PostUnary, PostTernary


# library of flatzinc predicates translated into numberjack constraints


def array_bool_and(x,y):
    if ((type(y) is int) and y != 0):
        return (Sum(x) == len(x))
    elif ((type(y) is int) and y == 0):
        return (Sum(x) != len(x))
    else:
        return (y == (Sum(x) == len(x)))


def array_bool_or(x,y):
    if ((type(y) is int) and y != 0):
        if len(x) == 2:
            return Or(x)
        else:
            return Disjunction(x)
    else:
        return (y == Disjunction(x))

def array_bool_xor(x):
    return ((Sum(x) % 2) == 1)

# cost functions direct access through MiniZinc
# only feasible with toulbar2 
def cost_function_unary(var, costs, costVar):
    return PostUnary(var,costs)

def cost_function_binary(var1, var2, costs, costVar):
    return PostBinary(var1, var2, costs)

def cost_function_ternary(var1, var2, var3, costs, costVar):
    return PostTernary(var1, var2, var3, costs)


def array_int_element(x, y, z):
    # Buggy Workaround, produces invalid values in some optimization cases.
    #aux = Variable(x.lb-1, x.ub-1, "somevar_minus1")
    #return [(z == Element([Variable(e,e,str(e)) if type(e) is int else e for e in y], aux)), (x >= 1), (x <= len(y)), (aux == x - 1)]
    #return [(z == Element([Variable(e,e,str(e)) if type(e) is int else e for e in y], x - 1)), int_le(1,x), int_le(x,len(y))]
    u = set()
    for e in y:
        u = u | set([e] if type(e) is int else range(e.lb, e.ub + 1))
    return [int_le(1,x), int_le(x,len(y)), set_in(z, u)] + [((z == e) | (x != i+1)) for i, e in enumerate(y)]

def array_var_int_element(x,y,z):
    return (array_int_element(x,y,z))

def array_bool_element(x,y,z):
    return (array_int_element(x,y,z))

def array_var_bool_element(x,y,z):
    return (array_var_int_element(x,y,z))

def bool2int(x, y):
    return (x == y)

def bool_and(x, y, z):
    return (And(x, y) if ((type(z) is int) and (z != 0)) else (z == And(x, y)))

def bool_clause(x, y):
    if len(x)>0 and len(y)==0:
      return Disjunction(x)
    elif len(x)==0 and len(y)>0:
      return Disjunction([(e == 0) for e in y])
    elif len(x)>0 and len(y)>0:
      return (Disjunction(x) | Disjunction([(e == 0) for e in y]))
    else:
      return []

def bool_le(x, y):
    return ((x == 0) | (y != 0))

def bool_le_reif(x, y, z):
    return [((x != 0) | (z != 0)), ((y != 0) | (z != 0)), ((x == 0) | (y != 0) | (z == 0))]

def bool_lt(x, y):
    return [(x == 0), (y != 0)]

def bool_lt_reif(x, y, z):
    return [((x == 0) | (z == 0)), ((y != 0) | (z == 0)), ((x != 0) | (y == 0) | (z != 0))]

def bool_not(x, y):
    return [((x == 0) | (y == 0)), ((x != 0) | (y != 0))]

def bool_or(x, y, z):
    return (z == (x | y ))

def bool_xor(x, y, z):
    return (z == (x != y))

def int_eq(x,y):
        '''
    if (type(y) is int) and issubclass(type(x), Expression) and x.is_var() and y >= x.lb and y <= x.ub:
        x.domain_ = None
        x.lb = y
        x.ub = y
        return []
    else:
        '''
        return (x == y)

def int_eq_reif(x,y,z):
    return (z == (x == y))

def bool_eq(x, y):
    return (int_eq(x,y))

def bool_eq_reif(x, y, z):
    return (int_eq_reif(x, y, z))

def int_le(x,y):
        '''
    if (type(y) is int) and issubclass(type(x), Expression) and x.is_var() and y >= x.lb: # and (x.domain_ is None)
        x.ub = min(x.ub, y)
        return []
    elif (type(x) is int) and issubclass(type(y), Expression) and y.is_var() and x <= y.ub: # and (y.domain_ is None)
        y.lb = max(y.lb, x)
        return []
    else:
        '''
        return (x <= y)

def int_le_reif(x,y,z):
    return (z == (x <= y))

def int_lt(x,y):
        '''
    if (type(y) is int) and issubclass(type(x), Expression) and x.is_var() and y > x.lb: # and (x.domain_ is None)
        x.ub = min(x.ub, y-1)
        return []
    elif (type(x) is int) and issubclass(type(y), Expression) and y.is_var() and x < y.ub: # and (y.domain_ is None)
        y.lb = max(y.lb, x+1)
        return []
    else:
        '''
        return (x < y)

def int_lt_reif(x,y,z):
    return (z == (x < y))

def int_ne(x,y):
    return (x != y)

def int_ne_reif(x,y,z):
    return (z == (x != y))

def int_lin_eq(coef,vars,res):
    if ((len(coef) == 1) and (coef[0] == 1)):
       return int_eq(vars[0],res)
    else:
       return (res == Sum(vars,coef))

def bool_lin_eq(coef,vars,res):
    return (int_lin_eq(coef,vars,res))

def int_lin_eq_reif(coef,vars,res,z):
    if ((len(coef) == 1) and (coef[0] == 1)):
       return int_eq_reif(vars[0],res,z)
    else:
       return (z == (res == Sum(vars, coef)))

def int_lin_le(coef,vars,res):
    if ((len(coef) == 1) and (coef[0] == 1)):
       return int_le(vars[0],res)
    else:
       return (res >= Sum(vars,coef))

def bool_lin_le(coef,vars,res):
    return (int_lin_le(coef,vars,res))

def int_lin_le_reif(coef,vars,res,z):
    if ((len(coef) == 1) and (coef[0] == 1)):
       return int_le_reif(vars[0],res,z)
    else:
       return (z == (res >= Sum(vars,coef)))

def int_lin_lt(coef,vars,res):
    if ((len(coef) == 1) and (coef[0] == 1)):
       return int_lt(vars[0],res)
    else:
       return (res > Sum(vars,coef))

def int_lin_lt_reif(coef,vars,res,z):
    if ((len(coef) == 1) and (coef[0] == 1)):
       return int_lt_reif(vars[0],res,z)
    else:
       return (z == (res > Sum(vars,coef)))

def int_lin_ne(coef,vars,res):
    if ((len(coef) == 1) and (coef[0] == 1)):
       return int_ne(vars[0],res)
    else:
       return (res != Sum(vars,coef))

def int_lin_ne_reif(coef,vars,res,z):
    if ((len(coef) == 1) and (coef[0] == 1)):
       return int_ne_reif(vars[0],res,z)
    else:
       return (z == (res != Sum(vars,coef)))

def int_abs(x,y):
    return (y == Abs(x))

def int_div(x,y,z):
    if (x is y): return (z == 1)
    return (z == (x / y))

def int_min(x,y,z):
    if (x is y): return (z == x)
    return ((z == Min([Variable(x,x,str(x)),y])) if (type(x) is int) else ((z == Min([x,Variable(y,y,str(y))])) if (type(y) is int) else (z == Min([x,y]))))

def int_max(x,y,z):
    if (x is y): return (z == x)
    return ((z == Max([Variable(x,x,str(x)),y])) if (type(x) is int) else ((z == Max([x,Variable(y,y,str(y))])) if (type(y) is int) else (z == Max([x,y]))))

def int_mod(x,y,z):
    if (x is y): return (z == 0)
    return (z == (x % y))

def int_plus(x,y,z):
    return (z == (x + y))

def int_times(x,y,z):
    return (z == (x * y))

def set_in(x,dom):
    if (type(x) is int):
        if (not(x in dom)):
            return [Variable(x,x,str(x)) != x]
        else:
            return []
#    return (Disjunction([(x == v) for v in dom]))
    return [(x != v) for v in range(x.get_min(),1+x.get_max()) if (not(v in dom))]
'''
    if (x.domain_ is not None):
        olddom = set(x.domain_)
        newdom = set(dom)
        x.domain_ = list(newdom & olddom)
    else:
        x.domain_ = [e for e in dom if e >= x.lb and e <= x.ub]
    x.domain_.sort()
    x.lb = x.domain_[0]
    x.ub = x.domain_[-1]
    return []
'''

def set_in_reif(x,dom,z):
    return (z == Disjunction([(x == v) for v in dom]))

# specific global constraints for numberjack

def all_different_int(x):
#    x = set(x_)
    if len(x) < 2:  # Some models specified alldiff on 1 variable
        return x
    return (AllDiff([Variable(e,e,str(e)) if type(e) is int else e for e in x]))

def lex_less_int(x,y):
    if len(x) == 1 and len(y) == 1:
        return x[0] < y[0]
    return LessLex(x, y)

def lex_lesseq_int(x,y):
    if len(x) == 1 and len(y) == 1:
        return x[0] <= y[0]
    return LeqLex(x, y)

def lex_less_bool(x,y):
    return (lex_less_int(x,y))

def lex_lesseq_bool(x,y):
    return (lex_lesseq_int(x,y))

def minimum_int(x,y):
#    y = set(y_)
    if(len(y)==1):
        return (x == y[0])
    else:
        return (x == Min(y))

def maximum_int(x,y):
#    y = set(y_)
    if(len(y)==1):
        return (x == y[0])
    else:
        return (x == Max(y))

def table_int(x,t):
    return (Table([Variable(e,e,str(e)) if type(e) is int else e for e in x],[tuple([t[i * len(x) + j] for j in range(len(x))]) for i in range(len(t) / len(x))]))

def table_bool(x,t):
    return (table_int(x, t))


def total_seconds(td):
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 1e6) / 1e6


def time_remaining(tcutoff):
    return max(tcutoff - total_seconds(datetime.datetime.now() - start_time), 0.0)


def run_solve(model, output_vars, param):
    # make it fail on purpose
 #   x = Variable(-1000000, 10000000, 'x')
#    model.add(x == model.variables[0])
#    param['solver'] = 'Mistral2'
    load_time = datetime.datetime.now()
    encoding = NJEncodings[param['encoding']] if param['encoding'] else None
    solver = model.load(param['solver'], encoding=encoding)
    solver.setVerbosity(param['verbose'])
    time_limit = max(int(param['tcutoff'] - total_seconds(datetime.datetime.now() - load_time)), 1)
    solver.setTimeLimit(time_limit)
    solver.setHeuristic(param['var'], param['val'], param['rand'])
    solver.setThreadCount(param['threads'])
    if param['solver'] == 'Toulbar2':
        solver.setOption('lds',param['lds'])
        solver.setOption('lcLevel',param['lcLevel'])
        solver.setOption('deadEndElimination',param['dee'])
        solver.setOption('btdMode',param['btd'])
        solver.setOption('splitClusterMaxSize',param['rds'])
##        uncomment the following lines to save the problem in wcsp format
#        solver.setOption('nopre')
#        solver.setOption('lcLevel',0)
#        solver.setOption("dumpWCSP",2)
    if param['solver'] == 'Mistral':
        solver.solveAndRestart(param['restart'], param['base'], param['factor'])
    else:
        solver.solve()
    return solver, output_vars


def solve_main(param):
    model, output_vars = get_model()
    return run_solve(model, output_vars, param)
