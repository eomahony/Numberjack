#! /usr/bin/env python
import datetime
from Numberjack import *


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


def array_int_element(x, y, z):
    # Buggy Workaround, produces invalid values in some optimization cases.
    # aux = Variable(x.lb-1, x.ub-1, "somevar_minus1")
    # return [(z == Element([Variable(e,e,str(e)) if type(e) is int else e for e in y], aux)), (x >= 1), (x <= len(y)), (aux == x - 1)]

    u = set()
    for e in y:
        u = u | set([e] if type(e) is int else range(e.lb, e.ub + 1))
    return [(x >= 1), (x <= len(y)), set_in(z, u)] + [((z == (Variable(e,e,str(e)) if type(e) is int else e)) | (x != i+1)) for i, e in enumerate(y)]

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
    return (Disjunction(x) | Disjunction([(e == 0) for e in y]))

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
    return (x == y)

def int_eq_reif(x,y,z):
    return [((x != y) | (z != 0)), ((x == y) | (z == 0))]

def bool_eq(x, y):
    return (int_eq(x,y))

def bool_eq_reif(x, y, z):
    return (int_eq_reif(x, y, z))

def int_le(x,y):
    return (x <= y)

def int_le_reif(x,y,z):
    return [((x > y) | (z != 0)), ((x <= y) | (z == 0))]

def int_lt(x,y):
    return (x < y)

def int_lt_reif(x,y,z):
    return [((x >= y) | (z != 0)), ((x < y) | (z == 0))]

def int_ne(x,y):
    return (x != y)

def int_ne_reif(x,y,z):
    return [((x == y) | (z != 0)), ((x != y) | (z == 0))]

def int_lin_eq(coef,vars,res):
    return (res == Sum(vars,coef))

def bool_lin_eq(coef,vars,res):
    return (int_lin_eq(coef,vars,res))

def int_lin_eq_reif(coef,vars,res,z):
    return (z == (res == Sum(vars, coef)))

def int_lin_le(coef,vars,res):
    return (res >= Sum(vars,coef))

def bool_lin_le(coef,vars,res):
    return (int_lin_le(coef,vars,res))

def int_lin_le_reif(coef,vars,res,z):
    return (z == (res >= Sum(vars,coef)))

def int_lin_lt(coef,vars,res):
    return (res > Sum(vars,coef))

def int_lin_lt_reif(coef,vars,res,z):
    return (z == (res > Sum(vars,coef)))

def int_lin_ne(coef,vars,res):
    return (res != Sum(vars,coef))

def int_lin_ne_reif(coef,vars,res,z):
    return (z == (res != Sum(vars,coef)))

def int_abs(x,y):
    return (y == Abs(x))

def int_div(x,y,z):
    return (z == (x / y))

def int_min(x,y,z):
    return ((z == Min([Variable(x,x,str(x)),y])) if (type(x) is int) else ((z == Min([x,Variable(y,y,str(y))])) if (type(y) is int) else (z == Min([x,y]))))

def int_max(x,y,z):
    return ((z == Max([Variable(x,x,str(x)),y])) if (type(x) is int) else ((z == Max([x,Variable(y,y,str(y))])) if (type(y) is int) else (z == Max([x,y]))))

def int_mod(x,y,z):
    return (z == (x % y))

def int_plus(x,y,z):
    return (z == (x + y))

def int_times(x,y,z):
    return (z == (x * y))

def set_in(x,dom):
#    return (Disjunction([(x == v) for v in dom]))
    return [(x != v) for v in range(x.get_min(),1+x.get_max()) if (not(v in dom))]

def set_in_reif(x,dom,z):
    return (z == Disjunction([(x == v) for v in dom]))

# specific global constraints for numberjack

def all_different_int(x):
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
    if(len(y)==1):
        return (x == y[0])
    else:
        return (x == Min(y))

def maximum_int(x,y):
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
    load_time = datetime.datetime.now()
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    time_limit = max(int(param['tcutoff'] - total_seconds(datetime.datetime.now() - load_time)), 1)
    solver.setTimeLimit(time_limit)
    solver.setHeuristic(param['var'], param['val'], param['rand'])
    if param['solver'] == 'Gurobi':
        solver.setThreadCount(param['threads'])
    if param['solver'] == 'Mistral':
        solver.solveAndRestart(param['restart'], param['base'], param['factor'])
    else:
        solver.solve()
    return solver, output_vars


def solve_main(param):
    model, output_vars = get_model()
    return run_solve(model, output_vars, param)
