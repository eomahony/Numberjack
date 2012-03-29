from Numberjack import *


class MIPParser(object):
    def __init__(self, model_path, data_path=None, solver_name=None):
        self.model_path = model_path
        self.data_path = data_path
        self.solver_name = solver_name
        self.model = None
        self.vars = {}

    def parse_model(self):
        # Import Glpk solver so that we have support for gmpl
        import OsiGlpk as Osi

        parts = self.model_path.split('.')
        extensions = parts[-2:len(parts)]
        filetype = extensions[0] if extensions[1] == 'gz' else extensions[1]

        parser = Osi.Solver(Model())
        if filetype == 'mod':
            parser.load_gmpl(self.model_path, self.data_path)
        elif filetype == 'mps':
            parser.load_mps(self.model_path, 'mps')
        elif filetype == 'lp':
            parser.load_lp(self.model_path, 0)
        else:
            print "Unknown filetype, Supported formats are [gmpl, mps, mps.gz, lp], exiting"
            return None
        prse = parser.solver

        self.model = Model()

        n = prse.num_expression()
        for i in range(n):
            expr = self.getNJExp(prse.get_expression(i), 0)
            self.model.add(expr)
        return self.model

    def solve(self):
        if self.model == None:
            self.parse_model()
        solver = self.model.load(self.solver_name)
        return (solver.solve(), self.model)

    def getNJExp(self, mexp, ident):
        exp_type = mexp.get_type()
        if exp_type == 'var':
            return self.getNJVar(mexp, ident)
        else:
            return self.getNJPred(mexp)

    def ifwhole(self, val):
        return val if val % 1 != 0 else int(val)

    def getNJPred(self, mexp):
        exp_type = mexp.get_type()

        children = []
        for i in range(mexp.get_arity()):
            children.append(self.getNJExp(mexp.get_child(i), 0))

        if exp_type == "le":
            return children[0] <= self.ifwhole(mexp.get_parameter(0))
        elif exp_type == "ge":
            return children[0] >= self.ifwhole(mexp.get_parameter(0))
        elif exp_type == "sum":
            weights = []
            for i in range(mexp.get_arity()):
                weights.append(self.ifwhole(mexp.get_parameter(i)))
            return Sum(children, weights)
        elif exp_type == "minimise":
            return Minimise(children[0])
        else:
            print "Error: Failed to parse expression:", type
            exit(1)

    def getNJVar(self, mexp, ident):
        ident = mexp.getVariableId()
        if ident not in self.vars:
            lb = mexp.get_min() if mexp.is_continuous() else int(mexp.get_min)
            ub = mexp.get_max() if mexp.is_continuous() else int(mexp.get_max)
            self.vars[ident] = Variable(lb, ub, mexp.name())
        return self.vars[ident]


if __name__ == "__main__":
    import os
    import sys
    if len(sys.argv) >= 2:
        filename = sys.argv[1]
        datafile = sys.argv[2] if len(sys.argv) == 3 else None

        parser = MIPParser(filename, datafile, 'OsiCbc')
        (solved, model) = parser.solve()
        if solved:
            print [(v.name(), v.get_value()) for v in model.variables]
        else:
            print "UNSAT"