from __future__ import print_function, division
from Numberjack import *
import xml.etree.cElementTree as ET
import sys


sys.setrecursionlimit(100000)


# Mapping of the XCSP functional notation to Numberjack predicates.
# Have to use lambda with the operators instead of Mul etc so that in the case
# of (2 * x), Mul([x, 2]) gets created instead of Mul([2, x]).
functional_map = {
    "abs": Abs,
    "add": lambda x_y: x_y[0] + x_y[1],
    "and": And,
    "div": lambda x_y1: x_y1[0] / x_y1[1],
    "eq": lambda x_y2: x_y2[0] == x_y2[1],
    "ge": lambda x_y3: x_y3[0] >= x_y3[1],
    "gt": lambda x_y4: x_y4[0] > x_y4[1],
    "le": lambda x_y5: x_y5[0] <= x_y5[1],
    "lt": lambda x_y6: x_y6[0] < x_y6[1],
    "max": Max,
    "min": Min,
    "mod": Mod,
    "mul": lambda x_y7: x_y7[0] * x_y7[1],
    "ne": lambda x_y8: x_y8[0] != x_y8[1],
    "neg": Neg,
    # "not": Not,
    "or": Or,
    "sub": lambda x_y9: x_y9[0] - x_y9[1],
    "iff": lambda x_y10: Or([x_y10[0] == 0, x_y10[1]]),
}

global_map = {
    "global:allDifferent": AllDiff,
    "global:element": lambda i_X_v: i_X_v[2] == Element(i_X_v[1], i_X_v[0]),
    "global:weightedSum": lambda X_W: Sum(X_W[0], X_W[1]),
}


class XCSPPredicate(object):

    def __init__(self, parameter_order, pred_str):
        self.parameter_order = parameter_order
        self.predicate = None
        self.children = []
        self.parse_functional(pred_str)

    def parse_functional(self, pred_str):
        # Split out the predicate name
        pred_str = pred_str.strip()
        bits = pred_str.split("(", 1)
        if len(bits) == 1:
            raise XCSPParserError("Expected predicate %s" % pred_str)
        pred_name = bits[0]
        remaining = bits[1][:-1]
        if pred_name not in functional_map:
            raise XCSPParserUnsupportedError("Unknown predicate name %s" % pred_name)

        self.pred_name = pred_name
        self.predicate = functional_map[pred_name]

        # Parse arguments
        def parse_arg(arg_str):
            if "(" in arg_str:
                return XCSPPredicate(self.parameter_order, arg_str)
            return arg_str

        def split_children(s):
            level = 0
            current = ""
            for c in s:
                if c == "(":
                    level += 1
                    current += c
                elif c == ")":
                    level -= 1
                    current += c
                elif c == "," and level == 0:
                    yield current
                    current = ""
                else:
                    current += c
            if level != 0:
                raise XCSPParserError("Error parsing predicate children, unbalanced parenthesis. %s" % s)
            yield current

        self.children = [parse_arg(arg) for arg in split_children(remaining)]

    def get_expr(self, args):
        arg_map = dict(list(zip(self.parameter_order, args)))
        return self.get_expr_with_arg_map(arg_map)

    def get_expr_with_arg_map(self, arg_map):
        def get_child(c):
            if isinstance(c, XCSPPredicate):
                return c.get_expr_with_arg_map(arg_map)
            elif isnumeric(c):
                return int(c)
            elif isinstance(c, str):
                try:
                    return arg_map[c]
                except Exception as e:
                    print(self.pred_name, self.predicate, self.children, self.parameter_order)
                    print("\n".join("%s:%s" % (k, str(v)) for k, v in arg_map.items()))
                    print(str(e))
                    raise e
            else:
                raise XCSPParserError("Error unknown child %s" % str(c))

        return self.predicate([get_child(c) for c in self.children])


class XCSPRelation(object):

    def __init__(self, semantics):
        self.semantics = "support" if "support" in semantics else "conflict"
        self.tuples = []

    def get_expr(self, args):
        if len(args) == 1:  # Benchmarks like si2-BVG define a Table over one variable, for some reason.
            x = args[0]
            if self.semantics == "conflict":
                return [x != t[0] for t in self.tuples]
            else:
                return Disjunction([x == t[0] for t in self.tuples])
        return Table(args, self.tuples, type=self.semantics)


class XCSPParserError(exceptions.Exception):

    def __init__(self, msg="XCSP parser error."):
        self.msg = msg

    def __str__(self):
        return "Error: %s" % self.msg


class XCSPParserUnsupportedError(XCSPParserError):

    def __init__(self, msg="Unsupported constraint."):
        self.msg = msg

    def __str__(self):
        return "Error: %s" % self.msg


class XCSPParser(object):

    def __init__(self, filename):
        self.filename = filename
        self.model = Model()
        self.tree = ET.parse(filename)
        self.root = self.tree.getroot()
        self.domains = {}
        self.variables = []
        self.variable_map = {}
        self.pred_and_rel = {}

        self.parse_domains()
        self.parse_variables()
        self.parse_relations()
        self.parse_predicates()
        self.parse_constraints()

    def parse_domains(self):

        def parse_domain_string(s):
            domain = []
            bits = s.strip().split(" ")
            for bit in bits:
                if ".." in bit:
                    l, u = list(map(int, bit.split("..")))
                    if len(bits) == 1:
                        # Domain is just specified by a single range
                        return [l, u]
                    else:
                        domain.extend(list(range(l, u + 1)))
                else:
                    domain.append(int(bit))
            return [domain]

        domains = self.root.find('domains')
        for d in domains:
            name = d.attrib['name']
            self.domains[name] = parse_domain_string(d.text)
            # print name, self.domains[name]

    def parse_variables(self):
        variables = self.root.find('variables')
        for v in variables:
            name = v.attrib['name']
            domain = self.domains[v.attrib['domain']]
            args = domain + [name]
            var = Variable(*args)
            self.variables.append(var)
            self.variable_map[name] = var

    def parse_relations(self):
        relations = self.root.find('relations')
        if relations:
            for r in relations:
                name = r.attrib['name']
                semantics = r.attrib['semantics']

                relation = XCSPRelation(semantics)
                if len(list(r)) > 0:
                    raise XCSPParserUnsupportedError("Only abridged notation for relations is supported for now.")
                if r.text is not None and len(r.text.strip()) > 0:
                    for t_str in r.text.split("|"):
                        bits = [s.strip() for s in t_str.split(" ") if len(s.strip()) > 0]
                        relation.tuples.append(list(map(int, bits)))
                # print semantics, relation.tuples
                self.pred_and_rel[name] = relation

    def parse_predicates(self):
        predicates = self.root.find('predicates')
        if predicates:
            for p in predicates:
                name = p.attrib['name']
                parameters = p.find('parameters')
                param_order = parameters.text.replace("int", "").split()
                expression = p.find('expression')
                functional = expression.find('functional')
                if functional is None:
                    raise XCSPParserUnsupportedError("Only functional predicate definitions are currently supported.")
                self.pred_and_rel[name] = XCSPPredicate(param_order, functional.text)

    def parse_constraints(self):
        def build_param(param):
            if isinstance(param, str):
                if isnumeric(param):
                    return int(param)
                return self.variable_map[param]
            else:
                raise XCSPParserError("Unknown parameter" % param)

        def parse_elem_param(s):
            s = s.strip()
            ind1 = s.find("[")
            ind2 = s.rfind("]")
            ind = s[:ind1].strip()
            xs = s[ind1 + 1:ind2].strip().split(" ")
            v = s[ind2 + 1:].strip()
            return ind, xs, v

        def parse_wsum_param(s):
            import re
            s = s.strip()
            W, X = [], []
            exp = re.compile(r"\{\s*(?P<coef>[-]?\d+)\s+(?P<varname>[\w\._]+)\s*\}")
            for match in exp.finditer(s):
                d = match.groupdict()
                W.append(int(d['coef']) if 'coef' in d else 1)
                X.append(str(d['varname']))
            return W, X

        constraints = self.root.find('constraints')
        for c in constraints:
            constraint = None
            reference = c.attrib['reference']

            if reference in self.pred_and_rel:
                predicate = self.pred_and_rel[reference]
                parameters = c.find('parameters')
                parameter_str = parameters.text if parameters is not None else c.attrib['scope']
                args = [build_param(p) for p in parameter_str.split(" ") if len(p.strip()) > 0]
                try:
                    constraint = predicate.get_expr(args)
                except Exception as e:
                    print(reference, parameter_str, args)
                    raise e

            elif "global:" in reference:
                if reference not in global_map:
                    raise XCSPParserUnsupportedError("Unknown global %s" % reference)

                if reference == "global:allDifferent":
                    parameters = c.find('parameters')
                    parameter_str = parameters.text if parameters else c.attrib['scope']
                    args = [build_param(p) for p in parameter_str.split(" ") if len(p.strip()) > 0]
                    constraint = global_map[reference](args)

                elif reference == "global:element":
                    parameter_str = c.find('parameters').text
                    ind, X, v = parse_elem_param(parameter_str)
                    X = [build_param(x) for x in X]
                    X = VarArray([Variable([x], str(x)) if isinstance(x, int) else x for x in X])  # Mistral requires a VarArray for Element
                    ind = build_param(ind)
                    v = build_param(v)

                    # XCSP indexes are 1-based :(
                    constraint = global_map[reference]([ind - 1, X, v])

                elif reference == "global:weightedSum":
                    parameters = c.find('parameters')
                    op_str = parameters[0].tag
                    lhs, rhs = list(parameters.itertext())
                    W, var_names = parse_wsum_param(lhs)
                    X = [build_param(x) for x in var_names]
                    rhs = build_param(rhs)
                    constraint = functional_map[op_str]([Sum(X, W), rhs])

                else:
                    raise XCSPParserUnsupportedError("Unknown global constraint %s" % reference)

            else:
                raise XCSPParserUnsupportedError("Unknown constraint predicate/relation: %s" % reference)

            self.model += constraint


def isnumeric(s):
    try:
        int(s)
    except (ValueError, TypeError):
        return False
    return True


if __name__ == '__main__':
    import datetime
    import time
    import os

    def usage():
        print("Usage: python %s -solver Mistral|MiniSat|... -xcsp xcspfilename.xml" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    default = {'solver': '', 'verbose': 0, 'tcutoff': 3600, 'xcsp': '', 'encoding': ''}
    param = input(default)
    filename = os.path.abspath(param['xcsp'])
    encoding = NJEncodings[param['encoding']] if param['encoding'] else None

    if not os.path.isfile(filename):
        print("Error: the file '%s' does not exist." % filename, file=sys.stderr)
        usage()
    if not param['solver']:
        print("Error: Please sepcify a solver.", file=sys.stderr)
        usage()

    t = datetime.datetime.now()
    c = time.clock()
    parser = XCSPParser(filename)
    print("c Time to parse: %.2f %.2f" % ((datetime.datetime.now() - t).total_seconds(), (time.clock() - c)))
    model, variables = parser.model, parser.variables
    # print "\n".join(str(v) for v in variables)
    # print model
    # sys.exit(0)
    # t = datetime.datetime.now()
    # model.preprocess()
    # print "c Time to preprocess: %.2f" % (datetime.datetime.now() - t).total_seconds()
    print("c Loading model")
    t = datetime.datetime.now()
    c = time.clock()
    s = model.load(param['solver'], encoding=encoding)
    print("c Time to load model: %.2f %.2f" % ((datetime.datetime.now() - t).total_seconds(), (time.clock() - c)))
    s.setVerbosity(param['verbose'])
    s.setTimeLimit(int(param['tcutoff'] - time.clock()))
    print("c Solve")
    t = datetime.datetime.now()
    c = time.clock()
    s.solve()
    print("c Time to solve: %.2f %.2f" % ((datetime.datetime.now() - t).total_seconds(), (time.clock() - c)))
    print("c Nodes %d" % s.getNodes())
    print("c Failures %d" % s.getFailures())
    print("c SolveTime %.4f" % s.getTime())
    if s.is_sat():
        print("s SATISFIABLE")
        print("v", " ".join(str(v.get_value()) for v in variables))
    elif s.is_unsat():
        print("s UNSATISFIABLE")
    else:
        print("s UNKNOWN")
