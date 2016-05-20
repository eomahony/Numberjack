# -*- coding: utf-8 -*-
#"""
#  Numberjack is a constraint satisfaction and optimisation library
#  Copyright (C) 2009 Cork Constraint Computation Center, UCC
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#  The authors can be contacted electronically at
#  numberjack.support@gmail.com

import types

class XCSPOutput():
    """
    This module handles outputting to XCSP 2.1 format
    Should be 2008 CSP Solver Competition rules compliant
    """

    # This maps the Numberjack global names to the XCSP ones
    global_map = {
        "AllDiff" : "allDifferent",
        "among"   : "among",
        "AtLeast" : "atleast",
        "AtMost"  : "atmost",
        #"NOT_SUPPORTED" : "cumulative",  # we need support for O1 E1 H1
        "cycle"   : "cycle",
        #"NOT_SUPPORTED" : "diffn",
        #"NOT_SUPPORTED" : "disjunctive", # we need support for O1 E1 H1
        "Element" : "element",
        "Gcc"     : "global_cardinality",
        "global_cardinality_with_costs" : "global_cardinality_with_costs",
        "minimum_weight_alldifferent" : "minimum_weight_alldifferent",
        "not_all_equal" : "not_all_equal",
        "nvalue"  : "nvalue",
        "nvalues" : "nvalues"
    }
    """
    Layout of the <parameters> section of global constraints
    S1     : Variable or Constant
    VARS   : List of variables and or constants
    VALS   : List of weighted variables
    E1     : Variable or Constant
    RELOP  : Relational operator (<eq/>)
    MATRIX : n*3 matrix [ {1 1 1} {1 1 0} {1 1 3} ]
    ORTHOTOPES: Dunno
    ORTHOTOPE: Dunno
    TASKS: Dunno
    """
    global_format = {
        "AllDiff": ["VARS"],           "among": ["S1","VARS","VALS"],
        "AtLeast": ["S1","VARS","E1"], "AtMost": ["S1","VARS","E1"],
        "cumulative": ["TASKS","E1"],  "cycle": ["S1","VALS"],
        "diffn": ["ORTHOTOPES"],       "disjunctive": ["TASKS"],
        "Element": ["VARS","E1"],      "Gcc": ["VARS","VALS"],
        "global_cardinality_with_costs":["VARS","VALS","MATRIX","E1"],
        "minimum_weight_alldifferent": ["VARS","MATRIX","E1"],
        "not_all_equal": ["VARS"],     "nvalue": ["S1","VARS"],
        "nvalues": ["VARS","RELOP","E1"]
    }

    # This maps integer expressions to boolean equivalents
    int2bool = {
        "eq" : "iff( #a, #b)",
        "ne" : "xor( #a, #b)",
        "ge" : "or( and(not(#a), not(#b)), #a)",
        "gt" : "and( #a, not(#b))",
        "le" : "or( not(#a), #b)",  # Implication
        "lt" : "not( or(#a, #b))"
    }

    def __init__(self, model):
        self.__model = model
        self.__model.close() # Close the model

        self.__domains = []
        self.__variables = []
        self.__relations = []
        self.__predicates = []
        self.__constraints = []
        self.__strings = []

        # These two record predicates and domains already specified
        self.pred_map = {}
        self.dom_map = {}

        # maxConstraintArity attribute must be printed in <presentation />
        self.__maxConstraintArity = 0

        # These are so that intermediate vars can be handled in a normalised form
        self.__local_var_idx = 0
        self.__local_rel_idx = 0
        self.__local_pred_idx = 0
        self.__local_con_idx = 0

        self.njvar_mapping = {}

    def output(self, path):
        """
        Outputs the model to a given file.
        """
        self.load_variables()
        self.load_expressions()

        outfile = open(path, "w")
        outfile.write("<instance>\n" +
                      "<presentation " +
                      "maxConstraintArity=\"%d\" " % self.__maxConstraintArity +
                      "format=\"XCSP 2.1\" />\n")

        self.output_domains(outfile)
        self.output_variables(outfile)
        self.output_relations(outfile)
        self.output_expressions(outfile)

        outfile.write("\n</instance>\n")
        outfile.close()

    def load_variables(self):
        """
        Load in all the variables of the problem and create their
        domains and variable definitions in XML format
        """
        for var in self.__model.variables:
            domain = var.get_domain()
            self.njvar_mapping[var] = self.create_variable(domain[0],
                                                domain[-1], domain,
                                                var.ident, var.name())

    def load_expressions(self):
        """
        Loads all the expressions from the problem into XML format.
        """
        for con in self.__model.get_exprs():
            if self.global_map.get(con.get_operator()) != None:
                self.handle_global_constraint(con)
            elif con.get_operator() == "Table":
                self.handle_relations(con)
            else:
                pred = self.get_pred(con)
                if pred != None:
                    self.__predicates.append(pred)
                self.__constraints.append(self.get_con(con))

    def get_pred(self, expr):
        """
        Handles the given predicate. Loads it in and stores it in XML format
        """
        self.__extract_param_counter = 0
        functional_expr = self.get_expr(expr)

        if self.pred_map.get(functional_expr) is None:
            # We haven"t seen this predicate before:
            params = self.get_params(expr)
            count = 0
            for p in set(params):
                count += 1

            pred = "\n\t<predicate name=\"P%d\">\n" % self.__local_pred_idx
            pred += "\t\t<parameters>\n\t\t\t"
            pred += "int X"+" int X".join([str(i) for i in range(count)])
            pred += "\n\t\t</parameters>\n"

            pred += "\t\t<expression>\n"
            pred += "\t\t\t<functional>\n"
            pred += "\t\t\t\t%s\n" % functional_expr
            pred += "\t\t\t</functional>\n"
            pred += "\t\t</expression>\n"
            pred += "\t</predicate>\n"

            self.pred_map[functional_expr] = (pred, self.__local_pred_idx)
            self.__local_pred_idx += 1
            return pred

    def get_con(self, expr):
        """
        This creates the constraint nescessary for the given expression.
        It adds it to the list of constraints.
        """
        params = []
        self.extract_params_for_constraint(expr, params)
        parameters, dvars = [], []
        for par in params:
            if par[0] == "VAR":
                if par[1] not in dvars:
                    dvars.append("V%d" % par[1])
                    parameters.append("V%d" % par[1])
            else:
                parameters.append("%d" % par[1])

        arity = len(dvars)
        if arity > self.__maxConstraintArity:
            self.__maxConstraintArity = arity

        pred_string = self.get_expr(expr)
        pred_ident  = self.pred_map.get(pred_string)[1]

        con  = "\n\t<constraint "
        con += "name=\"C%d\" arity=\"%d\" scope=\"%s\" reference=\"P%d\">\n" % (
                    self.__local_con_idx,
                    arity,
                    " ".join([str(i) for i in dvars]),
                    pred_ident )
        con += "\t\t<parameters>\n"
        con += "\t\t\t" + " ".join(parameters) + "\n"
        con += "\t\t</parameters>\n"
        con += "\t</constraint>\n"

        self.__local_con_idx += 1
        return con

    def get_params(self, expr):
        """
        Returns the list of normalised parameters for the given expression
        """
        params = []
        self.__extract_param_counter = 0
        self.extract_params(expr, params)
        return params

    def get_expr(self, expr):
        """
        Returns the expression in string form, normalised parameters
        """
        self.__dVars = {}
        self.__extract_param_counter = 0
        return self.extract_expr(expr)

    def extract_params(self, expr, params):
        """
        Extracts the normalised parameters from the given expression
        """
        if type(expr) == int:
            params.append(("INT", self.__extract_param_counter))
            self.__extract_param_counter += 1
        elif type(expr) == bytes:
            if expr in self.__strings:
                params.append(("INT", self.__strings.index(expr)))
            else:
                self.__strings.append(expr)
                params.append(("INT", self.__strings.index(expr)))
            self.__extract_param_counter += 1
        elif expr.get_operator() == "Sum":
            params.append(("VAR", expr.ident))
            self.__extract_param_counter += 1
        elif expr.is_var():
            params.append(("VAR", expr.ident))
            self.__extract_param_counter += 1
        elif len(expr.get_children()) >= 1:
            for ex in expr.get_children():
                self.extract_params(ex, params)

    def extract_expr(self, expr):
        """
        Extracts a string representation of the given expression
        with normalised parameteres.
        """
        extracted_expr = None
        if type(expr) == int:
            self.__extract_param_counter += 1
            return "X%d" % (self.__extract_param_counter -1)
        elif type(expr) == bytes:
            self.__extract_param_counter += 1
            return "X%d" % (self.__extract_param_counter -1)
        elif expr.get_operator() == "Sum":
            if(len(expr.parameters[0]) == 2 and
               expr.parameters[0][0] in [-1,1] and
               expr.parameters[0][1] in [-1,1]):
                if expr.parameters[0][1] == 1:
                    expr.operator = "add"
                else:
                    expr.operator = "sub"
                extracted_expr = self.extract_expr(expr)
            else:
                self.__extract_param_counter += 1
                extracted_expr = "X%d" % (self.__extract_param_counter - 1)
        elif self.global_map.get(expr.get_operator()) != None:
            self.__extract_param_counter += 1
            extracted_expr = "X%d" % (self.__extract_param_counter -1)
        elif expr.is_var():
            if self.__dVars.get(expr.ident) == None:
                self.__dVars[expr.ident] = self.__extract_param_counter
                self.__extract_param_counter += 1
                extracted_expr = "X%d" % (self.__extract_param_counter -1)
            else:
                extracted_expr = "X%d" % (self.__dVars.get(expr.ident))
        # If this expression should be taking ints but is being given bools:
        elif(expr.get_operator() in self.int2bool and
             expr.get_children()[0].get_operator() in self.int2bool and
             expr.get_children()[1].get_operator() in self.int2bool):
            extracted_expr = self.int2bool[expr.get_operator()].replace("#a",
                   self.extract_expr(expr.get_children()[0])).replace("#b",
                   self.extract_expr(expr.get_children()[1]))
        elif len(expr.get_children()) >= 1:
            extracted_expr = "%s(" % (expr.get_operator())
            subexpr = []
            for ex in expr.get_children():
                subexpr.append(self.extract_expr(ex))
            extracted_expr += ",".join(subexpr) +")"
        return extracted_expr

    def extract_params_for_constraint(self, expr, params):
        """
        Extracts the parameters in terms of their variable identifiers and
        wheather they are constants.
        """
        if type(expr) == int:
            params.append(("INT", expr))
        elif type(expr) == bytes:
            pass
        elif expr.get_operator() == "Sum":
            params.append(("VAR", self.handle_weighted_sum(expr)))
        elif self.global_map.get(expr.get_operator()) != None:
            params.append(("VAR", self.handle_global_constraint(expr)))
        elif expr.is_var():
            params.append(("VAR", expr.ident))
        elif len(expr.get_children()) >= 1:
            for ex in expr.get_children():
                self.extract_params_for_constraint(ex, params)

    def handle_global_constraint(self, expr):
        """
        This adds a global constraint to the list of constraints
        """
        if self.global_map.get(expr.get_operator()) is "NOT_SUPPORTED":
            print("Global operator, %s, is not supported, yet" % expr.get_operator())
        else:
            conformat = self.global_format.get(expr.get_operator())
            params = []
            for exp in expr.get_children():
                self.extract_params_for_constraint(exp, params)

            # Create intermediate variable for storing result of Element
            index = None
            if expr.get_operator() is "Element":
                index = self.create_variable(1, len(expr.get_children())-1,
                                             Name="Result of Element C%d" %
                                                self.__local_con_idx)
                print("C%d created aux element var" % (self.__local_con_idx), index)
                params.append(("VAR", index))
            dparams = sorted(set(params))
            scope = [str(par[1]) for par in dparams if par[0] == "VAR"]

            arity = len(scope)
            if arity > self.__maxConstraintArity:
                self.__maxConstraintArity = arity

            con_string = "\n\t<constraint name=\"C%d\" " % self.__local_con_idx
            con_string += "arity=\"%d\" scope=\"V%s\" reference=\"global:%s\">\n" % (
                arity,
                " V".join(scope),
                self.global_map.get(expr.get_operator())
            )
            self.__local_con_idx += 1

            con_string += "\t\t<parameters>\n"
            con_string += "\t\t\t"

            if "ORTHOTOPES" in conformat:
                con_string += "\t\t\t[\n"
                for row in expr.get_children()[0]:
                    con_string += "\t\t\t\t[ "
                    for orthotope in row:
                        con_string += "{ V%d V%d V%d }" % orthotope
                    con_string += " ]\n"
                con_string += "\t\t\t]\n"

            if "S1" in conformat:
                val = expr.get_children()[0]
                if isinstance(val, int):
                    con_string += "\t\t\t%d\n" % val
                else:
                    con_string += "\t\t\tV%d\n" % val.ident

            if "VARS" in conformat:
                if expr.get_operator() is "Element":
                    con_string += "V%d\n\t\t\t" % index
                    del(params[-1])
                con_string += "["
                for par in params:
                    if par[0] == "VAR":
                        con_string += " V%d" % par[1]
                    else:
                        con_string += " %d" % par[1]
                con_string += " ]\n"

            if "VALS" in conformat:
                weights = expr.parameters[1]
                vals = expr.get_children()
                con_string += "\t\t\t[ "
                for i in range(len(weights)):
                    con_string += "{ %d V%d }" % (weights[i], vals[i].ident)
                con_string += " ]\n"

            if "TASKS" in conformat:
                con_string += "\t\t\t[ "
                for task in expr.parameters[0]:
                    origin   = task.get_lb()
                    duration = task.duration
                    end      = task.get_ub()
                    height   = task#.something
                    con_string += "{ O%d D%d E%d H%d}" (origin,duration,end,height)
                con_string += " ]\n"

            if "MATRIX" in conformat:
                con_string += "\t\t\t[ "
                for row in expr.parameters[0]: # (int,int,int)
                    con_string += "{ %d %d %d }" % row
                con_string += " ]\n"

            if "RELOP" in conformat:
                con_string += "\t\t\t<%s/>\n" % expr.parameters[0][0]

            if "E1" in conformat:
                val = expr.get_children()[-1]
                if type(val) == int :
                    con_string += "\t\t\t%d\n" % val
                else:
                    con_string += "\t\t\tVV%d\n" % val.ident

            con_string += "\t\t</parameters>\n\t</constraint>\n"
            self.__constraints.append(con_string)

            return index

    def handle_weighted_sum(self, expr):
        """
        This handles a weighted sum. It posts the sum constraint and links
        back to a predicate tree with an intermediate variable that is created.
        """
        # Get the children
        children = [self.handle_extra_expression(child)
                    for child in expr.get_children()]
        params = []
        for i in range(len(children)):
            params.append((expr.parameters[0][i], children[i][0]))

        dchildren = []
        for child in children:
            if child[0] not in dchildren:
                dchildren.append(child[0])

        arity     = len(dchildren)
        scope     = "V" + " V".join([str(i) for i in sorted(dchildren)])

        # Create the intermediate variable V{idx} to store result
        # Sum(...+ ((-1) * (V{idx})) )
        weights = expr.parameters[0]
        offset = expr.parameters[1]
        lb, ub = offset, offset
        for i in range(len(children)):
            if weights[i] < 0:
                lb += weights[i]*children[i][2]
                ub += weights[i]*children[i][1]
            else:
                lb += weights[i]*children[i][1]
                ub += weights[i]*children[i][2]

        idx = self.create_variable(lb, ub, Name="Result of WSum:C%d" % self.__local_con_idx)
        params.append((-1, idx))
        scope += " V%d" % idx
        arity += 1

        if arity > self.__maxConstraintArity:
            self.__maxConstraintArity = arity

        con_string = "\n\t<constraint name=\"C%d\" " % self.__local_con_idx
        con_string += "arity=\"%d\" scope=\"%s\" " % ( arity, scope )
        con_string += "reference=\"global:weightedSum\">\n"
        self.__local_con_idx += 1

        con_string += "\t\t<parameters>\n\t\t\t"
        con_string += "[ "
        for coef, var in params:
            con_string += "{%d V%d} " % (coef, var)
        con_string += "]\n"
        con_string += "\t\t\t<eq/>\n"
        con_string += "\t\t\t%d\n" % -(offset)
        con_string += "\t\t</parameters>\n"
        con_string += "\t</constraint>\n"

        # Post the constraint
        self.__constraints.append(con_string)

        return idx

    def handle_extra_expression(self, expr):
        """
        Handles an expression that needs to be represented by an
        intermediate variable. Mainly for stuff under a weighted sum.
        """
        if type(expr) == int:
            return expr
        elif expr.is_var():
            return (expr.ident, expr.get_domain_tuple()[0], expr.get_domain_tuple()[1])
        elif expr.get_operator() == "Sum":
            self.handle_weighted_sum(expr)
        elif self.global_map.get(expr.get_operator()) != None:
            self.handle_global_constraint(expr)
        else:
            #if expr.get_children()[0].is_var():
            #    var = expr.get_children()[0]
            #else:
            #    var = self.handle_extra_expression(expr)

            # Create a variable
            v_lb, v_ub, domain = expr.get_children()[0].get_domain_tuple()
            coef = expr.get_children()[1]
            lb = min(v_lb * coef, v_ub * coef)
            ub = max(v_lb * coef, v_ub * coef)
            var_idx = self.create_variable(lb, ub)

            # Create a predicate
            function = self.get_expr(expr)
            params   = self.get_params(expr)

            dparams = []
            for par in params:
                if par[0] == "VAR":
                    if par[1] not in dparams:
                        dparams.append(par[1])

            arity = len(dparams)
            scope = "V" + " V".join([str(i) for i in sorted(dparams)])

            if arity > self.__maxConstraintArity:
                self.__maxConstraintArity = 2

            if self.pred_map.get(function) is None:
                pred_idx = self.__local_pred_idx
                self.__local_pred_idx += 1
                pred_string = "\n\t<predicate name=\"P%d\">\n" % pred_idx
                pred_string += "\n\t\t<parameters>\n"
                pred_string += "\t\t\t\t%s\n" % params
                pred_string += "\t\t</parameters>\n"

                pred_string += "\n\t\t<expression>\n"
                pred_string += "\t\t\t<functional>\n"
                pred_string += "\t\t\t\t%s\n" % function
                pred_string += "\t\t\t</functional>\n"
                pred_string += "\t\t</expression>\n"
                pred_string += "\t</predicate>\n"

                # Create a constraint
                self.__predicates.append(pred_string)
                self.pred_map[function] = (pred_string, pred_idx)

            else:
                pred_idx = self.pred_map.get(function)[1]

            con_string = "\n\t<constraint name=\"C%d\"" % self.__local_con_idx
            con_string += " arity=\"%d\" scope=\"%s\" reference=\"P%d\">\n" % (
                arity,
                scope,
                pred_idx
            )
            self.__local_con_idx += 1

            con_string += "\t\t<parameters>\n"
            con_string += "\t\t\t" + scope + " %d" % coef + "\n"
            con_string += "\t\t</parameters>\n"
            con_string += "\t</constraint>\n"

            self.__constraints.append(con_string)
            return (var_idx, lb, ub)

    def handle_relations(self, expr):
        """
        Handles conflict and support relations
        adds relation and constraint to XCSP document
        """
        tuples = expr.parameters[0]
        semantics = expr.parameters[1]
        scope = expr.get_children()

        if semantics is "conflict":
            semantics = "conflicts"
        elif semantics is "support":
            semantics = "supports"

        relation  = "\t<relation name=\"R%d\"" % self.__local_rel_idx
        relation += " arity=\"%d\"" % len(scope)
        relation += " nbTuples=\"%d\"" % len(tuples)
        relation += " semantics=\"%s\">\n\t\t" % semantics
        relation += "|".join([" ".join([str(t) for t in tuple]) for tuple in tuples])
        relation += "\n\t</relation>\n"
        self.__relations.append(relation)

        constraint  = "\n\t<constraint name=\"C%d\" arity=\"%d\" scope=\"%s\" reference=\"R%d\" />\n" % (
            self.__local_con_idx,
            len(scope),
            " ".join(["V%d" % var.ident for var in scope]),
            self.__local_rel_idx )

        self.__constraints.append(constraint)

        if self.__maxConstraintArity < len(scope):
            self.__maxConstraintArity = len(scope)

        self.__local_rel_idx += 1
        self.__local_con_idx += 1

    def output_domains(self, outfile):
        """
        Outputs all the domains to the file
        """
        outfile.write("\n<domains nbDomains=\"%d\">" % (
            len(self.__domains)))
        for dom in self.__domains:
            outfile.write(dom)
        outfile.write("</domains>\n")

    def output_variables(self, outfile):
        """
        Outputs all the created variables to the file
        """
        outfile.write("\n<variables nbVariables=\"%d\">\n" % (
            len(self.__variables)))
        for var in self.__variables:
            outfile.write(var)
        outfile.write("</variables>\n")

    def output_relations(self, outfile):
        """
        Outputs all the relations to the file
        """
        if len(self.__relations)>0:
            outfile.write("\n<relations nbRelations=\"%d\">\n" % (
                    len(self.__relations)))
            for rel in self.__relations:
                outfile.write(rel)
            outfile.write("</relations>\n")

    def output_expressions(self, outfile):
        """
        Print all the constraints and predicates to the file
        """
        if len(self.__predicates)>0:
            outfile.write("\n<predicates nbPredicates=\"%d\">" % len(self.__predicates))
            for pred in self.__predicates:
                outfile.write(pred)
            outfile.write("</predicates>\n")

        outfile.write("\n<constraints nbConstraints=\"%d\">" % len(self.__constraints))
        for con in self.__constraints:
            outfile.write(con)
        outfile.write("</constraints>\n")

    def create_variable(self, lb, ub, domain=None, index=None, Name=None):
        """
        Creates a variable in the variable list for the XML. If the domain
        has already been specified it uses the previously specified domain.
        Otherwise it also creates a domain, automatically reducing sequencial
        numbers in the domain to intervals
        """

        if domain is None:
            domain = [lb,ub]

        dom_string = ""

        if index is None:
            idx = self.__local_var_idx
            self.__local_var_idx += 1
        else:
            idx = index
            self.__local_var_idx = index + 1

        d_idx = len(self.__domains)

        nbvalues = 0
        dom_list = []
        if len(domain) == 2:
            nbvalues = len(list(range(lb, ub+1)))
            if lb != ub:
                dom_list.append("%d..%d" % (lb, ub))
            else:
                dom_list.append(str(ub))
        else:
            clean_domain = []
            for i in set(domain):
                if i not in clean_domain:
                    clean_domain.append(i)
            domain = sorted(clean_domain)
            nbvalues = len(domain)
            for interval in change_to_intervals(domain):
                if interval[0] != interval[1]:
                    dom_list.append("%d..%d" % (interval[0], interval[1]))
                else:
                    dom_list.append("%d" % interval[0])

        dom_string += "\n\t<domain name=\"D%d\" nbValues=\"%d\">\n\t\t" % (
                d_idx, nbvalues)

        dom_list_string = " ".join(dom_list)
        dom_string += dom_list_string

        dom_string += "\n\t</domain>\n"

        if self.dom_map.get(dom_list_string) is None:
            self.dom_map[dom_list_string] = d_idx
            self.__domains.append(dom_string)
        else:
            d_idx = self.dom_map.get(dom_list_string)

        var_string = "\t<variable name=\"V%d\" domain=\"D%d\"/>" % (idx, d_idx)
        if Name is not None and Name != "x":
            var_string += "\t<!-- %s -->" % Name
        var_string += "\n"

        self.__variables.append(var_string)
        return idx

def change_to_intervals(domain):
    """
    Represents the given domain as a list of intervals
    """
    intervals = []
    interval_start = domain[0]
    for i in range(0, len(domain) - 1):
        if domain[i] + 1 != domain[i+1]:
            interval_end = domain[i]
            intervals.append( (interval_start, interval_end) )
            interval_start = domain[i+1]
    interval_end = domain[-1]
    intervals.append( (interval_start, interval_end) )
    return intervals
