
# translate flatzinc format to python/Numberjack 

BEGIN {
	MAXINT = 10000000;
	CPUTIME = 900;

	print "#! /usr/bin/env python"
	print "from Numberjack import *"
	system("cat fzn2py.py");
	print "";
	print "model = Model()";
	parameter = 1;
	error = 0;
}

# { print $0 } # for debugging purposes only

/::_output_/ {
	gsub("::_output_","");
	$0 = $0 " ::_output_";
}

/^var / {
	parameter = 0;
}

/ of var / {
	parameter = 0;
}

parameter {
	sub(".*: ","");
	if ($2 == "=") print $0; # it might be a predicate instead if no =
}

/^var bool:/{
	print $3 " = Variable('" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "model.add(" $3 " == " $5 ")";
}

/^var int:/{
	print $3 " = Variable(-" MAXINT "," MAXINT ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "model.add(" $3 " == " $5 ")";
}

/^var range[(][-]*[0-9]+,1[+][-]*[0-9]+[)]:/{
	sub("range[(]","",$2);
	sub("[)]:","",$2);
	sub(",1[+]",",",$2);
	print $3 " = Variable(" $2 ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "model.add(" $3 " == " $5 ")";
}

/^var {[-]*[0-9]+(,[-]*[0-9]+)*}:/{
	sub("{","[",$2);
	sub("}:","]",$2);
	print $3 " = Variable(" $2 ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "model.add(" $3 " == " $5 ")";
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var bool:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	name = $6;
	print name " = VarArray(" isup ",'" name "')";
	if (match($0,"::_output_")) {
		for (i=0; i<isup; i++) {
			output[ name "[" i "]"] = 1;
		}
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "model.add(" name "[" i "] == " $0 "[" i "])";
		}
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var int:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	name = $6;
	print name " = VarArray(" isup ",-" MAXINT "," MAXINT ",'" name "')";
	if (match($0,"::_output_")) {
		for (i=0; i<isup; i++) {
			output[ name "[" i "]"] = 1;
		}
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "model.add(" name "[" i "] == " $0 "[" i "])";
		}
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var range[(][-]*[0-9]+,1[+][-]*[0-9]+[)]:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	sub("range[(]","",$5);
	sub("[)]:","",$5);
	sub(",1[+]",",",$5);
	name = $6;
	print name " = VarArray(" isup "," $5 ",'" name "')";
	if (match($0,"::_output_")) {
		for (i=0; i<isup; i++) {
			output[ name "[" i "]"] = 1;
		}
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "model.add(" name "[" i "] == " $0 "[" i "])";
		}
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var {[-]*[0-9]+(,[-]*[0-9]+)*}:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	sub("{","[",$5);
	sub("}:","]",$5);
	name = $6;
	print name " = [Variable(" $5 ",'" name "_" i "') for i in range(" isup ")]";
	if (match($0,"::_output_")) {
		for (i=0; i<isup; i++) {
			output[ name "[" i "]"] = 1;
		}
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "model.add(" name "[" i "] == " $0 "[" i "])";
		}
	}
}

/^constraint /{
	$1 = "model.add(";
    $NF = $NF " )";
	print $0;
}

/^solve / && / minimize /{
	print "model.add(Minimize(" $NF "))";
	objective = $NF;
}

/^solve / && / maximize /{
	print "model.add(Maximize(" $NF "))";
	objective = $NF;
}

END {
	if (!error) {

	print "solvers = ['Mistral', 'SCIP', 'MiniSat', 'toulbar2', 'Gurobi']";
	print "default = dict([('solver', 'Mistral'), ('verbose', 1), ('tcutoff', 900), ('var', 'DomainOverWDegree'), ('val', 'Lex'), ('rand', 2)])";
    print "param = input(default)";
    print "solver = model.load(param['solver'])";
    print "solver.setVerbosity(param['verbose'])";
    print "solver.setTimeLimit(param['tcutoff'])";
    print "solver.setHeuristic(param['var'], param['val'], param['rand'])";
    print "if param['solver'] == 'Mistral':";
    print "    solver.solveAndRestart(GEOMETRIC, 256, 1.3)";
    print "else:";
    print "    solver.solve()";

	print "sat_result = 's SATISFIABLE' if solver.is_sat() else 's UNKNOWN'"
	print "if solver.is_unsat():"
	print "    sat_result = 's UNSATISFIABLE'"
	print "print sat_result"
	if (objective){
		print "if solver.is_sat():"
		print "    print 'c Optimal', 1 if solver.is_opt() else 0";
	}
	print "print 'c SolveTime', solver.getTime()";
	print "print 'c Nodes', solver.getNodes()";
	if (objective) print "print 'c Objective'," objective ".get_value()";
	for (e in output) {
		name = e;
		sub("[]]","+1]",name);
		print "print 'constraint int_eq(" name ",'," e ".get_value(),');'";
	}
	}
}
