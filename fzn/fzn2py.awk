
# translate flatzinc format to python/Numberjack 

BEGIN {
	MAXINT = 10000000;
	print "#! /usr/bin/env python"
	print "import time"
	print "from Numberjack import *"
	system("cat " MZNNJ_DIR "/fzn2py.py");
	print "";
	print "model = Model()";	
	parameter = 1;
	error = 0;
}

# { print $0 } # for debugging purposes only

/::_output_var/ {
	gsub("::_output_var","");
	$0 = $0 " ::_output_";
}

/::_output_array/ {
	regexp = "::_output_array[(][[](range[(][0-9]+,1[+][0-9]+[)][ ,]*)+[]][)]";
	if (match($0,regexp)) {
	  array = substr($0, RSTART, RLENGTH);
	  sub(regexp,"");
	  gsub("[[]","",array);
	  gsub("[]]","",array);
	  gsub(" ","",array);
	  gsub("range[(]","",array);
	  gsub(",1[+]","..",array);
	  gsub("[)]","",array);
	  nbdim = array;
	  gsub("[^,]","",nbdim);
	  dim = length(nbdim)+1;
	  sub("::_output_array", "array" dim "d",array);
	  $0 = $0 " ::_output_" array;
	} else {
		print "ERROR WRONG OUTPUT ARRAY DEFINITION",$0;
		error = 3;
		exit(2);
	}
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
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
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
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
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
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
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
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
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
	print "solvers = ['Mistral', 'SCIP', 'MiniSat', 'Toulbar2', 'Gurobi']";
	print "default = dict([('solver', 'Mistral'), ('verbose', 0), ('tcutoff', 900), ('var', 'DomainOverWDegree'), ('val', 'Lex'), ('rand', 2), ('threads', 1)])";
	print "param = input(default)";
	print "solver = model.load(param['solver'])";
	print "solver.setVerbosity(param['verbose'])";
	print "solver.setTimeLimit(param['tcutoff'] - int(time.clock()+0.5))";
	print "solver.setHeuristic(param['var'], param['val'], param['rand'])";
	print "if param['solver'] == 'Gurobi':";
	print "    solver.setThreadCount(param['threads'])";
	print "if param['solver'] == 'Mistral':";
	print "    solver.solveAndRestart(GEOMETRIC, 256, 1.3)";
	print "else:";
	print "    solver.solve()";
	print "if solver.is_sat():"
	n = asorti(output,varnames);
	for (i=1; i<=n; i++) {
		e = varnames[i];
		if (output[e]==1) {
			if (e != objective && e != "objective" && e != "obj") print "    print '" e " = '," e ".get_value(),';'";
			else print "    print '" e " = ', (solver.getOptimum() if param['solver'] == 'Toulbar2' else " e ".get_value()),';'";
		} else {
			print "    print '" e " = " outputstring[e] ",',str(" e "),');'";
		}
	}
	print "    print '----------'";
	print "if solver.is_unsat():"
	print "    print '=====UNSATISFIABLE====='"
	if (objective){
		print "elif solver.is_opt():"
		print "    print '=========='"
	} else {
		print "elif solver.is_sat():"
		print "    print '=========='"
	}
        print "else:"
	print "    print '=====UNKNOWN====='"
	print "print '% SolveTime', solver.getTime()";
	print "print '% Nodes', solver.getNodes()";
	if (objective) print "if solver.is_sat(): print '% Objective', (solver.getOptimum() if param['solver'] == 'Toulbar2' else " objective ".get_value())";
	}
}
