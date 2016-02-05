
# translate flatzinc format to python/Numberjack 

BEGIN {
	MAXINT = 10000000;
	objective_type = "";
	system("cat " MZNNJ_DIR "/fzn2py.py");
	print "";
	print "def get_model():"
	print "    model = Model()";	
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
	if ($2 == "=") print "    " $0; # it might be a predicate instead if no =
}

/^var bool:/{
	print "    " $3 " = Variable('" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    model.add(" $3 " == " $5 ")";
}

/^var int:/{
	print "    " $3 " = Variable(-" MAXINT "," MAXINT ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    model.add(" $3 " == " $5 ")";
}

/^var range[(][-]*[0-9]+,1[+][-]*[0-9]+[)]:/{
	sub("range[(]","",$2);
	sub("[)]:","",$2);
	sub(",1[+]",",",$2);
	print "    " $3 " = Variable(" $2 ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    model.add(" $3 " == " $5 ")";
}

/^var {[-]*[0-9]+(,[-]*[0-9]+)*}:/{
	sub("{","[",$2);
	sub("}:","]",$2);
	print "    " $3 " = Variable(" $2 ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    model.add(" $3 " == " $5 ")";
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
	print "    " name " = VarArray(" isup ",'" name "')";
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "    model.add(" name "[" i "] == " $0 "[" i "])";
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
	print "#    " name " = VarArray(" isup ",-" MAXINT "," MAXINT ",'" name "')";
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "#    model.add(" name "[" i "] == " $0 "[" i "])";
		}
        print "    " name " = VarArray(" $0 ")";
	} else  if (match($7,"::var_is_introduced")) {
		sub(".*= ","");
		sub("::_output_.*","");
		sub("::var_is_introduced","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "#    model.add(" name "[" i "] == " $0 "[" i "])";
		}
        print "    " name " = VarArray(" $0 ")";
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
	print "    " name " = VarArray(" isup "," $5 ",'" name "')";
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "    model.add(" name "[" i "] == " $0 "[" i "])";
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
	print "    " name " = [Variable(" $5 ",'" name "_" i "') for i in range(" isup ")]";
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "    model.add(" name "[" i "] == " $0 "[" i "])";
		}
	}
}

/^constraint /{
	$1 = "    model.add(";
    	$NF = $NF " )";
	print $0;
}

/^solve / && / minimize /{
	print "    model.add(Minimize(" $NF "))";
	objective = $NF;
	objective_type = "Minimize";
}

/^solve / && / maximize /{
	print "    model.add(Maximize(" $NF "))";
	objective = $NF;
	objective_type = "Maximize";
}

function trimvarindex(s) {
	ind = index(s, "[");
	if(ind > 0) {
		return substr(s, 1, ind - 1);
	}
	return s;
}

END {
	if (!error) {
	
	# Create an array vars_array, keyed with the output variable names so that
	# we output a set of unique variable names. We will also trim indexed vars.
	if(objective) vars_array[trimvarindex(objective)] = 1;

	n = asorti(output,varnames);
	for (i=1; i<=n; i++) {
		e = varnames[i];
		if(e != objective) {
			vars_array[trimvarindex(e)] = 1;
		}
	}
	output_vars = sep = "";
	for (x in vars_array) {
	    output_vars = output_vars sep x;
	    sep = ", ";
	}
	print "    output_vars = (" output_vars ")";
	print "    return model, output_vars";
	print "\n";

	if(objective){
		print "def solve_dichotomic(param):";
		print "    model, output_vars = get_model()";
		if(output_vars) print "    " output_vars " = output_vars";
		print "    lb = reallb = " objective ".lb";
		print "    ub = realub = " objective ".ub";
		print "    best_sol = (None, output_vars)";
		print "    dichotomic_sat = dichotomic_opt = False";
		if(objective_type == "Minimize"){
			print "    while lb < ub - 1 and time_remaining(param['tcutoff']) > param['dichtcutoff']:";
		} else if(objective_type == "Maximize"){
			print "    while lb + 1 < ub and time_remaining(param['tcutoff']) > param['dichtcutoff']:";
		}
		print "        newobj = (lb + ub) / 2";
		print "        # print lb, ub, newobj";
		print "        dummymodel, output_vars = get_model()";
		if(output_vars) print "        " output_vars " = output_vars";
		if(objective_type == "Minimize"){
			print "        dummymodel.add(" objective " > reallb)";
			print "        dummymodel.add(" objective " <= newobj)";
		} else if(objective_type == "Maximize"){
			print "        dummymodel.add(" objective " <= realub)";
			print "        dummymodel.add(" objective " > newobj)";
		}
		print "        dichparam = dict(param)";
		print "        dichparam['tcutoff'] = param['dichtcutoff']";
		print "        solver, output_vars = run_solve(dummymodel, output_vars, dichparam)";
		print "";
		if(objective_type == "Minimize"){
			print "        if solver.is_opt():";
			print "            lb = reallb = " objective ".get_value() - 1";
			print "        if solver.is_sat():";
			print "            ub = " objective ".get_value()";
			print "            best_sol = solver, output_vars";
			print "            dichotomic_sat = True";
			print "        elif solver.is_unsat():";
			print "            lb = reallb = newobj";
			print "        else:";
			print "            lb = newobj";
			print "    if reallb < ub - 1:";
		} else if(objective_type == "Maximize"){
			print "        if solver.is_opt():";
			print "            ub = realub = " objective ".get_value() + 1";
			print "        if solver.is_sat():";
			print "            lb = " objective ".get_value()";
			print "            best_sol = solver, output_vars";
			print "            dichotomic_sat = True";
			print "        elif solver.is_unsat():";
			print "            ub = realub = newobj";
			print "        else:";
			print "            ub = newobj";print "";
			print "    if realub > lb + 1:";
		}
		print "        dummymodel, output_vars = get_model()";
		if(output_vars) print "        " output_vars " = output_vars";
		if(objective_type == "Minimize"){
			print "        dummymodel.add(" objective " > reallb)";
			print "        dummymodel.add(" objective " <= ub)";
		} else if(objective_type == "Maximize"){
			print "        dummymodel.add(" objective " <= realub)";
			print "        dummymodel.add(" objective " > lb)";
		}
		print "        tcutoff = time_remaining(param['tcutoff'])";
		print "        if tcutoff > 1.0:";
		print "            dichparam = dict(param)";
		print "            dichparam['tcutoff'] = tcutoff";
		print "            solver, output_vars = run_solve(dummymodel, output_vars, dichparam)";
		print "            if solver.is_sat():";
	    print "                best_sol = solver, output_vars";
		print "    else:";
		print "        dichotomic_opt = True";
		print "";
		print "    if not solver.is_sat() and dichotomic_sat:";
		print "        best_sol[0].is_sat = lambda: True";
		print "        best_sol[0].is_unsat = lambda: False";
		print "        if dichotomic_opt:";
		print "            best_sol[0].is_opt = lambda: True";
		print "    return best_sol";
	}

	print "\n";
	print "start_time = datetime.datetime.now()\n\n";
	print "if __name__ == '__main__':";
	print "    solvers = ['Mistral', 'SCIP', 'MiniSat', 'Toulbar2', 'Gurobi']";
	print "    default = dict([('solver', 'Mistral'), ('verbose', 0), ('tcutoff', 900), ('var', 'DomainOverWDegree'), ('val', 'Lex'), ('rand', 2), ('threads', 1), ('restart', GEOMETRIC), ('base', 256), ('factor', 1.3), ('lcLevel', 4), ('lds', 0), ('dee',0), ('btd',0), ('rds',0), ('dichotomic', 0), ('dichtcutoff', 5), ('encoding', '')])";
	print "    param = input(default)";
	if(objective){
		print "    if param['dichotomic'] == 1:";
	    print "        solver, output_vars = solve_dichotomic(param)";
	    print "    else:";
	    print "        solver, output_vars = solve_main(param)";
	} else {
		print "    solver, output_vars = solve_main(param)";
	}
	if(output_vars) print "    " output_vars " = output_vars";
	print "";
	print "    if not solver:";
    print "        print '=====UNKNOWN====='";
    print "        sys.exit(0)";
    print "";
	print "    if solver.is_sat():"
	for (i=1; i<=n; i++) {
		e = varnames[i];
		if (output[e]==1) {
			if (e != objective && e != "objective" && e != "obj") print "        print '" e " = '," e ".get_value(),';'";
			else print "        print '" e " = ', (solver.getOptimum() if param['solver'] == 'Toulbar2' else " e ".get_value()),';'";
		} else {
			print "        print '" e " = " outputstring[e] ",',str(" e "),');'";
		}
	}
	print "        print '----------'";
	if (objective){
		print "        if solver.is_opt():"
		print "            print '=========='"
	} else {
		print "        print '=========='"
	}
	print "    elif solver.is_unsat():"
	print "        print '=====UNSATISFIABLE====='"
    print "    else:"
	print "        print '=====UNKNOWN====='"
	print "    print '% SolveTime', solver.getTime()";
	print "    print '% Nodes', solver.getNodes()";
	print "    print '% Failures', solver.getFailures()";
	if (objective) print "    if solver.is_sat(): print '% Objective', (solver.getOptimum() if param['solver'] == 'Toulbar2' else " objective ".get_value())";
	}
}
