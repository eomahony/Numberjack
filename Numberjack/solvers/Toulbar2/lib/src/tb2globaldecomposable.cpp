#include "tb2globaldecomposable.hpp"

/// DECOMPOSABLE COST FUNCTION /////////////////////////////////////////

DecomposableGlobalCostFunction::DecomposableGlobalCostFunction() : arity(0), scope(NULL), label("empty") {
    ToulBar2::Berge_Dec=1;
}

DecomposableGlobalCostFunction::DecomposableGlobalCostFunction(unsigned int _arity, int* _scope) : arity(_arity), label("empty"){
    scope = new int[arity];
    for (unsigned int variable = 0 ; variable < _arity ; ++variable) {
        scope[variable] = _scope[variable];
    }
    ToulBar2::Berge_Dec=1;
}

DecomposableGlobalCostFunction::~DecomposableGlobalCostFunction() {
    delete [] scope;
}

DecomposableGlobalCostFunction* 
DecomposableGlobalCostFunction::FactoryDGCF(string type, unsigned int _arity, int* _scope, istream &file) {
    //cout << "Creating a " << type << " global cost function " << endl;
    if (type == "wamong")               return new WeightedAmong(_arity,_scope,file);
    if (type == "wvaramong")            return new WeightedVarAmong(_arity,_scope,file);
    if (type == "wsum")                 return new WeightedSum(_arity,_scope,file);
    if (type == "wvarsum")              return new WeightedVarSum(_arity,_scope,file);
    if (type == "woverlap")             return new WeightedOverlap(_arity,_scope,file);

    if (type == "walldifferent" || type == "walldiff")      return new WeightedAllDifferent(_arity,_scope,file);
    if (type == "wgcc")                 return new WeightedGcc(_arity,_scope,file);
    if (type == "wregular")             return new WeightedRegular(_arity,_scope,file);
    if (type == "wsame")                return new WeightedSame(_arity,_scope,file);
    if (type == "wsamegcc")             return new WeightedSameGcc(_arity,_scope,file);

    cout << type << " unknown decomposable global cost function" << endl;
    return 0;
}

void 
DecomposableGlobalCostFunction::color(int i) {
    switch (i) {
    case 1 : cout << "\033[41m"; break;
    case 2 : cout << "\033[42m"; break;
    case 3 : cout << "\033[43m"; break;
    case 4 : cout << "\033[44m"; break;
    case 5 : cout << "\033[45m"; break;
    case 6 : cout << "\033[46m"; break;
    case 7 : cout << "\033[47m"; break;
    case 8 : cout << "\033[40m\033[37m"; break;
    default : cout << "\033[0m"; break;
    };
}

/// WEIGHTED AMONG /////////////////////////////////////////////////////

WeightedAmong::WeightedAmong() : DecomposableGlobalCostFunction() {}

WeightedAmong::WeightedAmong(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {
}

WeightedAmong::WeightedAmong(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    unsigned int nbValue;
    file >> nbValue;
    for (unsigned int value = 0 ; value < nbValue ; ++value) {
        int valueRead;
        file >> valueRead;
        values.insert(valueRead);
    }
    file >> lb >> ub;
}

WeightedAmong::~WeightedAmong() {
    values.clear();
}

void
WeightedAmong::addToCostFunctionNetwork(WCSP* wcsp) {
    bool VERBOSE = false;
    bool VVERBOSE = false;
    int nbVariableCFN = wcsp->numberOfVariables();
    //cout << nbVariableCFN << endl;

    // -- new variables : counters -- //
    int addVariablesIndex[arity+1];
    for (int newVariable = 0 ; newVariable <= arity ; newVariable++) {
        string varname="WAmong" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
        addVariablesIndex[newVariable] = wcsp->makeEnumeratedVariable(varname,0,newVariable);
        if (VERBOSE) {color(5); cout << "new variable " << addVariablesIndex[newVariable] << "("<< ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getDomainInitSize()<< ")"; color(0); cout << endl;}
    }

    Cost top = wcsp->getUb();
    // -- ternary constraints : partial sum -- //
    for (int variable = 0 ; variable < arity ; ++variable) {
        int indexCi = addVariablesIndex[variable];
        int indexCj = addVariablesIndex[variable+1];
        int indexXi  = scope[variable];
        if (VERBOSE) {color(5); cout << indexCi << "--" << indexXi << "--" << indexCj; color(0); cout << endl;}
        EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi);
        EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj);
        EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
        wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
        wcsp->getListSuccessors()->at(indexXi).push_back(indexCj);



        unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
        vector<Cost> ternaryCosts(tableSize,top);

        for (unsigned long valuePosition = 0; valuePosition < varXi->getDomainInitSize() ; valuePosition++) {
            Value value =  varXi->toValue(valuePosition);
            for (unsigned long counterPosition = 0; counterPosition < varCi->getDomainInitSize() ; counterPosition++) {
                long counter =   varCi->toValue(counterPosition);

                int nextCounter = counter;
                if (values.find(value) != values.end()) {
                    nextCounter++;
                }
                if (VVERBOSE) 		cout << counter << "(" << counterPosition << ")" << " && " << value << "(" << valuePosition << ")" << " ==> " << nextCounter << "\t";
                //	unsigned long position =  (counter) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize())
                //							+ (value) 			* (varCj->getDomainInitSize())
                //							+ (nextCounter);
                unsigned long position =  (counterPosition) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize())
										        +  (valuePosition) 			* (varCj->getDomainInitSize())
										        +  (nextCounter);
                if (VVERBOSE)  cout << position << "/" << tableSize << endl;
                ternaryCosts[position] = 0;


            }
        }
        wcsp->postTernaryConstraint(indexCi, indexXi, indexCj, ternaryCosts);
    }

    // -- unary constraints : final variable -- //
    if (VERBOSE) {color(5); cout << "post unary constraint on " << addVariablesIndex[arity]; color(0); cout << endl;}
    vector<Cost> unaryCosts(arity+1,top);
    for (int i = 0 ; i <= arity ; i++) {
        int gap = max( 0 , max( int(lb - i), int(i - ub) ) );
        if (semantics == "hard") {
            if (((unsigned int) i) >= lb && ((unsigned int) i)  <= ub) unaryCosts[i] = 0;
            else unaryCosts[i] = min(top,baseCost);
        }
        if (semantics == "lin" || semantics == "var")  unaryCosts[i] = min (top, baseCost * gap);
        if (semantics == "quad") unaryCosts[i] = min (top,baseCost * gap * gap);
        if (VVERBOSE) cout << i << " => " << unaryCosts[i] << endl;
    }
    wcsp->postUnary(addVariablesIndex[arity],unaryCosts);
}

Cost 
WeightedAmong::evaluate(int* tuple) {
    int occurency = 0;
    for (int var = 0 ; var < arity ; var++) {
        if(values.find(tuple[var]) != values.end()) occurency++;
    }
    int gap = max( 0 , max( int(lb - occurency), int(occurency - ub) ) );
    if (gap) {
        if (semantics == "hard") return baseCost;
        if (semantics == "lin" || semantics == "var")  return baseCost * gap;
        if (semantics == "quad") return baseCost * gap * gap;
    }
    return 0;
}

void 
WeightedAmong::display() {
    cout << "WAmong (" << arity << ") : ";
    for (int variable = 0 ; variable < arity ; ++variable) {
        cout << scope[variable] << " ";
    }
    cout << endl;
    cout << "sem : " << semantics << " " << baseCost << endl;
    cout << "val : ";
    for (set<int>::iterator value = values.begin() ; value != values.end() ; value++) {
        cout << *value << " ";
    }
    cout << endl;
    cout << "bounds [" << lb << ":" << ub << "]" << endl;
}

/// WEIGHTED REGULAR ///////////////////////////////////////////////////

WeightedRegular::WeightedRegular() : DecomposableGlobalCostFunction(), automaton(0) {}

WeightedRegular::WeightedRegular(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope),automaton(0) {}

WeightedRegular::WeightedRegular(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    automaton = new WFA(file);
}

WeightedRegular::~WeightedRegular() {
    delete automaton;
}

void
WeightedRegular::addToCostFunctionNetwork(WCSP* wcsp) {
    ToulBar2::Berge_Dec=1;
    //automaton->display();
    Cost top = wcsp->getUb();
    int unsigned current_var_number = wcsp->numberOfVariables();
    int unsigned q0 = current_var_number;
    if ( ToulBar2::verbose > 1 ) 	{
        cout << "DEBUG>> wregular found : initial number of variables before creation = " <<  wcsp->numberOfVariables() <<endl;
        cout << "DEBUG>> wregular Automatum Total number of states: " <<  automaton->getNbStates() <<endl;
        cout << "DEBUG>> wregular Initial states number: " << automaton->getInitialStates().size()  <<endl;
        cout << "DEBUG>> wregular add new variable from " << q0 <<" to "<<  current_var_number+arity+1 <<endl;
    }
    if(  current_var_number > 0 ) {
        int unsigned domsize  = automaton->getNbStates()-1;
        string varname = "WR" + to_string(current_var_number);
        if ( ToulBar2::verbose > 1 ) cout << "DEBUG>> wregular q0 index "<< q0 << " domain = " << domsize+1 << endl;
        wcsp->makeEnumeratedVariable(varname,(Value) 0, (Value) domsize);			// add q0 variable
        if ( ToulBar2::verbose > 1 ) cout << "wregular add varname =" << varname <<"=> var index "<<  wcsp->numberOfVariables() << " domain size = " << domsize+1 << endl;
    } else { exit(EXIT_FAILURE);
    }
    //################################################initial state ##############################################
    if (automaton->getInitialStates().size() > 0 ) {
        vector<Cost> initial_states_costs(automaton->getNbStates(),top);
        list< pair<int,Cost> > initialStates = automaton->getInitialStates();
        for (list<pair<int,Cost> >::iterator it = initialStates.begin(); it !=  initialStates.end() ; ++it){
            pair<int,Cost> initial = *it;
            //cout << initial.first << " " << initial.second << endl;
            initial_states_costs[initial.first]=initial.second;
        }
        wcsp->postUnary(q0,initial_states_costs);
        if ( ToulBar2::verbose > 1 ) {
            cout << "DEBUG>> wregular initial state (q0) vector size ( nbre value) = "<<  initial_states_costs.size() << endl;
            cout << "DEBUG>> wregular var q0= "<< q0 <<" number of constraint wregular initialisation ==>" << wcsp->numberOfConstraints()  << endl;
        }
    }
    //################################################accepting state ##############################################
    for( int v = 1 ; v < arity+1 ; v++ )  {
        int unsigned domsize = automaton->getNbStates()-1;
        string varname = to_string(v+q0);

        DEBONLY(int theindex =) wcsp->makeEnumeratedVariable(varname,(Value) 0,(Value) domsize);			// add qi variable
        assert(theindex == v+ (int) current_var_number);
        if ( ToulBar2::verbose > 1 ) cout << "DEBUG>> wregular add varname =" << varname <<"=> rank "<<  wcsp->numberOfVariables() << " domain = " << domsize+1 << endl;
    }
    int unsigned q_last = wcsp->numberOfVariables() -1 ;
    if ( ToulBar2::verbose > 1 ) cout << "DEBUG>> wregular Final number of variables : " << wcsp->numberOfVariables() << endl;
    vector<Cost>final_states_costs(automaton->getNbStates(),top);

    list< pair<int,Cost> > acceptingStates = automaton->getAcceptingStates();
    if (acceptingStates.size()>=0) {

        for (list<pair<int,Cost> >::iterator it = acceptingStates.begin(); it !=  acceptingStates.end() ; ++it ) {
            pair<int,Cost> accept = *it;
            int unsigned t_index = accept.first;
            Cost ucost = accept.second;

            EnumeratedVariable* Qv  = (EnumeratedVariable *) wcsp->getVar(q_last); // get domaine size of last qi var
            unsigned long DomVar=Qv->getDomainInitSize();

            if(t_index < DomVar){
                final_states_costs[t_index]=ucost;
            } else {
                cout <<"wregular tuple error " << t_index << "out of domain" << DomVar << endl;
                exit(EXIT_FAILURE);
            }
        }
        wcsp->postUnary(q_last,final_states_costs);

        if ( ToulBar2::verbose > 1 )  {
            cout << "DEBUG>> wregular last q varname = "<<  q_last <<endl;
        }
    }
    /*
				//################################################### lecture des transition ????
				if ( ToulBar2::verbose > 1 ) 
				cout << "DEBUG>>wregular final number of Unary constraint Post after q0 and qi post : " << numberOfConstraints()  << endl;
				//==================
				// transition stat reading 
				//==================
				int nb_transition;
				vector<unsigned int> VQi;
				vector<unsigned int> VQj;
				vector<unsigned int> VXi;
				vector<Cost> transition_costs;
				file >> nb_transition;
				if ( ToulBar2::verbose > 1 ) cout << "DEBUG>> wregular transitions number :  " <<  nb_transition <<endl;
				for( int s = 0 ; s < nb_transition ; s++)
				{
					int qi;
					int xi;
					int qj;
					Cost transition_COST;
					file >> qi;
					file >> xi;
					file >> qj;
					file >> transition_COST;
					VQi.push_back(qi);
					VXi.push_back(xi);
					VQj.push_back(qj);
					transition_costs.push_back(transition_COST);

					if ( ToulBar2::verbose > 1 ) {
					cout << "DEBUG>> wregular read transition table qi =" << qi << " xi =" << xi << " qj =" << qj << " cost=" << transition_COST << endl;
					cout << "DEBUG>> wregular scope xi " << scopeIndex[xi] <<endl;
					}
				}
     */
    //##################################################ajout ternaire#############################
    for ( int q = 0 ; q < arity ; q++) {
        int qi = q0+q;
        int	xi = scope[q] ;
        int qj = qi+1;
        if ( ToulBar2::verbose > 1 ) cout << "DEBUG>>wregular  post ternary on  qi =" << qi << " xi =" << xi << " qj =" << qj << endl;
        // poiner on qi , xi , qj varibale
        EnumeratedVariable* Qi = (EnumeratedVariable *) wcsp->getVar(qi); //current qi variable;
        EnumeratedVariable* Xi = (EnumeratedVariable *) wcsp->getVar(xi); //current Xi variable;
        EnumeratedVariable* Qj = (EnumeratedVariable *) wcsp->getVar(qj); //current qj variable;
        // domain definition
        unsigned long DomQi=Qi->getDomainInitSize();
        unsigned long DomQj=Qi->getDomainInitSize();
        unsigned long DomXi=Xi->getDomainInitSize();
        unsigned long Domsize = long (Qi->getDomainInitSize() * Xi->getDomainInitSize() * Qj->getDomainInitSize());

        vector<Cost>tmp_ternary_costs(Domsize,top);
        list<WTransition*> transitions = automaton->getTransitions();
        for(list<WTransition*>::iterator it = transitions.begin() ; it != transitions.end() ; ++it)  {
            WTransition* transition = *it;
            int start = transition->start;
            int end = transition->end;
            unsigned int symbol = transition->symbol;
            int positionSymbol = 0;
            for (EnumeratedVariable::iterator iter = Xi->begin(); iter != Xi->end(); ++iter) {
                if (symbol == Xi->toIndex(*iter)) break;
                positionSymbol++;
            }
            int weight = transition->weight;

            if(symbol < DomXi) {
                unsigned long cindex= start*DomXi*DomQj + positionSymbol*DomQj + end;
                tmp_ternary_costs[cindex]=weight;
            }
            if(ToulBar2::verbose > 1) {
                cout << "DEBUG>> wregular init cost vector for ternary rank= " << q << endl;
                cout << "DEBUG>> wregular add ternary table qi =" << qi << " xi =" << xi << " qj =" << qj << endl;
                cout <<"DEBUG>> wregular Ternary const DOMAIN SIZE = " << Domsize <<" Dom qi ="<< DomQi << " Dom Xi=" << DomXi << " Dom Qj =" << DomQj <<" -------"<<endl;
                cout << "DEBUG>> wregular initial COST SIZE = " << tmp_ternary_costs.size() << endl;
                cout << "DEBUG>> wregular number of constraint before ternary cost function : " << wcsp->numberOfConstraints()  << endl;
            }
        }
        wcsp->postTernaryConstraint(qi, xi, qj, tmp_ternary_costs);
        wcsp->getListSuccessors()->at(qi).push_back(xi);
        wcsp->getListSuccessors()->at(xi).push_back(qj);
    }
    if ( ToulBar2::verbose >=1 ) cout << "DEBUG>> wregular Total number of constraints after wregular posting " << wcsp->numberOfConstraints()  << endl;
}

void 
WeightedRegular::display() {
    cout << "WRegular (" << arity << ") : ";
    for (int variable = 0 ; variable < arity ; ++variable) {
        cout << scope[variable] << " ";
    }
    cout << endl;
    if (automaton) {
        automaton->display();
    }
    else {
        cout << "no automaton associated" << endl;
    }
}

/// WEIGHTED SUM ///////////////////////////////////////////////////////

WeightedSum::WeightedSum() : DecomposableGlobalCostFunction() {}

WeightedSum::WeightedSum(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}

WeightedSum::WeightedSum(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    file >> comparator >> rightRes;
}

WeightedSum::~WeightedSum() {}

void 
WeightedSum::addToCostFunctionNetwork(WCSP* wcsp) {
    int nbVariableCFN = wcsp->numberOfVariables();
    //cout << nbVariableCFN << endl;

    // -- new variables : counters -- //
    int addVariablesIndex[arity+1];
    int cumulDOWN = 0;
    int cumulUP = 0;
    for (int newVariable = 0 ; newVariable <= arity ; newVariable++) {
        string varname="WSum" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
        addVariablesIndex[newVariable] = wcsp->makeEnumeratedVariable(varname,cumulDOWN,cumulUP);
        //cout << "\033[45m" << "new variable \033[0m" << addVariablesIndex[newVariable] << "("<< ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getInf()<< ":" << ":" << ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getSup() << ")" << "\033[0m" << endl;
        if (newVariable < arity) {
            cumulDOWN += ((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getInf();
            cumulUP   += ((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getSup();
            //cout << newVariable << " ["<<((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getInf()<<":"<<((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getSup()<<"] " << endl;
            //cout << cumulDOWN << "<->" << cumulUP << endl;
        }
    }

    Cost top = wcsp->getUb();
    // -- ternary constraints : partial sum -- //
    for (int variable = 0 ; variable < arity ; ++variable) {
        int indexCi = addVariablesIndex[variable];
        int indexCj = addVariablesIndex[variable+1];
        int indexXi = scope[variable];
        EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi);
        EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj);
        EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
        wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
        wcsp->getListSuccessors()->at(indexXi).push_back(indexCj);

        //cout << "\033[46m" << "post ternary constraint on " << indexCi <<"("<<varCi->getDomainInitSize()<<")" << ","  << indexXi <<"("<<varXi->getDomainInitSize()<<")" << ","  << indexCj <<"("<<varCj->getDomainInitSize()<<")" << "\033[0m" << endl;

        unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
        vector<Cost> ternaryCosts(tableSize,top);

        for (unsigned long valueXPosition = 0; valueXPosition < varXi->getDomainInitSize() ; valueXPosition++) {
            long value =  varXi->toValue(valueXPosition);
            for (unsigned long valueCiPosition = 0; valueCiPosition < varCi->getDomainInitSize() ; valueCiPosition++) {
                long counter_i =  varCi->toValue(valueCiPosition);
                for (unsigned long valueCjPosition = 0; valueCjPosition < varCj->getDomainInitSize() ; valueCjPosition++) {
                    long counter_j =  varCj->toValue(valueCjPosition);

                    if (counter_j == (counter_i + value)) {
                        //cout << counter_i << " + " << value << " = " << counter_j << endl;
                        unsigned long position =  (valueCiPosition) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize())
												         + (valueXPosition) 		* (varCj->getDomainInitSize())
												         + (valueCjPosition);
                        ternaryCosts[position] = 0;
                    }
                }
            }
        }
        wcsp->postTernaryConstraint(indexCi, indexXi, indexCj, ternaryCosts);
    }

    // -- unary constraints : final variable -- //
    //cout << "\033[46m" << "post unary constraint on " << addVariablesIndex[arity] << "\033[0m" << endl;
    int sizeLast = cumulUP - cumulDOWN + 1;
    //cout << cumulDOWN << " -> " << cumulUP << " :: " << sizeLast << endl;
    vector<Cost> unaryCosts(sizeLast ,top);
    int positionVar = 0;

    for (int i = cumulDOWN ; i <= cumulUP ; i++) {
        bool compFound = false;
        int gap = 0;

        if (comparator == "==") {
            compFound = true;
            if (i != rightRes) gap = (i < rightRes)?(rightRes - i):(i - rightRes);
        }
        if (comparator == "!=") {
            compFound = true;
            if (i == rightRes) gap=1;
        }
        if (comparator == "<" || comparator == "<=") {
            compFound = true;
            int newRightRes = rightRes;
            if (comparator == "<") newRightRes--;
            if (i > newRightRes)  gap = i - newRightRes;
        }
        if (comparator == ">" || comparator == ">=") {
            compFound = true;
            int newRightRes = rightRes;
            if (comparator == ">") newRightRes++;
            if (i < newRightRes)  gap = newRightRes - i;
        }

        if (semantics == "hard")  unaryCosts[positionVar] = (gap)?baseCost:0;
        if (semantics == "lin" || semantics == "var")   unaryCosts[positionVar] = gap*baseCost;
        if (semantics == "quad")  unaryCosts[positionVar] = gap*gap*baseCost;
        //cout << positionVar << " (=) " << i << " ==> " << unaryCosts[positionVar] << endl;

        if (!compFound) {
            cout << "comparator " << comparator << " not handle yet" << endl;
        }
        positionVar++;
    }
    //EnumeratedVariable* lastVar = (EnumeratedVariable *) wcsp->getVar(addVariablesIndex[arity]);
    //cout << addVariablesIndex[arity] << " " << lastVar->getDomainInitSize() << " [" <<  lastVar->getInf() << ":" << lastVar->getSup()<< "]" << endl;
    wcsp->postUnary(addVariablesIndex[arity],unaryCosts);
    //cout << "after  adding to CFN" << endl;

}

Cost 
WeightedSum::evaluate(int* tuple) {
    int sum = 0;
    for (int var = 0 ; var < arity ; var++) {
        sum += tuple[var];
    }
    if (comparator == "==") {
        int gap = (sum < rightRes)?  rightRes - sum: sum - rightRes;
        if (semantics == "hard") return min(gap*baseCost,baseCost);
        if (semantics == "lin" || semantics == "var")  return gap*baseCost;
        if (semantics == "quad") return gap*gap*baseCost;
    }
    if (comparator == "!=") {
        if (sum != rightRes) return baseCost;
    }
    if (comparator == "<" || comparator == "<=") {
        int newRightRes = rightRes; if (comparator == "<") newRightRes--;
        int gap = max(0,sum - newRightRes);
        if (semantics == "hard") return min(gap*baseCost,baseCost);
        if (semantics == "lin" || semantics == "var")  return gap*baseCost;
        if (semantics == "quad") return gap*gap*baseCost;
    }
    if (comparator == ">" || comparator == ">=") {
        int newRightRes = rightRes; if (comparator == ">") newRightRes++;
        int gap = max(0,newRightRes - sum);
        if (semantics == "hard") return min(gap*baseCost,baseCost);
        if (semantics == "lin" || semantics == "var")  return gap*baseCost;
        if (semantics == "quad") return gap*gap*baseCost;
    }
    return 0;
}

void 
WeightedSum::display() {
    cout << "WSum (" << arity << ") : ";
    for (int variable = 0 ; variable < arity ; ++variable) {
        cout << scope[variable] << " ";
    }
    cout << endl;
    cout << comparator << " " << rightRes << endl;
    cout << semantics << " " << baseCost << endl;
}

/// WEIGHTED VAR SUM ///////////////////////////////////////////////////////

WeightedVarSum::WeightedVarSum() : DecomposableGlobalCostFunction() {}

WeightedVarSum::WeightedVarSum(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}

WeightedVarSum::WeightedVarSum(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    file >> comparator;
}

WeightedVarSum::~WeightedVarSum() {}

void 
WeightedVarSum::addToCostFunctionNetwork(WCSP* wcsp) {
    int nbVariableCFN = wcsp->numberOfVariables();

    /// -- new variables : counters -- ///
    int addVariablesIndex[arity+1];
    int cumulDOWN = 0;
    int cumulUP = 0;
    for (int newVariable = 0 ; newVariable <= arity-1 ; newVariable++) {
        string varname="WVarSum" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
        addVariablesIndex[newVariable] = wcsp->makeEnumeratedVariable(varname,cumulDOWN,cumulUP);
        //cout << "\033[45m" << "new variable \033[0m" << addVariablesIndex[newVariable] << "("<< ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getInf()<< ":" << ":" << ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getSup() << ")" << "\033[0m" << endl;
        if (newVariable < arity) {
            cumulDOWN += ((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getInf();
            cumulUP   += ((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getSup();
            //cout << newVariable << " ["<<((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getInf()<<":"<<((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getSup()<<"] " << endl;
            //cout << cumulDOWN << "<->" << cumulUP << endl;
        }
    }

    Cost top = wcsp->getUb();
    /// -- ternary constraints : partial sum -- ///
    for (int variable = 0 ; variable < arity-1 ; ++variable) {
        int indexCi = addVariablesIndex[variable];
        int indexCj = addVariablesIndex[variable+1];
        int indexXi = scope[variable];
        EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi);
        EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj);
        EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
        wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
        wcsp->getListSuccessors()->at(indexXi).push_back(indexCj);

        //cout << "\033[46m" << "post ternary constraint on " <<
        //		indexCi <<"("<<varCi->getDomainInitSize()<<" [" << varCi->getInf() <<  ":" <<varCi->getSup()<< "])" << ","  <<
        //		indexXi <<"("<<varXi->getDomainInitSize()<<" [" << varXi->getInf() <<  ":" <<varXi->getSup()<< "])" << ","  <<
        //		indexCj <<"("<<varCj->getDomainInitSize()<<" [" << varCj->getInf() <<  ":" <<varCj->getSup()<< "])" << ","  <<
        //		"\033[0m" << endl;


        unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
        vector<Cost> ternaryCosts(tableSize,top);
        //		cout << "TableSize = " << tableSize << endl;

        EnumeratedVariable::iterator iterXi = varXi->begin();
        for (unsigned long valueXPosition = 0; valueXPosition < varXi->getDomainInitSize() ; valueXPosition++) {
            long value =  *iterXi;
            EnumeratedVariable::iterator iterCi = varCi->begin();
            for (unsigned long valueCiPosition = 0; valueCiPosition < varCi->getDomainInitSize() ; valueCiPosition++) {
                long counter_i =  *iterCi;
                EnumeratedVariable::iterator iterCj = varCj->begin();
                for (unsigned long valueCjPosition = 0; valueCjPosition < varCj->getDomainInitSize() ; valueCjPosition++) {
                    long counter_j =  *iterCj;
                    if (counter_j == (counter_i + value)) {

                        unsigned long position =  (valueCiPosition) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize())
												         + (valueXPosition) 		* (varCj->getDomainInitSize())
												         + (valueCjPosition);
                        //cout << counter_i << " + " << value << " = " << counter_j <<  " ==> "  << position << endl;
                        ternaryCosts[position] = 0;
                    }
                    ++iterCj;
                }
                ++iterCi;
            }
            ++iterXi;
        }
        //cout << "here" << endl;
        //for (int i = 0 ; i < tableSize ; i++) { cout << i << " => tableCost[] = " <<   ternaryCosts[i] << endl; }
        wcsp->postTernaryConstraint(indexCi, indexXi, indexCj, ternaryCosts);
        //cout << "here bis" << endl;
    }

    /// -- binary constraints : final variable -- ///
    //cout << "\033[46m" << "post binary constraint on " << addVariablesIndex[arity-1]  << " and " << scope[arity-1] << "\033[0m" << endl;
    EnumeratedVariable* lastCounter = (EnumeratedVariable *) wcsp->getVar(addVariablesIndex[arity-1]);
    EnumeratedVariable* lastVariable = (EnumeratedVariable *) wcsp->getVar(scope[arity-1]);
    unsigned int lastCounterSize = lastCounter->getDomainInitSize();
    unsigned int lastVariableSize = lastVariable->getDomainInitSize();
    unsigned int costSize = lastCounterSize * lastVariableSize;
    //cout << lastCounterSize << "*" << lastVariableSize << " => " << costSize << endl;
    vector<Cost> binaryCosts(costSize ,top);

    for (unsigned int positionLastCounter = 0 ; positionLastCounter < lastCounterSize ; positionLastCounter++) {
        for (unsigned int positionLastVariable = 0; positionLastVariable < lastVariableSize ; positionLastVariable++) {
            int valueCounter  = lastCounter->toValue(positionLastCounter);
            int valueVariable = lastVariable->toValue(positionLastVariable);
            unsigned int positionArray = positionLastCounter * lastVariableSize + positionLastVariable;
            //cout << valueCounter << " " << comparator << " " << valueVariable << " :: " ;
            int gap = 0;
            if (comparator == "==") {
                gap = (valueCounter < valueVariable)?(valueVariable-valueCounter):(valueCounter-valueVariable);
            }
            if (comparator == "!=") {
                gap = (valueCounter==valueVariable)?1:0;
            }
            if (comparator == "<" || comparator == "<=" ) {
                int modif = (comparator == "<")?1:0;
                gap = (valueCounter>valueVariable-modif)?(valueCounter-valueVariable+modif):0;
            }
            if (comparator == ">" || comparator == ">=" ) {
                int modif = (comparator == ">")?1:0;
                gap = (valueVariable>valueCounter-modif)?(valueVariable-valueCounter+modif):0;
            }

            if (semantics == "hard")  binaryCosts[positionArray] = (gap)?baseCost:0;
            if (semantics == "lin" || semantics == "var")   binaryCosts[positionArray] = (gap*baseCost >= top)?top:gap*baseCost;
            if (semantics == "quad")  binaryCosts[positionArray] = (gap*gap*baseCost >= top)?top:gap*gap*baseCost;
            //cout << valueCounter << "," << valueVariable << "," << gap << " : " << binaryCosts[positionArray] << endl;
        }
    }
    //cout << "here" << endl;
    wcsp->postBinaryConstraint(addVariablesIndex[arity-1],scope[arity-1],binaryCosts);
    //cout << "here bis" << endl;
}


void 
WeightedVarSum::display() {
    cout << "WVarSum (" << arity << ") : ";
    for (int variable = 0 ; variable < arity-1 ; ++variable) {
        cout << scope[variable] << " ";
    }
    cout << comparator << " " << scope[arity-1];
    cout << " (" << semantics << " " << baseCost << ")" << endl;
}

/// WEIGHTED OVERLAP ///////////////////////////////////////////////////

WeightedOverlap::WeightedOverlap() : DecomposableGlobalCostFunction() {}

WeightedOverlap::WeightedOverlap(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}

WeightedOverlap::WeightedOverlap(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    file >> comparator >> rightRes;
    //display();
}

WeightedOverlap::~WeightedOverlap() {}

void 
WeightedOverlap::addToCostFunctionNetwork(WCSP* wcsp) {
    int nbVariableCFN = wcsp->numberOfVariables();

    // -- new variables : counters OVERLAP -- //
    int addVariablesOverlap[arity/2];
    for (int newVariable = 0 ; newVariable < arity/2 ; newVariable++) {
        string varname="WOVERL_OVER_" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
        addVariablesOverlap[newVariable] = wcsp->makeEnumeratedVariable(varname,0,1);
        //cout << "add overlap " << newVariable << endl;
    }

    Cost top = wcsp->getUb();
    // -- ternary -- //
    for (int newVariable = 0 ; newVariable < arity/2 ; newVariable++) {
        int X = scope[newVariable];
        int Y = scope[newVariable+arity/2];
        int O = addVariablesOverlap[newVariable];
        EnumeratedVariable* varX = (EnumeratedVariable *) wcsp->getVar(X);
        EnumeratedVariable* varY = (EnumeratedVariable *) wcsp->getVar(Y);
        EnumeratedVariable* varO = (EnumeratedVariable *) wcsp->getVar(O);
        wcsp->getListSuccessors()->at(X).push_back(O);
        wcsp->getListSuccessors()->at(Y).push_back(O);

        unsigned long tableSize = long (varX->getDomainInitSize() * varY->getDomainInitSize() * varO->getDomainInitSize());
        vector<Cost> ternaryCosts(tableSize,top);

        //cout << X << "--" << Y  << "--" << O << endl;
        for (unsigned long valueXPosition = 0; valueXPosition < varX->getDomainInitSize() ; valueXPosition++) {
            long vX =  varX->toValue(valueXPosition);
            //for (unsigned int vY = 0 ; vY < varY->getDomainInitSize() ; vY++) {
            for (unsigned long valueYPosition = 0; valueYPosition < varY->getDomainInitSize() ; valueYPosition++) {
                long vY =  varY->toValue(valueYPosition);

                unsigned int vO = 0;
                if (vX == vY && vY != 0)  vO = 1;
                //cout << "X=" << vX << " Y=" << vY << " O=" << vO << endl;
                unsigned long position =  (valueXPosition) 		* (varY->getDomainInitSize()*varO->getDomainInitSize())
										        + (valueYPosition) 		* (varO->getDomainInitSize())
										        + (vO);
                ternaryCosts[position] = 0;
            }
        }
        wcsp->postTernaryConstraint(X,Y,O,ternaryCosts);
    }

    // -- new variables : counters Among -- //
    int addVariablesAmong[arity/2+1];
    for (int newVariable = 0 ; newVariable < arity/2+1 ; newVariable++) {
        string varname="WOVERL_AMONG_" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
        addVariablesAmong[newVariable] = wcsp->makeEnumeratedVariable(varname,0,newVariable);
        //cout << "add among " << newVariable << endl;
    }

    // -- ternary -- //
    for (int newVariable = 0 ; newVariable < arity/2 ; newVariable++) {
        int indexCi = addVariablesAmong[newVariable];
        int indexCj = addVariablesAmong[newVariable+1];
        int indexXi  = addVariablesOverlap[newVariable];

        //cout << indexCi << " -- " << indexXi << " -- "  << indexCj << endl;
        EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi);
        EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj);
        EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
        wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
        wcsp->getListSuccessors()->at(indexXi).push_back(indexCj);

        unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
        vector<Cost> ternaryCosts(tableSize,top);

        for (unsigned long value = 0; value < varXi->getDomainInitSize() ; value++) {
            for (unsigned long counter = 0; counter < varCi->getDomainInitSize() ; counter++) {
                int nextCounter = counter;
                if (value == 1) {
                    nextCounter++;
                }
                unsigned long position =  (counter) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize())
											        + (value) 			* (varCj->getDomainInitSize())
											        + (nextCounter);
                ternaryCosts[position] = 0;
            }
        }

        wcsp->postTernaryConstraint(indexCi,indexXi,indexCj,ternaryCosts);
    }

    // -- unary constraints : final variable -- //
    //cout << "\033[45m" << "post unary constraint on " << addVariablesIndex[arity] << "\033[0m" << endl;
    vector<Cost> unaryCosts(arity/2+1 ,top);
    if (comparator == "==") {
        for (int i = 0 ; i <= arity/2 ; i++) {
            if (i == rightRes) unaryCosts[i] = 0;
            else {
                int gap = (i < rightRes)?  rightRes - i: i - rightRes;
                if (semantics == "hard") unaryCosts[i] = baseCost;
                if (semantics == "lin" || semantics == "var")  unaryCosts[i] = gap*baseCost;
                if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
            }
        }
    } else
        if (comparator == "!=") {
            for (int i = 0 ; i <= arity/2 ; i++) {
                if (i != rightRes) unaryCosts[i] = 0;
                else {
                    unaryCosts[i] = baseCost;
                }
            }
        } else
            if (comparator == "<" || comparator == "<=") {
                int newRightRes = rightRes;
                if (comparator == "<") newRightRes--;

                for (int i = 0 ; i <= arity/2 ; i++) {
                    if (i <= newRightRes) unaryCosts[i] = 0;
                    else {
                        int gap = max(0,i - newRightRes);
                        if (semantics == "hard") unaryCosts[i] = baseCost;
                        if (semantics == "lin" || semantics == "var")  unaryCosts[i] = gap*baseCost;
                        if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
                    }
                }
            } else
                if (comparator == ">" || comparator == ">=") {
                    int newRightRes = rightRes;
                    if (comparator == ">") newRightRes++;

                    for (int i = 0 ; i <= arity/2 ; i++) {
                        if (i >= newRightRes) unaryCosts[i] = 0;
                        else {
                            int gap = max(0,newRightRes - i);
                            if (semantics == "hard") unaryCosts[i] = baseCost;
                            if (semantics == "lin" || semantics == "var")  unaryCosts[i] = gap*baseCost;
                            if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
                        }
                    }
                }
                else cout << "comparator " << comparator << " not handle yet" << endl;
    //cout << "control " << addVariablesAmong[arity/2] << endl;
    wcsp->postUnary(addVariablesAmong[arity/2],unaryCosts);
}

Cost 
WeightedOverlap::evaluate(int* tuple) {
    int occurency = 0;
    for (int var = 0 ; var < arity/2 ; var++) {
        if(tuple[var] && tuple[var] == tuple[var+arity/2]) {
            //cout << var << " && " << var+arity/2;
            occurency++;
        }
    }
    //cout << " => " << occurency << " " << comparator << " " << rightRes << endl;
    if (comparator == "==") {
        int gap = (occurency < rightRes)?  rightRes - occurency: occurency - rightRes;
        if (semantics == "hard") return min(gap*baseCost,baseCost);
        if (semantics == "lin" || semantics == "var")  return gap*baseCost;
        if (semantics == "quad") return gap*gap*baseCost;
    }
    if (comparator == "!=") {
        if (occurency != rightRes) return baseCost;
    }
    if (comparator == "<" || comparator == "<=") {
        int newRightRes = rightRes; if (comparator == "<") newRightRes--;
        int gap = max(0,occurency - newRightRes);
        if (semantics == "hard") return min(gap*baseCost,baseCost);
        if (semantics == "lin" || semantics == "var")  return gap*baseCost;
        if (semantics == "quad") return gap*gap*baseCost;
    }
    if (comparator == ">" || comparator == ">=") {
        int newRightRes = rightRes; if (comparator == ">") newRightRes++;
        int gap = max(0,newRightRes - occurency);
        if (semantics == "hard") return min(gap*baseCost,baseCost);
        if (semantics == "lin" || semantics == "var")  return gap*baseCost;
        if (semantics == "quad") return gap*gap*baseCost;
    }
    return 0;
}

void 
WeightedOverlap::display() {
    cout << "WOverlap (" << arity << ") : ";
    for (int variable = 0 ; variable < arity ; ++variable) {
        cout << scope[variable] << " ";
    }
    cout << endl;
    cout << semantics << " " << baseCost << endl;
    int i = 0;
    cout << "{ " ;
    for ( ; i < arity/2 ; i++) cout << scope[i]<< " ";
    cout << "}" << endl;
    cout << "{ " ;
    for ( ; i < arity ; i++) cout << scope[i]<< " ";
    cout << "}" << endl;
    cout << comparator << " " << rightRes << endl;
}

////////////////////////////////////////////////////////////////////////
// EXPERIMENTAL CONSTRAINTS                                           //
////////////////////////////////////////////////////////////////////////

/// WEIGHTED VAMONG /////////////////////////////////////////////////////

WeightedVarAmong::WeightedVarAmong() : DecomposableGlobalCostFunction() {}

WeightedVarAmong::WeightedVarAmong(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {
}

WeightedVarAmong::WeightedVarAmong(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    unsigned int nbValue;
    file >> nbValue;
    for (unsigned int value = 0 ; value < nbValue ; ++value) {
        int valueRead;
        file >> valueRead;
        values.insert(valueRead);
    }
    //file >> lb >> ub;
    index=scope[arity-1];
}

WeightedVarAmong::~WeightedVarAmong() {
    values.clear();
}

//TODO writing the other semantics
void
WeightedVarAmong::addToCostFunctionNetwork(WCSP* wcsp) {
    bool VERBOSE = false;
    int nbVariableCFN = wcsp->numberOfVariables();
    Cost top = wcsp->getUb();

    // -- new variables : counters -- //
    int addVariablesIndex[arity];
    for (int newVariable = 0 ; newVariable < arity ; newVariable++) {
        string varname="WAmong" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
        addVariablesIndex[newVariable] = wcsp->makeEnumeratedVariable(varname,0,newVariable);
        if (VERBOSE) cout << "\033[45m" << "new variable " << addVariablesIndex[newVariable] << "("<< ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getDomainInitSize()<< ")" << "\033[0m" << endl;
    }

    // -- ternary constraints : partial sum -- //
    for (int variable = 0 ; variable < arity-1 ; ++variable) {
        int indexCi = addVariablesIndex[variable];
        int indexCj = addVariablesIndex[variable+1];
        int indexXi  = scope[variable];
        if (VERBOSE) cout << "\033[45m" << indexCi << "--" << indexXi << "--" << indexCj << "\033[0m" << endl;
        EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi);
        EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj);
        EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
        wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
        wcsp->getListSuccessors()->at(indexXi).push_back(indexCj);

        unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
        vector<Cost> ternaryCosts(tableSize,top);

        for (unsigned long valuePosition = 0; valuePosition < varXi->getDomainInitSize() ; valuePosition++) {
            long value =  varXi->toValue(valuePosition);
            for (unsigned long counterPosition = 0; counterPosition < varCi->getDomainInitSize() ; counterPosition++) {
                long counter =   varCi->toValue(counterPosition);

                int nextCounter = counter;
                if (values.find(value) != values.end()) {
                    nextCounter++;
                }
                unsigned long position =  (counterPosition) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize())
											        + (valuePosition) 			* (varCj->getDomainInitSize())
											        + (nextCounter);
                ternaryCosts[position] = 0;
            }
        }
        wcsp->postTernaryConstraint(indexCi, indexXi, indexCj, ternaryCosts);
    }

    // -- binary constraints : final variable -- //
    if (semantics != "hard") { color(1) ; cout << "WARNING :: only hard semantics can be considered"; color(-1) ; cout << endl; }
    if (VERBOSE) cout << "\033[45m" << "post binary constraint on " << addVariablesIndex[arity-1] << " and " << scope[arity - 1] << "\033[0m" << endl;
    int indexCount = addVariablesIndex[arity-1];
    int indexLast  = scope[arity - 1];
    EnumeratedVariable* varCount = (EnumeratedVariable *) wcsp->getVar(indexCount);
    EnumeratedVariable* varLast  = (EnumeratedVariable *) wcsp->getVar(indexLast);
    wcsp->getListSuccessors()->at(indexCount).push_back(indexLast);
    unsigned long tableSize = long (varCount->getDomainInitSize() * varLast->getDomainInitSize());
    vector<Cost> binaryCosts(tableSize,top);

    for (unsigned long valuePosition = 0; valuePosition < varLast->getDomainInitSize() ; valuePosition++) {
        long value = varLast->toValue(valuePosition);

        for (unsigned long counterPosition = 0; counterPosition < varCount->getDomainInitSize() ; counterPosition++) {
            long counter =   varCount->toValue(counterPosition);

            //cout << value << " == " << counter << endl;

            unsigned long position = (counterPosition) * (varLast->getDomainInitSize()) + valuePosition;
            if (counter==value)
                binaryCosts[position] = 0;
            else
                binaryCosts[position] = baseCost;
        }
    }
    wcsp->postBinaryConstraint(indexCount, indexLast, binaryCosts);
}

void 
WeightedVarAmong::display() {
    cout << "WVarAmong (" << arity << ") : ";
    for (int variable = 0 ; variable < arity-1 ; ++variable) {
        cout << scope[variable] << " ";
    }
    cout << " == " << scope[arity - 1 ] << endl;
    cout << "sem : " << semantics << " " << baseCost << endl;
    cout << "val : ";
    for (set<int>::iterator value = values.begin() ; value != values.end() ; value++) {
        cout << *value << " ";
    }
    cout << endl;
}

/// WEIGHTED ALLDIFFERENT //////////////////////////////////////////////

WeightedAllDifferent::WeightedAllDifferent() : DecomposableGlobalCostFunction() {}

WeightedAllDifferent::WeightedAllDifferent(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}

WeightedAllDifferent::WeightedAllDifferent(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    //display();
}

WeightedAllDifferent::~WeightedAllDifferent() {}

void
WeightedAllDifferent::addToCostFunctionNetwork(WCSP* wcsp) {
    // Counting the number of value
    int inf = ((EnumeratedVariable *) wcsp->getVar(scope[0]))->getInf();
    int sup = ((EnumeratedVariable *) wcsp->getVar(scope[0]))->getSup();
    for (int variable = 0 ; variable < arity ; ++variable) {
        int tinf = ((EnumeratedVariable *) wcsp->getVar(scope[variable]))->getInf();
        int tsup = ((EnumeratedVariable *) wcsp->getVar(scope[variable]))->getSup();
        inf = min ( inf , tinf );
        sup = max (sup , tsup );
    }

    // Adding WeightedAmong over each variable
    for (int value = inf ; value <= sup ; value++) {
        WeightedAmong* wamong = new WeightedAmong(arity,scope);
        wamong->setSemantics(semantics);
        wamong->setBaseCost(baseCost);
        wamong->addValue(value);
        wamong->setBounds(0,1);
        wamong->addToCostFunctionNetwork(wcsp);
    }
}

void 
WeightedAllDifferent::display() {
    cout << "WeightedAllDifferent (" << arity << ") : ";
    for (int variable = 0 ; variable < arity ; ++variable) {
        cout << scope[variable] << " ";
    }
    cout << endl;
    cout << semantics << " " << baseCost << endl;
}


/// WEIGHTED GCC ///////////////////////////////////////////////////////

WeightedGcc::WeightedGcc() : DecomposableGlobalCostFunction() {}

WeightedGcc::WeightedGcc(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}

WeightedGcc::WeightedGcc(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    int nbValueToWatch;
    file >> nbValueToWatch;
    for (int idvalue = 0 ; idvalue < nbValueToWatch ; idvalue++) {
        Value value;
        unsigned int lb,ub;
        file >> value >> lb >> ub;
        setBounds(value,lb,ub);
    }
    //display();

}

WeightedGcc::~WeightedGcc() {}

void 
WeightedGcc::setBounds(Value value, unsigned int lb, unsigned int ub) {
    map<Value, pair<unsigned int, unsigned int> >::iterator it;
    it = bounds.find(value);
    if (it != bounds.end()) {
        cerr << "WeightedGcc::setBounds | Value " << value << " is already watch" << endl;
        exit(128);
    }
    bounds[value] = make_pair(lb,ub);
}

void 
WeightedGcc::addToCostFunctionNetwork(WCSP* wcsp) {
    //	int nbcounters = bounds.size();
    //	int counter = 0;
    //	int counters[nbcounters];
    //	int clb[nbcounters];
    //	int cub[nbcounters];
    //	int cscope[nbcounters];
    for (map<Value,pair <unsigned int,unsigned int> >::iterator it = bounds.begin(); it != bounds.end() ; ++it) {
        pair<Value,pair <unsigned int,unsigned int> > bound = *it;

        //Adding a wamong
        Value value = bound.first;
        unsigned int lb = (bound.second).first;
        //		clb[counter] = lb;
        unsigned int ub = (bound.second).second;
        //		cub[counter] = ub;
        WeightedAmong* wamong = new WeightedAmong(arity,scope);
        wamong->setSemantics(semantics);
        wamong->setBaseCost(baseCost);
        wamong->addValue(value);
        wamong->setBounds(lb,ub);
        wamong->addToCostFunctionNetwork(wcsp);
        //		counters[counter] = wcsp->numberOfVariables() - 1;
        //		counter++;
    }
    //	if (semantics == "hard") rec_sum_counters(wcsp, cscope, 0, 0, 0, counters, clb, cub, nbcounters, 0);
}

// Special additional constraint propagation for hard decomposed GCC
// add constraints on end-counters of wamong's decomposed GCC to enforce LB and UB for any subset of values
void WeightedGcc::rec_sum_counters(WCSP *wcsp, int *cscope, int carity, int totlb, int totub, int *counters, int *clb, int *cub, int nb, int rec)
{
    if (rec==nb) {
        if (carity==2) {
            vector<Cost> costs(wcsp->getDomainInitSize(cscope[0]) * wcsp->getDomainInitSize(cscope[1]), wcsp->getUb());
            for (int u=wcsp->getInf(cscope[0]); u<=wcsp->getSup(cscope[0]); u++) {
                for (int v=wcsp->getInf(cscope[1]); v<=wcsp->getSup(cscope[1]); v++) {
                    if (u+v >= totlb && u+v <= min(totub,arity)) {
                        costs[wcsp->toIndex(cscope[0],u) * wcsp->getDomainInitSize(cscope[1]) + wcsp->toIndex(cscope[1],v)] = MIN_COST;
                    }
                }
            }
            wcsp->postBinaryConstraint(cscope[0], cscope[1], costs);
        } else if (carity==3) {
            vector<Cost> costs(wcsp->getDomainInitSize(cscope[0]) * wcsp->getDomainInitSize(cscope[1]) * wcsp->getDomainInitSize(cscope[2]), wcsp->getUb());
            for (int u=wcsp->getInf(cscope[0]); u<=wcsp->getSup(cscope[0]); u++) {
                for (int v=wcsp->getInf(cscope[1]); v<=wcsp->getSup(cscope[1]); v++) {
                    for (int w=wcsp->getInf(cscope[2]); w<=wcsp->getSup(cscope[2]); w++) {
                        if (u+v+w >= totlb && u+v+w <= min(totub,arity)) {
                            costs[wcsp->toIndex(cscope[0],u) * wcsp->getDomainInitSize(cscope[1]) * wcsp->getDomainInitSize(cscope[2]) + wcsp->toIndex(cscope[1],v) * wcsp->getDomainInitSize(cscope[2]) + wcsp->toIndex(cscope[2],w)] = MIN_COST;
                        }
                    }
                }
            }
            wcsp->postTernaryConstraint(cscope[0], cscope[1], cscope[2], costs);
        } else if (carity>3) {
            wcsp->postWSum(cscope,carity,"hard",wcsp->getUb(),">=",totlb);
            wcsp->postWSum(cscope,carity,"hard",wcsp->getUb(),"<=",min(totub,arity));
        }
    } else {
        // try without variable at position rec
        rec_sum_counters(wcsp,cscope,carity,totlb,totub,counters,clb,cub,nb,rec+1);
        // try with variable at position rec
        cscope[carity]=counters[rec];
        rec_sum_counters(wcsp,cscope,carity+1,totlb+clb[rec],totub+cub[rec],counters,clb,cub,nb,rec+1);
    }
}

void 
WeightedGcc::display() {
    cout << "WeightedGcc (" << arity << ") : ";
    for (int variable = 0 ; variable < arity ; ++variable) {
        cout << scope[variable] << " ";
    }
    cout << endl;
    cout << semantics << " " << baseCost << endl;
    for (map<Value,pair <unsigned int,unsigned int> >::iterator it = bounds.begin(); it != bounds.end() ; ++it) {
        pair<Value,pair <unsigned int,unsigned int> > bound = *it;
        cout << bound.first << " [" << (bound.second).first << ":" << (bound.second).second << "]" << endl;
    }
}



/// WEIGHTED SAME //////////////////////////////////////////////////////

WeightedSame::WeightedSame() : DecomposableGlobalCostFunction() {}
WeightedSame::WeightedSame(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}
WeightedSame::WeightedSame(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    if (_arity % 2 == 1 ) {
        cerr << "WeightedSame::Constructor the scope must be even" << endl;
        exit(128);
    }
    //display();
}
WeightedSame::~WeightedSame() {}

void 
WeightedSame::addToCostFunctionNetwork(WCSP* wcsp) {
    Cost top = wcsp->getUb();

    // Counting the number of value
    int inf = ((EnumeratedVariable *) wcsp->getVar(scope[0]))->getInf();
    int sup = ((EnumeratedVariable *) wcsp->getVar(scope[0]))->getSup();
    for (int variable = 0 ; variable < arity ; ++variable) {
        int tinf = ((EnumeratedVariable *) wcsp->getVar(scope[variable]))->getInf();
        int tsup = ((EnumeratedVariable *) wcsp->getVar(scope[variable]))->getSup();
        inf = min ( inf , tinf );
        sup = max (sup , tsup );
    }
    int nbValue = sup - inf + 1;

    // Creating new counter variables
    int positionVar = 0;
    int** newVariable = new int*[nbValue];
    for (int value = inf ; value <= sup ; value++) {
        newVariable[positionVar] = new int[2];
        string varnamel = "WSame_ValueLeft=" + to_string(value);
        string varnamer = "WSame_ValueRight=" + to_string(value);
        newVariable[positionVar][0] = wcsp->makeEnumeratedVariable(varnamel,0,arity/2);
        newVariable[positionVar][1] = wcsp->makeEnumeratedVariable(varnamer,0,arity/2);
        positionVar++;
    }

    // Adding WeightedAmong over each variable
    positionVar = 0;
    for (int value = inf ; value <= sup ; value++) {
        int* newScopeL = new int[arity/2+1];
        int* newScopeR = new int[arity/2+1];
        newScopeL[arity/2] = newVariable[positionVar][0];
        newScopeR[arity/2] = newVariable[positionVar][1];
        for (int variable = 0 ; variable < arity/2 ; ++variable) {
            newScopeL[variable] = scope[variable];
            newScopeR[variable] = scope[variable+arity/2];
        }

        WeightedVarAmong* wamongL = new WeightedVarAmong(arity/2+1,newScopeL);
        wamongL->setSemantics("hard");
        wamongL->setBaseCost(top);
        wamongL->addValue(value);
        wamongL->addToCostFunctionNetwork(wcsp);
        delete[] newScopeL;

        WeightedVarAmong* wamongR = new WeightedVarAmong(arity/2+1,newScopeR);
        wamongR->setSemantics("hard");
        wamongR->setBaseCost(top);
        wamongR->addValue(value);
        wamongR->addToCostFunctionNetwork(wcsp);
        delete[] newScopeR;

        positionVar++;
    }

    // Adding Binary constraints
    for (int value = 0 ; value < nbValue ; value++) {
        EnumeratedVariable* left 	= (EnumeratedVariable *) wcsp->getVar(newVariable[value][0]);
        EnumeratedVariable* right 	= (EnumeratedVariable *) wcsp->getVar(newVariable[value][1]);
        unsigned long tableSize = long (left->getDomainInitSize() * right->getDomainInitSize());
        vector<Cost> binaryCosts(tableSize,top);
        //cout << "Binary = " << newVariable[value][0] << " " <<  newVariable[value][1] << endl;
        for (unsigned long valueL = 0; valueL < left->getDomainInitSize() ; valueL++) {
            for (unsigned long valueR = 0; valueR < right->getDomainInitSize() ; valueR++) {
                unsigned long position = (valueR) * (left->getDomainInitSize()) + valueL;

                int gap =  (valueL - valueR);
                if (gap < 0) gap *=-1;
                Cost currentCost = 0;
                if (gap && semantics == "hard") currentCost = baseCost;
                if (semantics == "lin" || semantics == "var") 		currentCost = baseCost*gap;
                if (semantics == "quad") 		currentCost = baseCost*gap*gap;
                binaryCosts[position] = currentCost;

                //cout << valueL << " - " << valueR << " ("<< gap << ") => " << currentCost << endl;

            }
        }

        wcsp->postBinaryConstraint(newVariable[value][0], newVariable[value][1], binaryCosts);
    }

}

void 
WeightedSame::display() {
    cout << "WeightedSame (" << arity << ") : ";
    for (int variable = 0 ; variable < arity ; ++variable) {
        cout << scope[variable] << " ";
        if (variable == arity/2-1) cout << "| ";
    }
    cout << endl;
    cout << semantics << " " << baseCost << endl;
}


/// WEIGHTED SAMEGCC////////////////////////////////////////////////////

WeightedSameGcc::WeightedSameGcc() : DecomposableGlobalCostFunction() {}
WeightedSameGcc::WeightedSameGcc(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}
WeightedSameGcc::WeightedSameGcc(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
    file >> semantics >> baseCost;
    file >> nbValueToWatch;
    for (int idvalue = 0 ; idvalue < nbValueToWatch ; idvalue++) {
        Value value;
        unsigned int lb,ub;
        file >> value >> lb >> ub;
        setBounds(value,lb,ub);
    }
    if (_arity % 2 == 1 ) {
        cerr << "WeightedSameGcc::Constructor the scope must be even" << endl;
        exit(128);
    }
    //display();
}
WeightedSameGcc::~WeightedSameGcc() {}

void 
WeightedSameGcc::setBounds(Value value, unsigned int lb, unsigned int ub) {
    map<Value, pair<unsigned int, unsigned int> >::iterator it;
    it = bounds.find(value);
    if (it != bounds.end()) {
        cerr << "WeightedSameGcc::setBounds | Value " << value << " is already watch" << endl;
        exit(128);
    }
    bounds[value] = make_pair(lb,ub);
}

void 
WeightedSameGcc::addToCostFunctionNetwork(WCSP* wcsp) {
    Cost top = wcsp->getUb();

    // Counting the number of value
    int inf = ((EnumeratedVariable *) wcsp->getVar(scope[0]))->getInf();
    int sup = ((EnumeratedVariable *) wcsp->getVar(scope[0]))->getSup();
    for (int variable = 0 ; variable < arity ; ++variable) {
        int tinf = ((EnumeratedVariable *) wcsp->getVar(scope[variable]))->getInf();
        int tsup = ((EnumeratedVariable *) wcsp->getVar(scope[variable]))->getSup();
        inf = min ( inf , tinf );
        sup = max (sup , tsup );
    }
    int nbValue = sup - inf + 1;

    // Creating new counter variables
    int positionVar = 0;
    int** newVariable = new int*[nbValue];
    for (int value = inf ; value <= sup ; value++) {
        newVariable[positionVar] = new int[2];
        string varnamel = "WSame_ValueLeft=" + to_string(value);
        string varnamer = "WSame_ValueRight=" + to_string(value);
        newVariable[positionVar][0] = wcsp->makeEnumeratedVariable(varnamel,0,arity/2);
        newVariable[positionVar][1] = wcsp->makeEnumeratedVariable(varnamer,0,arity/2);
        positionVar++;
    }

    // Adding WeightedAmong over each variable
    positionVar = 0;
    for (int value = inf ; value <= sup ; value++) {
        int* newScopeL = new int[arity/2+1];
        int* newScopeR = new int[arity/2+1];
        newScopeL[arity/2] = newVariable[positionVar][0];
        newScopeR[arity/2] = newVariable[positionVar][1];
        for (int variable = 0 ; variable < arity/2 ; ++variable) {
            newScopeL[variable] = scope[variable];
            newScopeR[variable] = scope[variable+arity/2];
        }

        WeightedVarAmong* wamongL = new WeightedVarAmong(arity/2+1,newScopeL);
        wamongL->setSemantics("hard");
        wamongL->setBaseCost(top);
        wamongL->addValue(value);
        wamongL->addToCostFunctionNetwork(wcsp);
        delete[] newScopeL;

        WeightedVarAmong* wamongR = new WeightedVarAmong(arity/2+1,newScopeR);
        wamongR->setSemantics("hard");
        wamongR->setBaseCost(top);
        wamongR->addValue(value);
        wamongR->addToCostFunctionNetwork(wcsp);
        delete[] newScopeR;

        positionVar++;
    }

    // Adding Binary constraints (SAME PART)
    for (int value = 0 ; value < nbValue ; value++) {
        EnumeratedVariable* left 	= (EnumeratedVariable *) wcsp->getVar(newVariable[value][0]);
        EnumeratedVariable* right 	= (EnumeratedVariable *) wcsp->getVar(newVariable[value][1]);
        unsigned long tableSize = long (left->getDomainInitSize() * right->getDomainInitSize());
        vector<Cost> binaryCosts(tableSize,top);
        //cout << "Binary = " << newVariable[value][0] << " " <<  newVariable[value][1] << endl;
        for (unsigned long valueL = 0; valueL < left->getDomainInitSize() ; valueL++) {
            for (unsigned long valueR = 0; valueR < right->getDomainInitSize() ; valueR++) {
                unsigned long position = (valueR) * (left->getDomainInitSize()) + valueL;

                int gap =  (valueL - valueR);
                if (gap < 0) gap *=-1;
                Cost currentCost = 0;
                if (gap && semantics == "hard") currentCost = baseCost;
                if (semantics == "lin" || semantics == "dec") 		currentCost = baseCost*gap;
                if (semantics == "quad") 		currentCost = baseCost*gap*gap;
                binaryCosts[position] = currentCost;

                //cout << valueL << " - " << valueR << " ("<< gap << ") => " << currentCost << endl;

            }
        }

        wcsp->postBinaryConstraint(newVariable[value][0], newVariable[value][1], binaryCosts);
    }

    positionVar = 0;
    // Adding Unary Constraints (GCC PART)
    for (int value = inf ; value <= sup ; value ++) {

        map<Value, pair<unsigned int, unsigned int> >::iterator it;
        it = bounds.find(value);
        if (it != bounds.end()) {
            pair<Value, pair<unsigned int, unsigned int> > bound = *it;
            unsigned int lb = (bound.second).first;
            unsigned int ub = (bound.second).second;

            { //LEFT
                vector<Cost> unaryCosts(arity+1 ,0);
                for (int count=0 ; count <= arity ; count++) {
                    Cost currentCost = 0;
                    int gap = max( 0 , max( int(lb - count), int(count - ub) ) );
                    if (gap && semantics == "hard") currentCost = baseCost;
                    if (semantics == "lin" || semantics == "dec") 		currentCost = baseCost*gap;
                    if (semantics == "quad") 		currentCost = baseCost*gap*gap;
                    unaryCosts[count] = currentCost;
                }
                wcsp->postUnary(newVariable[positionVar][0],unaryCosts);
            }
            { //RIGHT
                vector<Cost> unaryCosts(arity+1 ,0);
                for (int count=0 ; count <= arity ; count++) {
                    Cost currentCost = 0;
                    int gap = max( 0 , max( int(lb - count), int(count - ub) ) );
                    if (gap && semantics == "hard") currentCost = baseCost;
                    if (semantics == "lin" || semantics == "dec") 		currentCost = baseCost*gap;
                    if (semantics == "quad") 		currentCost = baseCost*gap*gap;
                    unaryCosts[count] = currentCost;
                }
                wcsp->postUnary(newVariable[positionVar][1],unaryCosts);
            }
        }
        positionVar++;
    }

    delete[] newVariable;

}

void 
WeightedSameGcc::display() {
    cout << "WeightedSameGcc (" << arity << ") : ";
    for (int variable = 0 ; variable < arity ; ++variable) {
        cout << scope[variable] << " ";
        if (variable == arity/2-1) cout << "| ";
    }
    cout << endl;
    cout << semantics << " " << baseCost << endl;
    for (map<Value,pair <unsigned int,unsigned int> >::iterator it = bounds.begin(); it != bounds.end() ; ++it) {
        pair<Value,pair <unsigned int,unsigned int> > bound = *it;
        cout << bound.first << " [" << (bound.second).first << ":" << (bound.second).second << "]" << endl;
    }
}


/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

