
#include "tb2naryconstr.hpp"
#include "tb2vac.hpp"
#include "tb2clusters.hpp"






NaryConstraint::NaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval) : AbstractNaryConstraint(wcsp, scope_in, arity_in),
        nonassigned(arity_in, &wcsp->getStore()->storeInt)
{
    int i;
    Char* tbuf = new Char [arity_in+1];
    tbuf[arity_in] = '\0';

    for(i=0;i<arity_in;i++) {
        conflictWeights.push_back(0);
        int domsize = scope_in[i]->getDomainInitSize();
        tbuf[i] = CHAR_FIRST;
        if(domsize + CHAR_FIRST > MAX_CHAR) {
            cerr << "Nary constraints overflow. Try undefine NARYCHAR in makefile." << endl;
            exit(EXIT_FAILURE);
        }
    }
    iterTuple = String(tbuf);
    evalTuple = String(tbuf);
    delete [] tbuf;

    Cost Top = wcsp->getUb();
    default_cost = defval;
    if(default_cost > Top) default_cost = Top;
    store_top = default_cost < Top;
}



NaryConstraint::NaryConstraint(WCSP *wcsp) : AbstractNaryConstraint(wcsp), nonassigned(0, &wcsp->getStore()->storeInt)
{
}


Cost NaryConstraint::evalsubstr( String& s, Constraint* ctr )
{
    int count = 0;

    for(int i=0;i<arity();i++) {
        int ind = ctr->getIndex( getVar(i) );
        if(ind >= 0) { evalTuple[i] = s[ind]; count++; }
    }
    assert(count <= arity());

    Cost cost;
    if(count == arity()) cost = eval( evalTuple );
    else cost = MIN_COST;

    return cost;
}

// USED ONLY DURING SEARCH
void NaryConstraint::assign(int varIndex) {
    if (connected(varIndex)) {
        deconnect(varIndex);
        nonassigned = nonassigned - 1;

        if (size()<=4 && universal()) { // check if it is satisfied (clause)
            //	  cout << "nary cost function satisfied: " << this << endl;
            deconnect();
        }

        if(nonassigned <= 3) {
            //	  cout << "Assign var " << *getVar(varIndex) << "  in  " << *this;
            deconnect();
            projectNary();
        }
    }
}


void NaryConstraint::projectNaryTernary(TernaryConstraint* xyz)
{
    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) xyz->setCluster( cluster );
    EnumeratedVariable* x = (EnumeratedVariable*) xyz->getVar(0);
    EnumeratedVariable* y = (EnumeratedVariable*) xyz->getVar(1);
    EnumeratedVariable* z = (EnumeratedVariable*) xyz->getVar(2);
    TernaryConstraint* ctr = x->getConstr(y,z);
    if(ctr && td && (ctr->getCluster() != getCluster())) {
        TernaryConstraint* ctr_ = x->getConstr(y,z,getCluster());
        if(ctr_) ctr = ctr_;
    }
    if (ToulBar2::verbose >= 1) {
        cout << "project nary to ternary (" << x->wcspIndex << "," << y->wcspIndex << "," << z->wcspIndex << ") ";
        if(td) cout << "   cluster nary: " << getCluster() << endl; else cout << endl;
        if(ctr) cout << "ctr exists" << endl;
    }
    if(!ctr || (ctr && td && cluster != ctr->getCluster())) {
        xyz->fillElimConstrBinaries();
        xyz->reconnect();
        if(ctr) xyz->setDuplicate();
    } else {
        ctr->addCosts(xyz);
        xyz = ctr;
    }
    xyz->propagate();
    assert( !td || (xyz->getCluster() == xyz->xy->getCluster() &&  xyz->getCluster() == xyz->xz->getCluster() &&  xyz->getCluster() == xyz->yz->getCluster()) );
}

void NaryConstraint::projectNaryBinary(BinaryConstraint* xy)
{
    TreeDecomposition* td = wcsp->getTreeDec();
    EnumeratedVariable* x = (EnumeratedVariable*) xy->getVar(0);
    EnumeratedVariable* y = (EnumeratedVariable*) xy->getVar(1);

    if (ToulBar2::verbose >= 1) cout << "project nary to binary (" << x->wcspIndex << "," << y->wcspIndex << ")" << endl;

    BinaryConstraint* ctr = NULL;
    if (td) ctr = x->getConstr(y, getCluster());
    if (!ctr) ctr = x->getConstr(y);

    if((ctr && !td) || (ctr && td && (getCluster() == ctr->getCluster())))
    {
        if (ToulBar2::verbose >= 2) cout << " exists -> fusion" << endl;
        ctr->addCosts(xy);
        xy = ctr;
    }
    else {
        if(td) {
            if(ctr) xy->setDuplicate();
            xy->setCluster( getCluster() );
        }
    }
    if (x->unassigned() && y->unassigned()) xy->reconnect();
    xy->propagate();

    if (ToulBar2::verbose >= 2) cout << " and the result: " << *xy << endl;
}



// USED ONLY DURING SEARCH to project the nary constraint
void NaryConstraint::projectNary()
{
    wcsp->revise(this);
    int indexs[3];
    EnumeratedVariable* unassigned[3] = {NULL,NULL,NULL};
    Char* tbuf = new Char [arity_ + 1];
    String t;
    bool flag = false;

    int i,nunassigned = 0;
    for(i=0;i<arity_;i++) {
        EnumeratedVariable* var = (EnumeratedVariable*) getVar(i);
        if(getVar(i)->unassigned()) {
            unassigned[nunassigned] = var;
            indexs[nunassigned] = i;
            tbuf[i] = CHAR_FIRST;
            nunassigned++;
        } else tbuf[i] = var->toIndex(var->getValue()) + CHAR_FIRST;
    }
    tbuf[arity_] =  '\0';
    t = tbuf;
    delete [] tbuf;

    EnumeratedVariable* x = unassigned[0];
    EnumeratedVariable* y = unassigned[1];
    EnumeratedVariable* z = unassigned[2];

    assert(nunassigned <= 3);
    if(nunassigned == 3) {
        TernaryConstraint* xyz = wcsp->newTernaryConstr(x,y,z,this);
        wcsp->elimTernOrderInc();
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
                    Value xval = *iterx;
                    Value yval = *itery;
                    Value zval = *iterz;
                    t[indexs[0]] =  x->toIndex(xval) + CHAR_FIRST;
                    t[indexs[1]] =  y->toIndex(yval) + CHAR_FIRST;
                    t[indexs[2]] =  z->toIndex(zval) + CHAR_FIRST;
                    Cost curcost = eval(t);
                    if (curcost > MIN_COST) flag = true;
                    xyz->setcost(x,y,z,xval,yval,zval,curcost);
                }}}
        if(flag) projectNaryTernary(xyz);
        //else cout << "ternary empty!" << endl;
    }
    else if(nunassigned == 2) {
        BinaryConstraint* xy = wcsp->newBinaryConstr(x,y,this);
        wcsp->elimBinOrderInc();
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                Value xval = *iterx;
                Value yval = *itery;
                t[indexs[0]] =  x->toIndex(xval) + CHAR_FIRST;
                t[indexs[1]] =  y->toIndex(yval) + CHAR_FIRST;
                Cost curcost = eval(t);
                if (curcost > MIN_COST) flag = true;
                xy->setcost(xval,yval,curcost);
                if (ToulBar2::verbose >= 5) {
                    Cout << t;
                    cout << " " << curcost << endl;
                }
            }}
        if(flag)projectNaryBinary(xy);
        //else cout << "binary empty!" << endl;
    }
    else if(nunassigned == 1) {
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            Value xval = *iterx;
            t[indexs[0]] =  x->toIndex(xval) + CHAR_FIRST;
            Cost c = eval(t);
            if(c > MIN_COST) {
                if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
                    TreeDecomposition* td = wcsp->getTreeDec();
                    if(td) td->addDelta(cluster,x,xval,c);
                }
                x->project(xval, c);
            }
        }
        x->findSupport();
    }
    else {
        Cost c = eval(t);
        projectLB(c);
    }
}


void NaryConstraint::projectFromZero(int index)
{
    return; // TO BE REMOVED !!!

    int i;
    int nsup = 0;
    Cost minc = MAX_COST;
    Cost c;
    EnumeratedVariable* var = (EnumeratedVariable*) getVar(index);
    for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
        c = var->getCost(*iter);
        if(c == MIN_COST) { if(nsup) return; nsup++; }
        else if(c < minc) minc = c;

    }
    for(i=0;i<arity_;i++) {
        if(i != index) {
            var = (EnumeratedVariable*) getVar(i);
            if(getVar(i)->unassigned()) {
                nsup = 0;
                for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
                    c = var->getCost(*iter);
                    if(c == MIN_COST) { if(nsup) return; nsup++; }
                    else if(c < minc) minc = c;
                }
            }
        }
    }
    int a = arity();
    Char* cht = new Char [a + 1];
    cht[a] = '\0';
    for(i=0;i<a;i++) {
        var = (EnumeratedVariable*) getVar(i);
        cht[i] = var->toIndex(var->getSup()) + CHAR_FIRST;
    }
    String t(cht);
    delete [] cht;
    c = eval(t);
    if(c < minc) minc = c;
    if(c > MIN_COST) starrule(t,minc);
}



void NaryConstraint::starrule(String& t, Cost minc)
{
    int i;
    for(i=0;i<arity_;i++) {
        EnumeratedVariable* var = (EnumeratedVariable*) getVar(i);
        if(getVar(i)->unassigned()) {
            for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
                Cost c = var->getCost(*iter);
                if(c > MIN_COST) var->extend(*iter, minc);
            }
        }

    }
    setTuple( t, eval(t) - minc );
    projectLB(minc);
}




void NaryConstraint::firstlex()
{
    it_values.clear();
    EnumeratedVariable* var;
    for(int i=0;i<arity_;i++) {
        var = (EnumeratedVariable*) getVar(i);
        it_values.push_back( var->begin() );
    }
}

bool NaryConstraint::nextlex( String& t, Cost& c)
{
    int i;
    int a = arity();
    EnumeratedVariable* var = (EnumeratedVariable*) getVar(0);
    if(it_values[0] == var->end()) return false;

    for(i=0;i<a;i++) {
        var = (EnumeratedVariable*) getVar(i);
        iterTuple[i] = var->toIndex(*it_values[i]) + CHAR_FIRST;
    }
    t = iterTuple;
    c = eval(t);

    // and now increment
    bool finished = false;
    i = a-1;
    while(!finished) {
        var = (EnumeratedVariable*) getVar(i);
        ++it_values[i];
        finished = it_values[i] != var->end();
        if(!finished) {
            if(i>0) {
                it_values[i] = var->begin();
                i--;
            } else finished = true;
        }
    }
    return true;
}





// ********************* NaryConstraintMap *********************
// _____________________________________________________________



NaryConstraintMap::NaryConstraintMap(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval)
: NaryConstraint(wcsp, scope_in, arity_in, defval)
{
    pf = new TUPLES;
    filters = NULL;
    // Cannot propagate here because cost tuples are not known yet
}


NaryConstraintMap::NaryConstraintMap(WCSP *wcsp)
: NaryConstraint(wcsp)
{
    pf = new TUPLES;
    filters = NULL;
}


NaryConstraintMap::~NaryConstraintMap()
{
    if(pf) delete pf;
}






void NaryConstraintMap::resetFilters()
{
    if(filters) {
        delete filters;
        filters = NULL;
    }
}

void NaryConstraintMap::fillFilters()
{
    if(filters == NULL) {
        filters = new set<Constraint*>;
        for(int i=0;i<arity();i++) {
            EnumeratedVariable* v = (EnumeratedVariable*) getVar(i);
            for (ConstraintList::iterator iter=v->getConstrs()->begin(); iter != v->getConstrs()->end(); ++iter) {
                Constraint* ctr = (*iter).constr;
                if(scopeIncluded(ctr)) filters->insert( ctr );
            }
        }
    }
}

Cost NaryConstraintMap::eval( String& tin, EnumeratedVariable** scope_in ) {
    String t(tin);
    if(scope_in) {
        for(int i = 0; i < arity_; i++) {
            int pos = getIndex(scope_in[i]);
            t[pos] = tin[i];
        }
    }
    return eval(t);
}

Cost NaryConstraintMap::eval( String& s ) {
    TUPLES& f = *pf;
    TUPLES::iterator  it = f.find(s);
    if(it != f.end()) return it->second;
    else return default_cost;
}


/// Set new default cost to df (df <= Top), keep existing costs SMALLER than this default cost in a Map
void NaryConstraintMap::keepAllowedTuples( Cost df )
{
    assert( CUT(wcsp->getUb(), df) );
    TUPLES* pfnew = new TUPLES;

    String t;
    Cost c;
    firstlex();
    while(nextlex(t,c)) {
        if(c < df) (*pfnew)[t] = c;
    }
    default_cost = df;

    delete pf;
    pf = pfnew;
}




bool NaryConstraintMap::consistent( String& t ) {
    int a = arity();
    bool ok = true;
    for(int i=0;i<a && ok;i++) {
        EnumeratedVariable* var = (EnumeratedVariable*) getVar(i);
        ok = ok && var->canbe( var->toValue(t[i] - CHAR_FIRST) );
    }
    return ok;
}


void NaryConstraintMap::first()
{
    tuple_it = pf->begin();
}

bool NaryConstraintMap::next( String& t, Cost& c)
{
    bool ok = false;
    while(!ok && (tuple_it != pf->end())) {
        t = tuple_it->first;
        c = tuple_it->second;
        ok = consistent(t);
        tuple_it++;
    }
    if(!ok && tuple_it == pf->end()) return false;
    else return true;
}

void NaryConstraintMap::first(EnumeratedVariable * vx, EnumeratedVariable * vz)
{
    it_values.clear();
    EnumeratedVariable* var;
    it_values.push_back( vx->begin() );
    for(int i=0;i<arity_;i++) {
        var = (EnumeratedVariable*) getVar(i);
        if(var != vx && var != vz) it_values.push_back( var->begin() );
    }
    it_values.push_back( vz->begin() );
}

bool NaryConstraintMap::separability(EnumeratedVariable * vx, EnumeratedVariable * vz){
    int a = arity(), i = a-1, k;
    bool finished = false;
    bool neweq = true;
    bool sev = true; // false if vx and vz are not separable
    Cost diff = 0;
    Cost c1, c;
    vector<int> ordre(a,-1);
    EnumeratedVariable* var;
    EnumeratedVariable *scope_in[a];
    String t1(iterTuple), t(iterTuple);

    EnumeratedVariable::iterator itvzfirst = vz->begin();
    EnumeratedVariable::iterator itvznext = itvzfirst;
    if(itvznext !=vz->end()) ++itvznext;
    first(vx,vz);

    scope_in[0] = vx; //.push_back( alpha );
    ordre[a-1] = 0;
    k = 1;
    if(ToulBar2::verbose >= 1) cout << "[ " << vx->getName() << " ";
    for (int j = 0; j < a; j++) {
        var = (EnumeratedVariable*) getVar(j);
        if(var != vx && var != vz) {
            scope_in[k] = var;
            ordre[k-1] = k;
            if(ToulBar2::verbose >= 1) cout << var->getName() << " ";
            ++k;
        }
    }

    scope_in[a-1] = vz;//.push_back( beta );
    ordre[a-2] = a-1;
    //if(it_values[ordre[a-2]] != vz->end()) ++it_values[ordre[a-2]];
    if(ToulBar2::verbose >= 1) cout << vz->getName() << " ] ? ";// << endl;

    while( sev &&  itvznext != vz->end()){
        first(vx,vz);
        while(sev && it_values[ordre[0]] != scope_in[ordre[0]]->end()){
            if(it_values[ordre[a-2]] == itvzfirst/*vz->begin()*/) ++it_values[ordre[a-2]];
            //if(it_values[ordre[a-2]] != vz->end()){
            for(int j = 0; j < a; j++) {
                t[j] = scope_in[j]->toIndex(*it_values[j]) + CHAR_FIRST; // PROBLEME CAR it_values[ordre[a-2]] = vz->end()
                t1[j] = scope_in[j]->toIndex(*it_values[j]) + CHAR_FIRST;
            }
            t1[a-1] = scope_in[a-1]->toIndex(*itvzfirst/*vz->begin()*/) + CHAR_FIRST;

            // etude de la difference
            c1 = eval(t1,scope_in);
            c = eval(t,scope_in);
            if(!universe(c, c1, wcsp->getUb())){

                if(ToulBar2::verbose >= 3) {if(!neweq) cout << "= \n"; else cout << endl;}
                if(ToulBar2::verbose >= 3){
                    cout << "C";
                    for(unsigned int j=0;j<t.size();j++) {
                        cout << t[j] - CHAR_FIRST << "." ;
                    }
                    cout << " - C";
                    for(unsigned int j=0;j<t1.size();j++) {
                        cout << t1[j] - CHAR_FIRST << "." ;
                    }
                    cout << " = " << c << " - " << c1;
                }
                if(neweq) {diff = squareminus(c,c1,wcsp->getUb()); neweq = false; }
                else sev = (diff == squareminus(c,c1,wcsp->getUb()));
                if(ToulBar2::verbose >= 3) cout << " = " << squareminus(c,c1,wcsp->getUb())  << endl;
            }
            else if(ToulBar2::verbose >= 3) cout << "universe\n";
            finished = false;
            i = a-1;
            while(sev && !finished) {
                var = scope_in[ordre[i]];
                ++it_values[ordre[i]];
                finished = it_values[ordre[i]] != var->end();
                if(!finished) {
                    //it_values[ordre[i]] = var->begin();
                    if(i>0) {
                        it_values[ordre[i]] = var->begin();
                        if(i == a-1) neweq = true;
                        i--;
                    }else finished = true;
                }
            }
            if(it_values[ordre[a-2]] == itvzfirst/*vz->begin()*/ && it_values[ordre[a-2]] != vz->end()) ++it_values[ordre[a-2]];
            //}
        }
        ++itvzfirst;
        ++itvznext;
        if(ToulBar2::verbose >= 3) cout << "--\n";
    }
    return sev;
}

void NaryConstraintMap::separate(EnumeratedVariable *vx, EnumeratedVariable *vz)
{
    Cost cost,minCost = MAX_COST;
    int a = arity();
    String t(a,'0'),tX(a-1,'0'),tZ(a-1,'0');
    int index, k;
    EnumeratedVariable* var;
    TernaryConstraint *existX = NULL, *existZ = NULL;
    Constraint *naryz,  *naryx;
    EnumeratedVariable * scope_in[a];
    EnumeratedVariable * subscopeX[a-1];
    EnumeratedVariable * subscopeZ[a-1];
    int scopeX[a-1];
    int scopeZ[a-1];

    first(vx,vz);

    // initialisation des scopes
    scope_in[0] = vx;
    scopeX[0] = vx->wcspIndex;
    subscopeX[0] = vx;

    k = 1;

    for (int j = 0; j < a; j++) {
        var = (EnumeratedVariable*) getVar(j);
        if(var != vx && var != vz) {
            scope_in[k] = var;

            subscopeX[k] = var;
            scopeX[k] = var->wcspIndex;

            subscopeZ[k-1] = var;
            scopeZ[k-1] = var->wcspIndex;

            ++k;
        }
    }

    scope_in[a-1] = vz;

    scopeZ[a-2] = vz->wcspIndex;
    subscopeZ[a-2] = vz;

    // creation de la nouvelle contrainte
    if(a == 4){
        int size = subscopeX[0]->getDomainInitSize()*subscopeX[1]->getDomainInitSize()*subscopeX[2]->getDomainInitSize();
        vector<Cost> costs(size,0);
        existX = subscopeX[0]->getConstr(subscopeX[1],subscopeX[2]);
        naryx = wcsp->newTernaryConstr(subscopeX[0],subscopeX[1],subscopeX[2],costs);
    }
    else{
        index = wcsp->postNaryConstraintBegin(scopeX,a-1,wcsp->getUb());
        naryx = wcsp->getCtr(index);
    }

    minCost = MAX_COST;
    int i;
    if(ToulBar2::verbose >= 3) {
        cout << "\n[ ";
        for(int j = 0; j < a-1; j++)
            cout << scopeX[j] << " ";
        cout << " ]"  << endl;
    }
    while(it_values[0] != scope_in[0]->end()){
        for(int j = 0; j < a-1; j++) {
            t[j] = scope_in[j]->toIndex(*it_values[j]) + CHAR_FIRST;
            tX[j] = scope_in[j]->toIndex(*it_values[j]) + CHAR_FIRST;
        }
        t[a-1] = scope_in[a-1]->toIndex(*it_values[a-1]) + CHAR_FIRST;

        cost = eval(t, scope_in);
        if(cost < minCost) minCost = cost;
        int finished = false;
        i = a-1;
        while( !finished) {

            ++it_values[i];
            finished = it_values[i] != scope_in[i]->end();
            if(!finished) {
                if(i == a-1) {
                    if(ToulBar2::verbose >= 3){
                        for(unsigned int j=0;j<tX.size();j++) {
                            cout << tX[j] - CHAR_FIRST<< " ";
                        }
                        cout << "  " << minCost << endl;
                    }
                    naryx->setTuple(tX, minCost,subscopeX);
                    minCost = MAX_COST;
                }
                if(i>0) {
                    it_values[i] = scope_in[i]->begin();
                    i--;
                }else finished = true;
            }
        }

    }


    if(ToulBar2::verbose >= 3) cout << "-----------------------------------" << endl;
    if(a == 4){

        int size = subscopeZ[0]->getDomainInitSize()*subscopeZ[1]->getDomainInitSize()*subscopeZ[2]->getDomainInitSize();
        vector<Cost> costs(size,0);
        existZ = subscopeZ[0]->getConstr(subscopeZ[1],subscopeZ[2]);
        naryz = wcsp->newTernaryConstr(subscopeZ[0],subscopeZ[1],subscopeZ[2],costs);
    }
    else {
        index = wcsp->postNaryConstraintBegin(scopeZ,a-1,wcsp->getUb());
        naryz = wcsp->getCtr(index);
    }

    if(ToulBar2::verbose >= 3){
        cout << "[ ";
        for(int j = 0; j < a-1; j++)
            cout << scopeZ[j] << " ";
        cout << " ]"  << endl;
    }
    first(vx,vz);

    Cost diffcost;
    Cost cX;
    while(it_values[1] != scope_in[1]->end()){
        for(int j = 1; j < a; j++) {
            t[j] = scope_in[j]->toIndex(*it_values[j]) + CHAR_FIRST;
            tZ[j-1] = scope_in[j]->toIndex(*it_values[j]) + CHAR_FIRST;
        }
        //t[0] = scope_in[0]->toIndex(*it_values[0]) + CHAR_FIRST;
        EnumeratedVariable::iterator it = scope_in[0]->begin();
        do{
            t[0] = scope_in[0]->toIndex(*it) + CHAR_FIRST;
            cost = eval(t, scope_in);
            cX = naryx->evalsubstr(t,naryx);
            diffcost = squareminus(cost,cX, wcsp->getUb());
            ++it;
            //cout << t[0]-CHAR_FIRST << endl;
        }while(it != scope_in[0]->end() && cost >= wcsp->getUb() && cX >= wcsp->getUb());

        cost = eval(t, scope_in);
        Cost cX = naryx->evalsubstr(t,naryx);
        diffcost = squareminus(cost,cX, wcsp->getUb());
        if(ToulBar2::verbose >= 3){
            for(unsigned int j=0;j<tZ.size();j++) {
                cout << tZ[j] - CHAR_FIRST<< " ";
            }

            cout << "  " <<  diffcost << endl;
        }
        naryz->setTuple(tZ, diffcost,subscopeZ);
        int finished = false;
        i = a-1;
        while( !finished) {

            ++it_values[i];
            finished = it_values[i] != scope_in[i]->end();
            if(!finished) {
                if(i>1) {
                    it_values[i] = scope_in[i]->begin();
                    i--;
                }else finished = true;
            }
        }

    }

    if(ToulBar2::verbose >= 3) cout << "-----------------------------------" << endl;
    assert(verifySeparate(naryx,naryz));

    if(!naryx->universal()){
        if(ToulBar2::verbose == 1) {
            cout << "\n[ ";
            for(int j = 0; j < a-1; j++)
                cout << scopeX[j] << " ";
            cout << " ]"  << endl;
        }
        if(a == 4){
            if(existX && !existX->universal()){
                ((TernaryConstraint*) naryx)->addCosts(existX);
                existX->deconnect();
            }
        }
        naryx->reconnect();
        naryx->propagate();
    } else naryx->deconnect();
    if(!naryz->universal()){
        if(ToulBar2::verbose == 1){
            cout << "[ ";
            for(int j = 0; j < a-1; j++)
                cout << scopeZ[j] << " ";
            cout << " ]"  << endl;
        }
        if(a == 4){
            if(existZ && !existZ->universal()){
                ((TernaryConstraint*) naryz)->addCosts(existZ);
                existZ->deconnect();
            }
        }
        naryz->reconnect();
        naryz->propagate();
    } else naryz->deconnect();
    deconnect();
}

void NaryConstraintMap::permute( EnumeratedVariable** scope_in )
{
    TUPLES* pf_old = pf;
    pf = new TUPLES;

    TUPLES::iterator it = pf_old->begin();
    while(it != pf_old->end()) {
        String s(it->first);
        setTuple(s, it->second, scope_in );
        it++;
    }

    delete pf_old;
    for(int i=0; i<arity_; i++)
    {
        map<int,int>::iterator it_pos =  scope_inv.find(scope_in[i]->wcspIndex);
        int i_old = it_pos->second;
        it_pos->second = i;

        scope_inv[ scope_in[i]->wcspIndex ] = i;
        scope[i] = scope_in[i];

        DLink<ConstraintLink>* l = links[i];
        links[i] = links[i_old];
        links[i_old] = l;

        links[i]->content.scopeIndex = i;
        l->content.scopeIndex = i_old;
    }
}


// for adding a tuple in f
// scope_in contains the order of the values in String tin
void NaryConstraintMap::setTuple( String& tin, Cost c, EnumeratedVariable** scope_in )
{
    String t(tin);
    if(scope_in) {
        for(int i = 0; i < arity_; i++) {
            int pos = getIndex(scope_in[i]);
            t[pos] = tin[i];
        }
    }
    (*pf)[t] = c;
}

void NaryConstraintMap::addtoTuple( String& tin, Cost c, EnumeratedVariable** scope_in )
{
    String t(tin);
    if(scope_in) {
        for(int i = 0; i < arity_; i++) {
            int pos = getIndex(scope_in[i]);
            t[pos] = tin[i];
        }
    }
    (*pf)[t] += c;
}

void NaryConstraintMap::setInfiniteCost(Cost ub)
{
    Cost mult_ub = ((ub < (MAX_COST / MEDIUM_COST))?(max(LARGE_COST, ub * MEDIUM_COST)):ub);
    for (TUPLES::iterator it = pf->begin(); it != pf->end(); ++it) {
        Cost c =  it->second;
        if (CUT(c, ub)) it->second = mult_ub;
    }
    if (CUT(default_cost, ub)) default_cost = mult_ub;
}

void NaryConstraintMap::insertSum( String& t1, Cost c1, Constraint* ctr1, String t2, Cost c2, Constraint* ctr2, bool bFilters )
{
    Cost Top = wcsp->getUb();
    if(c1 >= Top) return;
    if(c2 >= Top) return;
    Cost csum = c1 + c2;

    Char* t = new Char [arity_+1];

    for(int i = 0; i < arity_; i++) {
        EnumeratedVariable* v = scope[i];
        int pos = i;
        int pos1 = ctr1->getIndex(v);
        int pos2 = ctr2->getIndex(v);

        if((pos1 >= 0) && (pos2 >= 0)) {
            if(t1[pos1] != t2[pos2])  { delete [] t; return; }
            t[pos] = t1[pos1];
        }
        else if(pos1 >= 0) t[pos] = t1[pos1];
        else if(pos2 >= 0) t[pos] = t2[pos2];

        Cost unaryc = v->getCost( v->toValue(t[pos] - CHAR_FIRST)  );
        if(unaryc >= Top) return;
        csum += unaryc;
        if(csum >= Top) return;
    }
    t[arity_] = '\0';
    String tstr(t);

    if(bFilters && filters && (default_cost >= Top) ) {
        set<Constraint*>::iterator it = filters->begin();
        while(it != filters->end()) {
            Constraint* ctr = *it;
            if(ctr->connected()) {
                Cost c = ctr->evalsubstr(tstr, (Constraint*)this );
                if(c >= Top) return;
                csum += c;
            }
            if(csum >= Top) return;
            ++it;
        }
    }

    (*pf)[tstr] = c1 + c2;
    delete [] t;
}

// THIS CODE IS NEVER USED!!!
void NaryConstraintMap::sum( NaryConstraintMap* nary )
{
    deconnect(true);

    map<int,int> snew;
    set_union( scope_inv.begin(), scope_inv.end(),
            nary->scope_inv.begin(), nary->scope_inv.end(),
            inserter(snew, snew.begin()) );

    arity_ = snew.size();
    EnumeratedVariable** scope1 = scope;
    DLink<ConstraintLink>** links1 = links;
    scope = new EnumeratedVariable* [arity_];
    links = new DLink<ConstraintLink>* [arity_];

    int i = 0;
    map<int,int>::iterator its = snew.begin();
    while(its != snew.end()) {
        EnumeratedVariable* var = (EnumeratedVariable*) wcsp->getVar(its->first);
        its->second = i;
        scope[i] =  var;
        int index1 = getIndex(var);
        if(index1 >= 0) {
            links[i] = links1[index1];
            ConstraintLink e = {this, i};
            links[i]->content = e;
        }
        else links[i] = nary->links[ nary->getIndex(var) ];

        i++;
        its++;
    }

    TUPLES& f1 = *pf;
    TUPLES& f2 = *nary->pf;
    TUPLES::iterator  it1 = f1.begin();
    TUPLES& f = * new TUPLES;
    pf = &f;

    String t1,t2;
    Cost c1,c2;
    while(it1 != f1.end()) {
        t1 = it1->first;
        c1 =  it1->second;
        TUPLES::iterator  it2 = f2.begin();
        while(it2 != f2.end()) {
            t2 = it2->first;
            c2 =  it2->second;
            insertSum(t1, c1, this, t2, c2, nary);
            it2++;
        }
        it1++;
    }

    scope_inv = snew;
    delete [] scope1;
    delete [] links1;

    reconnect();
}

// Projection of variable x of the nary constraint
// complexity O(2|f|)
// this function is independent of the search
void NaryConstraintMap::project( EnumeratedVariable* x )
{
    int xindex = getIndex(x);
    if(xindex < 0) return;
    assert(x->getDegree() == 1);
    String t,tnext,tproj;
    Cost c;
    Value val;
    Cost Top = wcsp->getUb();
    TUPLES& f = *pf;
    TUPLES fproj;
    TUPLES::iterator  it;
    // First part of the projection: complexity O(|f|) we swap positions between the projected variable and the last variable
    while(!f.empty()) {
        it = f.begin();
        t = it->first;
        c =  it->second;
        assert( x->getDegree() == 1);
        val = x->toValue(t[xindex] - CHAR_FIRST);
        c += ((x->canbe(val))?(x->getCost( val )):Top);
        if(c > Top) c = Top;
        String tswap(t);
        Char a = tswap[arity_-1];
        tswap[arity_-1] = tswap[xindex];
        tswap[xindex] = a;
        fproj[tswap] = c;
        f.erase(t);
    }

    if(xindex < arity_-1) {
        // swap of links
        DLink<ConstraintLink>* linkx = links[xindex];
        links[xindex] = links[arity_-1];
        links[arity_-1] = linkx;
        //swap of scope array
        scope[ xindex ] = scope[arity_-1];
        scope[ arity_-1 ] = x;
        // update of links indexs
        links[xindex]->content.scopeIndex = xindex;
        links[arity_-1]->content.scopeIndex = arity_-1;
        scope_inv[scope[xindex]->wcspIndex] = xindex;
    }

    // Second part of the projection: complexity O(|f|) as the projected variable is in the last position,
    // it is sufficient to look for tuples with the same arity-1 prefix. If there are less than d (domain of
    // the projected variable) tuples, we have also to perform the minimum with default_cost
    // this is only true when the tuples are LEXICOGRAPHICALY ordered
    Cost negcost = 0;
    set<Value> markValue;
    it = fproj.begin();
    if(it != fproj.end()) {
        t = it->first;
        c = it->second;
        val = x->toValue(t[arity_-1] - CHAR_FIRST);
        bool end = false;
        unsigned int ntuples = ((x->canbe(val))?1:0);
        markValue.clear();
        markValue.insert(val);

        while(!end) {
            it++;
            end = (it == fproj.end());
            bool sameprefix = false;

            Cost cnext = MAX_COST;
            if(!end) {
                tnext = it->first;
                cnext = it->second;
                sameprefix = (t.compare(0,arity_-1,tnext,0,arity_-1) == 0);
                //cout << "<" << t << "," << c << ">   <" << tnext << "," << cnext << ">       : " << t.compare(0,arity_-1,tnext,0,arity_-1) << endl;
            }
            if(!sameprefix) {
                if(ntuples < x->getDomainSize()) {
                    for (EnumeratedVariable::iterator itv = x->begin(); itv != x->end(); ++itv) {
                        if (markValue.find(*itv) == markValue.end()) {
                            if (ToulBar2::isZ) {
                                c = wcsp->LogSumExp(c, default_cost + x->getCost(*itv));
                            } else if(default_cost + x->getCost(*itv) < c) c = default_cost + x->getCost(*itv);
                        }
                    }
                }
                if (ToulBar2::isZ && c < negcost) negcost = c;
                if(c != default_cost || ToulBar2::isZ) f[ t.substr(0,arity_-1) ] = c;
                t = tnext;
                c = cnext;
                val = x->toValue(t[arity_-1] - CHAR_FIRST);
                ntuples = ((x->canbe(val))?1:0);
                markValue.clear();
                markValue.insert(val);
            } else {
                val = x->toValue(t[arity_-1] - CHAR_FIRST);
                markValue.insert(val);
                if (x->canbe(val)) ntuples++;
                if (ToulBar2::isZ) {
                    c = wcsp->LogSumExp(c, cnext);
                } else if(cnext < c) c = cnext;
            }
        }
    }
    assert(negcost <= 0);
    if (negcost < 0) {
        for (it = f.begin(); it != f.end(); ++it) {
            t = (*it).first;
            c = (*it).second;
            f[t] = c - negcost;
        }
        default_cost -= negcost;
        assert(default_cost >= MIN_COST); // THIS TEST COULD BE REMOVED BUT NEGATIVE DEFAULT COST EFFECTS ARE UNSPECIFIED
        wcsp->decreaseLb(negcost);
    }
    if(x->unassigned()) {
        x->deconnect(links[arity_-1]);
        nonassigned = nonassigned - 1;
    }
    scope_inv.erase(x->wcspIndex);
    arity_--;
    iterTuple.resize(arity_);
    evalTuple.resize(arity_);
}


// Projects out all variables except x,y,z
// and gives the result at fproj
void NaryConstraintMap::projectxyz( EnumeratedVariable* x,
        EnumeratedVariable* y,
        EnumeratedVariable* z,
        TUPLES& fproj)
{
    assert(CUT(default_cost,wcsp->getUb()));

    Char   stxyz[4] = {CHAR_FIRST, CHAR_FIRST, CHAR_FIRST, '\0'};
    String txyz(stxyz);
    String t;
    Cost c;
    TUPLES::iterator  itproj;

    // compute in one pass of all tuples the projection
    first();
    while(next(t,c)) {
        txyz[0] = t[ getIndex(x) ];
        txyz[1] = t[ getIndex(y) ];
        txyz[2] = t[ getIndex(z) ];

        itproj = fproj.find(txyz);
        if(itproj != fproj.end()) {
            if(c < itproj->second) fproj[txyz] = c;
        } else {
            fproj[txyz] = c;
        }
    }

    // finially we substract the projection from the initial function
    first();
    while(next(t,c)) {
        txyz[0] = t[ getIndex(x) ];
        txyz[1] = t[ getIndex(y) ];
        txyz[2] = t[ getIndex(z) ];
        itproj = fproj.find(txyz);
        if(itproj != fproj.end()) { assert(CUT(c, itproj->second)); (*pf)[t] -= itproj->second; }
        else assert(false);
    }
}



// Projects out all variables except x,y
// and gives the result at fproj
void NaryConstraintMap::projectxy( EnumeratedVariable* x,
        EnumeratedVariable* y,
        TUPLES& fproj)
{
    assert(CUT(default_cost,wcsp->getUb()));

    Char   stxy[3] = {CHAR_FIRST, CHAR_FIRST, '\0'};
    String txy(stxy);
    String t;
    Cost c;
    TUPLES::iterator  itproj;

    // compute in one pass of all tuples the projection
    first();
    while(next(t,c)) {
        txy[0] = t[ getIndex(x) ];
        txy[1] = t[ getIndex(y) ];

        itproj = fproj.find(txy);
        if(itproj != fproj.end()) {
            if(c < itproj->second) fproj[txy] = c;
        } else {
            fproj[txy] = c;
        }
    }

    // finally we substract the projection from the initial function
    first();
    while(next(t,c)) {
        txy[0] = t[ getIndex(x) ];
        txy[1] = t[ getIndex(y) ];
        itproj = fproj.find(txy);
        if(itproj != fproj.end()) {
            assert(CUT(c, itproj->second));
            //		  if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
            (*pf)[t] -= itproj->second;
            //		  }
        } else assert(false);
    }
}


void NaryConstraintMap::preproject3()
{
    assert(connected());
    assert(CUT(default_cost,wcsp->getUb()));

    for(int i = 0; i < arity_ - 2; i++) {
        EnumeratedVariable* x = scope[i];
        EnumeratedVariable* y = scope[i+1];
        EnumeratedVariable* z = scope[i+2];

        TUPLES fproj;
        projectxyz(x,y,z,fproj);

        String t;
        vector<Cost> xyz;
        unsigned int a,b,c;
        unsigned int sizex = x->getDomainInitSize();
        unsigned int sizey = y->getDomainInitSize();
        unsigned int sizez = z->getDomainInitSize();

        for (a = 0; a < sizex; a++)
            for (b = 0; b < sizey; b++)
                for (c = 0; c < sizez; c++) xyz.push_back(default_cost);

        TUPLES::iterator it =  fproj.begin();
        while(it != fproj.end()) {
            t = it->first;
            a = t[0] - CHAR_FIRST;
            b = t[1] - CHAR_FIRST;
            c = t[2] - CHAR_FIRST;
            xyz[ a * sizey * sizez + b * sizez + c ]	= it->second;
            it++;
        }
        if(fproj.size() > 0 || default_cost > MIN_COST) wcsp->postTernaryConstraint(x->wcspIndex, y->wcspIndex, z->wcspIndex,xyz);
        if (deconnected()) return;
    }
}


inline bool cmp_pairvars(pair<EnumeratedVariable* ,EnumeratedVariable* > pv1, pair<EnumeratedVariable* ,EnumeratedVariable* > pv2) 
{ 
    return (pv1.first->wcspIndex < pv2.first->wcspIndex ||
            (pv1.first->wcspIndex == pv2.first->wcspIndex &&
                    pv1.second->wcspIndex < pv2.second->wcspIndex));
}
void NaryConstraintMap::preprojectall2()
{
    assert(connected());
    assert(CUT(default_cost,wcsp->getUb()));
    //  cout << "PREPROJECT " << *this << endl;

    // vector< pair<EnumeratedVariable* ,EnumeratedVariable* > > listxy;
    //  vector< pair<EnumeratedVariable* ,EnumeratedVariable* > > rndlistxy;

    // for(int i = 0; i < arity_; i++) {
    // for(int j = i+1; j < arity_; j++) {
    // 	pair<EnumeratedVariable* ,EnumeratedVariable* > p = make_pair(((scope[i]->wcspIndex < scope[j]->wcspIndex)?scope[i]:scope[j]),
    // 																  ((scope[i]->wcspIndex < scope[j]->wcspIndex)?scope[j]:scope[i]));
    // 	listxy.push_back(p);
    // }}
    // stable_sort(listxy.begin(), listxy.end(), cmp_pairvars);

    // while (listxy.size() > 0) {
    // 	int pos = myrand() % listxy.size();
    // 	//	cout << pos;
    // 	vector< pair<EnumeratedVariable* ,EnumeratedVariable* > >::iterator it = listxy.begin();
    // 	while (pos > 0) {
    // 	  ++it;
    // 	  --pos;
    // 	}
    // 	//	cout << "," << (*it).first << "," << (*it).second << endl;
    // 	rndlistxy.push_back(*it);
    // 	listxy.erase(it);
    // }

    // for(int i = 0; i < arity_; i++) {
    //   for(int j = i+1; j < arity_; j++) {
    // for (vector< pair<EnumeratedVariable* ,EnumeratedVariable* > >::iterator it = listxy.begin(); it != listxy.end(); ++it) {
    // for (vector< pair<EnumeratedVariable* ,EnumeratedVariable* > >::iterator it = rndlistxy.begin(); it != rndlistxy.end(); ++it) {

    TSCOPE scopeinv;
    getScope(scopeinv);
    for(TSCOPE::iterator it1 = scopeinv.begin(); it1 != scopeinv.end(); ++it1) {
        TSCOPE::iterator it2 = it1;
        for(++it2; it2 != scopeinv.end(); ++it2) {
            // int i = (*it).first;
            // int j = (*it).second;
            // EnumeratedVariable* x = (*it).first; //scope[i];
            // EnumeratedVariable* y = (*it).second; //scope[j];
            EnumeratedVariable* x = (EnumeratedVariable *) wcsp->getVar((*it1).first);
            EnumeratedVariable* y = (EnumeratedVariable *) wcsp->getVar((*it2).first);
            //	   cout << "try " << x->getName() << "," << y->getName() << ".." << endl;

            TUPLES fproj;
            projectxy(x,y,fproj);

            String t;
            vector<Cost> xy;
            unsigned int a,b;
            unsigned int sizex = x->getDomainInitSize();
            unsigned int sizey = y->getDomainInitSize();

            for (a = 0; a < sizex; a++)
                for (b = 0; b < sizey; b++)
                    xy.push_back(default_cost);

            TUPLES::iterator it =  fproj.begin();
            while(it != fproj.end()) {
                t = it->first;
                a = t[0] - CHAR_FIRST;
                b = t[1] - CHAR_FIRST;
                xy[ a * sizey + b ]	= it->second;
                it++;
            }
            if(fproj.size() > 0 || default_cost > MIN_COST){
                //int res =
                wcsp->postBinaryConstraint(x->wcspIndex, y->wcspIndex, xy);
                //		 if (!wcsp->getCtr(res)->universal()) cout << ".. succeed" << endl;
            }
            if (deconnected()) return;
        }
    }
}

double NaryConstraintMap::computeTightness()
{
    int count = 0;
    double sum = 0;
    Cost costs[pf->size()];
    TUPLES::iterator  it = pf->begin();
    while(it != pf->end()) {
        Cost c =  it->second;
        sum += to_double(min(wcsp->getUb(), c));
        costs[count] = min(wcsp->getUb(), c);
        count++;
        it++;
    }
    Long psize = getDomainSizeProduct();
    if (psize >= LONGLONG_MAX) {
        tight = to_double(min(wcsp->getUb(), default_cost));
    } else {
        if (ToulBar2::weightedTightness == 2) {
            tight = to_double(stochastic_selection<Cost>(costs, 0, count-1, count / 2)); // TO BE IMPPROVED!!!
        } else {
            if (psize > count) {
                tight =  ((double) default_cost * (psize - count) + sum) / (double) psize;
            } else {
                tight =  sum / (double) count;
            }
        }
    }
    return tight;
}

void NaryConstraintMap::print(ostream& os)
{
    TUPLES& f = *pf;
    os << endl << this << " f(";

    int unassigned_ = 0;
    Long totaltuples = 1;
    for(int i = 0; i < arity_;i++) {
        if(scope[i]->unassigned()) unassigned_++;
        os << scope[i]->wcspIndex;
        if(i < arity_-1) os << ",";
        totaltuples = totaltuples * scope[i]->getDomainInitSize();
    }
    os << ")    ";
    if (ToulBar2::weightedDegree) {
        os << "/" << getConflictWeight();
        for(int i = 0; i < arity_;i++) {
            os << "," << conflictWeights[i];
        }
    }
    os << " |f| = " << f.size() << " / " << totaltuples;
    os << "   default_cost: " << default_cost;
    os << "   arity: " << arity();
    os << "   unassigned: " << (int) nonassigned << "/" << unassigned_ << "         ";

    //	assert(nonassigned == unassigned_); // not valid when used with assignLS

    /*TSCOPE::iterator it = scope_inv.begin();
	while(it != scope_inv.end()) {
		os << "(" << it->first << ",idx: " << it->second << ") ";
		++it;
	}*/
    os << endl;


    if (ToulBar2::verbose >= 4) {
        os << "tuples: {";
        TUPLES::iterator  it = f.begin();
        while(it != f.end()) {
            String t = it->first;
            Cost c =  it->second;
            it++;
            os << "<";
            for(unsigned int i=0;i<t.size();i++) {
                os << t[i] - CHAR_FIRST;
                if (i<t.size()-1) os << " ";
            }
            os << "," << c << ">";
            if(it != f.end()) os << " ";
        }
        os << "} " << endl;
    }
}

void NaryConstraintMap::dump(ostream& os, bool original)
{
    if (original) {
        TUPLES& f = *pf;
        os << arity_;
        for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
        os << " " << default_cost << " " << f.size() << endl;

        TUPLES::iterator  it = f.begin();
        while(it != f.end()) {
            String t = it->first;
            Cost c =  it->second;
            it++;
            for(unsigned int i=0;i<t.size();i++) {
                os << t[i] - CHAR_FIRST << " ";
            }
            os << c << endl;
        }
    } else {
        os << nonassigned;
        for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
        String tuple;
        Cost cost;
        Long nbtuples = 0;
        first();
        while (next(tuple,cost)) {
            nbtuples++;
        }
        os << " " << default_cost << " " << nbtuples << endl;
        first();
        while (next(tuple,cost)) {
            for(unsigned int i=0;i<tuple.size();i++) {
                if (scope[i]->unassigned()) os << scope[i]->toCurrentIndex(scope[i]->toValue(tuple[i] - CHAR_FIRST)) << " ";
            }
            os << ((original)?cost:min(wcsp->getUb(),cost)) << endl;
        }
    }
}


TrieNode::TrieNode() {
    ptrs = NULL;
    letters = NULL;
    word = NULL;
    c = MIN_COST;
}


void TrieNode::iniLeaf(Char *suffix) {
    leaf = true;
    word = new Char[Strlen(suffix)+1];
    if (word == 0) exit(-1);
    Strcpy(word,suffix);
}

void TrieNode::iniNonLeaf(Char ch) {
    ptrs = new TrieNode*;
    letters = new Char[2];
    if (ptrs == 0 || letters == 0) exit(1);
    leaf = false;
    endOfWord = false;
    *ptrs = 0;
    *letters = ch;
    *(letters+1) = '\0';
}




Trie::Trie(Char* word, Cost c) : notFound(-1) {
    root = new TrieNode();
    root->iniNonLeaf(*word);
    TrieNode* lf = createLeaf(*word,word+1,root); // to avoid later tests;
    lf->c = c;
}

void Trie::printTrie(int depth, TrieNode *p, Char *prefix) {
    register int i;             // assumption: the root is not a leaf
    if (p->leaf) {              // and it is not null;
        TrieNode *lf = (TrieNode*) p;
        for (i = 1; i <= depth; i++) cout << "   ";
        cout << " >>" << prefix << "|" << lf->word << "   cost: " << lf->c << endl;

    }
    else {
        for (i = Strlen(p->letters)-1; i >= 0; i--)
            if (p->ptrs[i] != 0) {             // add the letter
                prefix[depth] = p->letters[i]; // corresponding to
                prefix[depth+1] = '\0';        // position i to prefix;
                printTrie(depth+1,p->ptrs[i],prefix);
            }
        if (p->endOfWord) {
            prefix[depth] = '\0';
            for (i = 1; i <= depth+1; i++) cout << "   ";
            cout << ">>>" << prefix << "\n";
        }
    }
}

int Trie::position(TrieNode *p, Char ch) {
    unsigned int i;
    for (i = 0; i < Strlen(p->letters) && p->letters[i] != ch; i++);
    if (i < Strlen(p->letters)) return i;
    else return notFound;
}

TrieNode* Trie::find(const Char *word) {
    TrieNode *p = root;
    int pos;
    while (true)
        if (p->leaf) {                      			   // node p is a leaf
            if (Strcmp(word,p->word) == 0)  return p;     // suffix of word should be found;
            else return NULL;
        }
        else if (*word == '\0')              			   // the end of word has
            if (p->endOfWord) return p;    			   // the endOfWord marker
            else return NULL;              			   // in node p set to true;
        else if ((pos = position(p,*word)) != notFound &&
                p->ptrs[pos] != 0) {        			   // continue
            p = p->ptrs[pos];               			   // path, if possible,
            word++;
        }
        else return NULL;                   			   // otherwise failure;
}

void Trie::addCell(Char ch, TrieNode *p, int stop) {
    int i, len = Strlen(p->letters);
    Char *s = p->letters;
    TrieNode **tmp = p->ptrs;
    p->letters = new Char[len+2];
    p->ptrs    = new TrieNode*[len+1];
    if (p->letters == 0 || p->ptrs == 0) { exit(1); }
    for (i = 0; i < len+1; i++) p->ptrs[i] = 0;
    if (stop < len)                // if ch does not follow all letters in p,
        for (i = len; i >= stop+1; i--) { // copy from tmp letters > ch;
            p->ptrs[i]    = tmp[i-1];
            p->letters[i] = s[i-1];
        }
    p->letters[stop] = ch;
    for (i = stop-1; i >= 0; i--) {           // and letters < ch;
        p->ptrs[i]    = tmp[i];
        p->letters[i] = s[i];
    }
    p->letters[len+1] = '\0';
    delete [] s;
}

TrieNode* Trie::createLeaf(Char ch, Char *suffix, TrieNode *p) {
    int pos = position(p,ch);
    if (pos == notFound) {
        for (pos = 0; (pos < (int)Strlen(p->letters)) && (p->letters[pos] < ch); pos++);
        addCell(ch,p,pos);
    }
    TrieNode* tn = new TrieNode();
    p->ptrs[pos] = (TrieNode*) tn;
    tn->iniLeaf(suffix);
    return tn;
}

void Trie::insert(Char *word, Cost c) {
    TrieNode *p = root;
    TrieNode *lf, *newlf;
    int offset;
    unsigned int pos;
    Char *hold = word;
    while (true) {
        if (*word == '\0') {            								// if the end of word reached,
            if (p->endOfWord) cout << "Duplicate entry1 " << hold << endl;
            else p->endOfWord = true;  								// set endOfWord to true;
            return;
        }                               								// if position in p indicated
        pos = position(p,*word);
        if ((int)pos == notFound) {          								 // by the first letter of word
            newlf = createLeaf(*word,word+1,p);								 // does not exist, create
            newlf->c = c;
            return;                    								 // a leaf and store in it the
        }                               								 // unprocessed suffix of word;
        else if ((int)pos != notFound && p->ptrs[pos]->leaf) {
            lf = (TrieNode*) p->ptrs[pos];      						 // hold this leaf;
            if (Strcmp(lf->word,word+1) == 0) {
                cout << "Duplicate entry2 " << hold << endl;
                return;
            }
            offset = 0;
            do {
                pos = position(p,word[offset]);						 // word == "ABC", leaf = "ABCDEF" => leaf = "DEF";
                if ((int)Strlen(word) == offset+1) {
                    p->ptrs[pos] = new TrieNode();
                    p->ptrs[pos]->iniNonLeaf(word[offset]);
                    p->ptrs[pos]->endOfWord = true;
                    newlf = createLeaf(lf->word[offset],lf->word + offset+1,p->ptrs[pos]);
                    newlf->c = c;
                    return;
                }                 															// word == "ABCDE", leaf = "ABC" => leaf = "DEF";
                else if ((int)Strlen(lf->word) == offset) {
                    p->ptrs[pos] = new TrieNode();
                    p->ptrs[pos]->iniNonLeaf(word[offset+1]);
                    p->ptrs[pos]->endOfWord = true;
                    newlf = createLeaf(word[offset+1],word+offset+2,p->ptrs[pos]);
                    newlf->c = c;
                    return;
                }
                p->ptrs[pos] = new TrieNode();
                p->ptrs[pos]->iniNonLeaf(word[offset+1]);
                p = p->ptrs[pos];
                offset++;
            } while (word[offset] == lf->word[offset-1]);
            offset--;
            // word = "ABCDEF", leaf = "ABCPQR" => leaf('D') = "EF", leaf('P') = "QR";
            // check whether there is a suffix left:
            // word = "ABCD", leaf = "ABCPQR" => leaf('D') = null, leaf('P') = "QR";
            Char emptystr[1];
            emptystr[0] = '\0';
            Char * s = &emptystr[0];
            if ((int)Strlen(word) > offset+2) s = word+offset+2;
            newlf = createLeaf(word[offset+1],s,p);
            newlf->c = c;
            // check whether there is a suffix left:
            // word = "ABCDEF", leaf = "ABCP" => leaf('D') = "EF", leaf('P') = null;
            if ((int)Strlen(lf->word) > offset+1) s = lf->word+offset+1; else s = &emptystr[0];
            newlf = createLeaf(lf->word[offset],s,p);
            newlf->c = lf->c;

            delete [] lf->word;
            delete lf;
            return;
        }
        else {
            p = p->ptrs[pos];
            word++;
        }
    }
}

void Trie::printTrie()
{
    *prefix = '\0';
    printTrie(0,root,prefix);
    cout << "__________________________" << endl;
}




NaryConstrie::NaryConstrie(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval)
: NaryConstraint(wcsp, scope_in, arity_in, defval)
{
    f = NULL;
}


NaryConstrie::NaryConstrie(WCSP *wcsp)
: NaryConstraint(wcsp)
{
    f = NULL;
}



NaryConstrie::~NaryConstrie()
{
    if(f) delete f;
}

// for adding a tuple in f
// scope_in contains the order of the values in String tin
void NaryConstrie::setTuple( String& tin, Cost c, EnumeratedVariable** scope_in )
{
    String t(tin);
    if(scope_in) {  for(int i = 0; i < arity_; i++) t[getIndex(scope_in[i])] = tin[i];  }
    Char tch[80];
    Strcpy(tch,t.c_str());
    if(!f) { f = new Trie( tch, c ); }
    else f->insert(tch, c);
}

void NaryConstrie::addtoTuple( String& tin, Cost c, EnumeratedVariable** scope_in ) {
    String t(tin);
    Cost csum = eval(t) + c;
    if(scope_in) {  for(int i = 0; i < arity_; i++) t[getIndex(scope_in[i])] = tin[i];  }
    Char tch[80];
    Strcpy(tch,t.c_str());
    if(!f) { f = new Trie( tch, csum ); }
    else f->insert(tch, csum);

}


Cost NaryConstrie::eval( String& s ) {
    Cost c = default_cost;
    if(f) {
        TrieNode* tn = f->find(s.c_str());
        if(tn) c = tn->c;
    }
    return c;
}



/*void NaryConstrie::assign(int varIndex) {
	int i;
    if (connected(varIndex)) {
        deconnect(varIndex);
	    nonassigned = nonassigned - 1;

	   if(nonassigned == 0) {
			Char* t = new Char [arity_ + 1];
			for(i = 0; i < arity_;i++) t[i] = CHAR_FIRST + scope[i]->toIndex(scope[i]->getValue());
			t[i] = '\0';
			deconnect();
	   	    projectLB(eval(String(t)));
			delete [] t;
	   }
    }
}*/


void NaryConstrie::print(ostream& os) {
    int unassigned_ = 0;
    Long totaltuples = 1;
    os << endl << this << " f(";
    for(int i = 0; i < arity_;i++) {
        if(scope[i]->unassigned()) unassigned_++;
        os << scope[i]->wcspIndex;
        if(i < arity_-1) os << ",";
        totaltuples = totaltuples * scope[i]->getDomainInitSize();
    }
    os << ")    ";
    os << "   default_cost: " << default_cost;
    os << "   arity: " << arity();
    os << "   unassigned: " << (int) nonassigned << "/" << unassigned_ << "         ";

    if(f) f->printTrie();
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

