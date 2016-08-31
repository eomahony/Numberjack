/*
 * **************** Data-structure to manage Cluster Tree Decomposition *******************
 *
 */

#include "tb2clusters.hpp"
#include "tb2naryconstr.hpp"
#include "tb2pedigree.hpp"
#include "tb2haplotype.hpp"

#include <list>
#include <algorithm>

/*
 * Comparison between cluster sons
 *
 */

int Cluster::clusterCounter = 0;

bool CmpClusterStructBasic::operator() (const Cluster *lhs, const Cluster *rhs) const
{
    return lhs && rhs && (lhs->getIndex() < rhs->getIndex());
}
bool CmpClusterStruct::operator() (const Cluster *lhs, const Cluster *rhs) const
{
    return lhs && rhs && (lhs->sepSize() < rhs->sepSize() || (lhs->sepSize() == rhs->sepSize() && (lhs->getNbVarsTree() < rhs->getNbVarsTree() || (lhs->getNbVarsTree() == rhs->getNbVarsTree() && lhs->getIndex() < rhs->getIndex()))));
}

/*
 * Separator class derived from NaryConstraint
 *
 */

Separator::Separator(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in) : AbstractNaryConstraint(wcsp, scope_in, arity_in),
        cluster(NULL),
        nonassigned(arity_in),
        isUsed(false),
        lbPrevious(MIN_COST),
        optPrevious(false)
{
    Char* tbuf = new Char[arity_in+1];
    tbuf[arity_in] = '\0';
    for(int i=0;i<arity_in;i++) {
        tbuf[i] = CHAR_FIRST;
        unsigned int domsize = scope_in[i]->getDomainInitSize();
        vars.insert(  scope_in[i]->wcspIndex );
        if(domsize + CHAR_FIRST > (unsigned int)std::numeric_limits<Char>::max()) {
            cerr << "Nary constraints overflow. Try undefine NARYCHAR in makefile." << endl;
            exit(EXIT_FAILURE);
        }
    }
    t = String(tbuf);
    delete [] tbuf;

    linkSep.content = this;

    // initial "delayed" propagation
    if (arity_ == 0) {
        queueSep();
    } else {
        for(int i=0;i<arity_;i++) {
            if (getVar(i)->assigned()) assign(i);
        }
    }
}

Separator::Separator(WCSP *wcsp) : AbstractNaryConstraint(wcsp),
        nonassigned(0),
        isUsed(false),
        lbPrevious(MIN_COST),
        optPrevious(false)
{
}

void Separator::setup(Cluster* cluster_in)
{
    cluster = cluster_in;
    AbstractNaryConstraint::cluster = cluster_in->getParent()->getId();
    delta.clear();
    TVars::iterator it = vars.begin();
    while(it != vars.end()) {
        EnumeratedVariable* var = (EnumeratedVariable*) cluster->getWCSP()->getVar(*it);
        delta.push_back( vector<StoreCost>(var->getDomainInitSize(), StoreCost(MIN_COST)) );
        ++it;
    }

    int nvars = cluster->getNbVars();
    if(!nvars) return;

    Char* sbuf = new Char [cluster->getNbVars()+1];
    int i = 0;
    int nproper = 0;
    it = cluster->beginVars();
    while(it != cluster->endVars()) {
        if (!cluster->isSepVar(*it)) nproper++;
        sbuf[i] = CHAR_FIRST;
        ++it;
        i++;
    }
    sbuf[nproper] = '\0';
    s = String(sbuf);
    delete [] sbuf;
}

void Separator::assign(int varIndex)
{
    if (connected(varIndex)) {
        deconnect(varIndex);
        nonassigned = nonassigned - 1;
        assert(nonassigned >= 0);
        if(nonassigned == 0) {
            assert(!cluster || cluster->isActive());
            queueSep();
        }
    }
}

void Separator::propagate()
{
    if (ToulBar2::verbose >= 3) cout << this << " propagate C" << cluster->getId() << " " << nonassigned << " " << cluster->getParent()->getId() << " " << connected() << endl;
    for(int i=0;connected() && i<arity_;i++) {
        if (getVar(i)->assigned()) assign(i);
    }
    if(nonassigned == 0 && wcsp->getTreeDec()->isInCurrentClusterSubTree(cluster->getParent()->getId())) {
        wcsp->revise(this);
        if(ToulBar2::allSolutions){
            Cost res = MIN_COST;
            BigInteger nb = 0.;
            getSg(res,nb);
            if (nb == 0.)  {
                if (ToulBar2::verbose >= 1) cout << "use #good " << this << endl;
                THROWCONTRADICTION;
            }
            else
                if(nb>0.)
                    unqueueSep();
        }
        else{
            Cost clb = MIN_COST;
            Cost cub = MAX_COST;
            get(clb,cub,&cluster->open);
            bool opt = (clb == cub);
            if (cluster->isActive()) {
                Cost lbpropa = cluster->getLbRec();
                Cost lb = clb - lbpropa;
                if(opt || lb>MIN_COST) {
                    if (ToulBar2::verbose >= 1) cout << "nogood C" << cluster->getId() << " used in advance (lbpropa=" << lbpropa << " ,lb+=" << lb << ")" << endl;
                    assert(lb >= MIN_COST);
                    if (opt) unqueueSep();
                    // atomic operations:
                    isUsed = true;
                    cluster->deactivate();
                    assert(cluster->getParent()->getId() == Constraint::cluster);
                    cluster->getParent()->increaseLb(lbpropa);
                    if (lb>MIN_COST) projectLB(lb); // project into global lb and into parent cluster
                    lbPrevious = clb;
                    optPrevious = opt;
                    // end of atomic operations.
                }
            } else if (isUsed && cluster->getParent()->isActive()) {
                if (clb > lbPrevious || (opt==true && optPrevious==false)) {
                    if (ToulBar2::verbose >= 1) cout << "nogood C" << cluster->getId() << " used in advance (lbPrevious=" << lbPrevious << " ,optPrevious=" << optPrevious << " ,clb=" << clb << " ,opt=" << opt << ")" << endl;
                    if (opt) unqueueSep();
                    // atomic operations:
                    if (clb > lbPrevious) projectLB(clb - lbPrevious); // project into global lb and into parent cluster
                    lbPrevious = clb;
                    optPrevious = opt;
                    // end of atomic operations.
                }
            }
        }
    }
}

void Separator::set( Cost clb, Cost cub, Solver::OpenList **open )
{
    assert(clb <= cub);
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    Cost deltares = MIN_COST;
    if (ToulBar2::verbose >= 1) cout << "( ";
    TVars::iterator it = vars.begin();
    while(it != vars.end()) {
        assert(wcsp->assigned(*it));
        Value val = wcsp->getValue(*it);
        if (ToulBar2::verbose >= 1) cout << "(" << *it << "," << val << ") ";
        t[i] = val + CHAR_FIRST;
        deltares += delta[i][val];
        ++it;
        i++;
    }
    if (ToulBar2::verbose >= 1) cout << ")";
    assert(clb < cub || clb + deltares >= MIN_COST);
    //assert(nogoods.find(String(t)) == nogoods.end() || nogoods[String(t)].second <= MAX(MIN_COST, c + deltares));
    TNoGoods::iterator itng = nogoods.find(t);
    if (ToulBar2::verbose >= 3) {
        cout << " <C" << cluster->getId() << ",";
        Cout << t;
        cout << "," << MAX(MIN_COST, clb + deltares) << "," << MAX(MIN_COST, cub + deltares) << ">" << endl;
    }
    if (open) {
        if (*open) {
            // open node list already found => the corresponding nogood has been created before
            assert(itng != nogoods.end());
            assert(*open == &itng->second.third);
            itng->second.first = MAX(itng->second.first, clb + deltares);
            itng->second.second = MIN(itng->second.second, MAX(MIN_COST, cub + ((cub < MAX_COST)?deltares:MIN_COST)));
            if (ToulBar2::verbose >= 1) cout << " Learn nogood " << itng->second.first << ", cub= " <<  itng->second.second << ", delta= " << deltares << " on cluster " << cluster->getId() << endl;
        } else {
            assert(itng == nogoods.end());
            nogoods[t] = make_triplet(MAX(MIN_COST, clb + deltares), MAX(MIN_COST, cub + ((cub < MAX_COST)?deltares:MIN_COST)), Solver::OpenList());
            if (ToulBar2::verbose >= 1) cout << " Learn nogood " << nogoods[t].first << ", cub= " <<  nogoods[t].second << ", delta= " << deltares << " on cluster " << cluster->getId() << endl;
            *open = &nogoods[t].third;
        }
    } else {
        if (itng == nogoods.end()) {
            nogoods[t] = make_triplet(MAX(MIN_COST, clb + deltares), MAX(MIN_COST, cub + ((cub < MAX_COST)?deltares:MIN_COST)), Solver::OpenList());
            if (ToulBar2::verbose >= 1) cout << " Learn nogood " << nogoods[t].first << ", cub= " <<  nogoods[t].second << ", delta= " << deltares << " on cluster " << cluster->getId() << endl;
        } else {
            itng->second.first = MAX(itng->second.first, clb + deltares);
            itng->second.second = MIN(itng->second.second, MAX(MIN_COST, cub + ((cub < MAX_COST)?deltares:MIN_COST)));
            if (ToulBar2::verbose >= 1) cout << " Learn nogood " << itng->second.first << ", cub= " <<  itng->second.second << ", delta= " << deltares << " on cluster " << cluster->getId() << endl;
        }
    }
}


void Separator::setSg( Cost c, BigInteger nb )
{
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    Cost deltares = MIN_COST;
    if (ToulBar2::verbose >= 1)	cout << "( ";
    TVars::iterator it = vars.begin();
    while(it != vars.end()) {
        assert(wcsp->assigned(*it));
        Value val = wcsp->getValue(*it);
        if (ToulBar2::verbose >= 1)	cout << "(" << *it << "," << val << ") ";
        t[i] = val + CHAR_FIRST;
        deltares += delta[i][val];
        ++it;
        i++;
    }
    assert(c + deltares >= MIN_COST);
    if (ToulBar2::verbose >= 1) cout << ") Learn #good with " << nb << " solutions" <<endl;// /" << cluster->getVarsTree().size() << endl;
    sgoods[t] = TPairSG(MAX(MIN_COST, c + deltares), nb);
}

Cost Separator::getCurrentDelta()
{
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    Cost sumdelta = MIN_COST;
    TVars::iterator it = vars.begin();
    while(it != vars.end()) {
        if (wcsp->assigned(*it)) {
            Value val = wcsp->getValue(*it);
            sumdelta += delta[i][val];
        } else {
            EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar(*it);
            if (wcsp->td->isDeltaModified(x->wcspIndex)) {
                Cost del = -MAX_COST;
                for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) {
                    Value val = *itx;
                    // Cost unaryc = x->getCost(val);
                    // Could use delta[i][val]-unaryc for pure RDS with only one separator per variable
                    if(del < delta[i][val]) del = delta[i][val];
                }
                assert(del > -MAX_COST);
                sumdelta += del;
            }
        }
        ++it;
        i++;
    }
    return sumdelta;
}

bool Separator::get( Cost& clb, Cost& cub, Solver::OpenList **open)
{
    int i = 0;
    clb = MIN_COST;
    cub = MIN_COST;

    if (ToulBar2::verbose >= 1) cout << "( ";
    TVars::iterator it = vars.begin();
    while(it != vars.end()) {
        assert(cluster->getWCSP()->assigned(*it));
        Value val = cluster->getWCSP()->getValue(*it);
        if (ToulBar2::verbose >= 1) cout << "(" << *it << "," << val << ") ";
        t[i] = val + CHAR_FIRST;	 // build the tuple
        clb -= delta[i][val];      // delta structure
        cub -= delta[i][val];      // delta structure
        ++it;
        i++;
    }
    TNoGoods::iterator itng = nogoods.find(t);
    if(itng != nogoods.end()) {
        TPairNG &p = itng->second; // it is crucial here to get a reference to the data triplet object instead of a copy, otherwise open node list would be copied
        if (ToulBar2::verbose >= 1) cout << ") Use nogood " << p.first << ", delta=" << clb << " (cub=" << p.second << ") on cluster " << cluster->getId() << " (active=" << cluster->isActive() << ")" << endl;
        assert(p.first < p.second || clb + p.first >= MIN_COST);
        clb += p.first;
        cub += p.second;
        cub = MAX(MIN_COST,cub);
        cluster->setUb(cub);
        if (open) *open = &p.third;
        if (ToulBar2::btdMode >= 2) {
            Cost lbrds = cluster->getLbRDS();
            assert(clb < cub || clb >= lbrds);
            clb = MAX(lbrds, clb);
        } else {
            clb = MAX(MIN_COST,clb);
        }
        return true;
    } else {
        clb = (ToulBar2::btdMode >= 2)?cluster->getLbRDS():MIN_COST;
        cub = MAX_COST;
        cluster->setUb(MAX_COST);
        if (open) *open = NULL;
        if (ToulBar2::verbose >= 1) cout << ") NOT FOUND for cluster " <<  cluster->getId() << endl;
        return false;
    }
}

BigInteger Separator::getSg( Cost& res, BigInteger& nb )
{
    int i = 0;
    res = MIN_COST;
    if (ToulBar2::verbose >= 1) cout << "( ";
    TVars::iterator it = vars.begin();
    while(it != vars.end()) {
        assert(cluster->getWCSP()->assigned(*it));
        Value val = cluster->getWCSP()->getValue(*it);
        if (ToulBar2::verbose >= 1) cout << "(" << *it << "," << val << ") ";
        t[i] = val + CHAR_FIRST;	 // build the tuple
        res -= delta[i][val];      // delta structure
        ++it;
        i++;
    }
    TSGoods::iterator itsg = sgoods.find(t);
    if(itsg != sgoods.end()) {
        TPairSG p = itsg->second;
        if (ToulBar2::verbose >= 1)	cout << ") Use #good  with nb = " << p.second << "solutions on cluster " << cluster->getId() << endl;
        /*		assert(res + p.first >= MIN_COST);
		res += p.first;*/
        nb = p.second;
        /*res = MAX(MIN_COST,res);*/
        return nb;
    } else {
        /*res = MIN_COST;*/
        if (ToulBar2::verbose >= 1)	cout << ") NOT FOUND for cluster " <<  cluster->getId() << endl;
        return nb=-1;
    }
}

bool Separator::solGet(TAssign& a, String& sol)
{
    int i = 0;
    TVars::iterator it = vars.begin();
    while(it != vars.end()) {
        Value val = a[*it];
        t[i] = val + CHAR_FIRST;	 // build the tuple
        ++it;
        i++;
    }
    TPairSol p;
    TSols::iterator itsol = solutions.find(t);
    if(itsol != solutions.end()) {
        p = itsol->second;
        sol = p.second;

        if (ToulBar2::verbose >= 1) {
            cout << "asking  solution  sep:";
            Cout << t;
            cout << "  cost: " << p.first << endl;
            Cout << "  sol: " << sol << endl;
        }

        return true;
    }
    return false;
}

void Separator::solRec(Cost ub)
{
    WCSP* wcsp = cluster->getWCSP();

    Cost deltares = MIN_COST;
    int i = 0;
    TVars::iterator it = vars.begin();
    while(it != vars.end()) {
        assert(wcsp->assigned(*it));
        Value val = wcsp->getValue(*it);
        t[i] = val + CHAR_FIRST;	 // build the tuple
        deltares += delta[i][val];
        ++it;
        i++;
    }

    //  	TPairSol p;
    //  	TSols::iterator itsol = solutions.find(t);
    //  	if(itsol != solutions.end()) {
    //  		p = itsol->second;
    //  	    assert(p.first > ub + deltares);
    //  	}

    wcsp->restoreSolution(cluster);

    i = 0;
    it = cluster->beginVars();
    while(it != cluster->endVars()) {
        assert(wcsp->assigned(*it));
        if (!cluster->isSepVar(*it)) {
            Value val = wcsp->getValue(*it);
            s[i] = val + CHAR_FIRST;
            i++;
        }
        ++it;
    }

    solutions[t] = TPairSol(ub + deltares,s);

    if (ToulBar2::verbose >= 1) {
        cout << "recording solution  " << " cost: " << ub << " + delta: " << deltares;
        Cout << " sol: " << s <<  " sep: " << t << endl;
    }
}

void Separator::resetLb()
{
    TNoGoods::iterator it = nogoods.begin();
    while(it != nogoods.end()) {
        (it->second).first = MIN_COST;
        (it->second).third = Solver::OpenList();
        ++it;
    }
}

void Separator::resetUb()
{
    TNoGoods::iterator it = nogoods.begin();
    while(it != nogoods.end()) {
        (it->second).second = MAX_COST;
        (it->second).third = Solver::OpenList();
        ++it;
    }
}

void Separator::print(ostream& os)
{
    os << this << " nogoods(";
    Double totaltuples = 1;
    for(int i = 0; i < arity_;i++) {
        os << scope[i]->getName();
        if(i < arity_-1) os << ",";
        totaltuples = totaltuples * scope[i]->getDomainInitSize();
    }
    os << ")    ";
    os << " |nogoods| = " << nogoods.size() << " / " << totaltuples << " min:" << min_element(nogoods.begin(), nogoods.end(), nogoods.value_comp())->second.first << " (" << cluster->getNbBacktracksClusterTree() << " bt)";
    if (ToulBar2::verbose >= 4) {
        os << "nogoods: {";
        TNoGoods::iterator  it = nogoods.begin();
        while(it != nogoods.end()) {
            TPairNG p = it->second;
            os << "<";
            for(unsigned int i=0;i<it->first.size();i++) {
                os << it->first[i] - CHAR_FIRST;
                if (i<it->first.size()-1) os << " ";
            }
            os << "," << p.first << ">";
            if(it != nogoods.end()) os << " ";
            it++;
        }
        os << "} " << endl;
    }
    os << endl;
}


/*
 * Cluster class
 *
 */

Cluster::Cluster(TreeDecomposition *tdin) : td(tdin), wcsp(tdin->getWCSP()), id(-1), parent(NULL), sep(NULL),
        lb(MIN_COST), ub(MAX_COST), lbRDS(MIN_COST),
        active(true),
        countElimVars(1),
        cp(NULL), open(NULL), hbfsGlobalLimit(LONGLONG_MAX), hbfsLimit(LONGLONG_MAX), nbBacktracks(0) {
    instance = clusterCounter++;
}

Cluster::~Cluster() {
    delete cp;
}

void Cluster::addVar( Variable* x ) { vars.insert(x->wcspIndex); }
void Cluster::removeVar( Variable* x ) { vars.erase(x->wcspIndex); }

void Cluster::addVars( TVars& morevars ) {
    set_union( vars.begin(), vars.end(),
            morevars.begin(), morevars.end(),
            inserter(vars, vars.begin()) );
}

void Cluster::addCtr( Constraint* c ) { ctrs.insert(c); }

void Cluster::addEdge( Cluster* c ) { edges.insert(c); }

TClusters::iterator Cluster::removeEdge( TClusters::iterator it ) {
    TClusters::iterator itaux = it;
    ++it;
    edges.erase(itaux);
    return it;
}

void Cluster::removeEdge( Cluster* c )  {
    TClusters::iterator it = edges.find(c);
    if(it != edges.end()) edges.erase(it);
}

void Cluster::addEdges( TClusters& cls )
{
    set_union( edges.begin(), edges.end(),
            cls.begin(), cls.end(),
            inserter(edges, edges.begin()) );
}

void Cluster::addCtrs( TCtrs& ctrsin )
{
    set_union( ctrs.begin(), ctrs.end(),
            ctrsin.begin(), ctrsin.end(),
            inserter(ctrs, ctrs.begin()) );
}

void Cluster::sum( TCtrs& c1, TCtrs& c2, TCtrs& ctout )
{
    set_union( c1.begin(), c1.end(),
            c2.begin(), c2.end(),
            inserter(ctout, ctout.begin()) );
}


TCtrs Cluster::getCtrsTree()
{
    TCtrs ctrsTree;
    for(TClusters::iterator it=descendants.begin(); it!=descendants.end(); ++it)
    {
        Cluster * c = *it;
        c->getCtrs();
        sum(ctrsTree,c->getCtrs(),ctrsTree);

    }
    return ctrsTree;
}

void Cluster::deconnectDiff(TCtrs listCtrsTot,TCtrs listCtrs)
{

    TCtrs listDiff;
    set_difference(listCtrsTot.begin(),listCtrsTot.end(),listCtrs.begin(),listCtrs.end(),inserter(listDiff,listDiff.begin()));
    for(TCtrs::iterator itctr=listDiff.begin(); itctr!=listDiff.end(); ++itctr)
    {
        Constraint * ctr= *itctr;
        ctr->deconnect();
    }
}


void Cluster::deconnectSep()
{
    if(!sep) return;
    TVars::iterator its = beginSep();
    while(its != endSep()) {
        Variable *x = wcsp->getVar(*its);
        ConstraintList* xctrs = x->getConstrs();
        for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
            Cluster* ctrc = td->getCluster(ctr->getCluster()); // warning! return parent cluster if sep
            if (!(ctr->isSep() && isDescendant(ctrc))) {
                // keep descendant separators connected
                //			    if (ctr->isSep()) cout << "deconnect separator parent " << ctr->cluster << " " << *ctr << endl;
                ctr->deconnect();
            }
        }
        x->assign( x->getSupport() );
        ++its;
    }
}

void Cluster::resetLbRec()
{
    if (sepSize() > 0) sep->resetLb();
    if (this != td->getRoot()) open = NULL;
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
        (*iter)->resetLbRec();
    }
}

void Cluster::resetUbRec(Cluster *root)
{
    if(!sep || sepSize()==0 || !root->sep || root->sepSize()==0) return;
    TVars inter;
    td->intersection(sep->getVars(), root->sep->getVars(), inter);
    if (inter.size() > 0) sep->resetUb();
    if (this != td->getRoot()) open = NULL;
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
        (*iter)->resetUbRec(root);
    }
}

Cost Cluster::getLbRec() const
{
    assert(isActive());
    Cost res = lb;
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
        if((*iter)->isActive()) res += (*iter)->getLbRec();
    }
    return res;
}

Cost Cluster::getLbRecRDS()
{
    assert(isActive());
    Cost res = lb;
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
        if((*iter)->isActive()) {
            Cost propa = (*iter)->getLbRecRDS();
            Cost rds = (*iter)->getLbRDS();
            res += MAX(propa,rds);
        }
    }
    return res;
}

void Cluster::reactivate()
{
    assert(!isActive());
    active = true;
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
        assert(!(*iter)->isActive());
        if (!(*iter)->sep->used()) (*iter)->reactivate();
    }
}

void Cluster::deactivate()
{
    if (isActive()) {
        if (ToulBar2::verbose >= 1) cout << "deactive cluster " << getId() << endl;
        active = false;
        for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
            (*iter)->deactivate();
        }
    }
}

void Cluster::setWCSP2Cluster()
{
    for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
        if(!isVar(i)) {
            EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar(i);
            for (ConstraintList::iterator it=x->getConstrs()->begin(); it != x->getConstrs()->end(); ++it) {
                Constraint* ctr = (*it).constr;
                ctr->deconnect();
            }
        }
    }
    for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
        if(!isVar(i)) {
            EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar(i);
            x->assign( x->getSupport() );
        }
    }
}

void Cluster::getElimVarOrder(vector<int> &elimVarOrder)
{
    for (TClusters::reverse_iterator iter = edges.rbegin(); iter != edges.rend(); ++iter) {
        Cluster* cluster = *iter;
        cluster->getElimVarOrder(elimVarOrder);
    }
    for (TVars::reverse_iterator itp = vars.rbegin(); itp != vars.rend(); ++itp) {
        if (!isSepVar(*itp)) {
            elimVarOrder.push_back(*itp);
        }
    }
}

// side-effect: remember last solution
void Cluster::getSolution( TAssign& sol )
{
    TVars::iterator it;
    if(parent == NULL || this == td->getRootRDS()) {
        if(vars.size() == 0) {
        } else {
            it = beginVars();
            while(it != endVars()) {
                assert(wcsp->assigned(*it));
                sol[*it] = wcsp->getValue(*it);
                //				cout << *it << " := " << sol[*it] << endl;
                ++it;
            }
        }
    }
    String s;
    if(sep) {
#ifndef NDEBUG
        bool found = sep->solGet(sol, s);
        assert(found);
#else
        sep->solGet(sol, s);
#endif
        int i = 0;
        it = beginVars();
        while(it != endVars()) {
            if (!isSepVar(*it)) {
                sol[*it] = ((EnumeratedVariable*) wcsp->getVar(*it))->toValue(s[i] - CHAR_FIRST);
                //				cout << *it << " := " << sol[*it] << endl;
                if (!ToulBar2::verifyOpt) wcsp->setBestValue(*it, sol[*it]);
                i++;
            }
            ++it;
        }
    }
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
        Cluster* cluster = *iter;
        cluster->getSolution(sol);
    }
}

bool Cluster::isEdge( Cluster* c )
{
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
        Cluster* cluster = *iter;
        if(c == cluster) return true;
    }
    return false;
}

void Cluster::setup()
{
    if(sep) sep->setup(this);
    if (ToulBar2::hbfs) {
        if (cp) delete cp;
        cp = new Solver::CPStore();
    }
}

void Cluster::accelerateDescendants()
{
    quickdescendants.clear();
    for(int j = 0; j < td->getNbOfClusters(); j++) {
        quickdescendants.push_back( descendants.find( td->getCluster(j) ) != descendants.end() );
    }
}

void Cluster::print()
{
    //cout << "(" << id << ",n:" << getNbVars() << ",lb:" << getLb() << ") ";

    cout << "cluster " << getId();

    cout << " vars {";
    TVars::iterator itp = beginVars();
    while(itp != endVars()) {
        if (!isSepVar(*itp)) {
            cout << wcsp->getVar(*itp)->getName() << ",";
            //cout << *itp << "C" << wcsp->getVar(*itp)->getCluster() << ",";
            assert(wcsp->getVar(*itp)->getCluster()==-1 || wcsp->getVar(*itp)->getCluster() == getId());
        }
        ++itp;
    }
    cout << "\b}";

    if(sep) {
        cout << " U sep {";
        TVars::iterator its = beginSep();
        while(its != endSep()) {
            cout << wcsp->getVar(*its)->getName();
            ++its;
            if(its != endSep()) cout << ",";
        }
        cout << "}";
    }

    if (!edges.empty()) {
        cout << " sons {";
        if (sortedEdges.size() == edges.size()) {
            TClusters::iterator itc = beginSortedEdges();
            while(itc != endSortedEdges()) {
                cout << (*itc)->getId();
                ++itc;
                if(itc != endSortedEdges()) cout << ",";
            }
        } else {
            TClusters::iterator itc = beginEdges();
            while(itc != endEdges()) {
                cout << (*itc)->getId();
                ++itc;
                if(itc != endEdges()) cout << ",";
            }
        }
        cout << "}";
    }

    /*
  	cout << " ctrs {";
  	TCtrs::iterator itctr = beginCtrs();
  	while(itctr != endCtrs()) {
  	  Constraint* ctr = *itctr;
  	  cout << "( ";
  	  for(int i=0;i<ctr->arity();i++) cout << ctr->getVar(i)->wcspIndex << " ";
  	  cout << ">C" << ctr->getCluster() << ")";
  	  ++itctr;
  	}
  	cout << "}";
	cout << " descendants {";
	TClusters::iterator itd = beginDescendants();
	while(itd != endDescendants()) {
      cout << (*itd)->getId();
      ++itd;
      if(itd != endDescendants()) cout << ",";
	}
	cout << "}";
     */

    cout << endl;
}

void Cluster::dump()
{
    //cout << "(" << id << ",n:" << getNbVars() << ",lb:" << getLb() << ") ";

    char clusterVarsFilename[128];
    char sepVarsFilename[128];
    char sonsFilename[128];
    char fatherFilename[128];
    char sepSizeFilename[128];

    sprintf(clusterVarsFilename,"%s.info/%d.vars",getWCSP()->getName().c_str(),getId());
    sprintf(sepVarsFilename,"%s.info/%d.sep",getWCSP()->getName().c_str(),getId());
    sprintf(sonsFilename,"%s.info/%d.sons",getWCSP()->getName().c_str(),getId());
    sprintf(fatherFilename,"%s.info/%d.father",getWCSP()->getName().c_str(),getId());
    sprintf(sepSizeFilename,"%s.info/%d.sepsize",getWCSP()->getName().c_str(),getId());

    ofstream clusterVarsFile(clusterVarsFilename);
    ofstream sepVarsFile(sepVarsFilename);
    ofstream sonsFile(sonsFilename);
    ofstream fatherFile(fatherFilename);
    ofstream sepSizeFile(sepSizeFilename);

    if (parent) { fatherFile << parent->getId(); }
    else { fatherFile << "-1"; }
    fatherFile.close();

    long double separatorSize = 1.0;

    if(sep) {
        TVars::iterator its = beginSep();
        while(its != endSep()) {
            clusterVarsFile << wcsp->getVar(*its)->getName() << " ";
            sepVarsFile << wcsp->getVar(*its)->getName() << " ";
            separatorSize *= (1.0 + wcsp->getVar(*its)->getSup()- wcsp->getVar(*its)->getInf());
            ++its;
        }
    }
    sepSizeFile << separatorSize;

    TVars::iterator itp = beginVars();
    while(itp != endVars()) {
        if (!isSepVar(*itp)) {
            clusterVarsFile << wcsp->getVar(*itp)->getName() << " ";
            assert(wcsp->getVar(*itp)->getCluster()==-1 || wcsp->getVar(*itp)->getCluster() == getId());
        }
        ++itp;
    }

    if (getNbVars() == 0)  clusterVarsFile << " ";

    if (!edges.empty()) {
        TClusters::iterator itc = beginEdges();
        while(itc != endEdges()) {
            sonsFile << (*itc)->getId();
            ++itc;
            if(itc != endEdges()) sonsFile << " ";
        }
    } 

    clusterVarsFile.close();
    sepVarsFile.close();
    sonsFile.close();
    sepSizeFile.close();
}

void Cluster::cartProduct(BigInteger& prodCart)
{
    for(TVars::iterator it = varsTree.begin(); it!= varsTree.end();it++)
    {
        Variable * x = (Variable *) wcsp->getVar(*it);
        prodCart*=x->getDomainSize();
    }
}


/*
 * Tree Decomposition class
 *
 */

TreeDecomposition::TreeDecomposition(WCSP* wcsp_in) :
          wcsp(wcsp_in),
          rootRDS(NULL),
          currentCluster(-1),
          deltaModified(vector<StoreInt>(wcsp_in->numberOfVariables(), StoreInt(false)))
{
}

bool TreeDecomposition::isInCurrentClusterSubTree(int idc)
{
    if (idc<0) return false;
    Cluster* ci = getCurrentCluster();
    Cluster* cj = getCluster(idc);
    assert(ci->isActive());
    return ci->isDescendant(cj);
}

bool TreeDecomposition::isActiveAndInCurrentClusterSubTree(int idc)
{
    if (idc<0) return false;
    Cluster* ci = getCurrentCluster();
    Cluster* cj = getCluster(idc);
    assert(ci->isActive());
    if(!cj->isActive()) return false;
    else return ci->isDescendant(cj);
}

void TreeDecomposition::fusion( Cluster* ci, Cluster* cj )
{
    if(!ci) return;
    if(!cj) return;

    if (ToulBar2::verbose >= 1) cout << "fusion: " << ci->getId() << " " << cj->getId() << endl;

    ci->addVars(cj->getVars());
    ci->addCtrs(cj->getCtrs());
    ci->addEdges(cj->getEdges());
    TClusters::iterator itk =  cj->beginEdges();
    while(itk != cj->endEdges()) {
        Cluster* ck = *itk;
        ++itk;
        ck->removeEdge(cj);
        ck->addEdge(ci);
    }
    ci->removeEdge(ci);
    clusters[ cj->getId() ] = NULL;
    if (ToulBar2::verbose >= 1) { cout << "fusion ci " <<  ci->getId() << ",  cj " <<  cj->getId() << endl; ci->print(); }
    delete cj;
}

bool TreeDecomposition::treeFusion( )
{
    bool done = false;
    for(int j= clusters.size()-1; j >= 0; j--) {
        if(!clusters[j]) continue;
        Cluster* cj = clusters[j];
        if (ToulBar2::verbose >= 3) { cout << "fusion testing "; cj->print(); }

        TClusters::iterator it =  cj->beginEdges();
        while(it != cj->endEdges()) {
            Cluster* c = *it;
            assert(c == clusters[c->getId()]);

            if((c->getId() < cj->getId()) &&
                    (included(c->getVars(), cj->getVars()) ||
                            included(cj->getVars(), c->getVars()) )){
                c->addVars(cj->getVars());
                c->addCtrs(cj->getCtrs());
                c->addEdges(cj->getEdges());
                TClusters::iterator itk =  cj->beginEdges();
                while(itk != cj->endEdges()) {
                    Cluster* ck = *itk;
                    ck->removeEdge(cj);
                    ck->addEdge(c);
                    ++itk;
                }
                c->removeEdge(c);
                clusters[ cj->getId() ] = NULL;
                if (ToulBar2::verbose >= 1) { cout << "fusion ci " <<  c->getId() << ",  cj " <<  cj->getId() << endl; c->print(); }
                delete cj;
                //					done = true;
                break;
            }
            ++it;
        }
    }
    return done;
}

void TreeDecomposition::treeFusions()
{
    while(treeFusion());

    TClusters visited;
    //	Cluster* croot = getBiggerCluster(visited);
    //	heuristicFusionRec(croot, croot);

    int treewidth = 0;
    TClusters sclu;
    for(unsigned int i=0; i < clusters.size(); i++) {
        if(clusters[i])	{
            Cluster* c = clusters[i];
            sclu.insert( c );
            if(c->getNbVars() > treewidth) treewidth = c->getNbVars();
        }
    }
    int i = 0;
    clusters.clear();
    TClusters::iterator it = sclu.begin();
    while(it != sclu.end()) {
        Cluster* c = *it;
        c->setId(i++);
        clusters.push_back(*it);
        ++it;
    }
    if(ToulBar2::verbose >= 2) cout << "Tree decomposition width  : " << treewidth - 1 << endl;
}

void TreeDecomposition::pathFusions(vector<int> &order)
{
    vector<Cluster*> rds;
    int size = clusters.size();
    vector<bool> connected;

    // detect singleton variables
    for(int i=0; i < size; i++) {
        bool isconnected = (clusters[i]->getNbVars() > 1);
        for(int j=0; j < i; j++) {
            if(clusters[j]->isVar(order[i])) {
                isconnected = true;
                break;
            }
        }
        connected.push_back(isconnected);
    }

    for(int i=0; i < size; i++) {
        assert(clusters[i] && clusters[i]->isVar(order[i]));
        Cluster* c = new Cluster( this );
        c->addVar(wcsp->getVar( order[i] ));
        if (connected[i]) {
            for(int j=0; j < i; j++) {
                if(clusters[j]->isVar(order[i])) {
                    for(int l=j+1; l < i; l++) {
                        if (connected[l]) rds[l]->addVar(wcsp->getVar( order[j] ));
                    }
                    c->addVar(wcsp->getVar( order[j] ));
                }
            }
            if (rds.size() > 0) {
                int last = rds.size()-1;
                while (last >= 0 && !connected[last]) last--;
                if (last >= 0) {
                    c->addEdge(rds[last]);
                    rds[last]->addEdge(c);
                }
            }
        }
        rds.push_back(c);
    }

    // fusion on included clusters
    for(int i=0; i < size; i++) {
        if (i < size-1 && included(rds[i]->getVars(), rds[i+1]->getVars())) {
            rds[i+1]->removeEdge(rds[i]);
            rds[i]->removeEdge(rds[i+1]);
            if (rds[i]->getEdges().size() > 0) {
                assert(rds[i]->getEdges().size() == 1);
                TClusters::iterator it = rds[i]->beginEdges();
                rds[i+1]->addEdge(*it);
                (*it)->removeEdge(rds[i]);
                (*it)->addEdge(rds[i+1]);
            }
            delete rds[i];
            rds[i] = NULL;
        }
    }
    clusters.clear();
    for(int i=0; i < size; i++) {
        Cluster* c = rds[i];
        if (c) {
            c->setId(clusters.size());
            clusters.push_back(c);
        }
    }
}

// Minimize tree height when the separators are included
bool TreeDecomposition::reduceHeight( Cluster* c, Cluster* cparent )
{
    assert(c != cparent);
    TClusters::iterator itj;
    itj =  c->beginEdges();
    while(itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj;
        if (cj != cparent) {
            TVars cjsep;
            intersection(c->getVars(), cj->getVars(), cjsep);
            if(cparent && included(cjsep, cparent->getVars())) {
                // replacing the cluster higher in the tree
                //				cout << "move " << cj->getId() << " from " << c->getId() << " to " << cparent->getId() << endl;
                c->removeEdge(cj);
                cparent->addEdge(cj);
                cj->removeEdge(c);
                cj->addEdge(cparent);
                return true;
            } else if(!cparent && cjsep.size() == 0) {
                c->removeEdge(cj);
                cj->removeEdge(c);
                return true;
            }
        }
    }
    itj =  c->beginEdges();
    while(itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj;
        if (cj != cparent) {
            if(reduceHeight(cj, c)) return true;
        }
    }
    return false;
}

int TreeDecomposition::getNextUnassignedVar(TVars *vars)
{
    return *(vars->begin());
}

int TreeDecomposition::getVarMinDomainDivMaxWeightedDegree(TVars *vars)
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (TVars::iterator iter = vars->begin(); iter!= vars->end(); ++iter) {
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize-1, domsize/2).cost;
        }
        double heuristic = (double) domsize / (double) (wcsp->getWeightedDegree(*iter) + 1 + unarymediancost);
        if (varIndex < 0 || heuristic < best - epsilon * best
                || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;
}

void TreeDecomposition::splitClusterRec( Cluster* c,  Cluster* father, unsigned int maxsize )
{
    TVars cvars = c->getVars();
    //  cout << c->getId() << " " << cvars.size() << endl;
    TVars csep;
    if (father) {
        intersection(father->getVars(), c->getVars(), csep);
    }
    TVars cproper;
    difference(cvars, csep, cproper);
    if (cproper.size() > maxsize && (!father || c->getEdges().size() != 1)) {
        Cluster* cprev = NULL;
        TClusters cedges = c->getEdges();
        if (father) cedges.erase(father);
        while (cproper.size() > 0) {
            TVars cnewvars;
            int varIndex = ((ToulBar2::Static_variable_ordering)?getNextUnassignedVar(&cproper):getVarMinDomainDivMaxWeightedDegree(&cproper));
            for (unsigned int i = 0; i < maxsize && varIndex >= 0; i++) {
                cnewvars.insert(varIndex);
                cproper.erase(varIndex);
                varIndex = ((ToulBar2::Static_variable_ordering)?getNextUnassignedVar(&cproper):getVarMinDomainDivMaxWeightedDegree(&cproper));
            }
            //	  TVars::iterator it = cproper.begin();
            //	  for (unsigned int i = 0; i < maxsize && it != cproper.end(); i++) {
            //		  cnewvars.insert(*it);
            //		  ++it;
            //	  }
            //	  TVars cpropernew;
            //	  difference(cproper, cnewvars, cpropernew);
            //	  cproper = cpropernew;
            if (!cprev) {
                c->getVars().clear();
                c->addVars(csep);
                c->addVars(cnewvars);
                c->getEdges().clear();
                if (father) c->addEdge(father);
                cprev = c;
            } else {
                Cluster* cnew = new Cluster( this );
                cnew->setId(clusters.size());
                clusters.push_back( cnew );
                cnew->addVars(cprev->getVars());
                cnew->addVars(cnewvars);
                cnew->addEdge(cprev);
                cprev->addEdge(cnew);
                cprev = cnew;
            }
        }
        assert(cprev->getEdges().size() == 1);
        father = *(cprev->beginEdges());
        cprev->addEdges(cedges);
        TClusters::iterator itj =  cprev->beginEdges();
        while(itj != cprev->endEdges()) {
            Cluster* cj = *itj;
            if (cj != father) {
                cj->removeEdge(c);
                cj->addEdge(cprev);
            }
            ++itj;
        }
        c = cprev;
    }
    TClusters::iterator itj =  c->beginEdges();
    while(itj != c->endEdges()) {
        Cluster* cj = *itj;
        if (cj != father) splitClusterRec(cj, c, maxsize);
        ++itj;
    }
}

TVars TreeDecomposition::boostingVarElimRec( Cluster* c,  Cluster* father,  Cluster* grandfather, unsigned int maxsize )
{
    TVars addedVarBySons;
    TClusters::iterator itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj; // warning! must be done before going inside boostingVarElimRec as it can delete current cluster/iterator by the following removeEdge(c) operation
        if (cj != father) {
            TVars cjaddedvars;
            cjaddedvars = boostingVarElimRec(cj, c, father, maxsize);
            sum(addedVarBySons, cjaddedvars, addedVarBySons);
        }
    }
    if (father && c->getEdges().size() == 1) {
        TVars fathersep;
        if (grandfather) {
            intersection(father->getVars(), grandfather->getVars(), fathersep);
        }
        TVars cvars;
        sum(fathersep, addedVarBySons, fathersep);
        difference(c->getVars(), fathersep, cvars);
        if (cvars.size() <= maxsize) {
            //	  	  cout << c->getId() << " which has " << cvars.size() << " vars (except whose from " << ((grandfather)?grandfather->getId():-1) << ") is merged into " << father->getId() << endl;
            TVars csep;
            intersection(c->getVars(), father->getVars(), csep);
            TVars cproper;
            difference(c->getVars(), csep, cproper);
            father->addVars(cproper);
            father->removeEdge(c);
            clusters.back()->setId(c->getId());
            clusters[ c->getId() ] = clusters.back();
            clusters.pop_back();
            sum(addedVarBySons, cproper, addedVarBySons);
            delete c;
        }
    }
    return addedVarBySons;
}

void TreeDecomposition::mergeClusterRec( Cluster* c,  Cluster* father, unsigned int maxsepsize, unsigned int minpropervar )
{
    TClusters::iterator itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj; // warning! must be done before going inside boostingVarElimRec as it can delete current cluster/iterator by the following removeEdge(c) operation
        if (cj != father) {
            mergeClusterRec(cj, c, maxsepsize, minpropervar);
        }
    }
    if (father) {
        TVars csep;
        intersection(c->getVars(), father->getVars(), csep);
        assert(csep.size()>0);
        if ((csep.size() > maxsepsize) || (c->getVars().size() - csep.size() < minpropervar)) {
            father->addVars(c->getVars());
            father->addEdges(c->getEdges());
            TClusters::iterator itk =  c->beginEdges();
            while(itk != c->endEdges()) {
                Cluster* ck = *itk;
                ++itk;
                ck->removeEdge(c);
                ck->addEdge(father);
            }
            father->removeEdge(father);
            father->removeEdge(c);
            clusters.back()->setId(c->getId());
            clusters[ c->getId() ] = clusters.back();
            clusters.pop_back();
            delete c;
        }
    }
}

// Specific code for cluster fusion based on heuristic criteria
void TreeDecomposition::heuristicFusionRec( Cluster* c, Cluster* noc )
{
    TClusters::iterator it =  c->beginEdges();
    while(it != c->endEdges()) {
        Cluster* cj = *it;
        ++it;
        if(cj == c) continue;
        if(cj == noc) continue;
        heuristicFusionRec( cj, c );
    }

    it =  c->beginEdges();
    while(it != c->endEdges()) {
        Cluster* cj = *it;
        ++it;
        if(cj == c) continue;
        if(cj == noc) continue;
        TVars varsum;
        TVars varinter;
        sum(c->getVars(), cj->getVars(), varsum);
        intersection(c->getVars(), cj->getVars(), varinter);

        int dif = 2;
        bool bf1 = (varinter.size() > 2) && (varsum.size() <= (unsigned int) c->getNbVars() + dif);
        bool bf2 = (varinter.size() > 2) && (varsum.size() <= (unsigned int) cj->getNbVars() + dif);
        bool bf3 = (varinter.size() > 100);
        if(bf1 || bf2 || bf3) {
            fusion(c,cj);
        }
    }
}

Cluster* TreeDecomposition::getBiggerCluster( TClusters& visited )
{
    Cluster* cmax = NULL;
    int maxsize = 0;
    for(unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        if(!c) continue;
        if(visited.find(c) == visited.end()) {
            if(c->getNbVars() > maxsize) {
                maxsize = c->getNbVars();
                cmax = c;
                if (ToulBar2::btdMode == 3) break;
            }
        }
    }
    return cmax;
}

int TreeDecomposition::height( Cluster* r, Cluster* father )
{
    int maxh = 0;
    TClusters::iterator it = r->beginEdges();
    while(it != r->endEdges()) {
        Cluster* adjr = *it;
        if(adjr != father) {
            int h = height(adjr,r);
            if(h > maxh) maxh = h;
        }
        ++it;
    }
    TVars rsep;
    intersection(r->getVars(), father->getVars(), rsep);
    return maxh + r->getNbVars() - rsep.size();
}

int TreeDecomposition::height( Cluster* r )
{
    int maxh = 0;
    TClusters::iterator it = r->beginEdges();
    while(it != r->endEdges()) {
        int h = height(*it,r);
        if(h > maxh) maxh = h;
        ++it;
    }
    return maxh + r->getNbVars();
}

void TreeDecomposition::makeDescendants( Cluster* c )
{
    c->getDescendants().insert(c);
    sum(c->getVarsTree(), c->getVars(), c->getVarsTree());
    TClusters::iterator itj = c->beginEdges();
    while(itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj;
        makeDescendants(cj);
        clusterSum(c->getDescendants(), cj->getDescendants(), c->getDescendants());
        sum(c->getVarsTree(), cj->getVarsTree(), c->getVarsTree());
    }
}

void TreeDecomposition::makeRootedRec( Cluster* c,  TClusters& visited )
{
    TClusters::iterator itj =  c->beginEdges();
    while(itj != c->endEdges()) {
        Cluster* cj = *itj;
        cj->removeEdge(c);
        cj->setParent(c);
        visited.insert(cj);

        TVars cjsep;
        intersection(c->getVars(), cj->getVars(), cjsep);

        //------- Add the constraint separator
        int i = 0;
        int arity = cjsep.size();
        EnumeratedVariable** scopeVars = new EnumeratedVariable* [arity];
        TVars::iterator it = cjsep.begin();
        while(it != cjsep.end()) {
            scopeVars[i] = (EnumeratedVariable *) wcsp->getVar(*it);
            ++it;
            i++;
        }
        cj->setSep( new Separator(wcsp,scopeVars,arity) );
        if(ToulBar2::approximateCountingBTD)
            cj->addCtr(cj->getSep());
        delete [] scopeVars;
        //-------

        makeRootedRec( cj, visited );
        ++itj;
    }
}

int TreeDecomposition::makeRooted()
{
    TClusters visited;
    roots.clear();
    Cluster* root;

    bool selected = false;
    while(visited.size() < clusters.size()) {
        if(!selected && ToulBar2::btdRootCluster >= 0 && ToulBar2::btdRootCluster < (int)clusters.size()) {
            root = getCluster(ToulBar2::btdRootCluster);
            selected = true;
        } else root = getBiggerCluster(visited);

        roots.push_back(root);
        visited.insert(root);
        while(reduceHeight(root,NULL));
        if (ToulBar2::splitClusterMaxSize >= 1) splitClusterRec(root, NULL, ToulBar2::splitClusterMaxSize);
        if (ToulBar2::maxSeparatorSize>=0||ToulBar2::minProperVarSize>=2) mergeClusterRec(root, NULL, ToulBar2::maxSeparatorSize, ToulBar2::minProperVarSize);
        if (ToulBar2::boostingBTD && ToulBar2::elimDegree >= 1) boostingVarElimRec(root, NULL, NULL, ToulBar2::elimDegree);
        while(reduceHeight(root,NULL));
        makeRootedRec(root, visited);
        makeDescendants(root);
    }

    // if it is a forest
    if(roots.size() > 1) {
        root = new Cluster( this );
        root->setId(clusters.size());
        clusters.push_back( root );

        for (list<Cluster*>::iterator iter = roots.begin(); iter!= roots.end(); ++iter) {
            Cluster* oneroot = *iter;

            EnumeratedVariable** scopeVars = new EnumeratedVariable* [1];
            oneroot->setSep( new Separator(wcsp,scopeVars,0) );
            if (oneroot->getNbVars() <= 1 && oneroot->getDescendants().size() == 1) {
                oneroot->getSep()->unqueueSep();
            }
            root->addEdge(oneroot);
            oneroot->setParent(root);
            root->getDescendants().insert( root );
            clusterSum(root->getDescendants(), oneroot->getDescendants(), root->getDescendants());
        }
        roots.clear();
        roots.push_back( root );
    }

    for(unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        c->accelerateDescendants();
        if(c->getSep()) c->getSep()->setSep();
        if (!ToulBar2::approximateCountingBTD)
        {
            int posx = 0;
            TVars::iterator itv = c->beginVars();
            while(itv != c->endVars()) {
                Variable* var = wcsp->getVar(*itv);
                if(!c->isSepVar( var->wcspIndex )) var->setCluster(c->getId());
                else {
                    var->setSep();
                    var->addCluster(c->getId(), posx++);  // we add the cluster and also the position of the variable for the delta structure
                }
                ++itv;
            }
        }
        c->setup();
    }
    rootRDS = NULL;
    root->sortEdgesRec();

    int treewidth = 0;
    for(unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        if(c->getNbVars() > treewidth) treewidth = c->getNbVars();
    }
    if(ToulBar2::verbose >= 0) cout << "Tree decomposition width  : " << treewidth - 1 << endl;
    //	if (treewidth > 30) {cout << "Sorry, cannot perform exact haplotype reconstruction! Please, try approximate methods." << endl; exit(EXIT_FAILURE);} // Warning! QTLmap specific exit

    return height(root);
}

void TreeDecomposition::buildFromOrder()
{
    vector<int> order;
    ((WCSP *) wcsp)->elimOrderFile2Vector(ToulBar2::varOrder, order);
    if (!ToulBar2::varOrder) {
        int n = wcsp->numberOfVariables();
        for(int i=n-1;i>=n/2;--i) {
            int tmp  = order[n-i-1];
            order[n-i-1] = order[i];
            order[i] = tmp;
        }
    }

    if(clusters.size() > 0) {
        for(unsigned int i=0;i<clusters.size();i++) {
            Cluster* c = clusters[i];
            if(c) delete c;
        }
    }
    clusters.clear();

    for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
        Cluster* c = new Cluster( this );
        c->setId(i);
        c->addVar( wcsp->getVar( order[i] ) );
        clusters.push_back( c );
    }
    ConstraintSet usedctrs;

    for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
        Variable* x = wcsp->getVar( order[i] );
        Cluster* c  = clusters[i];

        ConstraintList* xctrs = x->getConstrs();
        for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
            bool used = usedctrs.find( ctr ) != usedctrs.end();
            if(!used) {
                usedctrs.insert( ctr );
                c->addCtr(ctr);
                for(int k=0; k < ctr->arity(); k++) if (ctr->getVar(k)->unassigned()) c->addVar( ctr->getVar(k) );
            }
        }

        for(unsigned int j=i+1;j<wcsp->numberOfVariables();j++)
        {
            if(c->isVar(order[j])) {
                Cluster* cj  = clusters[j];
                TVars::iterator it = c->beginVars();
                while(it != c->endVars()) {
                    cj->addVar( wcsp->getVar(*it) );
                    ++it;
                }
                cj->removeVar(x);
                c->addEdge( cj );
                cj->addEdge( c );
                break;
            }
        }
    }
    buildFromOrderNext(order);
}

void TreeDecomposition::buildFromOrderForApprox()
{

    vector<int> order;
    bool firstComponent = true;
    int sizepart = 0;				//number of parts in the built partition
    ConstraintSet totalusedctrs;	// constraints already in a part
    vector<int> degreeinusedctr;	// number of constraints not adding for each variable
    //	int nbcstr = 0;					//
    double time;

    ((WCSP *) wcsp)->elimOrderFile2Vector(ToulBar2::varOrder, order);
    if (!ToulBar2::varOrder) {
        int n = wcsp->numberOfVariables();
        for(int i=n-1;i>=n/2;--i) {
            int tmp  = order[n-i-1];
            order[n-i-1] = order[i];
            order[i] = tmp;
        }
    }

    if(clusters.size() > 0) {
        for(unsigned int i=0;i<clusters.size();i++) {
            Cluster* c = clusters[i];
            if(c) delete c;
        }
    }
    clusters.clear();
    for(unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
        Variable* x = wcsp->getVar( i );
        //		degree.push_back(x->getTrueDegree());
        degreeinusedctr.push_back(x->getDegree());
    }
    time =cpuTime();
    while (totalusedctrs.size() < wcsp->numberOfConnectedConstraints() )//&& nbparties<4)
    {
        ConstraintSet currentusedctrs;						// liste des contraintes contenues dans la partie courante
        TVars currentusedvars;									// liste des variables contenues dans la partie courante
        TVars inusedvars;										// liste des variables qui n'ont pas encore ete etudiee dans la partie courante
        vector<Variable *> currentRevElimOrder;					// liste des variables dans l'ordre inverse construit
        sizepart++;
        for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
            if (wcsp->unassigned(i)) {
                if (wcsp->getDegree(i) == 0) {
                    if (firstComponent) {
                        currentRevElimOrder.push_back(wcsp->getVar(i));
                        currentusedvars.insert(i);
                    }
                }
                else{
                    if(degreeinusedctr[i]>0)
                        inusedvars.insert( i );
                }
            }
        }

        maxchord(sizepart, order, totalusedctrs, inusedvars, currentusedvars, currentRevElimOrder, currentusedctrs);

        // insert into tree decomposition

        // supprime les variables qui n'ont aucune contraintes dans la partition
        for(vector<Variable *>::iterator it=currentRevElimOrder.begin(); it != currentRevElimOrder.end();) {
            if (currentusedvars.find((*it)->wcspIndex) == currentusedvars.end()) {
                it = currentRevElimOrder.erase(it);
            } else ++it;
        }


        if(sizepart == 1) cout << endl;
        cout << "part " << sizepart << " : " << currentRevElimOrder.size() << " variables and " << currentusedctrs.size() << " constraints (really added)\n";
        if(ToulBar2::debug >=1 ||ToulBar2::verbose >=3)//affichage
        {
            cout << "\tVariables : ";
            for(vector<Variable *>::iterator it=currentRevElimOrder.begin(); it != currentRevElimOrder.end();it++) {
                cout << (*it)->wcspIndex << " " ;
            }
            cout << endl;
            cout << "\tContraintes : ";
            for(ConstraintSet::iterator it=currentusedctrs.begin(); it != currentusedctrs.end();it++) {
                cout << "[";
                for(int k=0; k < (*it)->arity(); k++) {
                    cout << (*it)->getVar(k)->wcspIndex;
                    if(k!= (*it)->arity()-1) cout << " ";
                }
                cout << "] ";
            }
            cout << endl;


        }

        insert(sizepart,currentRevElimOrder,currentusedctrs);
        firstComponent = false;
    }
    time = cpuTime() - time;
    cout << "--> number of parts : " << sizepart <<endl;
    cout << "--> time : " << time << " seconds. "<< endl << endl;
    buildFromOrderNext(order);

}

void TreeDecomposition::buildFromOrderNext(vector<int> &order)
{

    if (ToulBar2::verbose >= 2) {
        cout << "----- Before fusions process: " << endl;
        for(unsigned int i=0; i < clusters.size(); i++) {
            if(!clusters[i]) continue;
            Cluster* c = clusters[i];
            c->print();
        }
        cout << "----- fusions process starting... " << endl;
    }

    if (ToulBar2::btdMode == 3) pathFusions(order);
    else treeFusions();

    if (ToulBar2::verbose >= 2) {
        cout << "----- After fusions process: " << endl;
        for(unsigned int i=0; i < clusters.size(); i++) {
            if(!clusters[i]) continue;
            Cluster* c = clusters[i];
            c->print();
        }
        cout << "----- fusions process ended... " << endl;
    }

    for(unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        c->getDescendants().clear();
    }

    int h = makeRooted();
    if(ToulBar2::verbose >= 0) cout << "Tree decomposition height : " << h << endl;
    if (!ToulBar2::approximateCountingBTD)
    {
        // assign constraints to clusters and check for duplicate ternary constraints
        for (unsigned int i=0; i<wcsp->numberOfConstraints(); i++) {
            Constraint* ctr = wcsp->getCtr(i);
            ctr->assignCluster();
            if (ctr->connected() && !ctr->isSep()) {
                if(ctr->arity() == 3) {
                    TernaryConstraint* tctr = (TernaryConstraint*) ctr;
                    tctr->xy->setCluster( tctr->getCluster() );
                    tctr->xz->setCluster( tctr->getCluster() );
                    tctr->yz->setCluster( tctr->getCluster() );
                }
            }
        }
        for (int i=0; i<wcsp->elimBinOrder; i++) if (wcsp->elimBinConstrs[i]->connected()) {
            Constraint* ctr = wcsp->elimBinConstrs[i];
            ctr->assignCluster();
            if (ctr->connected() && !ctr->isSep()) {
                if(ctr->arity() == 3) {
                    TernaryConstraint* tctr = (TernaryConstraint*) ctr;
                    tctr->xy->setCluster( tctr->getCluster() );
                    tctr->xz->setCluster( tctr->getCluster() );
                    tctr->yz->setCluster( tctr->getCluster() );
                }
            }
        }
        for (int i=0; i<wcsp->elimTernOrder; i++) if (wcsp->elimTernConstrs[i]->connected()) {
            Constraint* ctr = wcsp->elimTernConstrs[i];
            ctr->assignCluster();
            if (ctr->connected() && !ctr->isSep()) {
                if(ctr->arity() == 3) {
                    TernaryConstraint* tctr = (TernaryConstraint*) ctr;
                    tctr->xy->setCluster( tctr->getCluster() );
                    tctr->xz->setCluster( tctr->getCluster() );
                    tctr->yz->setCluster( tctr->getCluster() );
                }
            }
        }

        for (unsigned int i=0; i<wcsp->numberOfConstraints(); i++) {
            Constraint* ctr = wcsp->getCtr(i);
            if (ctr->connected() && !ctr->isSep()) {
                if(ctr->arity() == 3) {
                    TernaryConstraint* tctr = (TernaryConstraint*) ctr;
                    tctr->setDuplicates();
                    assert(tctr->xy->getCluster() == tctr->getCluster() &&
                            tctr->xz->getCluster() == tctr->getCluster() &&
                            tctr->yz->getCluster() == tctr->getCluster() );
                }
            }
        }
        for (int i=0; i<wcsp->elimTernOrder; i++) if (wcsp->elimTernConstrs[i]->connected()) {
            Constraint* ctr = wcsp->elimTernConstrs[i];
            if (ctr->connected() && !ctr->isSep()) {
                if(ctr->arity() == 3) {
                    TernaryConstraint* tctr = (TernaryConstraint*) ctr;
                    tctr->setDuplicates();
                    assert(tctr->xy->getCluster() == tctr->getCluster() &&
                            tctr->xz->getCluster() == tctr->getCluster() &&
                            tctr->yz->getCluster() == tctr->getCluster() );
                }
            }
        }
    }
    if(ToulBar2::verbose >= 0) cout << "Number of clusters        : " << clusters.size() << endl;
    if(ToulBar2::debug >= 1 || ToulBar2::verbose >= 1) print();
    if (ToulBar2::dumpWCSP) dump();
    assert(verify());

}


void  TreeDecomposition::maxchord(int sizepart, vector<int> &order, ConstraintSet &totalusedctrs, TVars &inusedvars, TVars &currentusedvars, vector<Variable *> &currentRevElimOrder,ConstraintSet &currentusedctrs)
{
    vector<TVars>  listeVars(wcsp->numberOfVariables());	// liste des voisins d'ordre superieur de chaque variable
    int nbcstr = 0;
    double time, timetot = 0;
    while (inusedvars.size() > 0)
    {
        int maxsize = -1;
        Variable* maxvar = NULL;	/* next variable */

        //Choose the nex variable
        for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
            Variable* x = wcsp->getVar( order[i] );
            if (inusedvars.find( x->wcspIndex ) != inusedvars.end()) {
                int size = listeVars[x->wcspIndex].size();
                if (size > maxsize) {
                    maxsize = size;
                    maxvar = x;
                }
            }
        }

        if (maxvar) {
            //				cout << "Variable choisie: " << maxvar->wcspIndex << " ";
            ConstraintList* xctrs = maxvar->getConstrs();
            for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it)
            {
                Constraint* ctr = (*it).constr;
                bool used = totalusedctrs.find( ctr ) != totalusedctrs.end();
                if( !used )
                {
                    TVars scopectr;
                    TVars sc;
                    for(int k=0; k < ctr->arity(); k++) {
                        Variable *x = ctr->getVar(k);
                        if( x->wcspIndex!=maxvar->wcspIndex && wcsp->unassigned(x->wcspIndex) )
                        {
                            sc.insert(x->wcspIndex);
                            if( inusedvars.find(x->wcspIndex)!=inusedvars.end() )
                                scopectr.insert(x->wcspIndex);
                        }
                    }
                    if( scopectr.size()==0 )
                    { //all edges of the ctr are in the sub graph => the cstr is added in this current part
                        if( included(sc,listeVars[maxvar->wcspIndex]) )
                        {
                            ConstraintSet subctr;
                            nbcstr++;
                            currentusedctrs.insert(ctr);
                            totalusedctrs.insert(ctr);
                            time = cpuTime();
                            subctr = ctr->subConstraint();
                            ctrSum(totalusedctrs, subctr, totalusedctrs);
                            ctrSum(currentusedctrs, subctr, currentusedctrs);
                            time = time -cpuTime();
                            timetot += time;
                            sum(currentusedvars, sc, currentusedvars);
                            currentusedvars.insert(maxvar->wcspIndex);

                        }
                    }

                    for(TVars::iterator i=scopectr.begin(); i != scopectr.end();++i)
                    {
                        int vari=wcsp->getVar(*i)->wcspIndex;
                        int varj=maxvar->wcspIndex;
                        if( included(listeVars[vari] , listeVars[varj]) ){
                            listeVars[(*i)].insert(varj);
                            //--degree[(*i)];
                        }
                    }
                }
            }
            currentRevElimOrder.push_back(maxvar);
            inusedvars.erase(maxvar->wcspIndex);
        }
    }
}

void TreeDecomposition::insert(int sizepart, vector <Variable *> currentRevElimOrder, ConstraintSet currentusedctrs )
{
    int firstCluster = clusters.size();
    for(unsigned int i=0;i<currentRevElimOrder.size();i++) {
        Cluster* c = new Cluster( this );
        c->setId(clusters.size());
        c->addVar( currentRevElimOrder[currentRevElimOrder.size()-i-1] );
        clusters.push_back( c );
    }
    ConstraintSet usedctrs;

    for(unsigned int i=0;i<currentRevElimOrder.size();i++) {
        Cluster* c  = clusters[firstCluster+i];
        c->setPart(sizepart);
        Variable* x = currentRevElimOrder[currentRevElimOrder.size()-i-1];

        ConstraintList* xctrs = x->getConstrs();
        for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
            bool used = usedctrs.find( ctr ) != usedctrs.end();
            if(!used) {
                if (currentusedctrs.find( ctr ) != currentusedctrs.end()) {
                    usedctrs.insert( ctr );
                    c->addCtr(ctr);
                    for(int k=0; k < ctr->arity(); k++) {
                        if (ctr->getVar(k)->unassigned()) {
                            //assert(currentusedvars.find(ctr->getVar(k)->wcspIndex) != currentusedvars.end());
                            c->addVar( ctr->getVar(k) );
                        }
                    }
                }
            }
        }

        for(unsigned int j=i+1;j<currentRevElimOrder.size();j++)
        {
            if(c->isVar(currentRevElimOrder[currentRevElimOrder.size()-j-1]->wcspIndex)) {
                Cluster* cj  = clusters[firstCluster+j];
                TVars::iterator it = c->beginVars();
                while(it != c->endVars()) {
                    cj->addVar( wcsp->getVar(*it) );
                    ++it;
                }
                cj->removeVar(x);
                c->addEdge( cj );
                cj->addEdge( c );
                break;
            }
        }
    }
}

void TreeDecomposition::getElimVarOrder(vector<int> &elimVarOrder)
{
    getRoot()->getElimVarOrder(elimVarOrder);
}

void TreeDecomposition::addDelta(int cyid, EnumeratedVariable *x, Value value, Cost cost)
{
    Cluster* cy = getCluster( cyid );
    Cluster* cx = getCluster( x->getCluster() );
    if(! cy->isDescendant( cx ) ) {
        int ckid,posx;
        assert(x->clusters.size() > 0);
        if (cost != MIN_COST && !deltaModified[x->wcspIndex]) deltaModified[x->wcspIndex] = true;
        x->beginCluster();
        while( x->nextCluster(ckid,posx) ) {
            Cluster* ck = getCluster( ckid );
            if(ck->isDescendant(cy)) {
                if (ToulBar2::verbose >= 2) cout << "add delta " << cost << " to var " << x->wcspIndex << " (cluster " << cx->getId() << ") value " << value << " from subtree " << ck->getId() << " (cluster " << cyid << ")" << endl;
                ck->addDelta(posx, value, cost);
            }
        }
    }
}

// warning! variables are not assigned to the current new solution
// use assignment "a" instead
void TreeDecomposition::newSolution( Cost lb )
{   
    wcsp->setUb(lb);

    TAssign a;

    Cluster* root = getRoot();
    wcsp->restoreSolution(root);
    root->getSolution( a );

    if (ToulBar2::elimDegree>0 && root->getNbVars() == 0) {
        // recorded solutions in clusters containing a single variable eliminated in preprocessing may be wrong due to variable elimination in preprocessing; must be recovered after complete assignment and restoreSolution
        for (unsigned int i=0; i< wcsp->numberOfVariables(); i++) {
            if (wcsp->enumerated(i)) {
                EnumeratedVariable *x = (EnumeratedVariable *) wcsp->getVar(i);
                x->assignWhenEliminated( a[i] );
            }
        }
        wcsp->restoreSolution();
        for (unsigned int i=0; i< wcsp->numberOfVariables(); i++) {
            if (wcsp->enumerated(i)) {
                a[i] = wcsp->getValue(i);
            }
        }
    }
    if (!ToulBar2::allSolutions && !ToulBar2::isZ) wcsp->setSolution(&a);

    if (ToulBar2::showSolutions) {
        TAssign::iterator it = a.begin();
        while(it != a.end()) {
            Value v = it->second;
            cout << v << " ";
            ++it;
        }
        cout << endl;
    }

    if (ToulBar2::writeSolution) {
        ofstream file(ToulBar2::writeSolution);
        if (!file) {
            cerr << "Could not write file " << "solution" << endl;
            exit(EXIT_FAILURE);
        }
        TAssign::iterator it = a.begin();
        while(it != a.end()) {
            Value v = it->second;
            file << v << " ";
            ++it;
        }
        file << endl;
    }
    if (ToulBar2::maxsateval) {
        cout << "o " << lb << endl;
    }
    if(ToulBar2::xmlflag) {
        cout << "o " << lb << endl;
        wcsp->solution_XML();
    }
    else if(ToulBar2::uai || ToulBar2::uaieval) {
        wcsp->solution_UAI(lb);
    }
    // warning: cannot read solution from variable assignments
    // else if(ToulBar2::pedigree){
    // 	ToulBar2::pedigree->printSol(wcsp);
    // }
    // else if(ToulBar2::haplotype){
    //   ToulBar2::haplotype->printSol(wcsp);
    // }

    if (ToulBar2::newsolution) (*ToulBar2::newsolution)(wcsp->getIndex(), wcsp->getSolver());
}

void TreeDecomposition::intersection( TVars& v1, TVars& v2, TVars& vout )
{
    set_intersection( v1.begin(), v1.end(),
            v2.begin(), v2.end(),
            inserter(vout, vout.begin()) );
}

void TreeDecomposition::difference( TVars& v1, TVars& v2, TVars& vout )
{
    set_difference( v1.begin(), v1.end(),
            v2.begin(), v2.end(),
            inserter(vout, vout.begin()) );
}

void TreeDecomposition::sum( TVars& v1, TVars& v2, TVars& vout )
{
    set_union( v1.begin(), v1.end(),
            v2.begin(), v2.end(),
            inserter(vout, vout.begin()) );
}

bool TreeDecomposition::included( TVars& v1, TVars& v2 )
{
    TVars vout;
    set_union( v1.begin(), v1.end(),
            v2.begin(), v2.end(),
            inserter(vout, vout.begin()) );
    return vout.size() == v2.size();
}

void TreeDecomposition::clusterSum( TClusters& v1, TClusters& v2, TClusters& vout )
{
    set_union( v1.begin(), v1.end(),
            v2.begin(), v2.end(),
            inserter(vout, vout.begin()) );
}

void TreeDecomposition::ctrSum( TCtrs& v1, TCtrs& v2, TCtrs& vout )
{
    set_union( v1.begin(), v1.end(),
            v2.begin(), v2.end(),
            inserter(vout, vout.begin()) );
}

bool TreeDecomposition::verify()
{
    if (!ToulBar2::approximateCountingBTD) {
        for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
            Variable* x = wcsp->getVar( i );
            if(x->assigned()) continue;

            Cluster* ci  = clusters[x->getCluster()];
            if(!ci->isVar(x->wcspIndex) || ci->isSepVar(x->wcspIndex)) {
                cout << "cluster: " << ci->getId() << " , var " << x->wcspIndex << endl;
                return false;
            }
            //  	    ConstraintList* xctrs = x->getConstrs();
            //  	    for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
            //              Constraint* ctr = (*it).constr;
            //  			Cluster* cj  = clusters[ctr->getCluster()];
            //              int arity = ctr->arity();
            //              for(i=0;i<arity;i++) {
            //          		Variable* x = ctr->getVar(i);

            //              }
            //  	    }
        }
    }
    return true;
}

void TreeDecomposition::printStats( Cluster* c )
{
    if(!c) return;
    c->printStats();
    TClusters::iterator itj =  c->beginEdges();
    while(itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj;
        printStats(cj);
    }
}

void TreeDecomposition::print( Cluster* c, int recnum )
{
    if(!c) {
        //  		for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
        //  			Variable* x = wcsp->getVar(i);
        //  			x->beginCluster();
        //  			int c,posx;
        //  			cout << x->wcspIndex << " appears in sep {";
        //  			while(x->nextCluster(c,posx)) {
        //  				cout << c << " ";
        //  			}
        //  			cout << "}" << endl;
        //  		}
        if(roots.empty()) return;
        c = * roots.begin();
    }

    for (int i=0; i<recnum; i++) cout << "  ";
    c->print();


    TClusters::iterator ita = c->beginSortedEdges();
    while(ita != c->endSortedEdges()) {
        print( *ita, recnum+1 );
        ++ita;
    }
}

void TreeDecomposition::dump(Cluster* c)
{
    if(!c) {
        char tmpName[256];
        sprintf(tmpName,"%s.info",getWCSP()->getName().c_str());
#ifdef WIN32
        mkdir(tmpName);
#else
        mkdir(tmpName,0777);
#endif

        sprintf(tmpName,"%s.info/root",getWCSP()->getName().c_str());

        ofstream rootFile(tmpName);
        if(roots.empty()) {
            rootFile.close();
            return;
        }
        c = * roots.begin();
        rootFile << c->getId();
        rootFile.close();
    }

    c->dump();

    TClusters::iterator ita = c->beginEdges();
    while(ita != c->endEdges()) {
        dump(*ita);
        ++ita;
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

