/** \file tb2clusters.hpp
 *  \brief Cluster Tree Decomposition data-structures.
 *
 */

#ifndef TB2CLUSTERS_HPP_
#define TB2CLUSTERS_HPP_

#include "tb2wcsp.hpp"
#include "tb2solver.hpp"
#include "tb2enumvar.hpp"
#include "tb2naryconstr.hpp"

#include <set>
#include <list>


class Cluster;

typedef set<int>	       TVars;
typedef set<Constraint*>   TCtrs;
//typedef map<int,Value>     TAssign;


// sort cluster sons by mean separator size first and by number of variables in their subtree next
struct CmpClusterStruct {
    bool operator() (const Cluster *lhs, const Cluster *rhs) const;
};
typedef set<Cluster*>       TClusters;
typedef set<Cluster*, CmpClusterStruct>       TClustersSorted;

typedef triplet<Cost, Cost, Solver::OpenList >     TPairNG;
typedef pair<Cost,String>   TPairSol;

typedef map<String, TPairNG>  TNoGoods;
typedef map<String, TPairSol> TSols;

// for solution counting :
typedef pair<Cost,BigInteger>	TPairSG;
typedef map<String,TPairSG> TSGoods;






class Separator : public AbstractNaryConstraint
{
private:
    Cluster*					  cluster;
    TVars   					  vars;
    vector< vector<StoreCost> >   delta;    // structure to record the costs that leave the cluster
    StoreInt nonassigned;       			// number of non assigned variables during search
    StoreInt isUsed;
    StoreCost lbPrevious;
    StoreInt optPrevious;

    TNoGoods  					  nogoods;
    TSGoods						  sgoods;	// for solution counting
    TSols  						  solutions;
    DLink<Separator *>            linkSep; // link to insert the separator in PendingSeparator list

    String t;    // temporary buffer for a separator tuple
    String s;    // temporary buffer for a solution tuple

public:
    Separator(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in);
    Separator(WCSP *wcsp);

    void setup(Cluster* cluster_in);

    TVars& getVars() { return vars; }
    int getNbVars() { return vars.size();  }
    bool  is(int i) { return vars.find(i) != vars.end(); }
    TVars::iterator  begin() { return vars.begin(); }
    TVars::iterator  end()   { return vars.end(); }

    void set( Cost clb, Cost cub, Solver::OpenList **open = NULL );
    bool get( Cost& clb, Cost& cub, Solver::OpenList **open = NULL );

    void setSg( Cost c, BigInteger nb );
    BigInteger getSg( Cost& res, BigInteger& nb );

    void solRec( Cost ub );
    bool solGet( TAssign& a, String& sol );

    void resetLb();
    void resetUb();

    void queueSep() { wcsp->queueSeparator(&linkSep); }
    void unqueueSep() { wcsp->unqueueSeparator(&linkSep); }

    void addDelta( unsigned int posvar, Value value, Cost cost ) {
        assert( posvar < vars.size() );
        delta[posvar][value] += cost;
    }
    Cost getCurrentDelta(); // separator variables may be unassigned

    bool used() { return isUsed; }

    void assign(int varIndex);
    void propagate();

    // NaryConstraint methods not used by Separator class
    double computeTightness() { return 0; }
    bool   verify() {return true;}
    void   increase(int index) {}
    void   decrease(int index) {}
    void   remove(int index) {}
    //set<Constraint*> subConstraint(){TCtrs s; return s;}
    void   print(ostream& os);
};


class Cluster
{
private:
    TreeDecomposition*  td;
    WCSP*				  wcsp;
    int                 id; // corresponding to the vector index of the cluster in the tree decomposition
    TVars				  vars; // contains all variables inside a cluster including separator variables
    TCtrs			      ctrs; // intermediate usage by bucket elimination (DO NOT USE)
    TClusters           edges; // adjacent clusters (includes parent cluster before makeRooted is done)
    TClustersSorted     sortedEdges; // cluster sons are sorted after makeRooted is done

    Cluster*			  parent; // parent cluster
    TClusters           descendants; // set of cluster descendants (including itself)
    TVars               varsTree; // set of variables in cluster descendants (including itself)
    vector<bool>        quickdescendants;

    Separator* 	      sep; // associated separator with parent cluster
    StoreCost           lb; // current cluster lower bound deduced by propagation
    Cost                ub; // current cluster best known solution cost
    Cost				  lbRDS; // global cluster lower bound found by RDS
    StoreInt			  active; // unactive if a nogood including this cluster has been used by propagation

    StoreBigInteger     countElimVars;
    int 			  	  num_part; //for approximation: number of the corresponding partition

public:
    Cluster (TreeDecomposition* tdin);
    ~Cluster();

    void          setup();

    int           getId() const { return id; }
    void          setId(int iid) { id=iid; }

    WCSP* 		getWCSP() { return wcsp; }

    Separator*    getSep() { return sep; }
    void 			setSep( Separator* sepin ) { sep = sepin; }
    int			sepSize() const { if(sep) return sep->arity(); else return 0; }
    void          deconnectSep(); // deconnect all the constraints on separator variables and assigns separator variables to their support value
    void 			deconnectDiff(TCtrs listCtrsTot,TCtrs listCtrs);
    bool 			isSepVar( int i ) { if(!sep) return false; return sep->is(i); }

    Solver::CPStore*      cp;         // choice point cache for open nodes related to this cluster
    Solver::OpenList*     open;       // list of open nodes related to this cluster
    Long          hbfsGlobalLimit;    // global limit on number of backtracks for hybrid search on the subproblem rooted to this cluster
    Long          hbfsLimit;          // local limit on number of backtracks for hybrid search on this cluster only
    Long          nbBacktracks;       // current number of backtracks related to this cluster
    Long          getNbBacktracksClusterTree() const {Long res = nbBacktracks; for (TClusters::const_iterator iter = beginEdges(); iter != endEdges(); ++iter) res += (*iter)->getNbBacktracksClusterTree(); return res;}
    vector<Cluster *>   sons;         // copy of edges allowing sorting

    bool 			isVar( int i ) { TVars::iterator it = vars.find(i); return it != vars.end(); }
    int			getNbVars() const { return vars.size(); }
    TVars&		getVars() { return vars; }
    int           getNbVarsTree() const { return varsTree.size(); }
    TVars&		getVarsTree() { return varsTree; }
    TCtrs 		getCtrsTree();
    void 			addVars( TVars& vars );
    void 			addVar( Variable* x );
    void 			removeVar( Variable* x );

    void 			setParent(Cluster* p) { parent = p; }
    Cluster*		getParent() { return parent; }
    TClusters&	getEdges() { return edges; }
    void 			addEdges( TClusters& cls );
    void 			addEdge( Cluster* c );
    void 			removeEdge( Cluster* c );
    TClusters::iterator removeEdge( TClusters::iterator it );

    TClusters&	getDescendants() { return descendants; }
    bool  	    isEdge( Cluster* c );
    void          accelerateDescendants();
    bool			isDescendant( Cluster* c ) { return quickdescendants[c->getId()]; }

    TCtrs&		getCtrs() { return ctrs; }
    void 			addCtrs( TCtrs& ctrsin );
    void 			addCtr( Constraint* c );
    void 			sum( TCtrs& c1, TCtrs& c2, TCtrs& ctout );

    bool			isActive() const { int a = active; return a == 1; }
    void 			deactivate();
    void 			reactivate();

    Cost			getLb()  { return lb; }
    void			setLb(Cost c)  { lb = c; }
    void 			increaseLb( Cost addToLb ) { lb += addToLb; }
    Cost          getUb() const { return ub; }
    void          setUb(Cost c) {ub = c;}
    Cost		    getLbRDS() { Cost delta = getCurrentDelta(); return MAX(lbRDS - delta, MIN_COST); }
    void			setLbRDS(Cost c)  {assert(!sep || sep->getCurrentDelta()==MIN_COST); lbRDS = c; }
    Cost	        getLbRec() const;
    Cost	        getLbRecRDS();

    void          addDelta( int posvar, Value value, Cost cost ) {assert(sep); sep->addDelta(posvar,value,cost); }
    Cost          getCurrentDelta() {return (sep)?sep->getCurrentDelta():MIN_COST;}

    void          nogoodRec( Cost clb, Cost cub, Solver::OpenList **open = NULL ) { if(sep) sep->set(clb,cub,open); }
    bool          nogoodGet( Cost &clb, Cost &cub, Solver::OpenList **open = NULL ) { return sep->get(clb,cub,open); }

    void          resetLbRec();
    void          resetUbRec(Cluster *rootCluster);

    void			sgoodRec( Cost c, BigInteger nb) { if(sep) sep->setSg(c,nb);}
    BigInteger		sgoodGet( ){Cost c = MIN_COST; BigInteger nb; sep->getSg(c,nb); return nb; }
    BigInteger 		getCount() {return countElimVars;}
    void 			multCount(unsigned int s) {countElimVars = countElimVars * s;}
    void 			cartProduct(BigInteger& cartProd);
    int			getPart(){return num_part;}
    void			setPart(int num){num_part=num;}

    void          solutionRec(Cost c) { setUb(c); if(sep) sep->solRec(c); }
    void 		    getSolution( TAssign& sol ); // updates sol by the recorded solution found for a separator assignment also given in sol

    void 			setWCSP2Cluster();   // sets the WCSP to the cluster problem, deconnecting the rest
    void          getElimVarOrder(vector<int> &elimVarOrder);

    TVars::iterator beginVars() { return vars.begin(); }
    TVars::iterator endVars()   { return vars.end(); }
    TVars::iterator beginVarsTree() { return varsTree.begin(); }
    TVars::iterator endVarsTree()   { return varsTree.end(); }
    TVars::iterator beginSep() { return sep->begin(); }
    TVars::iterator endSep()   { return sep->end(); }
    TCtrs::iterator beginCtrs() { return ctrs.begin(); }
    TCtrs::iterator endCtrs()   { return ctrs.end(); }
    TClusters::iterator beginEdges() const { return edges.begin(); }
    TClusters::iterator endEdges() const { return edges.end(); }
    TClusters::iterator beginDescendants() { return descendants.begin(); }
    TClusters::iterator endDescendants()   { return descendants.end(); }

    void sortEdgesRec() {for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) (*iter)->sortEdgesRec(); TClustersSorted tmpset(edges.begin(), edges.end(), CmpClusterStruct()); sortedEdges = tmpset;}
    TClusters::iterator beginSortedEdges() const { return sortedEdges.begin(); }
    TClusters::iterator endSortedEdges() const { return sortedEdges.end(); }

    void print();
    void dump();
    void printStats() { if(!sep) return; sep->print(cout); }

    void printStatsRec() {
        TClusters::iterator it = beginSortedEdges();
        while(it != endSortedEdges()) {
            (*it)->sep->print(cout);
            (*it)->printStatsRec();
            ++it;
        }
    }
};


class TreeDecomposition  {
private:
    WCSP*			  wcsp;
    vector<Cluster*>  clusters;
    list<Cluster*> 	  roots; // intermediate list used by makeRooted method, only one root at the end

    Cluster*          rootRDS; // root cluster of the current RDS iteration

    StoreInt   		  currentCluster; // used to restrict local propagation (NC) and boosting by variable elimination to the current cluster's subtree
    vector<StoreInt>  deltaModified; // accelerator avoiding unecessary checks to delta structure if it is empty (Boolean value)

public:

    TreeDecomposition(WCSP* wcsp_in);

    WCSP* 		getWCSP() { return wcsp; }

    int			getNbOfClusters() { return clusters.size(); }
    Cluster*	getCluster( int i ) { assert( 0 <= i && i < (int) clusters.size() ); return clusters[i]; }
    Cluster*   	var2Cluster( int v );

    void setCurrentCluster(Cluster *c) {currentCluster = c->getId();}
    Cluster* getCurrentCluster() {return getCluster(currentCluster);}

    bool isInCurrentClusterSubTree(int idc);
    bool isActiveAndInCurrentClusterSubTree(int idc);

    //main function to build a cluster tree/path decomposition:
    // - builds a bucket for each variable following a given variable elimination order
    // - builds a tree/path decomposition from the buckets
    // - associate constraints to clusters, with special treatment for ternary constraints (duplicate flag)
    void buildFromOrder();
    void buildFromOrderNext(vector<int> &order);
    void getElimVarOrder(vector<int> &elimVarOrder);
    void treeFusions();                  	// merges all redundant clusters
    bool treeFusion();          		    // one fusion step
    void pathFusions(vector<int> &order);   // builds a path decomposition of clusters from a given order

    void buildFromOrderForApprox();	//builds a decomposition for approximation solution counting
    void maxchord(int sizepart, vector<int> &order, set<Constraint*> &totalusedctrs, TVars &inusedvars, TVars &currentusedvars, vector<Variable *> &currentRevElimOrder,set<Constraint*> &currentusedctrs);
    void insert(int sizepart, vector <Variable *> currentRevElimOrder, set<Constraint *> currentusedctrs );

    void fusion( Cluster* ci, Cluster* cj );
    bool reduceHeight( Cluster* c, Cluster *father );
    int getNextUnassignedVar(TVars *vars);
    int getVarMinDomainDivMaxWeightedDegree(TVars *vars);
    void splitClusterRec( Cluster* c,  Cluster* father, unsigned int maxsize );
    TVars boostingVarElimRec( Cluster* c,  Cluster* father,  Cluster* grandfather, unsigned int maxsize );
    void mergeClusterRec( Cluster* c,  Cluster* father, unsigned int maxsepsize, unsigned int minpropervar );
    void heuristicFusionRec( Cluster* c, Cluster* noc );

    void makeDescendants( Cluster* c );
    bool isDescendant( Variable* x, Variable* y ) { return getCluster(x->getCluster())->isDescendant(getCluster(y->getCluster())); }

    int      makeRooted(); // defines a rooted cluster tree decomposition from an undirected one
    void     makeRootedRec( Cluster* c,  TClusters& visited );
    Cluster* getBiggerCluster( TClusters& visited );
    Cluster* getRoot() {return roots.front();}
    Cluster* getRootRDS() {return rootRDS; }
    void setRootRDS(Cluster* rdsroot) { rootRDS = rdsroot; }

    int height( Cluster* r );
    int height( Cluster* r, Cluster* father );

    void intersection( TVars& v1, TVars& v2, TVars& vout );
    void difference( TVars& v1, TVars& v2, TVars& vout );
    void sum( TVars& v1, TVars& v2, TVars& vout );
    bool included( TVars& v1, TVars& v2 ); 	// true if v1 is included in v2
    void clusterSum( TClusters& v1, TClusters& v2, TClusters& vout );
    void ctrSum( TCtrs& v1, TCtrs& v2, TCtrs& vout );

    bool isDeltaModified(int varIndex) { return deltaModified[varIndex]; }
    Cost getLbRecRDS() {
        Cluster* c = getCluster(currentCluster);
        Cost res = c->getLbRecRDS();
        return MAX(res,c->getLbRDS());
    }
    void addDelta(int c, EnumeratedVariable *x, Value value, Cost cost);
    void newSolution( Cost lb );

    bool verify();

    void print( Cluster* c = NULL, int recnum = 0);
    void printStats( Cluster* c = NULL );
    void dump(Cluster* c = NULL);
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

