#ifndef TB2NARYCONSTR_HPP_
#define TB2NARYCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

#include <set>


class NaryConstraint : public AbstractNaryConstraint
{
public:

    Cost default_cost;          // default cost returned when tuple t is not found in TUPLES (used by function eval(t)
    bool store_top; 		    // this is true when default_cost < getUb() meaning that tuples with cost greater than ub must be stored
    StoreInt nonassigned;       // nonassigned variables during search, must be backtrackable (storeint) !

    String iterTuple;
    String evalTuple;

    vector<EnumeratedVariable::iterator> it_values;

    vector<Long> conflictWeights;
    Long getConflictWeight() const {return Constraint::getConflictWeight();} 
    Long getConflictWeight(int varIndex) const {assert(varIndex>=0);assert(varIndex<arity());return conflictWeights[varIndex]+Constraint::getConflictWeight();} 
    void incConflictWeight(Constraint *from) {
        //assert(fromElim1==NULL);
        //assert(fromElim2==NULL);
        if (from==this) {
            Constraint::incConflictWeight(1);
        } else if (deconnected()) {
            for (int i=0; i<from->arity(); i++) {
                int index = getIndex(from->getVar(i));
                if (index>=0) { // the last conflict constraint may be derived from two binary constraints (boosting search), each one derived from an n-ary constraint with a scope which does not include parameter constraint from
                    assert(index < arity());
                    conflictWeights[index]++;
                }
            }
        }
    }

    void firstlex();
    bool nextlex( String& t, Cost& c);

    NaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
    NaryConstraint(WCSP *wcsp);

    virtual int size() = 0;
    virtual void setTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL ) = 0;

    virtual void addtoTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL ) = 0;

    virtual void setDefCost( Cost df ) = 0;

    virtual void setTuple( int* tin, Cost c, EnumeratedVariable** scope_in )
    {
        Char* buf = new Char [arity_ + 1];
        for(int i=0;i<arity_;i++) buf[i] = tin[i]+CHAR_FIRST;
        buf[arity_] = '\0';
        String str = String(buf);
        setTuple( str, c, scope_in );
        delete [] buf;
    }

    virtual void addtoTuple( int* tin, Cost c, EnumeratedVariable** scope_in )
    {
        Char* buf = new Char [arity_ + 1];
        for(int i=0;i<arity_;i++) buf[i] = tin[i]+CHAR_FIRST;
        buf[arity_] = '\0';
        String str = String(buf);
        addtoTuple( str, c, scope_in );
        delete [] buf;
    }

    void reconnect() {
        if (deconnected()) {
            nonassigned = arity();
            AbstractNaryConstraint::reconnect();
        }
    }

    virtual Cost eval( String& s ) = 0;
    Cost evalsubstr( String& s, Constraint* ctr );

    void assign(int varIndex);

    void projectNary();
    void projectNaryTernary(TernaryConstraint* xyz);
    void projectNaryBinary(BinaryConstraint* xy);

    void propagate() {
        for(int i=0;connected() && i<arity_;i++) {
            if (getVar(i)->assigned()) assign(i);
        }
    };

    virtual void project( EnumeratedVariable* x ) = 0;

    bool   verify() {return true;}
    void   increase(int index) {}
    void   decrease(int index) {}
    void  remove(int index) {}

    void starrule(String& t, Cost minc);
    void projectFromZero(int index);

    virtual void print(ostream& os) {}
};


class NaryConstraintMap : public NaryConstraint
{
    typedef map<String,Cost> TUPLES;
    TUPLES* pf;


public:

    NaryConstraintMap(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
    NaryConstraintMap(WCSP *wcsp);
    virtual ~NaryConstraintMap();

    bool extension() const {return true;}

    TUPLES* getpf() {return pf;}
    int size() {return pf->size();}
    Long space() const {return (Long) pf->size() * (sizeof(Cost) + arity()*sizeof(Char));}
    bool consistent( String& t );
    Cost eval( String& s );
    Cost eval( String& s, EnumeratedVariable** scope_in );

    Cost getDefCost() { return default_cost; }
    void setDefCost( Cost df ) { default_cost = df; }
    void keepAllowedTuples( Cost df );

    set<Constraint*>* filters;
    void resetFilters();
    void fillFilters();

    void project( EnumeratedVariable* x );
    void sum( NaryConstraintMap* nary );
    double computeTightness();

    TUPLES::iterator  tuple_it;

    void first();
    bool next( String& t, Cost& c);

    void first(EnumeratedVariable* a, EnumeratedVariable* b);
    bool separability(EnumeratedVariable* a, EnumeratedVariable* b);
    void separate(EnumeratedVariable *a, EnumeratedVariable *c);

    void setTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL );
    void addtoTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL );
    void setInfiniteCost(Cost ub);
    void insertSum( String& t1, Cost c1, Constraint* ctr1, String t2, Cost c2, Constraint* ctr2, bool bFilters = false );
    void permute( EnumeratedVariable** scope_in );

    void projectxy( EnumeratedVariable* x, EnumeratedVariable* y, TUPLES& fproj);
    void projectxyz( EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, TUPLES& fproj);
    void preproject3();
    void preprojectall2();

    void fillRandom();
    void print(ostream& os);
    void dump(ostream& os, bool original = true);

};



class Trie;

class TrieNode {

public:
    TrieNode();

    void iniLeaf(Char *w);
    void iniNonLeaf(Char ch);

    Cost c;

private:
    bool leaf, endOfWord;
    Char *letters;
    Char *word;

    TrieNode **ptrs;
    friend class Trie;
};


class Trie {
public:
    Trie() : notFound(-1) {}
    Trie(Char*, Cost c);
    void insert(Char*, Cost c);
    TrieNode* find(const Char*);

    void printTrie();

private:
    TrieNode *root;
    const int notFound;
    Char prefix[80];
    int  position(TrieNode*,Char);
    void addCell(Char,TrieNode*,int);
    TrieNode* createLeaf(Char,Char*,TrieNode*);
    void printTrie(int,TrieNode*,Char*);
};





class NaryConstrie : public NaryConstraint
{

public:
    Trie* f;


    NaryConstrie(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
    NaryConstrie(WCSP *wcsp);
    virtual ~NaryConstrie();

    int size() {exit(EXIT_FAILURE);return 0;} // not implemented!!!

    void setTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL );
    void addtoTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL );

    void project( EnumeratedVariable* x ) {};

    double computeTightness() { return 0; }

    Cost eval( String& s );

    void print(ostream& os);

};



#endif /*TB2NARYCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

