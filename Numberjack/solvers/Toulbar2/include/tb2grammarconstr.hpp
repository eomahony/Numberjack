/** \file tb2grammarconstr.hpp
 *  \brief Dynamic programming based global cost function : sgrammar_dp
 */

#ifndef TB2GRAMMARCCONSTR_HPP_
#define TB2GRAMMARCCONSTR_HPP_

#include "tb2dpglobalconstr.hpp"
#include "tb2grammarutils.hpp" 
#include <vector>
#include <fstream>
#include <string>
using namespace std;

class GrammarConstraint : public DPGlobalConstraint {
private:

    // dimension: i x j x |N|
    Cost ***f;
    Cost ***up;
    bool ***marked;

    Cost ***curf;

    // dimension: i x |sigma|
    Cost **u;

    // grammar, assuming in CNF

    /*struct Rule {
        int from;        
        int weight;
        int to[2];
    };

    int nNonTerminals, nTerminals, startSymbol;
    vector<Rule> nonTerm2nonTerm;
    vector<Rule> nonTerm2term;*/
    WCNFCFG cfg;
    Cost top;        

    template<class T>
    void resizeTable(T*** &table) {
        table = new T**[arity() + 1];
        for (int i = 0; i < arity() + 1; i++) {
            table[i] = new T*[arity() + 1];
            for (int j = 0; j < arity() + 1; j++) {
                table[i][j] = new T[cfg.getNumNonTerminals()];
            }
        }
    }

    template<class T>
    void deleteTable(T*** &table) {
        for (int i = 0; i < arity() + 1; i++) {
            for (int j = 0; j < arity() + 1; j++) {
                delete[] table[i][j];
            }
            delete[] table[i];
        }
        delete[] table;
        table = NULL;
    }

    void recomputeTable(Cost*** table, Cost*** upTable = NULL);    
    void recompute();

    Cost unary(int ch, int var, Value v);

protected:
    Cost minCostOriginal();
    Cost minCostOriginal(int var, Value val, bool changed);
    Result minCost(int var, Value val, bool changed);

public:

    static const int WEIGHTED = 1;
    static const int VAR = 0;

    GrammarConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity);
    virtual ~GrammarConstraint();

    Cost eval(const String& s);

    void read(istream & file);
    //void setDefaultViolationCost(Cost cost) {if(configuring) def = cost;}                
    //void setViolationMeasure(int measure) {if(configuring) mode = measure;}            
    WeightedCNFCFG* getGrammar() {return &cfg;}     
    void initMemoization();   

    string getName() {
        return "sgrammar";
    }
    void dump(ostream& os, bool original = true);
};

#endif /* TB2GRAMMARCCONSTR_HPP_ */


/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

