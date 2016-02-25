#ifndef GRAMMARUTILS
#define GRAMMARUTILS

#include <vector>
#include <iostream>
#include <map>
#include "tb2types.hpp"

using namespace std;

class WeightedCNFCFG {
public:
    virtual void setStartSymbol(int start) = 0;    
    virtual void setNumNonTerminals(int num) = 0;    
    virtual void setNumTerminals(int num) = 0;        
    virtual void addProduction(int A, int B, int C, int weight) = 0;    
    virtual void addProduction(int A, int v, int weight) = 0;
};

struct WCNFRule {
   int from;        
   int weight;
   int to[2];
   WCNFRule(int A=0, int B=0, int C=0, int w=0):from(A), weight(w) {
       to[0] = B; to[1] = C;
   }
   
   bool operator==(const WCNFRule &rule) const {
       return (from == rule.from) && (to[0] == rule.to[0]) && (to[1] == rule.to[1]);
   }
   
   bool operator<(const WCNFRule &rule) const {
       if (from < rule.from) return true;
       if ((from == rule.from) && (to[0] < rule.to[0])) return true;
       if ((from == rule.from) && (to[0] == rule.to[0]) && (to[1] < rule.to[1])) return true;
       if ((from == rule.from) && (to[0] == rule.to[0]) && (to[1] == rule.to[1]) && (weight < rule.weight) ) return true;
       return false;
   }
};

class WCNFCFG : public WeightedCNFCFG { 
    
// Assumption: if A->a and B-> a, then A = B    
    
private:
    
    int nNonTerminals, nTerminals, startSymbol;
    vector<WCNFRule> nonTermProd; // for rules: A -> BC        
    vector<WCNFRule> termProd; // for rules: A -> a
    
    map<Value, int> valIndex;
    map<int, Value> indexValue;
    
public:
            
    WCNFCFG(int nNonTerminals=0)
            :nNonTerminals(nNonTerminals), nTerminals(0), startSymbol(0) {}
    
    ~WCNFCFG() {}
    
    void setStartSymbol(int start) { // Only one start symbol
        startSymbol = start;
    }
    
    void setNumNonTerminals(int num) { 
        nNonTerminals = num;
    }
    
    void setNumTerminals(int num) { 
        //nTerminals = num;
    }
    
    inline int getStartSymbol() { return startSymbol; }
    
    inline int getNumNonTerminals() { return nNonTerminals;}
    
    inline int getNumTerminals() { return nTerminals; }
    
    inline int toIndex(Value v) {return ((valIndex.find(v)==valIndex.end())?-1:valIndex[v]);}
    inline int toValue(int index) {return indexValue[index];}
            
    void addProduction(int A, int B, int C, int weight) {
        WCNFRule rule(A, B, C, weight);
        nonTermProd.push_back(rule);        
    }
    
    void addProduction(int A, int v, int weight) {                
        WCNFRule rule(A, v, -1, weight);        
        termProd.push_back(rule);
        valIndex[v] = nTerminals;
        indexValue[nTerminals] = v;
        nTerminals++;
    }
    
    void addRedundantRuleTo(Value v) {  
        int A = nNonTerminals;        
        addProduction(A, v, 0);
        nNonTerminals++;                
    }
                    
    void addVariableMeasure(int violationCost); // Convert to a weighted CNF s.t. it includes variable-based violation
    
    /* iterator for rules like A->v*/   
    typedef vector<WCNFRule>::iterator TermProdIterator;    
    inline TermProdIterator beginTermProd() {return termProd.begin();}
    inline TermProdIterator endTermProd() {return termProd.end();}
    
    /* iterator for rules like A->BC*/   
    typedef vector<WCNFRule>::iterator NonTermProdIterator;
    inline NonTermProdIterator beginNonTermProd() {return nonTermProd.begin();}
    inline NonTermProdIterator endNonTermProd() {return nonTermProd.end();}
    
    void print(ostream &ofs);
    
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

