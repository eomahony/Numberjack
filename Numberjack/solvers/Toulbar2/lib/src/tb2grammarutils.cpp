#include "tb2grammarutils.hpp"   
#include <set>
#include <algorithm>

using namespace std;

void WCNFCFG::addVariableMeasure(int violationCost) {
            
    // Assume A_i -> a_i   
    // 1) for all prod. rules C->A_iB, add C->A_jB with weight = w[C->A_jB] + vioCost, j != i
    // 2) for all prod. rules C->BA_i, add C->BA_j with weight = w[C->A_jB] + vioCost, j != i
    set<int> nonTerms;
    for (vector<WCNFRule>::iterator p = termProd.begin();p != termProd.end();++p) {         
        nonTerms.insert(p->from);
    }                
    for (int i=0;i<2;i++) {
        vector<WCNFRule> prods(nonTermProd);
        for (vector<WCNFRule>::iterator p = prods.begin();p != prods.end();++p) {                 
            if (nonTerms.find(p->to[i]) != nonTerms.end()) {                           
                for (set<int>::iterator A = nonTerms.begin();A != nonTerms.end();++A) {         
                    if (*A != p->to[i]) {
                        WCNFRule rule = *p;
                        rule.to[i] = *A;
                        rule.weight += violationCost;
                        nonTermProd.push_back(rule);
                    }
                }
            }
         }                        
    }
    
    sort(nonTermProd.begin(), nonTermProd.end());
    nonTermProd.erase(unique(nonTermProd.begin(), nonTermProd.end()), nonTermProd.end());
                
}

void WCNFCFG::print(ostream &ofs) {
     for (vector<WCNFRule>::iterator p = nonTermProd.begin(); p != nonTermProd.end(); ++p) {
         ofs << p->from << "->" << p->to[0] << " " << p->to[1] << ": " << p->weight << "\n";
     }   
      for (vector<WCNFRule>::iterator p = termProd.begin(); p != termProd.end(); ++p) {
         ofs << p->from << "->" << p->to[0] << ": " << p->weight << "\n";
     }  
}

void WCNFCFG::dump(ostream& os, bool original)
{
    assert(original); //TODO: case original is false
    os << nNonTerminals << " " << nTerminals << " " << startSymbol << endl;
    os << nonTermProd.size() + termProd.size() << endl;
    for (vector<WCNFRule>::iterator p = nonTermProd.begin(); p != nonTermProd.end(); ++p) {
        if (p->weight == 0) os << "1 " << p->from << " " << p->to[0] << " " << p->to[1] << endl;
        else os << "3 " << p->from << " " << p->to[0] << " " << p->to[1] << " " << p->weight << endl;
    }
    for (vector<WCNFRule>::iterator p = termProd.begin(); p != termProd.end(); ++p) {
        if (p->weight == 0) os << "0 " << p->from << " " << p->to[0] << endl;
        else os << "2 " << p->from << " " << p->to[0] << " " << p->weight << endl;
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

