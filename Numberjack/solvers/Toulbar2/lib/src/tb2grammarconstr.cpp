#include "tb2grammarconstr.hpp"
#include <vector>
#include <fstream>
#include <string>
#include <stack>
#include <set>
using namespace std;

GrammarConstraint::GrammarConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity) : DPGlobalConstraint(wcsp, scope, arity) {
    modeEnum["var"] = GrammarConstraint::VAR;
    modeEnum["weight"] = GrammarConstraint::WEIGHTED;    
}

GrammarConstraint::~GrammarConstraint(void) {
    deleteTable(f);
    deleteTable(curf);
    deleteTable(marked);
    deleteTable(curf);

    for (int i = 0; i < arity(); i++) {
        delete[] u[i];
    }
    delete[] u;
}

void GrammarConstraint::read(istream & file) {

    string str;
    file >> str >> def;

    /*if (str == "var") mode = VAR;
    else mode = WEIGHTED;*/
    this->setSemantics(str);

    // input the grammar, assuming already in CNF with	    
    // 1) the terminal and non-terminal are separated    
    // 2) if A-> a and B->a, then A = B

    // enter the number of non-terminals, terminals, and start symbol
    int nNonTerminals, nTerminals, startSymbol;
    file >> nNonTerminals >> nTerminals >> startSymbol;

    cfg.setNumNonTerminals(nNonTerminals);
    cfg.setNumTerminals(nTerminals);
    cfg.setStartSymbol(startSymbol);

    // enter the number of rules
    int nRules;
    file >> nRules;

    // Rule format:
    // for rule A->BC: 1 A B C <weight>
    // for rule A->a : 0 A a <weight>

    int type = 0;
    for (int i = 0; i < nRules; i++) {
        file >> type;
        switch (type) {
        case 0:
        {
            int A, v;
            file >> A >> v;
            cfg.addProduction(A, v, 0);
            break;
        }
        case 1:
        {
            int A, B, C;
            file >> A >> B >> C;
            cfg.addProduction(A, B, C, 0);
            break;
        }
        case 2:
        {
            int A, v, w;
            file >> A >> v >> w;
            cfg.addProduction(A, v, w);
            break;
        }
        case 3:
        {
            int A, B, C, w;
            file >> A >> B >> C >> w;
            cfg.addProduction(A, B, C, w);
            break;
        }
        default:
            cerr << "Error occur in reading grammar()" << endl;
            exit(1);
        }
    }

}

void GrammarConstraint::dump(ostream& os, bool original)
{
    if (original) {
        os << arity_;
        for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
    } else {
        os << nonassigned;
        for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
    }
    os << " -1 sgrammardp " << ((mode == VAR)?"var":"weight") << " " << def << endl;
    cfg.dump(os, original);
}

void GrammarConstraint::initMemoization() {

    if (mode == VAR) {
        set<int> allValues;
        for (int i = 0; i < arity(); i++) {        
            EnumeratedVariable *x = scope[i];
            for (EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
                allValues.insert(*it);
            }
        }      

        for (set<int>::iterator it=allValues.begin();it != allValues.end();++it) {
            if (cfg.toIndex(*it) == -1) cfg.addRedundantRuleTo(*it);
        }        

        cfg.addVariableMeasure(def);

    }    

    top = max(wcsp->getUb(), MAX_COST);

    //Create tables
    resizeTable(f);
    resizeTable(up);
    resizeTable(curf);
    resizeTable(marked);

    u = new Cost*[arity()+1];
    for (int i = 0; i < arity(); i++) {
        u[i] = new Cost[cfg.getNumNonTerminals()];
    }
}

Cost GrammarConstraint::minCostOriginal() {
    int n = arity();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < cfg.getNumTerminals(); j++) {
            u[i][j] = top;
            EnumeratedVariable *x = scope[i];
            for (EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
                if (u[i][j] > unary(cfg.toValue(j), i, *it))
                    u[i][j] = unary(cfg.toValue(j), i, *it);
            }
        }
    }

    recomputeTable(curf);

    int minCost = curf[0][n - 1][cfg.getStartSymbol()];

    return minCost;

}

Cost GrammarConstraint::minCostOriginal(int var, Value val, bool changed) {    
    return minCost(var, val, changed).first;
}

Cost GrammarConstraint::eval(const String& s) {
    int n = arity();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < cfg.getNumTerminals(); j++) {
            u[i][j] = unary(cfg.toValue(j), i, s[i] - CHAR_FIRST);
        }
    }

    recomputeTable(curf);
    int minCost = curf[0][n - 1][cfg.getStartSymbol()];           

    return minCost - projectedCost;
}

void GrammarConstraint::recompute() {    
    int n = arity();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < cfg.getNumTerminals(); j++) {
            u[i][j] = top;
            EnumeratedVariable *x = scope[i];
            for (EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
                if (u[i][j] > unary(cfg.toValue(j), i, *it)) {
                    u[i][j] = unary(cfg.toValue(j), i, *it);                                        
                }                
            }
        }
    }

    recomputeTable(f, up);

}

DPGlobalConstraint::Result GrammarConstraint::minCost(int var, Value val, bool changed) {
    if (changed) recompute();    

    int n = arity();
    Cost minCost = wcsp->getUb();

    for (WCNFCFG::TermProdIterator r = cfg.beginTermProd();
            r != cfg.endTermProd(); r++) {
        if ((r->to[0] == val) && marked[var][var][r->from]) {
            minCost = min(minCost,
                    unary(r->to[0], var, val) + r->weight - up[var][var][r->from] + f[0][n - 1][cfg.getStartSymbol()]);     
        }
    }            

    return DPGlobalConstraint::Result(minCost, NULL);
}

void GrammarConstraint::recomputeTable(Cost*** table, Cost*** upTable) {
    int n = arity();


    for (int i = 0; i < n; i++) {
        for (int A = 0; A < cfg.getNumNonTerminals(); A++) {
            table[i][i][A] = top;
        }
        /*for (vector<Rule>::iterator r = nonTerm2term.begin(); r != nonTerm2term.end(); r++) {

            if (table[i][i][r->from] > u[i][r->to[0]] + r->weight) {
                table[i][i][r->from] = u[i][r->to[0]] + r->weight;
            }

        }*/
        for (WCNFCFG::TermProdIterator r = cfg.beginTermProd();r != cfg.endTermProd(); ++r) {
            if (table[i][i][r->from] > u[i][cfg.toIndex(r->to[0])] + r->weight) {
                table[i][i][r->from] = u[i][cfg.toIndex(r->to[0])] + r->weight;
            }
        }
    }

    for (int len = 2; len <= n; len++) {
        for (int i = 0; i < n - len + 1; i++) {
            int j = i + len - 1;
            for (int A = 0; A < cfg.getNumNonTerminals(); A++) {
                table[i][j][A] = top;
            }
            /*for (vector<Rule>::iterator r = nonTerm2nonTerm.begin(); r != nonTerm2nonTerm.end(); r++) {
                for (int k = i; k < j; k++) {
                    int tmp = table[i][k][r->to[0]] + table[k + 1][j][r->to[1]] + r->weight;
                    if (table[i][j][r->from] > tmp) {
                        table[i][j][r->from] = tmp;
                    }
                }
            }*/
            for (WCNFCFG::NonTermProdIterator r = cfg.beginNonTermProd();r != cfg.endNonTermProd(); ++r) {                
                for (int k = i; k < j; k++) {
                    Cost tmp = table[i][k][r->to[0]] + table[k + 1][j][r->to[1]] + r->weight;                    
                    table[i][j][r->from] = min(table[i][j][r->from], tmp);                                        
                }
            }
        }
    }

    if (upTable != NULL) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int A = 0; A < cfg.getNumNonTerminals(); A++) {
                    marked[i][j][A] = false;
                    upTable[i][j][A] = -top;
                }
            }
        }
        upTable[0][n - 1][cfg.getStartSymbol()] = table[0][n - 1][cfg.getStartSymbol()];
        marked[0][n - 1][cfg.getStartSymbol()] = true;        
        for (int len = n; len >= 2; len--) {
            for (int i = 0; i < n - len + 1; i++) {
                int j = i + len - 1;
                //for (vector<Rule>::iterator r = nonTerm2nonTerm.begin(); r != nonTerm2nonTerm.end(); r++) {
                for (WCNFCFG::NonTermProdIterator r = cfg.beginNonTermProd();r != cfg.endNonTermProd(); ++r) {                
                    if (marked[i][j][r->from]) 
                    {
                        for (int k = i; k < j; k++) {
                            Cost tmp = table[i][k][r->to[0]] + table[k + 1][j][r->to[1]] + r->weight;
                            //if (tmp <= upTable[i][j][r->from])                                 
                            {
                                marked[i][k][r->to[0]] = true;
                                tmp = upTable[i][j][r->from] - table[k + 1][j][r->to[1]] - r->weight;
                                upTable[i][k][r->to[0]] = max(upTable[i][k][r->to[0]], tmp);                                

                                marked[k + 1][j][r->to[1]] = true;
                                tmp = upTable[i][j][r->from] - table[i][k][r->to[0]] - r->weight;
                                upTable[k + 1][j][r->to[1]] = max(upTable[k + 1][j][r->to[1]], tmp);                                

                            }
                        }
                    }
                }
            }
        }
    }
}

Cost GrammarConstraint::unary(int ch, int var, Value v) {       
    EnumeratedVariable *x = scope[var];
    Cost ucost = (v == ch) ? (-deltaCost[var][x->toIndex(v)]) : top;
    return ucost;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

