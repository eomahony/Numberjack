#include "tb2amongconstr.hpp"
#include <algorithm>
#include <sstream>

using namespace std;

AmongConstraint::AmongConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity) 
: DPGlobalConstraint(wcsp, scope, arity)
, f(NULL)
, invf(NULL)
, curf(NULL)
, minBarU(NULL)
, minU(NULL)
, ub(0), lb(0) {}

AmongConstraint::~AmongConstraint() {
    deleteTable(f);
    deleteTable(curf);
    deleteTable(invf);
}

void AmongConstraint::read(istream & file) {
    string str;
    file >> str >> def;

    if (str != "var") {
        cout << "Error in reading samong()\n";
        exit(1);
    }

    file >> lb >> ub;

    int nVal;
    file >> nVal;
    for (int i=0;i<nVal;i++) {
        Value tmp;
        file >> tmp;
        V.insert(tmp);
    }

}

void AmongConstraint::dump(ostream& os, bool original)
{
    assert(original); //TODO: case original is false
    if (original) {
        os << arity_;
        for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
    } else {
        os << nonassigned;
        for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
    }
    os << " -1 samong var " << def << " " << lb << " " << ub << endl;
    os << V.size();
    for (set<Value>::iterator it = V.begin(); it != V.end(); ++it) {
        os << " " << *it;
    }
    os << endl;
}

void AmongConstraint::initMemoization() {

    if (lb > ub) {
        cout << "Error in samong()\n";
        exit(1);
    }

    int n = arity();

    resizeTable(f, n+1, ub+1);
    resizeTable(invf, n+1, ub+1);    
    resizeTable(curf, n+1, ub+1);           

    minBarU = new UnaryTableCell[n+1]; 
    minU = new UnaryTableCell[n+1]; 

    top = max(MAX_COST, wcsp->getUb());
}

Cost AmongConstraint::minCostOriginal() {

    int n = arity();

    minBarU[0].val = minU[0].val = -1;
    for (int i = 1; i <= n; i++) {
        int minu, minbaru;
        minu = minbaru = wcsp->getUb();
        EnumeratedVariable *x = (EnumeratedVariable*)getVar(i-1);
        for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
            int uCost(0), baruCost(def);
            if (V.find(*v) == V.end()) {
                uCost = def; baruCost = 0;
            }
            minu = min(minu, uCost);
            minbaru = min(minbaru, baruCost);
        }
        minBarU[i].val = minbaru;
        minU[i].val = minu;
    }

    recomputeTable(curf);

    int minCost = wcsp->getUb();
    for (int j=lb;j<=ub;j++) {
        minCost = min(minCost, curf[n][j].val);
    }

    return minCost;
}

Cost AmongConstraint::minCostOriginal(int var, Value val, bool changed) {    
    Result result = minCost(var, val, changed);
    delete[] result.second;
    return result.first;
}

Cost AmongConstraint::evalOriginal(const String& s)
{
    int n = arity();
    int count = 0;
    for (int i=0;i<n;i++) {
        if (V.find(s[i] - CHAR_FIRST) != V.end()) count++;
    }
    return max(0, max(lb - count, count-ub))*def;
}

void AmongConstraint::recompute() {
    int n = arity();
    minBarU[0].val = minU[0].val = -1;
    for (int i = 1; i <= n; i++) {
        minBarU[i].val = computeMinBarU(i-1);
        minU[i].val = computeMinU(i-1);
    }

    recomputeTable(f, invf);

}

DPGlobalConstraint::Result AmongConstraint::minCost(int var, Value val, bool changed) {

    if (changed) recompute();

    Cost minCost = wcsp->getUb();
    Cost ucost(0), barucost(0);
    if (V.find(val) == V.end()) {
        ucost = def;
    } else {
        barucost = def;
    }
    EnumeratedVariable *x = (EnumeratedVariable*)getVar(var);
    ucost -= deltaCost[var][x->toIndex(val)];
    barucost -= deltaCost[var][x->toIndex(val)];

    minCost = f[var][0].val + barucost + invf[var+1][0].val;
    for (int j=1;j<=ub;j++) {
        Cost tmpMinCost = min(f[var][j].val + barucost + invf[var+1][j].val,
                f[var][j-1].val + ucost + invf[var+1][j].val);
        minCost = min(tmpMinCost, minCost);
    }

    return DPGlobalConstraint::Result(minCost, NULL);
}

void AmongConstraint::recomputeTable(DPTableCell** table, DPTableCell** invTable, int startRow) {
    int n = arity();

    if (startRow == 0)
    {
        for (int j=0;j<=ub;j++) {
            table[0][j].val = j*def;
            table[0][j].source = -1;

            if (invTable != NULL) {
                invTable[n][j].val = (j < lb)?top:0;
                invTable[n][j].source = -1;
            }
        }
        startRow++;
    }

    for (int i=startRow;i<=n;i++) {
        table[i][0].val = table[i-1][0].val + minBarU[i].val;
        table[i][0].source = 0;
        for (int j=1;j<=ub;j++) {
            int choice1 = table[i-1][j].val + minBarU[i].val;
            int choice2 = table[i-1][j-1].val + minU[i].val;
            if (choice1 > choice2) {
                table[i][j].val = choice2;
                table[i][j].source = 2;
            } else {
                table[i][j].val = choice1;
                table[i][j].source = 1;
            }
        }
    }

    if (invTable != NULL) {
        for (int i=n-1;i>=0;i--) {
            for (int j=0;j<ub;j++) {
                invTable[i][j].val = 0;
                int choice1 = invTable[i+1][j].val + minBarU[i+1].val;
                int choice2 = invTable[i+1][j+1].val + minU[i+1].val;
                if (choice1 > choice2) {
                    invTable[i][j].val = choice2;
                    invTable[i][j].source = 2;
                } else {
                    invTable[i][j].val = choice1;
                    invTable[i][j].source = 1;
                }
            }
            invTable[i][ub].val = invTable[i+1][ub].val + minBarU[i+1].val;
            invTable[i][ub].source = 1;
        }
    }

}

Cost AmongConstraint::computeMinU(int var)
{
    Cost minCost = top;
    EnumeratedVariable *x = (EnumeratedVariable*)getVar(var);
    for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
        Cost ucost = 0;
        if (V.find(*v) == V.end()) {
            ucost = def;
        }
        minCost = min(minCost, ucost - deltaCost[var][x->toIndex(*v)]);
    }
    return minCost;
}

Cost AmongConstraint::computeMinBarU(int var)
{
    Cost minCost = top;
    EnumeratedVariable *x = (EnumeratedVariable*)getVar(var);
    for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
        Cost ucost = def;
        if (V.find(*v) == V.end()) {
            ucost = 0;
        }
        minCost = min(minCost, ucost - deltaCost[var][x->toIndex(*v)]);
    }
    return minCost;
} 


/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

