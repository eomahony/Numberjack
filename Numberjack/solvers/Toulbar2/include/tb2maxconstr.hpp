/** \file tb2maxconstr.hpp
 *    \brief Dynamic programming based global Max constraint
 *
 */

#ifndef TB2MAXCONSTR_HPP_
#define TB2MAXCONSTR_HPP_

#include "tb2dpglobalconstr.hpp"

#include <fstream>
#include <vector>
#include <map>

using namespace std;

class MaxConstraint : public DPGlobalConstraint {

private:

    struct Entry{
        Cost weight;
        int var, val;
        Entry() : weight(MIN_COST), var(0), val(0) {}
        Entry(int var, int val, Cost weight){
            this->weight = weight;
            this->var = var;
            this->val = val;
        }
        bool operator<(const Entry &y)const{
            return weight < y.weight;
        }
    };

    vector< map<Value, Cost> > weightMap; //weight of each value
    Cost top;

    // Data structure for computing the min. cost
    vector<Entry> sorted;
    vector< vector<int> > stack;
    vector<Cost> cost;
    vector<Cost> best;
    vector<int> last;
    vector<int> query;
    vector<int> link;
    vector<int> tree;
    vector< map<Value, Cost> > mincosts;   

    Cost unary(int var, int val);
    int ancestor(int i);
    void recompute();

    //pick one value out of the n domains with smallest weight, 
    //and choose the largest one
    Cost largest;
    void findLargest();		

protected:

    Cost minCostOriginal();
    Cost minCostOriginal(int var, Value val, bool changed);
    Result minCost(int var, Value val, bool changed);

    Cost evalOriginal( const String& s );

public:
    MaxConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity);    
    virtual ~MaxConstraint();

    void read(istream &file);
    //void setDefaultViolationCost(Cost cost) {if(configuring) def = cost;}                
    //void setViolationMeasure(int measure) {if(configuring) mode = measure;}                                                   
    void setAssignmentWeight(EnumeratedVariable* variable, Value value, Cost weight) {    				

        int varID = -1;
        for (int j=0;j<arity_ && varID == -1;j++) {
            if (getVar(j) == variable) varID = j;            
        }

        if (varID == -1) {
            cout << "Error in reading max\n";
            exit(0);
        }		
        weightMap[varID][value] = weight;    
        top = max(top, weight);

    }
    void initMemoization();

    string getName(){return "max";}
    void dump(ostream& os, bool original = true);
};

#endif /*TB2MAXCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

