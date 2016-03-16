/** \file tb2sameconstr.hpp
 *  \brief Flow based global cost function : ssame_flow
 */

#ifndef TB2SAMECONSTR_HPP_
#define TB2SAMECONSTR_HPP_

#include "tb2flowbasedconstr.hpp"

class SameConstraint : public FlowBasedGlobalConstraint {
private:
    //int def;
    void buildIndex();
    vector<int> group[2];
    int nDistinctDomainValues;
    pair<int,int> mapto(int varindex, Value val);
    //void checkRemoved(Graph &graph, vector<int> &rmv);
    size_t GetGraphAllocatedSize();
    void buildGraph(Graph &g);
    //void getDomainFromGraph(Graph &graph, int varindex, vector<int> &domain);
    //void augmentGraph(Graph &graph, int &cost, int varindex);
public:
    SameConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int
            arity_in);

    ~SameConstraint() {}

    Cost evalOriginal (String s);
    /*void addToGroup(int gp, Variable *var) {
			for (int i=0;i<arity_;i++) {
				if (getVar(i) == var) {
					group[gp][size[gp]] = i;
					size[gp]++;
					break;
				}
			}
		}
		void addToGroupX(Variable *var) {addToGroup(0, var);}
		void addToGroupY(Variable *var) {addToGroup(1, var);}
     */
    string getName() {return "same constraint";}
    void read(istream &file);
    void addVariablesToGroup(EnumeratedVariable* variable, int groupID) {

        for (int j=0;j<arity_;j++) {
            if (getVar(j) == variable) {
                group[groupID].push_back(j);
                break;
            }
        }

    }
    void organizeConfig() {
        for (int g=0;g<2;g++) sort(group[g].begin(), group[g].end());
    }

    void print(ostream& os);
    void dump(ostream& os, bool original = true);
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

