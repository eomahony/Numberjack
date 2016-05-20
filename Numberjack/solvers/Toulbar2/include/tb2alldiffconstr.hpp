/** \file tb2alldiffconstr.hpp
 *  \brief AllDifferent constraint with mu_var and mu_dec measure
 *
 */

#ifndef TB2ALLDIFFCONSTR_HPP_
#define TB2ALLDIFFCONSTR_HPP_

#include "tb2flowbasedconstr.hpp"
#include "tb2binconstr.hpp"
#include "tb2vacutils.hpp"

class AllDiffConstraint : public FlowBasedGlobalConstraint
{
private:
    void buildIndex();
    pair<int,int> mapto(int varindex, Value val) {
        return make_pair(varindex+1, mapval[val]);
    }
    //void getDomainFromGraph(Graph &graph, int varindex, vector<int> &domain);
    size_t GetGraphAllocatedSize();
    void buildGraph(Graph &g);
public:
    static const int DECBI = 2;
    static const int DEC = 1;
    static const int VAR = 0;

    string getName() {return "allDifferent";}

    AllDiffConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int
            arity_in);

    ~AllDiffConstraint() {
    }
    Cost evalOriginal (String s);

    void read(istream &file);
    void organizeConfig();

    void decompose();
    //void initStructure() {if (mode != DECBI) FlowBasedGlobalConstraint::init();}
    //void end() {if (mode != DECBI) FlowBasedGlobalConstraint::end();}
    //void findFullSupport2(int index, vector<int> &supportProvide, bool isEAC);

    void print(ostream& os);
    void dump(ostream& os, bool original = true);
};

#endif /*TB2ALLDIFFCONSTR_HPP_*/


/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

