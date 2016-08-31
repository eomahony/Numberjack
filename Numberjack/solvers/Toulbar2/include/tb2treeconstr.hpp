/** \file tb2treeconstr.hpp
 *  \brief Dynamic programming based global constraint : tree
 */

#ifndef TB2TREECONSTR_HPP_
#define TB2TREECONSTR_HPP_

#include "tb2dpglobalconstr.hpp"

#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <string>

using namespace std;

//\brief Special query for Tree global constraint
template <class T>
class RangeMinQuery {

private:
    vector<T> A;

    int n;
    vector<int> pow2array;
    vector<int> log2array;
    vector<vector<int> > M;

public:
    RangeMinQuery(): n(0) {}

    ~RangeMinQuery() {
    }

    T& operator[](int i) {return A[i];}

    void push_back(const T &t){A.push_back(t);}

    int size() {return A.size();}

    void clear() {A.clear();}

    void pre_compute(){
        if (n != (int) A.size())
        {
            pow2array.clear();
            log2array.clear();
            M.clear();

            n = A.size();
            pow2array.resize(n+1);
            log2array.resize(n+1);

            for (int i=0;i<n+1;i++) log2array[i] = -1;
            pow2array[0] = 1;
            log2array[1] = 0;
            for (int i=1;i<n+1;i++) {
                pow2array[i] = pow2array[i-1]*2;
                if (pow2array[i] < n+1) log2array[pow2array[i]] = i;
            }

            int logVal = 0;
            for (int i=1;i<n+1;i++) {
                if (log2array[i] == -1) log2array[i] = logVal;
                else logVal = log2array[i];
            }

            M.resize(n);
            for (int i=0;i<n;i++) {
                M[i].resize(n);
            }
        }

        for (int i=0;i<n;i++) M[i][0] = i;
        for (int j=1;pow2array[j]<=n;j++) {
            for (int i=0;i < n - pow2array[j] + 1;i++) {
                int minL = M[i][j-1];
                int minR = M[i + pow2array[j-1]][j-1];
                if (A[minL] < A[minR]) M[i][j] = minL;
                else M[i][j] = minR;
            }
        }
    }

    int query(int start, int end) {
        int logWidth = log2array[end - start + 1];
        int minL = M[start][logWidth];
        int minR = M[end - pow2array[logWidth] + 1][logWidth];
        return ((A[minL] < A[minR])?minL:minR);
    }

};

class TreeConstraint : public DPGlobalConstraint
{
private:

    int curTreeCost;

    struct Edge {
        int u;
        int v;
        Cost weight;
        Edge(int u, int v, Cost w): u(u), v(v), weight(w) {}
        bool operator< (const Edge &e) const {return weight < e.weight;}
    };

    int minTreeEdgeCost;
    int maxTreeEdgeCost;
    set<pair<int, int> > treeEdge;

    struct CCTreeNode;  // Forward declaration
    vector<CCTreeNode> nodeStore;
    //typedef vector<CCTreeNode>::iterator CCTreeNodePtr;
    typedef CCTreeNode* CCTreeNodePtr;

    struct CCTreeNode {
        int nodeIndex;
        int u;
        int v;
        Cost weight;
        int height;
        CCTreeNodePtr parent;
        CCTreeNodePtr left;
        CCTreeNodePtr right;
        CCTreeNode():nodeIndex(0), u(-1), v(-1), weight(MIN_COST), height(0), parent(NULL), left(NULL), right(NULL) {}
    };

    vector<CCTreeNodePtr> ccTree;
    vector<CCTreeNodePtr> inorder;
    vector<CCTreeNodePtr> inorderNodeHeight;
    vector<int> pos;
    CCTreeNodePtr ccTreeRoot;
    RangeMinQuery<int> RMQ;

    //CCTreeNodePtr PtrNULL() {return nodeStore.end();}
    CCTreeNodePtr PtrNULL() {return NULL;}
    CCTreeNodePtr createNewNode();

    void joinCCTrees(int u, int v, Cost weight);
    CCTreeNodePtr findRoot(CCTreeNodePtr node);
    void InorderTransveral(CCTreeNodePtr root);

    // disjoint data set
    vector<int> p;
    int findParent(int index, vector<int>& p);
    void unionSet(int u, int v, vector<int>& p);

    map<int, int> val2VarIndex;

    int recomputeCurMST();
    int recomputeMST(vector<Edge> &edgeList);

protected:

    Cost minCostOriginal();
    Cost minCostOriginal(int var, Value val, bool changed);
    Result minCost(int var, Value val, bool changed);

    // This is a hard constraint. SNIC and D(G)AC* are equivalent to AC

    void propagateStrongNIC() {
        propagateAC();
    }

    void propagateDAC() {
        if (ToulBar2::LcLevel == LC_DAC) propagateAC();
    }

    // No need to run anything for (weak) ED(G)AC*
    bool isEAC(int var, Value val) {return true;}
    void findFullSupportEAC(int var) {}

public:
    TreeConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity);
    virtual ~TreeConstraint();

    Cost eval(const String& s);

    void read(istream & file) {} //No parameter needed
    void initMemoization();
    string getName(){return "MST";}
    void dump(ostream& os, bool original = true);
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

