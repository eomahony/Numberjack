#include "tb2treeconstr.hpp"
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <stack>

#include <algorithm>

using namespace std;

TreeConstraint::TreeConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity) : DPGlobalConstraint(wcsp, scope, arity),
        curTreeCost(0), minTreeEdgeCost(0), maxTreeEdgeCost(0), ccTreeRoot(NULL)
{
}


TreeConstraint::~TreeConstraint(void) {}

void TreeConstraint::initMemoization() {

    //check if the assignments represent a graph
    // check max value <= max var. index
    //build value to local var index mapping

    int n = arity();

    for (int i=0;i<n;i++) {
        EnumeratedVariable *x = scope[i];
        int varId = -1;
        for (int j=0;j< (int) wcsp->numberOfVariables() && varId == -1;j++) {
            if (getVar(j) == x) varId = j;
        }
        if (varId < 0) {
            cerr << "variable " << x->getName() << " not found" << endl;
            exit(1);
        }
        val2VarIndex[varId] = i;
        /*for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
			if (*it >= maxVal) maxVal = *it;
		}*/			
    }

    for (int i=0;i<n;i++) {
        EnumeratedVariable *x = scope[i];
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
            if (val2VarIndex.find(*it) == val2VarIndex.end()) {
                cerr << "Error invalid MST()" << endl;
                exit(1);
            }
        }
    }

    p.resize(n);

}

Cost TreeConstraint::minCostOriginal() {			
    if (curTreeCost >= wcsp->getUb()) return wcsp->getUb();
    return 0;
}

Cost TreeConstraint::minCostOriginal(int var, Value val, bool changed) {	
    DPGlobalConstraint::Result result = minCost(var, val, changed);
    return result.first;
}

Cost TreeConstraint::eval(const String& s) {

    int n = arity();
    int root = -1;
    int nRoot = 0;

    vector<vector<int> > edgeList;
    edgeList.resize(n);

    for (int i = 0; i < n; i++) {
        int val = s[i] - CHAR_FIRST;
        if (val2VarIndex[val] == i) {
            root = i;
            nRoot++;
        } else {
            edgeList[val2VarIndex[val]].push_back(i);
        }
    }

    if (nRoot != 1) return wcsp->getUb();       

    bool isLoopBack = false;

    stack<int> dfsStack;
    vector<int> visited;
    dfsStack.push(root);

    visited.resize(n);
    fill(visited.begin(), visited.end(), 0);

    while (!dfsStack.empty() && !isLoopBack) {
        int cur = dfsStack.top();
        dfsStack.pop();
        visited[cur] = 1;
        for (vector<int>::iterator it = edgeList[cur].begin(); it != edgeList[cur].end(); it++) {
            if (visited[*it] == 1) {
                isLoopBack = true;
                break;
            }
            dfsStack.push(*it);
        }        
    }        

    if (isLoopBack) return wcsp->getUb();

    bool allvisit = true;
    for (int i = 0; i < n; i++) {
        if (!visited[i]) allvisit = false;
    }

    if (!allvisit) return wcsp->getUb();

    return 0;

}

DPGlobalConstraint::Result TreeConstraint::minCost(int var, Value val, bool changed) {

    if (changed) curTreeCost = recomputeCurMST();

    bool consistent = true;
    if (curTreeCost >= wcsp->getUb()) {
        consistent = false;
    } else if (treeEdge.find(make_pair(var, val2VarIndex[val])) == treeEdge.end() &&
            treeEdge.find(make_pair(val2VarIndex[val], var)) == treeEdge.end())
    {
        EnumeratedVariable *x = scope[var];
        if (x->getCost(val) + curTreeCost - maxTreeEdgeCost > wcsp->getUb()) {
            consistent = false;
        } else if (x->getCost(val) + curTreeCost - minTreeEdgeCost > wcsp->getUb()) {
            int u = var;
            int v = val2VarIndex[val];
            CCTreeNodePtr lca = ((pos[u]<pos[v])?inorder[RMQ.query(pos[u],pos[v])]:inorder[RMQ.query(pos[v],pos[u])]);
            int maxWeight = lca->weight;
            if (x->getCost(val) + curTreeCost - maxWeight > wcsp->getUb()) consistent = false;
        }
    }

    if (consistent)
        return DPGlobalConstraint::Result(0, NULL);
    else
        return DPGlobalConstraint::Result(wcsp->getUb(), NULL);

}


int TreeConstraint::recomputeCurMST() {
    int n = arity();
    vector<Edge> edgeList;
    for (int i=0;i<n;i++) {
        EnumeratedVariable *x = scope[i];
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
            if (i != val2VarIndex[*it]) {
                edgeList.push_back(Edge(val2VarIndex[*it], i, x->getCost(*it)));
            }
        }
    }
    return recomputeMST(edgeList);
}

int TreeConstraint::recomputeMST(vector<TreeConstraint::Edge> &edgeList) {

    int n = arity();
    int treeCost = 0;

    minTreeEdgeCost = INT_MAX;
    maxTreeEdgeCost = 0;
    treeEdge.clear();
    inorder.clear();
    inorderNodeHeight.clear();
    ccTree.clear();
    nodeStore.clear();
    RMQ.clear();

    nodeStore.reserve(3*n);

    ccTreeRoot = PtrNULL();

    for (int i=0;i<n;i++) {
        p[i] = i;
        CCTreeNodePtr leaf = createNewNode();
        leaf->nodeIndex = i;
        ccTree.push_back(leaf);
    }

    sort(edgeList.begin(),edgeList.end());
    for (vector<Edge>::iterator e = edgeList.begin();e != edgeList.end();e++) {
        if (findParent(e->u, p) != findParent(e->v, p)) {
            unionSet(e->u, e->v, p);
            treeCost += e->weight;
            if (minTreeEdgeCost > e->weight) minTreeEdgeCost = e->weight;
            if (maxTreeEdgeCost < e->weight) maxTreeEdgeCost = e->weight;
            treeEdge.insert(make_pair(e->u, e->v));
            joinCCTrees(e->u, e->v, e->weight);
        }
    }

    int root = findParent(0, p);
    for (int i=1;i<n;i++) {
        if (findParent(i, p) != root) {treeCost = wcsp->getUb();}
    }

    if (treeCost < wcsp->getUb()) {
        for (vector<CCTreeNodePtr>::iterator it = ccTree.begin(); it != ccTree.end();it++) {
            if ((*it)->parent == PtrNULL()) ccTreeRoot = *it;
        }
        InorderTransveral(ccTreeRoot);
        pos.resize(ccTree.size());
        for (vector<CCTreeNodePtr>::iterator it = inorder.begin(); it != inorder.end();it++) {
            pos[(*it)->nodeIndex] = it - inorder.begin();
        }
        RMQ.pre_compute();
    }

    return treeCost;
}

int TreeConstraint::findParent(int index, vector<int>& p) {		
    while (index != p[index]) index = p[index];
    return index;
}

void TreeConstraint::unionSet(int u, int v, vector<int>& p) {
    int uRoot = findParent(u, p);
    int vRoot = findParent(v, p);

    p[uRoot] = vRoot;
}

void TreeConstraint::joinCCTrees(int u, int v, Cost weight) {

    CCTreeNodePtr uRoot = findRoot(ccTree[u]);
    CCTreeNodePtr vRoot = findRoot(ccTree[v]);

    CCTreeNodePtr newRootNode = createNewNode();
    newRootNode->nodeIndex = ccTree.size();
    newRootNode->u = u;
    newRootNode->v = v;
    newRootNode->weight = weight;
    newRootNode->height = max(uRoot->height, vRoot->height) + 1;
    newRootNode->left = uRoot;
    newRootNode->right = vRoot;

    uRoot->parent = newRootNode;
    vRoot->parent = newRootNode;

    ccTree.push_back(newRootNode);

}

TreeConstraint::CCTreeNodePtr TreeConstraint::findRoot(TreeConstraint::CCTreeNodePtr node) {

    if (node->parent == PtrNULL()) return node;
    if (node->parent == node) return node;
    CCTreeNodePtr parent = findRoot(node->parent);
    node->parent = parent;
    return parent;
}

void TreeConstraint::InorderTransveral(TreeConstraint::CCTreeNodePtr root) {

    if (root != PtrNULL()) {
        InorderTransveral(root->left);
        inorder.push_back(root);
        RMQ.push_back(-(root->height));
        InorderTransveral(root->right);
    }
}

TreeConstraint::CCTreeNodePtr TreeConstraint::createNewNode() {
    CCTreeNode newNode;
    newNode.parent = newNode.left = newNode.right = PtrNULL();
    nodeStore.push_back(newNode);            
    return &(nodeStore.back());
}

void TreeConstraint::dump(ostream& os, bool original)
{
    if (original) {
        os << arity_;
        for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
    } else {
        os << nonassigned;
        for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
    }
    os << " -1 smstdp" << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

