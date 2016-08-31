/** \file tb2graph.hpp
 *  \brief Multiple-edged Graph using the adjacent list structure for modelling the flow model
 * 
 */

#ifndef TB2GRAPH
#define TB2GRAPH

#include "tb2types.hpp"
#include "tb2store.hpp"
#include "tb2btlist.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <list>
#include <algorithm>
#include <assert.h>

//#define adj first
//#define weight second
#define INF (MAX_COST >> 5)
#define Free(pt) if(pt) {delete[] pt;pt=NULL;} 
#define NO_TAG (INT_MAX >> 1)

using namespace std;

// BTList Wrapper for easier usage
template <typename T>
class DLinkStore {
private:
    int blkSize;
    StoreInt curEmpty;
    StoreInt curUsingBlkIndex;
    vector<DLink<T>*> blockStore;
public:
    DLinkStore(int blkSize_)
    : blkSize(blkSize_)
    , curEmpty(0)
    , curUsingBlkIndex(0)
    {
        blockStore.push_back(new DLink<T>[blkSize]);
    }
    ~DLinkStore() {
        for (vector<DLink<int>*>::iterator it = blockStore.begin();
                it != blockStore.end(); it++) {
            delete[] *it;
        }
        blockStore.clear();
    }

    DLink<T>* allocate(const T &ele) {
        if (curEmpty >= blkSize) {
            curEmpty = 0;
            curUsingBlkIndex = curUsingBlkIndex + 1;
            if (curUsingBlkIndex >= (int) blockStore.size()) {
                blockStore.push_back(new DLink<T>[blkSize]);
            }
        }
        DLink<int>* container = &(blockStore[curUsingBlkIndex][curEmpty]);
        container->content = ele;
        curEmpty = curEmpty + 1;
        return container;
    }
};

template <typename T>
class BTListWrapper {
private:
    BTList<T> list;
    DLinkStore<T>* dlinkStore;
public:
    typedef typename BTList<T>::iterator iterator;

    BTListWrapper(DLinkStore<int>* dlinkStore)
    : list(&Store::storeIndexList), dlinkStore(dlinkStore) {}

    ~BTListWrapper() {}

    size_t size() {
        return list.getSize();
    }
    bool empty() {
        return list.getSize() == 0;
    }
    void push_back(const T &ele) {
        DLink<T>* container = dlinkStore->allocate(ele);
        list.push_back(container, true);
    }
    void remove(const T &ele) {
        DLink<T>* target = NULL;
        for(typename BTList<T>::iterator it = list.begin();it != list.end()
        && (target == NULL);++it) {
            if (*it == ele) {
                target = it.getElt();
            }
        }
        if (target != NULL) {
            list.erase(target, true);
        }
    }
    void erase(iterator &it) {
        list.erase(it.getElt(), true);
    }
    iterator begin() {return list.begin();}
    iterator end() {return list.end();}
};

class Graph {

private:

    // structure representing an edge
    struct List_Node {
        // entities need to backtrack
        StoreCost weight; // the weight
        StoreCost cap;    //the capacity, if cap = 0, the edge is set to "deleted"
        // entities need not backtrack
        int adj;     //the node connecting to
        int tag;     // the label of the edge
        int rEdgeIndex; // the pointer to the opposite edge
        List_Node (int depth, int a = -1, Cost w = 0, Cost c = 0, int t = NO_TAG, int rIndex = -1)
        : weight(w)
        , cap(c)
        , adj(a), tag(t), rEdgeIndex(rIndex) {}
    };

    // structure representing a node
    typedef int EdgePtr;
    struct Vertex {
        // entities need to backtrack
        vector<BTListWrapper<EdgePtr>*> edgeList;
        BTListWrapper<int> neighbor;

        Vertex (int n, int depth, DLinkStore<int>* dLinkStore)
        : edgeList(n)
        , neighbor(dLinkStore)
        {
            for (int i=0;i<n;i++) edgeList[i] = new BTListWrapper<int>(dLinkStore);
        }

        ~Vertex()
        {
            for (unsigned int i=0;i<edgeList.size();i++) delete edgeList[i];
        }
    };

    // native adjacent lists, storing all edges possibly appears during search
    vector<vector<List_Node*> > adjlist;

    // additional structure for speeding up traveral
    vector<Vertex*> vertexList;

    // potentials
    //vector<StoreCost> potential;

    // pre-allocated temporary structure
    vector<int> p;
    vector<int> counter;
    vector<Cost> d;

    // the number of node in the graph
    int gsize;

    // for backtractable structure
    int depth;
    DLinkStore<int> intDLinkStore;

    // do not allow copy
    Graph(const Graph &g);

public:
    // constructor
    Graph(int n, int depth);

    // destructor
    ~Graph();

    // connect the node from u to v with weight w, capacity c and a tag tag
    // if tag = -1, multiple edge allowed
    // if tag != -1, multiple edge not allowed
    // if addReverse = true, the reverse edge (v,u) with weight = -w and c = 0 will be automatically added
    void addEdge(int u, int v, Cost w, Cost capacity = 1, int tag = NO_TAG,
            bool addReverse = true) {
        addEdgeInternal(u, v, w, capacity, tag, addReverse, -1);
    }
    int addEdgeInternal(int u, int v, Cost w, Cost capacity, int tag ,
            bool addReverse, int index);

    // remove the edge from u to v with tag
    // return true if success
    // if tag = -1, remove the first edge in the adj. list
    bool removeEdge(int u, int v, int tag = NO_TAG);

    // modify the weight of the edge from u to v with tag
    // return true if success
    // if tag = -1, modify the weight of the first edge in the adj. list
    bool modifyCost(int u, int v, Cost cost, int tag = NO_TAG);

    // increase the weight of the edge from u to v with tag
    // return true if success
    // if tag = -1, increase the weight of the first edge in the adj. list
    bool increaseCost(int u, int v, Cost cost, int tag = NO_TAG);

    // retur the weight list froom u to v with tag
    // if tag = -1, return all edges
    vector<Cost> getWeight(int u, int v, int tag = NO_TAG);

    // Check if any edges going from u to v exists
    bool edgeExist(int u, int v);

    // retur the minimum weight froom all edges from u to v with tag
    // if tag = -1, consider all edges from u to v
    Cost getMinWeight(int u, int v, int tag = NO_TAG);

    // return the number of nodes in the graph
    int size() const {return gsize;}

    // add a flow of flow value flowval starting from u to v following the
    // shorteat path
    // if v and u is connected, the flow is also added to the edge (v,u)
    void addFlow(int u, int v, Cost flowval);

    // compute the cost of the maximum flow from s to t
    pair<int, Cost> minCostFlow(int s, int t);

    // add the flow from s to t following the shortest path which capacity is
    // bounded by the capacaity between s and t (if any)
    // if can_change = true, augment the graph to the residual graph
    // return :
    // a pair <cost,exist>, where cost is the cost of the augmenting path,
    // exist returns true if such a path from s to t exists.
    pair<Cost,bool> augment(int s, int t, bool can_change);
    pair<Cost,bool> augment(int s, int t, bool can_change, vector<pair<int,int> > &edges) {
        pair<Cost, bool> result = augment(s, t, can_change);
        if (result.second && !can_change) {
            edges.clear();
            int u = t;
            while (p[u] != u) {
                edges.push_back(make_pair(p[u], u));
                u = p[u];
            }
        }
        return result;
    }

    // remove all negative cycles in the graph by augmentation
    // cost is modified when all negative cycles are removed.
    void removeNegativeCycles(StoreCost &cost);

    // shortest path algorithms (using Bellmanford)
    void shortest_path(list<int> &sources, bool &nevloop);

    void shortest_path(int source) {
        list<int> sources;
        bool nevLoop;

        sources.push_back(source);
        shortest_path(sources, nevLoop);

        if (nevLoop) {
            cout << "negative loop exists from " <<endl;
            print();
            exit(1);
        }
    }

    void shortest_path(int source, vector<Cost> &pathCost) {
        shortest_path(source);
        pathCost.resize(size());
        for (int i=0;i<size();i++) pathCost[i] = d[i];
    }

    // shortest path algorithm (using Dijkstra with reweighting)
    // Not used due to error in computing potentials after argumentation
    /*void shortest_path_with_potential(int source);
		void shortest_path_with_potential(int source, vector<Cost> &pathCost) {
			shortest_path_with_potential(source);
			pathCost.resize(size());
			for (int i=0;i<size();i++) pathCost[i] = d[i];	
		}*/

    // just for checking
    void print() {
        print(cout);
    }
    void print(ostream &os);
    void printPath(int s, int t);

    // iterate each in-coming edges
    class iterator;
    friend class iterator;
    class iterator {
        int node;
        vector<List_Node*>::iterator next_edge;
    public:
        iterator() : node(0) {}
        iterator(int _node, vector<List_Node*>::iterator _start) :
            node(_node), next_edge(_start) {}

        iterator &operator++() {    // Prefix form
            next_edge++;
            return *this;
        }

        iterator &operator--() {    // Prefix form
            next_edge--;
            return *this;
        }

        int adjNode() {
            return (*next_edge)->adj;
        }

        Cost weight() {
            return (*next_edge)->weight;
        }

        Cost capacity() {
            return (*next_edge)->cap;
        }

        int tag() {
            return (*next_edge)->tag;
        }

        // To see if you're at the end:
        bool operator==(const iterator &iter) const {return next_edge == iter.next_edge;}
        bool operator!=(const iterator &iter) const {return next_edge != iter.next_edge;}
    };
    iterator begin(int node) {
        return iterator(node, adjlist[node].begin());
    }
    iterator end(int node) {
        return iterator(node, adjlist[node].end());
    }

    // iterate each in-coming edges between two nodes
    class edge_iterator;
    friend class edge_iterator;
    class edge_iterator {
        BTListWrapper<int>::iterator next_edge;
        vector<List_Node*> &edgeInfo;
    public:
        //edge_iterator() {}
        edge_iterator(BTListWrapper<int>::iterator _start,
                vector<List_Node*> &_edgeInfo)
        : next_edge(_start)
        , edgeInfo(_edgeInfo) {}

        edge_iterator &operator++() {    // Prefix form
            ++next_edge;
            return *this;
        }

        int adjNode() {
            return edgeInfo[*next_edge]->adj;
        }

        Cost weight() {
            return edgeInfo[*next_edge]->weight;
        }

        Cost capacity() {
            return edgeInfo[*next_edge]->cap;
        }

        int tag() {
            return edgeInfo[*next_edge]->tag;
        }

        // To see if you're at the end:
        bool operator==(const edge_iterator &iter) const {return next_edge == iter.next_edge;}
        bool operator!=(const edge_iterator &iter) const {return next_edge != iter.next_edge;}
    };
    edge_iterator begin(int u, int v) {
        return edge_iterator(vertexList[u]->edgeList[v]->begin(), adjlist[u]);
    }
    edge_iterator end(int u, int v) {
        return edge_iterator(vertexList[u]->edgeList[v]->end(), adjlist[u]);
    }

    // iterate each neigbouring nodes
    class node_iterator;
    friend class node_iterator;
    class node_iterator {
        BTListWrapper<int>::iterator next_node;
    public:
        node_iterator(BTListWrapper<int>::iterator _start)
        : next_node(_start)
        {}

        node_iterator &operator++() {    // Prefix form
            ++next_node;
            return *this;
        }

        int operator*() {
            return *next_node;
        }

        // To see if you're at the end:
        bool operator==(const node_iterator &iter) const {return next_node == iter.next_node;}
        bool operator!=(const node_iterator &iter) const {return next_node != iter.next_node;}

    };
    node_iterator node_begin(int node) {
        return node_iterator(vertexList[node]->neighbor.begin());
    }
    node_iterator node_end(int node) {
        return node_iterator(vertexList[node]->neighbor.end());
    }



};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

