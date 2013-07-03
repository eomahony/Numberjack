/** \file tb2graph.hpp
 *  \brief Multiple-edged Graph using the adjacent list structure for modelling the flow model
 * 
 */

#ifndef TB2GRAPH
#define TB2GRAPH

#include "tb2types.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <assert.h>

#define adj first
#define weight second
#define INF (MAX_COST >> 5)
#define Free(pt) if(pt) {delete[] pt;pt=NULL;} 

using namespace std;

class Graph {
	private:
		// structure for a element in a adjacent list
		struct List_Node {
			int adj;     //the node connecting to  
			Cost weight; // the weight
			Cost cap;    //the capacity 
			int tag;     // the label of the edge
			List_Node () {}
			List_Node (int a, Cost w, Cost c, int t) : adj(a), weight(w), cap(c), tag(t) {}
			bool operator<(const List_Node &a) const {return weight < a.weight;}  
			bool operator== (const List_Node &a) const {
				return adj == a.adj &&
					cap == a.cap &&
					weight == a.weight;
			}
		};

		vector<List_Node> *adjlist;  
		int *p, *counter;
		Cost *d; 
		Cost *potential;
		// the number of node in the graph
		int gsize;
	public:
		// constructor
        Graph() : adjlist(NULL), p(NULL), counter(NULL), d(NULL), potential(NULL), gsize(0) {} 
		// copy constructor
		Graph(const Graph &g);
		// destructor 
		~Graph() {
			//Free(out);
			//Free(low);
			Free(counter);
			Free(p);
			Free(d);
			Free(adjlist);
			Free(potential);
		}
		// set the number of in the graph, node numbered from 0 to n-1
		void setSize(int n);
		// clear all the edges in the graph but node remain
		void clearEdge(); 

		// connect the node from u to v with weight w, capacity c and a tag tag
		// if tag = -1, multiple edge allowed
		// if tag != -1, multiple edge not allowed 
		void addEdge(int u, int v, Cost w, Cost capacity = 1, int tag = -1);
		//void addEdge(int u, int v, Cost w, int c) {addEdge(u, v, w, c, -1);} 
		//void addEdge(int u, int v, Cost w) {addEdge(u, v, w, 1, -1);} 

		// remove the edge from u to v with tag 
		// return true if success
		// if tag = -1, remove the first edge in the adj. list
		bool removeEdge(int u, int v, int tag = -1); 
		//bool removeEdge(int u, int v) {return removeEdge(u, v, -1);}

		// modify the weight of the edge from u to v with tag 
		// return true if success
		// if tag = -1, modify the weight of the first edge in the adj. list
		bool modifyCost(int u, int v, Cost cost, int tag = -1); 
		//bool modifyCost(int u, int v, Cost cost) {return modifyCost(u, v, cost, -1);}

		// increase the weight of the edge from u to v with tag 
		// return true if success
		// if tag = -1, increase the weight of the first edge in the adj. list
		bool increaseCost(int u, int v, Cost cost, int tag = -1); 
		//bool increaseCost(int u, int v, Cost cost) {return increaseCost(u, v, cost, -1);}

		// retur the weight list froom u to v with tag 
		// if tag = -1, return all edges
		vector<Cost> getWeight(int u, int v, int tag = -1);
		//vector<Cost> getWeight(int u, int v) {return getWeight(u,v,-1);}

		// return the deg of the node i
		int getDegree(int index) {return adjlist[index].size();}

		// return the adjacent list of the node i stored as arrays for List_Node	
		vector<List_Node>& operator[] (int index) {return adjlist[index];}

		// return the number of nodes in the graph
		int size() const {return gsize;}

		// add a flow of floe value flowval starting from u to v following the
		// shorteat path
		// if v and u is connected, the flow is also added to the edge (v,u)
		void addFlow(int u, int v, Cost flowval);

		// compute the cost of the maximum flow from s to t
		pair<int, Cost> minCostFlow(int s, int t, Cost maxValue = INF);

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
		void removeNegativeCycles(Cost &cost);	

		void shortest_path(int source);
		void shortest_path(int source, vector<Cost> &pathCost) {
			shortest_path(source);
			pathCost.resize(size());
			for (int i=0;i<size();i++) pathCost[i] = d[i];
		}
		void shortest_path_with_potential(int source);
		void shortest_path_with_potential(int source, vector<Cost> &pathCost) {
			shortest_path_with_potential(source);
			pathCost.resize(size());
			for (int i=0;i<size();i++) pathCost[i] = d[i];	
		}

		// just for checking
		void print();
		void printPath(int s, int t);

		// iterate each in-coming edges
		class iterator;
		friend class iterator;
		class iterator {
			int node;
			vector<List_Node>::iterator next_edge;
			public:
			iterator() {}
			iterator(int _node, vector<List_Node>::iterator _start) :
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
				return next_edge->adj;
			}

			Cost weight() {
				return next_edge->weight;
			}

			Cost capacity() {
				return next_edge->cap;
			}

			int tag() {
				return next_edge->tag;
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

};

#endif
