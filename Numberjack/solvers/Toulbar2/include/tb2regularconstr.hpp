#include "tb2flowbasedconstr.hpp"

struct DFA {
	vector<int> init;
	vector<int> final;
	set<int> symbol;
	vector<pair<int, int> > *transition;
	int nstate;
	DFA() {}
	void setSize(int size) {transition = new vector<pair<int,int> >[size]; nstate = size;}
	int size() {return nstate;}
	void addTransition(int start, int end, int ch) {
		transition[start].push_back(make_pair(ch, end));
	}
	int getNextState(int start, int ch) {
		int next = -1;
		for (vector<pair<int,int> >::iterator i = transition[start].begin();
		i != transition[start].end();i++) {
			if (i->first == ch) next = i->second;
		}
		return next;
	}
	set<int> getSymbolNeed(int start, int end) {
		set<int> require;
		for (vector<pair<int,int> >::iterator i = transition[start].begin();
		i != transition[start].end();i++) {
			if (i->second == end) require.insert(i->first);
		}
		return require;
	}
	void print() {
		cout << "start state : ";
		for (vector<int>::iterator i = init.begin();i != init.end();i++) cout << *i << " ";
		cout << endl;
		for (int s=0;s<nstate;s++) {
			for (vector<pair<int,int> >::iterator i = transition[s].begin();i != transition[s].end();i++) 
				cout << s << " -" << i->first << "-> " << i->second << endl;
		}
		cout << "end state : ";
		for (vector<int>::iterator i = final.begin();i != final.end();i++) cout << *i << " ";
		cout << endl;
	}
};

class RegularConstraint : public FlowBasedGlobalConstraint
{
	private:
		int subdef, insdef, deldef;
		DFA dfa;
		Graph tgraph;
		vector<Cost> fromSource;
		vector<Cost> toSink;
		vector<vector<Cost> > table;
		vector<set<int> > tempdomain;
		vector<vector<vector<pair<int, int> > > > mapedge;
		pair<int,int> mapto(int varindex, Value val) {return make_pair(0,0);}
		Cost constructFlow(Graph &graph) { 
			//return graph.augment(0,graph.size()-1,false).first; 
			computeShortestPath(graph, cost);
			return cost;
		}
		void checkRemoved(Graph &graph, Cost &cost, vector<int> &rmv);
		void augmentStructure(Graph &graph, Cost &cost, int varindex, map<Value, Cost> &delta);
		void findProjection(Graph &graph, Cost cost, int varindex, map<Value, Cost> &delta);
		void buildGraph(Graph &g);
		void buildGraphBasic(Graph &g, bool needRebuildIndex);
		void computeShortestPath(Graph &g, Cost &cost);
		Cost evalOriginal(String s) {return 0;}
	public:
		RegularConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in);
		~RegularConstraint() {}
		string getName() {return "regular constraint";}
		Cost eval(String s);
		void read(istream &file);
		virtual Cost getMinCost () {
			return constructFlow(graph);
		//	return cost;
		}
		//void findFullSupportEAC(int varindex);
		//void fillEAC2(int index);
		//bool isEAC(int index, Value a); 
		//void findFullSupport2(int index, vector<int> &supportProvide, bool isEAC);

        void print(ostream& os);
  //        void dump(ostream& os, bool original = true);
};

