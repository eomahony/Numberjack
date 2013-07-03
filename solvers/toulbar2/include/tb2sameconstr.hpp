#include "tb2flowbasedconstr.hpp"

class SameConstraint : public FlowBasedGlobalConstraint
{
	private:
		//int def;
		void buildIndex();
		vector<int> group[2];
		pair<int,int> mapto(int varindex, Value val); 		
		//void checkRemoved(Graph &graph, vector<int> &rmv);
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
		//void findProjection(Graph &graph, int varindex);
		//void findFullSupport(int varindex);
		//void findFullSupport2(int index, vector<int> &supportProvide, bool isEAC);

        void print(ostream& os);
        void dump(ostream& os, bool original = true);
};
