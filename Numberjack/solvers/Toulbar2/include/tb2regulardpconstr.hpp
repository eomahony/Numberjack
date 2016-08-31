/** \file tb2regulardpconstr.hpp
 *  \brief Dynamic programming based global cost function : sregular_dp
 */

#ifndef TB2REGULARDPCONSTR_HPP_
#define TB2REGULARDPCONSTR_HPP_

#include "tb2dpglobalconstr.hpp"
#include "tb2automaton.hpp"
#include <vector>
#include <fstream>
#include <string>
using namespace std;

class RegularDPConstraint : public DPGlobalConstraint {
private:

    struct DFA : public WeightedAutomaton {
        vector<int> init;
        vector<int> final;
        vector<int> symbol;
        map<int, int> symbolIndex;
        vector<pair<int, int> > *transition;
        vector<pair<int, int> > *invTransition;
        int nstate;

        DFA() {
        }

        void setNumStates(int size) {
            transition = new vector<pair<int, int> >[size];
            invTransition = new vector<pair<int, int> >[size];
            nstate = size;
        }

        void addInitialState(int begin) {
            init.push_back(begin);
        }

        void addFinalState(int end) {
            final.push_back(end);
        }

        int size() {
            return nstate;
        }

        void addTransition(int start, int ch, int end, int weight) {
            transition[start].push_back(make_pair(ch, end));
            invTransition[end].push_back(make_pair(ch, start));
            symbol.push_back(ch);
        }

        void finalize() {
            sort(symbol.begin(), symbol.end());
            symbol.erase(unique(symbol.begin(), symbol.end()), symbol.end());
            for (vector<int>::iterator i = symbol.begin(); i != symbol.end(); i++) {
                symbolIndex[*i] = i - symbol.begin();
            }
        }

        void dump(ostream& os, bool original) {
            assert(original); //TODO: case original is false
            os << nstate << endl;
            os << init.size();
            for (vector<int>::iterator i = init.begin(); i != init.end(); i++) os << " " << *i;
            os << endl;
            os << final.size();
            for (vector<int>::iterator i = final.begin(); i != final.end(); i++) os << " " << *i;
            os << endl;
            int nbtrans = 0;
            for (int s = 0; s < nstate; s++) nbtrans += transition[s].size();
            os << nbtrans << endl;
            for (int s = 0; s < nstate; s++) {
                for (vector<pair<int, int> >::iterator i = transition[s].begin(); i != transition[s].end(); i++)
                    os << s << " " << i->first << " " << i->second << endl;
            }
        }

        void print() {
            cout << "start state : ";
            for (vector<int>::iterator i = init.begin(); i != init.end(); i++) cout << *i << " ";
            cout << endl;
            for (int s = 0; s < nstate; s++) {
                for (vector<pair<int, int> >::iterator i = transition[s].begin(); i != transition[s].end(); i++)
                    cout << s << " -" << i->first << "-> " << i->second << endl;
            }
            cout << "end state : ";
            for (vector<int>::iterator i = final.begin(); i != final.end(); i++) cout << *i << " ";
            cout << endl;
        }
    };

    DFA dfa;

    template <class Source>
    struct TableCell {
        int val;
        Source source;
    };

    typedef TableCell<pair<int, Value> > DPTableCell;
    DPTableCell **f;
    DPTableCell **curf;
    DPTableCell **invf;

    typedef TableCell<Value> UnaryTableCell;
    UnaryTableCell **u;

    int top;

    template <class T>
    void resizeTable(T** &table, int width, int heigth) {
        table = new T*[width];
        for (int i = 0; i <= arity(); i++) {
            table[i] = new T[heigth];
        }
    }

    template <class T>
    void deleteTable(T** &table) {
        for (int i = 0; i <= arity(); i++) delete[] table[i];
        delete[] table;
        table = NULL;
    }

    void recomputeTable(DPTableCell** table, DPTableCell** invTable = NULL, int startRow = 0);
    void recompute();

    Cost unary(int ch, int var, Value v);

protected:
    Cost minCostOriginal();
    Cost minCostOriginal(int var, Value val, bool changed);
    Result minCost(int var, Value val, bool changed);
    
    void initMemoization();

public:
    RegularDPConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity);
    virtual ~RegularDPConstraint();

    Cost eval(const String& s);
    void read(istream & file);
    WeightedAutomaton* getWeightedAutomaton() {return &dfa;}
    string getName() {
        return "sregulardp";
    }

    void dump(ostream& os, bool original);
    void print(ostream& os);
};

#endif /* TB2REGULARDPCONSTR_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

