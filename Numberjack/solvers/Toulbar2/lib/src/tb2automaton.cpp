#include "tb2automaton.hpp"

WFA::WFA() : nbStates(0) {}

WFA::WFA(int _nbStates) : nbStates(_nbStates) {}
WFA::WFA(istream &file) {
    unsigned int _nbStates, nbTransitions, nbStatesInit, nbStatesAccept;
    file >> _nbStates;
    nbStates = _nbStates;

    //cout << this << endl;
    //cout << "nbStates   =" << this->getNbStates() << endl;
    file >> nbStatesInit;
    //cout << "nb initSt  = " << nbStatesInit << endl;
    for (unsigned int state = 0 ; state < nbStatesInit ; state++) {
        unsigned int init;
        Cost weight;
        file >> init >> weight;
        //cout << "reading INIT = " << init << " " << weight << endl;
        pair<int, Cost> initSt = make_pair(init,weight);
        //cout << initSt.first << " " << initSt.second << endl;
        initialStates.push_back(initSt);
    }
    file >> nbStatesAccept;
    //cout << "nb acceptSt = " << nbStatesAccept << endl;
    for (unsigned int state = 0 ; state < nbStatesAccept ; state++) {
        unsigned int accept;
        Cost weight;
        file >> accept >> weight;
        //cout << "reading ACCEPT = " << accept << " " << weight << endl;
        pair<int, Cost> acceptSt = make_pair(accept,weight);
        //cout << acceptSt.first << " " << acceptSt.second << endl;
        acceptingStates.push_back(acceptSt);
    }
    file >> nbTransitions;
    for (unsigned int transition = 0 ; transition < nbTransitions ; transition++) {
        unsigned int start, end, symbol;
        Cost weight;
        file >> start >> symbol >> end >> weight;
        //cout << "TRANS " << start << "x" <<  symbol << "-->" << end << " w= " << weight << endl;
        transitions.push_back(new WTransition(start,end,symbol,weight));
    }
}

WFA::WFA(int nbSymbols, string forbiddenPattern, Cost cost) {
    /// Preparing the WFA : nbStates, initialStates, acceptingStates ///
    nbStates = forbiddenPattern.length();
    initialStates.push_back(make_pair(0,0));
    for (unsigned int state = 0 ; state < nbStates ; state++) {
        acceptingStates.push_back(make_pair(state,0));
    }
    /// Computing transition set
    for (unsigned int currentState = 0 ; currentState < nbStates; currentState++) {
        for (unsigned int symbol = 0 ; symbol < (unsigned int) nbSymbols ; symbol++) {
            string res = forbiddenPattern.substr(0,currentState) + ((char) (symbol+48));
            int weight = (res == forbiddenPattern)?cost:0;
            int start = currentState;
            int end = 0;
            for (int receptionState = ((int) min(currentState+1,nbStates-1)) ; receptionState >0 ; receptionState--) {
                int stringStart  = (currentState+1-receptionState);
                int stringLenght = currentState+1 - stringStart;
                string subCurrent = res.substr(stringStart,stringLenght);
                string subTarget  = forbiddenPattern.substr(0,receptionState);
                if (subCurrent==subTarget) {
                    end = receptionState;
                    break;
                }
            }
            transitions.push_back(new WTransition(start,end,symbol,weight));
        }
    }
}

////////////////////////////////////////////////////////////////////////

void 
WFA::display()  {
    cout << "Number of states = " << nbStates << endl;
    cout << "Initial States : " << endl;
    for (list<pair<int,Cost> >::iterator it = initialStates.begin() ; it != initialStates.end() ; it++) {
        pair<int,int> initial = *it;
        cout << initial.first << "(" << initial.second << ")" << endl;
    }
    cout << "Accepting States : " << endl;
    for (list<pair<int,Cost> >::iterator it = acceptingStates.begin() ; it != acceptingStates.end() ; it++) {
        pair<int,int> accepting = *it;
        cout << accepting.first << "(" << accepting.second << ")" << endl;
    }
    cout << "Transition : " << endl;
    for (list<WTransition*>::iterator it = transitions.begin() ; it != transitions.end() ; it++) {
        WTransition* transition = *it;
        transition->display();
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

