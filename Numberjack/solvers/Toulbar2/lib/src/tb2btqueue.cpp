/*
 * ****** Propagation backtrackable queue with time stamping *******
 */

#include "tb2btqueue.hpp"
#include "tb2variable.hpp"

void BTQueue::push(DLink<Variable *> *elt) {
    if (!inBTList(elt)) {
        push_back(elt, true);
    }
}

void BTQueue::remove(DLink<Variable *> *elt) {
    if (inBTList(elt)) {
        erase(elt, true);
    }
}

Variable* BTQueue::pop_back() {
    assert(!empty());
    DLink<Variable *> *elt = BTList<Variable *>::pop_back(true);
    return elt->content;
}

Variable *BTQueue::pop_first()
{
    assert(!empty());
    iterator iter=begin();
    DLink<Variable *> *elt = iter.getElt();
    erase(elt, true);
    return elt->content;
}

void BTQueue::print(ostream& os)
{
    os << "Queue: ";
    iterator iter=begin();
    if(iter != end()) {
        Variable* var = iter.getElt()->content;
        os << var->getName() ;
        for (++iter; iter != end(); ++iter) {
            var = iter.getElt()->content;
            os << " " << var->getName();
        }
    }
    os << endl;
}  

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

