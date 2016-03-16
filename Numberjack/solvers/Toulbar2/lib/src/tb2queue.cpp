/*
 * ****** Propagation queue with time stamping *******
 */

#include "tb2queue.hpp"
#include "tb2variable.hpp"


void Queue::push(DLink<VariableWithTimeStamp> *elt, Long curTimeStamp) {
    if (elt->content.timeStamp < curTimeStamp) {
        elt->content.timeStamp = curTimeStamp;
        push_back(elt, false);
    }
}

void Queue::push(DLink<VariableWithTimeStamp> *elt, EventType incdec, Long curTimeStamp) {
    elt->content.incdec |= incdec;
    push(elt, curTimeStamp);
}

void Queue::remove(DLink<VariableWithTimeStamp> *elt) {
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    erase(elt, false);
}

Variable* Queue::pop() {
    assert(!empty());
    DLink<VariableWithTimeStamp> *elt = pop_back(false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable* Queue::pop(int *incdec) {
    assert(!empty());
    *incdec = (*rbegin()).incdec;
    return pop();
}



Variable *Queue::pop_min()
{
    assert(!empty());
    iterator iter=begin();
    DLink<VariableWithTimeStamp> *elt = iter.getElt();
    int pos = (*iter).var->getDACOrder();
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->getDACOrder() < pos) {
            elt = iter.getElt();
            pos = (*iter).var->getDACOrder();
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable *Queue::pop_min(int *incdec)
{
    assert(!empty());
    iterator iter=begin();
    DLink<VariableWithTimeStamp> *elt = iter.getElt();
    int pos = (*iter).var->getDACOrder();
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->getDACOrder() < pos) {
            elt = iter.getElt();
            pos = (*iter).var->getDACOrder();
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    *incdec = elt->content.incdec;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable *Queue::pop_max()
{
    assert(!empty());
    iterator iter=begin();
    DLink<VariableWithTimeStamp> *elt = iter.getElt();
    int pos = (*iter).var->getDACOrder();
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->getDACOrder() > pos) {
            elt = iter.getElt();
            pos = (*iter).var->getDACOrder();
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable *Queue::pop_max(int *incdec)
{
    assert(!empty());
    iterator iter=begin();
    DLink<VariableWithTimeStamp> *elt = iter.getElt();
    int pos = (*iter).var->getDACOrder();
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->getDACOrder() > pos) {
            elt = iter.getElt();
            pos = (*iter).var->getDACOrder();
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    *incdec = elt->content.incdec;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable *Queue::pop_first()
{
    assert(!empty());
    iterator iter=begin();
    DLink<VariableWithTimeStamp> *elt = iter.getElt();
    erase(elt, false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}


void Queue::print(ostream& os)
{
    os << "Queue: ";
    iterator iter=begin();
    if(iter != end()) {
        VariableWithTimeStamp vts = iter.getElt()->content;
        os << "<var:" << vts.var->getName() << ",node:" << vts.timeStamp << "> ";
        for (++iter; iter != end(); ++iter) {
            vts = iter.getElt()->content;
            os << "<var:" << vts.var->getName() << ",node:" << vts.timeStamp << "> ";
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

