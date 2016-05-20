/** \file tb2queue.hpp
 *  \brief Propagation queue with time stamping.
 * 
 */

#ifndef TB2QUEUE_HPP_
#define TB2QUEUE_HPP_

#include "tb2btlist.hpp"


typedef enum {NOTHING_EVENT=0, INCREASE_EVENT=1, DECREASE_EVENT=2} EventType;

struct VariableWithTimeStamp 
{
    Variable *var;
    Long timeStamp;
    int incdec;
};

class Queue : private BTList<VariableWithTimeStamp>
{  
    // make it private because we don't want copy nor assignment
    Queue(const Queue &s);
    Queue& operator=(const Queue &s);

public:
    Queue() : BTList<VariableWithTimeStamp>(NULL) {}

    int getSize() const {return BTList<VariableWithTimeStamp>::getSize();}
    bool empty() const {return BTList<VariableWithTimeStamp>::empty();}

    void clear() {BTList<VariableWithTimeStamp>::clear();}

    void push(DLink<VariableWithTimeStamp> *elt, Long curTimeStamp);   
    void push(DLink<VariableWithTimeStamp> *elt, EventType incdec, Long curTimeStamp);

    void remove(DLink<VariableWithTimeStamp> *elt);

    Variable *pop();
    Variable *pop(int *incdec);
    Variable *pop_min();
    Variable *pop_min(int *incdec);
    Variable *pop_max();
    Variable *pop_max(int *incdec);
    Variable *pop_first();


    void print(ostream& o);
};

#endif /*TB2QUEUE_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

