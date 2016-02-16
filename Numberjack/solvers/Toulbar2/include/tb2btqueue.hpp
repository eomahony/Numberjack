/** \file tb2btqueue.hpp
 *  \brief Backtrackable propagation queue.
 * 
 */

#ifndef TB2BTQUEUE_HPP_
#define TB2BTQUEUE_HPP_

#include "tb2btlist.hpp"

/*
 * A backtrackable queue
 */
class BTQueue : public BTList<Variable *>
{  
    // make it private because we don't want copy nor assignment
    BTQueue(const BTQueue &s);
    BTQueue& operator=(const BTQueue &s);

public:
    BTQueue(StoreStack<BTList<Variable *>, DLink<Variable *> *> *sv) : BTList<Variable *>(sv) {}

    int getSize() const {return BTList<Variable *>::getSize();}
    bool empty() const {return BTList<Variable *>::empty();}

    void clear() {BTList<Variable *>::clear();}

    void push(DLink<Variable *> *elt);   

    void remove(DLink<Variable *> *elt);

    Variable *pop_back();

    Variable *pop_first();

    void print(ostream& o);
};

#endif /*TB2BTQUEUE_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

