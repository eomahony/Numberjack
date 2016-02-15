/*
 * **************** Storable enumerated domain **********************
 */

#include "tb2domain.hpp"

/*
 * Constructors and misc.
 * 
 */

Domain::Domain(Value inf, Value sup, StoreStack<BTList<Value>, DLink<Value> *> *s) : BTList<Value>(s), initSize(sup - inf + 1), distanceToZero(inf)
{
    init(inf, sup);
}

Domain::Domain(Value *d, int dsize, StoreStack<BTList<Value>, DLink<Value> *> *s) : BTList<Value>(s), initSize(max(d,dsize)-min(d,dsize)+1), distanceToZero(min(d,dsize))
{
    assert( dsize >= 1 );
    assert( dsize <= MAX_DOMAIN_SIZE );
    qsort(d, dsize, sizeof(Value), cmpValue);
    init(d[0], d[dsize-1]);
    int i = 0;
    for (iterator iter = begin(); iter != end(); ++iter) {
        if (*iter < d[i]) BTList<Value>::erase(&all[toIndex(*iter)], false);
        else i++;
    }
}

void Domain::init(Value inf, Value sup)
{
    assert( sup - inf + 1 >= 1 );
    assert( sup - inf + 1 <= MAX_DOMAIN_SIZE );
#if defined(WCSPFORMATONLY) && !defined(NUMBERJACK)
    assert(distanceToZero == 0);
#endif    
    all = new DLink<Value>[sup-inf+1];
    for (int idx=0; idx<sup-inf+1; idx++) {
        all[idx].content = idx + inf;
        push_back(&all[idx], false);
    }
}

int cmpValue(const void *v1, const void *v2)
{
    if (*((int *) v1) < *((int *) v2)) return -1;
    else if (*((int *)v1) > *((int *) v2)) return 1;
    else return 0;
}

ostream& operator<<(ostream& os, Domain &l)
{
    os << "{";
    for (Domain::iterator iter = l.begin(); iter != l.end(); ++iter) {
        os << " " << *iter;
    }
    os << " }";
    return os;
}    

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

