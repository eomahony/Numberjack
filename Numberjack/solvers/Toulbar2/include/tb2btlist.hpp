/** \file tb2btlist.hpp
 *  \brief Backtrackable double-linked list.
 * 
 * Convention: 
 * 
 * elements can be inserted at the end of the list only
 * these insertions can be undone in the reverse order of their insertion
 * 
 * elements can be removed in any order
 * these removals can be undone in the reverse order of their removal.
 * 
 */

#ifndef TB2BTLIST_HPP_
#define TB2BTLIST_HPP_

#include "tb2store.hpp"

template <class T>
struct DLink
{
    bool removed;       // true if the corresponding element has been removed
    DLink *next;
    DLink *prev;
    T content;

public: DLink<T>() : removed(true), next(NULL), prev(NULL) {}
};

template <class T>
class BTList
{
    StoreStack<BTList,DLink<T> *> *storeUndo;
    int size;
    DLink<T> *head;
    DLink<T> *last;

public:
    BTList(StoreStack<BTList,DLink<T> *> *s) : storeUndo(s), size(0), head(NULL), last(NULL) {}

    int getSize() const {return size;}
    bool empty() const {return size == 0;}

    // Warning! clear() is not a backtrackable operation
    void clear() {size = 0; head = NULL; last = NULL;}


    bool inBTList(DLink<T> *elt) {
        for(iterator iter = begin(); iter != end(); ++iter) {
            if(elt == iter.getElt()) return !elt->removed;
        }
        return false;
    }


    void push_back(DLink<T> *elt, bool backtrack) {
        assert( !inBTList(elt) );
        size++;
        elt->removed = false;
        if (last != NULL) {
            last->next = elt;
            elt->prev = last;
        } else {
            head = elt;
            elt->prev = NULL;
        }
        last = elt;
        last->next = NULL;
        if (backtrack) storeUndo->store(this, NULL);
    }


    void undoPushBack() {
        assert(last != NULL);
        size--;
        last->removed = true;
        if (last->prev != NULL) {
            last = last->prev;
            last->next->prev = NULL;
            last->next = NULL;
        } else {
            head = NULL;
            last = NULL;
        }
    }

    void erase(DLink<T> *elt, bool backtrack) {
        assert(!elt->removed);
        size--;
        elt->removed = true;
        if (elt->prev != NULL) {
            assert(!elt->prev->removed);
            assert(elt->prev->next == elt);
            elt->prev->next = elt->next;
        } else head = elt->next;
        if (elt->next != NULL) {
            assert(!elt->next->removed);
            assert(elt->next->prev == elt);
            elt->next->prev = elt->prev;
        } else last = elt->prev;
        if (backtrack) {
            storeUndo->store(this, elt->prev);
            storeUndo->store(this, elt);
        }
    }

    void undoErase(DLink<T> *elt, DLink<T> *prev) {
        assert(elt->removed);
        size++;
        elt->removed = false;
        if (prev != NULL) {
            assert(!prev->removed);
            elt->prev = prev;
            elt->next = prev->next;
            if (prev->next != NULL) prev->next->prev = elt;
            else last = elt;
            prev->next = elt;
        } else {
            if (head != NULL) head->prev = elt;
            else last = elt;
            elt->prev = NULL;
            elt->next = head;
            head = elt;
        }
    }

    // deprecated method to be used with erase(..) storing just one element
    //    void undoErase(DLink<T> *elt) {
    //        assert(elt->removed);
    //        size++;
    //        elt->removed = false;
    //        if (elt->prev != NULL) {
    //            assert(!elt->prev->removed);
    //            assert(elt->prev->next == elt->next);
    //            elt->prev->next = elt;
    //        } else head = elt;
    //        if (elt->next != NULL) {
    //            assert(!elt->next->removed);
    //            assert(elt->next->prev == elt->prev);
    //            elt->next->prev = elt;
    //        } else last = elt;
    //    }

    DLink<T> *pop_back(bool backtrack) {
        assert(last != NULL);
        DLink<T> *oldlast = last;
        erase(last, backtrack);
        return oldlast;
    }


    class iterator
    {
        DLink<T> *elt;
    public:
        iterator() { elt = NULL; }
        iterator(DLink<T> *e) : elt(e) {}

        T operator*() const {
            assert(elt != NULL);
            return elt->content;
        }

        DLink<T> *getElt() const {return elt;}

        iterator &operator++() {    // Prefix form
            if (elt != NULL) {
                while (elt->next != NULL && elt->next->removed) {
                    elt = elt->next;
                }
                elt = elt->next;
            }
            assert(elt == NULL || !elt->removed);
            return *this;
        }

        iterator &operator--() {    // Prefix form
            if (elt != NULL) {
                while (elt->prev != NULL && elt->prev->removed) {
                    elt = elt->prev;
                }
                elt = elt->prev;
            }
            assert(elt == NULL || !elt->removed);
            return *this;
        }

        // To see if you're at the end:
        bool operator==(const iterator &iter) const {return elt == iter.elt;}
        bool operator!= (const iterator &iter) const {return elt != iter.elt;}
    };



    iterator begin() {return iterator(head);}
    iterator end() {return iterator(NULL);}    
    iterator rbegin() {return iterator(last);}
    iterator rend() {return end();}

};

typedef BTList<ConstraintLink> ConstraintList;
typedef BTList<Variable *> VariableList;
typedef BTList<Separator *> SeparatorList;

/*
 * For internal use only! Interaction between tb2store and tb2btlist
 * 
 */

template <class T, class V> template <class Q> void StoreStack<T,V>::restore(BTList<Q> **l, DLink<Q> **elt, ptrdiff_t &x)
{
    if (elt[x] == NULL) {
        l[x]->undoPushBack();
    } else {
        assert(l[x] == l[x-1]);
        l[x]->undoErase(elt[x],elt[x-1]);
        x--;
    }
}

#endif /*TB2BTLIST_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

