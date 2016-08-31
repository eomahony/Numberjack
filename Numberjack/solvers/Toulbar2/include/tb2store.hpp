/** \file tb2store.hpp
 *  \brief Generic storable data management.
 *
 *	\defgroup backtrack Backtrack management
 *  Used by backtrack search methods.
 *  Allows to copy / restore the current state using Store::store and Store::restore methods.
 *  All storable data modifications are trailed into specific stacks.
 *
 *  Trailing stacks are associated to each storable type:
 *  - Store::storeValue for storable domain values ::StoreValue (value supports, etc)
 *  - Store::storeCost for storable costs ::StoreCost (inside cost functions, etc)
 *  - Store::storeDomain for enumerated domains (to manage holes inside domains)
 *  - Store::storeConstraint for backtrackable lists of constraints
 *  - Store::storeVariable for backtrackable lists of variables
 *  - Store::storeSeparator for backtrackable lists of separators (see tree decomposition methods)
 *  - Store::storeBigInteger for very large integers ::StoreBigInteger used in solution counting methods
 *
 *  Memory for each stack is dynamically allocated by part of \f$2^x\f$ with \e x initialized to ::STORE_SIZE and increased when needed.
 *  \note storable data are not trailed at depth 0.
 *  \warning ::StoreInt uses Store::storeValue stack (it assumes Value is encoded as int!).
 *  \warning Current storable data management is not multi-threading safe! (Store is a static virtual class relying on StoreBasic<T> static members)
 */

#ifndef TB2STORE_HPP_
#define TB2STORE_HPP_

#include "tb2types.hpp"

#ifndef NUMBERJACK
#ifdef BOOST
#include <boost/version.hpp>
#if (BOOST_VERSION >= 105600)
#include <boost/type_index.hpp>
#else
#include <typeinfo>
#endif
#else
#include <typeinfo>
#endif
#endif

template<class T> class BTList;
template<class T> class DLink;
class Constraint;
class Variable;
class Separator;
class ConstraintLink;

/*
 * Storable stack
 *
 */
template<class T, class V>
class StoreStack
{
    T **pointers;
    V *content;
    ptrdiff_t index;
    ptrdiff_t indexMax;
    ptrdiff_t base;

    // make it private because we don't want copy nor assignment
    StoreStack(const StoreStack &s);
    StoreStack& operator=(const StoreStack &s);

public:
    StoreStack(int powbckmemory = STORE_SIZE)  {
        if (pow(2., powbckmemory) >= SIZE_MAX) {
            cerr << "command-line initial memory size parameter " << powbckmemory << " power of two too large!" << endl;
            exit(EXIT_FAILURE);
        }
        indexMax = (ptrdiff_t) pow(2., powbckmemory);
        pointers = new T *[indexMax];
        content = new V[indexMax];
        index = 0;
        base = 0;
        if (ToulBar2::verbose > 0) {
            cout << "c " << indexMax * (sizeof(V) + sizeof(T *)) << " Bytes allocated for "
#ifndef NUMBERJACK
#if (BOOST_VERSION >= 105600)
                << boost::typeindex::type_id<T>().pretty_name()
#else
                << typeid(T).name()
#endif
#endif
                 << " stack." << endl;
        }
    }

    ~StoreStack() {
        delete[] pointers;
        delete[] content;
    }

    void realloc() {
        T **newpointers = new T *[indexMax * 2];
        V *newcontent = new V[indexMax * 2];
        if (!newpointers || !newcontent) {
            cerr
#ifndef NUMBERJACK
#if (BOOST_VERSION >= 105600)
                << boost::typeindex::type_id<T>().pretty_name()
#else
                << typeid(T).name()
#endif
#endif
                    << " stack out of memory!" << endl;
            exit(EXIT_FAILURE);
        }
        std::copy(pointers,pointers+indexMax,newpointers);
        std::copy(content,content+indexMax,newcontent);

        delete[] pointers;
        delete[] content;
        pointers = newpointers;
        content = newcontent;
        indexMax *= 2;
        if (ToulBar2::verbose >= 0) {
            cout << "c " << indexMax * (sizeof(V) + sizeof(T *)) << " Bytes allocated for "
#ifndef NUMBERJACK
#if (BOOST_VERSION >= 105600)
                << boost::typeindex::type_id<T>().pretty_name()
#else
                << typeid(T).name()
#endif
#endif
                    << " stack." << endl;
        }
    }

    void store(T *x, V y) {
        if (index > 0) {
            index++;
            if (index >= indexMax) realloc();
            content[index] = y;
            pointers[index] = x;
        }
    }

    void store(T *x) {
        if (index > 0) {
            index++;
            if (index >= indexMax) realloc();
            content[index] = *x;
            pointers[index] = x;
        }
    }

    void store() {
        index++;
        if (index >= indexMax) realloc();
        pointers[index] = (T *) (intptr_t) base;
        base = index;
    }

    //	void restore(int **adr, int *val, ptrdiff_t x) {
    //		*adr[x] = val[x];
    //	}

    void restore(Value **adr, Value *val, ptrdiff_t x) {
        *adr[x] = val[x];
    }

#ifndef INT_COST
    void restore(Cost **adr, Cost *val, ptrdiff_t x) {
        *adr[x] = val[x];
    }
#endif

    void restore(BigInteger **adr, BigInteger *val, ptrdiff_t x) {
        *adr[x] = val[x];
    }
    template<class Q> void restore(BTList<Q> **l, DLink<Q> **elt, ptrdiff_t &x);

    void restore() {
        if (index > 0) { // do nothing if already at depth = 0
            ptrdiff_t x, y;

            x = index + 1;
            y = base;
            while (--x != y) {
                restore(pointers, content, x);
            }

            index = y - 1;
            base = (ptrdiff_t) pointers[y];
        }
    }
};



/*
 * Storable basic types
 */
template<class T>
class StoreBasic
{
    T v;

public:
    StoreBasic(T vv) : v(vv) {
    } ///< \warning allows conversion from T to StoreBasic<T>, which may loose the compiler when mixing T and StoreBasic<T> in the same expression: explicit cast needed e.g. in T::v1 + (T) StoreBasic<T>::v2

    operator T() const {
        return v;
    } ///< allows conversion from StoreBasic to T

    StoreBasic(const StoreBasic &elt) : v(elt.v) {}

    static void store() { mystore.store(); };
    static void restore() { mystore.restore(); };
    
    StoreBasic &operator=(const StoreBasic &elt) { ///< \note assignment has to be backtrackable
        if (&elt != this) {
            mystore.store(&v);
            v = elt.v;
        }
        return *this;
    }

    StoreBasic &operator=(const T vv) {
        mystore.store(&v);
        v = vv;
        return *this;
    }
    StoreBasic &operator+=(const T vv) {
        mystore.store(&v);
        v += vv;
        return *this;
    }
    StoreBasic &operator-=(const T vv) {
        mystore.store(&v);
        v -= vv;
        return *this;
    }

    static StoreStack<T, T> mystore;
};

template<class T>
StoreStack<T, T> StoreBasic<T>::mystore(STORE_SIZE);

typedef StoreBasic<Value> StoreValue;
typedef StoreValue StoreInt;
typedef StoreBasic<Cost> StoreCost;
typedef StoreBasic<BigInteger> StoreBigInteger;

/*
 * Container for all storable stacks
 */
class Store
{
protected:
    virtual ~Store() = 0; // Trick to avoid any instantiation of Store

public:
    static int depth;
    static StoreStack<BTList<Value> , DLink<Value> *> storeDomain;
    static StoreStack<BTList<ConstraintLink> , DLink<ConstraintLink> *> storeConstraint;
    static StoreStack<BTList<Variable *> , DLink<Variable *> *> storeVariable;
    static StoreStack<BTList<Separator *> , DLink<Separator *> *> storeSeparator;

    /// \return the current (backtrack / tree search) depth
    static int getDepth() {
        return depth;
    }

    /// makes a copy of the current state
    static void store() {
        depth++;
        StoreValue::store();
        StoreCost::store();
        StoreBigInteger::store();
        storeDomain.store();
        storeConstraint.store();
        storeVariable.store();
        storeSeparator.store();
    }

    /// restores the current state to the last copy
    static void restore() {
        depth--;
        StoreValue::restore();
        StoreCost::restore();
        StoreBigInteger::restore();
        storeDomain.restore();
        storeConstraint.restore();
        storeVariable.restore();
        storeSeparator.restore();
    }

    /// restore the current state to the copy made at depth \c newDepth
    static void restore(int newDepth) {
        assert(depth >= newDepth);
        while (depth > newDepth) restore();
    }
};

#define storeIndexList storeDomain

#endif /*TB2STORE_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

