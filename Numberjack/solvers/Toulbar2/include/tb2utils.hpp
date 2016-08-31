/** \file tb2utils.hpp
 *  \brief Miscelaneous usefull functions.
 *
 */

#ifndef TB2UTILS_HPP_
#define TB2UTILS_HPP_

// these includes are needed if compiled on new g++ versions (>4.0?)
#include <climits>
#include <cstdlib>
#include <cstring>
#include <libgen.h>
#include <stdint.h>

#ifdef ILOGLUE
#include <ilsolver/ilosolverint.h>
#else
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#endif

#include <limits>
#include <vector>
#include <map>
#include <sstream>
#include <set>
#include <queue>
using namespace std;

#ifdef NDEBUG
#define DEBONLY(x)
#else
#define DEBONLY(x) x
#endif

template<typename T1, typename T2, typename T3>
struct triplet
{
    T1 first;
    T2 second;
    T3 third;
};

template<typename T1, typename T2, typename T3>
triplet<T1,T2,T3> make_triplet(const T1 &m1, const T2 &m2, const T3 &m3)
{
    triplet<T1,T2,T3> ans;
    ans.first = m1;
    ans.second = m2;
    ans.third = m3;
    return ans;
}

// template<class T>
// T abs(T x) {
//     if (x < 0) return -(x);
//     else return x;
// }

// Warning! Already defined in STL
//template<class T>
//T min(T x, T y) {
//    if (x < y) return x;
//    else return y;
//}
//
//template<class T>
//T max(T x, T y) {
//    if (x > y) return x;
//    else return y;
//}

template<class T>
T min(T *array, int size)
{
    assert(size >= 1);
    T res = array[0];
    for (int i=1; i < size; i++) {
        if (array[i] < res) {
            res = array[i];
        }
    }
    return res;
}

template<class T>
T max(T *array, int size)
{
    assert(size >= 1);
    T res = array[0];
    for (int i=1; i < size; i++) {
        if (array[i] > res) {
            res = array[i];
        }
    }
    return res;
}

template <class T>
inline std::string to_string (const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

template <typename T>
void free_all( T & t ) {
    T tmp;
    t.swap( tmp );
}

#include "tb2system.hpp"

// Cormen et al, 1990. pages 152, 158, and 184

template<class T>
int partition(T A[], int p, int r) {
    T x = A[p];
    int i = p - 1;
    int j = r + 1;
    while (true) {
        do {
            j = j - 1;
        } while (A[j] > x);
        do {
            i = i + 1;
        } while (A[i] < x);
        if (i < j) {
            T tmp = A[i];
            A[i] = A[j];
            A[j] = tmp;
        } else return j;
    }
}

template<class T>
int stochastic_partition(T A[], int p, int r) {
    int i = (myrand()%(r-p+1)) + p;
    T tmp = A[p];
    A[p] = A[i];
    A[i] = tmp;
    return partition(A, p, r);
}

template<class T>
T stochastic_selection(T A[], int p, int r, int i) {
    if (p == r) return A[p];
    int q = stochastic_partition(A, p, r);
    int k = q - p + 1;
    if (i <= k) return stochastic_selection(A, p, q, i);
    else return stochastic_selection(A, q+1, r, i-k);
}

#endif /* TB2UTILS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

