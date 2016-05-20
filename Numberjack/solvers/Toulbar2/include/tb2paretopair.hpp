/** \file tb2paretopair.hpp
 *  \brief ParetoPair numbers with basic operations.
 * 
 */

#ifndef TB2PARETOPAIR_HPP_
#define TB2PARETOPAIR_HPP_

struct ParetoPair {
    int p; //
    int q; //

    ParetoPair() : p(0), q(0) {}

    ParetoPair(int p_, int q_) : p(p_), q(q_) {}

    ParetoPair(int p_) : p(p_), q(p_) {
        cerr << "warning! implicit conversion from int to ParetoPair" << endl;
        exit(EXIT_FAILURE);
    }

    double to_double() const {cerr << "to_double not implemented on Paretopair"; exit(EXIT_FAILURE);}

    ParetoPair(const ParetoPair &r) : p(r.p), q(r.q) {}

    ParetoPair &operator=(const ParetoPair &r) {
        p = r.p;
        q = r.q;
        return *this;
    }
    ParetoPair &operator+=(const ParetoPair &r) {
        p = (p +  r.p);
        q = (q + r.q);
        return *this;
    }
    ParetoPair &operator-=(const ParetoPair &r) {
        p = (p  - r.p);
        q = (q  - r.q); 
        return *this;
    }
    ParetoPair &operator*=(const ParetoPair &r) {
        p = (p  * r.p);
        q = (q  * r.q); 
        return *this;
    }
    ParetoPair &operator/=(const ParetoPair &r) {
        p = (p  / r.p);
        q = (q  / r.q); 
        return *this;
    }
    const ParetoPair operator-() const {return ParetoPair(-p,-q);}

    friend const ParetoPair operator+(const ParetoPair& left, const ParetoPair& right) {
        return ParetoPair(left.p + right.p, left.q + right.q);
    }

    friend const ParetoPair operator-(const ParetoPair& left, const ParetoPair& right) {
        return ParetoPair(max(0,left.p - right.p), max(0,left.q - right.q));
    }


    friend const ParetoPair operator*(const ParetoPair& left, const ParetoPair& right) {
        return ParetoPair(left.p * right.p, left.q * right.q);
    }

    friend const ParetoPair operator/(const ParetoPair& left, const ParetoPair& right) {
        return ParetoPair(left.p / right.p, left.q / right.q);
    }

    friend const ParetoPair operator%(const ParetoPair& left, const ParetoPair& right) {
        return ParetoPair(left.p % right.p, left.q % right.q);
    }

    friend bool operator==(const ParetoPair& left, const ParetoPair& right) {
        return (left.p == right.p) & (left.q == right.q);
    }

    friend bool operator!=(const ParetoPair& left, const ParetoPair& right) {
        return (left.p != right.p)  | (left.q != right.q);
    }

    friend bool operator<=(const ParetoPair& left, const ParetoPair& right) {
        return (left.p <= right.p) & (left.q <= right.q);
    }

    friend bool operator>=(const ParetoPair& left, const ParetoPair& right) {
        return (left.p >= right.p) & ( left.q >= right.q);
    }

    friend bool operator<(const ParetoPair& left, const ParetoPair& right) {
        return ((left.p <= right.p) & ( left.q < right.q)) |
                ((left.p < right.p) & ( left.q <= right.q));
    }

    friend bool operator>(const ParetoPair& left, const ParetoPair& right) {
        return ((left.p >= right.p) & ( left.q > right.q)) |
                ((left.p > right.p) & ( left.q >= right.q));
    }

    void print(ostream& os) const { os << '(' << p << ',' << q << ')'; }

    friend ostream& operator<<(ostream& os, const ParetoPair &r) {
        os << '(' << r.p << ',' << r.q << ')';
        return os;
    }

    friend istream& operator>>(istream& is, ParetoPair& r) {
        char c;
        do {is.get(c);} while (c!='(');
        is >> r.p;
        do {is.get(c);} while (c!=',');
        is >> r.q;
        do {is.get(c);} while (c!=')');
        return is;
    }

    //   friend istream& operator>>(istream& is, ParetoPair& r) {
    //  		is >> r.p;
    //  		r.q = r.p;		// READ ONLY INTEGER, NOT PARETOPAIR !!!!!!!!
    //  		return is;
    //  	}
};

const ParetoPair PARETOPAIR_MIN = ParetoPair(0,0);
const ParetoPair PARETOPAIR_1 = ParetoPair(1,1);
const ParetoPair PARETOPAIR_3 = ParetoPair(3,3);
const ParetoPair PARETOPAIR_100 = ParetoPair(100,100);
const ParetoPair PARETOPAIR_MAX = ParetoPair((INT_MAX/2)/3, (INT_MAX/2)/3);

inline double to_double(const ParetoPair r) {cerr << "to_double not implemented on Paretopair"; exit(EXIT_FAILURE);}
inline Long ceil(const ParetoPair r) {exit(EXIT_FAILURE);return 0;}
inline Long floor(const ParetoPair r){exit(EXIT_FAILURE);return 0;}
inline ParetoPair randomCost(ParetoPair min, ParetoPair max) {return ParetoPair(min.p + (myrand() % (max.p - min.p + 1)), min.q + (myrand() % (max.q - min.q + 1)));}
inline ParetoPair string2Cost(char *ptr) {int p=0,q=0; sscanf(ptr, "(%d,%d)", &p, &q); return ParetoPair(p,q);}

inline int cost2log2(int x)
{
    if (x==0) return -1;
    register int l2 = 0;
    x>>=1;
    for (; x != 0; x >>=1)
    {
        ++ l2;
    }
    return (l2);
}

inline int cost2log2glb(const ParetoPair &r) {return  cost2log2(min(r.p, r.q));}
inline int cost2log2gub(const ParetoPair &r) {return  cost2log2(max(r.p, r.q));}

inline ParetoPair MIN(ParetoPair a, ParetoPair b) {if (a <= b) return a; else if (b <= a) return b; else exit(EXIT_FAILURE);}
inline ParetoPair MAX(ParetoPair a, ParetoPair b) {if (a >= b) return a; else if (b >= a) return b; else exit(EXIT_FAILURE);}
inline ParetoPair GLB(ParetoPair a, ParetoPair b) {return ParetoPair(min(a.p,b.p),min(a.q,b.q));}
inline ParetoPair LUB(ParetoPair a, ParetoPair b) {return ParetoPair(max(a.p,b.p),max(a.q,b.q));}
inline bool GLB(ParetoPair *a, ParetoPair b) {if (!(b >= *a)) {*a = GLB(*a,b); return true;} else return false;}
inline bool LUB(ParetoPair *a, ParetoPair b) {if (!(b <= *a)) {*a = LUB(*a,b); return true;} else return false;}
inline bool GLBTEST(ParetoPair a, ParetoPair b) {return (!(b >= a));}
inline bool LUBTEST(ParetoPair a, ParetoPair b) {return (!(b <= a));}
inline bool DACTEST(ParetoPair a, ParetoPair b) {return (a.p==0 && b.p>0) || (a.q==0 && b.q>0);}
inline bool SUPPORTTEST(ParetoPair a, ParetoPair b) {return DACTEST(a,b);}
inline bool SUPPORTTEST(ParetoPair a) {return (a.p==0 || a.q==0);}
inline bool CUT(ParetoPair lb, ParetoPair ub) {return !(lb < ub);}
inline bool CSP(ParetoPair lb, ParetoPair ub) {ParetoPair r = ub - lb; return (r == ParetoPair(1,0)) || (r == ParetoPair(0,1));}
void initCosts(ParetoPair ub);

#endif /*TB2PARETOPAIR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

