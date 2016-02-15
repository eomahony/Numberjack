/** \file tb2integer.hpp
 *  \brief Unlimited precision integers with basic operations.
 *
 */

#ifndef TB2ENTIERS_HPP_
#define TB2ENTIERS_HPP_

#include <gmp.h>

/// Unlimited precision integers with basic operations.
/// \note relies on GNU GMP library.
struct BigInteger {
    mpz_t integer; ///< the number

    BigInteger() {
        mpz_init(integer);
    }
    /// allows conversion from a simple double
    BigInteger(double d_) {
        mpz_init(integer);
        mpz_set_d(integer, d_);
    }

    BigInteger(const BigInteger &i) {
        mpz_init(integer);
        mpz_set(integer, i.integer);
    }
    ~BigInteger() {
        mpz_clear(integer);
    }

    BigInteger &operator=(const BigInteger &i) {
        mpz_set(integer, i.integer);
        return *this;
    }
    BigInteger &operator+=(const BigInteger &i) {
        mpz_add(integer, integer, i.integer);
        return *this;
    }
    BigInteger &operator-=(const BigInteger &i) {
        mpz_sub(integer, integer, i.integer);
        return *this;
    }
    BigInteger &operator*=(const BigInteger &i) {
        mpz_mul(integer, integer, i.integer);
        return *this;
    }
    BigInteger &operator/=(const BigInteger &i) {
        assert(i.integer != 0);
        mpz_div(integer, integer, i.integer);
        return *this;
    }
    const BigInteger operator-() const {
        BigInteger i;
        mpz_neg(i.integer, integer);
        return i;
    }
    friend const BigInteger operator+(const BigInteger& left,
            const BigInteger& right) {
        BigInteger i;
        mpz_add(i.integer, left.integer, right.integer);
        return i;
    }
    friend const BigInteger operator-(const BigInteger& left,
            const BigInteger& right) {
        BigInteger i;
        mpz_sub(i.integer, left.integer, right.integer);
        return i;
    }
    friend const BigInteger operator*(const BigInteger& left,
            const BigInteger& right) {
        BigInteger i;
        mpz_mul(i.integer, left.integer, right.integer);
        return i;
    }
    friend const BigInteger operator/(const BigInteger& left,
            const BigInteger& right) {
        BigInteger i;
        assert(right != 0);
        mpz_div(i.integer, left.integer, right.integer);
        return i;
    }
    friend bool operator==(const BigInteger& left, const BigInteger& right) {
        return (mpz_cmp(left.integer, right.integer) == 0);
    }
    friend bool operator!=(const BigInteger& left, const BigInteger& right) {
        return (!(mpz_cmp(left.integer, right.integer) == 0));
    }
    friend bool operator<=(const BigInteger& left, const BigInteger& right) {
        return (mpz_cmp(left.integer, right.integer) <= 0);
    }
    friend bool operator>=(const BigInteger& left, const BigInteger& right) {
        return (mpz_cmp(left.integer, right.integer) >= 0);
    }
    friend bool operator<(const BigInteger& left, const BigInteger& right) {
        return (mpz_cmp(left.integer, right.integer) < 0);
    }
    friend bool operator>(const BigInteger& left, const BigInteger& right) {
        return (mpz_cmp(left.integer, right.integer) > 0);
    }

    void print(ostream& os) const {
        char*p = NULL;
        p = mpz_get_str(p, 10, integer);
        if (strlen(p) > 300)
            //if(strlen(p)-1>=6)
        {
            os << p[0] << '.';
            for (int i = 1; i <= 5; i++)
                os << p[i];
            if (strlen(p) - 1 < 10)
                os << "e+0" << strlen(p) - 1;
            else
                os << "e+" << strlen(p) - 1;
        } else
            os << mpz_get_d(integer);//p;
        //os << mpz_get_d(integer);
    }
    friend ostream& operator<<(ostream& os, const BigInteger &i) {
        i.print(os);
        return os;
    }
    friend istream& operator>>(istream& is, BigInteger& i) {
        //
        double p;
        is >> p;
        mpz_set_d(i.integer, p);
        return is;
    }
};

#endif


/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

