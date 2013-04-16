
#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include "SatWrapper.hpp"


Lit Lit_True;
Lit Lit_False;


void addAMOClauses(Lits literals, Lits auxiliary, SatWrapperSolver *solver, EncodingConfiguration::AMOEncoding amo_encoding){
    /*
        'literals' is the vector of literals for which we are constraining AMO to be true.
        'auxiliary' is a vector of aux variables used for encodings such as the ladder encoding.
            If the ladder encoding is used then 'auxiliary' should contain literals.size()-1
            literals, or else 0 in which case auxiliary variables will be created.
    */
    Lits lits;
    if(amo_encoding == EncodingConfiguration::Pairwise){
        for(unsigned int i=0; i<literals.size(); i++){
            for(unsigned int j=i+1; j<literals.size(); j++){
                lits.clear();
                lits.push_back(~(literals[i]));
                lits.push_back(~(literals[j]));
                solver->addClause(lits);
            }
        }
    } else if(amo_encoding == EncodingConfiguration::Ladder){
        if(auxiliary.size() == 0){
            for(unsigned int i=0; i<literals.size()-1; i++){
                Lit l = Lit(solver->create_atom(NULL,0));
                auxiliary.push_back(l);
                if(i > 0){
                    // chain the bounds x<=i -> x<=i+1
                    lits.clear();
                    lits.push_back(~(auxiliary[i-1]));
                    lits.push_back(l);
                    solver->addClause(lits);
                }
            }
        } else if(auxiliary.size() != literals.size() - 1){
            std::cerr << "ERROR: mismatch in size of auxiliary literals to liters when using ladder encoding." << std::endl;
            exit(1);
        }
        for(unsigned int i=0; i<literals.size(); ++i) {
            // channel x=i -> x>i-1
            if(i) {
                lits.clear();
                lits.push_back(~(literals[i]));
                lits.push_back(~(auxiliary[i-1]));
                solver->addClause(lits);
            }

            // channel x=i -> x<=i
            if(i<literals.size()-1) {
                lits.clear();
                lits.push_back(~literals[i]);
                lits.push_back(auxiliary[i]);
                solver->addClause(lits);
            }
        }
    } else {
        std::cerr << "ERROR unknown AMO specified: " << amo_encoding << std::endl;
        exit(1);
    }
}

void addAMOClauses(Lits literals, SatWrapperSolver *solver, EncodingConfiguration::AMOEncoding amo_encoding){
    Lits dummy;
    addAMOClauses(literals, dummy, solver, amo_encoding);
}

/**************************************************************
 ********************      ENCODINGS        *******************
 **************************************************************/

std::ostream& DomainEncoding::display(std::ostream& o) const {
    if(_values) {
        o << "{" << _values[0];
        for(int i=1; i<_size; ++i)
            o << ", " << _values[i];
        o << "}";
    } else o << "[" << _lower << ".." << _upper << "]";
    o << " (" << _size << ")";
    return o;
}

std::ostream& OffsetDomain::display(std::ostream& o) const {
    o << offset << " + ";
    return _dom_ptr->display(o);
}

std::ostream& FactorDomain::display(std::ostream& o) const {
    o << factor << " * ";
    return _dom_ptr->display(o);
}

std::ostream& EqDomain::display(std::ostream& o) const {
    o << value << (spin ? " == " : " != ");
    return _dom_ptr->display(o);
}

std::ostream& LeqDomain::display(std::ostream& o) const {
    o << bound << (spin ? " < " : " >= ");
    return _dom_ptr->display(o);
}

std::ostream& ConstantDomain::display(std::ostream& o) const {
    o << value;
    return o;
}

DomainEncoding::DomainEncoding(SatWrapper_Expression *o) : AbstractDomain(o) {
    _lower = 0;
    _upper = 1;
    _size = 2;
    _values = NULL;

    // not initialised
    _direct_encoding = -1;
    _order_encoding = -1;
}

DomainEncoding::DomainEncoding(SatWrapper_Expression *o, const int nval) : AbstractDomain(o) {
    _lower = 0;
    _upper = nval-1;
    _size = nval;
    _values = NULL;

    // not initialised
    _direct_encoding = -1;
    _order_encoding = -1;
}

DomainEncoding::DomainEncoding(SatWrapper_Expression *o, const int lb, const int ub) : AbstractDomain(o) {
    _lower = lb;
    _upper = ub;
    _size = ub-lb+1;
    _values = NULL;

    // not initialised
    _direct_encoding = -1;
    _order_encoding = -1;
}

DomainEncoding::DomainEncoding(SatWrapper_Expression *o, SatWrapperIntArray& vals) : AbstractDomain(o) {
    _size   = vals.size();
    _values = new int[_size];
    _lower = _upper = vals.get_item(0);
    for(int i=0; i<_size; ++ i) {
        _values[i] = vals.get_item(i);
        if(_upper < _values[i]) _upper = _values[i];
        if(_lower > _values[i]) _lower = _values[i];
    }
}

DomainEncoding::~DomainEncoding() {
    delete [] _values;
}

int DomainEncoding::contain(const int value) const {
    if(_lower > value || _upper < value) return false;
    else if(_size == 2) return value == _lower || value == _upper;
    else if (!_values ) return true;
    int x = get_index_p(value);
    return (_values[x] == value);
}

int DomainEncoding::next(const int value, const int index) const {
    int nxt = value;

//   std::cout << _lower << " < " << value << " < " << _upper << std::endl;
//   std::cout << "size: " << _size << std::endl;

    if(nxt < _upper) {
        if(_size == 2) return _upper;
        if(_values) {
            if(index >= 0) {
                if(index < _size) return _values[index+1];
            } else {
                int x = get_index_p(value);

// 	std::cout << "index: " << x << std::endl;


                if(x < _size) return _values[x+1];
            }
        } else return ++nxt;
    }
    return value;
}

Lit DomainEncoding::less_or_equal(const int value, const int index) const {
    if(!owner->encoding->order) {
        std::cerr << "ERROR call to Domain::less_or_equal when not using the order encoding." << std::endl;
    }
    if(_lower > value) return Lit_False;
    else if(_upper <= value) return Lit_True;
    else if(_size == 2) {
        if(value == _lower) return ~Lit(_direct_encoding);
        if(value == _upper) return Lit(_direct_encoding);
        return Lit_False;
    }
    else if(index >= 0 && index < _size-1) {
        return Lit(_order_encoding+index);
    } else if(!_values) return Lit(_order_encoding+value-_lower);

    // We need to search for the right variable, the index of 'value'
    // if it is in the domain, or the index of the largest element
    // smaller than 'value' otherwise
    return Lit(_order_encoding+get_index_p(value));
}

Lit DomainEncoding::equal(const int value, const int index) const {
    if(!owner->encoding->direct && _size > 2) {
        std::cerr << "ERROR call to Domain::equal when not using the direct encoding." << std::endl;
    }
    if(_lower > value || _upper < value) return Lit_False;
    else if(_size == 2) {
        if(_lower == value) return ~Lit(_direct_encoding);
        if(_upper == value) return Lit(_direct_encoding);
        return Lit_False;
    } else {
        if(index >= 0 && index < _size) {
            return Lit(_direct_encoding+index);
        }
        if(!_values) return Lit(_direct_encoding+value-_lower);

        // We need to search for the right variable,: the index of 'value'
        int x = get_index_p(value);
        if(_values[x] == value) return Lit(_direct_encoding+x);
    }
    return Lit_False;
}

void DomainEncoding::encode(SatWrapperSolver *solver) {

#ifdef _DEBUGWRAP
    std::cout << "encode x" << owner->_ident << " in ";
    display(std::cout);
    std::cout << std::endl;
#endif

    // FIXME: modify special case of _size=2 to handle different encoding configs
    if(_size == 2) _direct_encoding = solver->create_atom(this,SELF);
    else {
        std::vector<Lit> lits;

        if(owner->encoding->direct){
            // create one atom for each value in the domain (direct encoding)
            _direct_encoding = solver->create_atom(this,DIRECT);
            for(int i=1; i<_size; ++i) solver->create_atom(this,DIRECT);

            // at least one value x=1 or x=2 or ...
            for(int i=0; i<_size; ++i) lits.push_back(equal(getval(i),i));
            solver->addClause(lits);

            // If the order encoding is also used, then the chaining between
            // the two representations will enforce AMO.
            if(!owner->encoding->order){
                addAMOClauses(lits, solver, owner->encoding->amo_encoding);
            }
        }

        if(owner->encoding->order){
            // create one atom for each value in the domain but the last (order encoding)
            _order_encoding = solver->create_atom(this,ORDER);
            for(int i=1; i<_size-1; ++i){
                solver->create_atom(this,ORDER);
                // chain the bounds x<=i -> x<=i+1
                lits.clear();
                lits.push_back(~(less_or_equal(getval(i-1),i-1)));
                lits.push_back(less_or_equal(getval(i),i));
                solver->addClause(lits);
            }
        }

        if(owner->encoding->direct && owner->encoding->order){
            // Chain the direct encoding representation with the order encoding
            Lits direct_lits, order_lits;
            for(int i=0; i<_size; i++) direct_lits.push_back(Lit(_direct_encoding + i));
            for(int i=0; i<_size - 1; i++) order_lits.push_back(Lit(_order_encoding + i));
            addAMOClauses(direct_lits, order_lits, solver, EncodingConfiguration::Ladder);
        }

        solver->validate();
    }
}

void DomainEncoding::print_lit(Lit p, const int type) const {
    std::cout << "(x" << owner->_ident << " " ;
    if(type == SELF) {
        std::cout << "== " << (sign(p) ? _lower : _upper);
    } else if(type == DIRECT) {
        std::cout << (sign(p) ? "!= " : "== ") << (_values ? _values[((var(p))-_direct_encoding)] : ((var(p))-_direct_encoding+_lower));
    } else if(type == ORDER) {
        std::cout << (sign(p) ? "> " : "<= ") << (_values ? _values[((var(p))-_order_encoding)] : ((var(p))-_order_encoding+_lower));
    }
    std::cout << ")";
}

OffsetDomain::OffsetDomain(SatWrapper_Expression *o, AbstractDomain *d, const int os) : AbstractDomain(o) {
    _dom_ptr = d;
    offset = os;
}

int OffsetDomain::getval(int idx) const {
    return (_dom_ptr->getval(idx)+offset);
}
int OffsetDomain::getmin() const {
    return (_dom_ptr->getmin()+offset);
}
int OffsetDomain::getmax() const {
    return (_dom_ptr->getmax()+offset);
}

Lit OffsetDomain::less_or_equal(const int value, const int index) const {
    return _dom_ptr->less_or_equal(value-offset,index);
}
Lit OffsetDomain::equal(const int value, const int index) const {
    return _dom_ptr->equal(value-offset,index);
}

FactorDomain::FactorDomain(SatWrapper_Expression *o, AbstractDomain *d, const int f) : AbstractDomain(o) {
    _dom_ptr = d;
    factor = f;
}

int FactorDomain::getval(int idx) const {
    return ((factor < 0 ? _dom_ptr->getval(_dom_ptr->getsize()-idx-1) : _dom_ptr->getval(idx))*factor);
}
int FactorDomain::getmin() const {
    return (factor < 0 ? factor*(_dom_ptr->getmax()) : factor*(_dom_ptr->getmin()));
}
int FactorDomain::getmax() const {
    return (factor < 0 ? factor*(_dom_ptr->getmin()) : factor*(_dom_ptr->getmax()));
}

Lit FactorDomain::less_or_equal(const int value, const int index) const {
    if(factor>0) {

        int quotient = (value/factor);
        if(value < 0 && value%factor) --quotient;

        return _dom_ptr->less_or_equal(quotient,-1);
        //return _dom_ptr->less_or_equal(value/factor,-1);
    } else {

        int quotient = (value/factor)-1;
        if(value < 0 && value%factor) ++quotient;

        return ~(_dom_ptr->less_or_equal(quotient,-1));
        //return ~(_dom_ptr->less_or_equal((value/factor)-1+(value%factor != 0),-1));
    }
}

Lit FactorDomain::equal(const int value, const int index) const {
    if(value%factor) return Lit_False;
    return _dom_ptr->equal(value/factor,-1);
}

EqDomain::EqDomain(SatWrapper_Expression *o, AbstractDomain *d, const int v, const int s) : AbstractDomain(o) {
    _dom_ptr = d;
    value = v;
    spin = s;
}

Lit EqDomain::less_or_equal(const int v, const int index) const {
    // make sense iff v == 0
    if(v < 0) return Lit_False;
    if(v > 0) return Lit_True;
    // <= 0, i.e., false iff the pointed domain is not equal to 'value' (reverse if neg spin)
    if(!_dom_ptr->owner->encoding->direct){
        return (spin ? ~(_dom_ptr->less_or_equal(value,-1)) : _dom_ptr->less_or_equal(value,-1));
    }
    return (spin ? ~(_dom_ptr->equal(value,-1)) : _dom_ptr->equal(value,-1));
}
Lit EqDomain::equal(const int v, const int index) const {
    // make sense iff v == 0 or v == 1
    if(v < 0 || v > 1) return Lit_False;
    // find out if it is in fact an equality or an inequality
    return ((v == spin) ? _dom_ptr->equal(value,-1) : ~(_dom_ptr->equal(value,-1)));
}

LeqDomain::LeqDomain(SatWrapper_Expression *o, AbstractDomain *d, const int b, const int s) : AbstractDomain(o) {
    _dom_ptr = d;
    bound = b;
    spin = s;
}

Lit LeqDomain::less_or_equal(const int v, const int index) const {
    // make sense iff v == 0
    if(v < 0) return Lit_False;
    if(v > 0) return Lit_True;
    // <= 0, i.e., false iff the pointed domain is <= bound (reverse if neg spin)
    return (spin ? ~(_dom_ptr->less_or_equal(bound,-1)) : _dom_ptr->less_or_equal(bound,-1));
}
Lit LeqDomain::equal(const int v, const int index) const {
    // make sense iff v == 0 or v == 1
    if(v < 0 || v > 1) return Lit_False;
    // find out if it is in fact an equality or an inequality
    return ((v == spin) ? _dom_ptr->less_or_equal(bound,-1) : ~(_dom_ptr->less_or_equal(bound,-1)));
}

Lit ConstantDomain::less_or_equal(const int v, const int index) const {
    if(value <= v) return Lit_True;
    else return Lit_False;
}
Lit ConstantDomain::equal(const int v, const int index) const {
    if(value == v) return Lit_True;
    else return Lit_False;
}

bool processClause(std::vector<Lit>& cl_in, std::vector<Lit>& cl_out) {
    cl_out.clear();
    unsigned int i=0;
    while(i<cl_in.size()) {
        if(cl_in[i] == Lit_True) break;
        else if(cl_in[i] != Lit_False) cl_out.push_back(cl_in[i]);
        ++i;
    }
    return i<cl_in.size();
}


bool less_or_equal(int x, int y){
    return x <= y;
}

bool greater_or_equal(int x, int y){
    return x >= y;
}

bool not_equal(int x, int y){
    return x != y;
}

bool equal_abs_first(int x, int y){
    return abs(x) == y;
}

bool equal_abs_second(int x, int y){
    return x == abs(y);
}


// Adds support clauses for comparison_func(X + K, Y)
void supportClauseEncoder(SatWrapper_Expression *X,
                          SatWrapper_Expression *Y,
                          const int K,
                          SatWrapperSolver *solver,
                          EncodingConfiguration *encoding,
                          bool(*comparison_func)(int,int)){
    Lits lits;
    unsigned int i=0, j=0;
    int x, y;
    bool subsumed;

    if(!encoding->direct){
        std::cerr << "ERROR supportClauseEncoder does not support this encoding." << std::endl;
        exit(1);
    }
    
    for(i=0;i<X->getsize();i++){
        x = X->getval(i);
        subsumed = true;
        lits.clear();
        for(j=0;j<Y->getsize();j++){
            y = Y->getval(j);
            if(comparison_func(x + K, y)){
                lits.push_back(Y->equal(y, j));
            } else {
                subsumed = false;
            }
        }
        // If every value in the domain of Y supports X=x, then there is no need
        // to add this clause because it is subsumed by the ALO.
        if(!subsumed){
            lits.push_back(~(X->equal(x, i)));
            solver->addClause(lits);
        }
    }
}


// (X + Y) == Z
void additionEncoder(SatWrapper_Expression *X,
                     SatWrapper_Expression *Y,
                     SatWrapper_Expression *Z,
                     SatWrapperSolver *solver,
                     EncodingConfiguration *encoding) {
    Lits lits;
    int i, j, x, y, z, prev_x, prev_y, prev_z;

    if(encoding->direct) {
        for(i=0; i<X->getsize(); ++i) {
            x = X->getval(i);
            for(j=0; j<Y->getsize(); ++j) {
                y = Y->getval(j);
                
                lits.clear();
                lits.push_back(~(X->equal(x, i)));
                lits.push_back(~(Y->equal(y, j)));
                lits.push_back(Z->equal(x+y));
                solver->addClause(lits);
            }
        }
    } else if(encoding->order){ 
        // The set of values in Z that were used for adding conflicts on Z later.
        std::set<int> used_zvalues;

        for(i=0; i<X->getsize(); ++i) {
            x = X->getval(i);
            for(j=0; j<Y->getsize(); ++j) {
                y = Y->getval(j);
                z = x + y;
                used_zvalues.insert(z);
#ifdef _DEBUGWRAP
                std::cout << "order additionEncoder i:" << i << " x:" << x << " j:" << j << " y:" << y << " z:" << z << std::endl;
#endif
                
                lits.clear();
                lits.push_back(~(X->less_or_equal(x, i)));
                // if(i > 0) lits.push_back(X->less_or_equal(prev_x, i-1));
                lits.push_back(~(Y->less_or_equal(y, j)));
                // if(j > 0) lits.push_back(Y->less_or_equal(prev_y, j-1));
                lits.push_back(Z->less_or_equal(z));
                solver->addClause(lits);

                if(i > 0 || j > 0){
                    // lits.pop_back();
                    lits.clear();
                    if(i > 0) lits.push_back(X->less_or_equal(prev_x, i-1));
                    if(j > 0) lits.push_back(Y->less_or_equal(prev_y, j-1));
                    lits.push_back(~(Z->less_or_equal(prev_z)));
                    solver->addClause(lits);
                }

                prev_y = y;
                prev_z = z;
            }
            prev_x = x;
        }
        // Need to add clauses saying that Z cannot be any of the values which
        // do not occurr as a result of X + Y
        for(i=0; i<Z->getsize(); i++){
            z = Z->getval(i);
            if(used_zvalues.find(z) == used_zvalues.end()){
                lits.clear();
                if(i > 0) lits.push_back(Z->less_or_equal(prev_z, i - 1));
                lits.push_back(~(Z->less_or_equal(z, i)));
                solver->addClause(lits);
            }
            prev_z = z;
        }
    } else {
       std::cerr << "ERROR: additionEncoder not implemented for this encoding, exiting." << std::endl;
       exit(1);
    }
}



// (X * Y) == Z
// This is almost the exact same as additionEncoder; TODO parameterise the operator.
void multiplicationEncoder(SatWrapper_Expression *X,
                     SatWrapper_Expression *Y,
                     SatWrapper_Expression *Z,
                     SatWrapperSolver *solver,
                     EncodingConfiguration *encoding) {
    Lits lits;
    int i, j, x, y, z, prev_x, prev_y, prev_z;

    if(encoding->direct) {
        for(i=0; i<X->getsize(); ++i) {
            x = X->getval(i);
            for(j=0; j<Y->getsize(); ++j) {
                y = Y->getval(j);
                
                lits.clear();
                lits.push_back(~(X->equal(x, i)));
                lits.push_back(~(Y->equal(y, j)));
                lits.push_back(Z->equal(x * y));
                solver->addClause(lits);
            }
        }
    } else if(encoding->order){
        // The set of values in Z that were used for adding conflicts on Z later.
        std::set<int> used_zvalues;

        for(i=0; i<X->getsize(); ++i) {
            x = X->getval(i);
            for(j=0; j<Y->getsize(); ++j) {
                y = Y->getval(j);
                z = x * y;
                used_zvalues.insert(z);
#ifdef _DEBUGWRAP
                std::cout << "order multiplicationEncoder i:" << i << " x:" << x << " j:" << j << " y:" << y << " z:" << z << std::endl;
#endif
                
                lits.clear();
                lits.push_back(~(X->less_or_equal(x, i)));
                lits.push_back(~(Y->less_or_equal(y, j)));
                lits.push_back(Z->less_or_equal(z));
                solver->addClause(lits);

                if(i > 0 || j > 0){
                    lits.clear();
                    if(i > 0) lits.push_back(X->less_or_equal(prev_x, i-1));
                    if(j > 0) lits.push_back(Y->less_or_equal(prev_y, j-1));
                    lits.push_back(~(Z->less_or_equal(prev_z)));
                    solver->addClause(lits);
                }

                prev_y = y;
                prev_z = z;
            }
            prev_x = x;
        }
        // Need to add clauses saying that Z cannot be any of the values which
        // do not occurr as a result of X * Y
        for(i=0; i<Z->getsize(); i++){
            z = Z->getval(i);
            if(used_zvalues.find(z) == used_zvalues.end()){
                lits.clear();
                if(i > 0) lits.push_back(Z->less_or_equal(prev_z, i - 1));
                lits.push_back(~(Z->less_or_equal(z, i)));
                solver->addClause(lits);
            }
            prev_z = z;
        }
    } else {
       std::cerr << "ERROR: multiplicationEncoder not implemented for this encoding, exiting." << std::endl;
       exit(1);
    }
}

// (X != Y)
void disequalityEncoder(SatWrapper_Expression *X,
                        SatWrapper_Expression *Y,
                        SatWrapperSolver *solver,
                        EncodingConfiguration *encoding) {
    Lits lits;
    int i=0, j=0, x, y, prev_x, prev_y;

    if(encoding->direct && encoding->conflict) {
        while( i<X->getsize() && j<Y->getsize() ) {
            x = X->getval(i);
            y = Y->getval(j);
            if(x == y) {
                lits.clear();
                lits.push_back(~(X->equal(x,i)));
                lits.push_back(~(Y->equal(y,j)));
                solver->addClause(lits);
            }
            if( x <= y ) ++i;
            if( y <= x ) ++j;
        }

    } else if(encoding->order) {
        while( i<X->getsize() && j<Y->getsize() ) {
            x = X->getval(i);
            y = Y->getval(j);

#ifdef _DEBUGWRAP
                std::cout << "order disequalityEncoder i:" << i << " x:" << x << " j:" << j << " y:" << y << " Xsize:" << X->getsize() << " Ysize:" << Y->getsize() << std::endl;
#endif

            if(x == y) {
                lits.clear();
                if(i>0) lits.push_back(X->less_or_equal(prev_x, i-1));
                if(j>0) lits.push_back(Y->less_or_equal(prev_y, j-1));
                lits.push_back(~(X->less_or_equal(x,i)));
                lits.push_back(~(Y->less_or_equal(y,j)));
                solver->addClause(lits);
            }
            if( x <= y ){
                prev_x = x;
                i++;
            }
            if( y <= x ){
                prev_y = y;
                j++;
            }
        }
    } else if(encoding->direct && encoding->support) {
        // Encode support clauses (x_1 \/ x_2 \/ !y_3)
        supportClauseEncoder(X, Y, 0, solver, encoding, &not_equal);
        supportClauseEncoder(Y, X, 0, solver, encoding, &not_equal);
    } else {
        std::cerr << "ERROR: disequalityEncoder not implemented for this encoding, exiting." << std::endl;
        exit(1);
    }
}

// (X == Abs(y))
void absoluteEncoder(SatWrapper_Expression *X,
                     SatWrapper_Expression *Y,
                     SatWrapperSolver *solver,
                     EncodingConfiguration *encoding) {
    Lits lits;
    int i=0, j=0, x, y, prev_x, prev_y;

    if(encoding->direct && encoding->conflict) {
        for(i=0; i<X->getsize(); i++){
            x = X->getval(i);
            lits.clear();
            lits.push_back(~(X->equal(x, i)));
            for(j=0; j<Y->getsize(); j++){
                y = Y->getval(j);
                if(x != abs(y)){
                    lits.push_back(~(Y->equal(y, j)));
                    solver->addClause(lits);
                    lits.pop_back();
                }
            }
        }
    } else if(encoding->order) {
        for(i=0; i<X->getsize(); i++){
            x = X->getval(i);
            for(j=0; j<Y->getsize(); j++){
                y = Y->getval(j);
                if(x != abs(y)) {
                    lits.clear();
                    if(i>0) lits.push_back(X->less_or_equal(prev_x, i-1));
                    if(j>0) lits.push_back(Y->less_or_equal(prev_y, j-1));
                    lits.push_back(~(X->less_or_equal(x,i)));
                    lits.push_back(~(Y->less_or_equal(y,j)));
                    solver->addClause(lits);
                }
                prev_y = y;
            }
            prev_x = x;
        }
    } else if(encoding->direct && encoding->support) {
        // Encode support clauses (x_1 \/ x_2 \/ !y_3)
        supportClauseEncoder(X, Y, 0, solver, encoding, &equal_abs_second);
        supportClauseEncoder(Y, X, 0, solver, encoding, &equal_abs_first);
    } else {
        std::cerr << "ERROR: absoluteEncoder not implemented for this encoding, exiting." << std::endl;
        exit(1);
    }
}

// (X == Y)
void equalityEncoder(SatWrapper_Expression *X,
                     SatWrapper_Expression *Y,
                     SatWrapperSolver *solver,
                     EncodingConfiguration *encoding) {
    std::vector<Lit> lits;
    int i=0, j=0, x, y;

    // equalityEncoder has to ignore the specifications of whether to encode
    // the conflicts & supports. Both need to be included to enforce this.

    if(encoding->direct){
        while( i<X->getsize() && j<Y->getsize() ) {
            x = X->getval(i);
            y = Y->getval(j);

            if(x < y) {
                lits.clear();
                lits.push_back(~(X->equal(x,i)));
                solver->addClause(lits);
                i++;
            } else if(y < x) {
                lits.clear();
                lits.push_back(~(Y->equal(y,j)));
                solver->addClause(lits);
                j++;
            } else if(x == y) {
                lits.clear();
                lits.push_back(~(X->equal(x,i)));
                lits.push_back((Y->equal(y,j)));
                solver->addClause(lits);
                
                lits.clear();
                lits.push_back((X->equal(x,i)));
                lits.push_back(~(Y->equal(y,j)));
                solver->addClause(lits);
                i++;
                j++;
            }
        }
        
        // Need to add not equals for any domain values that are after the overlap.
        // Values before the overlap are added in the while loop above.
        while(i<X->getsize()){
            x = X->getval(i);
            lits.clear();
            lits.push_back(~(X->equal(x, i)));
            solver->addClause(lits);
            i++;
        }
        while(j<Y->getsize()){
            y = Y->getval(j);
            lits.clear();
            lits.push_back(~(Y->equal(y, j)));
            solver->addClause(lits);
            j++;
        }
    }
    if(encoding->order){
        while( i<X->getsize() && j<Y->getsize() ) {
            x = X->getval(i);
            y = Y->getval(j);

            if(x < y) {
                lits.clear();
                lits.push_back(~(X->less_or_equal(x,i)));
                if(i > 0){
                    lits.push_back((X->less_or_equal(X->getval(i-1),i-1)));
                }
                if(i < X->getsize()-1){
                    lits.push_back(~(X->less_or_equal(X->getval(i+1),i+1)));
                }
                solver->addClause(lits);
                i++;
            } else if(y < x) {
                lits.clear();
                lits.push_back(~(Y->less_or_equal(y,j)));
                if(j > 0){
                    lits.push_back((Y->less_or_equal(Y->getval(j-1),j-1)));
                }
                if(j < Y->getsize()-1){
                    lits.push_back(~(Y->less_or_equal(Y->getval(j+1),j+1)));
                }
                solver->addClause(lits);
                j++;
            } else if(x == y) {
                lits.clear();
                if(i > 0) lits.push_back(X->less_or_equal(X->getval(i-1), i-1));
                lits.push_back(~(X->less_or_equal(x,i)));
                lits.push_back((Y->less_or_equal(y,j)));
                solver->addClause(lits);
                
                lits.clear();
                lits.push_back((X->less_or_equal(x,i)));
                lits.push_back(~(Y->less_or_equal(y,j)));
                if(j > 0) lits.push_back(Y->less_or_equal(Y->getval(j-1), j-1));
                solver->addClause(lits);
                i++;
                j++;
            }
        }
    }

    if(!encoding->direct && !encoding->order) {
        std::cerr << "ERROR: no encoding representation has been specified for equalityEncoder, exiting." << std::endl;
        exit(1);
    }
}

// (X == Y) <-> Z
void equalityEncoder(SatWrapper_Expression *X,
                     SatWrapper_Expression *Y,
                     SatWrapper_Expression *Z,
                     SatWrapperSolver *solver,
                     EncodingConfiguration *encoding,
                     const bool spin) {
    unsigned int num_clauses = solver->clause_base.size();
    Lit l;

    equalityEncoder(X, Y, solver, encoding);
    if(encoding->direct) l = spin ? Z->equal(0) : ~(Z->equal(0));
    else if(encoding->order) l = spin ? Z->less_or_equal(0) : ~(Z->less_or_equal(0));
    while(num_clauses < solver->clause_base.size())
        solver->clause_base[num_clauses++].push_back(l);

    disequalityEncoder(X, Y, solver, encoding);
    if(encoding->direct) l = spin ? ~(Z->equal(0)) : Z->equal(0);
    else if(encoding->order) l = spin ? ~(Z->less_or_equal(0)) : Z->less_or_equal(0);
    while(num_clauses < solver->clause_base.size())
        solver->clause_base[num_clauses++].push_back(l);
}


// X+K <= Y
void precedenceEncoder(SatWrapper_Expression *X,
                       SatWrapper_Expression *Y,
                       const int K,
                       SatWrapperSolver *solver,
                       EncodingConfiguration *encoding) {
    Lits lits;
    int i=0, j=0, x=0, y=0;

    if(encoding->order){
        while(i<X->getsize() && j<Y->getsize()) {
            if(i<X->getsize()) x = X->getval(i);
            if(j<Y->getsize()) y = Y->getval(j);
            
#ifdef _DEBUGWRAP
            std::cout << "order precedenceEncoder(K:" << K << ") i:" << i << " x:" << x << " j:" << j << " y:" << y << std::endl;
#endif
            
            if((x + K) == y){
                lits.clear();
                lits.push_back(X->less_or_equal(x, i));
                lits.push_back(~(Y->less_or_equal(y, j)));
                solver->addClause(lits);
                i++; j++;
            } else if((x + K) < y){
                i++;
            } else if((x + K) > y){
                if(x == X->getmin()) {
                    lits.clear();
                    lits.push_back(~(Y->less_or_equal(y,j)));
                    solver->addClause(lits);
                    ++j;
                } else if(y == Y->getmax()) {
                    lits.clear();
                    lits.push_back(X->less_or_equal(X->getval(i-1),i-1));
                    solver->addClause(lits);
                    ++i;
                } else {
                    lits.clear();
                    lits.push_back(X->less_or_equal(X->getval(i-1),i-1));
                    lits.push_back(~(Y->less_or_equal(y,j)));
                    solver->addClause(lits);
                    ++j;
                }
            }
        }
        
    } else if (encoding->direct && encoding->conflict) {
        for(i=0; i<X->getsize(); i++){
            x = X->getval(i);
            for(j=0; j<Y->getsize(); j++){
                y = Y->getval(j);
                if((x + K) > y){
                    lits.clear();
                    lits.push_back(~(X->equal(x, i)));
                    lits.push_back(~(Y->equal(y, j)));
                    solver->addClause(lits);
                } else {
                    // Every value of Y after this will satisfy x + k <= y.
                    break;
                }
            }
        }

    } else if(encoding->direct && encoding->support) {
        supportClauseEncoder(X, Y, K, solver, encoding, &less_or_equal);
        supportClauseEncoder(Y, X, -K, solver, encoding, &greater_or_equal);
    } else {
        std::cerr << "ERROR: precedenceEncoder not implemented for this encoding, exiting." << std::endl;
        exit(1);
    }
}

// (X+K <= Y) <-> Z
void precedenceEncoder(SatWrapper_Expression *X,
                       SatWrapper_Expression *Y,
                       SatWrapper_Expression *Z,
                       const int K,
                       SatWrapperSolver *solver,
                       EncodingConfiguration *encoding) {
    unsigned int num_clauses = solver->clause_base.size();
    Lit l;

    precedenceEncoder(X, Y, K, solver, encoding);
    if(encoding->direct) l = Z->equal(0);
    else if(encoding->order) l = Z->less_or_equal(0);
    while(num_clauses < solver->clause_base.size())
        solver->clause_base[num_clauses++].push_back(l);

    precedenceEncoder(Y, X, 1-K, solver, encoding);
    if(encoding->direct) l = ~(Z->equal(0));
    else if(encoding->order) l = ~(Z->less_or_equal(0));
    while(num_clauses < solver->clause_base.size())
        solver->clause_base[num_clauses++].push_back(l);
}


/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

void SatWrapper_Expression::initialise() {
    _ident = -1;
    _solver = NULL;
    domain = NULL;
    encoding = NULL;
}

SatWrapper_Expression::SatWrapper_Expression() {
    initialise();
    domain = new DomainEncoding(this);
}

SatWrapper_Expression::SatWrapper_Expression(const int nval) {
    initialise();
    domain = new DomainEncoding(this,nval);

#ifdef _DEBUGWRAP
    std::cout << "creating variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif

}

SatWrapper_Expression::SatWrapper_Expression(const int lb, const int ub) {
    initialise();
    domain = new DomainEncoding(this,lb,ub);

#ifdef _DEBUGWRAP
    std::cout << "creating variable [" << lb << ".." << ub << "]" << std::endl;
#endif

}

SatWrapper_Expression::SatWrapper_Expression(SatWrapperIntArray& vals) {
    initialise();
    domain = new DomainEncoding(this,vals);

#ifdef _DEBUGWRAP
    std::cout << "creating variable [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

}

SatWrapper_Expression::~SatWrapper_Expression() {
    delete domain;
}

int SatWrapper_Expression::get_size() const {
    int i=0, domsize=0;
    for(i=0; i<getsize(); ++i)
        if(encoding->direct) domsize += (_solver->truth_value(equal(getval(i),i)) != l_False);
        else if(encoding->order) domsize += (_solver->truth_value(less_or_equal(getval(i),i)) != l_False);
    return domsize;
}

int SatWrapper_Expression::next(const int v) const {
    int nxt = v;

    do nxt = domain->next(nxt);
    while(nxt < domain->getmax() &&
          ((encoding->direct && _solver->truth_value(equal(nxt)) == l_False) ||
           (encoding->order && _solver->truth_value(less_or_equal(nxt)) == l_False)));

    return nxt;
}

int SatWrapper_Expression::get_min() const {
    int v;
    for(int i=0; i<getsize(); ++i) {
        v = getval(i);
        if(encoding->direct){
            if(_solver->truth_value(equal(v,i)) != l_False) return v;
        } else if(encoding->order){
            if(i < getsize() - 1 && _solver->truth_value(less_or_equal(v,i)) != l_False) return v;
        }
    }
    // Special case of last literal encoding the domain under the order encoding does not need to be created as it is implied.
    v = getval(getsize()-2);
    if(encoding->order && _solver->truth_value(less_or_equal(v, getsize()-2)) == l_False) return getval(getsize()-1);
    return getmin();
}

int SatWrapper_Expression::get_max() const {
    int v;
    // Special case of last literal encoding the domain under the order encoding does not need to be created as it is implied.
    v = getval(getsize()-2);
    if(encoding->order && _solver->truth_value(less_or_equal(v, getsize()-2)) == l_False) return getval(getsize()-1);

    for(int i=getsize()-1; i>=0; --i){
        v = getval(i);
        if(encoding->direct){
            if(_solver->truth_value(equal(v,i)) != l_False) return v;
        } else if(encoding->order){
            if(i < getsize() - 1 && _solver->truth_value(less_or_equal(v,i)) == l_False) return getval(i+1);
        }
    }
    return getmax();
}

bool SatWrapper_Expression::contain(const int v) const {
    return (_solver->truth_value(equal(v)) != l_False);
}

int SatWrapper_Expression::get_value() const {
    if(_solver->cp_model) {
        return _solver->cp_model[_ident];
    } else
        return get_min();
}

bool SatWrapper_Expression::has_been_added() const {
    return (_solver != NULL);
}

SatWrapper_Expression* SatWrapper_Expression::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, true);
        solver->_variables.push_back(this);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

#ifdef _DEBUGWRAP
        std::cout << "+ add x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
        if(encoding) std::cout << "Expression encoding:" << (*encoding) << std::endl;
#endif

        domain->encode(_solver);
    }

    return this;
}


SatWrapper_add::SatWrapper_add(SatWrapper_Expression *arg1, SatWrapper_Expression *arg2)
    : SatWrapper_binop(arg1, arg2) {
    int _lower = arg1->getmin()+arg2->getmin();
    int _upper = arg1->getmax()+arg2->getmax();
    domain = new DomainEncoding(this, _lower, _upper);
}

SatWrapper_add::SatWrapper_add(SatWrapper_Expression *arg1, const int arg2)
    : SatWrapper_binop(arg1, arg2) {
    domain = new OffsetDomain(this,arg1->domain, arg2);

#ifdef _DEBUGWRAP
    std::cout << "creating offset expression [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

}

SatWrapper_Expression* SatWrapper_add::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);

        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

#ifdef _DEBUGWRAP
        std::cout << "creating add expression x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

        if(top_level) {

            std::cerr << "add predicate at the top level not supported" << std::endl;
            exit(1);

        } else {

            _vars[0] = _vars[0]->add(_solver, false);

            if(_vars[1]) {

                _vars[1] = _vars[1]->add(_solver, false);
                domain->encode(_solver);


#ifdef _DEBUGWRAP
                std::cout << "encode add predicate x" << _ident << " == x"
                          << _vars[0]->_ident << " + x" << _vars[1]->_ident << std::endl;
#endif

                additionEncoder(_vars[0], _vars[1], this, _solver, encoding);

            }
        }

        _solver->validate();
    }

    return this;
}

SatWrapper_add::~SatWrapper_add() {}

SatWrapper_mul::SatWrapper_mul(SatWrapper_Expression *arg1, const int arg2)
    : SatWrapper_binop(arg1, arg2) {

    domain = new FactorDomain(this,arg1->domain, arg2);

}

void get_mul_bounds(int lb1, int ub1, int lb2, int ub2, int *lb, int *ub){
    std::vector<int> bounds;
    bounds.push_back(lb1 * lb2);
    bounds.push_back(lb1 * ub2);
    bounds.push_back(ub1 * lb2);
    bounds.push_back(ub1 * ub2);
    *lb = *(std::min_element(bounds.begin(), bounds.end()));
    *ub = *(std::max_element(bounds.begin(), bounds.end()));
}

void get_mul_bounds(SatWrapper_Expression *arg1, SatWrapper_Expression *arg2, int *lb, int *ub){
    get_mul_bounds(arg1->getmin(), arg1->getmax(), arg2->getmin(), arg2->getmax(), lb, ub);
}


SatWrapper_mul::SatWrapper_mul(SatWrapper_Expression *arg1, SatWrapper_Expression *arg2)
    : SatWrapper_binop(arg1, arg2) {

    // We create a range between the lb and ub of the possible multiplications
    // between arg1 and arg2. Could possibly change this to only list the
    // possible values later.
    int lb, ub;
    get_mul_bounds(arg1, arg2, &lb, &ub);  
    domain = new DomainEncoding(this, lb, ub);
}

SatWrapper_Expression* SatWrapper_mul::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

#ifdef _DEBUGWRAP
        std::cout << "creating mul expression x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

        if(top_level) {

            std::cerr << "mul predicate at the top level not supported" << std::endl;
            exit(1);

        } else {
            _vars[0] = _vars[0]->add(_solver, false);

            if(_vars[1]) {
                _vars[1] = _vars[1]->add(_solver, false);
                domain->encode(_solver);
                multiplicationEncoder(_vars[0], _vars[1], this, solver, encoding);
            }
        }

        _solver->validate();
    }

    return this;
}

SatWrapper_mul::~SatWrapper_mul() {}

SatWrapper_Abs::SatWrapper_Abs(SatWrapper_Expression *arg1)
    : SatWrapper_Expression() {

    // Need to extract just the set of absolute values from arg1's domain.
    std::set<int> values_set;
    for(int i=0; i<arg1->getsize(); i++){
        values_set.insert(abs(arg1->getval(i)));
    }
    SatWrapperIntArray values;
    for(std::set<int>::iterator it=values_set.begin(); it!=values_set.end(); it++){
        values.add(*it);
    }

    _var = arg1;
    domain = new DomainEncoding(this, values);
}

SatWrapper_Expression* SatWrapper_Abs::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

#ifdef _DEBUGWRAP
        std::cout << "creating abs expression x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

        if(top_level) {

            std::cerr << "Abs predicate at the top level not supported" << std::endl;
            exit(1);

        } else {
            domain->encode(solver);
            _var = _var->add(solver, false);

#ifdef _DEBUGWRAP
            std::cout << "encode abs predicate Abs(x" << _ident << ") == x"
                      << _var->_ident << std::endl;
#endif

            absoluteEncoder(this, _var, solver, encoding);
        }

        _solver->validate();
    }

    return this;
}

SatWrapper_Abs::~SatWrapper_Abs() {}

void SatWrapper_AllDiff::addVar(SatWrapper_Expression* v) {
    _vars.add(v);
}

SatWrapper_AllDiff::SatWrapper_AllDiff( SatWrapper_Expression* arg1, SatWrapper_Expression* arg2 )
    : SatWrapper_Expression() {
    addVar(arg1);
    addVar(arg2);
}

SatWrapper_AllDiff::SatWrapper_AllDiff( SatWrapperExpArray& vars )
    : SatWrapper_Expression() {
    _vars = vars;
}

SatWrapper_AllDiff::~SatWrapper_AllDiff() {
    for(unsigned int i=0; i<_clique.size(); ++i)
        delete _clique[i];
}

SatWrapper_Expression* SatWrapper_AllDiff::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

#ifdef _DEBUGWRAP
        std::cout << "creating alldiff expression" << std::endl;
#endif

        if(top_level) {

            int i, j, n=_vars.size();
            for(i=0; i<n; ++i)
                _vars.set_item(i, (_vars.get_item(i))->add(_solver,false));

            SatWrapper_Expression *exp;
            for(i=1; i<n; ++i)
                for(j=0; j<i; ++j) {
                    exp = new SatWrapper_ne(_vars.get_item(j), _vars.get_item(i));
                    _solver->add(exp);
                    _clique.push_back(exp);
                }

        } else {

            std::cerr << "Alldiff not in top level insupported at the moment" << std::endl;
            exit(1);

        }
    }

    return NULL;
}


SatWrapper_Sum::SatWrapper_Sum(SatWrapperExpArray& vars,
                               SatWrapperIntArray& weights,
                               const int offset)
    : SatWrapper_Expression() {
    _offset = offset;
    _vars = vars;
    _weights = weights;
    initialise();
}

SatWrapper_Sum::SatWrapper_Sum(SatWrapper_Expression *arg1,
                               SatWrapper_Expression *arg2,
                               SatWrapperIntArray& w,
                               const int offset)
    : SatWrapper_Expression() {
    _offset = offset;
    _vars.add(arg1);
    _vars.add(arg2);
    _weights = w;
    initialise();
}

SatWrapper_Sum::SatWrapper_Sum(SatWrapper_Expression *arg,
                               SatWrapperIntArray& w,
                               const int offset)
    : SatWrapper_Expression() {
    _self = NULL;
    _offset = offset;
    _vars.add(arg);
    _weights = w;
    initialise();
}

SatWrapper_Sum::SatWrapper_Sum()
    : SatWrapper_Expression() {
    _offset = 0;
}

void SatWrapper_Sum::initialise() {
    if(_vars.size() == 2) {
        int _lower = 0;
        int _upper = 0;

        int weight;
        for(int i = 0; i < _vars.size(); ++i) {
            weight = _weights.get_item(i);

            if( weight > 0 ) {
                _lower += (weight * _vars.get_item(i)->getmin());
                _upper += (weight * _vars.get_item(i)->getmax());
            } else {
                _upper += (weight * _vars.get_item(i)->getmin());
                _lower += (weight * _vars.get_item(i)->getmax());
            }

        }
        _lower += _offset;
        _upper += _offset;

        domain = new DomainEncoding(this,_lower, _upper);
    }

#ifdef _DEBUGWRAP
    std::cout << "Intermediate variable has values: " << getmin() << " to " << getmax() << std::endl;
#endif

}

SatWrapper_Sum::~SatWrapper_Sum() {

#ifdef _DEBUGWRAP
    std::cout << "delete sum" << std::endl;
#endif

    for(unsigned int i=0; i<_subsum.size(); ++i) {
        delete _subsum[i];
    }
}

void SatWrapper_Sum::addVar(SatWrapper_Expression* v) {
    _vars.add(v);
}

void SatWrapper_Sum::addWeight(const int w) {
    _weights.add(w);
}

void SatWrapper_Sum::set_rhs(const int k) {
    _offset = k;
}

SatWrapper_Expression* SatWrapper_Sum::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

#ifdef _DEBUGWRAP
        std::cout << "add sum expression x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

        if(!top_level) {

            if(_vars.size() == 1) {

                _vars.set_item(0, _vars.get_item(0)->add(_solver, false));

            } else {

                int w1, w2;
                SatWrapper_Expression *exp, *exp1, *exp2;
                for(int i=0; i+1<_vars.size(); i+=2) {
                    w1 = _weights.get_item(i);
                    w2 = _weights.get_item(i+1);

                    if(w1 != 1) exp1 = new SatWrapper_mul(_vars.get_item(i), w1);
                    else exp1 = _vars.get_item(i);

                    if(w2 != 1) exp2 = new SatWrapper_mul(_vars.get_item(i+1), w2);
                    else exp2 = _vars.get_item(i+1);

                    exp = new SatWrapper_add(exp1, exp2);
                    exp->add(_solver, false);

                    _vars.add(exp);
                    _weights.add(1);
                    _subsum.push_back(exp);
                }
            }

            domain = new OffsetDomain(this,_vars.get_item(_vars.size()-1)->domain, _offset);

        } else {
            std::cout << "Warning SUM constraint on top level not supported" << std::endl;
        }
    }

    return this;
}


/* Binary operators */

SatWrapper_binop::SatWrapper_binop(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_Expression() {
    _vars[0] = var1;
    _vars[1] = var2;

#ifdef _DEBUGWRAP
    std::cout << "creating binary operator" << std::endl;
#endif

}

SatWrapper_binop::SatWrapper_binop(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_Expression() {
    _vars[0] = var1;
    _vars[1] = NULL;
    _rhs = rhs;

#ifdef _DEBUGWRAP
    std::cout << "creating binary operator" << std::endl;
#endif

}


SatWrapper_binop::~SatWrapper_binop() {

#ifdef _DEBUGWRAP
    std::cout << "delete binary operator" << std::endl;
#endif

}

SatWrapper_or::SatWrapper_or(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_binop(var1,var2) {

#ifdef _DEBUGWRAP
    std::cout << "creating or predicate" << std::endl;
#endif

    domain = new DomainEncoding(this);

}

SatWrapper_or::SatWrapper_or(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_binop(var1,rhs) {

#ifdef _DEBUGWRAP
    std::cout << "creating or predicate" << std::endl;
#endif

    if(rhs) domain = new ConstantDomain(this,1);
    else domain = new OffsetDomain(this,var1->domain,0);
}


SatWrapper_or::~SatWrapper_or() {

#ifdef _DEBUGWRAP
    std::cout << "delete or" << std::endl;
#endif

}

SatWrapper_Expression* SatWrapper_or::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        std::vector<Lit> lits;
        if(top_level) {

            if(_vars[1]) {

                _vars[0] = _vars[0]->add(_solver, false);
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "add or constraint" << std::endl;
#endif

                lits.clear();
                for(int i=0; i<2; ++i)
                    lits.push_back(~(_vars[i]->equal(0)));
                _solver->addClause(lits);

            } else if(_rhs == 0) {

                _vars[0] = _vars[0]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "add or constraint" << std::endl;
#endif

                lits.clear();
                lits.push_back(~(_vars[0]->equal(0)));
                _solver->addClause(lits);
            }

        } else {

            _vars[0] = _vars[0]->add(_solver, false);

            if(_vars[1]) {

                domain->encode(_solver);
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "add or constraint" << std::endl;
#endif

                // x -> y or z
                lits.clear();
                lits.push_back(this->equal(0));
                for(int i=0; i<2; ++i) {
                    lits.push_back(~(_vars[i]->equal(0)));
                }
                _solver->addClause(lits);

                // y -> x  and z -> x
                for(int i=0; i<2; ++i) {
                    lits.clear();
                    lits.push_back(_vars[i]->equal(0));
                    lits.push_back(~(this->equal(0)));
                    _solver->addClause(lits);
                }
            }
        }

        _solver->validate();
    }

    return this;
}



SatWrapper_and::SatWrapper_and(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_binop(var1,var2) {

#ifdef _DEBUGWRAP
    std::cout << "creating or predicate" << std::endl;
#endif

    domain = new DomainEncoding(this);
}

SatWrapper_and::SatWrapper_and(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_binop(var1,rhs) {

#ifdef _DEBUGWRAP
    std::cout << "creating or predicate" << std::endl;
#endif

    if(rhs) domain = new OffsetDomain(this,var1->domain,0);
    else domain = new ConstantDomain(this,0);
}


SatWrapper_and::~SatWrapper_and() {

#ifdef _DEBUGWRAP
    std::cout << "delete or" << std::endl;
#endif

}


SatWrapper_Expression* SatWrapper_and::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        std::vector<Lit> lits;
        if(top_level) {

            if(_vars[1]) {

                ///////////// OLD BUGGY?? VERSION ///////////////
                // 	_vars[0] = _vars[0]->add(_solver, false);
                // 	_vars[1] = _vars[1]->add(_solver, false);

                // #ifdef _DEBUGWRAP
                // 	std::cout << "add and constraint" << std::endl;
                // #endif

                // 	for(int i=0; i<2; ++i) {
                // 	  lits.clear();
                // 	  lits.push_back(~(_vars[i]->equal(0)));
                // 	  _solver->addClause(lits);
                // 	}

                _vars[0] = _vars[0]->add(_solver, true);
                _vars[1] = _vars[1]->add(_solver, true);

#ifdef _DEBUGWRAP
                std::cout << "add and constraint" << std::endl;
#endif


            } else if(_rhs != 0) {

                ///////////// NEW VERSION ///////////////
                _vars[0] = _vars[0]->add(_solver, true);

#ifdef _DEBUGWRAP
                std::cout << "add and constraint" << std::endl;
#endif

                ///////////// OLD BUGGY?? VERSION ///////////////
                // 	_vars[0] = _vars[0]->add(_solver, false);

                // #ifdef _DEBUGWRAP
                // 	std::cout << "add and constraint" << std::endl;
                // #endif

                // 	lits.clear();
                // 	lits.push_back(~(_vars[0]->equal(0)));
                // 	_solver->addClause(lits);
            }

        } else {

            if(_vars[1]) {

                domain->encode(_solver);

                _vars[0] = _vars[0]->add(_solver, false);
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "add and constraint" << std::endl;
#endif

                // y and z -> x
                // x or ~y or ~z
                lits.clear();
                lits.push_back(~(this->equal(0)));
                for(int i=0; i<2; ++i) {
                    lits.push_back(_vars[i]->equal(0));
                }
                _solver->addClause(lits);

                // x -> y  and x -> z
                for(int i=0; i<2; ++i) {
                    lits.clear();
                    lits.push_back(this->equal(0));
                    lits.push_back(~(_vars[i]->equal(0)));
                    _solver->addClause(lits);
                }

            }
        }

        _solver->validate();
    }

    return this;
}


SatWrapper_eq::SatWrapper_eq(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_binop(var1,var2) {
    domain = new DomainEncoding(this);
}

SatWrapper_eq::SatWrapper_eq(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_binop(var1,rhs) {
    domain = new EqDomain(this,var1->domain, rhs, 1);
}


SatWrapper_eq::~SatWrapper_eq() {

#ifdef _DEBUGWRAP
    std::cout << "delete eq" << std::endl;
#endif

}

SatWrapper_Expression* SatWrapper_eq::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

#ifdef _DEBUGWRAP
        std::cout << "creating eq expression x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

        Lits lits;

        _vars[0] = _vars[0]->add(_solver, false);
        if(top_level) {
            if(_vars[1]) {
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "encode equality constraint x"
                          << _vars[0]->_ident << " == x" << _vars[1]->_ident << std::endl;
#endif

                equalityEncoder(_vars[0], _vars[1], _solver, encoding);

            } else {
#ifdef _DEBUGWRAP
                std::cout << "encode unary equality constraint x"
                          << _vars[0]->_ident << " == " << _rhs << std::endl;
#endif

                if(encoding->direct){
                    lits.clear();
                    lits.push_back(_vars[0]->equal(_rhs));
                    _solver->addClause(lits);
                }
                if(encoding->order){
                    lits.clear();
                    lits.push_back(_vars[0]->less_or_equal(_rhs));
                    _solver->addClause(lits);

                    lits.clear();
                    lits.push_back(_vars[0]->greater_than(_rhs - 1));
                    _solver->addClause(lits);
                }
            }
        } else {
            if(_vars[1]) {
                domain->encode(_solver);
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "encode equality predicate x" << this->_ident << " == (x"
                          << _vars[0]->_ident << " == x" << _vars[1]->_ident << ")" << std::endl;
#endif

                equalityEncoder(_vars[0], _vars[1], this, _solver, encoding, true);
            }
        }
        _solver->validate();
    }

    return this;
}


/* Disequality operator */

SatWrapper_ne::SatWrapper_ne(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_binop(var1,var2) {
    domain = new DomainEncoding(this);
}

SatWrapper_ne::SatWrapper_ne(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_binop(var1,rhs) {
    domain = new EqDomain(this,var1->domain, rhs, 0);
}

SatWrapper_ne::~SatWrapper_ne() {

#ifdef _DEBUGWRAP
    std::cout << "delete notequal" << std::endl;
#endif

}

SatWrapper_Expression* SatWrapper_ne::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        std::vector<Lit> lits;

        _vars[0] = _vars[0]->add(_solver, false);
        if(top_level) {

            if(_vars[1]) {
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "add notequal constraint" << std::endl;
#endif

                disequalityEncoder(_vars[0], _vars[1], _solver, encoding);
            } else  {
#ifdef _DEBUGWRAP
                std::cout << "add notequal" << std::endl;
#endif
                if(encoding->direct){
                    if(encoding->conflict){
                        lits.clear();
                        lits.push_back(~(_vars[0]->equal(_rhs)));
                        _solver->addClause(lits);
                    } else if(encoding->support){
                        SatWrapper_ConstantInt *RHS = new SatWrapper_ConstantInt(_rhs);
                        RHS->encoding = encoding;
                        supportClauseEncoder(RHS, _vars[0], 0, solver, encoding, &not_equal);
                        delete RHS;
                    } else {
                        std::cerr << "ERROR SatWrapper_ne constant not implemented for this encoding yet." << std::endl;
                        exit(1);
                    }
                } else if(encoding->order){
                    lits.clear();
                    lits.push_back(_vars[0]->less_or_equal(_rhs-1));
                    lits.push_back(~(_vars[0]->less_or_equal(_rhs)));
                    _solver->addClause(lits);
                } else {
                    std::cerr << "ERROR SatWrapper_ne constant not implemented for this encoding yet." << std::endl;
                    exit(1);
                }
            }

        } else {
            if(_vars[1]) {
                domain->encode(_solver);
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "add notequal predicate" << std::endl;
#endif
                equalityEncoder(_vars[0], _vars[1], this, _solver, encoding, false);
            }
        }
        _solver->validate();
    }

    return this;
}


/* Leq operator */

SatWrapper_le::SatWrapper_le(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "Creating Le constraint" << std::endl;
#endif

    domain = new DomainEncoding(this);

}

SatWrapper_le::SatWrapper_le(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_binop(var1,rhs) {

#ifdef _DEBUGWRAP
    std::cout << "Leq on constant constructor" << std::endl;
#endif

    domain = new LeqDomain(this,var1->domain, rhs, 1);
}


SatWrapper_le::~SatWrapper_le() {

#ifdef _DEBUGWRAP
    std::cout << "delete lessequal" << std::endl;
#endif

}

SatWrapper_Expression* SatWrapper_le::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        std::vector<Lit> lits;

        _vars[0] = _vars[0]->add(_solver, false);
        if(top_level) {

#ifdef _DEBUGWRAP
                std::cout << "Adding leq constraint" << std::endl;
#endif

            if(_vars[1]) {
                _vars[1] = _vars[1]->add(_solver, false);

                precedenceEncoder(_vars[0], _vars[1], 0, _solver, encoding);
            } else {
                _vars[0] = _vars[0]->add(_solver, false);

                // x0 <= r
                if(_rhs < _vars[0]->getmin()) {
                    std::cerr << "model is inconsistent, -exiting" <<std::endl;
                    exit(1);
                } else if(_rhs < _vars[0]->getmax()) {
                    if(encoding->order){
                        lits.clear();
                        lits.push_back(_vars[0]->less_or_equal(_rhs));
                        _solver->addClause(lits);
                    } else if(encoding->direct){
                        SatWrapper_ConstantInt *RHS = new SatWrapper_ConstantInt(_rhs);
                        RHS->encoding = encoding;
                        precedenceEncoder(_vars[0], RHS, 0, _solver, encoding);
                        delete RHS;
                    } else {
                        std::cerr << "ERROR: SatWrapper_le not implemented for this encoding, exiting." << std::endl;
                        exit(1);
                    }
                }
            }
        } else {

#ifdef _DEBUGWRAP
            std::cout << "Adding leq predicate" << std::endl;
#endif

            if(_vars[1]) {
                domain->encode(_solver);
                _vars[1] = _vars[1]->add(_solver, false);

                precedenceEncoder(_vars[0], _vars[1], this, 0, _solver, encoding);
            }
        }

        _solver->validate();
    }

    return this;
}


/* Geq operator */

SatWrapper_ge::SatWrapper_ge(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "Creating Ge constraint" << std::endl;
#endif

    domain = new DomainEncoding(this);
}

SatWrapper_ge::SatWrapper_ge(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_binop(var1,rhs) {

#ifdef _DEBUGWRAP
    std::cout << "Geq on constant constructor" << std::endl;
#endif

    domain = new LeqDomain(this,var1->domain, rhs-1, 0);
}

SatWrapper_ge::~SatWrapper_ge() {

#ifdef _DEBUGWRAP
    std::cout << "delete greaterequal" << std::endl;
#endif

}

SatWrapper_Expression* SatWrapper_ge::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        std::vector<Lit> lits;

        _vars[0] = _vars[0]->add(_solver, false);
        if(top_level) {

            if(_vars[1]) {

                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "Adding geq constraint" << std::endl;
#endif

                precedenceEncoder(_vars[1], _vars[0], 0, _solver, encoding);

            } else {

#ifdef _DEBUGWRAP
                std::cout << "Adding geq constraint" << std::endl;
#endif

                // x0 >= r
                if(_rhs > _vars[0]->getmax()) {
                    std::cerr << "model is inconsistent, -exiting" <<std::endl;
                    exit(1);
                } else if(_rhs > _vars[0]->getmin()) {
                    if(encoding->order){
                        lits.clear();
                        lits.push_back(_vars[0]->greater_than(_rhs-1));
                        _solver->addClause(lits);
                    } else if(encoding->direct){
                        SatWrapper_ConstantInt *RHS = new SatWrapper_ConstantInt(_rhs);
                        RHS->encoding = encoding;
                        precedenceEncoder(RHS, _vars[0], 0, _solver, encoding);
                        delete RHS;
                    } else {
                        std::cerr << "ERROR: SatWrapper_ge not implemented for this encoding, exiting." << std::endl;
                        exit(1);
                    }
                }
            }
        } else {

#ifdef _DEBUGWRAP
            std::cout << "Adding geq predicate" << std::endl;
#endif

            if(_vars[1]) {

                domain->encode(_solver);
                _vars[1] = _vars[1]->add(_solver, false);

                precedenceEncoder(_vars[1], _vars[0], this, 0, _solver, encoding);

            }
        }

        _solver->validate();
    }

    return this;
}



/* Lt object */

SatWrapper_lt::SatWrapper_lt(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_binop(var1,var2) {
    domain = new DomainEncoding(this);
}

SatWrapper_lt::SatWrapper_lt(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_binop(var1,rhs) {
    domain = new LeqDomain(this,var1->domain, rhs-1, 1);
}

SatWrapper_lt::~SatWrapper_lt() {

#ifdef _DEBUGWRAP
    std::cout << "delete lessthan" << std::endl;
#endif

}

SatWrapper_Expression* SatWrapper_lt::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        std::vector<Lit> lits;

        _vars[0] = _vars[0]->add(_solver, false);
        if(top_level) {
            if(_vars[1]) {
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "Adding lt constraint" << std::endl;
#endif

                precedenceEncoder(_vars[0], _vars[1], 1, _solver, encoding);
            } else {
#ifdef _DEBUGWRAP
                std::cout << "Adding lt constraint" << std::endl;
#endif

                // x0 < r
                if(_rhs <= _vars[0]->getmin()) {
                    std::cerr << "model is inconsistent, -exiting" <<std::endl;
                    exit(1);
                } else if(_rhs <= _vars[0]->getmax()) {
                    SatWrapper_ConstantInt *RHS = new SatWrapper_ConstantInt(_rhs);
                    RHS->encoding = encoding;
                    precedenceEncoder(_vars[0], RHS, 1, _solver, encoding);
                    delete RHS;
                }
            }

        } else {
#ifdef _DEBUGWRAP
            std::cout << "Adding lt predicate" << std::endl;
#endif

            if(_vars[1]) {
                domain->encode(_solver);
                _vars[1] = _vars[1]->add(_solver, false);

                precedenceEncoder(_vars[0], _vars[1], this, 1, _solver, encoding);
            }
        }

        _solver->validate();
    }

    return this;
}


/* Gt object */

SatWrapper_gt::SatWrapper_gt(SatWrapper_Expression *var1, SatWrapper_Expression *var2)
    : SatWrapper_binop(var1,var2) {
    domain = new DomainEncoding(this);
}

SatWrapper_gt::SatWrapper_gt(SatWrapper_Expression *var1, int rhs)
    : SatWrapper_binop(var1,rhs) {
    domain = new LeqDomain(this, var1->domain, rhs, 0);
}

SatWrapper_gt::~SatWrapper_gt() {
}

SatWrapper_Expression* SatWrapper_gt::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        std::vector<Lit> lits;

        _vars[0] = _vars[0]->add(_solver, false);
        if(top_level) {

            if(_vars[1]) {
                _vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
                std::cout << "Adding gt constraint" << std::endl;
#endif

                precedenceEncoder(_vars[1], _vars[0], 1, _solver, encoding);

            } else {

#ifdef _DEBUGWRAP
                std::cout << "Adding gt constraint" << std::endl;
#endif

                // x0 > r
                if(_rhs >= _vars[0]->getmax()) {
                    std::cerr << "model is inconsistent, -exiting" <<std::endl;
                    exit(1);
                } else if(_rhs >= _vars[0]->getmin()) {
                    if(encoding->order){
                        lits.clear();
                        lits.push_back(_vars[0]->greater_than(_rhs));
                        _solver->addClause(lits);
                    } else if(encoding->direct){
                        SatWrapper_ConstantInt *RHS = new SatWrapper_ConstantInt(_rhs);
                        RHS->encoding = encoding;
                        precedenceEncoder(RHS, _vars[0], 1, _solver, encoding);
                        delete RHS;
                    }
                }
            }

        } else {

#ifdef _DEBUGWRAP
            std::cout << "Adding gt predicate" << std::endl;
#endif

            if(_vars[1]) {
                domain->encode(_solver);
                _vars[1] = _vars[1]->add(_solver, false);

                precedenceEncoder(_vars[1], _vars[0], this, 1, _solver, encoding);

            }
        }

        _solver->validate();
    }

    return this;
}




/* Minimise object */

SatWrapper_Minimise::SatWrapper_Minimise(SatWrapper_Expression *var)
    : SatWrapper_Expression() {
    _obj = var;
}

SatWrapper_Minimise::~SatWrapper_Minimise() {
}

SatWrapper_Expression* SatWrapper_Minimise::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        _obj = _obj->add(_solver, false);

        if(top_level) {

#ifdef _DEBUGWRAP
            std::cout << "Adding minimise objective" << std::endl;
#endif

            _solver->minimise_obj = _obj;

        }
    }

    return this;
}

/* Maximise object */

SatWrapper_Maximise::SatWrapper_Maximise(SatWrapper_Expression *var)
    : SatWrapper_Expression() {
    _obj = var;
}

SatWrapper_Maximise::~SatWrapper_Maximise() {
}

SatWrapper_Expression* SatWrapper_Maximise::add(SatWrapperSolver *solver, bool top_level) {
    if(!has_been_added()) {
        _solver = solver;
        _ident = _solver->declare(this, false);
        
        // If the encoding hasn't been overwritten for this expression, then we take the default for the solver.
        if(!encoding) encoding = solver->encoding;

        _obj = _obj->add(_solver, false);

        if(top_level) {

#ifdef _DEBUGWRAP
            std::cout << "Adding minimise objective" << std::endl;
#endif

            _solver->maximise_obj = _obj;

        }
    }

    return this;
}




/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

SatWrapperSolver::SatWrapperSolver() {

#ifdef _DEBUGWRAP
    std::cout << "create a SatWrapper solver" << std::endl;
#endif

    minimise_obj = NULL;
    maximise_obj = NULL;
    encoding = NULL;
    cp_model = NULL;

    Lit dummy(create_atom(NULL,0),true);
    Lit_True = Lit(dummy);
    Lit_False = ~Lit(dummy);

    current = 0;
    clause_limit = -1;
}

SatWrapperSolver::~SatWrapperSolver() {

#ifdef _DEBUGWRAP
    std::cout << "delete wrapped solver" << std::endl;
#endif

}

int SatWrapperSolver::declare(SatWrapper_Expression *exp, bool type) {
    int id = _expressions.size();
    _expressions.push_back(exp);
    return id;
}

int SatWrapperSolver::get_cb_size() {
    return clause_base.size();
}

int SatWrapperSolver::create_atom(DomainEncoding* dom, const int type) {
    unsigned int id = _atom_to_domain.size();
    _atom_to_domain.push_back(dom);
    _atom_to_type.push_back(type);
    return id;
}

void SatWrapperSolver::addClause(std::vector<Lit>& cl) {

    if(clause_limit != -1) {
        // We are limiting clauses
        if((unsigned int)clause_limit < clause_base.size()) {
            //std::cerr << "Warning: Clause limit reached" << std::endl;
            cl.clear();
            return;
        }
    }

    std::vector<Lit> clause;
    if(!processClause(cl,clause)) {
        clause_base.push_back(clause);
#ifdef _DEBUGWRAP
        displayClause(clause);
#endif
    }
}

void SatWrapperSolver::validate() {
    while(current < clause_base.size()) {

#ifdef _DEBUGWRAP
        std::cout << "add ";
        displayClause(clause_base[current]);
#endif

        ++current;
    }
}

void SatWrapperSolver::store_solution() {

#ifdef _DEBUGWRAP
    std::cout << "store a new solution" << std::endl;
#endif

}

void SatWrapperSolver::add(SatWrapper_Expression* arg) {

#ifdef _DEBUGWRAP
    std::cout << "add an expression to the solver" << std::endl;
    std::cout << "Solver is using encoding: " << *encoding << std::endl;
#endif

    if(arg != NULL) arg->add(this, true);
}

lbool SatWrapperSolver::truth_value(Lit x) {
    return l_Undef;
}

void SatWrapperSolver::initialise(SatWrapperExpArray& arg) {

#ifdef _DEBUGWRAP
    std::cout << "initialise the solver" << std::endl;
#endif

}

void SatWrapperSolver::initialise() {

#ifdef _DEBUGWRAP
    std::cout << "initialise the solver" << std::endl;
#endif

}

int SatWrapperSolver::solveAndRestart(const int policy,
                                      const unsigned int base,
                                      const double factor,
                                      const double decay) {

#ifdef _DEBUGWRAP
    std::cout << "call solve & restart" << std::endl;
#endif

    return is_sat();
}

int SatWrapperSolver::solve() {

#ifdef _DEBUGWRAP
    std::cout << "call solve" << std::endl;
#endif

    return is_sat();
}

bool SatWrapperSolver::propagate() {

#ifdef _DEBUGWRAP
    std::cout << "call propagate" << std::endl;
#endif

    return true;
}

void SatWrapperSolver::reset(bool full) {

#ifdef _DEBUGWRAP
    std::cout << "call reset" << std::endl;
#endif

}

bool SatWrapperSolver::undo(const int nlevel) {

#ifdef _DEBUGWRAP
    std::cout << "call undo" << std::endl;
#endif

    return true;
}

bool SatWrapperSolver::branch_right() {

#ifdef _DEBUGWRAP
    std::cout << "call branch right" << std::endl;
#endif

    return true;
}

void SatWrapperSolver::deduce() {

#ifdef _DEBUGWRAP
    std::cout << "call branch right" << std::endl;
#endif

}

void SatWrapperSolver::save() {

#ifdef _DEBUGWRAP
    std::cout << "call save" << std::endl;
#endif

}

void SatWrapperSolver::post(const char* op, SatWrapper_Expression* x, int v) {

#ifdef _DEBUGWRAP
    std::cout << "call save" << std::endl;
#endif

}

int SatWrapperSolver::startNewSearch() {

#ifdef _DEBUGWRAP
    std::cout << "start a new interuptable search" << std::endl;
#endif

    return 1;
}

int SatWrapperSolver::getNextSolution() {

#ifdef _DEBUGWRAP
    std::cout << "seek next solution" << std::endl;
#endif

    return 1;
}

int SatWrapperSolver::sacPreprocess(const int type) {

#ifdef _DEBUGWRAP
    std::cout << "enforces singleton arc consistency" << std::endl;
#endif

    return 1;
}

void SatWrapperSolver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand) {

#ifdef _DEBUGWRAP
    std::cout << "set the variable/value ordering" << std::endl;
#endif

}

void SatWrapperSolver::setFailureLimit(const int cutoff) {

#ifdef _DEBUGWRAP
    std::cout << "set a cutoff on the number of fails" << std::endl;
#endif

}

void SatWrapperSolver::setNodeLimit(const int cutoff) {

#ifdef _DEBUGWRAP
    std::cout << "set a cutoff on the number of nodes" << std::endl;
#endif

}

void SatWrapperSolver::setTimeLimit(const double cutoff) {

#ifdef _DEBUGWRAP
    std::cout << "set a time cutoff" << std::endl;
#endif

}

void SatWrapperSolver::setVerbosity(const int degree) {

#ifdef _DEBUGWRAP
    std::cout << "set the verbosity level" << std::endl;
#endif

}

void SatWrapperSolver::setRandomized(const int degree) {

#ifdef _DEBUGWRAP
    std::cout << "set the type of randomization" << std::endl;
#endif

}

void SatWrapperSolver::setRandomSeed(const int seed) {

#ifdef _DEBUGWRAP
    std::cout << "set the random seed" << std::endl;
#endif
}

bool SatWrapperSolver::is_sat() {

#ifdef _DEBUGWRAP
    std::cout << "was a solution found?" << std::endl;
#endif

    return false;
}

bool SatWrapperSolver::is_unsat() {

#ifdef _DEBUGWRAP
    std::cout << "was unfeasibility proven?" << std::endl;
#endif

    return false;
}

bool SatWrapperSolver::is_opt() {

#ifdef _DEBUGWRAP
    std::cout << "was optimality proven?" << std::endl;
#endif

    return false;
}

void SatWrapperSolver::printStatistics() {

#ifdef _DEBUGWRAP
    std::cout << "print a bunch of statistics" << std::endl;
#endif

}

int SatWrapperSolver::getBacktracks() {

#ifdef _DEBUGWRAP
    std::cout << "return the number of backtracks" << std::endl;
#endif

    return 0;
}

int SatWrapperSolver::getNodes() {

#ifdef _DEBUGWRAP
    std::cout << "return the number of nodes" << std::endl;
#endif

    return 0;
}

int SatWrapperSolver::getFailures() {

#ifdef _DEBUGWRAP
    std::cout << "return the number of failures" << std::endl;
#endif

    return 0;
}

int SatWrapperSolver::getChecks() {

#ifdef _DEBUGWRAP
    std::cout << "return the number of checks" << std::endl;
#endif

    return 0;
}

int SatWrapperSolver::getPropags() {

#ifdef _DEBUGWRAP
    std::cout << "return the number of propags" << std::endl;
#endif

    return 0;
}

double SatWrapperSolver::getTime() {

#ifdef _DEBUGWRAP
    std::cout << "return the cpu time" << std::endl;
#endif

    return 0.0;
}

void SatWrapperSolver::displayClause(std::vector<Lit>& cl) {
    if(cl.size()) {
        for(unsigned int i=0; i<cl.size(); ++i) {
            std::cout << " " ;
            if(sign(cl[i])) std::cout << "~";
            std::cout << var(cl[i]);
        }
        std::cout << " / ";
        std::cout.flush();
        Lit p = cl[0];
        displayLiteral(p);
        for(unsigned int i=1; i<cl.size(); ++i) {
            std::cout << " or ";
            displayLiteral(cl[i]);
        }
        std::cout << std::endl;
    } else std::cout << "{}" ;
}

void SatWrapperSolver::displayLiteral(Lit p) {
    int x = var(p);
    if(x>=0) {
        if(_atom_to_domain[x] != NULL)
            _atom_to_domain[x]->print_lit(p,_atom_to_type[x]);
        else std::cout << "(aux)";
    } else std::cout << "false" ;
}

void SatWrapperSolver::output_cnf(const char *filename){
    std::ofstream f(filename);

    // CNF header line. Number of atoms - 1 because of the dummy literal.
    f << "p cnf " << (_atom_to_domain.size() - 1) << " " << clause_base.size() << std::endl;
    for(unsigned int i=0;i<clause_base.size();i++){
        std::vector<Lit> clause = clause_base[i];
        if(clause.size() > 0){
            for(unsigned int j=0;j<clause.size();j++){
                if(sign(clause[j])) f << "-";
                f << var(clause[j]) << " ";
            }
            f << "0" << std::endl;
        }
    }
    f.close();
}

void SatWrapperSolver::setClauseLimit(int limit) {
    //std::cout << "Clause limit set to " << limit << std::endl;
    this->clause_limit = limit;
}
