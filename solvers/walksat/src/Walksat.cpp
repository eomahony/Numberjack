
#include "Walksat.hpp"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

Lit Lit_False;
Lit Lit_True;


static inline double cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }



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

DomainEncoding::DomainEncoding(Walksat_Expression *o) : AbstractDomain(o) {
  _lower = 0;
  _upper = 1;
  _size = 2;
  _values = NULL;

  // not initialised
  _direct_encoding = -1;
  _order_encoding = -1;
}

DomainEncoding::DomainEncoding(Walksat_Expression *o, const int nval) : AbstractDomain(o) {
  _lower = 0;
  _upper = nval-1;
  _size = nval;
  _values = NULL;

  // not initialised
  _direct_encoding = -1;
  _order_encoding = -1;
}

DomainEncoding::DomainEncoding(Walksat_Expression *o, const int lb, const int ub) : AbstractDomain(o) {
  _lower = lb;
  _upper = ub;
  _size = ub-lb+1;
  _values = NULL;

  // not initialised
  _direct_encoding = -1;
  _order_encoding = -1;
}

DomainEncoding::DomainEncoding(Walksat_Expression *o, WalksatIntArray& vals) : AbstractDomain(o) {
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

Lit DomainEncoding::less_or_equal(const int value, const int index) const {
  if(_lower > value) return Lit_False;
  else if(_upper <= value) return Lit_True;
  else if(_size == 2) return ~Lit(_direct_encoding);
  else if(index >= 0) return Lit(_order_encoding+index);
  else if(!_values) return Lit(_order_encoding+value-_lower);

  // We need to search for the right variable, the index of 'value' 
  // if it is in the domain, or the index of the largest element
  // smaller than 'value' otherwise
  return Lit(_order_encoding+get_index_p(value));
} 

Lit DomainEncoding::equal(const int value, const int index) const {
  if(_lower > value || _upper < value) return Lit_False;
  else if(_size == 2) {
    if(_lower == value) return ~Lit(_direct_encoding);
    if(_upper == value) return Lit(_direct_encoding);
    return Lit_False;
  } else {
    if(index >= 0) {return Lit(_direct_encoding+index);}
    if(!_values) return Lit(_direct_encoding+value-_lower);

    // We need to search for the right variable,: the index of 'value' 
    int x = get_index_p(value);
    if(_values[x] == value) return Lit(_direct_encoding+x);
  }
  return Lit_False;
} 

void DomainEncoding::encode(WalksatSolver *solver)
{

#ifdef _DEBUGWRAP
  std::cout << "encode x" << owner->_ident << " in ";
  display(std::cout);
  std::cout << std::endl;
#endif

  if(_size == 2) _direct_encoding = solver->create_atom(this,SELF);
  else {
    // create one atom for each value in the domain (direct encoding)
    _direct_encoding = solver->create_atom(this,DIRECT);
    for(int i=1; i<_size; ++i) solver->create_atom(this,DIRECT);
    
    // create one atom for each value in the domain but the last (order encoding)
    _order_encoding = solver->create_atom(this,ORDER);
    for(int i=1; i<_size-1; ++i) solver->create_atom(this,ORDER);
  
    std::vector<Lit> lits;
     
    // at least one value x=1 or x=2 or ...
    for(int i=0; i<_size; ++i) lits.push_back(equal(getval(i),i));
    solver->addClause(lits);

    for(int i=0; i<_size; ++i) {
      // chain the bounds x<=i -> x<=i+1
      if(i) {
	lits.clear();
	lits.push_back(~(less_or_equal(getval(i-1),i-1)));
	lits.push_back(less_or_equal(getval(i),i));
	solver->addClause(lits);
      }

      // channel x=i -> x>i-1
      if(i) {
	lits.clear();
	lits.push_back(~(equal(getval(i),i)));
	lits.push_back(~(less_or_equal(getval(i-1),i-1)));
	solver->addClause(lits);
      }

      // channel x=i -> x<=i
      if(i<_size-1) {
	lits.clear();
	lits.push_back(~equal(getval(i),i));
	lits.push_back(less_or_equal(getval(i),i));
	solver->addClause(lits);
      }
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

OffsetDomain::OffsetDomain(Walksat_Expression *o, AbstractDomain *d, const int os) : AbstractDomain(o) {
  _dom_ptr = d;
  offset = os;
}

int OffsetDomain::getval(int idx) const { return (_dom_ptr->getval(idx)+offset); }
int OffsetDomain::getmin() const { return (_dom_ptr->getmin()+offset); }
int OffsetDomain::getmax() const { return (_dom_ptr->getmax()+offset); }

Lit OffsetDomain::less_or_equal(const int value, const int index) const {
  return _dom_ptr->less_or_equal(value-offset,index);
}
Lit OffsetDomain::equal(const int value, const int index) const {
  return _dom_ptr->equal(value-offset,index);
}

FactorDomain::FactorDomain(Walksat_Expression *o, AbstractDomain *d, const int f) : AbstractDomain(o) {
  _dom_ptr = d;
  factor = f;
}

int FactorDomain::getval(int idx) const { return ((factor < 0 ? _dom_ptr->getval(_dom_ptr->getsize()-idx-1) : _dom_ptr->getval(idx))*factor); }
int FactorDomain::getmin() const { return (factor < 0 ? factor*(_dom_ptr->getmax()) : factor*(_dom_ptr->getmin())); }
int FactorDomain::getmax() const { return (factor < 0 ? factor*(_dom_ptr->getmin()) : factor*(_dom_ptr->getmax())); }

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

EqDomain::EqDomain(Walksat_Expression *o, AbstractDomain *d, const int v, const int s) : AbstractDomain(o) {
  _dom_ptr = d;
  value = v;
  spin = s;
}

Lit EqDomain::less_or_equal(const int v, const int index) const {
  // make sense iff v == 0
  if(v < 0) return Lit_False;
  if(v > 0) return Lit_True;
  // <= 0, i.e., false iff the pointed domain is not equal to 'value' (reverse if neg spin)
  return (spin ? ~(_dom_ptr->equal(value,-1)) : _dom_ptr->equal(value,-1));
}
Lit EqDomain::equal(const int v, const int index) const {
  // make sense iff v == 0 or v == 1
  if(v < 0 || v > 1) return Lit_False;
  // find out if it is in fact an equality or an inequality
  return ((v == spin) ? _dom_ptr->equal(value,-1) : ~(_dom_ptr->equal(value,-1)));
}

LeqDomain::LeqDomain(Walksat_Expression *o, AbstractDomain *d, const int b, const int s) : AbstractDomain(o) {
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

Lit ConstantDomain::less_or_equal(const int v, const int index) const {if(value <= v) return Lit_True; else return Lit_False;}
Lit ConstantDomain::equal(const int v, const int index) const {if(value == v) return Lit_True; else return Lit_False;}

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

// (X + Y) == Z
void additionEncoder(Walksat_Expression *X,
		     Walksat_Expression *Y,
		     Walksat_Expression *Z,
		     WalksatSolver *solver) {
  std::vector<Lit> lits;
  int i, j, x, y;

  for(i=0; i<X->getsize(); ++i) {
    x = X->getval(i);
    for(j=0; j<Y->getsize(); ++j) {
      y = Y->getval(j);
      
      lits.clear();
      if(i) lits.push_back(X->less_or_equal(x-1));
      if(j) lits.push_back(Y->less_or_equal(y-1));
      if(i || j) {
	lits.push_back(~(Z->less_or_equal(x+y-1)));
	solver->addClause(lits);
      }

      lits.clear();
      if(i<X->getsize()) lits.push_back(~(X->less_or_equal(x)));
      if(j<Y->getsize())lits.push_back(~(Y->less_or_equal(y)));
      if(i<X->getsize() || j<Y->getsize()) {
	lits.push_back(Z->less_or_equal(x+y));
	solver->addClause(lits);
      }
    }
  }
}

// (X != Y)
void disequalityEncoder(Walksat_Expression *X,
			Walksat_Expression *Y,
			WalksatSolver *solver) {
  std::vector<Lit> lits;
  int i=0, j=0, x, y;
  while( i<X->getsize() && j<Y->getsize() )
    {
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
}

// (X == Y)
void equalityEncoder(Walksat_Expression *X,
		     Walksat_Expression *Y,
		     WalksatSolver *solver) {
  std::vector<Lit> lits;
  int i=0, j=0, x, y;
  while( i<X->getsize() && j<Y->getsize() )
    {
      x = X->getval(i);
      y = Y->getval(j);

      if(x == y) {
	lits.clear();
	lits.push_back(~(X->equal(x,i)));
	lits.push_back(Y->equal(y,j));
	solver->addClause(lits);

	lits.clear();
	lits.push_back(X->equal(x,i));
	lits.push_back(~(Y->equal(y,j)));
	solver->addClause(lits);
	++i;
	++j;
      }
      if( x < y ) {
	lits.clear();
	lits.push_back(~(X->equal(x,i)));
	solver->addClause(lits);
	++i;
      }
      if( y < x ) {
	lits.clear();
	lits.push_back(~(Y->equal(y,j)));
	solver->addClause(lits);
	++j;
      }
    }
}

// (X == Y) <-> Z
void equalityEncoder(Walksat_Expression *X,
		     Walksat_Expression *Y,
		     Walksat_Expression *Z,
		     WalksatSolver *solver,
		     const bool spin) {
  unsigned int num_clauses = solver->clause_base.size();
  equalityEncoder(X,Y,solver);
  while(num_clauses < solver->clause_base.size()) 
    solver->clause_base[num_clauses++].push_back((spin ? Z->equal(0) : ~(Z->equal(0))));
  disequalityEncoder(X,Y,solver);
  while(num_clauses < solver->clause_base.size())
    solver->clause_base[num_clauses++].push_back((spin ? ~(Z->equal(0)) : Z->equal(0)));
}


// X+K <= Y
void precedenceEncoder(Walksat_Expression *X,
		       Walksat_Expression *Y,
		       const int K,
		       WalksatSolver *solver) {
  std::vector<Lit> lits;
  int i=0, j=0, x=0, y=0;
  
  while( i<X->getsize() )
    {
      if(i<X->getsize()) x = X->getval(i);
      if(j<Y->getsize()) y = Y->getval(j);

      std::cout << i << " " << j << " | " << x << " " << y << std::endl;


      if(x+K == y) {
	lits.clear();
	lits.push_back(X->less_or_equal(x,i));
	lits.push_back(~(Y->less_or_equal(y,j)));
	solver->addClause(lits);
	++i;
	++j;
      } else if(x+K < y) {
	++i;
      } else if(x+K > y) {
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
}

// (X+K <= Y) <-> Z
void precedenceEncoder(Walksat_Expression *X,
		       Walksat_Expression *Y,
		       Walksat_Expression *Z,
		       const int K,
		       WalksatSolver *solver) {
  unsigned int num_clauses = solver->clause_base.size();
  precedenceEncoder(X,Y,K,solver);
  while(num_clauses < solver->clause_base.size()) 
    solver->clause_base[num_clauses++].push_back(Z->equal(0));
  precedenceEncoder(Y,X,1-K,solver);
  while(num_clauses < solver->clause_base.size())
    solver->clause_base[num_clauses++].push_back(~(Z->equal(0)));
}


/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

void Walksat_Expression::initialise() { 
  _ident = -1;
  _solver = NULL;
  domain = NULL;
}

Walksat_Expression::Walksat_Expression() {
  initialise();
  domain = new DomainEncoding(this);
}

Walksat_Expression::Walksat_Expression(const int nval) {
  initialise();
  domain = new DomainEncoding(this,nval);

#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif

}

Walksat_Expression::Walksat_Expression(const int lb, const int ub) {
  initialise();
  domain = new DomainEncoding(this,lb,ub);

#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << lb << ".." << ub << "]" << std::endl;
#endif

}

Walksat_Expression::Walksat_Expression(WalksatIntArray& vals) {
  initialise();
  domain = new DomainEncoding(this,vals);

#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

}

Walksat_Expression::~Walksat_Expression() {
  delete domain;
}

int Walksat_Expression::get_size() const
{
  int i=0, domsize=0;
  for(i=0; i<getsize(); ++i)
    domsize += (_solver->truth_value(equal(getval(i),i)) != l_False);
  return domsize;
}

int Walksat_Expression::next(const int v) const
{
  int nxt = v;
  while( ++nxt <= getmax() ) 
    if(_solver->truth_value(equal(nxt)) != l_False) break;
  if(nxt > getmax()) nxt = v;
  return nxt;
}

int Walksat_Expression::get_min() const
{
  for(int i=0; i<getsize(); ++i) {
    if(_solver->truth_value(equal(getval(i),i)) != l_False) return getval(i);
  }
  return getmin();
}

int Walksat_Expression::get_max() const
{
  for(int i=getsize()-1; i>=0; --i)
    if(_solver->truth_value(equal(getval(i),i)) != l_False) return getval(i);
  return getmax();
}

bool Walksat_Expression::contain(const int v) const {
  return (_solver->truth_value(equal(v)) != l_False);
}

int Walksat_Expression::get_value() const
{
  if(_solver->cp_model) {
    return _solver->cp_model[_ident];
  } else 
    return get_min();
}

bool Walksat_Expression::has_been_added() const {
  return (_solver != NULL);
}

Walksat_Expression* Walksat_Expression::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {

    _solver = solver;
    _ident = _solver->declare(this, true);
    solver->_variables.push_back(this);

#ifdef _DEBUGWRAP
    std::cout << "+ add x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif
    
    domain->encode(_solver);
  }

  return this;
}


Walksat_add::Walksat_add(Walksat_Expression *arg1, Walksat_Expression *arg2)
  : Walksat_binop(arg1, arg2) {
  int _lower = arg1->getmin()+arg2->getmin();
  int _upper = arg1->getmax()+arg2->getmax();
  domain = new DomainEncoding(this, _lower, _upper);
}

Walksat_add::Walksat_add(Walksat_Expression *arg1, const int arg2)
  : Walksat_binop(arg1, arg2) {
  domain = new OffsetDomain(this,arg1->domain, arg2);

#ifdef _DEBUGWRAP
  std::cout << "creating offset expression [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

}

Walksat_Expression* Walksat_add::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this, false);

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

	additionEncoder(_vars[0], _vars[1], this, _solver);

      }
    }

    _solver->validate();
  } 

  return this;
}

Walksat_add::~Walksat_add() {}

Walksat_mul::Walksat_mul(Walksat_Expression *arg1, const int arg2)
  : Walksat_binop(arg1, arg2) {

  domain = new FactorDomain(this,arg1->domain, arg2);

}

Walksat_mul::Walksat_mul(Walksat_Expression *arg1, Walksat_Expression *arg2)
  : Walksat_binop(arg1, arg2) {

  std::cerr << "c NOT SUPPORTED (multiplication) - exiting" << std::endl;
  exit(1);

}

Walksat_Expression* Walksat_mul::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

#ifdef _DEBUGWRAP
    std::cout << "creating mul expression x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

    if(top_level) {
      
      std::cerr << "mul predicate at the top level not supported" << std::endl;
      exit(1);
      
    } else {
      
      _vars[0] = _vars[0]->add(_solver, false);

      if(_vars[1]) {

	std::cerr << "c NOT SUPPORTED (multiplication) - exiting" << std::endl;
	exit(1);
	
      } 
    }

    _solver->validate();
  } 

  return this;
}

Walksat_mul::~Walksat_mul() {}

void Walksat_AllDiff::addVar(Walksat_Expression* v) {
  _vars.add(v);
}

Walksat_AllDiff::Walksat_AllDiff( Walksat_Expression* arg1, Walksat_Expression* arg2 ) 
  : Walksat_Expression() {
  addVar(arg1);
  addVar(arg2);
}

Walksat_AllDiff::Walksat_AllDiff( WalksatExpArray& vars ) 
  : Walksat_Expression() {
  _vars = vars;
}

Walksat_AllDiff::~Walksat_AllDiff() {
  for(unsigned int i=0; i<_clique.size(); ++i)
    delete _clique[i];
}

Walksat_Expression* Walksat_AllDiff::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

#ifdef _DEBUGWRAP
    std::cout << "creating alldiff expression" << std::endl;
#endif

    if(top_level){
      
      int i, j, n=_vars.size();  
      for(i=0; i<n; ++i) 
	_vars.set_item(i, (_vars.get_item(i))->add(_solver,false));

      Walksat_Expression *exp;
      for(i=1; i<n; ++i)
	for(j=0; j<i; ++j) {
	  exp = new Walksat_ne(_vars.get_item(j), _vars.get_item(i));
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


Walksat_Sum::Walksat_Sum(WalksatExpArray& vars, 
				 WalksatIntArray& weights, 
				 const int offset)
  : Walksat_Expression() {
  _offset = offset;
  _vars = vars;
  _weights = weights;
  initialise();
}

Walksat_Sum::Walksat_Sum(Walksat_Expression *arg1, 
				 Walksat_Expression *arg2, 
				 WalksatIntArray& w, 
				 const int offset)
  : Walksat_Expression() {
  _offset = offset;
  _vars.add(arg1);
  _vars.add(arg2);
  _weights = w;
  initialise();
}

Walksat_Sum::Walksat_Sum(Walksat_Expression *arg, 
				 WalksatIntArray& w, 
				 const int offset)
  : Walksat_Expression() {
  _self = NULL;
  _offset = offset;
  _vars.add(arg);
  _weights = w;
  initialise();
}

Walksat_Sum::Walksat_Sum()
  : Walksat_Expression() {
  _offset = 0;
}

void Walksat_Sum::initialise() {
  if(_vars.size() == 2) {
    int _lower = 0;
    int _upper = 0;

    int weight;
    for(int i = 0; i < _vars.size(); ++i){
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

Walksat_Sum::~Walksat_Sum() {

#ifdef _DEBUGWRAP
  std::cout << "delete sum" << std::endl;
#endif

  for(unsigned int i=0; i<_subsum.size(); ++i) {
    delete _subsum[i];
  }
}

void Walksat_Sum::addVar(Walksat_Expression* v) {
  _vars.add(v);
}

void Walksat_Sum::addWeight(const int w) {
  _weights.add(w);
}

void Walksat_Sum::set_rhs(const int k) {
  _offset = k;
}

Walksat_Expression* Walksat_Sum::add(WalksatSolver *solver, bool top_level){
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this, false);
    
#ifdef _DEBUGWRAP    
      std::cout << "add sum expression x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

    if(!top_level) {
      
      if(_vars.size() == 1) {

	_vars.set_item(0, _vars.get_item(0)->add(_solver, false));

      } else {

	int w1, w2;
	Walksat_Expression *exp, *exp1, *exp2;
	for(int i=0; i+1<_vars.size(); i+=2) {
	  w1 = _weights.get_item(i);
	  w2 = _weights.get_item(i+1);

	  if(w1 != 1) exp1 = new Walksat_mul(_vars.get_item(i), w1);
	  else exp1 = _vars.get_item(i);

	  if(w2 != 1) exp2 = new Walksat_mul(_vars.get_item(i+1), w2);
	  else exp2 = _vars.get_item(i+1);

	  exp = new Walksat_add(exp1, exp2);
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

Walksat_binop::Walksat_binop(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_Expression() {
  _vars[0] = var1;
  _vars[1] = var2;

#ifdef _DEBUGWRAP
  std::cout << "creating binary operator" << std::endl;
#endif

}

Walksat_binop::Walksat_binop(Walksat_Expression *var1, int rhs)
  : Walksat_Expression() {
  _vars[0] = var1;
  _vars[1] = NULL;
  _rhs = rhs;

#ifdef _DEBUGWRAP
  std::cout << "creating binary operator" << std::endl;
#endif

}


Walksat_binop::~Walksat_binop() {

#ifdef _DEBUGWRAP
  std::cout << "delete binary operator" << std::endl;
#endif

}

Walksat_or::Walksat_or(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_binop(var1,var2) {

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  domain = new DomainEncoding(this);

}

Walksat_or::Walksat_or(Walksat_Expression *var1, int rhs)
  : Walksat_binop(var1,rhs) {

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  if(rhs) domain = new ConstantDomain(this,1);
  else domain = new OffsetDomain(this,var1->domain,0);
}


Walksat_or::~Walksat_or() {

#ifdef _DEBUGWRAP
  std::cout << "delete or" << std::endl;
#endif

}

Walksat_Expression* Walksat_or::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

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



Walksat_and::Walksat_and(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_binop(var1,var2) {

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  domain = new DomainEncoding(this);
}

Walksat_and::Walksat_and(Walksat_Expression *var1, int rhs)
  : Walksat_binop(var1,rhs) {

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  if(rhs) domain = new OffsetDomain(this,var1->domain,0);
  else domain = new ConstantDomain(this,0);
}


Walksat_and::~Walksat_and() {

#ifdef _DEBUGWRAP
  std::cout << "delete or" << std::endl;
#endif

}

Walksat_Expression* Walksat_and::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

    std::vector<Lit> lits;
    if(top_level) {
      
      if(_vars[1]) {
	
	_vars[0] = _vars[0]->add(_solver, false);
	_vars[1] = _vars[1]->add(_solver, false);
	
#ifdef _DEBUGWRAP
	std::cout << "add or constraint" << std::endl;
#endif
	
	for(int i=0; i<2; ++i) {
	  lits.clear();
	  lits.push_back(~(_vars[i]->equal(0)));
	  _solver->addClause(lits);
	}
	
      } else if(_rhs != 0) {
	
	_vars[0] = _vars[0]->add(_solver, false);
	
#ifdef _DEBUGWRAP
	std::cout << "add or constraint" << std::endl;
#endif
	
	lits.clear();
	lits.push_back(~(_vars[0]->equal(0)));
	_solver->addClause(lits);
      }
      
    } else {
      
      if(_vars[1]) {
	
	domain->encode(_solver);
	
	_vars[0] = _vars[0]->add(_solver, false);
	_vars[1] = _vars[1]->add(_solver, false);
	
#ifdef _DEBUGWRAP
	std::cout << "add or constraint" << std::endl;
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


Walksat_eq::Walksat_eq(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_binop(var1,var2) {
  domain = new DomainEncoding(this);
}

Walksat_eq::Walksat_eq(Walksat_Expression *var1, int rhs)
  : Walksat_binop(var1,rhs) {
  domain = new EqDomain(this,var1->domain, rhs, 1);
}


Walksat_eq::~Walksat_eq() {

#ifdef _DEBUGWRAP
  std::cout << "delete eq" << std::endl;
#endif

}

Walksat_Expression* Walksat_eq::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

#ifdef _DEBUGWRAP
    std::cout << "creating eq expression x" << _ident << " [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

    std::vector<Lit> lits;

    _vars[0] = _vars[0]->add(_solver, false);
    if(top_level) {
            
      if(_vars[1]) {
  
	_vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
	std::cout << "encode equality constraint x" 
		  << _vars[0]->_ident << " == x" << _vars[1]->_ident << std::endl;
#endif

	equalityEncoder(_vars[0], _vars[1], _solver);

      } else {

#ifdef _DEBUGWRAP
	std::cout << "encode unary equality constraint x" 
		  << _vars[0]->_ident << " == " << _rhs << std::endl;
#endif

	lits.clear();
	lits.push_back(_vars[0]->equal(_rhs));
	_solver->addClause(lits);

      }
    } else {
      
      if(_vars[1]) {
  
	domain->encode(_solver);
	_vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
	std::cout << "encode equality predicate x" << this->_ident << " == (x"
		  << _vars[0]->_ident << " == x" << _vars[1]->_ident << ")" << std::endl;
#endif	

	equalityEncoder(_vars[0], _vars[1], this, _solver, true);

      }
    }

    _solver->validate();
  } 

  return this;
}


/* Disequality operator */

Walksat_ne::Walksat_ne(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_binop(var1,var2) {
  domain = new DomainEncoding(this);
}

Walksat_ne::Walksat_ne(Walksat_Expression *var1, int rhs)
  : Walksat_binop(var1,rhs) {
  domain = new EqDomain(this,var1->domain, rhs, 0);
}

Walksat_ne::~Walksat_ne() {

#ifdef _DEBUGWRAP
  std::cout << "delete notequal" << std::endl;
#endif

}

Walksat_Expression* Walksat_ne::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

    std::vector<Lit> lits;

    _vars[0] = _vars[0]->add(_solver, false);
    if(top_level){

      if(_vars[1]) {
	_vars[1] = _vars[1]->add(_solver, false);
	
#ifdef _DEBUGWRAP
	std::cout << "add notequal constraint" << std::endl;
#endif
	
	disequalityEncoder(_vars[0], _vars[1], _solver);

      } else  {

#ifdef _DEBUGWRAP
	std::cout << "add notequal" << std::endl;
#endif

	lits.clear();
	lits.push_back(~(_vars[0]->equal(_rhs)));
	_solver->addClause(lits);

      }
      
    } else {

      if(_vars[1]) {
	
	domain->encode(_solver);
	_vars[1] = _vars[1]->add(_solver, false);
	
#ifdef _DEBUGWRAP
	std::cout << "add notequal predicate" << std::endl;
#endif
	
	equalityEncoder(_vars[0], _vars[1], this, _solver, false);

      } 
    }

    _solver->validate();
  }

  return this;
}


/* Leq operator */

Walksat_le::Walksat_le(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_binop(var1,var2) { 
#ifdef _DEBUGWRAP
  std::cout << "Creating Le constraint" << std::endl;
#endif

  domain = new DomainEncoding(this);

}

Walksat_le::Walksat_le(Walksat_Expression *var1, int rhs)
  : Walksat_binop(var1,rhs) {

#ifdef _DEBUGWRAP
  std::cout << "Leq on constant constructor" << std::endl;
#endif

  domain = new LeqDomain(this,var1->domain, rhs, 1);
}


Walksat_le::~Walksat_le() {

#ifdef _DEBUGWRAP
  std::cout << "delete lessequal" << std::endl;
#endif

}

Walksat_Expression* Walksat_le::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

    std::vector<Lit> lits;

    _vars[0] = _vars[0]->add(_solver, false);
    if(top_level) {

      if(_vars[1]) {

	_vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding leq constraint" << std::endl;
#endif

      precedenceEncoder(_vars[0], _vars[1], 0, _solver);
      } else {
	_vars[0] = _vars[0]->add(_solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding leq constraint" << std::endl;
#endif
	
	// x0 <= r
	if(_rhs < _vars[0]->getmin()) {
	  std::cerr << "model is inconsistent, -exiting" <<std::endl;
	  exit(1);
	} else if(_rhs < _vars[0]->getmax()) {
	  lits.clear();
	  lits.push_back(_vars[0]->less_or_equal(_rhs));
	  _solver->addClause(lits);
	}
      }
    } else {
  
#ifdef _DEBUGWRAP    
      std::cout << "Adding leq predicate" << std::endl;
#endif

      if(_vars[1]) {

	domain->encode(_solver);
	_vars[1] = _vars[1]->add(_solver, false);

	precedenceEncoder(_vars[0], _vars[1], this, 0, _solver);
	
      } 
    } 

    _solver->validate();
  }

  return this;
}


/* Geq operator */

Walksat_ge::Walksat_ge(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_binop(var1,var2) {
#ifdef _DEBUGWRAP
  std::cout << "Creating Ge constraint" << std::endl;
#endif

  domain = new DomainEncoding(this);
}

Walksat_ge::Walksat_ge(Walksat_Expression *var1, int rhs)
  : Walksat_binop(var1,rhs) {

#ifdef _DEBUGWRAP
  std::cout << "Geq on constant constructor" << std::endl;
#endif

  domain = new LeqDomain(this,var1->domain, rhs-1, 0);
}

Walksat_ge::~Walksat_ge() {

#ifdef _DEBUGWRAP
  std::cout << "delete greaterequal" << std::endl;
#endif

}

Walksat_Expression* Walksat_ge::add(WalksatSolver *solver, bool top_level) { 
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

    std::vector<Lit> lits;

    _vars[0] = _vars[0]->add(_solver, false);
    if(top_level) {
      
      if(_vars[1]) {
	
	_vars[1] = _vars[1]->add(_solver, false);
	
#ifdef _DEBUGWRAP
	std::cout << "Adding geq constraint" << std::endl;
#endif
	
	precedenceEncoder(_vars[1], _vars[0], 0, _solver);

      } else {

#ifdef _DEBUGWRAP
      std::cout << "Adding geq constraint" << std::endl;
#endif
	
	// x0 >= r
	if(_rhs > _vars[0]->getmax()) {
	  std::cerr << "model is inconsistent, -exiting" <<std::endl;
	  exit(1);
	} else if(_rhs > _vars[0]->getmin()) {
	  lits.clear();
	  lits.push_back(_vars[0]->greater_than(_rhs-1));
	  _solver->addClause(lits);
	}
      }
    } else {

#ifdef _DEBUGWRAP    
      std::cout << "Adding geq predicate" << std::endl;
#endif

      if(_vars[1]) {

	domain->encode(_solver);
	_vars[1] = _vars[1]->add(_solver, false);

	precedenceEncoder(_vars[1], _vars[0], this, 0, _solver);

      } 
    } 
    
    _solver->validate();
  }

  return this;
}



/* Lt object */

Walksat_lt::Walksat_lt(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_binop(var1,var2) {
  domain = new DomainEncoding(this);
}

Walksat_lt::Walksat_lt(Walksat_Expression *var1, int rhs)
  : Walksat_binop(var1,rhs) {
  domain = new LeqDomain(this,var1->domain, rhs-1, 1);
}

Walksat_lt::~Walksat_lt() {

#ifdef _DEBUGWRAP
  std::cout << "delete lessthan" << std::endl;
#endif

}

Walksat_Expression* Walksat_lt::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

    std::vector<Lit> lits;

    _vars[0] = _vars[0]->add(_solver, false);
    if(top_level) {

      if(_vars[1]) {

	_vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
	std::cout << "Adding lt constraint" << std::endl;
#endif
      
	precedenceEncoder(_vars[0], _vars[1], 1, _solver);

      } else {

#ifdef _DEBUGWRAP
	std::cout << "Adding lt constraint" << std::endl;
#endif
	
	// x0 < r
	if(_rhs <= _vars[0]->getmin()) {
	  std::cerr << "model is inconsistent, -exiting" <<std::endl;
	  exit(1);
	} else if(_rhs <= _vars[0]->getmax()) {
	  lits.clear();
	  lits.push_back(_vars[0]->less_or_equal(_rhs-1));
	  _solver->addClause(lits);
	} 
      }

    } else {

#ifdef _DEBUGWRAP
	std::cout << "Adding lt predicate" << std::endl;
#endif	

      if(_vars[1]) {

	domain->encode(_solver);
	_vars[1] = _vars[1]->add(_solver, false);

	precedenceEncoder(_vars[0], _vars[1], this, 1, _solver);

      }
    } 

    _solver->validate();
  }

  return this;
}


/* Gt object */

Walksat_gt::Walksat_gt(Walksat_Expression *var1, Walksat_Expression *var2)
  : Walksat_binop(var1,var2) {
  domain = new DomainEncoding(this);
}

Walksat_gt::Walksat_gt(Walksat_Expression *var1, int rhs)
  : Walksat_binop(var1,rhs) {
  domain = new LeqDomain(this, var1->domain, rhs, 0);
}

Walksat_gt::~Walksat_gt() {
}

Walksat_Expression* Walksat_gt::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

    std::vector<Lit> lits;

    _vars[0] = _vars[0]->add(_solver, false);
    if(top_level) {

      if(_vars[1]) {
	_vars[1] = _vars[1]->add(_solver, false);

#ifdef _DEBUGWRAP
	std::cout << "Adding gt constraint" << std::endl;
#endif

	precedenceEncoder(_vars[1], _vars[0], 1, _solver);

      } else {

#ifdef _DEBUGWRAP
      std::cout << "Adding gt constraint" << std::endl;
#endif
	
	// x0 > r
	if(_rhs >= _vars[0]->getmax()) {
	  std::cerr << "model is inconsistent, -exiting" <<std::endl;
	  exit(1);
	} else if(_rhs >= _vars[0]->getmin()) {
	  lits.clear();
	  lits.push_back(_vars[0]->greater_than(_rhs));
	  _solver->addClause(lits);
	}
      }

    } else {

#ifdef _DEBUGWRAP
	std::cout << "Adding gt predicate" << std::endl;
#endif	

      if(_vars[1]) {
	domain->encode(_solver);
	_vars[1] = _vars[1]->add(_solver, false);

	precedenceEncoder(_vars[1], _vars[0], this, 1, _solver);

      }
    } 

    _solver->validate();
  }

  return this;
}




/* Minimise object */

Walksat_Minimise::Walksat_Minimise(Walksat_Expression *var)
  : Walksat_Expression() {
  _obj = var;
}

Walksat_Minimise::~Walksat_Minimise(){
}

Walksat_Expression* Walksat_Minimise::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

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

Walksat_Maximise::Walksat_Maximise(Walksat_Expression *var)
  : Walksat_Expression() {
  _obj = var;
}

Walksat_Maximise::~Walksat_Maximise(){
}

Walksat_Expression* Walksat_Maximise::add(WalksatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this, false);

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

WalksatSolver::WalksatSolver() 
{

#ifdef _DEBUGWRAP
  std::cout << "create a minisat solver" << std::endl;
#endif

  ////////////// Walksat Specific ////////////////
  STARTTIME = cpuTime();
  nbSolutions = 0;
  ////////////// Walksat Specific ////////////////

  minimise_obj = NULL;
  maximise_obj = NULL;
  cp_model = NULL;

  Lit dummy(0,true);
  Lit_True = Lit(dummy);
  Lit_False = ~Lit(dummy);

  current = 0;
  _atom_to_domain.push_back(NULL);
  _atom_to_type.push_back(0);
}

WalksatSolver::~WalksatSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "delete wrapped solver" << std::endl;
#endif

}

int WalksatSolver::declare(Walksat_Expression *exp, bool type) {
  int id = _expressions.size(); 
  _expressions.push_back(exp); 
  return id;
}

int WalksatSolver::create_atom(DomainEncoding* dom, const int type) {
  unsigned int id = _atom_to_domain.size(); 
  _atom_to_domain.push_back(dom); 
  _atom_to_type.push_back(type); 
  return id; 
}

void WalksatSolver::addClause(std::vector<Lit>& cl) { 
  std::vector<Lit> clause;
  if(!processClause(cl,clause))
    clause_base.push_back(clause);
}

void WalksatSolver::validate() {
}

void WalksatSolver::store_solution() {
#ifdef _DEBUGWRAP
  std::cout << "store a new solution" << std::endl;
#endif
  
  ++nbSolutions;
  if(!cp_model) cp_model = new int[_expressions.size()];
  for(unsigned int i=0; i<_variables.size(); ++i) {
    if(_variables[i]) {
      cp_model[_variables[i]->_ident] = _variables[i]->get_min();
    } else
      cp_model[_variables[i]->_ident] = 0;
  }
}

void WalksatSolver::add(Walksat_Expression* arg)
{

#ifdef _DEBUGWRAP
  std::cout << "add an expression to the solver" << std::endl;
#endif

  if(arg != NULL)
    arg->add(this, true);
}

void WalksatSolver::initialise(WalksatExpArray& arg)
{

  initialise();

}

void WalksatSolver::initialise()
{
  int i, j;
  
  int *storeptr = NULL;
  int freestore;
  int lit;
  //int simplified=0;

  std::vector<Lit> clause;

#ifdef _DEBUGWRAP
  std::cout << "initialise the solver" << std::endl;
#endif

//   for(i=0; (unsigned int)i<clause_base.size(); ++i)
//     {
//       displayClause(clause_base[i]);
//     }
  
  wsat.numatom = _atom_to_domain.size()-1;
  wsat.numclause = clause_base.size();
  
#ifdef DYNAMIC
  wsat.clause = (int **) malloc(sizeof(int *)*(wsat.numclause+1));
  wsat.size = (int *) malloc(sizeof(int)*(wsat.numclause+1));
  wsat.falsified = (int *) malloc(sizeof(int)*(wsat.numclause+1));
  wsat.lowfalse = (int *) malloc(sizeof(int)*(wsat.numclause+1));
  wsat.wherefalse = (int *) malloc(sizeof(int)*(wsat.numclause+1));
  wsat.numtruelit = (int *) malloc(sizeof(int)*(wsat.numclause+1));
#else
  if(wsat.numclause > MAXCLAUSE)                     
    {                                      
      fprintf(stderr,"ERROR - too many clauses\n"); 
      exit(-1);                              
    }                                        
#endif

  freestore = 0;
  wsat.numliterals = 0;
  for(i = 0;i < 2*MAXATOM+1;i++)
    wsat.numoccurence[i] = 0;


  for(i = 0;i < wsat.numclause;i++)
    {
      //if(!processClause(clause_base[i],clause)) {
      wsat.size[i] = -1;
      if (freestore < MAXLENGTH)
	{
	  storeptr = (int *) malloc( sizeof(int) * STOREBLOCK );
	  freestore = STOREBLOCK;
	  fprintf(stderr,"allocating memory...\n");
	}
      wsat.clause[i] = storeptr;
      if(clause_base[i].size() > MAXLENGTH)
	{
	  printf("ERROR - clause too long\n");
	  exit(-1);
	}
      else 
	{
	  wsat.size[i] = clause_base[i].size();
	  for(j = 0;j<wsat.size[i];j++)
	    {
	      lit = (sign(clause_base[i][j]) ? -1 : 1)*var(clause_base[i][j]);
	      *(storeptr++) = lit;
	      freestore--;
	      wsat.numliterals++;
	      wsat.numoccurence[lit+MAXATOM]++;
	    }
	}
      //} else ++simplified;
    }
  //wsat.numclause -= simplified;
  
  for(i = 0;i < 2*MAXATOM+1;i++)
    {
      if (freestore < wsat.numoccurence[i])
	{
	  storeptr = (int *) malloc( sizeof(int) * STOREBLOCK );
	  freestore = STOREBLOCK;
	  fprintf(stderr,"allocating memory...\n");
	}
      wsat.occurence[i] = storeptr;
      freestore -= wsat.numoccurence[i];
      storeptr += wsat.numoccurence[i];
      wsat.numoccurence[i] = 0;
    }
  
  for(i = 0;i < wsat.numclause;i++)
    {
      for(j = 0;j < wsat.size[i];j++)
	{
	  wsat.occurence[wsat.clause[i][j]+MAXATOM]
	    [wsat.numoccurence[wsat.clause[i][j]+MAXATOM]] = i;
	  wsat.numoccurence[wsat.clause[i][j]+MAXATOM]++;
	}
    } 
}

lbool WalksatSolver::truth_value(Lit x)
{
  return (wsat.solution[var(x)] == 1 ? l_True : l_False);
}

int WalksatSolver::solveAndRestart(const int policy, 
				   const unsigned int base, 
				   const double factor,
				   const double decay)
{
  return solve();
}

int WalksatSolver::solve()
{

#ifdef _DEBUGWRAP
  std::cout << "call solve" << std::endl;  
#endif 

#ifdef _DEBUGWRAP
  std::cout << "solve!" << std::endl;  
#endif 

  STARTTIME = cpuTime();

  wsat.initialize_statistics();
  wsat.print_statistics_header();
  wsat.walk_solve();

  return is_sat();
}

bool WalksatSolver::propagate()
{
  return true;
}

void WalksatSolver::reset(bool full) {
}

bool WalksatSolver::undo(const int nlevel)
{
  return true;
}

bool WalksatSolver::branch_right()
{
  return true;
}

void WalksatSolver::deduce()
{
}

void WalksatSolver::save()
{
}

void WalksatSolver::post(const char* op, Walksat_Expression* x, int v)
{
}

int WalksatSolver::startNewSearch()
{
  std::cout << "start a new interuptable search" << std::endl;
  return 0;
}

int WalksatSolver::getNextSolution()
{
  std::cout << "seek next solution" << std::endl;
  return 0;
}

int WalksatSolver::sacPreprocess(const int type)
{
  std::cout << "enforces singleton arc consistency" << std::endl;
  return 0;
}

void WalksatSolver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand)
{
  std::cout << "set the variable/value ordering" << std::endl;
}

void WalksatSolver::setFailureLimit(const int cutoff)
{
}

void WalksatSolver::setNodeLimit(const int cutoff)
{
}

void WalksatSolver::setTimeLimit(const double cutoff)
{
}

void WalksatSolver::setVerbosity(const int degree)
{
  std::cout << "set the verbosity" << std::endl;
}

void WalksatSolver::setRandomized(const int degree)
{
  std::cout << "set the type of randomization" << std::endl;
}

void WalksatSolver::setRandomSeed(const int seed)
{
  //setRandomSeed((double)seed);
}

bool WalksatSolver::is_sat()
{
  return wsat.numsuccesstry > 0;
}

bool WalksatSolver::is_unsat()
{
  return false;
}

bool WalksatSolver::is_opt()
{
  return false;
}

void WalksatSolver::printStatistics()
{
  wsat.print_statistics_final();
}

int WalksatSolver::getBacktracks()
{
  std::cout << "print the number of backtracks" << std::endl;
  return 0;
}

int WalksatSolver::getNodes()
{
  std::cout << "print the number of nodes" << std::endl;
  return 0;
}

int WalksatSolver::getFailures()
{
  std::cout << "print the number of failures" << std::endl;
  return 0;
}

int WalksatSolver::getChecks()
{
  std::cout << "print the number of checks" << std::endl;
  return 0;
}

int WalksatSolver::getPropags()
{
  std::cout << "print the number of propags" << std::endl;
  return 0;
}

double WalksatSolver::getTime()
{
  return cpuTime() - STARTTIME;
}

// int WalksatSolver::solveDimacs(const char *filename) {
//   return solve();
// }

void WalksatSolver::displayClause(std::vector<Lit>& cl) {
  if(cl.size()) {
    for(unsigned int i=0; i<cl.size(); ++i) {
      std::cout << " " ;
      if(sign(cl[i])) std::cout << "~";
      std::cout << (var(cl[i])+1);
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

void WalksatSolver::displayLiteral(Lit p) { 
  int x = var(p); 
  if(x)
    _atom_to_domain[x]->print_lit(p,_atom_to_type[x]); 
  else std::cout << "false" ;
}


