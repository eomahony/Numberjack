
#include <iostream>
#include "MiniSat.hpp"
#include <math.h>
#include <Sort.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

Lit Lit_False;
Lit Lit_True;


static inline double cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }

void printStats(Solver& S)
{
    double   cpu_time = cpuTime();
    uint64_t mem_used = 0;//memUsed();
    reportf("restarts              : %lld\n", (unsigned long long int)(S.starts));
    reportf("conflicts             : %-12lld   (%.0f /sec)\n", (unsigned long long int)(S.conflicts)   , S.conflicts   /cpu_time);
    reportf("decisions             : %-12lld   (%4.2f %% random) (%.0f /sec)\n", (unsigned long long int)(S.decisions), (float)S.rnd_decisions*100 / (float)S.decisions, S.decisions   /cpu_time);
    reportf("propagations          : %-12lld   (%.0f /sec)\n", (unsigned long long int)(S.propagations), S.propagations/cpu_time);
    reportf("conflict literals     : %-12lld   (%4.2f %% deleted)\n", (unsigned long long int)(S.tot_literals), (S.max_literals - S.tot_literals)*100 / (double)S.max_literals);
    if (mem_used != 0) reportf("Memory used           : %.2f MB\n", mem_used / 1048576.0);
    reportf("CPU time              : %g s\n", cpu_time);
}

#define CHUNK_LIMIT 1048576

class StreamBuffer {
    gzFile  in;
    char    buf[CHUNK_LIMIT];
    int     pos;
    int     size;

    void assureLookahead() {
        if (pos >= size) {
            pos  = 0;
            size = gzread(in, buf, sizeof(buf)); } }

public:
    StreamBuffer(gzFile i) : in(i), pos(0), size(0) {
        assureLookahead(); }

    int  operator *  () { return (pos >= size) ? EOF : buf[pos]; }
    void operator ++ () { pos++; assureLookahead(); }
};


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<class B>
static void skipWhitespace(B& in) {
    while ((*in >= 9 && *in <= 13) || *in == 32)
        ++in; }

template<class B>
static void skipLine(B& in) {
    for (;;){
        if (*in == EOF || *in == '\0') return;
        if (*in == '\n') { ++in; return; }
        ++in; } }

template<class B>
static int parseInt(B& in) {
    int     val = 0;
    bool    neg = false;
    skipWhitespace(in);
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '0' || *in > '9') reportf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val; }

template<class B>
static void readClause(B& in, SimpSolver& S, vec<Lit>& lits) {
    int     parsed_lit, var;
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        while (var >= S.nVars()) S.newVar();
        lits.push( (parsed_lit > 0) ? Lit(var) : ~Lit(var) );
    }
}

template<class B>
static bool match(B& in, const char* str) {
    for (; *str != 0; ++str, ++in)
        if (*str != *in)
            return false;
    return true;
}


template<class B>
static void parse_DIMACS_main(B& in, SimpSolver& S) {
    vec<Lit> lits;
    for (;;){
        skipWhitespace(in);
        if (*in == EOF) break;
        else if (*in == 'p'){
            if (match(in, "p cnf")){
                int vars    = parseInt(in);
                int clauses = parseInt(in);
                reportf("|  Number of variables:  %-12d                                         |\n", vars);
                reportf("|  Number of clauses:    %-12d                                         |\n", clauses);

                // SATRACE'06 hack
                if (clauses > 4000000)
                    S.eliminate(true);
            }else{
                reportf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
        } else if (*in == 'c' || *in == 'p')
            skipLine(in);
        else{
            readClause(in, S, lits);
            S.addClause(lits); }
    }
}

// Inserts problem into solver.
//
static void parse_DIMACS(gzFile input_stream, SimpSolver& S) {
    StreamBuffer in(input_stream);
    parse_DIMACS_main(in, S); }


//=================================================================================================

SimpSolver* solver_ptr;
static void SIGINT_handler(int signum) {
    reportf("\n"); reportf("*** INTERRUPTED ***\n");
    printStats(*solver_ptr);
    reportf("\n"); reportf("*** INTERRUPTED ***\n");
    exit(1); }




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

DomainEncoding::DomainEncoding(MiniSat_Expression *o) : AbstractDomain(o) {
  _lower = 0;
  _upper = 1;
  _size = 2;
  _values = NULL;

  // not initialised
  _direct_encoding = -1;
  _order_encoding = -1;
}

DomainEncoding::DomainEncoding(MiniSat_Expression *o, const int nval) : AbstractDomain(o) {
  _lower = 0;
  _upper = nval-1;
  _size = nval;
  _values = NULL;

  // not initialised
  _direct_encoding = -1;
  _order_encoding = -1;
}

DomainEncoding::DomainEncoding(MiniSat_Expression *o, const int lb, const int ub) : AbstractDomain(o) {
  _lower = lb;
  _upper = ub;
  _size = ub-lb+1;
  _values = NULL;

  // not initialised
  _direct_encoding = -1;
  _order_encoding = -1;
}

DomainEncoding::DomainEncoding(MiniSat_Expression *o, MiniSatIntArray& vals) : AbstractDomain(o) {
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

void DomainEncoding::encode(MiniSatSolver *solver)
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

OffsetDomain::OffsetDomain(MiniSat_Expression *o, AbstractDomain *d, const int os) : AbstractDomain(o) {
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

FactorDomain::FactorDomain(MiniSat_Expression *o, AbstractDomain *d, const int f) : AbstractDomain(o) {
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

EqDomain::EqDomain(MiniSat_Expression *o, AbstractDomain *d, const int v, const int s) : AbstractDomain(o) {
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

LeqDomain::LeqDomain(MiniSat_Expression *o, AbstractDomain *d, const int b, const int s) : AbstractDomain(o) {
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
void additionEncoder(MiniSat_Expression *X,
		     MiniSat_Expression *Y,
		     MiniSat_Expression *Z,
		     MiniSatSolver *solver) {
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
void disequalityEncoder(MiniSat_Expression *X,
			MiniSat_Expression *Y,
			MiniSatSolver *solver) {
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
void equalityEncoder(MiniSat_Expression *X,
		     MiniSat_Expression *Y,
		     MiniSatSolver *solver) {
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
void equalityEncoder(MiniSat_Expression *X,
		     MiniSat_Expression *Y,
		     MiniSat_Expression *Z,
		     MiniSatSolver *solver,
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
void precedenceEncoder(MiniSat_Expression *X,
		       MiniSat_Expression *Y,
		       const int K,
		       MiniSatSolver *solver) {
  std::vector<Lit> lits;
  int i=0, j=0, x=0, y=0;
  
  while( i<X->getsize() )
    {
      if(i<X->getsize()) x = X->getval(i);
      if(j<Y->getsize()) y = Y->getval(j);

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
void precedenceEncoder(MiniSat_Expression *X,
		       MiniSat_Expression *Y,
		       MiniSat_Expression *Z,
		       const int K,
		       MiniSatSolver *solver) {
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

void MiniSat_Expression::initialise() { 
  _ident = -1;
  _solver = NULL;
  domain = NULL;
}

MiniSat_Expression::MiniSat_Expression() {
  initialise();
  domain = new DomainEncoding(this);
}

MiniSat_Expression::MiniSat_Expression(const int nval) {
  initialise();
  domain = new DomainEncoding(this,nval);

#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif

}

MiniSat_Expression::MiniSat_Expression(const int lb, const int ub) {
  initialise();
  domain = new DomainEncoding(this,lb,ub);

#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << lb << ".." << ub << "]" << std::endl;
#endif

}

MiniSat_Expression::MiniSat_Expression(MiniSatIntArray& vals) {
  initialise();
  domain = new DomainEncoding(this,vals);

#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

}

MiniSat_Expression::~MiniSat_Expression() {
  delete domain;
}

int MiniSat_Expression::get_size() const
{
  int i=0, domsize=0;
  for(i=0; i<getsize(); ++i)
    domsize += (_solver->truth_value(equal(getval(i),i)) != l_False);
  return domsize;
}

int MiniSat_Expression::next(const int v) const
{
  int nxt = v;
  while( ++nxt <= getmax() ) 
    if(_solver->truth_value(equal(nxt)) != l_False) break;
  if(nxt > getmax()) nxt = v;
  return nxt;
}

int MiniSat_Expression::get_min() const
{
  for(int i=0; i<getsize(); ++i) {
    if(_solver->truth_value(equal(getval(i),i)) != l_False) return getval(i);
  }
  return getmin();
}

int MiniSat_Expression::get_max() const
{
  for(int i=getsize()-1; i>=0; --i)
    if(_solver->truth_value(equal(getval(i),i)) != l_False) return getval(i);
  return getmax();
}

bool MiniSat_Expression::contain(const int v) const {
  return (_solver->truth_value(equal(v)) != l_False);
}

int MiniSat_Expression::get_value() const
{
  if(_solver->cp_model) {
    return _solver->cp_model[_ident];
  } else 
    return get_min();
}

bool MiniSat_Expression::has_been_added() const {
  return (_solver != NULL);
}

MiniSat_Expression* MiniSat_Expression::add(MiniSatSolver *solver, bool top_level) {
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


MiniSat_add::MiniSat_add(MiniSat_Expression *arg1, MiniSat_Expression *arg2)
  : MiniSat_binop(arg1, arg2) {
  int _lower = arg1->getmin()+arg2->getmin();
  int _upper = arg1->getmax()+arg2->getmax();
  domain = new DomainEncoding(this, _lower, _upper);
}

MiniSat_add::MiniSat_add(MiniSat_Expression *arg1, const int arg2)
  : MiniSat_binop(arg1, arg2) {
  domain = new OffsetDomain(this,arg1->domain, arg2);

#ifdef _DEBUGWRAP
  std::cout << "creating offset expression [" << getmin() << ".." << getmax() << "]" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_add::add(MiniSatSolver *solver, bool top_level) {
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

MiniSat_add::~MiniSat_add() {}

MiniSat_mul::MiniSat_mul(MiniSat_Expression *arg1, const int arg2)
  : MiniSat_binop(arg1, arg2) {

  domain = new FactorDomain(this,arg1->domain, arg2);

}

MiniSat_mul::MiniSat_mul(MiniSat_Expression *arg1, MiniSat_Expression *arg2)
  : MiniSat_binop(arg1, arg2) {

  std::cerr << "c NOT SUPPORTED (multiplication) - exiting" << std::endl;
  exit(1);

}

MiniSat_Expression* MiniSat_mul::add(MiniSatSolver *solver, bool top_level) {
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

MiniSat_mul::~MiniSat_mul() {}

void MiniSat_AllDiff::addVar(MiniSat_Expression* v) {
  _vars.add(v);
}

MiniSat_AllDiff::MiniSat_AllDiff( MiniSat_Expression* arg1, MiniSat_Expression* arg2 ) 
  : MiniSat_Expression() {
  addVar(arg1);
  addVar(arg2);
}

MiniSat_AllDiff::MiniSat_AllDiff( MiniSatExpArray& vars ) 
  : MiniSat_Expression() {
  _vars = vars;
}

MiniSat_AllDiff::~MiniSat_AllDiff() {
  for(unsigned int i=0; i<_clique.size(); ++i)
    delete _clique[i];
}

MiniSat_Expression* MiniSat_AllDiff::add(MiniSatSolver *solver, bool top_level) {
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

      MiniSat_Expression *exp;
      for(i=1; i<n; ++i)
	for(j=0; j<i; ++j) {
	  exp = new MiniSat_ne(_vars.get_item(j), _vars.get_item(i));
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


MiniSat_Sum::MiniSat_Sum(MiniSatExpArray& vars, 
				 MiniSatIntArray& weights, 
				 const int offset)
  : MiniSat_Expression() {
  _offset = offset;
  _vars = vars;
  _weights = weights;
  initialise();
}

MiniSat_Sum::MiniSat_Sum(MiniSat_Expression *arg1, 
				 MiniSat_Expression *arg2, 
				 MiniSatIntArray& w, 
				 const int offset)
  : MiniSat_Expression() {
  _offset = offset;
  _vars.add(arg1);
  _vars.add(arg2);
  _weights = w;
  initialise();
}

MiniSat_Sum::MiniSat_Sum(MiniSat_Expression *arg, 
				 MiniSatIntArray& w, 
				 const int offset)
  : MiniSat_Expression() {
  _self = NULL;
  _offset = offset;
  _vars.add(arg);
  _weights = w;
  initialise();
}

MiniSat_Sum::MiniSat_Sum()
  : MiniSat_Expression() {
  _offset = 0;
}

void MiniSat_Sum::initialise() {
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

MiniSat_Sum::~MiniSat_Sum() {

#ifdef _DEBUGWRAP
  std::cout << "delete sum" << std::endl;
#endif

  for(unsigned int i=0; i<_subsum.size(); ++i) {
    delete _subsum[i];
  }
}

void MiniSat_Sum::addVar(MiniSat_Expression* v) {
  _vars.add(v);
}

void MiniSat_Sum::addWeight(const int w) {
  _weights.add(w);
}

void MiniSat_Sum::set_rhs(const int k) {
  _offset = k;
}

MiniSat_Expression* MiniSat_Sum::add(MiniSatSolver *solver, bool top_level){
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
	MiniSat_Expression *exp, *exp1, *exp2;
	for(int i=0; i+1<_vars.size(); i+=2) {
	  w1 = _weights.get_item(i);
	  w2 = _weights.get_item(i+1);

	  if(w1 != 1) exp1 = new MiniSat_mul(_vars.get_item(i), w1);
	  else exp1 = _vars.get_item(i);

	  if(w2 != 1) exp2 = new MiniSat_mul(_vars.get_item(i+1), w2);
	  else exp2 = _vars.get_item(i+1);

	  exp = new MiniSat_add(exp1, exp2);
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

MiniSat_binop::MiniSat_binop(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_Expression() {
  _vars[0] = var1;
  _vars[1] = var2;

#ifdef _DEBUGWRAP
  std::cout << "creating binary operator" << std::endl;
#endif

}

MiniSat_binop::MiniSat_binop(MiniSat_Expression *var1, int rhs)
  : MiniSat_Expression() {
  _vars[0] = var1;
  _vars[1] = NULL;
  _rhs = rhs;

#ifdef _DEBUGWRAP
  std::cout << "creating binary operator" << std::endl;
#endif

}


MiniSat_binop::~MiniSat_binop() {

#ifdef _DEBUGWRAP
  std::cout << "delete binary operator" << std::endl;
#endif

}

MiniSat_or::MiniSat_or(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2) {

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  domain = new DomainEncoding(this);

}

MiniSat_or::MiniSat_or(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs) {

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  if(rhs) domain = new ConstantDomain(this,1);
  else domain = new OffsetDomain(this,var1->domain,0);
}


MiniSat_or::~MiniSat_or() {

#ifdef _DEBUGWRAP
  std::cout << "delete or" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_or::add(MiniSatSolver *solver, bool top_level) {
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



MiniSat_and::MiniSat_and(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2) {

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  domain = new DomainEncoding(this);
}

MiniSat_and::MiniSat_and(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs) {

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  if(rhs) domain = new OffsetDomain(this,var1->domain,0);
  else domain = new ConstantDomain(this,0);
}


MiniSat_and::~MiniSat_and() {

#ifdef _DEBUGWRAP
  std::cout << "delete or" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_and::add(MiniSatSolver *solver, bool top_level) {
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


MiniSat_eq::MiniSat_eq(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2) {
  domain = new DomainEncoding(this);
}

MiniSat_eq::MiniSat_eq(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs) {
  domain = new EqDomain(this,var1->domain, rhs, 1);
}


MiniSat_eq::~MiniSat_eq() {

#ifdef _DEBUGWRAP
  std::cout << "delete eq" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_eq::add(MiniSatSolver *solver, bool top_level) {
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

MiniSat_ne::MiniSat_ne(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2) {
  domain = new DomainEncoding(this);
}

MiniSat_ne::MiniSat_ne(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs) {
  domain = new EqDomain(this,var1->domain, rhs, 0);
}

MiniSat_ne::~MiniSat_ne() {

#ifdef _DEBUGWRAP
  std::cout << "delete notequal" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_ne::add(MiniSatSolver *solver, bool top_level) {
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

MiniSat_le::MiniSat_le(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2) { 
#ifdef _DEBUGWRAP
  std::cout << "Creating Le constraint" << std::endl;
#endif

  domain = new DomainEncoding(this);

}

MiniSat_le::MiniSat_le(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs) {

#ifdef _DEBUGWRAP
  std::cout << "Leq on constant constructor" << std::endl;
#endif

  domain = new LeqDomain(this,var1->domain, rhs, 1);
}


MiniSat_le::~MiniSat_le() {

#ifdef _DEBUGWRAP
  std::cout << "delete lessequal" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_le::add(MiniSatSolver *solver, bool top_level) {
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

MiniSat_ge::MiniSat_ge(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2) {
#ifdef _DEBUGWRAP
  std::cout << "Creating Ge constraint" << std::endl;
#endif

  domain = new DomainEncoding(this);
}

MiniSat_ge::MiniSat_ge(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs) {

#ifdef _DEBUGWRAP
  std::cout << "Geq on constant constructor" << std::endl;
#endif

  domain = new LeqDomain(this,var1->domain, rhs-1, 0);
}

MiniSat_ge::~MiniSat_ge() {

#ifdef _DEBUGWRAP
  std::cout << "delete greaterequal" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_ge::add(MiniSatSolver *solver, bool top_level) { 
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

MiniSat_lt::MiniSat_lt(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2) {
  domain = new DomainEncoding(this);
}

MiniSat_lt::MiniSat_lt(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs) {
  domain = new LeqDomain(this,var1->domain, rhs-1, 1);
}

MiniSat_lt::~MiniSat_lt() {

#ifdef _DEBUGWRAP
  std::cout << "delete lessthan" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_lt::add(MiniSatSolver *solver, bool top_level) {
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

MiniSat_gt::MiniSat_gt(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2) {
  domain = new DomainEncoding(this);
}

MiniSat_gt::MiniSat_gt(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs) {
  domain = new LeqDomain(this, var1->domain, rhs, 0);
}

MiniSat_gt::~MiniSat_gt() {
}

MiniSat_Expression* MiniSat_gt::add(MiniSatSolver *solver, bool top_level) {
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

MiniSat_Minimise::MiniSat_Minimise(MiniSat_Expression *var)
  : MiniSat_Expression() {
  _obj = var;
}

MiniSat_Minimise::~MiniSat_Minimise(){
}

MiniSat_Expression* MiniSat_Minimise::add(MiniSatSolver *solver, bool top_level) {
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

MiniSat_Maximise::MiniSat_Maximise(MiniSat_Expression *var)
  : MiniSat_Expression() {
  _obj = var;
}

MiniSat_Maximise::~MiniSat_Maximise(){
}

MiniSat_Expression* MiniSat_Maximise::add(MiniSatSolver *solver, bool top_level) {
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

MiniSatSolver::MiniSatSolver() : SimpSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "create a minisat solver" << std::endl;
#endif

  ////////////// MiniSat Specific ////////////////
  STARTTIME = cpuTime();
  result = l_Undef;
  nbSolutions = 0;
  
  // search stuff
  conflict_clause = NULL;
  backtrack_level = 0;
  conflictC = 0;

  first_decision_level = -1;
  last_decision = lit_Undef;
  saved_level = -1;
  ////////////// MiniSat Specific ////////////////

  minimise_obj = NULL;
  maximise_obj = NULL;
  cp_model = NULL;

  vec<Lit> lits;
  Lit dummy(Solver::newVar(),true);
  Lit_True = Lit(dummy);
  Lit_False = ~Lit(dummy);

  lits.push(dummy);
  Solver::addClause(lits);

  current = 0;
  _atom_to_domain.push_back(NULL);
  _atom_to_type.push_back(0);
}

MiniSatSolver::~MiniSatSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "delete wrapped solver" << std::endl;
#endif

}

int MiniSatSolver::declare(MiniSat_Expression *exp, bool type) {
  int id = _expressions.size(); 
  _expressions.push_back(exp); 
  return id;
}

int MiniSatSolver::create_atom(DomainEncoding* dom, const int type) {
  unsigned int id = Solver::newVar();
  assert(id == _atom_to_domain.size()); 
  _atom_to_domain.push_back(dom); 
  _atom_to_type.push_back(type); 
  return id; 
}

void MiniSatSolver::addClause(std::vector<Lit>& cl) { 
  std::vector<Lit> clause;
  if(!processClause(cl,clause))
    clause_base.push_back(clause);
}

void MiniSatSolver::validate() {
  vec<Lit> cl;
  unsigned int i;
  //std::vector<Lit> clause;
  while(current < clause_base.size()) {
    //if(!processClause(clause_base[current],clause)) {
    cl.clear();
    for(i=0; i<clause_base[current].size(); ++i)
      cl.push(clause_base[current][i]);
    //displayClause(clause_base[current]);
    Solver::addClause(cl);
    //}
    ++current;
  }
}

void MiniSatSolver::store_solution() {
  
#ifdef _DEBUGWRAP
  std::cout << "store a new solution" << std::endl;
#endif
  
  ++nbSolutions;
  if(model.size() < nVars()) {
    model.growTo(nVars());
    for (int i = 0; i < nVars(); i++) model[i] = value(i);
  }
  if(!cp_model) cp_model = new int[_expressions.size()];
  for(unsigned int i=0; i<_variables.size(); ++i) {
    if(_variables[i]) {
      cp_model[_variables[i]->_ident] = _variables[i]->get_min();
    } else
      cp_model[_variables[i]->_ident] = 0;
  }
}

void MiniSatSolver::add(MiniSat_Expression* arg)
{

#ifdef _DEBUGWRAP
  std::cout << "add an expression to the solver" << std::endl;
#endif

  if(arg != NULL)
    arg->add(this, true);
}

void MiniSatSolver::initialise(MiniSatExpArray& arg)
{

#ifdef _DEBUGWRAP
  std::cout << "initialise the solver" << std::endl;
#endif

}

void MiniSatSolver::initialise()
{

#ifdef _DEBUGWRAP
  std::cout << "initialise the solver" << std::endl;
#endif

}

lbool MiniSatSolver::truth_value(Lit x)
{  
  if(model.size() > 0) {
    return modelValue(x);
  } else {
    return value(x);
  }
}

int MiniSatSolver::solveAndRestart(const int policy, 
				   const unsigned int base, 
				   const double factor,
				   const double decay)
{
  return solve();
}

int MiniSatSolver::solve()
{

#ifdef _DEBUGWRAP
  std::cout << "call solve" << std::endl;  
#endif 

#ifdef _DEBUGWRAP
  std::cout << "print to file" << std::endl;  
  toDimacs("dimacs.out");
#endif 

#ifdef _DEBUGWRAP
  std::cout << "solve!" << std::endl;  
#endif 

  start_time = getRunTime();
  saved_level = init_level;
  if(init_level < decisionLevel())
    init_level = decisionLevel();

  solver_ptr = this;
  signal(SIGINT,SIGINT_handler);

  if(minimise_obj) {
    vec<Lit> lits;
    int objective = minimise_obj->getmax();
    
    result = SimpSolver::solve(true,true);
    
    while(result == l_True && !(limitsExpired())) {
      cancelUntil(init_level);

      store_solution();
      objective = minimise_obj->get_value()-1;
      if(objective < minimise_obj->getmin()) break;

      if(verbosity >= 0) {
	std::cout << "c  new objective: " << objective+1 << std::endl;
      }
      
      lits.clear();
      lits.push(minimise_obj->less_or_equal(objective));
      Solver::addClause(lits);

      result = SimpSolver::solve(true,true);

      if(result == l_True) {
	++objective;
      }
    }
  }
  else if(maximise_obj) {
    vec<Lit> lits;
    int objective = maximise_obj->getmin();
    
    result = SimpSolver::solve(true,true);
    
    while(result == l_True && !(limitsExpired())) {
      cancelUntil(init_level);

      store_solution();
      objective = maximise_obj->get_value()+1;
      if(objective > maximise_obj->getmax()) break;

      if(verbosity >= 0) {
	std::cout << "c  new objective: " << objective-1 << std::endl;
      }
      
      lits.clear();
      lits.push(maximise_obj->greater_than(objective-1));
      Solver::addClause(lits);
      
      result = SimpSolver::solve(true,true);
      if(result == l_True) {
	++objective;
      } 
    }
  }
  else {
    result = SimpSolver::solve(true,true);
    if(result == l_True) {
      store_solution();
    }
  }

#ifdef _DEBUGWRAP
  std::cout << "print results" << std::endl; 
  if(is_sat())
    for(int i=0; i<nVars(); ++i) {
      int res = modelValue(Lit(i)).toInt();
      std::cout << "** " << i << " = " << (res == 1) 
		<< std::endl;
    }
  else
    std::cout << "unsatisfiable" << std:: endl;
#endif 

  return is_sat();
}

bool MiniSatSolver::propagate()
{
  conflict_clause = NULL;
  conflict_clause = SimpSolver::propagate();
  if(conflict_clause) {
    // CONFLICT
    conflicts++; 
    conflictC++;
    return false;
  } 
  return true;
}

void MiniSatSolver::reset(bool full) {
  learnt_clause.clear();
  backtrack_level = init_level;
  forced_decisions.clear();

  cancelUntil(backtrack_level);

  model.clear();
  delete [] cp_model;
  cp_model = NULL;

  init_level = saved_level;
  ok = true;
}

bool MiniSatSolver::undo(const int nlevel)
{
  int okay = true;
  backtrack_level = decisionLevel()-nlevel;
  if(backtrack_level < first_decision_level) okay = false;

  learnt_clause.clear();
  
  if(backtrack_level < 0) backtrack_level = 0;

  for(int i=decisionLevel()-1; i>backtrack_level; --i) forced_decisions.pop();
  last_decision = forced_decisions.last();
  forced_decisions.pop();
  cancelUntil(backtrack_level);

  return okay;
}

bool MiniSatSolver::branch_right()
{
  if (decisionLevel() == first_decision_level) return false;

  learnt_clause.clear();
  if(conflict_clause) 
    analyze(conflict_clause, learnt_clause, backtrack_level);
  else backtrack_level = decisionLevel()-1;

  for(int i=decisionLevel()-1; i>backtrack_level; --i) forced_decisions.pop();
  last_decision = forced_decisions.last();
  forced_decisions.pop();

  cancelUntil(backtrack_level);

  if (learnt_clause.size() > 0){
    if (learnt_clause.size() == 1){
      uncheckedEnqueue(learnt_clause[0]);
    }else{
      Clause* c = Clause_new(learnt_clause, true);
      learnts.push(c);
      attachClause(*c);
      claBumpActivity(*c);
      uncheckedEnqueue(learnt_clause[0], c);
    }

    varDecayActivity();
    claDecayActivity();

  } else deduce();

  return true;
}

void MiniSatSolver::deduce()
{
  uncheckedEnqueue(~(last_decision));
}


void MiniSatSolver::save()
{
  decisions++;
  newDecisionLevel();
}

void MiniSatSolver::post(const char* op, MiniSat_Expression* x, int v)
{
  if(op[1] == 't') {
    if(op[0] == 'g') ++v;
    else --v;
  }

  int lvl = decisionLevel();
  if(first_decision_level < 0) 
    first_decision_level = lvl-1;

  //learnt_clause.clear();
  Lit next = lit_Undef;


  switch(op[0]) {
  case 'e': next =  (x->equal(v)); break;
  case 'n': next = ~(x->equal(v)); break;
  case 'g': next =  (x->greater_than(v-1)); break;
  case 'l': next =  (x->less_or_equal(v )); break;
  }

  forced_decisions.push(next);

  uncheckedEnqueue(next);
}

int MiniSatSolver::startNewSearch()
{
  std::cout << "start a new interuptable search" << std::endl;
  return 0;
}

int MiniSatSolver::getNextSolution()
{
  std::cout << "seek next solution" << std::endl;
  return 0;
}

int MiniSatSolver::sacPreprocess(const int type)
{
  std::cout << "enforces singleton arc consistency" << std::endl;
  return 0;
}

void MiniSatSolver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand)
{
  std::cout << "set the variable/value ordering" << std::endl;
}

void MiniSatSolver::setFailureLimit(const int cutoff)
{
  fail_limit = cutoff;
}

void MiniSatSolver::setNodeLimit(const int cutoff)
{
  node_limit = cutoff;
}

void MiniSatSolver::setTimeLimit(const double cutoff)
{
  time_limit = cutoff;
}

void MiniSatSolver::setVerbosity(const int degree)
{
  //std::cout << "set the verbosity" << std::endl;
  verbosity = degree-1;
}

void MiniSatSolver::setRandomized(const int degree)
{
  //std::cout << "set the type of randomization" << std::endl;
  if(degree) SimpSolver::setRandomized();
}

void MiniSatSolver::setRandomSeed(const int seed)
{
  setRandomSeed((double)seed);
}

bool MiniSatSolver::is_sat()
{
  //return (result == 1 || nbSolutions);
  //return (model.size() > 0);
  return (cp_model != NULL);
}

bool MiniSatSolver::is_unsat()
{
  return (result == 0 && nbSolutions == 0);
}

bool MiniSatSolver::is_opt()
{
  return (result == 0 && nbSolutions);
}

void MiniSatSolver::printStatistics()
{
  //std::cout << "print a bunch of statistics" << std::endl;
  //printStats(minisolver);
  printStats(*this);
}

int MiniSatSolver::getBacktracks()
{
  std::cout << "print the number of backtracks" << std::endl;
  return 0;
}

int MiniSatSolver::getNodes()
{
  //std::cout << "print the number of nodes" << std::endl;
  //return NULL;
  return decisions;
}

int MiniSatSolver::getFailures()
{
  std::cout << "print the number of failures" << std::endl;
  return 0;
}

int MiniSatSolver::getChecks()
{
  std::cout << "print the number of checks" << std::endl;
  return 0;
}

int MiniSatSolver::getPropags()
{
  std::cout << "print the number of propags" << std::endl;
  return 0;
}

double MiniSatSolver::getTime()
{
  //std::cout << "print the cpu time" << std::endl;
  //return NULL;
  return cpuTime() - STARTTIME;
}

int MiniSatSolver::solveDimacs(const char *filename) {
  gzFile in = gzopen(filename, "rb");
  //parse_DIMACS(in, minisolver);
  parse_DIMACS(in, *this);
  
  return solve();
}

void MiniSatSolver::displayClause(std::vector<Lit>& cl) {
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

void MiniSatSolver::displayLiteral(Lit p) { 
  int x = var(p); 
  if(x)
    _atom_to_domain[x]->print_lit(p,_atom_to_type[x]); 
  else std::cout << "false" ;
}


