
#include <iostream>
#include "MiniSat.hpp"
#include <math.h>
//#include <vector>
#include <Sort.h>
//#include "SimpSolver.hpp"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

Lit Lit_False;
Lit Lit_True;


void printClause(vec<Lit>& cl)
{
  for(int i=0; i<cl.size(); ++i)
    std::cout << " " << (sign(cl[i]) ? "-" : "") << var(cl[i]);
  std::cout << std::endl;
}

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
 ********************     EXPRESSION        *******************
 **************************************************************/

void MiniSat_Expression::initialise()
{ 
  _ident = -1;
  _solver = NULL;
  _direct_encoding = NULL;
  _order_encoding = NULL;
  _values = NULL;
  //_encoding_type = NO;
  _encoding_type = (ORDER | DIRECT); // we force all variables to have both, 
  // otherwise the states of the literals do not exactly reflect the state of the cp var
  // ans this is bad for the get_min/get_max/get_size
  _lower = 0;
  _upper = 0; 
  _size  = 0;  
}

MiniSat_Expression::MiniSat_Expression()
{
  initialise();

#ifdef _DEBUGWRAP
  std::cout << "creating an empty expression" << std::endl;
#endif

}

MiniSat_Expression::MiniSat_Expression(const int nval)
{
  initialise();
  _upper = nval-1;
  _size  = nval;

#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif

}

MiniSat_Expression::MiniSat_Expression(const int lb, const int ub)
{
  initialise();
  _upper = ub;
  _lower = lb;
  _size  = (ub-lb+1);

#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << lb << ".." << ub << "]" << std::endl;
#endif

}

MiniSat_Expression::MiniSat_Expression(MiniSatIntArray& vals)
{
  initialise();
  _size   = vals.size();
  _values = new int[_size]; 
  _lower = _upper = vals.get_item(0);

  std::cout << "init: " << _lower << " " << _upper << std::endl;

  for(int i=0; i<_size; ++ i) {
    _values[i] = vals.get_item(i);

    std::cout << _values[i] << " ";

    if(_upper < _values[i]) _upper = _values[i];
    if(_lower > _values[i]) _lower = _values[i];
  }
  
  std::cout << std::endl;


#ifdef _DEBUGWRAP
  std::cout << "creating variable [" << _lower << ".." << _upper << "]" << std::endl;
#endif

}

// int MiniSat_Expression::get_lower() const {
//   return _lower;
// }

// int get_upper() const {

// }


int MiniSat_Expression::getval(const int i) const {
  if(_values) return _values[i];
  else return _lower+i;
}

int MiniSat_Expression::getNextEq(const int v) const {
  if(!_values) {
    if(v <= _lower) return _lower;
    else if(v <= _upper) return v;
  } else return _values[getind(v)];
  return v-1;
}

int MiniSat_Expression::getNext(const int v) const {
  if(!_values) {
    if(v < _lower) return _lower;
    else if(v < _upper) return v+1;
  } else return _values[getind(v+1)];
  return v;
}

int MiniSat_Expression::getPrevEq(const int v) const {
  if(!_values) {
    if(v >= _upper) return _upper;
    else if(v >= _lower) return v;
  } else {
    int x = _values[getind(v)];
    if(_values[x] == v) return v;
    else if(x) return _values[x-1];
  }
  return v-1;
}

int MiniSat_Expression::getPrev(const int v) const {
  if(!_values) {
    if(v >= _upper) return _upper;
    else if(v >= _lower) return v;
  } else {
    int x = _values[getind(v)];
    if(x) return _values[x-1];
  }
  return v-1;
}

int MiniSat_Expression::getind(const int v) const {
  if(!_values) return v-_lower;
  int lb = 0, ub = _size-1, x;
  while(lb < ub) {
    x = (lb+ub)/2;
    if(_values[x] == v) return x;
    if(_values[x] < v) lb = x+1;
    else ub = x-1;
  }
  return lb;
}


void MiniSat_Expression::setDirectEncoding()
{
#ifdef _DEBUGWRAP
  std::cout << "set encoding (exp): direct" << std::endl;
#endif

  _encoding_type |= DIRECT;
}

void MiniSat_Expression::setOrderEncoding()
{
#ifdef _DEBUGWRAP
  std::cout << "set encoding (exp): order" << std::endl;
#endif

  _encoding_type |= ORDER;
}

int MiniSat_Expression::get_size() const
{
  int i=0, lb=0, ub=0, domsize=0;
  lbool tv;
  if(_direct_encoding) {
    for(i=0; i<_size; ++i) {
      tv = _solver->truth_value(equal(getval(i)));
      if(tv == l_True) {domsize=1; break;}
      if(tv != l_False) ++domsize;
    }
  } else if(_order_encoding) {
    for(i=0; i<_size; ++i) {
      lb = getval(i);
      if(_solver->truth_value(less_or_equal(lb)) != l_False) 
	break;
    }
    for(i=_size-1; i>=0; --i) {
      ub = getval(i);
      if(_solver->truth_value(greater_than(ub-1)) != l_False) 
	break;
    }
    domsize = ub-lb+1;
  } 

  return domsize;
}

int MiniSat_Sum::get_size() const
{
  if(_vars.size() > 1)
    return _self->get_size();
  else
    return MiniSat_Expression::get_size();
}


int MiniSat_Expression::get_min() const
{
  int i=0;
  if(_direct_encoding) {
    for(i=0; i<_size; ++i) {
      if(_solver->truth_value(equal(getval(i))) != l_False)
	break;
    }
  } else if(_order_encoding) {
    for(i=0; i<_size; ++i)
      if(_solver->truth_value(less_or_equal(getval(i))) != l_False)
	break;
  } 
  return getval(i);
}

int MiniSat_Sum::get_min() const
{
  if(_vars.size() > 1)
    return _self->get_min();
  else
    return MiniSat_Expression::get_min();
}

int MiniSat_Expression::get_max() const
{
  int i=0;
  if(_direct_encoding) {
    for(i=_size-1; i>=0; --i)
      if(_solver->truth_value(equal(getval(i))) != l_False)
	break;
  } else if(_order_encoding) {
    for(i=_size-1; i>=0; --i)
      if(_solver->truth_value(less_or_equal(getval(i))) != l_False)
	break;
  }

  return getval(i);
}

int MiniSat_Sum::get_max() const
{
  if(_vars.size() > 1)
    return _self->get_max();
  else
    return MiniSat_Expression::get_max();
}


bool MiniSat_Expression::contain(const int v) const {
  bool isin = false;
  if(_direct_encoding) 
    isin = _solver->truth_value(equal(v)) != l_False;
  else if(_order_encoding) 
    isin = (_solver->truth_value(less_or_equal(v)) != l_False &&
	    _solver->truth_value(greater_than(v-1)) != l_False);
  return isin;
}

int MiniSat_Expression::get_value() const
{
  if(_solver->cp_model) {
    return _solver->cp_model[_ident];
  } else 
    return get_min();
}




// int MiniSat_Expression::get_value() const
// {
//   return _solver->cp_model[_ident];
// }


MiniSat_Expression::~MiniSat_Expression()
{
  if(_direct_encoding) {
    _direct_encoding += _lower;
    delete [] _direct_encoding;
  }
  if(_order_encoding) {
    _order_encoding += _lower;
    delete [] _order_encoding;
  }
  delete [] _values;
}

void MiniSat_Expression::channel(MiniSatSolver *solver)
{
#ifdef _DEBUGWRAP
  std::cout << "\tchannel encoding" << std::endl;
#endif
  vec<Lit> lits;
  int x;
  for(int i=0; i<_size; ++i) {
    x = getval(i);
    if(x>_lower) {

      std::cout << "channel: " << _ident << " = " << x << " -> " << _ident << " > " << (getval(i-1)) << std::endl;

      lits.clear();
      lits.push(~(equal(x)));
      lits.push(greater_than(getval(i-1)));

      printClause(lits);
      
      solver->addClause(lits);
    }
    if(x<_upper) {

      std::cout << "channel: " << _ident << " = " << x << " -> " << _ident << " <= " << x << std::endl;

      lits.clear();
      lits.push(~equal(x));
      lits.push(less_or_equal(x));

      printClause(lits);

      solver->addClause(lits);
    }
  }


//   for(int i=_lower; i<=_upper; ++i) {
//     // WE SHOULDN'T NEED A QUADRATIC NUMBER OF CLAUSES
//     //       // -xij \/ -yik for all k < j
//     //       for(j=_lower; j<i; ++j) {
//     // 	lits.clear();
//     // 	lits.push(~Lit(_direct_encoding[i]));
//     // 	lits.push(~Lit(_order_encoding[j]));
//     // 	solver->addClause(lits);
//     //       }
//     //       // -xij \/ yik for all k >= j
//     //       for(j=i; j<_upper; ++j) {
//     // 	lits.clear();
//     // 	lits.push(~Lit(_direct_encoding[i]));
//     // 	lits.push(Lit(_order_encoding[j]));
//     // 	solver->addClause(lits);
//     //       }
//     if(i>_lower) {
//       lits.clear();
//       lits.push(~Lit(_direct_encoding[i]));
//       lits.push(~Lit(_order_encoding[i-1]));
//       solver->addClause(lits);
//     }
//     if(i<_upper) {
//       lits.clear();
//       lits.push(~Lit(_direct_encoding[i]));
//       lits.push(Lit(_order_encoding[i]));
//       solver->addClause(lits);
//     }
//   }
}

void MiniSat_Expression::encode(MiniSatSolver *solver)
{
#ifdef _DEBUGWRAP
  std::cout << "encode variable x" << _ident << std::endl;
#endif

  assert(_encoding_type == (DIRECT | ORDER));

 
  if(!_encoding_type)
    _encoding_type |= ORDER;


  vec<Lit> lits;
  int i, x, y, val;
  if( _encoding_type & DIRECT ) {
    if( !_direct_encoding ) {

#ifdef _DEBUGWRAP
      std::cout << "\tdirect encoding" << std::endl;
#endif
      _direct_encoding = new Var[_upper-_lower+1];
      _direct_encoding -= _lower;


      std::cout << "direct: " ;
      
      for(i=0; i<_size; ++i) {
	val = getval(i);
	_direct_encoding[val] = solver->newVar(this, (val << 1));
      }

      for(i=0; i<_size; ++i) {
	lits.push(equal(getval(i)));

	std::cout << _ident << " = " << getval(i) << " or " ;

	assert(var(lits[i]) == var(equal(getval(i))));
	assert(solver->_lit_to_var[var(lits[i])] == this);
	assert(((solver->_lit_to_val[var(lits[i])]) >> 1) == getval(i));
	assert(((solver->_lit_to_val[var(lits[i])]) % 2) == 0);
	assert(sign(lits[i]) == false);

      }

      std::cout << std::endl;

      solver->addClause(lits);
      
      printClause(lits);


//       // those are optional, but can provides stronger filtering
//       for(i=_lower; i<_upper; ++i)
// 	for(j=i+1; j<=_upper; ++j)
// 	  if(_direct_encoding[i] >= 0 && _direct_encoding[j] >= 0)
// 	    {
// 	      lits.clear();
// 	      lits.push(~Lit(_direct_encoding[i]));
// 	      lits.push(~Lit(_direct_encoding[j]));
// 	      solver->addClause(lits);
// 	    }
      if(_order_encoding)
	this->channel(_solver);
    }
  }

  if( _encoding_type & ORDER ) {
    if( !_order_encoding) {

#ifdef _DEBUGWRAP
      std::cout << "\torder encoding" << std::endl;
#endif

      _order_encoding = new Var[_upper-_lower];
      _order_encoding -= _lower;
      
      for(i=0; i<_size-1; ++i) {
	val = getval(i);
	_order_encoding[val] = solver->newVar(this, (val << 1)+1);
      }
      
      for(i=1; i<_size-1; ++i) {

	x = getval(i);
	y = getval(i-1);


	std::cout << "order: " << _ident << " > " << y << " or " << _ident << " <= " << x << std::endl;

	lits.clear();
	lits.push(greater_than(y));
	lits.push(less_or_equal(x));

	printClause(lits);

	solver->addClause(lits);


// 	assert(solver->_lit_to_var[var(lits[1])] == this);
// 	assert(((solver->_lit_to_val[var(lits[1])]) >> 1) == x);
// 	assert(((solver->_lit_to_val[var(lits[1])]) % 2) == 1);
// 	assert(sign(lits[1]) == false);

// 	assert(solver->_lit_to_var[var(lits[0])] == this);
// 	assert(((solver->_lit_to_val[var(lits[0])]) >> 1) == y);
// 	assert(((solver->_lit_to_val[var(lits[0])]) % 2) == 1);
// 	assert(sign(lits[0]) == true);
      }

      if(_direct_encoding)
	this->channel(_solver);
    }
  } 

  //std::cout << "end encode" << std::endl;
}

Lit MiniSat_Expression::less_or_equal(const int value) const {
  if(value < _lower) {
    return Lit_False;
  } else if(value >= _upper) {
    return Lit_True;
  } else {
    if(!_order_encoding) {
      std::cerr << "Warning: order encoding not defined" << std::endl;
      exit(1);
    } else {
      //int v=value;
      //while(v>=_lower && _order_encoding[v] < 0) --v;
      return Lit(_order_encoding[value]);
    }
  }
}

Lit MiniSat_Expression::greater_than(const int value) const {
  if(value < _lower) {
    return Lit_True;
  } else if(value >= _upper) {
    return Lit_False;
  } else {
    if(!_order_encoding) {
      std::cerr << "Warning: order encoding not defined" << std::endl;
      exit(1);
    } else {
      //int v=value;
      //while(v<_upper && _order_encoding[v] < 0) ++v;
      return ~Lit(_order_encoding[value]);
    }
  }
}

Lit MiniSat_Expression::equal(const int value) const {
  if(value < _lower || value > _upper) {
    return Lit_False;
  } else {
    if(!_direct_encoding) {
      std::cerr << "Warning: direct encoding not defined" << std::endl;
      exit(1);
    } else {
      if(_direct_encoding[value] < 0) {
	return Lit_False;
      } else {
	return Lit(_direct_encoding[value]);
      }
    }
  }
}

// Lit MiniSat_Expression::not_equal(const int value) {
//   if(value < _lower || value > _upper) {
//     return Lit_True;
//   } else {
//     if(!_direct_encoding) {
//       std::cerr << "Warning: direct encoding not defined" << std::endl;
//       exit(1);
//     } else {
//       if(_direct_encoding[value] < 0) {
// 	return Lit_True;
//       } else {
// 	return ~Lit(_order_encoding[value]);
//       }
//     }
//   }
// }

bool MiniSat_Expression::has_been_added() const
{
  return (_solver != NULL);
}

MiniSat_Expression* MiniSat_Expression::add(MiniSatSolver *solver, bool top_level){

#ifdef _DEBUGWRAP
  std::cout << "+ add variable [" << _lower << ".." << _upper << "]" << std::endl;
#endif

  if(!has_been_added()) {

#ifdef _DEBUGWRAP
    std::cout << "+ add variable to minisat" << std::endl;    
#endif

    _solver = solver;
    _ident = _solver->declare(this);

  }
  encode(_solver);

  return this;
}

void MiniSat_binop::setDirectEncoding()
{
#ifdef _DEBUGWRAP
  std::cout << "set encoding (bin op): direct" << std::endl;
#endif

  _encoding_type |= DIRECT;

  if(!_vars[1]) // for views
    _vars[0]->setDirectEncoding();
}

void MiniSat_binop::setOrderEncoding()
{
#ifdef _DEBUGWRAP
  std::cout << "set encoding (bin op): order" << std::endl;
#endif

  _encoding_type |= ORDER;

  if(!_vars[1]) // for views
    _vars[0]->setOrderEncoding();
}


int MiniSat_add::get_min() const {
  return this->MiniSat_Expression::get_min();
}

int MiniSat_add::getval(const int i) const {
  if(_vars[1]) return this->MiniSat_Expression::getval(i);
  return (_vars[0]->getval(i)+_rhs);
}

Lit MiniSat_add::less_or_equal(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::less_or_equal(value);
  else return _vars[0]->less_or_equal(value-_rhs);
}

Lit MiniSat_add::greater_than(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::greater_than(value);
  else return _vars[0]->greater_than(value-_rhs);
}

Lit MiniSat_add::equal(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::equal(value);
  else return _vars[0]->equal(value-_rhs);
}

// Lit MiniSat_add::not_equal(const int value) {
//   if(_vars[1]) return this->MiniSat_Expression::not_equal(value);
//   else return _vars[0]->not_equal(value-_rhs);
// }

MiniSat_add::MiniSat_add(MiniSat_Expression *arg1, MiniSat_Expression *arg2)
  : MiniSat_binop(arg1, arg2)
{
  initialise();
  _lower = arg1->_lower+arg2->_lower;
  _upper = arg1->_upper+arg2->_upper;
  _size = _upper-_lower+1;

#ifdef _DEBUGWRAP
  std::cout << "creating add expression [" << _lower << ".." << _upper << "]" << std::endl;
#endif

  arg1->setOrderEncoding();
  arg2->setOrderEncoding();
  this->setOrderEncoding();

}

MiniSat_add::MiniSat_add(MiniSat_Expression *arg1, const int arg2)
  : MiniSat_binop(arg1, arg2)
{
  initialise();
  _lower = arg1->_lower+arg2;
  _upper = arg1->_upper+arg2;
  _size = arg1->_size;

#ifdef _DEBUGWRAP
  std::cout << "creating offset expression [" << _lower << ".." << _upper << "]" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_add::add(MiniSatSolver *solver, bool top_level) {
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    if(top_level) {
      
      std::cerr << "add predicate at the top level not supported" << std::endl;
      exit(1);
      
    } else {

      _vars[0] = _vars[0]->add(solver, false);            

      if(_vars[1]) {

	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	this->encode(solver);

#ifdef _DEBUGWRAP
      std::cout << "add add predicate" << std::endl;
#endif

	vec<Lit> lits;
	int i, j, x, y;
	Lit p;
	bool valid;

	for(i=0; i<_vars[0]->_size; ++i) {
	  x = _vars[0]->getval(i);
	  for(j=0; j<_vars[1]->_size; ++j) {
	    y = _vars[1]->getval(j);

	    lits.clear();
	    valid = false;

	    p = _vars[0]->less_or_equal(x-1);
	    if(p != Lit_False) {
	      lits.push(p);
	      valid |= (p == Lit_True);
	    }

	    p = _vars[1]->less_or_equal(y-1);
	    if(p != Lit_False) {
	      lits.push(p);
	      valid |= (p == Lit_True);
	    }

	    p = this->greater_than(x+y-1);
	    if(p != Lit_False) {
	      lits.push(p);
	      valid |= (p == Lit_True);
	    }

	    if(!valid) {
	      solver->addClause(lits);
	    }

	    valid = false;
	    lits.clear();

	    p = _vars[0]->greater_than(x);
	    if(p != Lit_False) {
	      lits.push(p);
	      valid |= (p == Lit_True);
	    }

	    p = _vars[1]->greater_than(y);
	    if(p != Lit_False) {
	      lits.push(p);
	      valid |= (p == Lit_True);
	    }

	    p = this->less_or_equal(x+y);
	    if(p != Lit_False) {
	      lits.push(p);
	      valid |= (p == Lit_True);
	    }

	    if(!valid) {
	      solver->addClause(lits);
	    }
	    
	  }
	}
      }
    }
  } 

  return this;
}

MiniSat_add::~MiniSat_add() {}

int MiniSat_mul::get_min() const {
  return this->MiniSat_Expression::get_min();
}

int MiniSat_mul::getval(const int i) const {
  if(_vars[1]) return this->MiniSat_Expression::getval(i);
  return (_vars[0]->getval(i)*_rhs);
}

Lit MiniSat_mul::less_or_equal(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::less_or_equal(value);
  else {
    
    //        std::cout << _rhs << " * X <= " << value << " <-> X ";

    if(_rhs>=0) {

//       std::cout << "<= " << ((value/_rhs)// -(value%_rhs != 0 //&& value<0
// // 					   )
// 			     ) << std::endl;
      return _vars[0]->less_or_equal((value/_rhs)// -(value%_rhs != 0 //&& value<0
// 						   )
				     );

    } else {
      
//       std::cout << ">= " << ((value/_rhs)+(value%_rhs != 0//  && value>0
// 					   )) << std::endl;
      
      return _vars[0]->greater_than((value/_rhs)-1+(value%_rhs != 0//  && value>0
						  ));
    }
  }
}

Lit MiniSat_mul::greater_than(const int value) const {
  if(!_vars[1]) {

//     std::cout << _rhs << " * X > " << value << " <-> X " ;

    if(_rhs>=0) {

//       std::cout << "> " << ((value/_rhs)) << std::endl;

      return _vars[0]->greater_than((value/_rhs));
    } else {

//       std::cout << "< " << ((value/_rhs)+(value%_rhs != 0)) << std::endl;

      return _vars[0]->less_or_equal((value/_rhs)-1+(value%_rhs != 0));
    }
  }
  return this->MiniSat_Expression::greater_than(value);
}

Lit MiniSat_mul::equal(const int value) const {
  if(!_vars[1]) {
    if(value%_rhs) return Lit_False;
    return _vars[0]->equal(value/_rhs);
  }
  return this->MiniSat_Expression::equal(value);
}

// Var MiniSat_mul::getDirectVar(const int value) {
//   if(_vars[1]) return MiniSat_Expression::getDirectVar(value);
//   else {

//     std::cout <<value<< "%" << _rhs << " = "<< (value%_rhs) << std::endl;

//     if(value%_rhs) return -1;
//     return _vars[0]->getDirectVar(value/_rhs);
//   }
// }

// Var MiniSat_mul::getInterVar(const int value) {
//   if(_vars[1]) return MiniSat_Expression::getInterVar(value);
//   else {

//     std::cout <<value<< "%" << _rhs << " = "<< (value%_rhs) << std::endl;

//     if(value%_rhs) return -1;
//     return _vars[0]->getInterVar(value/_rhs);
//   }
// }

MiniSat_mul::MiniSat_mul(MiniSat_Expression *arg1, const int arg2)
  : MiniSat_binop(arg1, arg2)
{
  initialise();
  if(arg2 > 0) {
    _lower = arg1->_lower*arg2;
    _upper = arg1->_upper*arg2;
  } else {
    _lower = arg1->_upper*arg2;
    _upper = arg1->_lower*arg2;
  }
  _size = arg1->_size;

#ifdef _DEBUGWRAP
  std::cout << "creating factor expression [" << _lower << ".." << _upper << "]" << std::endl;
#endif

}

MiniSat_mul::MiniSat_mul(MiniSat_Expression *arg1, MiniSat_Expression *arg2)
  : MiniSat_binop(arg1, arg2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating mul expression [" << _lower << ".." << _upper << "]" << std::endl;
#endif

  std::cerr << "c NOT SUPPORTED (multiplication) - exiting" << std::endl;
  exit(1);

}

MiniSat_Expression* MiniSat_mul::add(MiniSatSolver *solver, bool top_level) {
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    if(top_level) {
      
      std::cerr << "mul predicate at the top level not supported" << std::endl;
      exit(1);
      
    } else {
      
#ifdef _DEBUGWRAP
      std::cout << "add mul predicate" << std::endl;
#endif
      
      _vars[0] = _vars[0]->add(solver, false);

      if(_vars[1]) {

	std::cerr << "c NOT SUPPORTED (multiplication) - exiting" << std::endl;
	exit(1);

      } 
    }
  } 

  return this;
}

MiniSat_mul::~MiniSat_mul() {}



void MiniSat_AllDiff::addVar(MiniSat_Expression* v){
  _vars.add(v);
}

MiniSat_AllDiff::MiniSat_AllDiff( MiniSat_Expression* arg1, MiniSat_Expression* arg2 ) 
  : MiniSat_Expression()
{
  addVar(arg1);
  addVar(arg2);

#ifdef _DEBUGWRAP
  std::cout << "creating bin alldiff" << std::endl;
#endif

  arg1->setDirectEncoding();
  arg2->setDirectEncoding();

}

MiniSat_AllDiff::MiniSat_AllDiff( MiniSatExpArray& vars ) 
  : MiniSat_Expression() 
{

  _vars = vars;

#ifdef _DEBUGWRAP
  std::cout << "creating alldiff" << std::endl;
#endif

  int i, n=_vars.size();  
  for(i=0; i<n; ++i)
    (_vars.get_item(i))->setDirectEncoding();

}

MiniSat_AllDiff::~MiniSat_AllDiff()
{
  for(unsigned int i=0; i<_clique.size(); ++i)
    delete _clique[i];
}

MiniSat_Expression* MiniSat_AllDiff::add(MiniSatSolver *solver, bool top_level) {
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    if(top_level){
      
      int i, j, n=_vars.size();  
      for(i=0; i<n; ++i) 
	_vars.set_item(i, (_vars.get_item(i))->add(_solver,false));

#ifdef _DEBUGWRAP
      std::cout << "add alldiff constraint" << std::endl;
#endif

      MiniSat_Expression *exp;
      for(i=1; i<n; ++i)
	for(j=0; j<i; ++j) {
	  exp = new MiniSat_ne(_vars.get_item(j), _vars.get_item(i));
	  solver->add(exp);
	  _clique.push_back(exp);
	}
      
    } else {

      std::cerr << "Alldiff not in top level insupported at the moment" << std::endl;
      exit(1);

    }
  }
  
  return NULL;
}

Lit MiniSat_Sum::less_or_equal(const int value) const {
  int factor = _weights.get_item(0);  
  if(_vars.size() > 1) return _self->MiniSat_Expression::less_or_equal(value);
  else return _vars.get_item(0)->less_or_equal(((value-_offset)/factor)+((value-_offset)%factor != 0));
}

Lit MiniSat_Sum::greater_than(const int value) const {
  int factor = _weights.get_item(0);  
  if(_vars.size() > 1) return _self->MiniSat_Expression::greater_than(value);
  else return _vars.get_item(0)->greater_than(((value-_offset)/factor)+((value-_offset)%factor != 0));
}

Lit MiniSat_Sum::equal(const int value) const {
  if(_vars.size() > 1) return _self->MiniSat_Expression::equal(value);
  else {
    int factor = _weights.get_item(0);
    if((value-_offset)%factor) return Lit_False;
    return _vars.get_item(0)->equal((value-_offset)/factor);
  }
}

MiniSat_Sum::MiniSat_Sum(MiniSatExpArray& vars, 
			 MiniSatIntArray& weights, 
			 const int offset)
  : MiniSat_Expression() 
{
  _self = NULL;
  _offset = offset;
  _vars = vars;
  _weights = weights;
  initialise();
}

MiniSat_Sum::MiniSat_Sum(MiniSat_Expression *arg1, 
			 MiniSat_Expression *arg2, 
			 MiniSatIntArray& w, 
			 const int offset)
  : MiniSat_Expression() 
{
  _self = NULL;
  _offset = offset;
  _vars.add(arg1);
  _vars.add(arg2);
  _weights = w;
  initialise();
}

MiniSat_Sum::MiniSat_Sum(MiniSat_Expression *arg, 
			 MiniSatIntArray& w, 
			 const int offset)
  : MiniSat_Expression() 
{
  _self = NULL;
  _offset = offset;
  _vars.add(arg);
  _weights = w;
  initialise();
}

MiniSat_Sum::MiniSat_Sum()
  : MiniSat_Expression()
{
  _offset = 0;
}

void MiniSat_Sum::initialise() {

#ifdef _DEBUGWRAP
  std::cout << "creating sum: Size of parameters is " << _vars.size() << std::endl; 
#endif
    
  _lower = 0;
  _upper = 0;

  int weight;
  for(int i = 0; i < _vars.size(); ++i){
    weight = _weights.get_item(i);
    
    if( weight > 0 ) {
      _lower += (weight * _vars.get_item(i)->_lower);
      _upper += (weight * _vars.get_item(i)->_upper);
    } else {
      _upper += (weight * _vars.get_item(i)->_lower);
      _lower += (weight * _vars.get_item(i)->_upper);
    }

    _vars.get_item(i)->setDirectEncoding();
  }
  _lower += _offset;
  _upper += _offset;

  if(_vars.size() > 1)
    _size = (_upper - _lower +1);
  else
    _size = _vars.get_item(0)->_size;

#ifdef _DEBUGWRAP
  std::cout << "Intermediate variable has values: " << _lower << " to " << _upper << std::endl;
#endif

}

MiniSat_Sum::~MiniSat_Sum(){

#ifdef _DEBUGWRAP
  std::cout << "delete sum" << std::endl;
#endif

//   for(unsigned int i=0; i<_subsum.size()-1; ++i) {
    //     delete _subsum[i];
//   }

}

void MiniSat_Sum::addVar(MiniSat_Expression* v){
  _vars.add(v);
}

void MiniSat_Sum::addWeight(const int w){
  _weights.add(w);
}

void MiniSat_Sum::set_rhs(const int k){
  _offset = k;
}

MiniSat_Expression* MiniSat_Sum::add(MiniSatSolver *solver, bool top_level){
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);
    
    if(!top_level) {
      
      if(_vars.size() == 1) {

	_vars.set_item(0, _vars.get_item(0)->add(solver, false));
	_self = this;

#ifdef _DEBUGWRAP    
      std::cout << "add sum" << std::endl;
#endif

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
	  if(_encoding_type & DIRECT) exp->setDirectEncoding();
	  exp->add(solver, false);
	  _vars.add(exp);
	  _weights.add(1);
	  _subsum.push_back(exp);
	} 

#ifdef _DEBUGWRAP    
      std::cout << "add sum predicate" << std::endl;
#endif
	
	_self = _vars.get_item(_vars.size()-1);
	//_solver->_variables[_ident] = _self;
	_solver->_variables[_ident] = NULL;
      }

      
    } else {
      std::cout << "Warning SUM constraint on top level not supported" << std::endl;
    }
  }

  return _self;
}


/* Binary operators */

MiniSat_binop::MiniSat_binop(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_Expression()
{
  _vars[0] = var1;
  _vars[1] = var2;

#ifdef _DEBUGWRAP
  std::cout << "creating binary operator" << std::endl;
#endif

}

MiniSat_binop::MiniSat_binop(MiniSat_Expression *var1, int rhs)
  : MiniSat_Expression()
{
  _vars[0] = var1;
  _vars[1] = NULL;
  _rhs = rhs;

#ifdef _DEBUGWRAP
  std::cout << "creating binary operator" << std::endl;
#endif

}


MiniSat_binop::~MiniSat_binop(){

#ifdef _DEBUGWRAP
  std::cout << "delete binary operator" << std::endl;
#endif

}

MiniSat_or::MiniSat_or(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setDirectEncoding();
  _vars[1]->setDirectEncoding();

}

MiniSat_or::MiniSat_or(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs)
{

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setDirectEncoding();

}


MiniSat_or::~MiniSat_or(){

#ifdef _DEBUGWRAP
  std::cout << "delete or" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_or::add(MiniSatSolver *solver, bool top_level){
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    vec<Lit> lits;
    if(top_level) {
            
      if(_vars[1]) {
  
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	

#ifdef _DEBUGWRAP
	std::cout << "add or constraint" << std::endl;
#endif

	lits.clear();
	for(int i=0; i<2; ++i) 
	  lits.push(~(_vars[i]->equal(0)));
	solver->addClause(lits);
      
      } else if(_rhs == 0) {
	
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "add or constraint" << std::endl;
#endif

	lits.clear();
	lits.push(~(_vars[0]->equal(0)));
	solver->addClause(lits);
      }

    } else {

      this->setDirectEncoding();
      this->encode(solver);
    
      if(_vars[1]) {
  
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	

#ifdef _DEBUGWRAP
	std::cout << "add or constraint" << std::endl;
#endif

	// x -> y or z
	lits.clear();
	lits.push(this->equal(0));
	for(int i=0; i<2; ++i) {
	  lits.push(~(_vars[i]->equal(0)));
	}
	solver->addClause(lits);
      
	// y -> x  and z -> x
	for(int i=0; i<2; ++i) {
	  lits.clear();
	  lits.push(_vars[i]->equal(0));
	  lits.push(~(this->equal(0)));
	  solver->addClause(lits);
	}

      } else if(_rhs == 0) {
	
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "add or constraint" << std::endl;
#endif

	// x -> y 
	lits.clear();
	lits.push(this->equal(0));
	lits.push(~(_vars[0]->equal(0)));
	solver->addClause(lits);
      
	// y -> x 
	lits.clear();
	lits.push(_vars[0]->equal(0));
	lits.push(~(this->equal(0)));
	solver->addClause(lits);
      }
      
    }
  } 

  return this;
}



MiniSat_and::MiniSat_and(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setDirectEncoding();
  _vars[1]->setDirectEncoding();

}

MiniSat_and::MiniSat_and(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs)
{

#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setDirectEncoding();

}


MiniSat_and::~MiniSat_and(){

#ifdef _DEBUGWRAP
  std::cout << "delete or" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_and::add(MiniSatSolver *solver, bool top_level){
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    vec<Lit> lits;
    if(top_level) {
            
      if(_vars[1]) {
  
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	

#ifdef _DEBUGWRAP
	std::cout << "add or constraint" << std::endl;
#endif

	for(int i=0; i<2; ++i) {
	  lits.clear();
	  lits.push(~(_vars[i]->equal(0)));
	  solver->addClause(lits);
	}
      
      } else if(_rhs == 1) {
	
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "add or constraint" << std::endl;
#endif

	lits.clear();
	lits.push(~(_vars[0]->equal(0)));
	solver->addClause(lits);
      }

    } else {

      this->setDirectEncoding();
      this->encode(solver);

    
      if(_vars[1]) {
  
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	

#ifdef _DEBUGWRAP
	std::cout << "add or constraint" << std::endl;
#endif

	// y and z -> x
	// x or ~y or ~z
	lits.clear();
	lits.push(~(this->equal(0)));
	for(int i=0; i<2; ++i) {
	  lits.push(_vars[i]->equal(0));
	}
	solver->addClause(lits);
      
	// x -> y  and x -> z
	for(int i=0; i<2; ++i) {
	  lits.clear();
	  lits.push(this->equal(0));
	  lits.push(~(_vars[i]->equal(0)));
	  solver->addClause(lits);
	}

      } else if(_rhs == 1) {
	
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "add or constraint" << std::endl;
#endif

	// x -> y 
	lits.clear();
	lits.push(this->equal(0));
	lits.push(~(_vars[0]->equal(0)));
	solver->addClause(lits);
      
	// y -> x 
	lits.clear();
	lits.push(_vars[0]->equal(0));
	lits.push(~(this->equal(0)));
	solver->addClause(lits);
      }
      
    }
  } 

  return this;
}


MiniSat_eq::MiniSat_eq(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating equality" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setDirectEncoding();
  _vars[1]->setDirectEncoding();

}

MiniSat_eq::MiniSat_eq(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs)
{

#ifdef _DEBUGWRAP
  std::cout << "creating equality" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setDirectEncoding();

}


MiniSat_eq::~MiniSat_eq(){

#ifdef _DEBUGWRAP
  std::cout << "delete eq" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_eq::add(MiniSatSolver *solver, bool top_level){
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    vec<Lit> lits;
    if(top_level) {
            
      if(_vars[1]) {
  
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);


#ifdef _DEBUGWRAP
      std::cout << "add equality" << std::endl;
#endif

	int lb = std::max(_vars[0]->_lower, _vars[1]->_lower);
	int ub = std::min(_vars[0]->_upper, _vars[1]->_upper);

	Lit p;

	for(int x=0; x<2; ++x)
	  for(int i=_vars[x]->_lower; i<lb; ++i) {
	    lits.clear();
	    p = ~(_vars[x]->equal(i));
	    if(p != Lit_True) {
	      lits.push(p);
	      solver->addClause(lits);
	    }
	  }

	for(int x=0; x<2; ++x)
	  for(int i=ub+1; i<=_vars[x]->_upper; ++i) {
	    lits.clear();
	    p = ~(_vars[x]->equal(i));
	    if(p != Lit_True) {
	      lits.push(p);
	      solver->addClause(lits);
	    }
	  }

	bool valid;
	for(int i=lb; i<=ub; ++i) {
	  valid = false;
	  lits.clear();
	  p = ~(_vars[0]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }
	  
	  p = _vars[1]->equal(i);
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }

	  if(!valid) solver->addClause(lits);

	  valid = false;
	  lits.clear();
	  p = ~(_vars[1]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }
	  
	  p = _vars[0]->equal(i);
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }

	  if(!valid) solver->addClause(lits);
	}

      } else {
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "add equality" << std::endl;
#endif

	lits.clear();
	lits.push(_vars[0]->equal(_rhs));
	solver->addClause(lits);
      }
    } else {
      
      //std::cout << "add equal predicate (not supported)" << std::endl;
      //exit(1);

      if(_vars[1]) {
  
	this->setDirectEncoding();
	this->encode(solver);

	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	
	
#ifdef _DEBUGWRAP
	std::cout << "add equality predicate" << std::endl;
#endif
	
	int lb = std::max(_vars[0]->_lower, _vars[1]->_lower);
	int ub = std::min(_vars[0]->_upper, _vars[1]->_upper);

	Lit p;

	for(int x=0; x<2; ++x)
	  for(int i=_vars[x]->_lower; i<lb; ++i) {
	    // k not in z and yk -> -x
	    lits.clear();
	    p = ~(_vars[x]->equal(i));
	    if(p != Lit_True) {
	      lits.push(this->equal(0));
	      lits.push(p);
	      solver->addClause(lits);
	    }
	  }

	for(int x=0; x<2; ++x)
	  for(int i=ub+1; i<=_vars[x]->_upper; ++i) {
	    lits.clear();
	    p = ~(_vars[x]->equal(i));
	    if(p != Lit_True) {
	      lits.push(this->equal(0));
	      lits.push(p);
	      solver->addClause(lits);
	    }
	  }

	bool valid;
	for(int i=lb; i<=ub; ++i) {
	  // yi and zi -> x
	  // -yi or -zi or x
	  valid = false;
	  lits.clear();
	  lits.push(this->equal(1));

	  p = ~(_vars[0]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }
	  
	  p = ~(_vars[1]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }

	  if(!valid) solver->addClause(lits);

	  // -yi and zi -> -x
	  // yi or -zi or -x
	  valid = false;
	  lits.clear();
	  lits.push(this->equal(0));

	  p = _vars[0]->equal(i);
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }
	  
	  p = ~(_vars[1]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }

	  if(!valid) solver->addClause(lits);

	  // yi and -zi -> -x
	  // -yi or zi or -x
	  valid = false;
	  lits.clear();
	  lits.push(this->equal(0));

	  p = ~(_vars[0]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }
	  
	  p = _vars[1]->equal(i);
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }

	  if(!valid) solver->addClause(lits);
	}

      } else {
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "add equality" << std::endl;
#endif

// 	lits.clear();
// 	lits.push(this->equal(0));
// 	lits.push(_vars[0]->equal(_rhs));
// 	solver->addClause(lits);

// 	lits.clear();
// 	lits.push(this->equal(1));
// 	lits.push(~(_vars[0]->equal(_rhs)));
// 	solver->addClause(lits);
      }
      
    }
  } 

  return this;
}


Lit MiniSat_eq::less_or_equal(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::less_or_equal(value);
  else {
    if(value < 0) return Lit_False;
    if(value >= 1) return Lit_True;
    // the value is 0, eq is false iff vars[0] != _rhs
    return ~(_vars[0]->equal(_rhs));
  }
}

Lit MiniSat_eq::greater_than(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::greater_than(value);
  else {
    if(value < 0) return Lit_True;
    if(value >= 1) return Lit_False;
    // the value is 0, eq is true iff vars[0] == _rhs
    return _vars[0]->equal(_rhs);
  }
}

Lit MiniSat_eq::equal(const int value) const {
  if(!_vars[1]) {
    if(value < 0 || value > 1) return Lit_False;

    if(value == 0) {
      // the value is 0, eq is false iff vars[0] != _rhs
      return ~(_vars[0]->equal(_rhs));
    }
    if(value == 1) {
      // the value is 1, eq is true iff vars[0] == _rhs
      return _vars[0]->equal(_rhs);
    }
  }return this->MiniSat_Expression::equal(value);
}

/* Disequality operator */

MiniSat_ne::MiniSat_ne(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2)
{
  //_vars[0]->setDirectEncoding();
  //_vars[1]->setDirectEncoding();
  
  _vars[0]->setDirectEncoding();
  _vars[1]->setDirectEncoding();
  _lower = 0;
  _upper = 1;
  _size = 2;

}

MiniSat_ne::MiniSat_ne(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs)
{

  _vars[0]->setDirectEncoding();
  _lower = 0;
  _upper = 1;
  _size = 2;

}

MiniSat_ne::~MiniSat_ne(){

#ifdef _DEBUGWRAP
  std::cout << "delete notequal" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_ne::add(MiniSatSolver *solver, bool top_level){
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    vec<Lit> lits;
    if(top_level){

      if(_vars[1]) {
	
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	
#ifdef _DEBUGWRAP
	std::cout << "add notequal constraint" << std::endl;
#endif
	
	int lb = std::max(_vars[0]->_lower, _vars[1]->_lower);
	int ub = std::min(_vars[0]->_upper, _vars[1]->_upper);

	int val;
	
	for(val=lb; val<=ub; ++val) {
	  lits.clear();
	  lits.push(~(_vars[0]->equal(val)));
	lits.push(~(_vars[1]->equal(val)));
	solver->addClause(lits);
	}

      } else  {

	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
	std::cout << "add notequal" << std::endl;
#endif

	lits.clear();
	lits.push(~(_vars[0]->equal(_rhs)));
	solver->addClause(lits);

      }
      
    } else {

      if(_vars[1]) {
	
	this->setDirectEncoding();
	this->encode(solver);

	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	
#ifdef _DEBUGWRAP
	std::cout << "add notequal predicate" << std::endl;
#endif
	
	int lb = std::max(_vars[0]->_lower, _vars[1]->_lower);
	int ub = std::min(_vars[0]->_upper, _vars[1]->_upper);

	Lit p;


	for(int x=0; x<2; ++x)
	  for(int i=_vars[x]->_lower; i<lb; ++i) {
	    // k not in z and yk -> x
	    lits.clear();
	    p = ~(_vars[x]->equal(i));
	    if(p != Lit_True) {
	      lits.push(this->equal(1));
	      lits.push(p);
	      solver->addClause(lits);
	    }
	  }

	for(int x=0; x<2; ++x)
	  for(int i=ub+1; i<=_vars[x]->_upper; ++i) {
	    lits.clear();
	    p = ~(_vars[x]->equal(i));
	    if(p != Lit_True) {
	      lits.push(this->equal(1));
	      lits.push(p);
	      solver->addClause(lits);
	    }
	  }

	bool valid;
	for(int i=lb; i<=ub; ++i) {
	  // yi and zi -> -x
	  // -yi or -zi or -x
	  valid = false;
	  lits.clear();
	  lits.push(this->equal(0));

	  p = ~(_vars[0]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }
	  
	  p = ~(_vars[1]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }

	  if(!valid) solver->addClause(lits);

	  // -yi and zi -> x
	  // yi or -zi or x
	  valid = false;
	  lits.clear();
	  lits.push(this->equal(1));

	  p = _vars[0]->equal(i);
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }
	  
	  p = ~(_vars[1]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }

	  if(!valid) solver->addClause(lits);

	  // yi and -zi -> x
	  // -yi or zi or x
	  valid = false;
	  lits.clear();
	  lits.push(this->equal(1));

	  p = ~(_vars[0]->equal(i));
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }
	  
	  p = _vars[1]->equal(i);
	  if(p != Lit_False) {
	    valid |= (p == Lit_True);
	    lits.push(p);
	  }

	  if(!valid) solver->addClause(lits);
	}


      } else  {
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
	std::cout << "add notequal" << std::endl;
#endif

// 	lits.clear();
// 	lits.push(this->equal(1));
// 	lits.push(_vars[0]->equal(_rhs));
// 	solver->addClause(lits);

// 	lits.clear();
// 	lits.push(this->equal(0));
// 	lits.push(~(_vars[0]->equal(_rhs)));
// 	solver->addClause(lits);
      }

    }
  }

  return this;
}


Lit MiniSat_ne::less_or_equal(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::less_or_equal(value);
  else {
    if(value < 0) return Lit_False;
    if(value >= 1) return Lit_True;
    // the value is 0, ne is false iff vars[0] == _rhs
    return _vars[0]->equal(_rhs);
  }
}

Lit MiniSat_ne::greater_than(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::greater_than(value);
  else {
    if(value < 0) return Lit_True;
    if(value >= 1) return Lit_False;
    // the value is 0, ne is true iff vars[0] != _rhs
    return ~(_vars[0]->equal(_rhs));
  }
}

Lit MiniSat_ne::equal(const int value) const {
  if(!_vars[1]) {
    if(value < 0 || value > 1) return Lit_False;

    if(value == 0) {
      // the value is 0, ne is false iff vars[0] == _rhs
      return _vars[0]->equal(_rhs);
    }
    if(value == 1) {
      // the value is 1, ne is true iff vars[0] != _rhs
      return ~(_vars[0]->equal(_rhs));
    }
  }return this->MiniSat_Expression::equal(value);
}


/* Leq operator */

MiniSat_le::MiniSat_le(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2)
{ 
#ifdef _DEBUGWRAP
  std::cout << "Creating Le constraint" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setOrderEncoding();
  _vars[1]->setOrderEncoding();
}

MiniSat_le::MiniSat_le(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs)
{

#ifdef _DEBUGWRAP
  std::cout << "Leq on constant constructor" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setOrderEncoding();
}


MiniSat_le::~MiniSat_le(){

#ifdef _DEBUGWRAP
  std::cout << "delete lessequal" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_le::add(MiniSatSolver *solver, bool top_level) {
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this);

    if(top_level) {

      if(_vars[1]) {

	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding leq constraint" << std::endl;
#endif

	vec<Lit> lits;
	int i, j;
	
	// x1_j -> x0_j-1 
	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
	  if(j < _vars[0]->_lower) {
	    //x1_y is inconsistent
	    lits.clear();
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j < _vars[0]->_upper) { 
	    lits.clear();
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    lits.push(_vars[0]->less_or_equal(j));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[1]->_upper; i<_vars[0]->_upper; ++i) {
	  lits.clear();
	  lits.push(_vars[0]->less_or_equal(i));
	  solver->addClause(lits);
	}
      } else {
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding leq constraint" << std::endl;
#endif

	vec<Lit> lits;
	
	// x0 <= r
	if(_rhs < _vars[0]->_lower) {
	  std::cerr << "model is inconsistent, -exiting" <<std::endl;
	  exit(1);
	} else if(_rhs < _vars[0]->_upper) {
	  lits.clear();
	  lits.push(_vars[0]->less_or_equal(_rhs));
	  solver->addClause(lits);
	}
      }
    } else {
  
#ifdef _DEBUGWRAP    
      std::cout << "Adding leq predicate" << std::endl;
#endif

      if(_vars[1]) {

	this->setDirectEncoding();
	this->encode(solver);
	
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);

	vec<Lit> lits;
	int i, j;
	
	// x -> y <= z
	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
	  if(j < _vars[0]->_lower) {
	    //x1_y is inconsistent
	    lits.clear();
	    lits.push((this->equal(0)));
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j < _vars[0]->_upper) { 
	    lits.clear();
	    lits.push((this->equal(0)));
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    lits.push(_vars[0]->less_or_equal(j));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[1]->_upper; i<_vars[0]->_upper; ++i) {
	  lits.clear();
	  lits.push((this->equal(0)));
	  lits.push(_vars[0]->less_or_equal(i));
	  solver->addClause(lits);
	}

	// ~x -> z < y
	for(j=_vars[0]->_lower; j<_vars[0]->_upper; ++j) {
	  if(j-1 < _vars[1]->_lower) {
	    //x1_y is incomsistent
	    lits.clear();
	    lits.push(~(this->equal(0)));
	    lits.push(~(_vars[0]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j-1 < _vars[1]->_upper) { 
	    lits.clear();
	    lits.push(~(this->equal(0)));
	    lits.push(~(_vars[0]->less_or_equal(j)));
	    lits.push(_vars[1]->less_or_equal(j-1));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[0]->_upper-1; i<_vars[1]->_upper; ++i) {
	  lits.clear();
	  lits.push(~(this->equal(0)));
	  lits.push(_vars[1]->less_or_equal(i));
	  solver->addClause(lits);
	}


      } else {
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding leq predicate" << std::endl;
#endif

// 	vec<Lit> lits;
	
// 	// x -> y <= r
// 	if(_rhs < _vars[0]->_lower) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	} else if(_rhs < _vars[0]->_upper) {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	  lits.push(_vars[0]->less_or_equal(_rhs));
// 	} else {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	}
// 	solver->addClause(lits);

// 	// x -> r < y
// 	if(_rhs >= _vars[0]->_upper) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	} else if(_rhs >= _vars[0]->_lower) {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	  lits.push(_vars[0]->greater_than(_rhs));
// 	} else {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	}
// 	solver->addClause(lits);
      }



    } 
  }

  return this;
}


Lit MiniSat_le::less_or_equal(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::less_or_equal(value);
  else {
    if(value < 0) return Lit_False;
    if(value >= 1) return Lit_True;
    // the value is 0, le is false iff vars[0] > _rhs
    return _vars[0]->greater_than(_rhs);
  }
}

Lit MiniSat_le::greater_than(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::greater_than(value);
  else {
    if(value < 0) return Lit_True;
    if(value >= 1) return Lit_False;
    // the value is 0, le is true iff vars[0] <= _rhs
    return _vars[0]->less_or_equal(_rhs);
  }
}

Lit MiniSat_le::equal(const int value) const {
  if(!_vars[1]) {
    if(value < 0 || value > 1) return Lit_False;

    if(value == 0) {
      // the value is 0, le is false iff vars[0] > _rhs
      return _vars[0]->greater_than(_rhs);
    }
    if(value == 1) {
      // the value is 1, le is true iff vars[0] <= _rhs
      return _vars[0]->less_or_equal(_rhs);
    }
  }
  return this->MiniSat_Expression::equal(value);
}


/* Geq operator */

MiniSat_ge::MiniSat_ge(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating Ge constraint" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setOrderEncoding();
  _vars[1]->setOrderEncoding();
}

MiniSat_ge::MiniSat_ge(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs)
{

#ifdef _DEBUGWRAP
  std::cout << "Geq on constant constructor" << std::endl;
#endif

  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setOrderEncoding();

}

MiniSat_ge::~MiniSat_ge(){

#ifdef _DEBUGWRAP
  std::cout << "delete greaterequal" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_ge::add(MiniSatSolver *solver, bool top_level){ 
  if(!has_been_added()) {
    _solver = solver;
    _ident = _solver->declare(this);

    if(top_level) {

      if(_vars[1]) {

	// reverse less than or equal constraint
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);

	MiniSat_Expression *aux = _vars[1];
	_vars[1] = _vars[0];
	_vars[0] = aux;

#ifdef _DEBUGWRAP
      std::cout << "Adding geq constraint" << std::endl;
#endif

	vec<Lit> lits;
	int i, j;
	
	// x1_j -> x0_j-1 
	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
	  if(j < _vars[0]->_lower) {
	    //x1_y is inconsistent
	    lits.clear();
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j < _vars[0]->_upper) { 
	    lits.clear();
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    lits.push(_vars[0]->less_or_equal(j));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[1]->_upper; i<_vars[0]->_upper; ++i) {
	  lits.clear();
	  lits.push(_vars[0]->less_or_equal(i));
	  solver->addClause(lits);
	}
      } else {

	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding geq constraint" << std::endl;
#endif

	vec<Lit> lits;
	
	// x0 >= r
	if(_rhs > _vars[0]->_upper) {
	  std::cerr << "model is inconsistent, -exiting" <<std::endl;
	  exit(1);
	} else if(_rhs > _vars[0]->_lower) {
	  lits.clear();
	  lits.push(_vars[0]->greater_than(_rhs-1));
	  solver->addClause(lits);
	}
      }
    } else {

      //std::cout << "do not support this constraint yet" << std::endl;
      //exit(1);
  
#ifdef _DEBUGWRAP    
      std::cout << "Adding geq predicate" << std::endl;
#endif

      if(_vars[1]) {

	this->setDirectEncoding();
	this->encode(solver);

	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);

	MiniSat_Expression *aux = _vars[1];
	_vars[1] = _vars[0];
	_vars[0] = aux;

	vec<Lit> lits;
	int i, j;
	
	// x -> y <= z
	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
	  if(j < _vars[0]->_lower) {
	    //x1_y is inconsistent
	    lits.clear();
	    lits.push((this->equal(0)));
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j < _vars[0]->_upper) { 
	    lits.clear();
	    lits.push((this->equal(0)));
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    lits.push(_vars[0]->less_or_equal(j));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[1]->_upper; i<_vars[0]->_upper; ++i) {
	  lits.clear();
	  lits.push((this->equal(0)));
	  lits.push(_vars[0]->less_or_equal(i));
	  solver->addClause(lits);
	}

	// ~x -> z < y
	for(j=_vars[0]->_lower; j<_vars[0]->_upper; ++j) {
	  if(j-1 < _vars[1]->_lower) {
	    //x1_y is incomsistent
	    lits.clear();
	    lits.push(~(this->equal(0)));
	    lits.push(~(_vars[0]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j-1 < _vars[1]->_upper) { 
	    lits.clear();
	    lits.push(~(this->equal(0)));
	    lits.push(~(_vars[0]->less_or_equal(j)));
	    lits.push(_vars[1]->less_or_equal(j-1));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[0]->_upper-1; i<_vars[1]->_upper; ++i) {
	  lits.clear();
	  lits.push(~(this->equal(0)));
	  lits.push(_vars[1]->less_or_equal(i));
	  solver->addClause(lits);
	}


      } else {

	_vars[0] = _vars[0]->add(solver, false);
	//_vars[1] = NULL;

#ifdef _DEBUGWRAP
      std::cout << "Adding geq predicate" << std::endl;
#endif

// 	vec<Lit> lits;
	
// 	// x -> y <= r
// 	if(_rhs < _vars[0]->_lower) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	} else if(_rhs < _vars[0]->_upper) {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	  lits.push(_vars[0]->less_or_equal(_rhs));
// 	} else {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	}
// 	solver->addClause(lits);

// 	// x -> r < y
// 	if(_rhs >= _vars[0]->_upper) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	} else if(_rhs >= _vars[0]->_lower) {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	  lits.push(_vars[0]->greater_than(_rhs));
// 	} else {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	}
// 	solver->addClause(lits);
      }

    } 
  }

  return this;
}

Lit MiniSat_ge::less_or_equal(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::less_or_equal(value);
  else {
    if(value < 0) return Lit_False;
    if(value >= 1) return Lit_True;
    // the value is 0, ge is false iff vars[0] < _rhs
    return _vars[0]->less_or_equal(_rhs-1);
  }
}

Lit MiniSat_ge::greater_than(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::greater_than(value);
  else {
    if(value < 0) return Lit_True;
    if(value >= 1) return Lit_False;
    // the value is 0, ge is true iff vars[0] >= _rhs
    return _vars[0]->greater_than(_rhs-1);
  }
}

Lit MiniSat_ge::equal(const int value) const {
  if(!_vars[1]) {
    if(value < 0 || value > 1) return Lit_False;

    if(value == 0) {
      // the value is 0, ge is false iff vars[0] < _rhs
      return _vars[0]->less_or_equal(_rhs-1);
    }
    if(value == 1) {
      // the value is 1, ge is true iff vars[0] >= _rhs
      return _vars[0]->greater_than(_rhs-1);
    }
  }
  return this->MiniSat_Expression::equal(value);
}


/* Lt object */

MiniSat_lt::MiniSat_lt(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2)
{
  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setOrderEncoding();
  _vars[1]->setOrderEncoding();
}

MiniSat_lt::MiniSat_lt(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs)
{
  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setOrderEncoding();
}

MiniSat_lt::~MiniSat_lt(){

#ifdef _DEBUGWRAP
  std::cout << "delete lessthan" << std::endl;
#endif

}

MiniSat_Expression* MiniSat_lt::add(MiniSatSolver *solver, bool top_level) {
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    if(top_level) {

      if(_vars[1]) {
	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding lt constraint" << std::endl;
#endif

	vec<Lit> lits;
	int i, j, x, y;
	
	// x1_j -> x0_j-1 
	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
	  if(j-1 < _vars[0]->_lower) {
	    //x1_y is incomsistent
	    lits.clear();
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j-1 < _vars[0]->_upper) { 
	    lits.clear();
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    lits.push(_vars[0]->less_or_equal(j-1));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[1]->_upper-1; i<_vars[0]->_upper; ++i) {
	  lits.clear();
	  lits.push(_vars[0]->less_or_equal(i));
	  solver->addClause(lits);
	}


// 	//x1 <= j -> x0 <= j-1 
// 	std::cout << " X < Y " << std::endl;
// 	i=0;
// 	j=0;
// 	x = _vars[0]->getval(i);
// 	y = _vars[1]->getval(j);
// 	while(y <= x) {
	  
// 	}


// 	for(i=0; i<_vars[0]->_size; ++i) {



// 	for(i=0; i<_vars[1]->_size; ++i) {
// 	  x = _vars[1]->getval(i);
// 	  y = _vars[0]->getPrev(x);
// 	  std::cout << _vars[1]->_ident << " <= " << x << " -> " << _vars[0]->_ident << " <= " << y << std::endl;

// 	  if(y<x) {
// 	    lits.clear();
// 	    lits.push(_vars[1]->greater_than(x));
// 	    lits.push(_vars[0]->less_or_equal(y));
// 	  } else
// 	  for(j=0; j<_vars[1]->_size; ++j) { 


// 	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
// 	  if(j-1 < _vars[0]->_lower) {
// 	    //x1_y is incomsistent
// 	    lits.clear();
// 	    lits.push(~(_vars[1]->less_or_equal(j)));
// 	    solver->addClause(lits);
// 	  } else if(j-1 < _vars[0]->_upper) { 
// 	    lits.clear();
// 	    lits.push(~(_vars[1]->less_or_equal(j)));
// 	    lits.push(_vars[0]->less_or_equal(j-1));
// 	    solver->addClause(lits);
// 	  } 
// 	}
// 	for(i=_vars[1]->_upper-1; i<_vars[0]->_upper; ++i) {
// 	  lits.clear();
// 	  lits.push(_vars[0]->less_or_equal(i));
// 	  solver->addClause(lits);
// 	}

      } else {
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding lt constraint" << std::endl;
#endif

	vec<Lit> lits;
	
	// x0 < r
	if(_rhs <= _vars[0]->_lower) {
	  std::cerr << "model is inconsistent, -exiting" <<std::endl;
	  exit(1);
	} else if(_rhs <= _vars[0]->_upper) {
	  lits.clear();
	  lits.push(_vars[0]->less_or_equal(_rhs-1));
	  solver->addClause(lits);
	}
      }

    } else {

#ifdef _DEBUGWRAP
	std::cout << "Adding lt predicate" << std::endl;
#endif	

      if(_vars[1]) {

	this->setDirectEncoding();
	this->encode(solver);

	_vars[0] = _vars[0]->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);

	vec<Lit> lits;
	int i, j;

	// x -> y < z
	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
	  if(j-1 < _vars[0]->_lower) {
	    lits.clear();
	    lits.push((this->equal(0)));
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j-1 < _vars[0]->_upper) { 
	    lits.clear();
	    lits.push((this->equal(0)));
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    lits.push(_vars[0]->less_or_equal(j-1));
	    solver->addClause(lits);

	  } 
	}
	for(i=_vars[1]->_upper-1; i<_vars[0]->_upper; ++i) {
	  lits.clear();
	  lits.push((this->equal(0)));
	  lits.push(_vars[0]->less_or_equal(i));
	  solver->addClause(lits);
	}

	// ~x -> y >= z
	for(j=_vars[0]->_lower; j<_vars[0]->_upper; ++j) {
	  if(j < _vars[1]->_lower) {
	    //x1_y is inconsistent
	    lits.clear();
	    lits.push(~(this->equal(0)));
	    lits.push(~(_vars[0]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j < _vars[1]->_upper) { 
	    lits.clear();
	    lits.push(~(this->equal(0)));
	    lits.push(~(_vars[0]->less_or_equal(j)));
	    lits.push(_vars[1]->less_or_equal(j));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[0]->_upper; i<_vars[1]->_upper; ++i) {
	  lits.clear();
	  lits.push(~(this->equal(0)));
	  lits.push(_vars[1]->less_or_equal(i));
	  solver->addClause(lits);
	}

      } else {

	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
	std::cout << "Adding lt constraint" << std::endl;
#endif
	
// 	vec<Lit> lits;
	
// 	// x -> y < r
// 	if(_rhs <= _vars[0]->_lower) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	} else if(_rhs <= _vars[0]->_upper) {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	  lits.push(_vars[0]->less_or_equal(_rhs-1));
// 	} else {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	}
// 	solver->addClause(lits);

// 	// x -> r < y
// 	if(_rhs > _vars[0]->_upper) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	} else if(_rhs > _vars[0]->_lower) {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	  lits.push(_vars[0]->greater_than(_rhs-1));
// 	} else {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	}
// 	solver->addClause(lits);
      }

    } 
  }

  return this;
}


Lit MiniSat_lt::less_or_equal(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::less_or_equal(value);
  else {
    if(value < 0) return Lit_False;
    if(value >= 1) return Lit_True;
    // the value is 0, lt is false iff vars[0] >= _rhs
    return _vars[0]->greater_than(_rhs-1);
  }
}

Lit MiniSat_lt::greater_than(const int value) const {
  if(_vars[1]) return this->MiniSat_Expression::greater_than(value);
  else {
    if(value < 0) return Lit_True;
    if(value >= 1) return Lit_False;
    // the value is 0, lt is true iff vars[0] < _rhs
    return _vars[0]->less_or_equal(_rhs-1);
  }
}

Lit MiniSat_lt::equal(const int value) const {
  if(!_vars[1]) {
    if(value < 0 || value > 1) return Lit_False;

    if(value == 0) {
      // the value is 0, lt is false iff vars[0] >= _rhs
      return _vars[0]->greater_than(_rhs-1);
    }
    if(value == 1) {
      // the value is 1, lt is true iff vars[0] < _rhs
      return _vars[0]->less_or_equal(_rhs-1);
    }
  }return this->MiniSat_Expression::equal(value);
}


/* Gt object */

MiniSat_gt::MiniSat_gt(MiniSat_Expression *var1, MiniSat_Expression *var2)
  : MiniSat_binop(var1,var2)
{
  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setOrderEncoding();
  _vars[1]->setOrderEncoding();
}

MiniSat_gt::MiniSat_gt(MiniSat_Expression *var1, int rhs)
  : MiniSat_binop(var1,rhs)
{
  _lower = 0;
  _upper = 1;
  _size = 2;
  _vars[0]->setOrderEncoding();
}

MiniSat_gt::~MiniSat_gt(){
}

MiniSat_Expression* MiniSat_gt::add(MiniSatSolver *solver, bool top_level) {
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    if(top_level) {

      if(_vars[1]) {

	
	MiniSat_Expression *aux = _vars[1]->add(solver, false);
	_vars[1] = _vars[0]->add(solver, false);
	_vars[0] = aux;

#ifdef _DEBUGWRAP
      std::cout << "Adding gt constraint" << std::endl;
#endif

	vec<Lit> lits;
	int i, j;
	
	// x1_j -> x0_j-1 
	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
	  if(j-1 < _vars[0]->_lower) {
	    //x1_y is incomsistent
	    lits.clear();
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j-1 < _vars[0]->_upper) { 
	    lits.clear();
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    lits.push(_vars[0]->less_or_equal(j-1));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[1]->_upper-1; i<_vars[0]->_upper; ++i) {
	  lits.clear();
	  lits.push(_vars[0]->less_or_equal(i));
	  solver->addClause(lits);
	}

      } else {
	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
      std::cout << "Adding gt constraint" << std::endl;
#endif
      
	vec<Lit> lits;
	
	// x0 > r
	if(_rhs >= _vars[0]->_upper) {
	  std::cerr << "model is inconsistent, -exiting" <<std::endl;
	  exit(1);
	} else if(_rhs >= _vars[0]->_lower) {
	  lits.clear();
	  lits.push(_vars[0]->greater_than(_rhs));
	  solver->addClause(lits);
	}
      }

    } else {

//       this->setDirectEncoding();
//       this->encode(solver);

#ifdef _DEBUGWRAP
	std::cout << "Adding gt predicate" << std::endl;
#endif	

      if(_vars[1]) {
	
	this->setDirectEncoding();
	this->encode(solver);

	MiniSat_Expression *aux = _vars[1]->add(solver, false);
	_vars[1] = _vars[0]->add(solver, false);
	_vars[0] = aux;
	vec<Lit> lits;
	int i, j;

	// x -> y < z
	for(j=_vars[1]->_lower; j<_vars[1]->_upper; ++j) {
	  if(j-1 < _vars[0]->_lower) {
	    lits.clear();
	    lits.push((this->equal(0)));
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j-1 < _vars[0]->_upper) { 
	    lits.clear();
	    lits.push((this->equal(0)));
	    lits.push(~(_vars[1]->less_or_equal(j)));
	    lits.push(_vars[0]->less_or_equal(j-1));
	    solver->addClause(lits);

	  } 
	}
	for(i=_vars[1]->_upper-1; i<_vars[0]->_upper; ++i) {
	  lits.clear();
	  lits.push((this->equal(0)));
	  lits.push(_vars[0]->less_or_equal(i));
	  solver->addClause(lits);
	}

	// ~x -> y >= z
	for(j=_vars[0]->_lower; j<_vars[0]->_upper; ++j) {
	  if(j < _vars[1]->_lower) {
	    //x1_y is inconsistent
	    lits.clear();
	    lits.push(~(this->equal(0)));
	    lits.push(~(_vars[0]->less_or_equal(j)));
	    solver->addClause(lits);
	  } else if(j < _vars[1]->_upper) { 
	    lits.clear();
	    lits.push(~(this->equal(0)));
	    lits.push(~(_vars[0]->less_or_equal(j)));
	    lits.push(_vars[1]->less_or_equal(j));
	    solver->addClause(lits);
	  } 
	}
	for(i=_vars[0]->_upper; i<_vars[1]->_upper; ++i) {
	  lits.clear();
	  lits.push(~(this->equal(0)));
	  lits.push(_vars[1]->less_or_equal(i));
	  solver->addClause(lits);
	}

      } else {

	_vars[0] = _vars[0]->add(solver, false);

#ifdef _DEBUGWRAP
	std::cout << "Adding gt constraint" << std::endl;
#endif
	
// 	vec<Lit> lits;
	
// 	// x -> y > r
// 	if(_rhs >= _vars[0]->_upper) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	} else if(_rhs >= _vars[0]->_lower) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	  lits.push(_vars[0]->greater_than(_rhs));
// 	} else {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	}
// 	solver->addClause(lits);

// 	// y > r -> x
// 	if(_rhs > _vars[0]->_upper) {
// 	  lits.clear();
// 	  lits.push(this->equal(0));
// 	} else if(_rhs > _vars[0]->_lower) {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	  lits.push(~(_vars[0]->greater_than(_rhs)));
// 	} else {
// 	  lits.clear();
// 	  lits.push(~(this->equal(0)));
// 	}
// 	solver->addClause(lits);
      }

    } 
  }

  return this;
}


Lit MiniSat_gt::less_or_equal(const int value) const {
  if(!_vars[1]) {
    if(value < 0) return Lit_True;
    if(value >= 1) return Lit_False;
    // the value is 0, gt is false iff vars[0] <= _rhs
    return _vars[0]->less_or_equal(_rhs);
  }
  return this->MiniSat_Expression::less_or_equal(value);
}

Lit MiniSat_gt::greater_than(const int value) const {
  if(!_vars[1]) {
    if(value < 0) return Lit_False;
    if(value >= 1) return Lit_True;
    // the value is 0, gt is true iff vars[0] > _rhs
    return _vars[0]->greater_than(_rhs);
  }
  return this->MiniSat_Expression::greater_than(value);
}

Lit MiniSat_gt::equal(const int value) const {
  if(!_vars[1]) {
    if(value < 0 || value > 1) return Lit_False;
    
    if(value == 0) {
      // the value is 0, gt is false iff vars[0] <= _rhs
      return _vars[0]->less_or_equal(_rhs);
    }
    if(value == 1) {
      // the value is 1, gt is true iff vars[0] > _rhs
      return _vars[0]->greater_than(_rhs);
    }
  }
  return this->MiniSat_Expression::equal(value);
}

/* Minimise object */

MiniSat_Minimise::MiniSat_Minimise(MiniSat_Expression *var)
  : MiniSat_Expression()
{
  _obj = var;
  _obj->setOrderEncoding();
}

MiniSat_Minimise::~MiniSat_Minimise(){
}

MiniSat_Expression* MiniSat_Minimise::add(MiniSatSolver *solver, bool top_level) {
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    _obj = _obj->add(solver, false);

    if(top_level) {
      
#ifdef _DEBUGWRAP
      std::cout << "Adding minimise objective" << std::endl;
#endif
      
      solver->minimise_obj = _obj;
      
    } 
  }

  return this;
}

/* Maximise object */

MiniSat_Maximise::MiniSat_Maximise(MiniSat_Expression *var)
  : MiniSat_Expression()
{
  _obj = var;
  _obj->setOrderEncoding();
}

MiniSat_Maximise::~MiniSat_Maximise(){
}

MiniSat_Expression* MiniSat_Maximise::add(MiniSatSolver *solver, bool top_level) {
  if(!has_been_added()){
    _solver = solver;
    _ident = _solver->declare(this);

    _obj = _obj->add(solver, false);

    if(top_level) {
      
#ifdef _DEBUGWRAP
      std::cout << "Adding minimise objective" << std::endl;
#endif
      
      solver->maximise_obj = _obj;
      
    } 
  }

  return this;
}

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

MiniSatSolver::MiniSatSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "create a minisat solver" << std::endl;
#endif

  minimise_obj = NULL;
  maximise_obj = NULL;
  cp_model = NULL;

  conflict_clause = NULL;
  conflictC = 0;

  result = l_Undef;
  nbSolutions = 0;

  Var dummy = SimpSolver::newVar();
  Lit_True = Lit(dummy);
  Lit_False = ~Lit(dummy);

  vec<Lit> lits;
  lits.push(Lit_True);
  addClause(lits);

  _lit_to_var.push_back(NULL);
  _lit_to_val.push_back(0);

  first_decision_level = -1;

  STARTTIME = cpuTime();

  //minisolver = new SimpSolver();

  // Create the scip object
}

MiniSatSolver::~MiniSatSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "delete wrapped solver" << std::endl;
#endif

//   for(unsigned int i=0; i<_variables.size(); ++i) {
//     delete _variables[i];
//   }

}

int MiniSatSolver::declare(MiniSat_Expression *exp) {
  int id = _variables.size();
  _variables.push_back(exp);
  return id;
}

void MiniSatSolver::store_solution() {
  ++nbSolutions;
  if(model.size() < nVars()) {
    model.growTo(nVars());
    for (int i = 0; i < nVars(); i++) model[i] = value(i);
  }
  if(!cp_model) cp_model = new int[_variables.size()];
  for(unsigned int i=0; i<_variables.size(); ++i) {
    if(_variables[i])
      cp_model[i] = _variables[i]->get_min();
    else
      cp_model[i] = 0;
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

// bool MiniSatSolver::not_false(Lit x)
// {  
//   int res =  (modelValue(x) != l_False);
//   return res;
// }

Var MiniSatSolver::newVar(MiniSat_Expression* var, int val)
 {
   Var x = SimpSolver::newVar();
   //std::cout << x << " " << (_lit_to_var.size()) << std::endl;
   assert(_lit_to_var.size() == (unsigned int)x);
   _lit_to_var.push_back(var);
   _lit_to_val.push_back(val);

   return x;
 }

// int MiniSatSolver::nVars()
// {
//   return nVars();
// }

// void MiniSatSolver::addClause(vec<Lit>& cl)
//  {
//    printClause(cl);

//    Solver::addClause(cl);
//  }

int MiniSatSolver::solveAndRestart(const int policy, 
				   const unsigned int base, 
				   const double factor,
				   const double decay)
{
//   if(decay>0 && decay<1) 
//     var_decay = decay;
//   clause_decay = decay;
//   restart_first = base;
//   restart_inc = factor;
  return solve();
}

int MiniSatSolver::solve()
{

#ifdef _DEBUGWRAP
  std::cout << "call solve" << std::endl;  
#endif 

#ifdef _DEBUGWRAP
  std::cout << "print to file" << std::endl;  
  //if(verbosity > 1)  
  toDimacs("dimacs.out");
#endif 

#ifdef _DEBUGWRAP
  std::cout << "solve!" << std::endl;  
#endif 

 
  //solve(true, true);
  //printStats(*minisolver);

//   verbosity = 1;

  start_time = getRunTime();
  saved_level = init_level;
  if(init_level < decisionLevel())
    init_level = decisionLevel();

  //solver_ptr = &minisolver;
  solver_ptr = this;
  signal(SIGINT,SIGINT_handler);

  //std::cout << "c solve using MiniSat" << std::endl;


  if(minimise_obj) {
    vec<Lit> lits;
    int objective = minimise_obj->_upper;
    
    result = SimpSolver::solve(true,true);
    
    while(result == l_True && !(limitsExpired())) {
      cancelUntil(init_level);

      store_solution();
      objective = minimise_obj->get_value()-1;
      if(objective < minimise_obj->_lower) break;

      if(verbosity >= 0) {
	std::cout << "c  new objective: " << objective+1 << std::endl;
      }
      
//       for(; objective>minimise_obj->_lower; --objective) {

// 	if(!value(minimise_obj->less_or_equal(objective-1))) 
// 	  break;
//       }

//       //      std::cout << "current objective = " << objective << std::endl;

//       --objective;

      
      lits.clear();
      lits.push(minimise_obj->less_or_equal(objective));
      addClause(lits);

      result = SimpSolver::solve(true,true);



      if(result == l_False) {
	++objective;
      }
    }
    //if(is_sat()) result = true;
  }
  else if(maximise_obj) {
    vec<Lit> lits;
    int objective = maximise_obj->_lower;
    
    result = SimpSolver::solve(true,true);
    
    while(result == l_True && !(limitsExpired())) {
      cancelUntil(init_level);

      store_solution();
      objective = maximise_obj->get_value()+1;
      if(objective > maximise_obj->_upper) break;

      if(verbosity >= 0) {
	std::cout << "c  new objective: " << objective-1 << std::endl;
      }
      
//       for(; objective>minimise_obj->_lower; --objective) {

// 	if(!value(minimise_obj->less_or_equal(objective-1))) 
// 	  break;
//       }

//       //      std::cout << "current objective = " << objective << std::endl;

//       --objective;

      
      lits.clear();
      lits.push(maximise_obj->greater_than(objective-1));
      addClause(lits);

      result = SimpSolver::solve(true,true);
      if(result == l_False) {
	++objective;
      } 
    }
    //if(is_sat()) result = true;
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
	//<< " (x" << (_lit_to_var[i+1]->nbj_ident) << "/" <<  _lit_to_val[i+1] << ")" << std::endl;
    }
  else
    std::cout << "unsatisfiable" << std:: endl;
#endif 


  if(is_sat()) return true;
  return false;
  
  //return result;

  //std::cout << "finished" << std::endl;

  //return 0;
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

  //for(int i=0; i<learnts.size(); ++i)
  // removeClause(*(learnts[i]));
  //learnts.clear();

  model.clear();
  delete [] cp_model;
  cp_model = NULL;

  //std::cout << decisionLevel() << std::endl;
  init_level = saved_level;
  ok = true;
}

int MiniSat_Expression::next(int v)
{
  int nxt = v;

  int i=getind(v)+1;
  if(_direct_encoding) {
    for(; i<_size; ++i)
      if(_solver->truth_value(equal(getval(i))) != l_False)
	break;
  } else if(_order_encoding) {
    for(; i<_size; ++i)
      if(_solver->truth_value(less_or_equal(getval(i))) != l_False)
	break;
  }

  if(i < _size)
    nxt = getval(i); 
  return nxt;
}


bool MiniSatSolver::undo(const int nlevel)
{

  //std::cout << decisionLevel() << " -> " << decisionLevel()-nlevel << " / " << first_decision_level << std::endl;

  //if (decisionLevel() == first_decision_level) return false;

  int okay = true;
  backtrack_level = decisionLevel()-nlevel;
  if(backtrack_level < first_decision_level) okay = false;

  learnt_clause.clear();
  //backtrack_level = lvl-nlevel;
 
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
  //Lit p = forced_decisions.last();
  //forced_decisions.pop();
  //std::cout << "b lvl " << decisionLevel() << ": " << var(last_decision) << std::endl;
  uncheckedEnqueue(~(last_decision));
}

// void MiniSatSolver::learn()
// {
//   if (learnt_clause.size() > 0){
//     if (learnt_clause.size() == 1){
//       uncheckedEnqueue(learnt_clause[0]);
//     }else{
//       Clause* c = Clause_new(learnt_clause, true);
//       learnts.push(c);
//       attachClause(*c);
//       claBumpActivity(*c);
//       uncheckedEnqueue(learnt_clause[0], c);
//     }

//     varDecayActivity();
//     claDecayActivity();

//   }
// }

// void MiniSatSolver::branch_on(const char* op, MiniSat_Expression* x, int v)
// {
//   if(op[1] == 't') 
//     if(op[0] == 'g') ++v;
//     else --v;

//   decisions++;

//   Lit next = lit_Undef;
//   switch(op[0]) {
//   case 'e': next =  (x->equal(v)); break;
//   case 'n': next = ~(x->equal(v)); break;
//   case 'g': next =  (x->greater_than(v-1)); break;
//   case 'l': next =  (x->less_or_equal(v )); break;
//   }

//   newDecisionLevel();
//   uncheckedEnqueue(next);
// }

void MiniSatSolver::save()
{
  decisions++;
  newDecisionLevel();
  //std::cout << "saving: " << decisionLevel() << std::endl;
}

void MiniSatSolver::post(const char* op, MiniSat_Expression* x, int v)
{
  if(op[1] == 't') 
    if(op[0] == 'g') ++v;
    else --v;

  int lvl = decisionLevel();
  if(first_decision_level < 0) 
    first_decision_level = lvl-1;

  //learnt_clause.clear();
  Lit next = lit_Undef;


  switch(op[0]) {
  case 'e': {
//     for(i=0; i<x->_size; ++i) {
//       k = x->getval(i);
//       std::cout << v << " " << k << std::endl;
//       if(k == v) learnt_clause.push(x->equal(v));
//       else learnt_clause.push(~(x->equal(v)));
//       //if(k < v) learnt_clause.push(x->greater_than(v));
//       //else learnt_clause.push(x->less_or_equal(v));
//     }
    next =  (x->equal(v)); 
  } break;
  case 'n': next = ~(x->equal(v)); break;
  case 'g': next =  (x->greater_than(v-1)); break;
  case 'l': next =  (x->less_or_equal(v )); break;
  }

  forced_decisions.push(next);
  //std::cout << "a lvl " << lvl << ": " << var(next) << std::endl;
  

//   for(i=0; i<learnt_clause.size(); ++i)
//     uncheckedEnqueue(learnt_clause[i]);
//   learnt_clause.clear();

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



