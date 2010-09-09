
#include <iostream>
#include "MiniSat.hpp"
#include <math.h>
#include <Sort.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


//Lit Lit_False;
//Lit Lit_True;


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

// class StreamBuffer {
//     gzFile  in;
//     char    buf[CHUNK_LIMIT];
//     int     pos;
//     int     size;

//     void assureLookahead() {
//         if (pos >= size) {
//             pos  = 0;
//             size = gzread(in, buf, sizeof(buf)); } }

// public:
//     StreamBuffer(gzFile i) : in(i), pos(0), size(0) {
//         assureLookahead(); }

//     int  operator *  () { return (pos >= size) ? EOF : buf[pos]; }
//     void operator ++ () { pos++; assureLookahead(); }
// };


// //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// template<class B>
// static void skipWhitespace(B& in) {
//     while ((*in >= 9 && *in <= 13) || *in == 32)
//         ++in; }

// template<class B>
// static void skipLine(B& in) {
//     for (;;){
//         if (*in == EOF || *in == '\0') return;
//         if (*in == '\n') { ++in; return; }
//         ++in; } }

// template<class B>
// static int parseInt(B& in) {
//     int     val = 0;
//     bool    neg = false;
//     skipWhitespace(in);
//     if      (*in == '-') neg = true, ++in;
//     else if (*in == '+') ++in;
//     if (*in < '0' || *in > '9') reportf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
//     while (*in >= '0' && *in <= '9')
//         val = val*10 + (*in - '0'),
//         ++in;
//     return neg ? -val : val; }

// template<class B>
// static void readClause(B& in, SimpSolver& S, vec<Lit>& lits) {
//     int     parsed_lit, var;
//     lits.clear();
//     for (;;){
//         parsed_lit = parseInt(in);
//         if (parsed_lit == 0) break;
//         var = abs(parsed_lit)-1;
//         while (var >= S.nVars()) S.newVar();
//         lits.push( (parsed_lit > 0) ? Lit(var) : ~Lit(var) );
//     }
// }

// template<class B>
// static bool match(B& in, const char* str) {
//     for (; *str != 0; ++str, ++in)
//         if (*str != *in)
//             return false;
//     return true;
// }


// template<class B>
// static void parse_DIMACS_main(B& in, SimpSolver& S) {
//     vec<Lit> lits;
//     for (;;){
//         skipWhitespace(in);
//         if (*in == EOF) break;
//         else if (*in == 'p'){
//             if (match(in, "p cnf")){
//                 int vars    = parseInt(in);
//                 int clauses = parseInt(in);
//                 reportf("|  Number of variables:  %-12d                                         |\n", vars);
//                 reportf("|  Number of clauses:    %-12d                                         |\n", clauses);

//                 // SATRACE'06 hack
//                 if (clauses > 4000000)
//                     S.eliminate(true);
//             }else{
//                 reportf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
//             }
//         } else if (*in == 'c' || *in == 'p')
//             skipLine(in);
//         else{
//             readClause(in, S, lits);
//             S.addClause(lits); }
//     }
// }

// // Inserts problem into solver.
// //
// static void parse_DIMACS(gzFile input_stream, SimpSolver& S) {
//     StreamBuffer in(input_stream);
//     parse_DIMACS_main(in, S); }


//=================================================================================================

SimpSolver* solver_ptr;
static void SIGINT_handler(int signum) {
    reportf("\n"); reportf("*** INTERRUPTED ***\n");
    printStats(*solver_ptr);
    reportf("\n"); reportf("*** INTERRUPTED ***\n");
    exit(1); }


/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

MiniSatSolver::MiniSatSolver() : SatWrapperSolver(), SimpSolver()
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

  Solver::newVar();
  //vec<Lit> lits;
  //lits.push(Lit_True);
  //Solver::addClause(lits);
}

MiniSatSolver::~MiniSatSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "delete wrapped solver" << std::endl;
#endif

}

int MiniSatSolver::create_atom(DomainEncoding* dom, const int type) {
  unsigned int id = SatWrapperSolver::create_atom(dom, type);
  Solver::newVar();
  return id; 
}

void MiniSatSolver::validate() {
  vec<Lit> cl;
  unsigned int i;
  while(current < clause_base.size()) {
    cl.clear(false);
    for(i=0; i<clause_base[current].size(); ++i)
      cl.push(clause_base[current][i]);
    //displayClause(clause_base[current]);
    Solver::addClause(cl);
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

bool ret = is_sat();

 return ret;
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
  nbSolutions = 0;

  learnt_clause.clear();
  backtrack_level = init_level;
  //forced_decisions.clear();
  
//   //if(forced_decisions.size()) {
//   for(int i=decisionLevel()-1; i>backtrack_level; --i) {

//     std::cout << "UNDOpop "  << (forced_decisions.size()) << std::endl;


//   //forced_decisions.pop();
//   forced_decisions.pop_back();
//   }


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

  for(int i=decisionLevel()-1; i>backtrack_level; --i) {
    forced_decisions.pop_back();
  }

  //last_decision = forced_decisions.last();
  last_decision = forced_decisions.back();

  forced_decisions.pop_back();

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

  for(int i=decisionLevel()-1; i>backtrack_level; --i) {
    forced_decisions.pop_back();
  }

  //last_decision = forced_decisions.last();
  last_decision = forced_decisions.back();
  //forced_decisions.pop();
  forced_decisions.pop_back();

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
  //forced_decisions.push(lit_Undef);
  forced_decisions.push_back(lit_Undef);

}

void MiniSatSolver::post(const char* op, SatWrapper_Expression* x, int v)
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

  //forced_decisions.push(next);
  //forced_decisions.last() = next;
  forced_decisions.back() = next;

  uncheckedEnqueue(next);
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
  verbosity = degree-1;
}

void MiniSatSolver::setRandomized(const int degree)
{
  if(degree) SimpSolver::setRandomized();
}

void MiniSatSolver::setRandomSeed(const int seed)
{
  setRandomSeed((double)seed);
}

bool MiniSatSolver::is_sat()
{
  return (cp_model != NULL);
}

bool MiniSatSolver::is_unsat()
{
  return (result == l_False && nbSolutions == 0);
}

bool MiniSatSolver::is_opt()
{
  return (result == l_False && nbSolutions);
}

void MiniSatSolver::printStatistics()
{
  printStats(*this);
}

int MiniSatSolver::getBacktracks()
{
  return conflicts;
}

int MiniSatSolver::getNodes()
{
  return decisions;
}

int MiniSatSolver::getFailures()
{
  return conflicts;
}

int MiniSatSolver::getPropags()
{
   return propagations;
}

double MiniSatSolver::getTime()
{
  return cpuTime() - STARTTIME;
}




