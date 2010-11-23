
/*
  Mistral is a constraint satisfaction and optimisation library
  Copyright (C) 2003-2005  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/

/**********************************************
 * DVOs and Random number generator
 **********************************************/

#include <mistral_sol.h>
#include <mistral_csp.h>
#include <mistral_dvo.h>

#include <math.h>

using namespace Mistral;



BranchingStrategy::BranchingStrategy() {}
BranchingStrategy::~BranchingStrategy() {}
    
Contention::Contention(VarArray* scp, int w) {
  width = w;
  length = scp[0].size();
  VariableInt *x;

  int i, j, k;

  num_vars = width*length;
  min_dx = new int[num_vars];
  tmp_proba_ = new double*[num_vars];
  proba_ = new double*[num_vars];
  scope_ = new VariableInt*[num_vars];
  neighborhood = new Vector<int>[num_vars];
  
  for(i=0; i<width; ++i) {
    for(j=0; j<length; ++j) {
      x = scp[i][j].var_ptr_->getVariable();
      if(x && !x->isGround()) {
	active.push(i*length+j);
	neighborhood[i*length+j].init(0,width+length-2);
	k = x->max() - x->min() + 1;
	tmp_proba_[i*length+j] = new double[k];
	proba_[i*length+j] = new double[k];
	std::fill(proba_[i*length+j], 
		  proba_[i*length+j]+k, 
		  1.0/(double)(x->domsize()));
	min_dx[i*length+j] = x->min();
	tmp_proba_[i*length+j] -= min_dx[i*length+j];
	proba_[i*length+j] -= min_dx[i*length+j];
	scope_[i*length+j] = x;
      } else {
	scope_[i*length+j] = NULL;
      }
    }
  }

//   DomainIterator *valit;
   for(i=0; i<width; ++i) 
     for(j=0; j<length; ++j) if(scope_[i*length+j]) {
	for(k=0; k<length; ++k) 
	  if(k!=j && scope_[i*length+k] 
	     && scope_[i*length+k]->intersect(scope_[i*length+j]))
	    neighborhood[i*length+j].push(i*length+k);
	
	for(k=0; k<width; ++k) 
	  if(k!=i && scope_[k*length+j]
	     && scope_[k*length+j]->intersect(scope_[i*length+j]))
	    neighborhood[i*length+j].push(k*length+j);
       }


}
Contention::~Contention() {
  for(int i=0; i<num_vars; ++i) 
    if(scope_[i]) {
      proba_[i] += min_dx[i];
      tmp_proba_[i] += min_dx[i];
      delete [] proba_[i];
      delete [] tmp_proba_[i];
    }
  delete [] neighborhood;
  delete [] scope_;
}

void Contention::extract() {
  for(int i=0; i<num_vars; ++i) {
    if(scope_[i]) {
      delete scope_[i]->branch;
      scope_[i]->branch = new ValSelectorContention(scope_[i], this, i);
    }
  }
}

void Contention::compute_proba() {
  int i, j, k, l, v;
  double p, total;
  DomainIterator *valit;

  for(i=0; i<active.size; ++i) {
    k = active[i];
    total = 0.0;
    valit = scope_[k]->begin();
    do {
      v = *valit;
      p = 1.0;
      for(j=0; j<neighborhood[k].size; ++j) {
	l = neighborhood[k][j];
	if(scope_[l]->contain(v)) {
	  p *= (1.0-proba_[l][v]);
	} 
      }
      tmp_proba_[k][v] = p;
      total += tmp_proba_[k][v];
    } while( valit->next() );

    valit = scope_[k]->begin();
    do {
      v = *valit;
      tmp_proba_[k][v] /= total;
    } while( valit->next() );

  }

  double **tmp = proba_;
  proba_ = tmp_proba_;
  tmp_proba_ = tmp;

}

void Contention::reset() {
  int i, k;
  DomainIterator *valit;
  for(i=0; i<active.size; ++i) {
    k = active[i];
    valit = scope_[k]->begin();
    do proba_[k][*valit] = 1.0/(scope_[k]->domsize());
    while(valit->next());
  }
}

int Contention::get_best(const int id) {
  int v, val = scope_[id]->min();
  DomainIterator *valit;

//   std::cout << "branch on ";
//   scope_[id]->print(std::cout);
//   std::cout << std::endl;

  valit = scope_[id]->begin();
  do {
    v = *valit;
    if(proba_[id][val] < proba_[id][v]) 
      val = v;
  } while( valit->next() );

  //std::cout << "   ==> " << val << std::endl << std::endl;

  return val;
}

int Contention::update(const int id) {
  int j, k, v, val = scope_[id]->min();
  double p, best = 0.0;
  DomainIterator *valit;

//   std::cout << "branch on ";
//   scope_[id]->print(std::cout);
//   std::cout << std::endl;

  valit = scope_[id]->begin();
  do {
    v = *valit;
    p = 1.0;
    for(j=0; j<neighborhood[id].size; ++j) {
      k = neighborhood[id][j];
      if(scope_[k]->contain(v)) 
	p *= (1.0-(1.0/((double)(scope_[k]->domsize()))));
    }


    //std::cout << "\t" << v << " (" << p << ")" << std::endl; 

    if(best < p) {
      best = p;
      val = v;
    }
  } while( valit->next() );

  //std::cout << "   ==> " << val << std::endl << std::endl;

  return val;
}

std::ostream& Contention::print(std::ostream& os, const int id) const {
  DomainIterator *valit;
  valit = scope_[id]->begin();
  do {
    os 
      << " " << std::setw(10)
      << std::setprecision(3)
      << (proba_[id][*valit] >= 0.001 ? proba_[id][*valit] : 0);
  } while( valit->next() );
  os << std::endl;
  return os;
}

std::ostream& Contention::print(std::ostream& os) const {
  int i, j;
//   for(i=0; i<width; ++i) {
//     for(j=0; j<length; ++j) {
//       if(scope[i][j] && !scope[i][j]->isGround()) {

// 	os << setw(2) << i << "." 
// 	   << setw(2) << j << " ";
// 	scope[i][j]->print(os);

// 	valit = scope[i][j]->begin();
// 	do {
// 	  os << " " << (*valit) << " = " 
// 	     << setprecision(3) << proba[i][j][*valit];
// 	} while( valit->next() );
// 	os << std::endl;
//       }
//     }
//   }

  for(i=0; i<width; ++i) {
    for(j=0; j<length; ++j) {
      if(scope_[i*length+j]) {
	print(os, i*length+j);
      }
    }
  }
  
  return os; 
}

PredicateDisjunctive** _garbage_disjuncts;
void free_disjuncts() {
//   for(int i=0; i<_garbage_disjuncts.size; ++i)
//     delete [] _garbage_disjuncts[i];
}

Weighter::Weighter( Solver* s ) 
  : level(s->level)
{
  init_level = s->init_level;
}

WeighterDegree::WeighterDegree( Solver* s ) 
  : Weighter(s)
{

 #ifdef _WEIGHT_STATS

  int i;
  
  solver = s;
  numvars = s->length;
  choice = new int*[numvars+1];

  for(i=0; i<numvars; ++i) {
    choice[i] = new int[numvars];
    std::fill(choice[i], choice[i]+numvars, 0);
  }

  choice[numvars] = new int[numvars];
  for(i=0; i<numvars; ++i) 
    choice[numvars][i] = (numvars-i);


//   std::cout << "HERE" << std::endl;

//   int i, j;
//   VariableInt *aux;

//   numvars = s->length;
//   vars = new VariableInt*[numvars];
//   //vars = s->sequence;
//   order = new int[numvars];
//   mean_weight = 0.0;
//   for(i=0; i<numvars; ++i) {
//     order[i] = i;
//     vars[i] = s->sequence[i];
//     mean_weight += s->sequence[i]->weight; 
//     for(j=i; j && vars[j]->weight < vars[j-1]->weight; --j) {
//       aux = vars[j];
//       vars[j] = vars[j-1];
//       vars[j-1] = aux;

//       order[vars[j]->id] = j;
//       order[vars[j-1]->id] = j-1;
//     }
//   }
//   mean_weight /= numvars;

//   for(int i=1; i<numvars; ++i) {
//     assert( vars[i-1]->weight <= vars[i]->weight );
//   }

#endif

  threshold = (s->length)/2;

}


#ifdef _WEIGHT_STATS

void WeighterDegree::notifyChoice() 
{
  int idx = solver->decision.back()->id;
  ++choice[level-1][idx];
  ++numchoices;
}  


int **num_choices; 
int compar_qsort(const void *a, const void *b) {
  int x = *(int*)a;
  int y = *(int*)b;
  int i=0;

  while(num_choices[i][x] == num_choices[i][y]) ++i;
  if(num_choices[i][x] < num_choices[i][y]) return 1;
  else return -1;
}

int *count90;
int compar_qsort2(const void *a, const void *b) {
  int x = *(int*)a;
  int y = *(int*)b;
  if(count90[x] < count90[y]) return 1;
  else if(count90[x] > count90[y]) return -1;
  else return 0;
}

void WeighterDegree::init_choices() {
  numchoices = 0;
  for(int i=0; i<numvars; ++i) 
    std::fill(choice[i], choice[i]+numvars, 0);
}

void WeighterDegree::print_choices(std::ostream& o) {
  int i, j, current, threshold, vars[numvars], total[numvars], 
    count[numvars], cumul[numvars], vars90[numvars], ct90, cl90;
  BitSet is_in(0, numvars-1, BitSet::empt);
  BitSet is_in_90(0, numvars-1, BitSet::empt);
  for(i=0; i<numvars; ++i) {
    vars[i] = i;
    vars90[i] = i;
    total[i] = 0;
    count[i] = 0;
    cumul[i] = 0;
  }

  num_choices = choice;
  std::qsort(vars, numvars, sizeof(int), compar_qsort);

  for(i=0; i<numvars; ++i) {
    for(j=0; j<numvars; ++j) if(choice[i][j]) {
	++count[i];
	total[i] += choice[i][j];
	is_in.insert(j);
      }
    cumul[i] = is_in.size();
    if(total[i]) {
      o << i << " ";
      o << count[i] << " ";
      o << cumul[i] << " ";
      o << (((double)(count[i]))/((double)numvars)) << " ";
      o << (((double)(cumul[i]))/((double)numvars)) << " ";
      o << total[i] << " ";
      o << (((double)(total[i]))/((double)numchoices)) << " ";
      o << numvars << " " << numchoices << " ";

      count90 = choice[i];
      std::qsort(vars90, numvars, sizeof(int), compar_qsort2);
      threshold = (int)(0.9*((double)total[i]));
      current = 0;
      ct90 = 0;
      for(j=0; j<numvars && choice[i][vars90[j]]>0 && current<=threshold; ++j) {
	current += choice[i][vars90[j]];
	++ct90;
	is_in_90.insert(vars90[j]);
      }
      cl90 = is_in_90.size();
      o << ct90 << " ";
      o << cl90 << " ";
      o << (((double)(ct90))/((double)numvars)) << " ";
      o << (((double)(cl90))/((double)numvars)) << " ";
//       //for(j=0; j<numvars && count[vars90[j]]>0; ++j)
//       //std::cout << " " << (count[vars90[j]]) ;
//       //std::cout << std::endl;
//       for(j=0; j<numvars; ++j) {
// 	if(total[i] && choice[i][vars[j]]) 
// 	  o << (((double)(choice[i][vars[j]])) / ((double)(total[i]))) << " ";
// 	else o << "0 ";
//       }
      o << std::endl;
    }
  }

}

#endif 


void WeighterDegree::notifyFailure( Constraint *con )
{

//   //++con->weight;
//   int old_weight = con->weight;
  
  
//   //int new_weight = (int)(1000 + ceil((double)old_weight * 0.9));

//   int new_weight = (int)ceil(1.01 * (double)old_weight);

//   con->weight = new_weight;
 
//   //std::cout << old_weight << " -> " << new_weight << std::endl;

//   new_weight -= old_weight;


//   for(int i=0; i<con->arity; ++i)
//     con->scope[i]->weight += new_weight;
//   //std::cout << con->weight << std::endl;

// #ifdef _WEIGHT_STATS

//   std::cout << "THERE" << std::endl;

//   int id, j;
//   VariableInt *aux;
//   mean_weight = (mean_weight * numvars);// + con->arity) / numvars);

// #endif


  if( con->arity <= threshold ) {
    ++con->weight;
    for(int i=0; i<con->arity; ++i) {
      ++con->_scope[i]->weight;

// #ifdef _WEIGHT_STATS
      
//       id = con->_scope[i]->id;

//       if(id < numvars) {

// 	assert( vars[order[id]] == con->_scope[i] );

// 	mean_weight += 1.0;

// 	for(j=order[id]; (j < numvars-1) && (vars[j]->weight > vars[j+1]->weight); ++j) {

// // 	  for(int x=id; i<numvars; ++i) {
// // 	    std::cout << " " << vars[x]->weight;
// // 	  }
// // 	  std::cout << std::endl;

// 	  aux = vars[j];
// 	  vars[j] = vars[j+1];
// 	  vars[j+1] = aux;

// 	  order[vars[j]->id] = j;
// 	  order[vars[j+1]->id] = j+1;
// 	}
//       }

// #endif

    }
  }


// #ifdef _WEIGHT_STATS

//   mean_weight = (mean_weight / numvars);

//   for(int i=0; i<numvars; ++i) {
//     std::cout << " " << vars[i]->weight;
//   }
//   std::cout << std::endl;

//   double real_mean = vars[0]->weight;
//   for(int i=1; i<numvars; ++i) {

//     if(vars[i-1]->weight > vars[i]->weight)
//       std::cout << vars[i-1]->weight << " " << vars[i]->weight << std::endl;

//     assert( vars[i-1]->weight <= vars[i]->weight );
//     real_mean += vars[i]->weight;
//   }
  
//   real_mean /= numvars;
//   assert( (mean_weight - real_mean) < 0.001 );
//   assert( (real_mean - mean_weight) < 0.001 );

// #endif

//   Solver *sol = con->scope[0]->solver;

//   sol->printWeightProfile(std::cout, INT_MAX, 1);
//   std::cout << std::endl;

//   //std::cout << con->weight << std::endl;
}

WeighterLevelDegree::WeighterLevelDegree( Solver* s )
  : Weighter(s)
{
  lmax = (int)(log2((double)(s->length)));
}

void WeighterLevelDegree::notifyFailure( Constraint *con )
{
  if( lmax < level ) 
    lmax = level;
  con->weight += (1 + lmax - level);
  for(int i=0; i<con->arity; ++i)
    con->_scope[i]->weight += (1 + lmax - level);
}


// WeighterPruning::WeighterPruning( Solver* s ) 
//   : Weighter(s)
// {
// }

// void WeighterPruning::notifyFailure( Constraint *con )
// {
//   ++con->weight;
//   for(int i=0; i<con->arity; ++i)
//     ++con->scope[i]->weight;
// }


void WeighterImpact::print( std::ostream& o ) const
{
  for(unsigned int i=0; i<ilength; ++i) 
    {
      std::cout << "c ";
      variables[i]->printshort( std::cout );
      std::cout << ":\t";
      int v=variables[i]->min();
      do {
	std::cout << " " << std::setprecision(2) << decision_impact[i][v] << ":" << v; 
      } while( variables[i]->setNext(v) );
      std::cout << std::endl;
    }
}

WeighterImpact::~WeighterImpact( )
{
  for(unsigned int i=0; i<ilength; ++i) 
    {
      decision_impact[i] += varmin[i];
      decision_count [i] += varmin[i];

      delete [] decision_impact[i];
      delete [] decision_count [i];
    }
  delete [] decision_impact;
  delete [] decision_count;
  delete [] X;
  delete [] domain_size;
  delete [] varmin;
  delete [] isRange;
}

WeighterImpact::WeighterImpact( Solver* s ) 
  : Weighter(s), 
    decision(s->decision.stack_),
    first(s->future), 
    last(s->empty),
    //verylast(s->sequence+s->numvars), 
    variables(s->variables.stack_)
{

  //verylast = s->sequence+s->numvars; 

  unsigned int i, m;

  needUpdate = false;
  ilength = s->length;
  varmin = new int[ilength];
  isRange = new bool[ilength];
  decision_impact = new double*[ilength];
  decision_count  = new unsigned int*[ilength];
  X = new VariableInt*[ilength];
  domain_size = new unsigned int[ilength];
    
  for(i=0; i<ilength; ++i) 
    {
      isRange[i] = ( s->variables[i]->getType() == VariableInt::RANGE );
      
      if( isRange[i] ) 
	{
	  decision_impact[i] = new double[2];
	  decision_count [i] = new unsigned int[2];
	  
	  std::fill( decision_count [i], decision_count [i]+2, 1 );
	  std::fill( decision_impact[i], decision_impact[i]+2, 1 );
	  
	  varmin[i] = 0;


	  //std::cout << "Range: " << decision_impact[i] << std::endl;
	}
      else 
	{
	  m = (s->variables[i]->maxCapacity() - s->variables[i]->minCapacity());

	  decision_impact[i] = new double[m];
	  decision_count [i] = new unsigned int[m];
	  
	  std::fill( decision_count [i], decision_count [i]+m, 0 );
	  std::fill( decision_impact[i], decision_impact[i]+m, 1 );//0.0001 );

	  //std::cout << "Domain: " << decision_impact[i] << std::endl;

	  // dom/wdeg
	  //             1/0.361      490 NDS     9841 BTS/s        199956 CKS 0.04837 s


	  // impact/wdeg
	  // 1           1/0.361      670 NDS     9673 BTS/s        238509 CKS 0.06777 s
	  // 0.1         1/0.361      680 NDS     9698 BTS/s        242066 CKS 0.06856 s
	  // 0.01        1/0.361      645 NDS     9597 BTS/s        232364 CKS 0.06566 s
	  // 0.001       1/0.361      615 NDS     9481 BTS/s        224070 CKS 0.06334 s
	  // 0.0001      1/0.361      580 NDS     9362 BTS/s        215040 CKS 0.06042 s
	  // 0.00001     1/0.361      660 NDS     9201 BTS/s        240103 CKS 0.07017 s
	  // 0.000001    1/0.361     1007 NDS     9691 BTS/s        346201 CKS 0.10252 s
	  // 0.000000001 1/0.361     1673 NDS    10094 BTS/s        544769 CKS 0.16436 s


	  // sac + impact/wdeg
	  // 1           1/0.361     1280 NDS     5354 BTS/s        341568 CKS 0.10538 s
	  // 0.0001      1/0.361     1221 NDS     5455 BTS/s        319567 CKS 0.09264 s
	  // 0           1/0.361     1678 NDS     6865 BTS/s        463694 CKS 0.14020 s


	  // impact
	  // 1           1/0.361     1957 NDS     9534 BTS/s        661834 CKS 0.20370 s
	  // 0.0001      1/0.361     1004 NDS     9460 BTS/s        369642 CKS 0.10464 s
	  // 0.0000001   1/0.361     2601 NDS     9004 BTS/s        898515 CKS 0.28732 s


	  // sac + impact
	  // 0.0001      1/0.361     1526 NDS     7019 BTS/s        428557 CKS 0.11546 s






	  // dom/wdeg
	  //             1/0.68     25680 NDS     9001 BTS/s      13660881 CKS 2.84940 s

	  // probe (100x30) + dom/wdeg
	  //             1/0.68     24185 NDS     8760 BTS/s      11673523 CKS 2.62100 s

	  // probe (100x50) + dom/wdeg
	  //             1/0.68     25070 NDS     9751 BTS/s      11574613 CKS 2.44440 s

	  // probe (100x75) + dom/wdeg
	  //             1/0.68    26578 NDS    10049 BTS/s      11686517 CKS  2.5202 s

	  // probe (200x30) + dom/wdeg
	  //             1/0.68     27658 NDS     9003 BTS/s      12353985 CKS 2.80400 s

	  // probe (200x50) + dom/wdeg
	  //             1/0.68     30975 NDS    10208 BTS/s      12946528 CKS 2.79400 s

	  // impact/wdeg
	  // 0.0001      1/0.68     51859 NDS     8332 BTS/s      23709083 CKS 6.22020 s


	  // sac + impact/wdeg
	  // 0.0001      1/0.68     36246 NDS     8047 BTS/s      17021756 CKS 4.38840 s


	  // impact
	  // 1           1/0.68    174116 NDS     9196 BTS/s      73627047 CKS  18.928 s
	  // 0.0001      1/0.68     66433 NDS     8519 BTS/s      31399406 CKS  7.7944 s
	  // 


	  // sac + impact
	  // 0.0001     1/0.68      49709 NDS     8506 BTS/s      22890752 CKS  5.7344 s


	  // probe + impact/wdeg
	  // 0.0001     1/0.68      32999 NDS     8692 BTS/s      14480867 CKS  3.6554 s


	  varmin[i] = s->variables[i]->minCapacity();

	  decision_impact[i] -= varmin[i];
	  decision_count [i] -= varmin[i];
	}
    }

}

void WeighterImpact::notifyFailure( Constraint *con )
{

  if( level > init_level && needUpdate )
    {
      // as a result of the decision, the problem is emptied
      // the inverse impact therefore is 0.
      int x = decision[level]->id;
      int v;

      if( isRange[x] ) {
	
	//std::cout << "failure" << std::endl;

	v = ((VariableRange*)(decision[level]))->whichBound();

	//std::cout << "update " << (v ? "upper" : "lower") << " bound's impact" << std::endl;

      } else { 
	v = decision[level]->value();
      }
      decision_impact[x][v] = 
	(decision_impact[x][v] * decision_count[x][v]) / ++decision_count[x][v];
      needUpdate = false;
    }
}

void WeighterImpact::notifySuccess( )
{

  // Compute the reduction Sa/Sb
  double reduction = 1;
  int v;
  unsigned int x;

  for(x=0; x<nbvars; ++x)
    if( X[x]->domsize() < (int)(domain_size[x]) ) {
      reduction *= ((double)(X[x]->domsize()) / (double)(domain_size[x]));
    }

  x = decision[level-1]->id;

  if( isRange[x] ) {
    //std::cout << "success" << std::endl;

    v = ((VariableRange*)(decision[level-1]))->whichBound();

    //std::cout << "update " << (v ? "upper" : "lower") << " bound's impact" << std::endl;

  }
  else 
    v = decision[level-1]->value();

  decision_impact[x][v] = 
    (((decision_impact[x][v]) * (decision_count[x][v])) + reduction) / ++(decision_count[x][v]);

  needUpdate = false;
}

void WeighterImpact::notifyChoice( )
{
  VariableInt **var_iterator;
  nbvars=0;
  for( var_iterator = first; var_iterator != last; ++var_iterator ) {
    // store the var and its size
    X[nbvars] = (*var_iterator);
    domain_size[nbvars] = X[nbvars]->domsize();
    ++nbvars;
  }
  needUpdate =true;
} 


bool WeighterSAC::isNotSac() { 
  return ins;
}

WeighterSAC::~WeighterSAC() 
{
  delete [] SACdomain;
  delete [] isRange;
}

WeighterSAC::WeighterSAC(Solver* s) 
  : Weighter(s)
{
  complete = false;
  ins = true;
  solver = s;
  ilength = s->length;
  needPruning = false;
  SACdomain = new BitSet[ilength];
  isRange = new bool[ilength];
  sacvals = 0;
  domsize = 0;
  for(int i=0; i<ilength; ++i) 
    {
      isRange[i] = ( s->variables[i]->getType() == VariableInt::RANGE );

      if( isRange[i] )
	{
	  SACdomain[i].init( 0, 1, BitSet::empt );
	  domsize += 2;
	}
      else 
	{
	  SACdomain[i].init( s->variables[i]->min(),
			     s->variables[i]->max(),
			     BitSet::empt );
	  domsize += s->variables[i]->domsize();
	  if( s->variables[i]->isGround() ) {
	    ++sacvals;
	  }
	}
    }
}


void WeighterSAC::notifyChoice()
{
  int i, j;
  if( level == init_level && needPruning )
    {
      domsize = 0;
      if( complete ) 
	{
	  j=solver->decision[1]->id;
	  //sacvals = SACdomain[j].size();
	  sacvals = 0;
	  for(i=0; i<ilength; ++i)
	    {
	      domsize += solver->variables[i]->domsize();
	      //if( i != j ) 
	      SACdomain[i].clear();
	    }
	}
      else
	{
	  for(i=0; i<ilength; ++i)
	    domsize += solver->variables[i]->domsize();
	}
    }
  needPruning = false;
}


void WeighterSAC::notifyFailure( Constraint *con ) 
{
  int j, k;
  if( level > init_level ) 
    {
      ins =false;

      VariableInt **var_iterator = solver->sequence,
	**last = solver->future;
      
      while( var_iterator < last )
	{
	  j = (*var_iterator)->id;
	  if( isRange[j] )
	    k = ((VariableRange*)(*var_iterator))->whichBound();
	  else
	    k = (*var_iterator)->value();

	  ++var_iterator;
	  if( SACdomain[j].member(k) ) continue;

	  SACdomain[j].insert(k);
	  ins = true;
	  ++sacvals;
	}
      ins = (ins && (domsize > sacvals));
    }
  else if( level == init_level )
    {
      ++solver->FAILLIMIT;
      needPruning = true;
    }
}



WeighterISAC::WeighterISAC(Solver* s)
 : WeighterSAC(s)
{
}

WeighterISAC::~WeighterISAC()
{
}

void WeighterISAC::notifyFailure( Constraint *con )
{
 int j, k;
 if( level > init_level )
   {
     ins =false;

     for(int i=1; i<(level + 1); i++ )
       {
         j=solver->decision[i]->id;
         if( isRange[j] )
           k = ((VariableRange*)(solver->decision[i]))->whichBound();
         else
           k = solver->decision[i]->value();

         if( SACdomain[j].member(k) ) continue;

         SACdomain[j].insert(k);
         ins = true;
         ++sacvals;
       }
     ins = (ins && (domsize > sacvals));
   }
 else if( level == init_level )
   {
     ++solver->FAILLIMIT;
     needPruning = true;
   }
}


DVOSingletonAC::DVOSingletonAC(Solver* s, WeighterSAC* w) 
  : DVO(s), SACdomain(w->SACdomain), isRange(w->isRange)
{
}

inline VariableInt* DVOSingletonAC::select() 
{
  VariableInt *var, **var_iterator=first;
  int i;

  do {
    var = *var_iterator;
    i=var->id;
  } while( 
	  ((isRange[i] && SACdomain[i].empty())
	   ||
	   (!isRange[i] && var->included( SACdomain[i] )))
	  && 
	  ++var_iterator != last );
  
  return var;
}



WeighterRestartNogood::WeighterRestartNogood( Solver* s ) 
  : Weighter(s), decision(s->decision.stack_)
{
  // nBranches = new int[s->numvars+1];
  choices = new int[s->numvars+1];
  lvl = 0;
}

WeighterRestartNogood::~WeighterRestartNogood() 
{  
  delete [] choices;
  // delete [] nBranches;
}

void WeighterRestartNogood::notifyFailure( Constraint *con )
{

  lvl = level;

  //if( level > 1 ) {
  if( level > init_level ) {
    int failed = choices[level];
    while( path.back() != failed ) {
      path.pop();
    }
  } else {
    path.clear();
  }

//   if( level ) {
//     VariableInt *failed = decision[level];
//     while( path.back() != failed ) {
//       path.pop();
//       polarity.pop();
//     }
//   } else {
//     path.clear();
//     polarity.clear();
//   }

//   if( level ) {
//     Atom failed = decision[level]->id+1;
//     while( atom(path.back()) != failed ) 
//       path.pop();
//   } else {
//     path.clear();
//   }

// //   if( nBranches[path.size-1]++ ) {
// //     path.pop();
// //     ++nBranches[path.size-1];
// //   }

//   std::cout <<setw(2)<< level << "\t";
//   std::cout.flush();
//   for(int i=1; i<=level; ++i) {
//     std::cout << " " << setw(2) << (decision[i]->id+1);
//   }
//   std::cout << std::endl;
//   std::cout <<setw(2)<< (path.size) << "\t";
//   std::cout.flush();
//   for(int i=0; i<path.size; ++i) {
//     std::cout << " " << setw(2) << path[i]->id+1;
//   }
//   std::cout << std::endl;
// //   std::cout <<setw(2)<< (path.size) << "\t";
// //   std::cout.flush();
// //   for(int i=0; i<path.size; ++i) {
// //     std::cout << " " << setw(2) << nBranches[i];
// //   }
// //   std::cout << std::endl;
}

void WeighterRestartNogood::notifyChoice()
{
  
  //nBranches[path.size] = 0;
  //path.push( decision[level] );  
  //polarity.push( decision[level]->equal(0) ? 1 : -1 );

  if(decision[level]) {
    int id = decision[level]->id+1;
    choices[level] = ( decision[level]->equal(0) ? id : -id );
  } else {
    choices[level] = 0;
  }
  
  path.push( choices[level] );

//   std::cout << std::endl <<setw(2)<< level << "\t";
//   std::cout.flush();
//   for(int i=1; i<=level; ++i) {
//     std::cout << " " << setw(2) << (decision[i]->id+1);
//   }
//   std::cout << std::endl;
//   std::cout <<setw(2)<< (path.size) << "\t";
//   std::cout.flush();
//   for(int i=0; i<path.size; ++i) {
//     std::cout << " " << setw(2) << path[i]->id+1;
//   }
//   std::cout << std::endl;
// //   std::cout <<setw(2)<< (path.size) << "\t";
// //   std::cout.flush();
// //   for(int i=0; i<path.size; ++i) {
// //     std::cout << " " << setw(2) << nBranches[i];
// //   }
// //   std::cout << std::endl;
}


void WeighterRestartNogood::notifyRestart() 
{


//   std::cout << std::endl <<std::setw(3)<< lvl << "   ";
//   std::cout.flush();
//   for(int i=1; i<lvl; ++i) {
//     std::cout << " " << std::setw(2) << choices[i];
//   }
//   std::cout << std::endl;
//   std::cout << std::setw(3)<< (path.size) << "   ";
//   std::cout.flush();
//   for(int i=0; i<path.size; ++i) {
//     std::cout << " " << std::setw(2) << path[i];
//   }
//   std::cout << std::endl ;



  int i=1, j=0, n=lvl, m=path.size;


  //decision[n] = NULL;
  choices[n] = 0;

//   // find first disagreemnet between path and decision
//   while( i<n && decision[i] == path[j] ) {
//     ++i;
//     ++j;
//     // loop through all right 
//     while( j<m && path[j] != decision[i] ) {
//       std::cout << "      ";
//       for(int x=1; x<i; ++x)
// 	std::cout << " " << setw(2) << polarity[x]decision[x]->id+1;
//       std::cout << " " << setw(2) << path[j]->id+1 << std::endl;
//       ++j;
//     }
//   }
//   std::cout << std::endl;


  Vector< Literal > clause;
  while( j<m && path[j] != choices[i] )
    ++j;

  // find first disagreemnet between path and decision
  while( i<n && choices[i] == path[j] ) {
    ++i;
    ++j;
    // loop through all right 
    while( j<m && path[j] != choices[i] ) {
      clause.clear();
      //std::cout << "      ";
      for(int x=1; x<i; ++x) {
	clause.push(choices[x]);
	//std::cout << " " << setw(2) << choices[x];
      }
      clause.push( path[j] );
      //std::cout << " " << setw(2) << path[j] << std::endl;
      ++j;

      sat->addClause( sat->learnt, clause, sat->stats.learnt_avg_size );
    }
  }
  //std::cout << std::endl;


//   std::cout << "Assumption: ";
//   for(i=0; i<sat->assumptions.size; ++i)
//     std::cout << sat->polarity[sat->assumptions[i]] << " ";
//   std::cout << std::endl;

//   std::cout << "Assumption: ";
//   for(i=0; i<sat->assumptions.size; ++i) {
//     sat->X[sat->assumptions[i]]->print( std::cout ) ;
//     std::cout << " ";
//   }
//   std::cout << std::endl;


//   if( path.size )
//  exit(0);

  path.clear();
  lvl = 0;
}


WeighterRestartGenNogood::WeighterRestartGenNogood( Solver* s ) 
  : Weighter(s), decision(s->branching_decision.stack_)
{
  bad_choices = new Vector< Decision >[s->numvars+1];
  bad_choices -= s->init_level;
  
  depth = s->init_level;
}

WeighterRestartGenNogood::~WeighterRestartGenNogood() 
{  
  bad_choices += init_level;
  delete [] bad_choices;
}

void WeighterRestartGenNogood::forget(const int l) 
{
  base->forget(l);
}

void WeighterRestartGenNogood::reinit()
{
  for(int i=init_level+2; i<=depth; ++i) {
    bad_choices[i].clear();
  }
}

void WeighterRestartGenNogood::notifyChoice()
{
  //std::cout << " choice " << level << " " ;
  //decision[level].print(std::cout);
  //std::cout << std::endl;

  bad_choices[level+1].clear();
}

void WeighterRestartGenNogood::notifyFailure( Constraint *con )
{
//      std::cout << " failure " << level << " " ;

//    if(level > init_level+1) {
//      decision[level].print(std::cout);
//      std::cout << std::endl;
//    }

  bad_choices[level].push(decision[level]);
  depth = level;
  if(depth < init_level+2) {
    bad_choices[init_level+2].clear();
  }
}

void WeighterRestartGenNogood::notifyRestart() 
{
  
  //  std::cout << level << " " << init_level << " " << depth << " restart!!" << std::endl;
  //exit(1);

//   for(int i=init_level+2; i<=depth; ++i) {
//     std::cout << "bad (" << (bad_choices[i].size) << ") ";
//     for(int j=0; j<bad_choices[i].size; ++j) {
//       bad_choices[i][j].print(std::cout);
//       std::cout << " ";
//     }
//     if(i<depth) {
//       std::cout << "good: ";
//       decision[i].print(std::cout);
//     }
//     std::cout << std::endl;
//   }
  
  Vector< Decision > learnt;
  int d = init_level+2;
  while(d <= depth) {
    for(int i=0; i<bad_choices[d].size; ++i) {
      learnt.push(bad_choices[d][i]);
      learnt.back().revert();

//       std::cout << "add ";
//       for(int j=0; j<=(d-init_level-2); ++j) {
//    	learnt[j].print(std::cout);
//    	std::cout << " ";
//       }
//       std::cout << std::endl;
      
      base->add(learnt);

      learnt.pop();
    }
    learnt.push(decision[d]);
    learnt.back().revert();

    bad_choices[d].clear();
    ++d;
  }

  depth = init_level+2;

//   int i=1, j=0, n=lvl, m=path.size;
//   choices[n] = 0;
//   Vector< Literal > clause;

//   while( j<m && path[j] != choices[i] )
//     ++j;

//   // find first disagreemnet between path and decision
//   while( i<n && choices[i] == path[j] ) {
//     ++i;
//     ++j;
//     // loop through all right 
//     while( j<m && path[j] != choices[i] ) {
//       clause.clear();
//       for(int x=1; x<i; ++x) {
// 	clause.push(choices[x]);
//       }
//       clause.push( path[j] );
//       ++j;

//       //base->add( );
//     }
//   }
//   path.clear();
//   lvl = 0;

}



WeighterSwitchPromise::WeighterSwitchPromise( Solver* s, const int t ) 
  : Weighter(s), sequence(s->sequence)
{
  length = s->length;
  n_restart = 0;
  threshold = t;
}

WeighterSwitchPromise::~WeighterSwitchPromise() 
{  
}

void WeighterSwitchPromise::notifyRestart() 
{
  //std::cout << "switch restart? " << n_restart << " " << threshold << std::endl;

  if( ++n_restart == threshold ) {
    for(int i=0; i<length; ++i)
      ((ValSelectorSOSP*)(sequence[i]->branch))->mode = 1;
  }
}

WeighterSAT::WeighterSAT( Solver* s, ConstraintClauseBase *c ) 
  : Weighter(s), sat(c), first(s->future)
{
  //branch = LEFT;
  // dec = 0;
  //nright = 0;
//   s->binds( fail );
//   fail.setValue( 0 );
}

WeighterSAT::~WeighterSAT() 
{
}

void WeighterSAT::notifyFailure( Constraint *con )
{
  
  sat->analyzeConflict();

}

// //   if( branch == LEFT ) {
// //     if( con == sat ) {
      
// // #ifdef _DEBUGNOGOOD
      
// //       std::cout << "FAIL ON A DECISION" << std::endl;
      
// //       sat->printDecisions( std::cout, 0 );
// // #endif
      
// //       //sat->printDecisions( std::cout, 0 );

// //       //Literal p;

// //       sat->backjump();

// // 	 //sat->btLevel = sat->analyze( sat->conflict, p );
      
// // //       //     std::cout << p << " because: ";
// // //       //     sat->printClause( std::cout, sat->learnt.back() );
// // //       //     std::cout << std::endl;
      
// //  //       sat->reason[atom(p)] = sat->learnt.back();
// // //        //fail = p;
// // //        //fail.push(p);
       
// // //        std::cout << "c left  fail " << p  << " because ";
// // //        sat->printClause( std::cout, sat->learnt.back() );
// // //        std::cout << std::endl;


// //       branch = RIGHT;

// //     } 
// //   }
// // //   else {
    
// // // //     int q = sat->solver->decision[*(sat->level)+1]->id+1;
// // // //     std::cout << "c right fail " << q << " ";
// // // //     std::cout << std::endl;
// // // //  //    //sat->printDecisions( std::cout, 0 );
// // // // //     if( sat->reason[q] ) {
// // // // //       sat->resolve( sat->reason[q], q );
// // // // //     }    
// // //   }
// }

void WeighterSAT::notifyChoice()
{

//   std::cout << "Notify ";
//   first[0]->print( std::cout );
//   std::cout << std::endl;

//   int p = first[0]->id;

//   //if( first[0]->equal(0) ) p = -p;
//   //sat->makeDecision( p );
//   //sat->save();


//   std::cout << std::endl << "before the notification: " << std::endl;
//   sat->printDecisions( std::cout );

  //sat->notifyChoice();

  //std::cout << "after the notification: " << std::endl;
  //sat->printDecisions( std::cout );

  sat->is_choice = true;
  //sat->decayActivity();


//   //fail = p;
//   std::cout << "choice " << p << " " << *(sat->level) << std::endl;


//   //branch = LEFT;
//   //nright = 0;
//   //dec = 0;
}


/**********************************************
 * DVOs Dynamic Variable Ordering
 **********************************************/

DVO::DVO(Solver *s, const int l) 
  : variables(s->variables.stack_),
    first(s->future), 
    last(s->empty), 
    decision(s->decision.stack_), 
    level(s->level)
{
//   VariableInt **var_iterator;
//   for(var_iterator = first; var_iterator != last; ++var_iterator)
//     {
//     }
  limit = l;
  //verbosity = 0;
}

PredicateDisjunctive** DVO::get_disjuncts() { return _garbage_disjuncts; }

inline VariableInt* DVONoOrder::select() 
{
  return *first;
}

DVOLexicographic::~DVOLexicographic()
{
  delete [] sequence;
}

DVOLexicographic::DVOLexicographic(Solver* s) : DVO(s) 
{ 
  //last.init(s,0); 
  s->binds( lastIdx );
  lastIdx.setValue( 0 );
  
  if(s->length) {
    sequence = new VariableInt*[s->length];
    std::memcpy(sequence, s->sequence, (s->length)*sizeof(VariableInt*));
  }

//   for(int i=0; i<s->length; ++i)
//     {
//       sequence[i]->printshort( std::cout );
//       std::cout << std::endl;
//     }
}

inline VariableInt* DVOLexicographic::select() 
{
  int i=lastIdx;  
//   while( !variables[i]->isLinked() ) ++i;
//   lastIdx = i;
//   return variables[i];

  while( !sequence[i]->isLinked() ) {
    ++i;
  }
   lastIdx = i;
  return sequence[i];
}


DVOAntiLex::~DVOAntiLex()
{
  delete [] sequence;
}

DVOAntiLex::DVOAntiLex(Solver* s) : DVO(s) 
{ 

  s->binds( lastIdx );
  lastIdx.setValue( 0 );
  
  if(s->length) {
    sequence = new VariableInt*[s->length];
    for(int i=0; i<s->length; ++i)
      sequence[i] = s->sequence[s->length-i-1];
    //std::memcpy(sequence, s->sequence, (s->length)*sizeof(VariableInt*));
  }

//   for(int i=0; i<s->length; ++i)
//     {
//       sequence[i]->printshort( std::cout );
//       std::cout << std::endl;
//     }
}

inline VariableInt* DVOAntiLex::select() 
{
  int i=lastIdx;  
//   while( !variables[i]->isLinked() ) ++i;
//   lastIdx = i;
//   return variables[i];

  while( !sequence[i]->isLinked() ) {
    ++i;
  }
  lastIdx = i;
  return sequence[i];
}

OSP::~OSP() {
}

PredicateDisjunctive** OSP::get_disjuncts(Solver *s) {
  // collect disjuncts
  
  PredicateDisjunctive **disjunct = new PredicateDisjunctive*[s->variables.size];
  std::fill(disjunct, disjunct+(s->variables.size), (PredicateDisjunctive *)NULL);
  _garbage_disjuncts = disjunct;
  int i, j, n = s->constraints.size;
  Constraint **cons = s->constraints.stack_;
  VariableInt *x;
  for(i=0; i<n; ++i)
    if( cons[i]->arity == 3 ) {
      x = cons[i]->scope[2];
      j = x->id;
      disjunct[j] = (PredicateDisjunctive*)(cons[i]);
      if(promise == 1) {
	//std::cout << "PROMISE" << std::endl;
	delete x->branch;
	x->branch = new ValSelectorLOSP( x, disjunct[j] );
      } else if(promise == -1) {
	//std::cout << "ANTI" << std::endl;
	delete x->branch;
	x->branch = new ValSelectorMOSP( x, disjunct[j] );
      }
    }
  return disjunct;
}

DVO* OSP::extract( Solver* s )
{
  s->setLearner( Weighter::WDG );
  if( size > 1 ) {
    int i;
    switch(strategy) {
    case DOMAIN_O_NOT: {
      GenericSchedulingRandomDVO<VarSelectorOSP_Domain> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_Domain>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOMAIN_P_TWEIGHT: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DomainWeight> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DomainWeight>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_BOOLWEIGHT: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DoBoolWeight> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DoBoolWeight>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_TASKWEIGHT: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DoTaskWeight> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DoTaskWeight>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_BOOLTASKWEIGHT: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DoBoolTaskWeight> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DoBoolTaskWeight>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOMAIN_O_NOTTYPE: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DomainType> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DomainType>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_BOOLWEIGHTTYPE: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DoBoolWeightType> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DoBoolWeightType>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_TASKWEIGHTTYPE: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DoTaskWeightType> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DoTaskWeightType>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_BOOLTASKWEIGHTTYPE: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DoBoolTaskWeightType> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DoBoolTaskWeightType>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i)
	var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    default: {
      GenericSchedulingRandomDVO<VarSelectorOSP_DoWeakWeight> *var_heuristic = 
	new GenericSchedulingRandomDVO<VarSelectorOSP_DoWeakWeight>(s, size);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      for(i=0; i<=size; ++i) {
	var_heuristic->bests[i].disjuncts = disjunct;
	var_heuristic->bests[i].weak_ = strategy;
      }
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->current.weak_ = strategy;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    }
    // return new GenericRandomDVO<VarSelectorOSP>(s, size);
  } else {
    switch(strategy) {
    case DOMAIN_O_NOT: {
      GenericSchedulingDVO<VarSelectorOSP_Domain> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_Domain>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOMAIN_P_TWEIGHT: {
      GenericSchedulingDVO<VarSelectorOSP_DomainWeight> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DomainWeight>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_BOOLWEIGHT: {
      GenericSchedulingDVO<VarSelectorOSP_DoBoolWeight> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DoBoolWeight>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_TASKWEIGHT: {
      GenericSchedulingDVO<VarSelectorOSP_DoTaskWeight> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DoTaskWeight>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_BOOLTASKWEIGHT: {
      GenericSchedulingDVO<VarSelectorOSP_DoBoolTaskWeight> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DoBoolTaskWeight>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOMAIN_O_NOTTYPE: {
      GenericSchedulingDVO<VarSelectorOSP_DomainType> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DomainType>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_BOOLWEIGHTTYPE: {
      GenericSchedulingDVO<VarSelectorOSP_DoBoolWeightType> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DoBoolWeightType>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_TASKWEIGHTTYPE: {
      GenericSchedulingDVO<VarSelectorOSP_DoTaskWeightType> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DoTaskWeightType>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    case DOM_O_BOOLTASKWEIGHTTYPE: {
      GenericSchedulingDVO<VarSelectorOSP_DoBoolTaskWeightType> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DoBoolTaskWeightType>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    default: {
      GenericSchedulingDVO<VarSelectorOSP_DoWeakWeight> *var_heuristic = 
	new GenericSchedulingDVO<VarSelectorOSP_DoWeakWeight>(s);
      PredicateDisjunctive** disjunct = get_disjuncts(s);
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->best.weak_ = strategy;
      var_heuristic->current.weak_ = strategy;
      var_heuristic->the_disjuncts = disjunct;
      return var_heuristic;
    }
    }
  }
  return NULL;
}


void PFSP::get_disjuncts(Solver* s,
			 PredicateDisjunctive***& disjunct,
			 int*& sdegree) {
  int i, n = s->variables.size;
  VariableInt *x;
  MistralNode<Constraint*> *nd;

  disjunct = new PredicateDisjunctive**[n];
  sdegree = new int[n];

  //_garbage_disjuncts = disjunct;
  for(i=0; i<n; ++i) {
    x = s->variables[i];

    sdegree[i] = 0;
    disjunct[i] = new PredicateDisjunctive*[x->degree];

    if(x->getType() == VariableInt::BOOL) {
      nd = x->constraintsOnValue();
      while( nextNode(nd) ) {
	if(nd->elt->arity == 3) 
	  disjunct[i][sdegree[i]++] = ((PredicateDisjunctive*)(nd->elt));
      }

       if(promise == 1) {
	 x->branch = new ValSelectorPFSP( x, disjunct[i], sdegree[i] );
       }
       // else if(promise == -1)
// 	x->branch = new ValSelectorMOSP( x, disjunct[j] );
    }
  }
}

DVO* PFSP::extract( Solver* s )
{
  s->setLearner( Weighter::WDG );
  if( size > 1 ) {
    int i;
    GenericRandomDVO<VarSelectorPFSP> *var_heuristic = 
      new GenericRandomDVO<VarSelectorPFSP>(s, size);

    PredicateDisjunctive*** disjunct;
    int* sdegree;
    get_disjuncts(s, disjunct, sdegree);

    for(i=0; i<=size; ++i) {
      var_heuristic->bests[i].disjuncts = disjunct;
      var_heuristic->bests[i].sdegree = sdegree;
    }
    var_heuristic->current.disjuncts = disjunct;
    var_heuristic->current.sdegree = sdegree;
    return var_heuristic;
    
  } else {

    GenericDVO<VarSelectorPFSP> *var_heuristic = 
      new GenericDVO<VarSelectorPFSP>(s);

    PredicateDisjunctive*** disjunct;
    int* sdegree;
    get_disjuncts(s, disjunct, sdegree);

    var_heuristic->best.disjuncts = disjunct;
    var_heuristic->best.sdegree = sdegree;
    var_heuristic->current.disjuncts = disjunct;
    var_heuristic->current.sdegree = sdegree;
    return var_heuristic;
    
  }
  return NULL;
}


DVO* OSPSAT::extract( Solver* s )
{
  s->setLearner( Weighter::WDG );
  if( size > 1 ) {
    GenericRandomDVO<VarSelectorOSPSAT> *var_heuristic = new GenericRandomDVO<VarSelectorOSPSAT>(s, size);
      
      // collect disjuncts
      PredicateDisjunct **disjunct = new PredicateDisjunct*[s->variables.size];
      int i, j, n = s->constraints.size;
      Constraint **cons = s->constraints.stack_;
      VariableInt *x;
      for(i=0; i<n; ++i)
	if( cons[i]->arity == 3 ) {
	  x = cons[i]->scope[2];
	  j = x->id;
	  disjunct[j] = (PredicateDisjunct*)(cons[i]);
	  //x->branch = new ValSelectorOSP( x, disjunct[j] );
	}
      for(i=0; i<=size; ++i) {
	var_heuristic->bests[i].disjuncts = disjunct;
	var_heuristic->bests[i].activity = s->sat->activity;
      }
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->current.activity = s->sat->activity;

      return var_heuristic;
      // return new GenericRandomDVO<VarSelectorOSP>(s, size);
  } else {
      GenericDVO<VarSelectorOSPSAT> *var_heuristic = new GenericDVO<VarSelectorOSPSAT>(s);
      
      // collect disjuncts
      PredicateDisjunct **disjunct = new PredicateDisjunct*[s->variables.size];
      int i, j, n = s->constraints.size;
      Constraint **cons = s->constraints.stack_;
      VariableInt *x;
      for(i=0; i<n; ++i)
	if( cons[i]->arity == 3 ) {
	  x = cons[i]->scope[2];
	  j = x->id;
	  disjunct[j] = (PredicateDisjunct*)(cons[i]);
	  //x->branch = new ValSelectorOSP( x, disjunct[j] );
	}
      var_heuristic->best.disjuncts = disjunct;
      var_heuristic->best.activity = s->sat->activity;
      var_heuristic->current.disjuncts = disjunct;
      var_heuristic->current.activity = s->sat->activity;

      return var_heuristic;
    }
}


DVO* FPP::extract( Solver* s )
{
  s->setLearner( Weighter::WDG );
  if( size > 1 ) {
    GenericRandomDVO<VarSelectorFPP> *var_heuristic = new GenericRandomDVO<VarSelectorFPP>(s, size);
      
    // collect disjuncts
    PredicateLess **prec = new PredicateLess*[s->variables.size];
    PredicateUpperBound **inter = new PredicateUpperBound*[s->variables.size];
    int *isInterval = new int[s->variables.size];
    int i, j, n = s->constraints.size;
    Constraint **cons = s->constraints.stack_;
    VariableInt *x;
    for(i=0; i<n; ++i)
      if( cons[i]->arity == 3 ) {
	x = cons[i]->scope[2];
	j = x->id;
	prec[j] = (PredicateLess*)(cons[i]);
	isInterval[j] = 0;
      } else if( cons[i]->arity == 2 ) {
	x = cons[i]->scope[1];
	j = x->id;
	inter[j] = (PredicateUpperBound*)(cons[i]);
	isInterval[j] = 1;
      }
    for(i=0; i<=size; ++i) {
      var_heuristic->bests[i].precs = prec;
      var_heuristic->bests[i].inters = inter;
      var_heuristic->bests[i].isInterval = isInterval;
    }
    var_heuristic->current.precs = prec;
    var_heuristic->current.inters = inter;
    var_heuristic->current.isInterval = isInterval;
    
    return var_heuristic;
    // return new GenericRandomDVO<VarSelectorFPP>(s, size);
  } else {
    GenericDVO<VarSelectorFPP> *var_heuristic = new GenericDVO<VarSelectorFPP>(s);
    
    // collect disjuncts
    PredicateLess **prec = new PredicateLess*[s->variables.size];
    PredicateUpperBound **inter = new PredicateUpperBound*[s->variables.size];
    int *isInterval = new int[s->variables.size];
    int i, j, n = s->constraints.size;
    Constraint **cons = s->constraints.stack_;
    VariableInt *x;
    for(i=0; i<n; ++i)
      if( cons[i]->arity == 3 ) {
	x = cons[i]->scope[2];
	j = x->id;
	prec[j] = (PredicateLess*)(cons[i]);
      } else if( cons[i]->arity == 2 ) {
	x = cons[i]->scope[1];
	j = x->id;
	inter[j] = (PredicateUpperBound*)(cons[i]);
	isInterval[j] = 1;
      }
    var_heuristic->best.precs = prec;
    var_heuristic->current.precs = prec;
    var_heuristic->best.inters = inter;
    var_heuristic->current.inters = inter;
    var_heuristic->best.isInterval = isInterval;
    var_heuristic->current.isInterval = isInterval;
    
    
    return var_heuristic;
  }
}


DVO* DomOverWLDeg::extract( Solver* s )
{
  s->setLearner( Weighter::WLD );
  if( size > 1 )
    return new GenericRandomDVO<VarSelectorDomainOverWeight>(s, size);
  else 
    return new GenericDVO<VarSelectorDomainOverWeight>(s);
}


DVO* DomOverWDeg::extract( Solver* s )
{
  s->setLearner( Weighter::WDG );
  if( size > 1 )
    return new GenericRandomDVO<VarSelectorDomainOverWeight>(s, size);
  else 
    return new GenericDVO<VarSelectorDomainOverWeight>(s);
}



void linkImpact( Solver* s, double**** bdi, double*** cdi, const int size ) 
{
  WeighterImpact *w = ((WeighterImpact*)(s->setLearner( Weighter::IPT )));
  for(int i=0; i<s->length; ++i)
    {
      if( !s->sequence[i]->branch && 
	  s->sequence[i]->getType() != VariableInt::RANGE )
	//s->sequence[i]->branch = new ValSelectorWeight( s->sequence[i], w->decision_impact[i] );
	s->sequence[i]->branch = new ValSelectorRand
( s->sequence[i] );
      //s->sequence[i]->branch = new ValSelectorMin( s->sequence[i] );
    }
  for(int i=0; i<size; ++i) 
    *(bdi[i]) = w->decision_impact;
  *cdi = w->decision_impact;
  //  std::cout << "o impact " << (w->decision_impact) << std::endl;
}

DVO* Impact::extract( Solver* s )
{
  DVO *h;
  if( !size ) size=1;
  double ***bdi[size];
  double ***cdi;
  //GenericRandomDVO<VarSelectorImpact> *hr;
  if( size > 1 ) {
    GenericRandomDVO<VarSelectorImpact> *hr = 
      new GenericRandomDVO<VarSelectorImpact>(s, size);
    for(int i=0; i<size; ++i)
      bdi[i] = &(hr->bests[i].decision_impact);
    cdi = &(hr->current.decision_impact);
    h = hr;
  } else {
    GenericDVO<VarSelectorImpact> *hn = 
      new GenericDVO<VarSelectorImpact>(s);
    bdi[0] = &(hn->best.decision_impact);
    cdi = &(hn->current.decision_impact);
    h = hn;
  }
  linkImpact( s, bdi, cdi, size );

//   std::cout << "r impact " << (((GenericDVO<VarSelectorImpact> *)h)->best.decision_impact) << std::endl;
//   std::cout << "r impact " << (((GenericDVO<VarSelectorImpact> *)h)->current.decision_impact) << std::endl;

  return h;
}

DVO* ImpactOverDeg::extract( Solver* s )
{
  DVO *h;
  if( !size ) size=1;
  double ***bdi[size];
  double ***cdi;
  if( size > 1 ) {
    GenericRandomDVO<VarSelectorImpactOverDegree> *hr = 
      new GenericRandomDVO<VarSelectorImpactOverDegree>(s, size);
    for(int i=0; i<size; ++i)
      bdi[i] = &(hr->bests[i].decision_impact);
    cdi = &(hr->current.decision_impact);
    h = hr;
  } else {
    GenericDVO<VarSelectorImpactOverDegree> *hn = 
      new GenericDVO<VarSelectorImpactOverDegree>(s);
    bdi[0] = &(hn->best.decision_impact);
    cdi = &(hn->current.decision_impact);
    h = hn;
  }
  linkImpact( s, bdi, cdi, size );
  return h;
}

DVO* ImpactOverWLDeg::extract( Solver* s )
{
  DVO *h;
  if( !size ) size=1;
  double ***bdi[size];
  double ***cdi;
  if( size > 1 ) {
    GenericRandomDVO<VarSelectorImpactOverWeight> *hr = 
      new GenericRandomDVO<VarSelectorImpactOverWeight>(s, size);
    for(int i=0; i<size; ++i)
      bdi[i] = &(hr->bests[i].decision_impact);
    cdi = &(hr->current.decision_impact);
    h = hr;
  } else {
    GenericDVO<VarSelectorImpactOverWeight> *hn = 
      new GenericDVO<VarSelectorImpactOverWeight>(s);
    bdi[0] = &(hn->best.decision_impact);
    cdi = &(hn->current.decision_impact);
    h = hn;
  }
  linkImpact( s, bdi, cdi, size );
  return h;
}

DVO* ImpactOverWDeg::extract( Solver* s )
{
  DVO *h;
  if( !size ) size=1;
  double ***bdi[size];
  double ***cdi;
  if( size > 1 ) {
    GenericRandomDVO<VarSelectorImpactOverWeight> *hr = 
      new GenericRandomDVO<VarSelectorImpactOverWeight>(s, size);
    for(int i=0; i<size; ++i)
      bdi[i] = &(hr->bests[i].decision_impact);
    cdi = &(hr->current.decision_impact);
    h = hr;
  } else {
    GenericDVO<VarSelectorImpactOverWeight> *hn = 
      new GenericDVO<VarSelectorImpactOverWeight>(s);
    bdi[0] = &(hn->best.decision_impact);
    cdi = &(hn->current.decision_impact);
    h = hn;
  }
  linkImpact( s, bdi, cdi, size );
  return h;
}

DVO* Probedvo::extract( Solver* s )
{
  s->setLearner( ltype );

  switch( htype ) {
  case Weighter::NO : return new DVORandom(s);
  case Weighter::WDG : return new GenericDVO<VarSelectorDomainOverWeight>(s);
  case Weighter::WLD : return new GenericDVO<VarSelectorDomainOverWeight>(s);
  case Weighter::IPT : return new GenericDVO<VarSelectorImpact>(s);
  default : return NULL;
  }
}


DVORandom::DVORandom(Solver* s) : DVO(s) 
{
}


inline VariableInt* DVORandom::select() 
{
  //return first[(rand() % ( last - first ) )];
  return first[(randint(last - first))];
}


// GenericDVO::~GenericDVO()
// {
//   delete [] _garbage_disjuncts;
// }

// GenericRandomDVO::~GenericRandomDVO() 
// {
//   delete [] bests;
//   delete [] bestvars;
//   delete [] _garbage_disjuncts;
// }
// GenericRandomDVO::GenericRandomDVO(Solver* s, const int sz) : DVO(s) 
//   {
//     size = sz;
//     bests = new T[size+1];
//     bestvars = new VariableInt*[size+1];
//   }

// GenericRandomDVO::~GenericRandomDVO() 
//   {
//     delete bests;
//     delete bestvars;
//   }

/**********************************************
 * Objective Functions 
 **********************************************/

int ObjectiveFunction::update()
{
  int nscore = score();
  if( nscore < upper_bound )
    upper_bound = nscore;
  if(upper_bound == 0) {
    return OPT;
  }
  return UNKNOWN; 
}

MaximiseVar::MaximiseVar(Solver *s, VariableInt *x) 
{
  umore = new UnaryConstraintMore(s, x, x->min());
  maxX = x->max();
  upper_bound = (maxX - x->min());
}

MaximiseVar::~MaximiseVar() 
{
  delete umore;
}

int MaximiseVar::score() 
{
  return (maxX - umore->X->max());
}

int MaximiseVar::solution_score() 
{
  return umore->X->max();
}

int MaximiseVar::update() 
{
  int res = ObjectiveFunction::update();
  if( res == UNKNOWN ) {
    umore->bound = (maxX - upper_bound + 1);
    umore->activate();
  }

  return res; 
}

MinimiseVar::MinimiseVar(Solver *s, VariableInt *x) 
{
  uless = new UnaryConstraintLess(s, x, x->max());
  minX = x->min();
  upper_bound = x->max();
}

MinimiseVar::~MinimiseVar() 
{
  delete uless;
}

int MinimiseVar::score() 
{
  return uless->X->min();
}

int MinimiseVar::solution_score() 
{
  return uless->X->min();
}

int MinimiseVar::update() 
{
  int res = ObjectiveFunction::update();

  //if( res == UNKNOWN ) {
    uless->bound = (upper_bound - 1);
    uless->activate();
    //}
  return res; 
}


