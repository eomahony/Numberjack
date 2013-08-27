

#include <mistral_sat.hpp>
#include <stdlib.h>
#include <fstream>
 
using namespace std;


// class Watched {
  
//   Watched( const int nc2 ) {
//     for(int i=0; i<n
//   }

//   // for each pair in (n \choose 2), 
//   int

// }



class Job {

public:
  int initial_rank;

  int duration_on_first_machine;
  int duration_on_second_machine;

  int start_on_first_machine;
  int start_on_second_machine;

  Job(const int d1, const int d2=0) {
    duration_on_first_machine = d1;
    duration_on_second_machine = d2;

    start_on_first_machine = 0;
    start_on_second_machine = 0;

    initial_rank = 0;
  }

};


class Flowshop {

public:
  int target_makespan;
  int **sequ;

  Vector<Job> schedule;

  int get_size() { return schedule.size-1; } 

  Flowshop() { Job j(0); add(j); }
  virtual ~Flowshop() {
    int ntasks = get_size();
    for(int i=0; i<=ntasks; ++i) 
      delete [] sequ[i];
    delete [] sequ;
  }

  int get_makespan() {return schedule.back().start_on_second_machine + schedule.back().duration_on_second_machine; }

  bool is_minimal() {
    bool minimal = true;
    for(int i=get_size(); minimal && --i;) {
      minimal = !swap_safe_probe(i);
    }
    return minimal;
  }

  void set_time(const int i) {

    //std::cout << "set time " << i << std::endl;

    int m1_availability = 0;
    int m2_availability = 0;

    if(i) {
      Job last = schedule[i-1];
    
      m1_availability = last.start_on_first_machine + last.duration_on_first_machine;
      m2_availability = last.start_on_second_machine + last.duration_on_second_machine;
    }

    schedule[i].start_on_first_machine = m1_availability;
    schedule[i].start_on_second_machine = m1_availability+schedule[i].duration_on_first_machine;
    if(schedule[i].start_on_second_machine < m2_availability) schedule[i].start_on_second_machine = m2_availability;
  }

  void close() {

    int ntasks = get_size();
    sequ = new int*[ntasks+1];
    for(int i=0; i<=ntasks; ++i) 
      sequ[i] = new int[ntasks+1];


    int k=0;
    for(int i=1; i<ntasks; ++i) {
      for(int j=i+1; j<=ntasks; ++j) {
	sequ[i][j] = k;
	sequ[j][i] = k;
	
	++k;
      }
    }
  }

  void add(Job j) {
    j.initial_rank = schedule.size;
    schedule.add(j);
    set_time(schedule.size-1);
  }

  void swap(const int x) {
    Job j = schedule[x];
    schedule[x] = schedule[x+1];
    schedule[x+1] = j;
    for(unsigned int i=x; i<schedule.size; ++i) {
      set_time(i);
    }
  }


  bool swap_safe_probe(const int x) {
    bool res = false;
    if(schedule[x].initial_rank < schedule[x+1].initial_rank) {
      swap(x);
      if(get_makespan() <= target_makespan) res = true;
      swap(x);
    }
    return res;
  } 

  bool swap_probe(const int x) {
    if(schedule[x].initial_rank < schedule[x+1].initial_rank) {
      swap(x);
      if(get_makespan() <= target_makespan) return true;
      else swap(x);
    }
    return false;
  } 

  int get_seq_id(const int x, int **seq) {
    return seq[schedule[x].initial_rank][schedule[x+1].initial_rank];
  }

  int get_seq_id(const int x) {
    return sequ[schedule[x].initial_rank][schedule[x+1].initial_rank];
  }


  std::ostream& display(std::ostream& os) const {
    // for(unsigned int i=1; i<schedule.size; ++i) {
    //   os << schedule[i].start_on_first_machine << ":" << schedule[i].duration_on_first_machine << " ";
    // }
    // os << std::endl;
    // for(unsigned int i=1; i<schedule.size; ++i) {
    //   os << schedule[i].start_on_second_machine << ":" << schedule[i].duration_on_second_machine << " ";
    // }
    // os << " (" << target_makespan << ")" << std::endl;

    for(unsigned int i=1; i<schedule.size; ++i) {
      os << schedule[i].initial_rank << " ";
    }
    return os;
  }
  
};

std::ostream& operator<< (std::ostream& os, const Flowshop& x) {
  return x.display(os);
}


void extract_from_file( const char* filename, Flowshop& x, int rho ) {
  int dur;
  std::vector<int> duration;
  double obj;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> obj;
  x.target_makespan = (int)obj;
  
  if(rho>0) {
    x.target_makespan *= (100+rho);
    x.target_makespan /= 100;
  }

  while( true ) {
    infile >> dur;
    if(!infile.good()) break;
    duration.push_back(dur);
  }

  int nJobs = duration.size()/2;

  for(int i=0; i<nJobs; ++i) {
    Job j(duration[i], duration[i+nJobs]);
    x.add(j);
  }

  x.close();

  //exit(1);
}


void explore(Flowshop& F) {

  int ntasks = F.get_size();
  int max_depth = (ntasks * (ntasks-1) / 2);

  BitSet current(0, max_depth, BitSet::empt);

  // //int task_i[max_depth];
  // //int task_j[max_depth];
  // int *sequ[ntasks+1];
  // for(int i=0; i<=ntasks; ++i) 
  //   sequ[i] = new int[ntasks+1];


  // int k=0;
  // for(int i=1; i<ntasks; ++i) {
  //   for(int j=i+1; j<=ntasks; ++j) {

  //     //std::cout << i << "<" << j << " <-> " << k << std::endl;

  //     //task_i[k] = i;
  //     //task_j[k] = j;
      
  //     sequ[i][j] = k;
  //     sequ[j][i] = k;

  //     ++k;
  //   }
  // }

  Vector< BitSet > solutions;
  Vector< BitSet > top_nogoods;

  //std::cout << F << std::endl;

  //std::cout << current << std::endl;

  //Vector< int > local_opt; // 1 if 

  // stack of decisions (a decision is encoded as the rank of the swap)
  Vector< int > decisions;

  // the next decision (swap) to be computed
  int next_swap;

  // data
  int num_nodes = 0;

  // push a decision onto the stack
  decisions.add(0);
  //local_opt.add(1);

  // flag indicating if we stopped because we are in a local optima
  bool local_opt = true;

  while(true) {
    ++num_nodes;
    // if((++num_nodes%1000000)==0) {
    //   std::cout << ".";
    //   std::cout.flush();
    // }
    
    // for(unsigned int i=0; i<decisions.size; ++i)
    //   std::cout << " ";
    // std::cout // << decisions
    //   //<< " " << F  // << " " << local_opt.back()
    //   << " " << current << " " << top_nogoods
    //   << std::endl;

    // if(local_opt.back() != left_branch) {
    //   std::cout << "zarb (1)" << std::endl;
    //   exit(1);
    // }


    // below we look for the next possible move downward in the lattice.
 
   // we skip all the moves that we already tried
    next_swap = decisions.back()+1;
 
   while(true) {

     // we probe the swaps to see if they degrade the solution too much
     while(next_swap<ntasks && !F.swap_probe(next_swap)) {
       //std::cout << " (" << next_swap << ")";
       ++next_swap;
     }
     
      if(next_swap < ntasks) {
	// We found a non-degrading move!

	//current.add(F.get_seq_id(next_swap, sequ));
	current.add(F.get_seq_id(next_swap));
	bool different_all = true;
	unsigned int i=0;
	for(; different_all && i<solutions.size; ++i) {
	  
	  //std::cout << std::endl << current << " inter " << solutions[i] << std::endl;

	  different_all = (current != solutions[i]);
	}
	--i;
	//current.fast_remove(F.get_seq_id(next_swap, sequ));
	current.fast_remove(F.get_seq_id(next_swap));
	if(!different_all) {

	  // However the current yields an old solution

	  F.swap(next_swap);
	  //std::cout << " do not try " << next_swap << " (" << F.get_seq_id(next_swap, sequ) << ") because it is equal to " << solutions[i] << std::endl;
	  //std::cout << " {" << next_swap << "}";	  
	  //local_opt.back() = 0;

	  local_opt = false;

	  ++next_swap;
	} else break;
      } else break;
    }
     
   
   // stopping conditions 
   if(decisions.size > max_depth || next_swap >= ntasks) {

     // we backtrack
     decisions.pop();
     if(decisions.size) {
	
	// if(local_opt.back() != left_branch) {
	//   std::cout << "zarb (2)" << std::endl;
	//   exit(1);
	// }

	//if(local_opt.pop()) { 

	//local_opt.pop();

       // is the current sequence a local optima (if none of the children was one)
       BitSet sol(current);
       if(local_opt) { 
	 //sol.flip();
	 solutions.add(sol);
	 // for(unsigned int i=0; i<decisions.size; ++i)
	 //   std::cout << " ";
	 // std::cout << "NEW SOLUTION" << std::endl;
	 // // //std::cout << F << std::endl;
	 // // //std::cout << sol ;//<< std::endl;
	 // // //std::cout << solutions ;//<< std::endl;
       }


       //top_nogoods.add(sol);
       for(int i=top_nogoods.size-1; i>=0; --i) {
	 if(top_nogoods[i].includes(sol)) {

      
	   //std::cout << "copy " << top_nogoods.back() << " into " << top_nogoods[i] << std::endl;


	
	   if(i<top_nogoods.size-1)
	     top_nogoods[i].copy(top_nogoods.back());
	   --top_nogoods.size;
	 }
       }
       top_nogoods.add(sol);
       


	// continue to backtrack
	F.swap(decisions.back());

	// for(unsigned int i=0; i<decisions.size; ++i)
	//   std::cout << " ";
	// std::cout << current << " - " << F.get_seq_id(decisions.back(), sequ) << std::endl;
	
	//current.fast_remove(F.get_seq_id(decisions.back(), sequ));
	current.fast_remove(F.get_seq_id(decisions.back()));
	//local_opt.back() = 0;

	local_opt = false;

     } else {
       
       // there is no more possible moves at level 0, we exhausted the search tree

       break;
     }
    } else {
      //F.swap(next_swap);



      // for(unsigned int i=0; i<decisions.size; ++i)
      // 	std::cout << " ";
      // std::cout << current // << " swap " << next_swap << "&" << next_swap+1
      // 		<< " + "
      // 		<< F.get_seq_id(next_swap, sequ) << " | " ; //<< std::endl;


      // make the downward move
     //current.add(F.get_seq_id(next_swap, sequ));
     current.add(F.get_seq_id(next_swap));
      decisions.back() = next_swap;
      decisions.add(0);
      //local_opt.add(1);

      //std::cout << F << std::endl;

      local_opt = true;
    }

    //std::cout << std::endl;
  }

  std::cout << std::endl << solutions.size << " solutions " 
	    << num_nodes << " nodes " << std::endl;

}


bool not_included(const int pswap, BitSet& current_lb, Vector< BitSet>& solutions) {
  bool included = false;
  unsigned int i=0;
  for(; !included && i<solutions.size; ++i) {
    included = (solutions[i].contain(pswap) && 
		current_lb.included(solutions[i]));
    if(included)
      std::cout << " -> " << current_lb  << "+" << pswap << " in " << solutions[i] << std::endl;
  }

  
  // if(included) {
  //   std::cout << current_lb << " + " << pswap << " is included in " << solutions[i-1] << std::endl;
  //   exit(1);
  // }

  return !included;
}

int explore2(Flowshop& F) {

  int ntasks = F.get_size();
  int max_depth = (ntasks * (ntasks-1) / 2);

  BitSet current_lb(0, max_depth, BitSet::empt);
  BitSet current_ub(0, max_depth, BitSet::full);

  // stack of decisions (a decision is encoded as the rank of the swap)
  Vector< int > decisions;
  Vector< int > implications;
  Vector< int > num_implications;

  Vector< BitSet > solutions;
  solutions.initialise(2046);

  // the next decision (swap) to be computed
  int pswap, next_swap;

  // data
  int num_nodes = 0;

  // push a decision onto the stack
  decisions.add(0);
  num_implications.add(0);
  //local_opt.add(1);

  // flag indicating if we stopped because we are in a local optima
  bool local_opt = false;

  while(true) {
    ++num_nodes;
    // if((++num_nodes%1000000)==0) {
    //   std::cout << ".";
    //   std::cout.flush();
    // }
    
    // for(unsigned int i=0; i<decisions.size; ++i)
    //   std::cout << " ";
    // std::cout << current_lb << ".." << current_ub << " - " << decisions << " " 
    //   //<< implications << " " << num_implications //<< std::endl;
    // 	      << std::endl;


    
    // below we look for the next possible move downward in the lattice.
    
    // we skip all the moves that we already tried
    next_swap = decisions.back();

    // we probe the swaps to see if they degrade the solution too much
    while(++next_swap<ntasks && 
	  // it is an acceptable move if it corresponds to a swap in the upper bound
	  (!(current_ub.contain((pswap = F.get_seq_id(next_swap))) 
	     && F.swap_probe(next_swap)
	     //&& not_included(pswap,current_lb,solutions)
	     )
	   ) // and it does not degrade the criterion
	  );
  
    
    // stopping conditions 
    if(decisions.size > max_depth || next_swap >= ntasks) {
      
      // we backtrack
      decisions.pop();

      if(decisions.size) {
	
       // is the current sequence a local optima (if none of the children was one)
       if(local_opt) { 
	 //BitSet sol(current_lb);
	 //sol.flip();

	 // test if this is a "new" solution
	 bool not_included = true;
	 for(int i=0; not_included && i<solutions.size; ++i) {
	   not_included = !(current_lb.included(solutions[i]));
	 }

	 if(not_included) {
	   //std::cout << "c " << F << std::endl;

	   solutions.add(current_lb);
	   
	   //   for(unsigned int i=0; i<decisions.size; ++i)
	   //     std::cout << " ";
	   // //std::cout << sol ;//<< std::endl;
	   // //std::cout << solutions ;//<< std::endl;
	 }//  else {
	 //   std::cout << "INCLUDED!" << std::endl;
	 // }
       }
       
       
       // continue to backtrack
       int dec = decisions.back();
       F.swap(dec);
       
       // for(unsigned int i=0; i<decisions.size; ++i)
       // 	 std::cout << " ";
       // // std::cout << current_lb << ".." << current_ub // << " - " << F.get_seq_id(decisions.back()) 
       // // 		 << std::endl;
       // std::cout << " rem "
       // 		 << F.get_seq_id(dec) 
       // 		 << std::endl;
      
       
       pswap = F.get_seq_id(dec);
       current_lb.fast_remove(pswap);
	//local_opt.back() = 0;

       
       //std::cout << "restore ub [" ;

       num_implications.pop();
       int restored = num_implications.back();

       // for(int i=implications.size; i>restored; ) {
       // 	 std::cout << " " << implications[--i] ;
       // }

       // std::cout << " ] (" << restored << ")" << std::endl;

       while(implications.size > restored) {
	 current_ub.fast_add(implications.pop());
       }
       implications.add(pswap);
       ++num_implications.back();
       current_ub.fast_remove(pswap);

       
       local_opt = false;
       
      } else {
       
       // there is no more possible moves at level 0, we exhausted the search tree

       break;
     }
    } else {
      //F.swap(next_swap);

      // for(unsigned int i=0; i<decisions.size; ++i)
      //  	std::cout << " ";
      // std::cout //<< current_lb // 
      // 	// << " swap (" << next_swap 
      // 	// << ":" << F.schedule[next_swap+1].initial_rank << ")"
      // 	// << " & (" << next_swap+1 
      // 	// << ":" << F.schedule[next_swap].initial_rank << ") -" 
      // 	<< " add "
      // 	<< F.get_seq_id(next_swap) 
      // 	<< std::endl;
      
      // make the downward move
      pswap = F.get_seq_id(next_swap);
      current_lb.add(pswap);
      decisions.back() = next_swap;
      decisions.add(0);
      num_implications.add(implications.size);
      //local_opt.add(1);

      local_opt = true;
    }

    //std::cout << std::endl;
  }

  // std::cout << solutions.size << " solutions " 
  // 	    << num_nodes << " nodes " << std::endl;

  std::cout << "\nd NODES " << num_nodes << std::endl;


  // std::cout << implications.size << "/" << implications.capacity << std::endl;
  // std::cout << num_implications.size << "/" << num_implications.capacity << std::endl;
  // std::cout << decisions.size << "/" << decisions.capacity << std::endl;
  // std::cout << solutions.size << "/" << solutions.capacity << std::endl;

  return  solutions.size;
}


bool print_swap(Flowshop& F, const int next_swap) {
  std::cout << "try " << next_swap << ":" << F.get_seq_id(next_swap) << " " << F.swap_safe_probe(next_swap) << std::endl;
  return true;
}

//bool is_local_optimal(Flowshop& F, )


int explore3(Flowshop& F) {

  int ntasks = F.get_size();
  int max_depth = (ntasks * (ntasks-1) / 2);

  BitSet gamma_ub(0, max_depth, BitSet::full);
  BitSet gamma(0, max_depth, BitSet::empt);

  // stack of decisions (index of the swap)
  Vector< int > decisions;

  // stack of decisions (a decision is encoded as the rank of the swap)
  Vector< int > delta_gamma;

  // number of deductions at each level
  Vector< int > size_delta;

  Vector< BitSet > solutions;

  // the next decision (swap) to be computed
  int pswap, next_swap;

  // data
  int num_nodes = 0;

  // flag indicating if we stopped because we are in a local optima
  bool local_opt = false;

  // no deduction at level 0
  size_delta.add(0);

  // 
  decisions.add(0);

  while(true) {
    ++num_nodes;



    // if((++num_nodes%1000000)==0) {
    //   std::cout << ".";
    //   std::cout.flush();
    // }

    // for(unsigned int i=0; i<size_delta.size; ++i)
    //    std::cout << " ";
    //  std::cout // << decisions
    //    //<< " " << F  // << " " << local_opt.back()
    //    << " " << F << " " << gamma << " " << gamma_ub << " " << std::endl;
    // // int j=0;

    // std::cout << "[";
    // for(int k=0; k<size_delta[0]; ++k)
    //   std::cout << " -" << delta_gamma[j++];
    // std::cout << "]";
    // for(int i=1; i<size_delta.size; ++i) {
    //   std::cout << "[" << delta_gamma[j++];
    //   for(int k=0; k<size_delta[i]; ++k)
    // 	std::cout << " -" << delta_gamma[j++];
    //   std::cout << "]";
    // }
    // if(j<delta_gamma.size) {
    //   std::cout << "[" << delta_gamma[j++];
    //   while(j<delta_gamma.size) 
    // 	std::cout << " -" << delta_gamma[j++];
    //   std::cout << "]";
    // }
    // std::cout << " " << size_delta.size << std::endl;

    
    // for(unsigned int i=0; i<size_delta.size; ++i)
    //   std::cout << " ";
    // std::cout << gamma << std::endl;
    // // for(unsigned int j=0; j<solutions.size; ++j) {
    // //   BitSet s(0, max_depth, BitSet::empt);
    // //   s.copy(solutions[j]);
    // //   s.flip();
    // //   //// s.intersect_with(gamma_ub);

    // //   //s.copy(gamma);
      

    // //   for(unsigned int i=0; i<size_delta.size; ++i) 
    // // 	std::cout << " ";
    // //   std::cout << "   " << s << " \\ " << gamma ;//<< std::endl;
    // //   s.setminus_with(gamma);
    // //   std::cout << " = " << s;             


    // //   if(s.size() <= 1) {
    // // 	std::cout << " *** " << std::endl;
    // //   } else if(s.size() ) {
    // // 	std::cout << std::endl;
    // //   } else {
    // // 	std::cout << " *** " << std::endl;
    // // 	exit(1);
    // //   }
    // //   //if(s.size() <= 3) exit(1);
    // // }

    // // for(next_swap=0; ++next_swap<ntasks;)
    // //   std::cout << next_swap << " => " <<  F.get_seq_id(next_swap) << std::endl;
    // // exit(1);

    // //std::cout << F << std::endl; 

    next_swap = decisions.back();
    // we probe the swaps to see if they degrade the solution too much
    while(++next_swap<ntasks && 
	  //print_swap(F, next_swap) &&
	  // it is an acceptable move if it corresponds to a swap in the upper bound
	  (!(gamma_ub.contain((pswap = F.get_seq_id(next_swap))) 
	     
	     //&& not_included(pswap,gamma,solutions)

	     && F.swap_probe(next_swap)
	     
	     )
	   ) // and it does not degrade the criterion
	  );
  
    
    // stopping conditions 
    if(size_delta.size >= max_depth || next_swap >= ntasks) {
      
      if(!size_delta.empty()) {
	
       // is the current sequence a local optima (if none of the children was one)
       if(local_opt) { 

	 // bool not_included = true;
	 // for(int i=solutions.size; not_included && i--;) {
	 //   //for(int i=0; not_included && i<solutions.size; ++i) {
	 //    not_included = !(gamma.included(solutions[i]));
	 //  }
	 //  //not_included = (solutions.empty() || !(gamma.included(solutions.back())));

	 // if(not_included) {
	 //   //std::cout << "c " << F << std::endl;
	 //   solutions.add(gamma);
	 // }

	 if(F.is_minimal()) {
	   solutions.add(gamma);
	 }
       }
       
       
       // undo all the deduction at this level
       int deduction = delta_gamma.size-1;
       int restore = size_delta.pop();
       for(int i=0; i<restore; ++i) {
	 gamma_ub.fast_add(delta_gamma.pop());
       }

       
       if(!delta_gamma.empty()) {
	 // the last element in delta_gamma (before the deductions) 
	 // corresponds to the last decision, it becomes a deduction
	 pswap = delta_gamma.back(); //[deduction-restore];
	 gamma.fast_remove(pswap);
	 gamma_ub.fast_remove(pswap);
	 ++size_delta.back();
	 decisions.pop();
	 F.swap(decisions.back());
       }

       local_opt = false;
       
      } else {
       
       // there is no more possible moves at level 0, we exhausted the search tree

       break;
     }
    } else {

      // make the downward move
      decisions.back() = next_swap;
      pswap = F.get_seq_id(next_swap);
      gamma.add(pswap);
      size_delta.add(0);
      decisions.add(0);
      delta_gamma.add(pswap);

      local_opt = true;
    }

    //std::cout << std::endl;
  }

  // std::cout << solutions.size << " solutions " 
  // 	    << num_nodes << " nodes " << std::endl;

  std::cout << "\nd NODES " << num_nodes << std::endl;


  // std::cout << implications.size << "/" << implications.capacity << std::endl;
  // std::cout << num_implications.size << "/" << num_implications.capacity << std::endl;
  // std::cout << decisions.size << "/" << decisions.capacity << std::endl;
  // std::cout << solutions.size << "/" << solutions.capacity << std::endl;

  return  solutions.size;
}

bool set_next_swap(Flowshop& F, const int next_swap, int& pswap) {
  pswap = F.get_seq_id(next_swap);
  //std::cout << " try " <<  next_swap << "/" << pswap << std::endl;  
  return true;
}


bool print(int i) {
  std::cout << i << std::endl;
  return true;
}

void print_delta(Vector<int>& delta, Vector<int>& size_delta) {
  //delta[size_delta[i]] up delta[size_delta[i]-1] gives all the assumptions at level i  
  for(int i=1; i<size_delta.size; ++i) {
    std::cout << " l" << (i-1) << ":";
    for(int j=size_delta[i-1]; j<size_delta[i]; ++j)
      std::cout << " "<< delta[j] ;
  }
  std::cout << std::endl;
}


int explore4(Flowshop& F) {

  int ntasks = F.get_size();
  int max_depth = (ntasks * (ntasks-1) / 2);

  BitSet gamma_ub(0, max_depth, BitSet::full);
  BitSet gamma_lb(0, max_depth, BitSet::empt);


  // stack of decisions (index of the swap: 1 to n-1)
  Vector< int > decisions;

  // secondary stack of decision, keeps the pair of swapped element (element of gamma)
  Vector< int > delta;
  // number of decision per level, used for backtracking
  Vector< int > size_delta;

  //delta[size_delta[i]] up delta[size_delta[i]-1] gives all the assumptions at level i



  Vector< BitSet > solutions;

  // the next decision (swap) to be computed
  int pswap, next_swap;

  // data
  int num_nodes = 0;

  // flag indicating if we stopped because we are in a local optima
  bool local_opt = true;

  int assumption, last_level;


  // no deduction at level 0
  size_delta.add(0);
  size_delta.add(0);

  decisions.add(0);

  while(true) {
    ++num_nodes;

    for(int i=0; i<decisions.size; ++i) std::cout << " ";
    std::cout << F << " -- "<< gamma_lb << " <= X <= " << gamma_ub << std::endl;   
    // for(int i=0; i<decisions.size; ++i) std::cout << " ";
    // print_delta(delta, size_delta);

    // std::cout << decisions << std::endl;
    // for(int i=0; i<decisions.size; ++i) std::cout << " ";
    // std::cout << delta << std::endl;
    // for(int i=0; i<decisions.size; ++i) std::cout << " ";
    // //int k=0;
    // for(int i=2; i<size_delta.size; ++i) {
    //   for(int k=size_delta[i-1];k<size_delta[i]-1;++k) {
    // 	if(delta[k] < 10) std::cout << "  ";
    // 	else std::cout << "   ";
    //   }
    //   if(delta[size_delta[i]-1] < 10)
    // 	std::cout << " |";
    //   else
    // 	std::cout << " | ";
    // }
    // std::cout << std::endl << size_delta << std::endl;;


    next_swap = decisions.back();
    // we probe the swaps to see if they degrade the solution too much
    while(++next_swap<ntasks &&
	  set_next_swap(F, next_swap, pswap) &&
	  (
	   // it is an acceptable move if it corresponds to a swap in the upper bound
	   !gamma_ub.contain(pswap) ||
	   // and to a swap not in the lower bound
	   gamma_lb.contain(pswap) ||
	   // and it remains optimal
	   !F.swap_probe(next_swap)
	   )
	  );
      
    // left branch
    if(next_swap < ntasks) {
      // for(int i=0; i<decisions.size; ++i) std::cout << " ";
      // std::cout << " branch left (swap " << F.schedule[next_swap].initial_rank << "&" 
      // 		<< F.schedule[next_swap-1].initial_rank << ", i.e., rem " 
      // 		<< pswap << ") "<< std::endl;       

      // 
      decisions.back() = next_swap;
      decisions.add(0);

      //
      //gamma_lb.fast_remove(pswap);
      gamma_ub.fast_remove(pswap);
      delta.add(pswap);

      size_delta.add(delta.size);
      
      local_opt = true;

      // right branch
    } else {
      
      
      if(decisions.size < 2) {
	break;
	// std::cout << "here!!" << std::endl;
	// exit(1);
      }

      //if(size_delta.empty()) break;

      assumption = size_delta.pop();
      last_level = size_delta.back();
      //pswap = delta.back();

      // forget the last decision
      //gamma_ub.fast_add(delta.back());
      //gamma_lb.fast_add(delta.back());
      
      //forget all assumptions from the last level
      while(--assumption>last_level) {
	gamma_lb.fast_remove(delta[assumption]);
      }
      
      // make the deduction
      delta.size = assumption+1;
      gamma_ub.fast_add(delta.back());
      gamma_lb.fast_add(delta.back());
      ++size_delta.back();

      // for(int i=0; i<decisions.size; ++i) std::cout << " ";
      // std::cout << " branch right (add " << delta.back() << ")" << std::endl;

      decisions.pop();


      if(local_opt) {

	bool not_included = true;
	for(int i=solutions.size; not_included && i--;) {
	  //for(int i=0; not_included && i<solutions.size; ++i) {
	  not_included = (solutions[i].included(gamma_ub));
	}
	
	if(not_included) {
	  //std::cout << not_included << " | " << F << std::endl;
	  solutions.add(gamma_ub);
	}
      }
      local_opt = false;
      

      // std::cout << "swap " << decisions.back() << std::endl;
      F.swap(decisions.back());
      
      
    }

  }

  std::cout << "\nd NODES " << num_nodes << std::endl;


  return  solutions.size;
}




int main(int argc, char **argv)
{

  Flowshop f;
  
  // Job j1(10,5);
  // f.add(j1);
  // std::cout << f << std::endl << std::endl;
  
  // Job j2(2,6);
  // f.add(j2);
  // std::cout << f << std::endl << std::endl;
  
  // Job j3(3,3);
  // f.add(j3);
  // std::cout << f << std::endl << std::endl;
  
  // Job j4(7,2);
  // f.add(j4);
  // std::cout << f << std::endl << std::endl;


  // f.swap(1);
  // std::cout << f << std::endl << std::endl;

  // f.swap(2);
  // std::cout << f << std::endl << std::endl;

  // f.swap(3);
  // std::cout << f << std::endl << std::endl;


  //


  double start_stamp = get_run_time();

  const char* filename = argv[1];
  int rho = 0;
  if(argc>2) rho = atoi(argv[2]);

  extract_from_file( filename, f, rho );

  double read_stamp = get_run_time();
  
  //std::cout << f << std::endl;
  
  int num_solutions = explore3(f);

  double solve_stamp = get_run_time();

  std::cout << "d READTIME " 
	    << read_stamp-start_stamp << std::endl;

  std::cout << "d SOLVETIME " 
	    << solve_stamp-read_stamp << std::endl;

  std::cout << "d TOTALTIME " 
	    << solve_stamp-start_stamp << std::endl;

  std::cout << "d SOLUTIONS " 
	    << num_solutions << std::endl;



  std::cout << "\ns SATISFIABLE\nv " << num_solutions << std::endl;

  return 0;
}



