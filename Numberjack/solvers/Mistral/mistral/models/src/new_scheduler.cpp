 

#include "mistral_scd.h"


using namespace Mistral;


int main( int argc, char** argv )
{

  ParameterList params(argc, argv);
  usrand(params.Seed);

  StatisticList stats;
  stats.start();

  Instance jsp(params);
  
  if(params.Verbose >= 0) {
    jsp.print(std::cout);
    jsp.printStats(std::cout);
  }

  SchedulingModel *model;
  if(params.Objective == "makespan") {
    if(params.Verbose >= 0) {
      std::cout << "c Minimising Makespan" << std::endl;
    }
    if(params.Type == "now") model = new No_wait_Model(jsp, &params, -1, 0);
    else if(params.Type == "now2") {
      //params.Type = "now";
      model = new No_wait_Model(jsp, &params, -1, 1);
    }
    else model = new C_max_Model(jsp, &params, -1);
  } else if(params.Objective == "tardiness") {
    if(params.Verbose >= 0) {
      std::cout << "c Minimising Tardiness" << std::endl;
    }
    model = new L_sum_Model(jsp, &params, -1);
  } else if(params.Objective == "depth") {
    if(params.Verbose >= 0) {
      std::cout << "c Minimising Depth" << std::endl;
    }    
    model = new Depth_Model(jsp, &params, jsp.getMakespanUpperBound());
  } else if(params.Objective == "weight") {
    if(params.Verbose >= 0) {
      std::cout << "c Minimising Weight" << std::endl;
    }
    model = new DTP_Model(jsp, &params, 1000);
  } else {
    std::cout << "c unknown objective, exiting" << std::endl;
    exit(1);
  }

  SchedulingSolver solver(model, &params, &stats);
  usrand(params.Seed);

  //solver.print(std::cout);
  //exit(1);

  if(params.Verbose >= 0) {
    params.print(std::cout);  
    model->printStats(std::cout);  
    stats.print(std::cout, "INIT");  
  }

  //if
  //solver.jtl_presolve();

  //exit(1);

  if(solver.status == UNKNOWN) solver.dichotomic_search();
  else if( solver.status == SAT ) {
    std::cout << "c Solved while building!" << std::endl;
    exit(1);
    
  } else if( solver.status == UNSAT ) {
    std::cout << "c Found inconsistent while building!" << std::endl;
    exit(1);

  }
      
  if(params.Verbose >= 0) {
    stats.print(std::cout, "DS");
  }

  if(!stats.solved()) {
    if(params.Algorithm == "bnb")
      solver.branch_and_bound();
    else if(params.Algorithm == "lns")
      solver.large_neighborhood_search();
  }

  if(params.Verbose >= 0) {
    stats.print(std::cout, "");  
    std::cout << "s " << (stats.num_solutions ? "SATISFIABLE" : "UNSATISFIABLE") 
	      << " \nv 00" << std::endl;
  }

  int i = params.PrintSolution;
  //if(i > solver.pool->size()) i = 
  int j = solver.pool->size();
  while(j-- && i--) {

    //std::cout << i << " " << j << std::endl;

    (*(solver.pool))[j]->print(std::cout);
  }
  delete model;

}
  




