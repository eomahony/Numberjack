 

#include "mistral_scheduler.hpp"


using namespace Mistral;


int main( int argc, char** argv )
{

  ParameterList params(argc, argv);
  usrand(params.Seed);

  StatisticList stats;
  stats.start();

  Instance jsp(params);
  
  std::cout << std::endl;
  
	jsp.print(std::cout);
  
	jsp.printStats(std::cout);
  params.print(std::cout);


  SchedulingSolver *solver;
  if(params.Objective == "makespan") {
    std::cout << "c Minimising Makespan" << std::endl;
    if(params.Type == "now") solver = new No_wait_Model(jsp, &params, -1, 0);
    else if(params.Type == "now2") {
      //params.Type = "now";
      solver = new No_wait_Model(jsp, &params, -1, 1);
    }
    else solver = new C_max_Model(&jsp, &params, &stats);
  } else if(params.Objective == "tardiness") {
    std::cout << "c Minimising Tardiness" << std::endl;
    solver = new L_sum_Model(&jsp, &params, &stats);
  } // else if(params.Objective == "depth") {
  //   std::cout << "c Minimising Depth" << std::endl;
  //   solver = new Depth_Model(jsp, &params, jsp.getMakespanUpperBound());
  // } else if(params.Objective == "weight") {
  //   std::cout << "c Minimising Weight" << std::endl;
  //   solver = new DTP_Model(jsp, &params, 1000);
  // }
  else {
    std::cout << "c unknown objective, exiting" << std::endl;
    exit(1);
  }
	
	solver->setup();
	

  // SchedulingSolver solver(model, &params, &stats);
  usrand(params.Seed);
  
  solver->consolidate();

  //std::cout << solver << std::endl;
  
  solver->dichotomic_search();
  

  if(!stats.solved()) {
    if(params.Algorithm == "bnb")
      solver->branch_and_bound();
    //else if(params.Algorithm == "lns")
    //solver.large_neighborhood_search();
  }

  stats.print(std::cout, "");  
  std::cout << "s " << (stats.num_solutions ? "SATISFIABLE" : "UNSATISFIABLE") 
   	    << " \nv 00" << std::endl;


#ifdef _PROFILING
  std::cout << solver->statistics.total_propag_time << std::endl;
#endif


}
  




