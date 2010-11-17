 

#include "mistral_scd.h"


using namespace MistralScheduler;


int main( int argc, char** argv )
{
  ParameterList params(argc, argv);
  StatisticList stats;

  usrand( params.Seed );

  params.print(std::cout);

  Instance jsp(params);
  
  SchedulingModel *model;
  if(params.Objective == "makespan") {
    std::cout << "c Minimising Makespan" << std::endl;
    if(params.Type == "now") model = new No_wait_Model(jsp, -1);
    else model = new C_max_Model(jsp, -1);
  } else if(params.Objective == "tardiness") {
    std::cout << "c Minimising Tardiness" << std::endl;
    model = new L_sum_Model(jsp, -1);
  } else {
    std::cout << "c unknown objective, exiting" << std::endl;
    exit(1);
  }
  SchedulingSolver solver(model, params, stats);

  stats.print(std::cout, "INIT");  

  if(solver.status == UNKNOWN) solver.dichotomic_search();

  stats.print(std::cout, "DS");

  if(!stats.solved()) solver.branch_and_bound();

  stats.print(std::cout, "");
  
  std::cout << "s SATISFIABLE \nv 00" << std::endl;

  delete model;

}
  




