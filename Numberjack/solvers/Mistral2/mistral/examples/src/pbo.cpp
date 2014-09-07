

//#include <cmath>

#include <mistral_sat.hpp>
#include <stdlib.h>
#include <fstream>
 
using namespace std;



int main(int argc, char **argv)
{


  // double x = 2;
  // double k = -7;

  // double exp = pow(x,k);

  // std::cout << exp << std::endl;

  // exit(1);

  try 
    {
      // Define the command line object.
      SolverCmdLine cmd("Mistral (PBO)", ' ', "2.0");      
      cmd.parse(argc, argv);

  
  usrand(cmd.get_seed());
  
  Solver solver;
  
  cmd.set_parameters(solver);
  
  solver.parse_pbo(cmd.get_filename().c_str());


  
  solver.consolidate();
  
  
   
  BranchingHeuristic *strategy;
  


  if(solver.parameters.backjump) {
    if(cmd.get_value_ordering() == "minval")
      strategy = new GenericHeuristic< VSIDS<2>, MinValue >(&solver);
    else if(cmd.get_value_ordering() == "maxval")
      strategy = new GenericHeuristic< VSIDS<2>, MaxValue >(&solver);
    else if(cmd.get_value_ordering() == "random")
      strategy = new GenericHeuristic< VSIDS<2>, RandomMinMax >(&solver);
    else if(cmd.get_value_ordering() == "minweight")
      strategy = new GenericHeuristic< VSIDS<2>, BoolMinWeightValue >(&solver);
    else  if(cmd.get_value_ordering() == "maxweight")
      strategy = new GenericHeuristic< VSIDS<2>, BoolMaxWeightValue >(&solver);
    else if(cmd.get_value_ordering() == "minval+guided")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< MinValue > >(&solver);
    else if(cmd.get_value_ordering() == "maxval+guided")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< MaxValue > >(&solver);
    else  if(cmd.get_value_ordering() == "random+guided")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< RandomMinMax > >(&solver);
    else if(cmd.get_value_ordering() == "minweight+guided")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< BoolMinWeightValue > >(&solver);
    else if(cmd.get_value_ordering() == "maxweight+guided")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< BoolMaxWeightValue > >(&solver);
    else {
      std::cout << "This value ordering is not handled!\n";
      exit(1);
    }
  } else { 
    strategy = new GenericHeuristic< WDEG<2>, RandomMinMax >(&solver);
  }


  if(cmd.print_model())
    std::cout << solver << std::endl;

  

  RestartPolicy *policy = NULL;

  if(cmd.get_restart_policy() == "no")  
    policy = new NoRestart();
  else if(cmd.get_restart_policy() == "luby")  
    policy = new Luby(solver.parameters.restart_base);
  else if(cmd.get_restart_policy() == "geom")  
    policy = new Geometric(solver.parameters.restart_base,
					  solver.parameters.restart_factor);
  

  VarArray X;
  for(int i=0; i<solver.variables.size; ++i) {
    if(solver.variables[i].is_bool())
      X.add(solver.variables[i]);
  }

  /*
  solver.monitor_list << "objective = ";
  solver.monitor_list << solver.objective->objective;
  solver.monitor_list << " c13 = ";
  solver.monitor_list << solver.constraints[13];
  */

  Outcome result = solver.depth_first_search(X, strategy, policy);
  
  if(cmd.print_statistics())
    cout << solver.statistics ;
  
  if(result == SAT || result == OPT) {
    Solution sol(solver.variables);
    
    if(cmd.print_solution())
      cout << " c  " << solver.variables.size << " " << sol << endl;
  }

#ifdef _CHECK_NOGOOD
  solver.read_solution("sol.txt");
  solver.check_nogoods();
#endif

    } catch (TCLAP::ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
  
  return 0;
}



