#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <iomanip>

#include "flatzinc.hpp"
#include <set>
#include <map>
using namespace std;
#include <typeinfo>

#ifdef _VERIFICATION
void write_solution(FlatZinc::FlatZincModel *fm, string filename)
{
  if (fm->finished())
    {
      //cout << fm->verif_constraints.size() << filename <<endl;
      unsigned int __size=filename.size();
      for (__size; __size>0; __size--)
	if (filename[__size-1] == '/')
	  break;
      filename.insert(__size,"sol_" );

      ofstream myfile;
      myfile.open (filename.c_str());

      //cout << filename <<endl;
      for (unsigned int i = 0; i < fm->verif_constraints.size() ; i++)
	{
	  if (fm->verif_constraints[i].first != "int_in")
	    {
	      myfile << "constraint " << fm->verif_constraints[i].first << "(" ;
	      int size = fm->verif_constraints[i].second.size();
	      for (unsigned int j = 0; j <size ; j++)
		{
		  myfile << fm->verif_constraints[i].second[j].get_string() ;
		  if (j< (size -1))
		    myfile << " , ";
		}
	      myfile << ");" << endl;
	    }
	}
      myfile <<"solve satisfy;" << endl;
      myfile.close();
      std::cout <<" c DONE" << endl;
    }
}
#endif


int main(int argc, char *argv[])
{
#ifdef _FLATZINC_OUTPUT
  cout << "%";
#endif
  
  SolverCmdLine cmd("Mistral (fzn)", ' ', "2.0");      
  
  TCLAP::SwitchArg annotationArg("","follow_annotations","Uses the annotations", false);
  cmd.add( annotationArg );

  TCLAP::ValueArg<int> parityArg("","parity","Uses parity processing", false, 0, "int");
  cmd.add( parityArg );

  TCLAP::SwitchArg simple_rewriteArg("","simple_rewrite","Uses simple rewriting", false);
  cmd.add( simple_rewriteArg );

  cmd.parse(argc, argv);
  
  usrand(cmd.get_seed());
  
  Solver s;
  
  cmd.set_parameters(s);

  double cpu_time = get_run_time() ;

#ifdef _VERBOSE_PARSER
  std::cout << " " << s.parameters.prefix_comment << " Parse: ";
#endif

  FlatZinc::Printer p;
  FlatZinc::FlatZincModel *fm = 0L;

  fm = parse(cmd.get_filename(), s, p);

  if( !fm ) return 0;
  double parse_time = get_run_time() - cpu_time;
#ifdef _VERBOSE_PARSER
  std::cout << std::endl;
#endif
  cout << " " << s.parameters.prefix_statistics << " PARSETIME " << parse_time << std::endl;




  FlatZinc::SolutionPrinter *sp = new FlatZinc::SolutionPrinter(&p, fm, &s);
  s.add(sp);


  if(s.parameters.time_limit>0) std::cout << " " << s.parameters.prefix_statistics 
					  << " CUTOFF " << s.parameters.time_limit << std::endl;


  // set flatzinc model options

  fm->set_strategy(cmd.get_variable_ordering(), cmd.get_value_ordering(), cmd.get_randomization(), cmd.get_restart_policy());
  fm->set_display_model(cmd.print_model());
  fm->set_display_solution(cmd.print_solution());
  fm->set_annotations(annotationArg.getValue());
  //fm->set_annotations(false);
  fm->set_rewriting(cmd.use_rewrite());
  fm->set_simple_rewriting(simple_rewriteArg.getValue());
  fm->set_parity_processing(parityArg.getValue());
  fm->set_enumeration(cmd.enumerate_solutions());
  fm->encode_clauses();
  fm->run(cout , p);
  
  if(cmd.print_solution())
    fm->print_final(cout , p);

  if(cmd.print_statistics())
    s.statistics.print_full(std::cout);

#ifdef _VERIFICATION
  write_solution(fm, args.back());
#endif

  delete fm;
  delete sp;
  return 0;
}
