/*=============================================================================
 * Copyright (c) 2008 Olivier ROUSSEL (olivier.roussel <at> cril.univ-artois.fr)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *=============================================================================
 */

#include <iostream>

 #include "xml/XMLParser_libxml2.hh"
 #include "xml/MistralCallback.hh"
 #include "xml/CSPParserCallback.hh"
//#include <mistral_solver.hpp>

using namespace CSPXMLParser;
using namespace Mistral;
using namespace std;


const int nia = 27;
const char* int_ident[nia] = {
  "-base", 
  "-randomized", 
  "-seed",

  "-timelimit",
  "-nodelimit",
  
  "-model", 
  "-clique",
  "-maxclique",
  "-branch_extra_vars", 
  "-max_extra_size",
  "-sat",

  "-probe", 
  "-probe_iterations", 
  "-probe_limit", 
  "-probe_heuristic",
  "-probe_valheuristic",
  "-update_lweight", 
  "-update_impact",  

  "-sac", 
  "-lds", 

  "-dom_splitting",

  "-feature_extraction",
  "-verbose",
  "-all",
  "-val_h",
  "-sleep",
  "-h" 
};

const char* int_man[nia] = {
  "set restart base (default: dynamic)", 
  "randomize the list of variables before each restart / randomize the heuristic \n                       (0: no randomization, >0: randomize the sequence, k: randomize the heuristic (|k| best choices) | default: 1)", 
  "set random seed (default: 11041979)\n",

  "set time limit (default: no limit)",
  "set the node limit (default: no limit)\n", 

  "set model type (0: extra variables, 1: no extra variables, 3: dynamic | default: dynamic)", 
  "infers implied AllDifferent constraints (0: no, 1: some, 2: all | default: some)", 
  "set a limit on the number of infered implied AllDifferent constraints (default: 10000)", 
  "set extra variable tolerance (use extra vars for branching if there are less than this number | default: 12000)", 
  "set extra variable tolerance (use only exta vars whose domain is smaller than this number | default: 150)\n", 
  "switch to SAT solver when possible (0: no, 1: yes | default yes)",

  "set probing (0: no, 1: yes | default: no)", 
  "set number of probes (default: 100)", 
  "set probe limit (default: 30)", 
  "set heuristic used during probing (0: random, 1: wdeg, 3: wldeg, 8: impact | default: random)", 
  "set val heuristic used during probing (0: lex, 1: reverse lex, 2: random | default: lex)",
  "update level-based weights when probing (0: no, 1: yes | default: no)", 
  "update impacts weights when probing (0: no, 1: yes | default: no)\n", 

  "singleton arc consistency preprocessing (0: no, 1: imcomplete, 2: complete 3: dynamic | default: no)\n", 
  "limited dicrepancy search (amount of time used to initialised the oredering default: no)\n", 

  "use domain splitting on large domains (0: no, 1: yes, 2: dynamic | default: no)",

  "extract features (0: no, 1: yes | default: no)",
  "set verbosity (default: quiet)", 
  "finds all solutions",
  "sleep for x seconds before solving",
  "set val heuristic (0: lex, 1: reverse lex, 2: random | default: lex)",
  "display this help message" 
   };

			      
int int_param[nia];

const int nsa = 4;
const char* str_ident[nsa] = {
  "-name",

  "-heuristic", 

  "-restart",
  "-factor" 
};
const char* str_man[nsa] = {
  "set output file extension (default: no output in file)\n",

  "set heuristic (dom, dom/deg, dom/wdeg, dom/wldeg, impact, impact/wdeg, etc.. default: dom/wdeg)\n",
 
  "set restart policy (no, geom, luby, dyn | default: dyn)", 
  "set geometric restart factor (default: 1.333)" 
};
const char* str_param[nsa];


void outputHelpMessage()
{
  cout << "\nUsage: \t bin/solver [path to the xml file] <args>" << endl << endl;
  for(int i=0; i<nsa; ++i)
    cout  << setw(20) << str_ident[i] << ":  " << str_man[i] << endl;
  for(int i=0; i<nia; ++i)
    cout  << setw(20) << int_ident[i] << ":  " << int_man[i] << endl;
  cout << endl;
}


int main(int argc, char **argv)
{
  MistralCallback cb; // my interface between the parser and the solver
  

  try
    {
      // the parser that will decode the XML file
      XMLParser_libxml2<MistralCallback> parser(cb); 
      parser.setPreferredExpressionRepresentation(TREE);

      str_param[0] = "";
      // getCommandLine(int_ident, int_param, nia,
      // 		     str_ident, str_param, nsa,
      // 		     &(argv[1]), argc-1);

      if( argc < 2 || !strcmp( argv[argc-1], "-h" ) )
  	{
  	  outputHelpMessage();
  	}
      else
  	{
  	  cout << "c mistral" <<  " version 2.0" << endl;
	  
  	  parser.parse( argv[1] ); // parse the input file
	  
	  Outcome result = cb.solve();

	  cout << outcome2str(result) << endl;

	  cout << cb.solver.statistics << endl;

  	  //cb.print_outcome( cb.solve() );
  	  //cb.cleanup();
  	}

    }
  catch (exception &e)
    {
      cout.flush();
      cerr << "\n\tUnexpected exception :\n";
      cerr << "\t" << e.what() << endl;
      exit(1);
    }

  return 0;
}
