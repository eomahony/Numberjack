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

#include "xcsp/XMLParser_libxml2.hh"
#include "xcsp/MistralCallback.hh"
#include "xcsp/CSPParserCallback.hh"

using namespace CSPXMLParser;
using namespace Mistral;



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
      getCommandLine(int_ident, int_param, nia,
		     str_ident, str_param, nsa,
		     &(argv[1]), argc-1);

      if( argc < 2 || int_param[nia-1] != NOVAL || !strcmp( argv[argc-1], "-h" ) )
	{
	  outputHelpMessage();
	}
      else
	{
	  cb.setName         ( strcmp(str_param[0],"nil") ? str_param[0] : "default" );
	  
	  cb.setHeuristic    ( strcmp(str_param[1],"nil") ? str_param[1] : "dom/wldeg" );
	  
	  cb.setRestartPolicy( strcmp(str_param[2],"nil") ? str_param[2] : "dyn" );
	  cb.setRestartFactor( strcmp(str_param[3],"nil") ? atof(str_param[3]) :  (1+1.0/3.0) );
	  cb.setRestartBase  ( int_param[0 ] != NOVAL ? int_param[0 ] : -1 );
	  cb.setRandomize    ( int_param[1 ] != NOVAL ? int_param[1 ] : 1 );
	  cb.setRandomSeed   ( int_param[2 ] != NOVAL ? int_param[2 ] : 11041979 );
	  
	  cb.setTimeLimit( (double)(int_param[3] != NOVAL ? int_param[3] : -1 ) );
	  cb.setNodeLimit    ( int_param[4 ] != NOVAL ? int_param[4] : -1 );

	  cb.setModel        ( int_param[5 ] != NOVAL ? int_param[5 ] :  4 );
	  cb.setInferCliques ( int_param[6 ] != NOVAL ? int_param[6 ] :  2 );
	  cb.setCliquesLimit ( int_param[7 ] != NOVAL ? int_param[7 ] : 10000 );
	  cb.setEVarLimit    ( int_param[8 ] != NOVAL ? int_param[8 ] : 12000 );
	  cb.setEVarSize     ( int_param[9 ] != NOVAL ? int_param[9 ] : 150 );
	  cb.setRevertToSAT  ( int_param[10] != NOVAL ? int_param[10] : 1 );
	  
	  cb.setProbing      ( int_param[11] != NOVAL ? int_param[11] : 0 );
	  cb.setPIteration   ( int_param[12] != NOVAL ? int_param[12] : 100 );
	  cb.setPLimit       ( int_param[13] != NOVAL ? int_param[13] : 30 );
	  cb.setProbeHeuris  ( int_param[14] != NOVAL ? int_param[14] : 0 );
          cb.setPValHeur     ( int_param[15] != NOVAL ? int_param[15] : 0 );
	  cb.setUpdateLW     ( int_param[16] != NOVAL ? int_param[16] : 0 );
	  cb.setUpdateIP     ( int_param[17] != NOVAL ? int_param[17] : 0 );
	  
	  cb.setSAC          ( int_param[18] != NOVAL ? int_param[18] : 0 );
	  cb.setLDS          ( int_param[19] != NOVAL ? int_param[19] : 0 );

	  cb.setDomainSplit  ( int_param[20] != NOVAL ? int_param[20] : 0 );
	  
	  cb.setFeatureExt   ( int_param[21] != NOVAL ? int_param[21] : 0 );
	  cb.setVerbosity    ( int_param[22] != NOVAL ? int_param[22] : 0 );
	  cb.setAllSolution  ( int_param[23] != NOVAL ? int_param[23] : 0 );
          cb.setValHeur      ( int_param[24] != NOVAL ? int_param[24] : 0 );
	  
	  if( int_param[25] != NOVAL )
	    sleep( int_param[25] );

	  cout << "c mistral" << ( cb.use_sac ? "-option" : "-prime") <<  " version 1.550" << endl;
	  
	  parser.parse( argv[1] ); // parse the input file
	  
	  cb.print_outcome( cb.solve() );
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
