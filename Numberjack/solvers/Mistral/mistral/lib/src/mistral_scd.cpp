 

#include "mistral_scd.h"
#include <math.h>

using namespace Mistral;


#define INFTY 0xffffff
#define BIG 0xffff

#define DICHO 0
#define BNB   1
#define LNS   2


StatisticList::StatisticList() {
  best_solution_index     = 0;
  branch_and_bound_index  = 0;

  lower_bound             = INFTY;
  upper_bound             = 0;

  num_nogoods             = 0;
  avg_nogood_size         = 0;
  num_solutions           = 0;
  
  avg_distance            = 0;
  min_distance            = INFTY;
  max_distance            = 0;

  normalized_objective    = 1.0;

  real_start_time         = 0.0;
}

StatisticList::~StatisticList() {}

void StatisticList::start() {
  real_start_time = getRunTime();
}

// void StatisticList::stop() {
//   real_time = (getRunTime() - real_time);
// }

bool StatisticList::solved() {
  return (lower_bound == upper_bound);
}

double StatisticList::get_total_time() {
  double total_time  = 0.0;
  for(unsigned int i=0; i<time.size(); ++i)
    total_time += time[i];
  return total_time;
}

double StatisticList::get_lowerbound_time() {
  double total_time  = 0.0;
  for(unsigned int i=0; i<time.size(); ++i)
    if(outcome[i] == UNSAT)
      total_time += time[i];
  return total_time;
}

void StatisticList::add_info(const int objective, int tp) {

  //bool update_ub = true;
  DBG("Update statistics%s\n", "");


  time.push_back(solver->ENDTIME);
  soltime.push_back(solver->SOLTIME);
  nodes.push_back(solver->NODES);
  backtracks.push_back(solver->BACKTRACKS);
  fails.push_back(solver->FAILURES);
  propags.push_back(solver->PROPAGS);
  types.push_back(tp);

  outcome.push_back(solver->status);

  if(outcome.back() == SAT || outcome.back() == OPT) {
    ++num_solutions;
    best_solution_index = outcome.size()-1;
    upper_bound = objective;
    if(outcome.back() == OPT) lower_bound = objective;
  } else if(types.back() != LNS && outcome.back() == UNSAT) {
    lower_bound = objective+1;
  }
}



std::ostream& StatisticList::print(std::ostream& os, 
				   const char* prefix,
				   const int start, 
				   const int end) {

  int k, i=start, j=outcome.size();
  if(end >= 0) j=end;

  double total_time      = 0.0;
  double opt_time        = 0.0;
  double ub_time         = 0.0;
  double lb_time         = 0.0;
  double lost_time       = 0.0;
  double proof_time      = 0.0;
  double avg_cutoff_time = 0.0;
  long unsigned int total_nodes         = 0;
  long unsigned int total_backtracks    = 0;
  long unsigned int total_fails         = 0;
  long unsigned int total_propags       = 0;
  
  int nb_unknown = 0;
  for(k=i; k<j; ++k) {
//     if(k<=best_solution_index) {
//       if(outcome[k] != OPT)
// 	opt_time += time[k];
//       else
// 	opt_time += soltime[k];
//     }
    if(k==best_solution_index) opt_time += soltime[k];
    else opt_time += time[k];
    
    if(outcome[k] != OPT && outcome[k] != UNSAT) {
      lost_time += (time[k] - soltime[k]);
    }
    
    ub_time            += soltime[k];
    total_time         += time[k];
    total_nodes        += nodes[k];
    total_backtracks   += backtracks[k];
    total_fails        += fails[k];
    total_propags      += propags[k];

    if(types[k]==DICHO && outcome[k] == UNKNOWN)
      {
	++nb_unknown;
	avg_cutoff_time += time[k];
      }

  }

  if(nb_unknown)
    avg_cutoff_time /= (double)nb_unknown;
  proof_time = (total_time - opt_time);
  lb_time = (total_time - ub_time - lost_time);

  if(lb_time < 0.00001) lb_time = 0.0;
  if(lost_time < 0.00001) lost_time = 0.0;
  if(proof_time < 0.00001) proof_time = 0.0;


  avg_distance = 0;
  min_distance = NOVAL;
  max_distance = 0;
  for(i=1; (unsigned int)i<solver->pool->size(); ++i) {
    double dist = (*(solver->pool))[i-1]->distance((*(solver->pool))[i]);
    avg_distance += dist;
    if(dist>max_distance) max_distance = dist;
    if(dist<min_distance) min_distance = dist;
  }
  avg_distance /= solver->pool->size();


  int plength = 0;
  while(prefix[plength] != '\0') ++plength;
  
  os << "c =================[ statistics ]==================" << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "LOWERBOUND "    << std::right << std::setw(21) << lower_bound << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "UPPERBOUND "    << std::right << std::setw(21) << upper_bound << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "OBJECTIVE "     << std::right << std::setw(21) << upper_bound << std::endl 
    << "d " << prefix << std::left << std::setw(28-plength)  << "NORMOBJECTIVE " << std::right << std::setw(21) << normalized_objective << std::endl 
    << "d " << prefix << std::left << std::setw(28-plength)  << "REALTIME "      << std::right << std::setw(21) << (getRunTime() - real_start_time) << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "RUNTIME "       << std::right << std::setw(21) << total_time << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "OPTTIME "       << std::right << std::setw(21) << opt_time << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "PROOFTIME "     << std::right << std::setw(21) << proof_time << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "UBTIME "        << std::right << std::setw(21) << ub_time << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "LBTIME "        << std::right << std::setw(21) << lb_time << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "LOSTTIME "      << std::right << std::setw(21) << lost_time << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "NODES "         << std::right << std::setw(21) << total_nodes << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "BACKTRACKS "    << std::right << std::setw(21) << total_backtracks << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "FAILS "         << std::right << std::setw(21) << total_fails << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "PROPAGS "       << std::right << std::setw(21) << total_propags << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "NOGOODS "       << std::right << std::setw(21) << num_nogoods << std::endl;
  if(num_nogoods)
    std::cout << "d " << prefix << std::left << std::setw(28-plength)  << "NOGOODSIZE "    << std::right << std::setw(21) << avg_nogood_size/(double)num_nogoods << std::endl;
  std::cout << "d " << prefix << std::left << std::setw(28-plength)  << "NODES/s "       << std::right << std::setw(21) ;
  if(total_time > 0)
    std::cout << (int)((double)total_nodes/total_time);
  else 
    std::cout << "N/A";
  std::cout << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "BACKTRACKS/s "  << std::right << std::setw(21) ;
  if(total_time > 0)
    std::cout << (int)((double)total_backtracks/total_time);
  else 
    std::cout << "N/A";
  std::cout << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "FAILS/s "       << std::right << std::setw(21) ;
  if(total_time > 0)
    std::cout << (int)((double)total_fails/total_time);
  else 
    std::cout << "N/A";
  std::cout << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "PROPAGS/s "     << std::right << std::setw(21) ;
  if(total_time > 0)
    std::cout << (int)((double)total_propags/total_time);
  else 
    std::cout << "N/A";
  std::cout 
    << std::endl
    //<< "d " << prefix << std::left << std::setw(28-plength)  << "RESTARTS "      << std::right << std::setw(21) << total_restarts << std::endl
    
    << "d " << prefix << std::left << std::setw(28-plength)  << "SOLUTIONS "     << std::right << std::setw(21) << num_solutions << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "AVGCUTOFF "     << std::right << std::setw(21) << avg_cutoff_time << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "AVGDISTANCE "   << std::right << std::setw(21) << avg_distance << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "MINDISTANCE "   << std::right << std::setw(21) << min_distance << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "MAXDISTANCE "   << std::right << std::setw(21) << max_distance << std::endl
    << "d " << prefix << std::left << std::setw(28-plength)  << "OPTIMAL "       << std::right << std::setw(21) << (lower_bound == upper_bound) << std::endl
    << "c =================[ statistics ]==================" << std::endl
    << std::endl;

  return os;
}  


const char* ParameterList::int_ident[ParameterList::nia] = 
  {"-ub", "-lb", "-check", "-seed", "-cutoff", "-dichotomy", 
   "-base", "-randomized", "-verbose", "-optimise", "-nogood", 
   "-dyncutoff", "-nodes", "-hlimit", "-init", "-neighbor", 
   "-initstep", "-fixtasks", "-order", "-ngdt", "-print_sol"};

const char* ParameterList::str_ident[ParameterList::nsa] = 
  {"-heuristic", "-restart", "-factor", "-decay", "-type", 
   "-value", "-dvalue", "-ivalue", "-skew", "-objective", 
   "-algo", "-presolve"};


// ParameterList::ParameterList() {
// }

// ParameterList::ParameterList(const ParameterList& pl) {
//   initialise(pl);
// }

ParameterList::ParameterList(int length, char **commandline) {

  if( length < 2 )
    {
      std::cerr << "need a data file" << std::endl;
      exit( 0 );
    }

  
  data_file = commandline[1];

  int i=0;
  while(commandline[1][i] != '\0') ++i;
  while(commandline[1][i] != '/') --i;
  data_file_name = &(commandline[1][i+1]);

  getCommandLine(ParameterList::int_ident,
		 int_param,
		 ParameterList::nia,
		 ParameterList::str_ident,
		 str_param,
		 ParameterList::nsa,
		 &commandline[1],length-1);

  if(strcmp(str_param[4],"nil")) Type = str_param[4];
  else {
    Type = "jsp";
    std::cout << "c Warning: no type specified, treating the data as Taillard's jsp" << std::endl;
  }

  UBinit      = -1;
  LBinit      = -1;
  Checked     = true;
  Seed        = 12345;
  Cutoff      = 300;
  NodeCutoff  = 0;
  NodeBase    = 30;
  Dichotomy   = 128;
  Base        = 256;
  Randomized  = -1;
  Precision   = 100;
  Hlimit      = 20000;
  InitBound   = 1000;
  InitStep    = 1;
  Neighbor    = 2;

  Verbose     = 0;
  Optimise    = 3600;
  Rngd        = 2;

  Policy    = "geom";
  Factor    = 1.3;
  Decay     = 0.0;
  Value     = "guided";
  DValue    = "guided";
  IValue    = "promise";
  Skew      = -1.0;
  Objective = "makespan";
  Algorithm = "bnb";
  Presolve  = "default";

  Heuristic = "none";
  PolicyRestart = GEOMETRIC;
  FixTasks  = 0;
  NgdType   = 2;
  OrderTasks = 1;

  PrintSolution = 0;


  if(Type == "osp") {
    Objective = "makespan";
    if(Heuristic == "none")
      Heuristic = "osp-b";
  } else if(Type == "sds") {
    Objective = "makespan";
    if(Heuristic == "none")
      Heuristic = "osp-b";
  } else if(Type == "jtl") {
    Objective = "makespan";
    Presolve = "jtl";
    if(Heuristic == "none")
      Heuristic = "osp-t";
  } else if(Type == "now" || Type == "now2") {
    Objective = "makespan";
    Presolve = "default";
    if(Heuristic == "none")
      Heuristic = "osp-dw";
  } else if(Type == "jla") {
    Objective = "makespan";
    if(Heuristic == "none")
      Heuristic = "osp-t";
  } else if(Type == "jsp") {
    Objective = "makespan";
    if(Heuristic == "none")
      Heuristic = "osp-t";
  } else if(Type == "fsp") {
    Verbose = -1;
    OrderTasks = -1;
    Objective = "makespan";
    if(Heuristic == "none")
      Heuristic = "osp-t";
  } else if(Type == "jet") {
    Objective = "tardiness";
    if(Heuristic == "none")
      Heuristic = "osp-t";
  } else if(Type == "dyn") {
    Objective = "tardiness";
    if(Heuristic == "none")
      Heuristic = "osp-t";
  }

  if(int_param[0]  != NOVAL) UBinit      = int_param[0];
  if(int_param[1]  != NOVAL) LBinit      = int_param[1];
  if(int_param[2]  != NOVAL) Checked     = int_param[2];
  if(int_param[3]  != NOVAL) Seed        = int_param[3];
  if(int_param[4]  != NOVAL) Cutoff      = int_param[4];
  if(int_param[5]  != NOVAL) Dichotomy   = int_param[5];
  if(int_param[6]  != NOVAL) Base        = int_param[6];
  if(int_param[7]  != NOVAL) Randomized  = int_param[7]; 
  if(int_param[8]  != NOVAL) Verbose     = int_param[8];
  if(int_param[9]  != NOVAL) Optimise    = int_param[9]; 
  if(int_param[10] != NOVAL) Rngd        = int_param[10];
  if(int_param[13] != NOVAL) Hlimit      = int_param[13]; 
  if(int_param[14] != NOVAL) InitBound   = int_param[14]; 
  if(int_param[15] != NOVAL) Neighbor    = int_param[15]; 
  if(int_param[16] != NOVAL) InitStep    = int_param[16]; 
  if(int_param[17] != NOVAL) FixTasks    = int_param[17]; 
  if(int_param[18] != NOVAL) OrderTasks  = int_param[18]; 
  if(int_param[19] != NOVAL) NgdType     = int_param[19]; 
  if(int_param[20] != NOVAL) PrintSolution = int_param[20]; 

  if(strcmp(str_param[0 ],"nil")) Heuristic  = str_param[0];
  if(strcmp(str_param[1 ],"nil")) Policy     = str_param[1];
  if(strcmp(str_param[2 ],"nil")) Factor     = atof(str_param[2]);
  if(strcmp(str_param[3 ],"nil")) Decay      = atof(str_param[3]);
  if(strcmp(str_param[5 ],"nil")) Value      = str_param[5];
  if(strcmp(str_param[6 ],"nil")) DValue     = str_param[6];
  if(strcmp(str_param[7 ],"nil")) IValue     = str_param[7];
  if(strcmp(str_param[8 ],"nil")) Skew       = atof(str_param[8]);
  if(strcmp(str_param[9 ],"nil")) Objective  = str_param[9];
  if(strcmp(str_param[10],"nil")) Algorithm  = str_param[10];
  if(strcmp(str_param[11],"nil")) Presolve   = str_param[11];

}

void ParameterList::initialise(SchedulingSolver *s) {

  solver = s;
  if(int_param[11] != NOVAL) NodeBase    = int_param[11];

  // since the model is different for no-wait, we use a different factor
  unsigned long long int NodeFactor = 20000000;
  NodeCutoff = (NodeFactor * NodeBase);
  if(int_param[12] != NOVAL) NodeCutoff  = int_param[12]; 

  if(Policy == "luby")
    PolicyRestart = LUBY;
  else if(Policy == "geom")
    PolicyRestart = GEOMETRIC;
  else
    PolicyRestart = NO;
}

double ParameterList::getSkew() {

  if(Skew < 0) Skew = ((double)(solver->stats->upper_bound)/
		       ((double)solver->stats->lower_bound+1));
  if(Skew > 100.0) Skew = 100.0;
  if(Skew != 1.0) {
    if(Skew < 1.0) Skew = 1.0; 
    double rs = 2 + randreal()/10; //(randreal() - .5)/2.0;
    //Skew += rs;
    Skew *= rs;
  }

  return Skew;

}

std::ostream& ParameterList::print(std::ostream& os) {
  os << "c =================[ parameters ]==================" << std::endl;
  os << std::left << std::setw(30) << "c data file " << ":" << std::right << std::setw(20) << data_file_name << std::endl;
  os << std::left << std::setw(30) << "c type " << ":" << std::right << std::setw(20) << Type << std::endl;
  os << std::left << std::setw(30) << "c seed " << ":" << std::right << std::setw(20) << Seed << std::endl;
  os << std::left << std::setw(30) << "c greedy iterations " << ":" << std::right << std::setw(20) << InitBound << std::endl;
  os << std::left << std::setw(30) << "c skew " << ":" << std::right << std::setw(20) << // getSkew()
    Skew
     << std::endl;
  os << std::left << std::setw(30) << "c use initial probe " << ":" << std::right << std::setw(20) << (InitStep ? "yes" : "no") << std::endl;
  os << std::left << std::setw(30) << "c time cutoff " << ":" << std::right << std::setw(20) << Cutoff << std::endl;
  os << std::left << std::setw(30) << "c node cutoff " << ":" << std::right << std::setw(20) << NodeCutoff << std::endl;
  os << std::left << std::setw(30) << "c dichotomy " << ":" << std::right << std::setw(20) << (Dichotomy ? "yes" : "no") << std::endl;
  os << std::left << std::setw(30) << "c restart policy " << ":" << std::right << std::setw(20) << Policy << std::endl;
  os << std::left << std::setw(30) << "c base " << ":" << std::right << std::setw(20) << Base << std::endl;
  os << std::left << std::setw(30) << "c factor " << ":" << std::right << std::setw(20) << Factor << std::endl;
  os << std::left << std::setw(30) << "c heuristic " << ":" << std::right << std::setw(20) << Heuristic << " (" << abs(Randomized) << ")" << std::endl;
  os << std::left << std::setw(30) << "c value ordering (init step) " << ":" << std::right << std::setw(20) << IValue << std::endl;
  os << std::left << std::setw(30) << "c value ordering (dichotomy) " << ":" << std::right << std::setw(20) << DValue << std::endl;
  os << std::left << std::setw(30) << "c value ordering (optim) " << ":" << std::right << std::setw(20) << Value << std::endl;
  os << std::left << std::setw(30) << "c randomization " << ":" ;
     
  if(Randomized<-1) 
    os << std::right
       << std::setw(17) << "i-shuff random (" << std::setw(2) << -Randomized << ")" << std::endl;
  else if(Randomized>1) 
    os << std::right
       << std::setw(17) << "shuff & random (" << std::setw(2) << Randomized << ")" << std::endl;
  else if(Randomized==1) 
    os << std::right
       << std::setw(20) << "shuffled" << std::endl;
  else if(Randomized==0) 
    os << std::right
       << std::setw(20) << "not random" << std::endl;
  else 
    os << std::right
       << std::setw(20) << "init shuffle" << std::endl;
  
  os << "c =================[ parameters ]==================" << std::endl;
  
  return os;
}


Instance::Instance(ParameterList& params) {

  DBG("Build instance %s\n", params.data_file);

  dtp_nodes = 0;
  
  setup_time    = NULL;
  time_lag[0]   = NULL;
  time_lag[1]   = NULL;
  jsp_duedate   = NULL;
  jsp_latecost  = NULL;
  jsp_earlycost = NULL;
  jsp_floatcost = NULL;

  max_makespan  = INFTY;
  
  if(params.Type == "osp") {
    osp_readData( params.data_file );
  } else if(params.Type == "sds") {
    sds_readData( params.data_file );
  } else if(params.Type == "jtl") {
    jtl_readData( params.data_file );
  } else if(params.Type == "now" || params.Type == "now2") {
    now_readData( params.data_file );
  } else if(params.Type == "jla") {
    jla_readData( params.data_file );
  } else if(params.Type == "tsp") {
    tsp_readData( params.data_file );
  } else if(params.Type == "fsp") {
    fsp_readData( params.data_file );
  } else if(params.Type == "pfsp") {
    fsp_readData( params.data_file );
  } else if(params.Type == "jsp") {
    jsp_readData( params.data_file );
  } else if(params.Type == "jet") {
    jet_readData( params.data_file );
  } else if(params.Type == "dyn") {
    dyn_readData( params.data_file, params.Precision );
  } else if(params.Type == "dtp") {
    dtp_readData( params.data_file );
  } 
}

Instance::~Instance() {
  int i, j;

  if(hasSetupTime()) {
    for(i=0; i<nMachines(); ++i) {
      for(j=0; j<nJobs(); ++j) {
	delete [] setup_time[i][j];
      }
      delete [] setup_time[i];
    }
    delete [] setup_time;
  }

  if(hasTimeLag()) {
    for(i=0; i<nJobs(); ++i) {
      delete [] time_lag[0][i];
      delete [] time_lag[1][i];
    }
    delete [] time_lag[0];
    delete [] time_lag[1];
  }

  if(hasJobDueDate()) {
    delete [] jsp_duedate;
  }

  if(hasLateCost()) {
    delete [] jsp_latecost;
  }

  if(hasEarlyCost()) {
    delete [] jsp_earlycost;
  }
}


int straight_compar(const void *x, const void *y) {
  int a = *(int*)x;
  int b = *(int*)y;
  return(a==b ? 0 : (a<b ? -1 : 1));
}

void Instance::get_placement(const int x, const int y, Vector<int>& intervals) 
{  
  // we compute the allowed/forbidden start times of y relative to x
  // (like 'y in [x+30, x+70])

  // stores the start of each forbidden interval
  Vector<int> forbidden_starts;
  // stores the end of each forbidden interval
  Vector<int> forbidden_ends;

  // for each machine M the starting time of the task processed on M 
  // for job x
  int *timetable_x = new int[nMachines()];
  int *timetable_y = new int[nMachines()];

  // for each machine M the duration of the task processed on M
  // for job y
  int *duration_x = new int[nMachines()];
  int *duration_y = new int[nMachines()];


  int dur = 0;
  for(int k=0; k<nTasksInJob(x); ++k) {
    int t = getJobTask(x,k);
    for(int i=0; i<nMachines(t); ++i) {
      int m = getMachine(t,i);
      duration_x[m] = getDuration(t);
      timetable_x[m] = dur;
    }
    dur += getDuration(t);
  }
  dur = 0;
  for(int k=0; k<nTasksInJob(y); ++k) {
    int t = getJobTask(y,k);
    for(int i=0; i<nMachines(t); ++i) {
      int m = getMachine(t,i);
      duration_y[m] = getDuration(t);
      timetable_y[m] = dur;
    }
    dur += getDuration(t);
  }



  // compute all the forbidden intervals
  for(int k=0; k<nMachines(); ++k) {
    // forbidden interval for machine k
    // y cannot start in [st_x-duration_y+1, st_x+duration_x-1]
    // (since x is using machine k at this time)
    forbidden_starts.push(timetable_x[k]-timetable_y[k]-duration_y[k]+1);
    forbidden_ends.push(timetable_x[k]-timetable_y[k]+duration_x[k]);
  }


  // Now the cool part, we want to compute the 'real' intervals, since
  // some that we just computed can be merged if they overlap

  // we sort the intervals ends
  qsort(forbidden_starts.stack_, forbidden_starts.size, sizeof(int), straight_compar);
  qsort(forbidden_ends.stack_, forbidden_ends.size, sizeof(int), straight_compar);

  int i=0, j=0;
  int current = 0;    
  
  // now we can go forward from the earliest interval to latest
  while(i<forbidden_starts.size && j<forbidden_ends.size) {
    if(forbidden_starts[i]<forbidden_ends[j]) {
      ++i;

      // when we see the start of an interval, the current number
      // of machine that forbids this start time increases
      if(current++==0) {
	// if it goes from 0 to 1, it is the start of a forbidden interval
	intervals.push(forbidden_starts[i-1]-1);
      }
    } else if(forbidden_starts[i]>forbidden_ends[j]) {
      ++j;

      // when we see the end of an interval, the current number
      // of machine that forbids this start time decreases
      if(--current==0) {
	// if it goes from 1 to 0, it is the end of a forbidden interval
	intervals.push(forbidden_ends[j-1]);
      }
    } else {
      ++i;
      ++j;
    }
  }

  intervals.push(forbidden_ends.back());

//   int num__intervals = intervals.size/2+1;
//   int *max__intervals = new int[num__intervals];
//   int *min__intervals = new int[num__intervals];
//   min__intervals[0] = -NOVAL;
//   max__intervals[num__intervals-1] = NOVAL;
//   int k=0;
//   for(int i=0; i<num__intervals-1; ++i) {
//     max__intervals[i] = intervals[k++];
//     min__intervals[i+1] = intervals[k++];
//   }

}

void Instance::get_placement2(const int x, const int y) {
  // compute x's and y's lengths.
  int x_length=0, y_length=0;

  for(int k=0; k<nTasksInJob(x); ++k) {
    x_length += getDuration(getJobTask(x,k));
  }

  for(int k=0; k<nTasksInJob(y); ++k) {
    y_length += getDuration(getJobTask(y,k));
  }



  int *timetable_x_start = new int[nMachines()];
  int *timetable_x_end = new int[nMachines()];

  int *timetable_y_start = new int[nMachines()];
  int *timetable_y_end = new int[nMachines()];

  int *duration_x = new int[nMachines()];
  int *duration_y = new int[nMachines()];

  // the slack is the maximum right shift of y staying right of the given machine
  int *slack = new int[nMachines()];
  // the delay is the minimum right shift of y to get past the given machine
  int *delay = new int[nMachines()];

  //BitSet *timetable_y = new BitSet[nMachines()];

  
  int dur = 0;
  for(int k=0; k<nTasksInJob(x); ++k) {
    int t = getJobTask(x,k);
    for(int i=0; i<nMachines(t); ++i) {
      int m = getMachine(t,i);
      duration_x[m] = getDuration(t);
      timetable_x_start[m] = dur;
      timetable_x_end[m] = dur+getDuration(t);
    }
    dur += getDuration(t);
  }

  dur = -y_length;
  for(int k=0; k<nTasksInJob(y); ++k) {
    int t = getJobTask(y,k);
    for(int i=0; i<nMachines(t); ++i) {
      int m = getMachine(t,i);
      duration_y[m] = getDuration(t);
      timetable_y_start[m] = dur;
      timetable_y_end[m] = dur+getDuration(t);
    }
    dur += getDuration(t);
  }

  while(true) {
    int lb_shift=INFTY;
    int ub_shift=0;

    for(int i=0; i<nMachines(); ++i) {
      slack[i] = timetable_x_start[i]-timetable_y_end[i];
      delay[i] = timetable_x_end[i]-timetable_y_start[i];
      if(slack[i]<0) slack[i] = INFTY;

      if(slack[i]<lb_shift) lb_shift=slack[i];
      if(delay[i]>ub_shift) ub_shift=delay[i];
    }
  
    std::cout << "job" << x << ": " << x_length << " ";
    for(int k=0; k<nTasksInJob(x); ++k) {
      int t = getJobTask(x,k);
      for(int i=0; i<nMachines(t); ++i) {
	int m = getMachine(t,i);
	std::cout << "[" << std::setw(4) << timetable_x_start[m] << ":"
		  << std::setw(2) << getDuration(t) 
		  << " " << m
		  << "]";
      }
    }
    std::cout << std::endl;
    
    std::cout << "job" << y << ": " << y_length << " ";
    for(int k=0; k<nTasksInJob(y); ++k) {
      int t = getJobTask(y,k);
      for(int i=0; i<nMachines(t); ++i) {
	int m = getMachine(t,i);
	std::cout << "[" << std::setw(4) << timetable_y_start[m] << ":"
		  << std::setw(2) << getDuration(t) 
		  << " " << m
		  << "]";
      }
    }
    std::cout << std::endl;

    std::cout << "job" << y << ": " << y_length << " ";
    for(int k=0; k<nTasksInJob(y); ++k) {
      int t = getJobTask(y,k);
      for(int i=0; i<nMachines(t); ++i) {
	int m = getMachine(t,i);
	std::cout << "[" << std::setw(4) << slack[m] << "," << std::setw(4) << delay[m] << "]";
      }
    }
    std::cout << std::endl;

    //std::cout << min_forced_shift << " " << max_accepted_shift << " " << max_jump_shift << std::endl;
  

    
    // the shift should be, for each machine, either less than the slack, or more than the delay
    // therefore, we compute a forbidden interval, between the minimum slack, and the maximum delay
    if(lb_shift>=0) {
      std::cout << lb_shift << "]";
    }
    if(ub_shift>=0) {
      std::cout << "[" << ub_shift ;
    }

    std::cout << std::endl;

    
    for(int i=0; i<nMachines(); ++i) {
      timetable_y_start[i] += ub_shift;
      timetable_y_end[i] += ub_shift;
    }
  }

  exit(1);
}

int Instance::addTask(const int dur, const int job, const int machine) {
  int index = duration.size();
  if(job >= 0) addTaskToJob(index, job);
  if(machine >= 0) addTaskToMachine(index, machine);
  duration.push_back(dur);
  due_date.push_back(INFTY);
  release_date.push_back(0);
  
  return index;
}

void Instance::addTaskToJob(const unsigned int index, const unsigned int j) {
  if(tasks_in_job.size() <= j) tasks_in_job.resize(j+1);
  if(jobs_of_task.size() <= index) jobs_of_task.resize(index+1);
  if(task_rank_in_job.size() <= index) task_rank_in_job.resize(index+1);
  tasks_in_job[j].push_back(index);
  jobs_of_task[index].push_back(j);
  pair_ x(j, tasks_in_job[j].size()-1);
  task_rank_in_job[index].push_back(x);
}

void Instance::addTaskToMachine(const unsigned int index, const unsigned int j) {
  if(tasks_in_machine.size() <= j) tasks_in_machine.resize(j+1);
  if(machines_of_task.size() <= index) machines_of_task.resize(index+1);
  if(task_rank_in_machine.size() <= index) task_rank_in_machine.resize(index+1);
  tasks_in_machine[j].push_back(index);
  machines_of_task[index].push_back(j);
  pair_ x(j, tasks_in_machine[j].size()-1);
  task_rank_in_machine[index].push_back(x);
}

int Instance::getSetupTime(const int k, const int i, const int j) const {
  // get the rank of task i in machine k
  int ri = 0;
  for(unsigned int x=0; x<task_rank_in_machine[i].size(); ++x)
    if(task_rank_in_machine[i][x].element == k) {
      ri = task_rank_in_machine[i][x].rank;
      break;
    }

  int rj = 0;
  for(unsigned int x=0; x<task_rank_in_machine[j].size(); ++x)
    if(task_rank_in_machine[j][x].element == k) {
      rj = task_rank_in_machine[j][x].rank;
      break;
    }
  
  return setup_time[k][ri][rj];
}

std::ostream& Instance::print(std::ostream& os) {
  os << "c " << (nJobs()) << " jobs, " 
     << nMachines() << " machines ("
     << nTasks() << " tasks)" << std::endl;
  for(int i=0; i<nJobs(); ++i) {
    if(nTasksInJob(i) > 1) {
      os << "c ";
      for(int j=1; j<nTasksInJob(i); ++j)
	os << "  t" << tasks_in_job[i][j-1] << "+" << (duration[tasks_in_job[i][j-1]]) 
	   << " <= t" << tasks_in_job[i][j];
      os << std::endl;
    }
  }
  for(int i=0; i<nMachines(); ++i) {
    if(tasks_in_machine[i].size() > 0) {
      os << "c machine" << i << ": t" << tasks_in_machine[i][0];
      for(unsigned int j=1; j<tasks_in_machine[i].size(); ++j)
	os << ", t" << tasks_in_machine[i][j];
      os << std::endl;
    }
  }

//   for(int i=0; i<task_rank_in_machine.size(); ++i) {
//     std::cout << "task_" << i << " is"; 
//     for(int j=0; j<task_rank_in_machine[i].size(); ++j) {
//       pair_ p = task_rank_in_machine[i][j];
//       std::cout << " the " << p.rank << "/" << getMachine(i,j) << "th in machine_" << p.element;
//       if(task_rank_in_machine[i][j].element != machines_of_task[i][j]) {
// 	std::cout << "INCONSISTENCY" << std::endl;
// 	exit(1);
//       }
//     }
//     std::cout << std::endl;
//   }

  return os;
}

int Instance::nDisjuncts() const {
  int n_disjuncts = 0;
  for(int i=0; i<nMachines(); ++i) {
    n_disjuncts += (nTasksInMachine(i) * (nTasksInMachine(i)-1))/2;
  }
  return n_disjuncts;
}

int Instance::nPrecedences() const {
  int n_precedences = 0;
  for(int i=0; i<nJobs(); ++i) {
    n_precedences += (nTasksInJob(i)-1);
  }
  if(hasTimeLag()) {
    //n_precedences *= 2;
    for(int i=0; i<nJobs(); ++i) 
      for(int j=1; j<nTasksInJob(i); ++j) 
	if(getMaxLag(i,j-1) >= 0) 
	  ++n_precedences;
  }
  return n_precedences;
}

double Instance::getNormalizer() const {
  double normalizer=0.0, cost;
  int i, j, job_dur;
  for(i=0; i<nJobs(); ++i) {
    job_dur = 0;
    for(j=0; j<nTasksInJob(i); ++j) job_dur += getDuration(getJobTask(i,j));
    if(hasFloatCost())
      cost = getJobFloatCost(i);
    else if(hasEarlyCost() || hasLateCost()) {
      cost = 0.0;
      if(hasEarlyCost()) cost += getJobEarlyCost(i);
      if(hasLateCost()) cost += getJobLateCost(i);
    } else cost = 1.0;
    normalizer += ((double)job_dur * cost);
  }
  return normalizer;
}

std::ostream& Instance::printStats(std::ostream& os) {
  os << "c ==================[ instance ]===================" << std::endl
     << "d " << std::left << std::setw(28)  << "NUMTASKS "      << std::right << std::setw(21) << nTasks() << std::endl
     << "d " << std::left << std::setw(28)  << "NUMJOBS "       << std::right << std::setw(21) << nJobs() << std::endl
     << "d " << std::left << std::setw(28)  << "NUMMACHINES "   << std::right << std::setw(21) << nMachines() << std::endl
     << "d " << std::left << std::setw(28)  << "NUMDISJUNCTS "  << std::right << std::setw(21) << nDisjuncts() << std::endl
     << "d " << std::left << std::setw(28)  << "NUMPRECEDENCES "<< std::right << std::setw(21) << nPrecedences() << std::endl
//     << "d " << std::left << std::setw(28)  << "LBMAKESPAN "    << std::right << std::setw(21) << lb_C_max << std::endl
//     << "d " << std::left << std::setw(28)  << "UBMAKESPAN "    << std::right << std::setw(21) << ub_C_max << std::endl
     << "c ==================[ instance ]===================" << std::endl;
  return os;
}

int Instance::getMakespanLowerBound() {
  int mkp = 0, length;
  for(int i=0; i<nJobs(); ++i) {
    length = 0;
    for(int j=0; j<nTasksInJob(i); ++j)
      length += getDuration(getJobTask(i,j));
    if(mkp < length) mkp = length;
  }
  for(int i=0; i<nMachines(); ++i) {
    length = 0;
    for(int j=0; j<nTasksInMachine(i); ++j)
      length += getDuration(getMachineTask(i,j));
    if(mkp < length) mkp = length;
  }

  DBG("Get instance's makespan lb (%d)\n", mkp);

  return mkp;
}

int Instance::getMakespanUpperBound(const int iterations) {


  if(max_makespan < INFTY) return max_makespan;

  int best_makespan = INFTY;
  if(!hasTimeLag()) {
    int current_job[nJobs()];
    int current_job_bound[nJobs()];
    int current_machine_bound[nMachines()];
    int current_machine_task[nMachines()];
    int ranks[nTasks()];
    int random_jobs[nTasks()+1];
    int m[nMachines()];
    int k=0, i, j, t, n=0;

    std::fill(current_job, current_job+nJobs(), 0);
    while(k<nTasks()) {
      ranks[n++] = k;
      for(i=0; i<nJobs(); ++i)
	if(current_job[i]<nTasksInJob(i)) {
	  random_jobs[k++] = i;
	  ++current_job[i];
	  //std::cout << " " << i ;
	}
    }
    //std::cout << std::endl;
    ranks[n] = nTasks();

    int iter = iterations;
  
    while(iter--) {
      std::fill(current_job, current_job+nJobs(), 0);
      std::fill(current_job_bound, current_job_bound+nJobs(), 0);
      std::fill(current_machine_bound, current_machine_bound+nMachines(), 0);
      std::fill(current_machine_task, current_machine_task+nMachines(), -1);
  
      int makespan = 0;

      if(iter<iterations-1)
	for(i=0; i<n; ++i) {
	  for(t=ranks[i]; t<ranks[i+1]; ++t) {
	    j = randint(ranks[i+1]-t);
	    k = random_jobs[t];
	    random_jobs[t] = random_jobs[ranks[i]+j];
	    random_jobs[ranks[i]+j] = k;
	  }
	}
      //      for(i=0; i<nTasks(); ++i) {
      //	j = randint(nTasks()-i);
      // 	k = random_jobs[i];
      // 	random_jobs[i] = random_jobs[i+j];
      // 	random_jobs[i+j] = k;
      //       }
    
      for(i=0; i<nTasks(); ++i) {

// 	for(int jj=0; jj<nJobs(); ++jj) {
// 	  std::cout << " " << std::setw(4) << current_job_bound[jj];
// 	}
// 	std::cout << std::endl;

// 	for(int jj=0; jj<nMachines(); ++jj) {
// 	  std::cout << " " << std::setw(4) << current_machine_bound[jj];
// 	}
// 	std::cout << std::endl;


	// pick the next job
	j = random_jobs[i];


	// find out which task is that
	t = getJobTask(j, current_job[j]);

	DBG("pick task t%d\n", t);

	// find out which machine is that
	for(k=0; k<nMachines(t); ++k) {
	  m[k] = getMachine(t,k);
	  DBG("  -> uses machine%d\n", m[k]);
	}
      
// 	std::cout << "pick task " << t << "(job=" << j 
// 		  << ", mach=" << m[0] << ")" << std::endl;


	// find the current timepoint for this job
	int j_mkp = current_job_bound[j];

	// find the current timepoint for this machine
	int m_mkp = current_machine_bound[m[0]];
	DBG("m%d = %d\n", m[0], current_machine_bound[m[0]]);
	for(k=1; k<nMachines(t); ++k)
	  if(m_mkp < current_machine_bound[m[k]]) {
	    m_mkp = current_machine_bound[m[k]];
	    DBG("m%d = %d\n", m[k], current_machine_bound[m[k]]);
	  }

	// get the start time for this task
	int timepoint = (j_mkp < m_mkp ? m_mkp : j_mkp);

	//	std::cout << "earliest start time = " << timepoint << std::endl;

	// check its release date
	if(getReleaseDate(t) > timepoint) timepoint = getReleaseDate(t);

	DBG("timepoint = %d\n", timepoint);

	// add setup time, if any
	if(hasSetupTime()) {
	  int setup = 0;
	  int setup_mk;
	  for(k=0; k<nMachines(t); ++k) {
	    if(current_machine_task[m[k]] >= 0) {
	      setup_mk = getSetupTime(m[k], current_machine_task[m[k]], t);
	      if(setup < setup_mk) setup = setup_mk;
	      DBG("setup = %d\n", setup_mk);
	    }
	  }
	  timepoint += setup;
	}

	// get the finish time for this task
	timepoint += getDuration(t);

	// update machin and job bounds
	for(k=0; k<nMachines(t); ++k) {
	  current_machine_bound[m[k]] = timepoint;
	  current_machine_task[m[k]] = t;
	}
	//current_machine_bound[m] = timepoint;
	current_job_bound[j] = timepoint;

	// get the final makespan
	if(makespan < timepoint) makespan = timepoint;

	++current_job[j];
      }
      if(best_makespan > makespan) best_makespan = makespan;

      //exit(1);
      //std::cout << "\t" << makespan << " " << best_makespan << std::endl;
    }
  } else {

    best_makespan = 0;
    for(int t=0; t<nTasks(); ++t) best_makespan += getDuration(t);

  }

  DBG("Get instance's makespan ub (%d)\n", best_makespan);

  return best_makespan;
}

int Instance::getEarlinessTardinessLowerBound(const int c_max) {
  return 0;
}
int Instance::getEarlinessTardinessUpperBound(const int c_max) {
  int i, ti, sum_ub = 0;

  if(hasEarlyCost()) 
    for(i=0; i<nJobs(); ++i) {
      ti = getLastTaskofJob(i);
      sum_ub += ((getJobDueDate(i) - (getReleaseDate(ti) + getDuration(ti)))*getJobEarlyCost(i));
    }
    
  if(hasLateCost()) 
    for(i=0; i<nJobs(); ++i) {
      ti = getLastTaskofJob(i);
      sum_ub += ((c_max - getJobDueDate(i))*getJobLateCost(i));
    }

  return sum_ub;
}


void Instance::osp_readData( const char* filename ) {

  DBG("Read (osp)%s\n", "");

  int opt, lb, i=0, j, k, nJobs, nMachines, dur, bufsize=1000;
   char buf[bufsize];
   std::ifstream infile( filename, std::ios_base::in );
	
   do {
     infile.getline( buf, bufsize, '\n' );
   } while( buf[0] == '#' );
	
   while( buf[i] != ' ' ) ++i;
   buf[i] = '\0';
   lb = atoi( buf );
	
   while( buf[i] == '\0' || buf[i] == ' ' || buf[i] == '*' ) ++i;
   j = i;
   while( buf[i] != ' ' && buf[i] != '\n' && buf[i] != '\0' ) ++i;
   buf[i] = '\0';
   opt = atoi( &(buf[j]) );
	
   do {
     infile.get( buf[0] );
     if( buf[0] != '#' ) break;
     infile.getline( buf, bufsize, '\n' );
   } while( true );
   infile.unget();
	
   infile >> nJobs;
   infile >> nMachines;

   infile.getline( buf, bufsize, '\n' );
	
   do {
     infile.get( buf[0] );
     if( buf[0] != '#' ) break;
     infile.getline( buf, bufsize, '\n' );
   } while( true );
   infile.unget();
	
   k = 0;
   for(i=0; i<nJobs; ++i) {
     for(j=0; j<nMachines; ++j) {
       infile >> dur;
       addTask(dur, k, -1);

       addTaskToMachine(k, i);
       addTaskToMachine(k, nJobs+j);       

       ++k;
     }
   }
}

void Instance::sds_readData( const char* filename ) {

  DBG("Read (sds)%s\n", "");

  int i, j, nJobs, nMachines, nFamilies, dur, mach;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> nMachines;
  infile >> nJobs;
  infile >> nFamilies;

  int **family_matrix = new int*[nFamilies+1];
  for(i=0; i<=nFamilies; ++i)
    family_matrix[i] = new int[nFamilies];
	
  setup_time = new int**[nMachines];
  for(i=0; i<nMachines; ++i) {
    setup_time[i] = new int*[nJobs];
    for(j=0; j<nJobs; ++j) {
      setup_time[i][j] = new int[nJobs];
    }
  }
	
  int **family = new int*[nJobs];
  for(i=0; i<nJobs; ++i) 
    family[i] = new int[nMachines];

  for(i=0; i<nJobs; ++i) {
		
    infile >> j;
    assert(j==nMachines);
		
    for(j=0; j<nMachines; ++j) {

      infile >> dur;
      infile >> mach;
      --mach;

      addTask(dur, i, mach);
      //addTaskToMachine(k++, mach);
			
      infile >> family[i][j];
      --family[i][j];
    }
  }
	
  for(i=0; i<=nFamilies; ++i)
    for(j=0; j<nFamilies; ++j)
      infile >> family_matrix[i][j];
	
  for(int k=0; k<nMachines; ++k) {
    for(i=0; i<nJobs; ++i) {
      for(j=0; j<nJobs; ++j) {
	setup_time[k][i][j] = family_matrix[1+family[i][k]][family[j][k]] ;
	std::cout << " " << setup_time[k][i][j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  
  int k=0;
  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j)
      release_date[k++] = family_matrix[0][family[i][j]];
  }	

}
void Instance::jtl_readData( const char* filename ) {

  DBG("Read (jtl)%s\n", "");

  int i, j, dur, mach, nJobs, nMachines, opt;
  std::string tag;
  char c;
  std::ifstream infile( filename, std::ios_base::in );
  
  infile >> nJobs;
  infile >> nMachines;

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {

      infile >> mach;
      infile >> dur;

      addTask(dur, i, mach);
    }
  }


  infile >> tag;
  infile >> opt;
  infile.ignore( 100, '\n' );
  
  infile.get(c);
  assert( c == 'T' );
  infile.get(c);
  assert( c == 'L' );
  infile.get(c);
  assert( c == '=' );

  
  time_lag[0] = new int*[nJobs];
  time_lag[1] = new int*[nJobs];

  for(i=0; i<nJobs; ++i) {
    time_lag[0][i] = new int[nMachines];
    time_lag[1][i] = new int[nMachines];
    for(j=0; j<nMachines; ++j) {
      infile >> time_lag[0][i][j];
      infile >> time_lag[1][i][j];
    }
  }

}
void Instance::now_readData( const char* filename ) {

  DBG("Read (now)%s\n", "");

  int i, j, dur, mach, nJobs, nMachines;

  std::string tag;

  std::ifstream infile( filename, std::ios_base::in );

  infile >> nJobs;

  infile >> nMachines;

  for(i=0; i<nJobs; ++i) {

    for(j=0; j<nMachines; ++j) {

      infile >> mach;

      infile >> dur;

      addTask(dur, i, mach);

    }

  }

  time_lag[0] = new int*[nJobs];

  time_lag[1] = new int*[nJobs];

  for(i=0; i<nJobs; ++i) {

    time_lag[0][i] = new int[nMachines];

    time_lag[1][i] = new int[nMachines];

    for(j=0; j<nMachines; ++j) {

      time_lag[0][i][j] = 0;

      time_lag[1][i][j] = 0;

    }

  }

//   print(std::cout);
//   std::cout << std::endl;
 
//   Vector<int> intervals;

//   for(int i=0; i<nJobs; ++i) {
//     for(int j=i+1; j<nJobs; ++j) {
//       intervals.clear();
//       get_placement(i,j,intervals);
//       intervals.print(std::cout);
//     }
//     std::cout << std::endl;
//   }
//   //exit(1);
}
void Instance::jla_readData( const char* filename ) {

  DBG("Read (jla)%s\n", "");

  int i, j, dur, mach, nJobs, nMachines;
  std::ifstream infile( filename, std::ios_base::in );
  
  infile >> nJobs;
  infile >> nMachines;

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {

      infile >> mach;
      infile >> dur;

      addTask(dur, i, mach);
    }
  }

}
void Instance::tsp_readData( const char* filename ) {

  DBG("Read (tsp)%s\n", "");

}
void Instance::fsp_readData( const char* filename ) {

  DBG("Read (fsp)%s\n", "");

  int dur;
  std::vector<int> duration;
  double obj;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> obj;
  max_makespan = (int)obj;

  while( true ) {
    infile >> dur;
    if(!infile.good()) break;
    duration.push_back(dur);

    //std::cout << dur << std::endl;
  }

  int nJobs = duration.size()/2;

  for(int i=0; i<nJobs; ++i) {
    addTask(duration[i], i, 0);
    addTask(duration[i+nJobs], i, 1);
  }

  //this->print(std::cout);

  //exit(1);
}

void Instance::jsp_readData( const char* filename ) {

  DBG("Read (jsp)%s\n", "");

  int i, j, k, dur, mach;
  long int dump;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  int nJobs;
  int nMachines;
  infile >> nJobs;
  infile >> nMachines;
	
  infile >> dump;
  infile >> dump;
	
  infile >> dump;
  infile >> dump;
	
  infile >> tag;

  assert( tag == "Times" );

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      infile >> dur;
      addTask(dur, i, -1);
    }
  }
	
  infile >> tag;
  assert( tag == "Machines" );
  
  k = 0;
  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      infile >> mach;
      addTaskToMachine(k++, --mach);
    }
}

void Instance::jet_readData( const char* filename ) {

  DBG("Read (jet)%s\n", "");

  int i, j, dur, mach, nJobs, nMachines;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> nJobs;
  infile >> nMachines;
	
  jsp_duedate = new int[nJobs];
  jsp_earlycost = new int[nJobs];
  jsp_latecost = new int[nJobs];
	
  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
			
      infile >> mach;
      infile >> dur;
      addTask(dur, i, mach);

    }
		
    infile >> jsp_duedate[i];
    infile >> jsp_earlycost[i];
    infile >> jsp_latecost[i];
  }

}

void Instance::dyn_readData( const char* filename, const int precision ) {

  DBG("Read (dyn)%s\n", "");

  int i, j, k, dur, mach, nJobs, nMachines;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> nJobs;
  infile >> nMachines;
	
  jsp_duedate = new int[nJobs];
  jsp_earlycost = new int[nJobs];
  jsp_latecost = new int[nJobs];
  jsp_floatcost = new double[nJobs];

  for(i=0; i<nJobs; ++i) {
    k = release_date.size();

    for(j=0; j<nMachines; ++j) {
      
      infile >> mach;
      infile >> dur;
      
      if(mach != -1 && dur != -1) addTask(dur, i, mach);
    }		

    infile >> jsp_duedate[i];
    infile >> j;
    //infile >> jsp_latecost[i];
    infile >> jsp_floatcost[i];
    
    setReleaseDate(k, j);
    //jsp_earlycost[i] = jsp_latecost[i];
  }


  double min_cost = jsp_floatcost[0];
  double max_cost = jsp_floatcost[0];
  for(i=1; i<nJobs; ++i) {
    if(min_cost > jsp_floatcost[i]) min_cost = jsp_floatcost[i];
    if(min_cost < jsp_floatcost[i]) max_cost = jsp_floatcost[i];
  }
  for(i=0; i<nJobs; ++i) {
    int approx = (int)((jsp_floatcost[i]/max_cost)*precision);
    jsp_earlycost[i] = jsp_latecost[i] = approx;
    //std::cout << approx << std::endl;
  }

}

void Instance::dtp_readData( const char* filename ) {

  DBG("Read (dtp)%s\n", "");

  char comment = '#';
  std::ifstream infile( filename, std::ios_base::in );

  while(true) {
    comment = infile.get();
    if(comment == '#') {
      infile.ignore(1000, '\n');
    } else break;
  }
  infile.unget();
  
  std::string exp;
  int i, j, x, y, d;
  int step = 0;

  while(true) {
    infile >> exp;
    if(!infile.good())  break;

    if(step == 0) {
      std::vector< Term > clause;
      dtp_clauses.push_back(clause);
      dtp_weights.push_back(randint(100)+1);
    }

    //

    if(step%2 == 0) {
      i = 0;
      while(exp[++i] != '-');
      x = atoi(exp.substr(1,i-1).c_str());
      j = i+2;
      while(exp[++j] != '<');
      y = atoi(exp.substr(i+2,j-i-2).c_str());
      d = atoi(exp.substr(j+2,10).c_str());

      Term t;
      t.X = x;
      t.Y = y;
      t.k = d;

      if(t.X > dtp_nodes) dtp_nodes = t.X;
      if(t.Y > dtp_nodes) dtp_nodes = t.Y;

      dtp_clauses.back().push_back(t);

      //std::cout << "t" << x << " - t" << y  << " <= " << d << std::endl << std::endl;
    }

    step = (step+1) % 4;
  }

    ++dtp_nodes;

  // for(int i=0; i<dtp_clauses.size(); ++i) {
  //   for(int j=0; j<dtp_clauses[i].size(); ++j) {
  //     std::cout << "t" << dtp_clauses[i][j].X << " - t" << dtp_clauses[i][j].Y << " <= " << dtp_clauses[i][j].k << " OR ";
  //   }
  //   std::cout << std::endl;
  // }


}



SchedulingModel::SchedulingModel(Instance& inst, ParameterList *params, 
				 const int max_makespan) : CSP() {

  setup(inst, params, max_makespan);

}

std::ostream& SchedulingModel::printStats(std::ostream& os) {
  os << "c ====================[ model ]====================" << std::endl
     << "d " << std::left << std::setw(28)  << "LBMAKESPAN "    << std::right << std::setw(21) << lb_C_max << std::endl
     << "d " << std::left << std::setw(28)  << "UBMAKESPAN "    << std::right << std::setw(21) << ub_C_max << std::endl
     << "c ====================[ model ]====================" << std::endl;
  return os;
}


void SchedulingModel::setup(Instance& inst, 
			    ParameterList *params, 
			    const int max_makespan) {

  int i,j,k, lb, ub, ti, tj, rki, rkj, hi, hj, aux;
  
  lb_C_max = (params->LBinit<0 ? inst.getMakespanLowerBound() : params->LBinit);
  if(max_makespan < 0) ub_C_max = (params->UBinit<0 ? inst.getMakespanUpperBound(params->InitBound) : params->UBinit);
  else ub_C_max = max_makespan; 
  if(params->Objective == "tardiness") {
    int max_due_date = inst.getJobDueDate(0);
    for(int i=1; i<inst.nJobs(); ++i) {
      if(max_due_date < inst.getJobDueDate(i))
	max_due_date = inst.getJobDueDate(i);
    }
    if(max_due_date < ub_C_max) max_due_date = ub_C_max;
    ub_C_max += max_due_date;
  }

  //for(i=0; i<inst.nJobs(); ++i)
  //ub_C_max += inst.getDuration(inst.getLastTaskofJob(i));

  lb_L_sum = inst.getEarlinessTardinessLowerBound(ub_C_max);
  ub_L_sum = inst.getEarlinessTardinessUpperBound(ub_C_max);


  data = &inst;

  // create one variable per task
  for(i=0; i<inst.nTasks(); ++i) {
    lb = inst.getReleaseDate(i);
    ub = std::min(ub_C_max, inst.getDueDate(i)) - inst.getDuration(i);

    if(lb > ub) {
      std::cout << "INCONSISTENT" << std::endl;
      exit(1);
    }

    Variable t(lb, ub);
    tasks.add(t);
  }

  // precedence constraints
  for(i=0; i<inst.nJobs(); ++i) 
    for(j=1; j<inst.nTasksInJob(i); ++j) {
      ti = inst.getJobTask(i,j-1);
      tj = inst.getJobTask(i,j);
      add( Precedence(tasks[ti], 
		      (inst.getDuration(ti) + 
		       (inst.hasTimeLag() ? inst.getMinLag(i,j-1) : 0)), 
		      tasks[tj]) );
    }

  // time lags constraints
  if(inst.hasTimeLag()) {
    for(i=0; i<inst.nJobs(); ++i) 
      for(j=1; j<inst.nTasksInJob(i); ++j) if(inst.getMaxLag(i,j-1) >= 0) {
	  ti = inst.getJobTask(i,j-1);
	  tj = inst.getJobTask(i,j);
	  add( Precedence(tasks[tj], 
			  -(inst.getDuration(ti)+inst.getMaxLag(i,j-1)), 
			  tasks[ti]) );
	}
  }

  // mutual exclusion constraints
  for(k=0; k<inst.nMachines(); ++k) 
    for(i=0; i<inst.nTasksInMachine(k); ++i)
      for(j=i+1; j<inst.nTasksInMachine(k); ++j) {
	ti = inst.getMachineTask(k,i);
	tj = inst.getMachineTask(k,j);
	
	if(params->OrderTasks==1) {
	  rki = inst.getRankInJob(ti);
	  rkj = inst.getRankInJob(tj);
	  hi = inst.getHeadInJob(ti);
	  hj = inst.getHeadInJob(tj);
	  
	  if(rkj < rki || (rkj == rki && hj < hi))
	    {
	      aux = ti;
	      ti = tj;
	      tj = aux;
	    }
	} else if(params->OrderTasks==-1) {
      
	  aux = ti;
	  ti = tj;
	  tj = aux;

	}

//  	std::cout << "t" << ti << " <> t" << tj 
//  		  << " (" << inst.getDuration(ti) << "." 
// 		  << inst.getDuration(tj) << ")" 
// 		  << std::endl; 

	first_job.push(inst.getJob(ti,0));
	second_job.push(inst.getJob(tj,0));

	first_task_of_disjunct.push_back(ti);
	second_task_of_disjunct.push_back(tj);

	disjuncts.add( Disjunctive( tasks[ti],
				    
				    inst.getDuration(ti)
				    +(inst.hasSetupTime() ? inst.getSetupTime(k,ti,tj) : 0),
				    
				    tasks[tj],
				    
				    inst.getDuration(tj)
				    +(inst.hasSetupTime() ? inst.getSetupTime(k,tj,ti) : 0) ) );
      }
  add( disjuncts );

  //exit(1);

  Variable x_cmax(lb_C_max, ub_C_max);
  C_max = x_cmax;
  for(i=0; i<inst.nJobs(); ++i) {
    ti = inst.getLastTaskofJob(i);
    add(Precedence(tasks[ti], inst.getDuration(ti), C_max));
  }

  if(inst.hasJobDueDate()) {
    // early/late bools - whether the last tasks are early or late
    for(i=0; i<inst.nJobs(); ++i) {
      ti = inst.getLastTaskofJob(i);
      if(inst.hasEarlyCost())
	earlybool.add(tasks[ti] < (inst.getJobDueDate(i) - inst.getDuration(ti)));
      if(inst.hasLateCost())
	latebool.add(tasks[ti] > (inst.getJobDueDate(i) - inst.getDuration(ti)));
    }
    
    //DG Added for calculating cost
    VarArray sum_scope;
    int weights[2*inst.nJobs()+1];
    int sum_ub = 0;
    int n_weights = 0;

    k=0;
    if(inst.hasEarlyCost()) {
      for(i=0; i<inst.nJobs(); ++i) {
	ti = inst.getLastTaskofJob(i);
	Variable x(0, inst.getJobDueDate(i)-inst.getDuration(ti)); // amount of "earliness"
	add(x == (earlybool[i]*((tasks[ti]*-1) - (inst.getDuration(ti) - inst.getJobDueDate(i)))));
	sum_scope.add(x);
	weights[k++] = inst.getJobEarlyCost(i);
	sum_ub += ((inst.getJobDueDate(i) - (inst.getReleaseDate(ti) + inst.getDuration(ti)))*inst.getJobEarlyCost(i));
      }
      n_weights += inst.nJobs();
    }
    
    if(inst.hasLateCost()) {
      for(i=0; i<inst.nJobs(); ++i) {
	ti = inst.getLastTaskofJob(i);
	Variable x(0, ub_C_max-inst.getJobDueDate(i));
	add(x == (latebool[i]*(tasks[ti] + (inst.getDuration(ti) - inst.getJobDueDate(i)))));
	sum_scope.add(x);
	weights[k++] = inst.getJobLateCost(i);
	sum_ub += ((ub_C_max - inst.getJobDueDate(i))*inst.getJobLateCost(i));
      }
      n_weights += inst.nJobs();
    }
    
    weights[n_weights] = 0;
    
    Variable x_lsum(0, sum_ub);
    L_sum = x_lsum;
    add( L_sum == Sum(sum_scope, weights) );
  }


  // initialise last tasks
  for(i=0; i<inst.nJobs(); ++i) last_tasks.add(tasks[inst.getLastTaskofJob(i)]);

  // initialise searched vars
  if(data->hasJobDueDate()) {
    for(int i=0; i<earlybool.size(); ++i) SearchVars.add(earlybool[i]);
    for(int i=0; i< latebool.size(); ++i) SearchVars.add(latebool [i]);
    for(int i=0; i<data->nJobs(); ++i) SearchVars.add(tasks[data->getLastTaskofJob(i)]);
  }
  for(int i=0; i<disjuncts.size(); ++i) SearchVars.add(disjuncts[i]);
  if(params->FixTasks)
    for(int i=0; i<tasks.size(); ++i) SearchVars.add(tasks[i]);

  lb_Depth = 0;
  ub_Depth = disjuncts.size();

  if(params->Objective == "depth") {
    Variable depth(lb_Depth, ub_Depth);
    Depth = depth;
    VarArray scope;
    for(int i=0; i<disjuncts.size()/2; ++i) {
      add(disjuncts[i] == disjuncts[i+disjuncts.size()/2]);
      scope.add(disjuncts[i]);
    }
    add(depth == Sum(scope));
  }
}

SchedulingModel::~SchedulingModel() {}


void No_wait_Model::setup(Instance& inst, ParameterList *params, const int max_makespan) {

  //if(type==0) {
  int i,j,k, lb, ub, ti=0, tj=0;
  std::vector<int> offset;
  for(i=0; i<inst.nJobs(); ++i) {
    offset.push_back(0);
    for(j=0; j<inst.nTasksInJob(i)-1; ++j)
      offset.push_back(offset.back() + inst.getDuration(inst.getJobTask(i,j)));
  }

  lb_C_max = inst.getMakespanLowerBound();
  if(max_makespan < 0) {
    ub_C_max = 0;
    for(i=0; i<inst.nJobs(); ++i) {
      ub_C_max += (offset[inst.getLastTaskofJob(i)] + inst.getDuration(inst.getLastTaskofJob(i)));
    }
  } else ub_C_max = max_makespan;

  data = &inst;

  // create one variable per task
  for(i=0; i<inst.nJobs(); ++i) {

    lb = inst.getReleaseDate(inst.getJobTask(i,0));
    ub = std::min(ub_C_max, inst.getDueDate(inst.getLastTaskofJob(i))) 
      - (offset[inst.getLastTaskofJob(i)] + inst.getDuration(inst.getLastTaskofJob(i)));


    if(lb > ub) {
      std::cout << "INCONSISTENT" << std::endl;
      exit(1);
    }

    Variable t(lb, ub);
    tasks.add(t);
  }

  if(type==0) {

    // mutual exclusion constraints
    for(k=0; k<inst.nMachines(); ++k) 
      for(i=0; i<inst.nTasksInMachine(k); ++i)
	for(j=i+1; j<inst.nTasksInMachine(k); ++j) {
	  ti = inst.getMachineTask(k,i);
	  tj = inst.getMachineTask(k,j);
	  
	  first_task_of_disjunct.push_back(ti);
	  second_task_of_disjunct.push_back(tj);

	  disjuncts.add( Disjunctive( tasks[inst.getJob(ti,0)],
				      
				      inst.getDuration(ti)
				      +(inst.hasSetupTime() ? inst.getSetupTime(k,ti,tj) : 0)
				      +offset[ti]
				      -offset[tj],
				      
				      tasks[inst.getJob(tj,0)],
				      
				      inst.getDuration(tj)
				      +(inst.hasSetupTime() ? inst.getSetupTime(k,tj,ti) : 0)
				      +offset[tj]
				      -offset[ti] ) );
	}
    add( disjuncts );
  } else {
    
    //       print(std::cout);
    //   std::cout << std::endl;

    if(true) {
    Vector<int> intervals;
    
    for(int i=0; i<data->nJobs(); ++i) {
      for(int j=i+1; j<data->nJobs(); ++j) {
	intervals.clear();
	data->get_placement(i,j,intervals);
	
	int job_start = disjuncts.size();
	for(int k=0; k<intervals.size; k+=2) {

	  first_task_of_disjunct.push_back(i);
	  second_task_of_disjunct.push_back(j);

	  //disjuncts.add( Disjunctive(tasks[j], -intervals[k], tasks[i], intervals[k+1]) );
	  disjuncts.add( Disjunctive(tasks[i], intervals[k+1], tasks[j], -intervals[k]) );
	}
	int ds = disjuncts.size();
  	for(int k=job_start; k<ds; ++k) {
//  	  add(disjuncts[k-1] <= disjuncts[k]);

	  job_size.push(ds-job_start);
	  job_index.push(k-job_start);
 	}

      	//gen_disjuncts.add( GenDisjunctive(tasks[i], tasks[j], intervals) );

      }
      //std::cout << std::endl;
    }
    
//     job_size.print(std::cout);
//     std::cout << std::endl;

//     job_index.print(std::cout);
//     std::cout << std::endl;

    add(disjuncts);


    } else {


    Vector<int> intervals;
    
    for(int i=0; i<data->nJobs(); ++i) {
      for(int j=i+1; j<data->nJobs(); ++j) {
	intervals.clear();
	data->get_placement(i,j,intervals);
	//intervals.print(std::cout);
      	gen_disjuncts.add( GenDisjunctive(tasks[i], tasks[j], intervals) );

      }
      //std::cout << std::endl;
    }
    
    add(gen_disjuncts);
    //exit(1);
    }
  }



  Variable x_cmax(lb_C_max, ub_C_max);
  C_max = x_cmax;
  for(i=0; i<inst.nJobs(); ++i) {
    ti = inst.getLastTaskofJob(i);
    k = inst.getDuration(ti) + offset[ti];
    add(Precedence(tasks[i], k, C_max));
  }

  
  if(disjuncts.size()>0) {
    for(int i=0; i<disjuncts.size(); ++i) SearchVars.add(disjuncts[i]);
  } 
  if(gen_disjuncts.size()>0) {
    for(int i=0; i<gen_disjuncts.size(); ++i) SearchVars.add(gen_disjuncts[i]);
  }

  

}


void DTP_Model::setup(Instance& inst, ParameterList *params, const int max_makespan) {

  ub_C_max = max_makespan; 
  ub_Weight = 0;
  for(unsigned int i=0; i<inst.dtp_clauses.size(); ++i) {
    ub_Weight += inst.dtp_weights[i];
  }
  inst.dtp_weights.push_back(0);

  for(int i=0; i<inst.dtp_nodes; ++i) {
    Variable t(0, max_makespan);
    tasks.add(t);
  }

  // x-y <= k
  for(unsigned int i=0; i<inst.dtp_clauses.size(); ++i) {
    for(int j=0; j<2; ++j) {
      SearchVars.add( Precedence( tasks[inst.dtp_clauses[i][j].X], 
				  -inst.dtp_clauses[i][j].k, 
				  tasks[inst.dtp_clauses[i][j].Y] ) );
    }
    int k = SearchVars.size();
    disjuncts.add( SearchVars[k-2] || SearchVars[k-1] );
  }

  // for(int i=0; i<disjuncts.size(); ++i) {
  //   SearchVars.add(disjuncts[i]);
  // }
  
  Weight = Variable(0,ub_Weight);
  add(Weight == (-Sum(disjuncts, &(inst.dtp_weights[0])) + ub_Weight));
  //  Weight = Sum(disjuncts, &(inst.dtp_weights[0]));
  //add(Maximise(Weight));
}


int L_sum_Model::get_lb() {
  return lb_L_sum;
}

int L_sum_Model::get_ub() {
  return ub_L_sum;
}

int Depth_Model::get_lb() {
  return lb_Depth;
}

int Depth_Model::get_ub() {
  return ub_Depth;
}

int DTP_Model::get_lb() {
  return 0;
}

int DTP_Model::get_ub() {
  return ub_Weight;
}

int C_max_Model::get_lb() {
  return lb_C_max;
}

int C_max_Model::get_ub() {
  return ub_C_max;
}

VariableInt* L_sum_Model::get_objective_var() {
  return L_sum.getVariable();
}

VariableInt* Depth_Model::get_objective_var() {
  return Depth.getVariable();
}

VariableInt* DTP_Model::get_objective_var() {
  return Weight.getVariable();
}

VariableInt* C_max_Model::get_objective_var() {
  return C_max.getVariable();
}

int L_sum_Model::set_objective(const int obj) {
  //VariableInt *x_lsum = L_sum.getVariable();
  return (L_sum.getVariable()->setMax(obj) ? UNKNOWN : UNSAT);
}

int Depth_Model::set_objective(const int obj) {
//   //VariableInt *x_lsum = Depth.getVariable();
//   Depth.print(std::cout);
//   std::cout << std::endl;
//   Depth.getVariable()->print(std::cout);
//   std::cout << std::endl;
  return (Depth.getVariable()->setMax(obj) ? UNKNOWN : UNSAT);
}

int DTP_Model::set_objective(const int obj) {
  return (Weight.getVariable()->setMax(obj) ? UNKNOWN : UNSAT);
}

int C_max_Model::set_objective(const int obj) {
  VariableInt *x_cmax = C_max.getVariable();
  bool consistent = true;
  if(x_cmax) {
    consistent = x_cmax->setMax(obj);
  } else {
    int i, ti;
    VariableInt *t;
    for(i=0; consistent && i<data->nJobs(); ++i) {
      ti = data->getLastTaskofJob(i);
      t = tasks[ti].getVariable();
      if(t) {
	consistent &= t->setMax(obj - data->getDuration(ti));
      } else {
	if(tasks[ti].value() > (obj-data->getDuration(ti))) consistent = false;
      }
    }
  }
  return (consistent ? UNKNOWN : UNSAT);
}

int L_sum_Model::get_objective() {
  return L_sum.value();
}

int Depth_Model::get_objective() {
  return Depth.value();
}

int DTP_Model::get_objective() {
  return Weight.value();
}

double L_sum_Model::get_normalized_objective() {
  double total_cost = 0, cost;
  for(int i=0; i<latebool.size(); ++i) {
    cost = 0.0;
    if(earlybool[i].value()) {
      cost = (data->getJobDueDate(i) - (last_tasks[i].value()+data->getDuration(data->getLastTaskofJob(i))))*(data->hasFloatCost() ? data->getJobFloatCost(i) : data->getJobEarlyCost(i));
    } else if(latebool[i].value()) {
      cost = ((last_tasks[i].value()+data->getDuration(data->getLastTaskofJob(i))) - data->getJobDueDate(i))*(data->hasFloatCost() ? data->getJobFloatCost(i) : data->getJobLateCost(i));
    }
    total_cost += cost;
  }
  
  total_cost /= data->getNormalizer();

  return total_cost;
}

int C_max_Model::get_objective() {
  return C_max.value();
}


Solution::Solution(SchedulingModel *m, SchedulingSolver *s, int target) {
  model = m;
  solver = s;

  int i, n=solver->variables.size;

  earlybool_value = NULL;
  latebool_value = NULL;
  task_min = NULL;
  task_max = NULL;
  ltask_value = NULL;
  disjunct_value = NULL;
  search_value = NULL;
  all_value = new int[n];
  min_objective_value = m->get_objective();
  max_objective_value = (target < 0 ? min_objective_value : target);

  

  for(i=0; i<n; ++i) {
    all_value[i] = solver->variables[i]->value();
//     std::cout << "store " ;
//     solver->variables[i]->print(std::cout);
//     std::cout << std::endl;
  }


  n = m->disjuncts.size();
  if(n) {
    disjunct_value = new int[n];
    for(i=0; i<n; ++i) disjunct_value[i] = m->disjuncts[i].value();
  } 

  n = m->SearchVars.size();
  if(n) {
    search_value = new int[n];
    for(i=0; i<n; ++i) search_value[i] = m->SearchVars[i].value();
  }

  n = m->tasks.size();
  if(n) {
    task_min = new int[n];
    task_max = new int[n];
    for(i=0; i<n; ++i) {
      task_min[i] = m->tasks[i].min();
      task_max[i] = m->tasks[i].max();
    } 
  }

  if(model->data->hasJobDueDate()) {
    n = m->earlybool.size();
    if(n) {
      earlybool_value = new int[n];
      for(i=0; i<n; ++i) earlybool_value[i] = m->earlybool[i].value();
    } 
    n = m->latebool.size();
    if(n) {
      latebool_value = new int[n];
      for(i=0; i<n; ++i) latebool_value[i] = m->latebool[i].value();
    } 
    // n = m->tasks.size();
    // if(n) {
    //   task_min = new int[n];
    //   task_max = new int[n];
    //   for(i=0; i<n; ++i) {
    // 	task_min[i] = m->tasks[i].min();
    // 	task_max[i] = m->tasks[i].max();
    //   } 
    // }

    n = m->last_tasks.size();
    if(n) {
      ltask_value = new int[n];
      for(i=0; i<n; ++i) {
	ltask_value[i] = m->last_tasks[i].value();
      }
    }
  }
}


Solution::~Solution() {
  delete [] earlybool_value;
  delete [] latebool_value;
  delete [] task_min;
  delete [] task_max;
  delete [] ltask_value;
  delete [] disjunct_value;
  delete [] search_value;
  delete [] all_value;
}

void Solution::guide_search() {
  
//   std::cout << "guide search with ";
//   print(std::cout);
//   std::cout << std::endl;

//   //solver->setGuidedOrdering(model->disjuncts, disjunct_value, "spl");
//   if( ((SchedulingSolver*)solver)->params->Type == "now2" )
//     solver->setGuidedBoundsOrdering(model->SearchVars, search_value);
//   else
    
    solver->setGuidedOrdering(model->SearchVars, search_value);

  if(model->data->hasJobDueDate()) {
    solver->setGuidedOrdering(model->earlybool, earlybool_value);
    solver->setGuidedOrdering(model->latebool, latebool_value);
    solver->setGuidedOrdering(model->last_tasks, ltask_value, "nbd");
  }
}

void Solution::guide_search_bounds() {
  
//   std::cout << "guide bound search with ";
//   print(std::cout);
//   std::cout << std::endl;

//   //solver->setGuidedOrdering(model->disjuncts, disjunct_value, "spl");
//   if( ((SchedulingSolver*)solver)->params->Type == "now2" )
//     solver->setGuidedBoundsOrdering(model->SearchVars, search_value);
//   else
    
    solver->setGuidedOrdering(model->SearchVars, search_value);
  if(model->data->hasJobDueDate()) {
    solver->setGuidedOrdering(model->earlybool, earlybool_value);
    solver->setGuidedOrdering(model->latebool, latebool_value);
    solver->setGuidedOrdering(model->last_tasks, ltask_value, "nbd");
  }
}

double Solution::distance(Solution* s) {
  
  int i, //n = model->disjuncts.size();
    n = model->SearchVars.size();
  double dist = 0.0;

  for(i=0; i<n; ++i) {
    //dist += (double)(disjunct_value[i] != s->disjunct_value[i]);
    dist += (double)(search_value[i] != s->search_value[i]);
  } 
  if(model->data->hasJobDueDate()) {
    n = model->last_tasks.size();
    for(i=0; i<n; ++i) {dist += 
	((ltask_value[i] - s->ltask_value[i])*(ltask_value[i] - s->ltask_value[i]));
    }
  }
  
  return sqrt(dist);
}

std::ostream& Solution::print(std::ostream& os, std::string type) {
  int i, //n = model->disjuncts.size();
    n = model->SearchVars.size();

//   for(i=0; i<n; ++i) {
//     os << disjunct_value[i] ;
//   } 
//   os << std::endl;

  if(type == "fsp") {
    int j, m = model->data->nJobs(), k = 0;
    int *rank = new int[m];
    int *order = new int[m];
    std::fill(rank, rank+m, 0);

    for(i=0; i<m-1; ++i) {
      for(j=i+1; j<m; ++j) {
	//if(!disjunct_value[k++]) ++rank[i];
	if(!disjunct_value[k++]) ++rank[i];
	else ++rank[j];
      }
    }
    
    for(i=0; i<m; ++i) {
      order[rank[i]] = i;
    }

    os << "c ";

    for(i=0; i<m; ++i) {
      os << (order[i]+1) << " ";
    }

    //os << std::endl;

    delete [] rank;
    delete [] order;

  } else {

    os << min_objective_value << " " << max_objective_value << std::endl;
    n = model->tasks.size();
    for(i=0; i<n; ++i) {
      os << i << " " 
	 << model->data->getDuration(i) << " "
	 << model->data->nMachines(i) << " ";
      for(int j=0; j<model->data->nMachines(i); ++j)
	os << model->data->getMachine(i,j) << " ";
      os << model->data->nJobs(i) << " ";
      for(int j=0; j<model->data->nJobs(i); ++j)
	os << model->data->getJob(i,j) << " ";
      os << task_min[i] << " " << task_max[i] << std::endl ;
    }
    os << std::endl;


    // for(i=0; i<n; ++i) {
    //   os << std::setw(3) << model->SearchVars[i].getVariable()->id ;
    // }
    // os << std::endl;
    // for(i=0; i<n; ++i) {
    //   os << std::setw(3) << search_value[i] ;
    // } 
    // os << std::endl;
  }
  return os;
}


SchedulingSolver::SchedulingSolver(SchedulingModel* m, 
				   ParameterList* p,
				   StatisticList* s) 
  : Solver(*m,m->SearchVars), stats(s) 
{ 
  model = m; 
  params = p;
  stats = s;
  stats->solver = this;

  params->initialise(this);

  stats->lower_bound = model->get_lb();//lower_bound;
  stats->upper_bound = model->get_ub();//upper_bound;

  params->getSkew();
  
  nogoods = NULL;
  pool = new SolutionPool();

  addHeuristic( params->Heuristic, params->Randomized, params->IValue, params->Hlimit );
}

std::ostream& SchedulingSolver::print_weights(std::ostream& os) {
  int i, n=variables.size;
  for(i=0; i<n; ++i) {
    sequence[i]->print(os);
    os << " " << sequence[i]->weight;
    os << " " << std::setw(4) << heuristic->get_value(sequence[i]) << std::endl;
  }
  //os << std::endl;
  n=constraints.size;
  for(i=0; i<n; ++i) {
    constraints[i]->print(os);
    os << " " << std::setw(4) << constraints[i]->weight << std::endl;
  }
  //os << std::endl;
  return os;
}

void SchedulingSolver::decay_weights(const double decay) {
  int i, w, n=variables.size;
  for(i=0; i<n; ++i) {
    w = (int)(((double)(variables[i]->weight))*decay);
    if(w > variables[i]->degree) variables[i]->weight = w;
    else variables[i]->weight = variables[i]->degree;
  }
  n = constraints.size;
  for(i=0; i<n; ++i) {
    if(constraints[i]->arity < variables.size/2) {
      w = (int)(((double)(constraints[i]->weight))*decay);
      if(w > 1) constraints[i]->weight = w;
      else constraints[i]->weight = 1;
    }
  }
}

void SchedulingSolver::jtl_presolve() 
{
  int objective = stats->upper_bound, new_objective;
  
  setVerbosity(params->Verbose);

  setNodeLimit(50000);

  addHeuristic("osp-t", 1, "anti", 1);
  JOB Heuristic( model, false );
  add( Heuristic );

  model->set_objective(stats->upper_bound);
  addObjective();

  if(params->Verbose>=0) {  
  std::cout << "c =============[ greedy step for jtl ]=============" << std::endl;
  std::cout << std::left << std::setw(30) << "c node cutoff" << ":"  
	    << std::right << std::setw(20) << 50000 << std::endl;
  //for(int iteration=0; iteration<params->InitBound; ++iteration) {
  }

  ((JobByJob*)heuristic)->shuffle();
  solve_and_restart();
  //solve();
  
  if( status == SAT ) {
    new_objective = model->get_objective();
    
    if(objective>new_objective) {
      pool->add(new Solution(model, this, stats->upper_bound));
      objective = new_objective;
     
      if(params->Verbose>=0) {      
	//std::cout << std::left << std::setw(30) << "c iteration / objective" << ":" << std::right << std::setw(11) << (iteration+1) << " / " << std::setw(6) << objective << ")" << std::endl;
	std::cout << std::left << std::setw(30) << "c objective" << ":" << std::right << std::setw(20) << objective << std::endl;
      }
    }
    
  } else {
    if(params->Verbose>=0) {
      std::cout << std::left << std::setw(30) << "c jtl presolve " << ":" << std::right << std::setw(20) << "failed" << std::endl;
    }
    params->InitBound = 0;
  }
  
  reset(true);
  decay_weights(params->Decay);
  
  
  //exit(1);
  //}
  
  if(params->Verbose>=0) {
    std::cout << "c ===============[ end greedy step ]===============" << std::endl;   
    std::cout << std::endl;
  }

  //resetBacktrackLimit();
  resetNodeLimit();
  removeObjective();
  
  if(pool->size()) {
    addHeuristic( params->Heuristic, params->Randomized, params->DValue, params->Hlimit );
  } else {
    addHeuristic( params->Heuristic, params->Randomized, params->IValue, params->Hlimit );
  }
  
  stats->upper_bound = objective;
  
}

void SchedulingSolver::old_jtl_presolve() 
{
  int objective = stats->upper_bound, new_objective;
  
  setVerbosity(params->Verbose);

  setBacktrackLimit(1000);

  addHeuristic("osp-t", 1, "anti", 1);
  JOB Heuristic( model, true );
  add( Heuristic );

 //  model->set_objective(stats->upper_bound);
//   addObjective();
  
  std::cout << "c =============[ greedy step for jtl ]=============" << std::endl;
  std::cout << std::left << std::setw(30) << "c backtrack cutoff" << ":"  
	    << std::right << std::setw(20) << 1000 << std::endl;
  for(int iteration=0; iteration<params->InitBound; ++iteration) {
    
    ((JobByJob*)heuristic)->shuffle();
    //solve_and_restart();
    solve();
    
    if( status == SAT ) {
      new_objective = model->get_objective();
      
      if(objective>new_objective) {
	objective = new_objective;
	pool->add(new Solution(model, this));
	
	std::cout << std::left << std::setw(30) << "c iteration / objective" << ":" << std::right << std::setw(11) << (iteration+1) << " / " << std::setw(6) << objective << ")" << std::endl;
	//std::cout << std::left << std::setw(30) << "c objective" << ":" << std::right << std::setw(20) << objective << std::endl;
      }
      
    } else {
      std::cout << std::left << std::setw(30) << "c jtl presolve " << ":" << std::right << std::setw(20) << "failed" << std::endl;
      params->InitBound = 0;
    }
    
    reset(true);
    decay_weights(params->Decay);

    
    //exit(1);
  }
   
  std::cout << "c ===============[ end greedy step ]===============" << std::endl;   
  std::cout << std::endl;
  
  resetBacktrackLimit();
  //resetNodeLimit();
  //removeObjective();
  
  if(pool->size()) {
    addHeuristic( params->Heuristic, params->Randomized, params->DValue, params->Hlimit );
  } else {
    addHeuristic( params->Heuristic, params->Randomized, params->IValue, params->Hlimit );
  }
  
  stats->upper_bound = objective;
  
}


bool SchedulingSolver::probe_ub() 
{
  bool ret_value = false;
  int objective = stats->upper_bound;
  
  setVerbosity(params->Verbose);

  addHeuristic( params->Heuristic, params->Randomized, params->IValue, params->Hlimit );

  if(params->Verbose>=0) {
    std::cout << "c =============[ initial probing step ]============" << std::endl;
    std::cout << std::left << std::setw(30) << "c propag cutoff" << ":"  
	      << std::right << std::setw(20) << (params->NodeCutoff/100) << std::endl;
    setPropagsLimit(params->NodeCutoff/100);
  }

  status = model->set_objective(objective);

  if(status == UNKNOWN) {
    solve_and_restart(params->PolicyRestart, params->Base, params->Factor);
  }

  if( status == SAT ) {
    ret_value = true;
    objective = model->get_objective();
    stats->normalized_objective = model->get_normalized_objective();
    pool->add(new Solution(model, this));
    stats->upper_bound = objective;

    if(params->Verbose>=0) {	
      std::cout << std::left << std::setw(30) << "c solutions's objective" << ":" << std::right << std::setw(20) << objective << std::endl;
    }

  } else if( status == UNSAT ) {

    stats->lower_bound = objective+1;
    
  } else {
    if(params->Verbose>=0) {
      std::cout << std::left << std::setw(30) << "c probing " << ":" << std::right << std::setw(20) << "failed" << std::endl;
    }
  }

  reset(true);
  decay_weights(params->Decay);

  if(pool->size()) {
    addHeuristic( params->Heuristic, params->Randomized, params->DValue, params->Hlimit );
  } 

  stats->upper_bound = objective;

  resetPropagsLimit();
  resetFailureLimit();

  if(params->Verbose>=0) {
    std::cout << "c ===============[ end probing step ]==============" << std::endl;
  }

  return ret_value;
}

void SchedulingSolver::dichotomic_search()
{

  presolve();
  
  if(!probe_ub()) {
    if(params->Presolve == "jtl_old") {
      old_jtl_presolve();
    } else if(params->Presolve == "jtl") {
      jtl_presolve();
    }
  }
  
  int iteration = 1;
  
  int minfsble = stats->lower_bound;
  int maxfsble = stats->upper_bound;

  int objective = -1;
  int new_objective = -1;
  int ngd_stamp = 0;
  int lit_stamp = 0;
    
  
  setVerbosity(params->Verbose);
  setTimeLimit(params->Cutoff);
  if(params->Randomized > 0)
    setRandomized(params->Randomized);
  else if(params->Randomized < 0) {
    randomizeSequence();
  }

  if(level < init_level) presolve();

  nogoods = NULL;
  if(params->Rngd) nogoods = setRestartGenNogood();

  ////////// dichotomic search ///////////////
  if(status == UNKNOWN) {
    while( minfsble<maxfsble && 
	   iteration<params->Dichotomy
	   ) {
      
      double remaining_time = params->Optimise - stats->get_total_time();
      
      if(remaining_time < (2*params->NodeBase)) break;

      objective = (int)(floor(((double)minfsble + (double)maxfsble)/2));

      if(params->Verbose>=0) {
	std::cout << "c ============[ start dichotomic step ]============" << std::endl;
      }

      setPropagsLimit(params->NodeCutoff);

      if(params->Verbose>=0) {
	std::cout << std::left << std::setw(30) << "c current dichotomic range" << ":" 
		  << std::right << std::setw(6) << " " << std::setw(5) << minfsble 
		  << " to " << std::setw(5) << maxfsble << std::endl;
	std::cout << std::left << std::setw(30) << "c target objective" << ":"  
		  << std::right << std::setw(20) << objective << std::endl;
      }

      save();

      status = model->set_objective(objective);
      if(pool->size()) {
	if(params->DValue == "guided") {
	  
	  //std::cout << params->Type << std::endl;

// 	  if(params->Type == "now2") 
// 	    pool->getBestSolution()->guide_search_bounds();
// 	  else
	    pool->getBestSolution()->guide_search();
	}
      }
      
      if(nogoods) {
	ngd_stamp = (params->Rngd>1 ? nogoods->base->nogood.size : 0);
	lit_stamp = sUnaryCons.size;
      }

      if(status == UNKNOWN) {
	solve_and_restart(params->PolicyRestart, params->Base, params->Factor);
      }
      
      if( status == SAT ) {
	new_objective = model->get_objective();

	stats->normalized_objective = model->get_normalized_objective();
	
	maxfsble = new_objective;
	pool->add(new Solution(model, this, objective));
	
	
	
	if(nogoods) {
	  for(int i=ngd_stamp; i<nogoods->base->nogood.size; ++i)
	    stats->avg_nogood_size += (double)(nogoods->base->nogood[i]->size);
	  if(params->Rngd>1) stats->num_nogoods = nogoods->base->nogood.size;
	  else stats->num_nogoods += nogoods->base->nogood.size;
	}

	if(params->Verbose>=0) {
	  std::cout << std::left << std::setw(30) << "c solutions's objective" << ":" << std::right << std::setw(20) << new_objective << std::endl;
	}
	//pool->getBestSolution()->print(std::cout);

      } else {

	new_objective = objective;
	minfsble = objective+1;

	if(nogoods) {
	  nogoods->forget(ngd_stamp);
	  nogoods->reinit();
	}
      }

      stats->add_info(new_objective, DICHO);

      if(params->Verbose>=0) {
	printStatistics(std::cout, ((params->Verbose ? RUNTIME : 0) + ((params->Verbose || status != UNKNOWN)  ? BTS + PPGS : 0) + OUTCOME) );
     } 

      reset(true);
      decay_weights(params->Decay);

      if(pool->size() && (params->DValue != params->IValue)) {
	addHeuristic( params->Heuristic, params->Randomized, params->DValue, params->Hlimit );
      }
      
      if(params->Verbose>=0) {
	std::cout << "c =============[ end dichotomic step ]=============" << std::endl;
      }

      ++iteration;
    } 
  } else if( status == SAT ) {
    std::cout << "c Solved during preprocessing!" << std::endl;

  } else if( status == UNSAT ) {
    std::cout << "c Found inconsistent during preprocessing!" << std::endl;

  }

  if(params->Verbose>=0) {      
    std::cout << std::endl;
  }
}


void SchedulingSolver::all_solutions_search()
{

  SolutionPool *second_pool = new SolutionPool();

  presolve();
  
  Vector< Decision > literals;
  int iteration = 0;
  int minfsble = stats->lower_bound;
  int maxfsble = stats->upper_bound;

  setVerbosity(params->Verbose);
  setTimeLimit(params->Cutoff);
  if(params->Randomized > 0)
    setRandomized(params->Randomized);
  else if(params->Randomized < 0) {
    randomizeSequence();
  }

  nogoods = setRestartGenNogood();

  while(true) {

    int ds_iteration = 1;

    int objective = -1;
    int new_objective = -1;
    int ngd_stamp = 0;
    int lit_stamp = 0;
    
    bool solution_found = false;

    if(params->Verbose>=0)
      std::cout << "c ============[ start dichotomic step ]============" << std::endl;

    ////////// dichotomic search ///////////////
    if(status == UNKNOWN) {
      while( minfsble<maxfsble && 
	     ds_iteration<params->Dichotomy
	     ) {
	
	double remaining_time = params->Optimise - stats->get_total_time();
	
	if(remaining_time < (2*params->NodeBase)) break;

	objective = (int)(floor(((double)minfsble + (double)maxfsble)/2));

	setPropagsLimit(params->NodeCutoff);

	if(params->Verbose>=0) {
	  std::cout << std::left << std::setw(30) << "c current dichotomic range" << ":" 
		  << std::right << std::setw(6) << " " << std::setw(5) << minfsble 
		    << " to " << std::setw(5) << maxfsble << std::endl;
	  std::cout << std::left << std::setw(30) << "c target objective" << ":"  
		    << std::right << std::setw(20) << objective << std::endl;
	}

	save();
	
	status = model->set_objective(objective);
	if(pool->size()) {
	  if(params->DValue == "guided") {

	    //std::cout << params->Type << std::endl;

// 	    if(params->Type == "now2")
// 	      pool->getBestSolution()->guide_search_bounds();
// 	    else
	      pool->getBestSolution()->guide_search();
	  }
	}
	
	ngd_stamp = (params->Rngd>1 ? nogoods->base->nogood.size : 0);
	lit_stamp = sUnaryCons.size;
	
	if(status == UNKNOWN) {
	  
// 	  std::cout << "nogood base:" << std::endl;
// 	  nogoods->print(std::cout);

	  for(int i=0; i<literals.size; ++i) {
	    sUnaryCons.push(literals[i]);
	  }

	  solve_and_restart(params->PolicyRestart, params->Base, params->Factor);
	}
	
	nogoods->forget(ngd_stamp);
	nogoods->reinit();
	
	if( status == SAT ) {

	  //std::cout << "c makespan " << model->C_max.value() << std::endl;

	  solution_found = true;
	  new_objective = model->get_objective();
	  
	  stats->normalized_objective = model->get_normalized_objective();
	  
	  maxfsble = new_objective;
	  pool->add(new Solution(model, this));
	  
	  if(params->Verbose>=0) {
	    std::cout << std::left << std::setw(30) << "c solutions's objective" << ":" << std::right << std::setw(20) << new_objective << std::endl;
	  }

	} else {
	  
	  new_objective = objective;
	  minfsble = objective+1;
	  
	}

	stats->add_info(new_objective, DICHO);
	
	//printStatistics(std::cout, ((params->Verbose ? RUNTIME : 0) + ((params->Verbose || status != UNKNOWN)  ? BTS + PPGS : 0) + OUTCOME) );
	
	
	reset(true);
	decay_weights(params->Decay);

	if(pool->size() && (params->DValue != params->IValue)) {
	  addHeuristic( params->Heuristic, params->Randomized, params->DValue, params->Hlimit );
	}
      
      
	++ds_iteration;
      } 
    } else if( status == SAT ) {
      std::cout << "c Solved during preprocessing!" << std::endl;
      
    } else if( status == UNSAT ) {
      std::cout << "c Found inconsistent during preprocessing!" << std::endl;
      
    }
    if(params->Verbose>=0) {
      std::cout << "c =============[ end dichotomic step ]=============" << std::endl;
    }

    //std::cout << status << " =?= " << UNKNOWN << std::endl;

    if(!solution_found) break;

    // std::cout << maxfsble << ": ";
    //  pool->getBestSolution()->print(std::cout, "fsp");
    //  std::cout << " ";
    // pool->getBestSolution()->print(std::cout, "jsp");
    // std::cout << std::endl;
    // // add the nogood corresponding to that solution

    second_pool->add(pool->getBestSolution());

    //std::cout << maxfsble << std::endl;

    if(maxfsble>0) {

      //int choice = params->NgdType;
      //if(choice == 2) choice = (maxfsble<=(model->disjuncts.size()/4));
      
      //for(int k=0; k<2; ++k) {
	Vector< Decision > learnt;
	int k=1;
	for(int i=0; i<model->disjuncts.size()/2; ++i) {
	  if(pool->getBestSolution()->disjunct_value[i] == k) {
	    Decision d('n', k, model->disjuncts[i].getVariable());
	    learnt.push(d);
	    // 	  d.print(std::cout);
	    // 	  std::cout << " ";
	  }
	}
	//std::cout << std::endl;
	if(learnt.size>1)
	  nogoods->base->add(learnt);
	else
	  literals.push(learnt[0]);
	//}

      //learnt[0].propagate();
      stats->upper_bound = maxfsble = model->disjuncts.size();

      //exit(1);
    }

    ++iteration;

    if(maxfsble == 0) break;
    //if(iteration>2)
    //break;
  }
      
  stats->num_solutions = iteration;
  std::cout << "d ITERATIONS  " << iteration << std::endl;


//   for(int i=0; i<second_pool->size(); ++i) {
//     for(int j=i+1; j<second_pool->size(); ++j) {
//       std::cout << "dist=" << (*second_pool)[i]->distance((*second_pool)[j]) << " (";
//       (*second_pool)[i]->print(std::cout, "fsp");
//       std::cout << " / " ;
//       (*second_pool)[j]->print(std::cout, "fsp");
//       std::cout << ")" << std::endl ;

//     }
//   }

}


void SolutionGuidedSearch::execute() 
  { 

    //std::cout << "HERE" << std::endl;

    SchedulingSolver *ss = (SchedulingSolver*)solver;
    pool->add(new Solution(ss->model, ss));
    if(pool->size()) pool->getBestSolution()->guide_search();
    StoreStats::execute();
  }


int SchedulingSolver::virtual_iterative_dfs()
    {
      SimpleUnaryConstraint last_decision;
      
      while( status == UNKNOWN ) {
	
	if( filtering() ) {
	  if(verbosity>2) {
	    VariableInt *t;
	    for(int i=0; i<model->data->nJobs(); ++i) {
	      std::cout << "j" << std::left << std::setw(2) << i << std::right ; //<< " ";
	      //std::cout << " (" << model->data->nTasksInJob(i) << ") ";
	      //std::cout.flush();

	      //if(params->Type != "now" && params->Type != "now2") {
	      if(model->tasks.size() > model->data->nJobs()) {
		for(int j=0; j<model->data->nTasksInJob(i); ++j) {
		  t = model->tasks[model->data->getJobTask(i,j)].getVariable();
		  std::cout << " " << t->id+1 << "[" << t->min() << ".." << t->max() << "]" ;
		}
	      } else {
		t = model->tasks[i].getVariable();
		std::cout << " " << t->id+1 << "[" << t->min() << ".." << t->max() << "]" ;
	      }
	      std::cout << std::endl;
	    }
	  }

	  if( future == empty ) {
	    solutionFound(init_level);
	  } else {

	    newNode();
	   
	    if(verbosity>2) {
	      VariableInt *d = decision.back();
	      MistralNode<Constraint*> *nd = d->constraintsOnValue();
	      
	      while( nextNode(nd) ) 
		if(nd->elt->arity == 3) {
		  PredicateDisjunctive *p = (PredicateDisjunctive *)(nd->elt);
		  //std::cout << "HERE: " << p;
		  std::cout << "c";
		  for(int k=0; k<=level; ++k) std::cout << " ";
		  p->print(std::cout);
		  std::cout << std::endl;
		}
	    }

	  }
	} else {
	  if( level <= init_level ) {
	    
#ifdef _DEBUGSEARCH
	    if(verbosity > 2) {
	      std::cout << "c UNSAT!" << std::endl;
	    }
#endif
	    
	    status = UNSAT;
	  } else if( limitsExpired() ) {
	    
#ifdef _DEBUGSEARCH
	    if(verbosity > 2) {
	      std::cout << "c";
	      for(int k=0; k<=level; ++k) std::cout << " ";
	      SimpleUnaryConstraint d = branching_decision[level];
	      d.revert();
	      d.print(std::cout);
	      std::cout << " (limit expired at level " << level << ")" << std::endl;
	    }
#endif
	    
	    status = LIMITOUT;
	  } else {
	    last_decision = branching_decision[level];

#ifdef _DEBUGSEARCH
	    if(verbosity > 2) {
	      if( level > backtrackLevel+1 ) {
		std::cout << "c";
		for(int k=0; k<=level; ++k) std::cout << " ";
		std::cout << " backjump to level " << backtrackLevel << std::endl;
	      }
	    }
#endif
	    
	    backtrackTo( backtrackLevel );	    
	    last_decision.deduce();

#ifdef _DEBUGSEARCH
	    if(verbosity > 2) {
	      std::cout << "c";
	      for(int k=0; k<=level; ++k) std::cout << " ";
	      last_decision.print( std::cout );
	      std::cout << std::endl;
	    }
#endif

	  }
	}
      }
      return status;
    }

void SchedulingSolver::branch_and_bound()
{
  int ngd_stamp = 0;
  int lit_stamp = 0;
  //resetNodeLimit();
  resetPropagsLimit();

  save();
  model->set_objective(stats->upper_bound-1);
  addObjective();

  setVerbosity(params->Verbose);
  //setRandomSeed( params->Seed );

  double time_limit = (params->Optimise - stats->get_total_time());

  if(time_limit > 0) {
    setTimeLimit( time_limit ); 
    addHeuristic( params->Heuristic, params->Randomized, params->Value, params->Hlimit );
    if(params->Value == "guided") {
      function = new SolutionGuidedSearch( this, pool, stats );
      if(pool->size()) pool->getBestSolution()->guide_search();
    } else {
      function = new StoreStats( this, stats );
    }
    
    if(params->Verbose>=0) {
      std::cout << "c ============[ start branch & bound ]=============" << std::endl;
      std::cout << std::left << std::setw(30) << "c current range" << ":" 
		<< std::right << std::setw(6) << " " << std::setw(5) << stats->lower_bound 
		<< " to " << std::setw(5) << goal->upper_bound << std::endl;
      std::cout << std::left << std::setw(30) << "c run for " << ":"
		<< std::right << std::setw(19) << (time_limit) << "s" << std::endl;
    }

    if(status == UNKNOWN) {

      if(nogoods) {
	ngd_stamp = nogoods->base->nogood.size;
	lit_stamp = sUnaryCons.size;
      }

      solve_and_restart(params->PolicyRestart, params->Base, params->Factor);

      if(nogoods) {
	for(int i=ngd_stamp; i<nogoods->base->nogood.size; ++i)
	  stats->avg_nogood_size += (double)(nogoods->base->nogood[i]->size);
	if(params->Rngd>1) stats->num_nogoods = nogoods->base->nogood.size;
	else stats->num_nogoods += nogoods->base->nogood.size;
      }
    }

    stats->add_info(goal->upper_bound, BNB);

    if(params->Verbose>=0) {    
      printStatistics(std::cout, ((params->Verbose ? RUNTIME : 0) + (params->Verbose ? BTS + PPGS : 0) + OUTCOME) );
    }
    
    reset(true);

    if(params->Verbose>=0) {
      std::cout << "c =============[ end branch & bound ]==============" << std::endl;
    }
  }
}

void StoreStats::execute()
{ 
  stats->normalized_objective = ss->model->get_normalized_objective();
}  

Vector<VariableInt*> unstable;

void SchedulingSolver::extract_stable(List& neighbors, 
				      Vector<VariableInt*>& stable)
{
  neighbors.clear();
  neighbors.random_fill(params->Neighbor);

  stable.clear();
  unstable.clear();
  int i, j, k=0, n=neighbors.capacity, m=0;
  
  for(i=0; i<n; ++i) if(!neighbors.member(i)) {
      m = model->data->nTasksInMachine(i);
      for(j=0; j<(m*(m-1)/2); ++j) {
	stable.push(model->disjuncts[k++].getVariable());
      }
    } else {
      k += (m*(m-1)/2);
//       m = model->data->nTasksInMachine(i);
//       for(j=0; j<(m*(m-1)/2); ++j) {
// 	unstable.push(model->disjuncts[k++].getVariable());
//       }
    }
}

void SchedulingSolver::large_neighborhood_search() {
   Vector< VariableInt* > stable;
   List neighbors;
   neighbors.init(0, model->data->nMachines()-1);
   int objective = INFTY;

   save();
   model->set_objective(stats->upper_bound);
   addObjective();

   while(true) {
     // select a subset of variables that will stay stable
     extract_stable(neighbors,stable);
    
     // optimise the remaining part
     objective = stats->upper_bound;
     repair(pool->getBestSolution(), stable);
     if(stats->upper_bound <= objective)
       pool->add(new Solution(model, this));
     reset(true);
   }
}


void SchedulingSolver::repair(Solution *sol, Vector<VariableInt*>& stable)  
{
  save();
  for(int i=0; i<stable.size; ++i) {
    Decision branch('e', (*sol)[stable[i]], stable[i]);
    branch.propagate();
  }
  solve_and_restart() ;
  stats->add_info(goal->upper_bound, LNS);
}


void SchedulingSolver::print_solution(std::ostream& os, std::string type)
{
  pool->getBestSolution()->print(os, type);
  os << std::endl;
}
