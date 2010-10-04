 

#include "new_scheduler.hpp"


using namespace MistralScheduler;


#define INFTY 0xffffff
#define BIG 0xffff


StatisticList::StatisticList() {
  best_solution_index     = 0;
  branch_and_bound_index  = 0;

  lower_bound         = INFTY;
  upper_bound         = 0;

  num_nogoods         = 0;
  num_solutions       = 0;
  
  avg_distance        = 0;
  min_distance        = INFTY;
  max_distance        = 0;

  /*
  total_time          = 0.0;
  total_nodes         = 0;
  total_backtracks    = 0;
  total_fails         = 0;
  total_propags       = 0;
  */
}

StatisticList::~StatisticList() {}

int StatisticList::get_total_time() {
  double total_time  = 0.0;
  for(unsigned int i=0; i<time.size(); ++i)
    total_time += time[i];
  return total_time;
}

void StatisticList::add_info(SchedulingSolver *s, const int objective) {

  DBG("Update statistics%s\n", "");

  time.push_back(s->ENDTIME);
  nodes.push_back(s->NODES);
  backtracks.push_back(s->BACKTRACKS);
  fails.push_back(s->FAILURES);
  propags.push_back(s->PROPAGS);

  outcome.push_back(s->status);

  /*
  total_time         += time.back();
  total_nodes        += nodes.back();
  total_backtracks   += backtracks.back();
  total_fails        += fails.back();
  total_propags      += propags.back();
  */

  if(outcome.back() == SAT) {
    ++num_solutions;
    best_solution_index = outcome.size()-1;
    upper_bound = objective;
  } else if(outcome.back() == UNSAT) {
    lower_bound = objective+1;
  }

}



std::ostream& StatisticList::print(std::ostream& os, 
				   const char* prefix,
				   const int start, 
				   const int end) {

  int k, i=start, j=outcome.size();
  if(end >= 0) j=end;

  double total_time  = 0.0;
  double opt_time    = 0.0;
  double proof_time  = 0.0;
  long unsigned int total_nodes         = 0;
  long unsigned int total_backtracks    = 0;
  long unsigned int total_fails         = 0;
  long unsigned int total_propags       = 0;

  for(k=i; k<j; ++k) {
    if(k<=best_solution_index) opt_time += time[k];
    total_time         += time[k];
    total_nodes        += nodes[k];
    total_backtracks   += backtracks[k];
    total_fails        += fails[k];
    total_propags      += propags[k];
  }
  proof_time = (total_time - opt_time);


  int plength = 0;
  while(prefix[plength] != '\0') ++plength;
  
  os << "c =================[ statistics ]==================" << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "LOWERBOUND "    << std::right << std::setw(21) << lower_bound << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "OBJECTIVE "     << std::right << std::setw(21) << upper_bound << std::endl 
     << "d " << prefix << std::left << std::setw(28-plength)  << "RUNTIME "       << std::right << std::setw(21) << total_time << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "OPTTIME "       << std::right << std::setw(21) << opt_time << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "PROOFTIME "       << std::right << std::setw(21) << proof_time << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "NODES "         << std::right << std::setw(21) << total_nodes << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "BACKTRACKS "    << std::right << std::setw(21) << total_backtracks << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "FAILS "         << std::right << std::setw(21) << total_fails << std::endl
     << "d " << prefix << std::left << std::setw(28-plength)  << "PROPAGS "       << std::right << std::setw(21) << total_propags << std::endl
    //<< "d " << prefix << std::left << std::setw(28-plength)  << "RESTARTS "      << std::right << std::setw(21) << total_restarts << std::endl
    
     << "d " << prefix << std::left << std::setw(28-plength)  << "SOLUTIONS "     << std::right << std::setw(21) << num_solutions << std::endl
    //<< std::left << std::setw(28) << "d AVGDISTANCE "     << std::right << std::setw(21) << avg_distance << std::endl
    //<< std::left << std::setw(28) << "d MINDISTANCE "     << std::right << std::setw(21) << min_distance << std::endl
    //<< std::left << std::setw(28) << "d MAXDISTANCE "     << std::right << std::setw(21) << max_distance << std::endl;
     << "d " << prefix << std::left << std::setw(28-plength)  << "OPTIMAL " << std::right << std::setw(21) << (lower_bound == upper_bound) << std::endl
     << "c =================[ statistics ]==================" << std::endl
     << std::endl;

  return os;
}  


// void StatisticList::add_solution(SolutionPool *s) {
//   best_solution_index = 
// }


const char* ParameterList::int_ident[ParameterList::nia] = 
  {"-ub", "-lb", "-check", "-seed", "-cutoff", "-dichotomy", 
   "-base", "-randomized", "-verbose", "-optimise", "-nogood"};

const char* ParameterList::str_ident[ParameterList::nsa] = 
  {"-heuristic", "-restart", "-factor", "-decay", "-type", 
   "-value", "-dvalue", "-ivalue", "-skew", "-objective"};

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
  Cutoff      = 30;
  Dichotomy   = 128;
  Base        = 256;
  Randomized  = 1;

  Verbose     = 1;
  Optimise    = 3600;
  Rngd        = 1;

  Policy    = "geom";
  Factor    = 1.3;
  Decay     = 0.0;
  Value     = "guided";
  DValue    = "guided";
  IValue    = "promise";
  Skew      = -1.0;
  Objective = "makespan";

  Heuristic = "osp-t";
  PolicyRestart = GEOMETRIC;



  if(Type == "osp") {
    Cutoff      = 30;
    Objective = "makespan";
    Heuristic = "osp-b";
  } else if(Type == "sds") {
    Cutoff      = 30;
    Objective = "makespan";
    Heuristic = "osp-b";
  } else if(Type == "jtl") {
    Cutoff      = 30;
    Objective = "makespan";
    Heuristic = "osp-t";
  } else if(Type == "now") {
    Cutoff      = 30;
    Objective = "makespan";
    Heuristic = "osp-dw";
  } else if(Type == "jla") {
    Cutoff      = 300;
    Objective = "makespan";
    Heuristic = "osp-t";
  } else if(Type == "jsp") {
    Cutoff      = 300;
    Objective = "makespan";
    Heuristic = "osp-t";
  } else if(Type == "jet") {
    Cutoff      = 30;
    Objective = "tardiness";
    Heuristic = "osp-t";
  } else if(Type == "dyn") {
    Cutoff      = 30;
    Objective = "tardiness";
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

  if(strcmp(str_param[0],"nil")) Heuristic  = str_param[0];
  if(strcmp(str_param[1],"nil")) Policy     = str_param[1];
  if(strcmp(str_param[2],"nil")) Factor     = atof(str_param[2]);
  if(strcmp(str_param[3],"nil")) Decay      = atof(str_param[3]);
  if(strcmp(str_param[5],"nil")) Value      = str_param[5];
  if(strcmp(str_param[6],"nil")) DValue     = str_param[6];
  if(strcmp(str_param[7],"nil")) IValue     = str_param[7];
  if(strcmp(str_param[8],"nil")) Skew       = atof(str_param[8]);
  if(strcmp(str_param[9],"nil")) Objective  = str_param[9];


  if(Policy == "luby")
    PolicyRestart = LUBY;
  else if(Policy == "geom")
    PolicyRestart = GEOMETRIC;
  else
    PolicyRestart = NO;
}

std::ostream& ParameterList::print(std::ostream& os) {
  os << std::endl << "c =================[ parameters ]==================" << std::endl;
  os << std::left << std::setw(30) << "c data file " << ":" << std::right << std::setw(20) << data_file_name << std::endl;
  os << std::left << std::setw(30) << "c type " << ":" << std::right << std::setw(20) << Type << std::endl;
  os << std::left << std::setw(30) << "c seed " << ":" << std::right << std::setw(20) << Seed << std::endl;
  os << std::left << std::setw(30) << "c cutoff " << ":" << std::right << std::setw(20) << Cutoff << std::endl;
  os << std::left << std::setw(30) << "c dichotomy " << ":" << std::right << std::setw(20) << (Dichotomy ? "yes" : "no") << std::endl;
  os << std::left << std::setw(30) << "c restart policy " << ":" << std::right << std::setw(20) << Policy << std::endl;
  os << std::left << std::setw(30) << "c base " << ":" << std::right << std::setw(20) << Base << std::endl;
  os << std::left << std::setw(30) << "c factor " << ":" << std::right << std::setw(20) << Factor << std::endl;
  os << std::left << std::setw(30) << "c heuristic " << ":" << std::right << std::setw(20) << Heuristic << " (" << abs(Randomized) << ")" << std::endl;
  os << std::left << std::setw(30) << "c value ordering (init step) " << ":" << std::right << std::setw(20) << IValue << std::endl;
  os << std::left << std::setw(30) << "c value ordering (dichotomy) " << ":" << std::right << std::setw(20) << DValue << std::endl;
  os << std::left << std::setw(30) << "c value ordering (optim) " << ":" << std::right << std::setw(20) << Value << std::endl;
  os << "c =================[ parameters ]==================" << std::endl;
  
  return os;
}


Instance::Instance(const ParameterList& params) {

  DBG("Build instance %s\n", params.data_file);

  setup_time = NULL;
  time_lag[0] = NULL;
  time_lag[1] = NULL;
  jsp_duedate = NULL;
  
  if(params.Type == "osp") {
    osp_readData( params.data_file );
  } else if(params.Type == "sds") {
    sds_readData( params.data_file );
  } else if(params.Type == "jtl") {
    jtl_readData( params.data_file );
  } else if(params.Type == "now") {
    jtl_readData( params.data_file );
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
    dyn_readData( params.data_file );
  }

  // close();

}

Instance::~Instance() {
}

// void Instance::close() {
//   // remove unused machines:
//   for()
  
// }


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
  
  //std::cout << k << " " << i << " " << j << " -> ";

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

  //std::cout << k << " " << ri << " " << rj << " : " << setup_time[k][ri][rj] << std::endl;

  //std::cout << setup_time << std::endl;
  //std::cout << setup_time[k] << std::endl;
  //std::cout << setup_time[k][ri] << std::endl;
  //std::cout << setup_time[k][ri][rj] << std::endl;
  return setup_time[k][ri][rj];}

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
  int best_makespan = INFTY;
  if(!hasTimeLag()) {
    int current_job[nJobs()];
    int current_job_bound[nJobs()];
    int current_machine_bound[nMachines()];
    int current_machine_task[nMachines()];
    int random_jobs[nTasks()];
    int m[nMachines()];
    int k=0;

    for(int i=0; i<nJobs(); ++i)
      for(int j=0; j<nTasksInJob(i); ++j)
	random_jobs[k++] = i;

    int iter = iterations;
  
    while(iter--) {
      std::fill(current_job, current_job+nJobs(), 0);
      std::fill(current_job_bound, current_job_bound+nJobs(), 0);
      std::fill(current_machine_bound, current_machine_bound+nMachines(), 0);
      std::fill(current_machine_task, current_machine_task+nMachines(), -1);
  
      int makespan = 0;

      for(int i=0; i<nTasks(); ++i) {
	int j = randint(nTasks()-i);
	int k = random_jobs[i];
	random_jobs[i] = random_jobs[i+j];
	random_jobs[i+j] = k;
      }
    
      for(int i=0; i<nTasks(); ++i) {
	// pick the next job
	int j = random_jobs[i];
	// find out which task is that
	int t = getJobTask(j, current_job[j]);

	DBG("pick task t%d\n", t);

	// find out which machine is that
	for(int j=0; j<nMachines(t); ++j) {
	  m[j] = getMachine(t,j);
	  DBG("  -> uses machine%d\n", m[j]);
	}
      
	// find the current timepoint for this job
	int j_mkp = current_job_bound[j];

	// find the current timepoint for this machine
	int m_mkp = current_machine_bound[m[0]];
	DBG("m%d = %d\n", m[0], current_machine_bound[m[0]]);
	for(int j=1; j<nMachines(t); ++j)
	  if(m_mkp < current_machine_bound[m[j]]) {
	    m_mkp = current_machine_bound[m[j]];
	    DBG("m%d = %d\n", m[j], current_machine_bound[m[j]]);
	  }

	// get the start time for this task
	int timepoint = (j_mkp < m_mkp ? m_mkp : j_mkp);

	// check its release date
	if(getReleaseDate(t) > timepoint) timepoint = getReleaseDate(t);

	DBG("timepoint = %d\n", timepoint);

	// add setup time, if any
	if(hasSetupTime()) {
	  int setup = 0;
	  int setup_mj;
	  for(int j=0; j<nMachines(t); ++j) {
	    if(current_machine_task[m[j]] >= 0) {
	      setup_mj = getSetupTime(m[j], current_machine_task[m[j]], t);
	      if(setup < setup_mj) setup = setup_mj;
	      DBG("setup = %d\n", setup_mj);
	    }
	  }
	  timepoint += setup;
	}

	// get the finish time for this task
	timepoint += getDuration(t);

	// update machin and job bounds
	for(int j=0; j<nMachines(t); ++j) {
	  current_machine_bound[m[j]] = timepoint;
	  current_machine_task[m[j]] = t;
	}
	//current_machine_bound[m] = timepoint;
	current_job_bound[j] = timepoint;

	// get the final makespan
	if(makespan < timepoint) makespan = timepoint;

	++current_job[j];
      }
      if(best_makespan > makespan) best_makespan = makespan;
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
    for(i=0; i<nJobs; ++i) 
      for(j=0; j<nJobs; ++j) {
	setup_time[k][i][j] = family_matrix[1+family[i][k]][family[j][k]] ;
      }
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

void Instance::dyn_readData( const char* filename ) {

  DBG("Read (dyn)%s\n", "");

  int i, j, k, dur, mach, nJobs, nMachines;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> nJobs;
  infile >> nMachines;
	
  jsp_duedate = new int[nJobs];
  jsp_earlycost = new int[nJobs];
  jsp_latecost = new int[nJobs];
  

  double *float_weights = new double[nJobs];
  
  for(i=0; i<nJobs; ++i) {
    k = release_date.size();

    for(j=0; j<nMachines; ++j) {
      
      infile >> mach;
      infile >> dur;
      
      if(mach != -1 && dur != -1) addTask(dur, i, mach);
    }		



    infile >> jsp_duedate[i];
    infile >> j;
    infile >> jsp_latecost[i];
    infile >> float_weights[i];

    std::cout << (k) << " " << release_date.size() << std::endl;
    setReleaseDate(k, j);
    jsp_earlycost[i] = jsp_latecost[i];
  }
}



SchedulingModel::SchedulingModel(Instance& inst, const int max_makespan) : CSP() {

  setup(inst, max_makespan);

}

void SchedulingModel::setup(Instance& inst, const int max_makespan) {

  int i,j,k, lb, ub, ti, tj;
  
  lb_C_max = inst.getMakespanLowerBound();
  if(max_makespan < 0) ub_C_max = inst.getMakespanUpperBound(1);
  else ub_C_max = max_makespan;

  for(i=0; i<inst.nJobs(); ++i)
    ub_C_max += inst.getDuration(inst.getLastTaskofJob(i));

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

	disjuncts.add( Disjunctive( tasks[ti],
				    
				    inst.getDuration(ti)
				    +(inst.hasSetupTime() ? inst.getSetupTime(k,ti,tj) : 0),
				    
				    tasks[tj],
				    
				    inst.getDuration(tj)
				    +(inst.hasSetupTime() ? inst.getSetupTime(k,tj,ti) : 0) ) );
      }
  add( disjuncts );

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
    
    Variable x(0, sum_ub);
    L_sum = x;
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
}

SchedulingModel::~SchedulingModel() {}


void No_wait_Model::setup(Instance& inst, const int max_makespan) {

  //std::cout << "nowait setup" << std::endl;
  //exit(1);

  int i,j,k, lb, ub, ti=0, tj=0;
  std::vector<int> offset;
  //**offset=new int*[inst.nJobs()];
  for(i=0; i<inst.nJobs(); ++i) {
    //offset[i] = new int[inst.nTasksInJob(i)+1];
    //offset[i][0] = 0;
    //for(j=0; j<inst.nTasksInJob(i); ++j)
    //offset[i][j+1] = (offset[i][j] + inst.getDuration(inst.getJobTask(i,j)));

    offset.push_back(0);
    for(j=0; j<inst.nTasksInJob(i)-1; ++j)
      offset.push_back(offset.back() + inst.getDuration(inst.getJobTask(i,j)));
  }



  lb_C_max = inst.getMakespanLowerBound();
  if(max_makespan < 0) {
    ub_C_max = 0;
    for(i=0; i<inst.nJobs(); ++i) {

      //std::cout << offset.size() << " " << i << " " << ((i+1)*inst.nTasksInJob(i)-1) << std::endl;

      //ub_C_max += (offset[(i+1)*inst.nTasksInJob(i)-1] + inst.getDuration(inst.getLastTaskofJob(i)));
      ub_C_max += (offset[inst.getLastTaskofJob(i)] + inst.getDuration(inst.getLastTaskofJob(i)));
    }
  } else ub_C_max = max_makespan;


  //exit(1);

  //lb_L_sum = inst.getEarlinessTardinessLowerBound(ub_C_max);
  //ub_L_sum = inst.getEarlinessTardinessUpperBound(ub_C_max);

  data = &inst;

  // create one variable per task
  for(i=0; i<inst.nJobs(); ++i) {

    lb = inst.getReleaseDate(inst.getJobTask(i,0));
    ub = std::min(ub_C_max, inst.getDueDate(inst.getLastTaskofJob(i))) 
      //- (offset[(i+1)*inst.nTasksInJob(i)-1] + inst.getDuration(inst.getLastTaskofJob(i)));
      - (offset[inst.getLastTaskofJob(i)] + inst.getDuration(inst.getLastTaskofJob(i)));


    if(lb > ub) {
      std::cout << "INCONSISTENT" << std::endl;
      exit(1);
    }

    Variable t(lb, ub);
    tasks.add(t);
  }

  // mutual exclusion constraints
  for(k=0; k<inst.nMachines(); ++k) 
    for(i=0; i<inst.nTasksInMachine(k); ++i)
      for(j=i+1; j<inst.nTasksInMachine(k); ++j) {
	ti = inst.getMachineTask(k,i);
	tj = inst.getMachineTask(k,j);

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

  Variable x_cmax(lb_C_max, ub_C_max);
  C_max = x_cmax;
  for(i=0; i<inst.nJobs(); ++i) {
    ti = inst.getLastTaskofJob(i);
    k = inst.getDuration(ti) + offset[ti];
    add(Precedence(tasks[i], k, C_max));
  }
  for(int i=0; i<disjuncts.size(); ++i) SearchVars.add(disjuncts[i]);

//   for(i=0; i<inst.nJobs(); ++i) {
//     for(j=0; j<inst.nTasksInJob(i); ++j) {
//       ti = inst.getJobTask(i,j);
//       std::cout << offset[ti] << "+t" << ti << "+" << inst.getDuration(ti) << " ";
//     }
//     std::cout << std::endl;
//   }

//   //exit(1);

}


int L_sum_Model::get_lb() {
  return lb_L_sum;
}

int L_sum_Model::get_ub() {
  return ub_L_sum;
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
// void L_sum_Model::add_objective() {
//   add(Minimise(L_sum));
// }

VariableInt* C_max_Model::get_objective_var() {
  return C_max.getVariable();
}
// void C_max_Model::add_objective() {
//   add(Minimise(C_max));
// }

int L_sum_Model::set_objective(const int obj) {
  //return UNKNOWN;

  std::cout << (L_sum.getVariable()) << std::endl;
  L_sum.getVariable()->print(std::cout);
  std::cout << " <= " << obj << std::endl;

  return (L_sum.getVariable()->setMax(obj) ? UNKNOWN : UNSAT);
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
  //return -1;

  L_sum.getVariable()->print(std::cout);
  std::cout << std::endl;

  int total_cost = 0;
  for(int i=0; i<latebool.size(); ++i) {
    if(earlybool[i].getVariable()) earlybool[i].getVariable()->print(std::cout) ;
    else std::cout << (earlybool[i].value());
    std::cout << " ";
    if(latebool[i].getVariable()) latebool[i].getVariable()->print(std::cout);
    else std::cout << (latebool[i].value());
    std::cout << " ";
    if(last_tasks[i].getVariable()) last_tasks[i].getVariable()->print(std::cout);
    else std::cout << "-";
    std::cout << " " << data->getJobDueDate(i) << " " 
	      << data->getJobEarlyCost(i) << " " 
	      << data->getJobLateCost(i) << " ";

    int cost = 0;
    if(earlybool[i].value()) {
      cost = (data->getJobDueDate(i) - (last_tasks[i].value()+data->getDuration(data->getLastTaskofJob(i))))*data->getJobEarlyCost(i);
    } else if(latebool[i].value()) {
      cost = ((last_tasks[i].value()+data->getDuration(data->getLastTaskofJob(i))) - data->getJobDueDate(i))*data->getJobLateCost(i);
    }

    std::cout << cost << std::endl;

    total_cost += cost;
  }

  std::cout << total_cost << std::endl;

  return L_sum.value();
}
int C_max_Model::get_objective() {
//   int obj=0, aux, ti;
//   for(int i=0; i<data->nJobs(); ++i) {
//     ti = data->getLastTaskofJob(i);
//     aux = (tasks[ti].value() + data->getDuration(ti));
//     if(aux > obj) obj = aux;
//   }
//   return obj;
  return C_max.value();
}


Solution::Solution(SchedulingModel *m, SchedulingSolver *s) {
  model = m;
  solver = s;

  earlybool_value = NULL;
  latebool_value = NULL;
  task_min = NULL;
  task_max = NULL;
  ltask_value = NULL;
  disjunct_value = NULL;

  int i, n;

  n = m->disjuncts.size();
  if(n) {
    disjunct_value = new int[n];
    for(i=0; i<n; ++i) disjunct_value[i] = m->disjuncts[i].value();
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
    n = m->tasks.size();
    if(n) {
      task_min = new int[n];
      task_max = new int[n];
      for(i=0; i<n; ++i) {
	task_min[i] = m->tasks[i].min();
	task_max[i] = m->tasks[i].max();
      } 
    }

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
}

void Solution::guide_search() {
  solver->setGuidedOrdering(model->disjuncts, disjunct_value);
  if(model->data->hasJobDueDate()) {
    solver->setGuidedOrdering(model->earlybool, earlybool_value);
    solver->setGuidedOrdering(model->latebool, latebool_value);
    solver->setGuidedOrdering(model->last_tasks, ltask_value, "2nd");
  }
}

std::ostream& SchedulingSolver::print_weights(std::ostream& os) {
  int n=100; //length;
  for(int i=0; i<n; ++i) {
    os << " " << std::setw(4) << variables[i]->weight;
  }
  os << std::endl;
  return os;
}

void SchedulingSolver::decay_weights(const double decay) {
  int w, n=variables.size;
  for(int i=0; i<n; ++i) {
    w = (int)(((double)(variables[i]->weight))*decay);
    if(w > variables[i]->degree) variables[i]->weight = w;
    else variables[i]->weight = variables[i]->degree;
  }
}


void SchedulingSolver::dichotomic_search()
{
  
  int iteration = 0;
  
  int minfsble = lower_bound;
  int maxfsble = upper_bound;
  int objective = -1;
  int new_objective = -1;
  int Solved = true;
  
  presolve();
  
  setVerbosity(params.Verbose);
  setTimeLimit(params.Cutoff);
  
  
  ////////// dichotomic search ///////////////
  if(status == UNKNOWN) 
    while( minfsble < maxfsble && iteration < params.Dichotomy ) {
      setRandomSeed( params.Seed );
      
      objective = (int)(floor(((double)minfsble + (double)maxfsble)/2));
      
      std::cout << "c ============[ start dichotomic step ]============" << std::endl;
      std::cout << std::left << std::setw(30) << "c current dichotomic range" << ":" 
		<< std::right << std::setw(6) << " " << std::setw(5) << minfsble 
		<< " to " << std::setw(5) << maxfsble << std::endl;
      std::cout << std::left << std::setw(30) << "c target objective" << ":"  
		<< std::right << std::setw(20) << objective << std::endl;
      //std::cout << "c ===========[ start dichotomic search ]===========" << std::endl;
	  
      save();
      status = model->set_objective(objective);

      if(pool->size()) pool->getBestSolution()->guide_search();
      
      if(status == UNKNOWN) {
	solve_and_restart(params.PolicyRestart, params.Base, params.Factor, params.Decay);
      }
      
      if( status == SAT ) {
	
	new_objective = model->get_objective();
	maxfsble = new_objective;
	pool->add(new Solution(model, this));
	
	if(pool->size()) {
	  addHeuristic( params.Heuristic, params.Randomized, params.DValue );
	}
	
	std::cout << std::left << std::setw(30) << "c solutions's objective" << ":" << std::right << std::setw(20) << new_objective << std::endl;
	
      } else {
	
	new_objective = objective;
	minfsble = objective+1;
	if( status != UNSAT ) Solved = false;
	else lower_bound = objective;
      }
      
      stats.add_info(this, new_objective);
      
      printStatistics(std::cout, (RUNTIME + BTS + PPGS + OUTCOME) );
      
      reset(true);
      decay_weights(params.Decay);
      
      std::cout << "c =============[ end dichotomic step ]=============" << std::endl;
      
      ++iteration;
    } 
      
  std::cout << std::endl;
  
  upper_bound = maxfsble;

  Solved = (upper_bound == lower_bound);
  
}


void SchedulingSolver::branch_and_bound()
{
  int Solved = false;

  addObjective();

  setRandomSeed( params.Seed );
  setTimeLimit( params.Optimise - stats.get_total_time() ); 
  function = new SolutionGuidedSearch( this, pool );

  save();
  model->set_objective(upper_bound-1);
     
  std::cout << "c ============[ start branch & bound ]============" << std::endl;
  std::cout << std::left << std::setw(30) << "c current range" << ":" 
	    << std::right << std::setw(6) << " " << std::setw(5) << lower_bound 
	    << " to " << std::setw(5) << upper_bound << std::endl;

  if(pool->size()) pool->getBestSolution()->guide_search();
  
  if(status == UNKNOWN) {
    solve_and_restart(params.PolicyRestart, params.Base, params.Factor, params.Decay);
  }
  if( status == OPT ) {
    lower_bound = goal->upper_bound+1;
  }

  upper_bound = this->goal->upper_bound+1;
  stats.add_info(this, upper_bound);
      
  printStatistics(std::cout, (RUNTIME + BTS + PPGS + OUTCOME) );
      
  reset(true);
  std::cout << "c =============[ end branch & bound ]=============" << std::endl;
  
  Solved = (upper_bound == lower_bound);

}
  

int main( int argc, char** argv )
{
  ParameterList params(argc, argv);
  StatisticList stats;

  usrand( params.Seed );

  params.print(std::cout);

  Instance jsp(params);
  
  jsp.print(std::cout);

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
  
  solver.print(std::cout);

  if(solver.status == UNKNOWN) solver.dichotomic_search();
  if(solver.status == UNKNOWN) solver.branch_and_bound();

  stats.print(std::cout, "DS");
  
}
  




