 

#include "mistral_scheduler.hpp"
#include <math.h>
#include <assert.h>

using namespace Mistral;


//#define INFTY 0xffffff
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
	real_start_time = get_run_time();
}

// void StatisticList::stop() {
//   real_time = (get_run_time() - real_time);
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

	double runtime = (solver->statistics.end_time - solver->statistics.start_time);
	//std::cout << "ADDINFO RUNTIME=" << runtime << std::endl;

	time.push_back(runtime);
	soltime.push_back(runtime);
	nodes.push_back(solver->statistics.num_nodes);
	backtracks.push_back(solver->statistics.num_backtracks);
	fails.push_back(solver->statistics.num_failures);
	propags.push_back(solver->statistics.num_propagations);
	types.push_back(tp);

	//std::cout << outcome2str(solver->statistics.outcome) << std::endl;

	outcome.push_back(solver->statistics.outcome);


  

	if(tp == BNB) {
		if(outcome.back() == OPT) {
			lower_bound = upper_bound = objective;
			best_solution_index = outcome.size()-1;
		} else if(outcome.back() == UNSAT) {
			lower_bound = objective+1;
		}
	} else {
		if(outcome.back() == SAT || outcome.back() == OPT) {
			++num_solutions;
			best_solution_index = outcome.size()-1;
			upper_bound = objective;
			if(outcome.back() == OPT) lower_bound = objective;
		} else if((types.back() != LNS) && (outcome.back() == UNSAT)) {
			lower_bound = objective+1;
		}
	}

	//std::cout << outcome2str(outcome.back()) << " ==> [" << lower_bound << ".." << upper_bound << "]" << std::endl;

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
  
	os << "\n c +==============[ statistics ]===============+" << std::endl
		<< " d " << prefix << std::left << std::setw(25-plength)  << " LOWERBOUND "    << std::right << std::setw(19) << lower_bound << std::endl
			<< " d " << prefix << std::left << std::setw(25-plength)  << " UPPERBOUND "    << std::right << std::setw(19) << upper_bound << std::endl
				<< " d " << prefix << std::left << std::setw(25-plength)  << " OBJECTIVE "     << std::right << std::setw(19) << upper_bound << std::endl 
					<< " d " << prefix << std::left << std::setw(25-plength)  << " NORMOBJECTIVE " << std::right << std::setw(19) << normalized_objective << std::endl 
						<< " d " << prefix << std::left << std::setw(25-plength)  << " REALTIME "      << std::right << std::setw(19) << (get_run_time() - real_start_time) << std::endl
							<< " d " << prefix << std::left << std::setw(25-plength)  << " RUNTIME "       << std::right << std::setw(19) << total_time << std::endl
								<< " d " << prefix << std::left << std::setw(25-plength)  << " OPTTIME "       << std::right << std::setw(19) << opt_time << std::endl
									<< " d " << prefix << std::left << std::setw(25-plength)  << " PROOFTIME "     << std::right << std::setw(19) << proof_time << std::endl
										<< " d " << prefix << std::left << std::setw(25-plength)  << " UBTIME "        << std::right << std::setw(19) << ub_time << std::endl
											<< " d " << prefix << std::left << std::setw(25-plength)  << " LBTIME "        << std::right << std::setw(19) << lb_time << std::endl
												<< " d " << prefix << std::left << std::setw(25-plength)  << " LOSTTIME "      << std::right << std::setw(19) << lost_time << std::endl
													<< " d " << prefix << std::left << std::setw(25-plength)  << " NODES "         << std::right << std::setw(19) << total_nodes << std::endl
														<< " d " << prefix << std::left << std::setw(25-plength)  << " BACKTRACKS "    << std::right << std::setw(19) << total_backtracks << std::endl
															<< " d " << prefix << std::left << std::setw(25-plength)  << " FAILS "         << std::right << std::setw(19) << total_fails << std::endl
																<< " d " << prefix << std::left << std::setw(25-plength)  << " PROPAGS "       << std::right << std::setw(19) << total_propags << std::endl
																	<< " d " << prefix << std::left << std::setw(25-plength)  << " NOGOODS "       << std::right << std::setw(19) << num_nogoods << std::endl;
	if(num_nogoods)
		std::cout << " d " << prefix << std::left << std::setw(25-plength)  << " NOGOODSIZE "    << std::right << std::setw(19) << avg_nogood_size/(double)num_nogoods << std::endl;
	std::cout << " d " << prefix << std::left << std::setw(25-plength)  << " NODES/s "       << std::right << std::setw(19) ;
	if(total_time > 0)
		std::cout << (int)((double)total_nodes/total_time);
	else 
		std::cout << "N/A";
	std::cout << std::endl
		<< " d " << prefix << std::left << std::setw(25-plength)  << " BACKTRACKS/s "  << std::right << std::setw(19) ;
	if(total_time > 0)
		std::cout << (int)((double)total_backtracks/total_time);
	else 
		std::cout << "N/A";
	std::cout << std::endl
		<< " d " << prefix << std::left << std::setw(25-plength)  << " FAILS/s "       << std::right << std::setw(19) ;
	if(total_time > 0)
		std::cout << (int)((double)total_fails/total_time);
	else 
		std::cout << "N/A";
	std::cout << std::endl
		<< " d " << prefix << std::left << std::setw(25-plength)  << " PROPAGS/s "     << std::right << std::setw(19) ;
	if(total_time > 0)
		std::cout << (int)((double)total_propags/total_time);
	else 
		std::cout << "N/A";
	std::cout 
		<< std::endl
			//<< " d " << prefix << std::left << std::setw(25-plength)  << "RESTARTS "      << std::right << std::setw(19) << total_restarts << std::endl
    
			<< " d " << prefix << std::left << std::setw(25-plength)  << " SOLUTIONS "     << std::right << std::setw(19) << num_solutions << std::endl
				<< " d " << prefix << std::left << std::setw(25-plength)  << " AVGCUTOFF "     << std::right << std::setw(19) << avg_cutoff_time << std::endl
					<< " d " << prefix << std::left << std::setw(25-plength)  << " AVGDISTANCE "   << std::right << std::setw(19) << avg_distance << std::endl
						<< " d " << prefix << std::left << std::setw(25-plength)  << " MINDISTANCE "   << std::right << std::setw(19) << min_distance << std::endl
							<< " d " << prefix << std::left << std::setw(25-plength)  << " MAXDISTANCE "   << std::right << std::setw(19) << max_distance << std::endl
								<< " d " << prefix << std::left << std::setw(25-plength)  << " OPTIMAL "       << std::right << std::setw(19) << (lower_bound == upper_bound) << std::endl
									<< " c +==============[ statistics ]===============+" << std::endl
										<< std::endl;

	return os;
}  


const char* ParameterList::int_ident[ParameterList::nia] = 
	{"-ub", "-lb", "-check", "-seed", "-cutoff", "-dichotomy", 
"-base", "-randomized", "-verbose", "-optimise", "-nogood", 
"-dyncutoff", "-nodes", "-hlimit", "-init", "-neighbor", 
"-initstep", "-fixtasks", "-order", "-ngdt"};

const char* ParameterList::str_ident[ParameterList::nsa] = 
	{"-heuristic", "-restart", "-factor", "-decay", "-type", 
"-value", "-dvalue", "-ivalue", "-objective", 
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

	get_command_line(ParameterList::int_ident,
	int_param,
	ParameterList::nia,
	ParameterList::str_ident,
	str_param,
	ParameterList::nsa,
	&commandline[1],length-1);

	if(strcmp(str_param[4],"nil")) Type = str_param[4];
	else {
		Type = "jsp";
		std::cout << " c Warning: no type specified, treating the data as Taillard's jsp" << std::endl;
	}

	UBinit      = -1;
	LBinit      = -1;
	Checked     = true;
	Seed        = 12345;
	Cutoff      = 300;
	NodeCutoff  = 0;
	NodeBase    = 20;
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
	Objective = "makespan";
	Algorithm = "bnb";
	Presolve  = "default";

	Heuristic = "none";
	PolicyRestart = GEOMETRIC;
	FixTasks  = 0;
	NgdType   = 2;
	OrderTasks = 1;



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

	if(strcmp(str_param[0 ],"nil")) Heuristic  = str_param[0];
	if(strcmp(str_param[1 ],"nil")) Policy     = str_param[1];
	if(strcmp(str_param[2 ],"nil")) Factor     = atof(str_param[2]);
	if(strcmp(str_param[3 ],"nil")) Decay      = atof(str_param[3]);
	if(strcmp(str_param[5 ],"nil")) Value      = str_param[5];
	if(strcmp(str_param[6 ],"nil")) DValue     = str_param[6];
	if(strcmp(str_param[7 ],"nil")) IValue     = str_param[7];
	if(strcmp(str_param[8 ],"nil")) Objective  = str_param[8];
	if(strcmp(str_param[9 ],"nil")) Algorithm  = str_param[9 ];
	if(strcmp(str_param[10],"nil")) Presolve   = str_param[10];

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
		PolicyRestart = NORESTART;
}


std::ostream& ParameterList::print(std::ostream& os) {
	os << " c +==============[ parameters ]===============+" << std::endl;
	os << std::left << std::setw(30) << " c | data file " << ":" << std::right << std::setw(15) << data_file_name << " |" << std::endl;
	os << std::left << std::setw(30) << " c | type " << ":" << std::right << std::setw(15) << Type << " |" << std::endl;
	os << std::left << std::setw(30) << " c | seed " << ":" << std::right << std::setw(15) << Seed << " |" << std::endl;
	os << std::left << std::setw(30) << " c | greedy iterations " << ":" << std::right << std::setw(15) << InitBound << " |" << std::endl;
	os << std::left << std::setw(30) << " c | use initial probe " << ":" << std::right << std::setw(15) << (InitStep ? "yes" : "no") << " |" << std::endl;
	os << std::left << std::setw(30) << " c | time cutoff " << ":" << std::right << std::setw(15) << Cutoff << " |" << std::endl;
	os << std::left << std::setw(30) << " c | node cutoff " << ":" << std::right << std::setw(15) << NodeCutoff << " |" << std::endl;
	os << std::left << std::setw(30) << " c | dichotomy " << ":" << std::right << std::setw(15) << (Dichotomy ? "yes" : "no") << " |" << std::endl;
	os << std::left << std::setw(30) << " c | restart policy " << ":" << std::right << std::setw(15) << Policy << " |" << std::endl;
	os << std::left << std::setw(30) << " c | base " << ":" << std::right << std::setw(15) << Base << " |" << std::endl;
	os << std::left << std::setw(30) << " c | factor " << ":" << std::right << std::setw(15) << Factor << " |" << std::endl;
	os << std::left << std::setw(30) << " c | heuristic " << ":" << std::right << std::setw(15) << Heuristic << " (" << abs(Randomized) << ")" << " |" << std::endl;
	os << std::left << std::setw(30) << " c | value ord. (init step) " << ":" << std::right << std::setw(15) << IValue << " |" << std::endl;
	os << std::left << std::setw(30) << " c | value ord. (dichotomy) " << ":" << std::right << std::setw(15) << DValue << " |" << std::endl;
	os << std::left << std::setw(30) << " c | value ord. (optim) " << ":" << std::right << std::setw(15) << Value << " |" << std::endl;
	os << std::left << std::setw(30) << " c | randomization " << ":" ;
     
	if(Randomized<-1) 
		os << std::right
			<< std::setw(12) << "i-shuff random (" << std::setw(2) << -Randomized << ")" << " |" << std::endl;
	else if(Randomized>1) 
		os << std::right
			<< std::setw(12) << "shuff & random (" << std::setw(2) << Randomized << ")" << " |" << std::endl;
	else if(Randomized==1) 
		os << std::right
			<< std::setw(15) << "shuffled" << " |" << std::endl;
	else if(Randomized==0) 
		os << std::right
			<< std::setw(15) << "not random" << " |" << std::endl;
	else 
		os << std::right
			<< std::setw(15) << "init shuffle" << " |" << std::endl;
  
	os << " c +==============[ parameters ]===============+" << std::endl;
  
	return os;
}


Instance::Instance(ParameterList& params) {

	DBG("Build instance %s\n", params.data_file);

	dtp_nodes = 0;
  
	setup_time    = NULL;
	time_lag[0]   = NULL;
	time_lag[1]   = NULL;
	jsp_duedate   = NULL;
	// jsp_latecost  = NULL;
	// jsp_earlycost = NULL;
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

	// if(hasLateCost()) {
	// 	delete [] jsp_latecost;
	// }
	//
	// if(hasEarlyCost()) {
	// 	delete [] jsp_earlycost;
	// }
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
		forbidden_starts.add(timetable_x[k]-timetable_y[k]-duration_y[k]+1);
		forbidden_ends.add(timetable_x[k]-timetable_y[k]+duration_x[k]);
	}


	// Now the cool part, we want to compute the 'real' intervals, since
	// some that we just computed can be merged if they overlap

	// we sort the intervals ends
	qsort(forbidden_starts.stack_, forbidden_starts.size, sizeof(int), straight_compar);
	qsort(forbidden_ends.stack_, forbidden_ends.size, sizeof(int), straight_compar);

	unsigned int i=0, j=0;
	int current = 0;    
  
	// now we can go forward from the earliest interval to latest
	while(i<forbidden_starts.size && j<forbidden_ends.size) {
		if(forbidden_starts[i]<forbidden_ends[j]) {
			++i;

			// when we see the start of an interval, the current number
			// of machine that forbids this start time increases
			if(current++==0) {
				// if it goes from 0 to 1, it is the start of a forbidden interval
				intervals.add(forbidden_starts[i-1]-1);
			}
		} else if(forbidden_starts[i]>forbidden_ends[j]) {
			++j;

			// when we see the end of an interval, the current number
			// of machine that forbids this start time decreases
			if(--current==0) {
				// if it goes from 1 to 0, it is the end of a forbidden interval
				intervals.add(forbidden_ends[j-1]);
			}
		} else {
			++i;
			++j;
		}
	}

	intervals.add(forbidden_ends.back());

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
	os << " c " << (nJobs()) << " jobs, " 
		<< nMachines() << " machines ("
			<< nTasks() << " tasks)" << std::endl;
	for(int i=0; i<nJobs(); ++i) {
		if(nTasksInJob(i) > 1) {
			os << " c ";
			for(int j=1; j<nTasksInJob(i); ++j)
				os << "  t" << tasks_in_job[i][j-1] << "+" << (duration[tasks_in_job[i][j-1]]) 
					<< " <= t" << tasks_in_job[i][j];
			os << std::endl;
		}
	}
	for(int i=0; i<nMachines(); ++i) {
		if(tasks_in_machine[i].size() > 0) {
			os << " c machine" << i << ": t" << tasks_in_machine[i][0];
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
	os << " c +===============[ instance ]================+" << std::endl
		<< " d " << std::left << std::setw(25)  << " NUMTASKS "      << std::right << std::setw(19) << nTasks() << std::endl
			<< " d " << std::left << std::setw(25)  << " NUMJOBS "       << std::right << std::setw(19) << nJobs() << std::endl
				<< " d " << std::left << std::setw(25)  << " NUMMACHINES "   << std::right << std::setw(19) << nMachines() << std::endl
					<< " d " << std::left << std::setw(25)  << " NUMDISJUNCTS "  << std::right << std::setw(19) << nDisjuncts() << std::endl
						<< " d " << std::left << std::setw(25)  << " NUMPRECEDENCES "<< std::right << std::setw(19) << nPrecedences() << std::endl
							//     << " d " << std::left << std::setw(25)  << "LBMAKESPAN "    << std::right << std::setw(19) << lb_C_max << std::endl
							//     << " d " << std::left << std::setw(25)  << "UBMAKESPAN "    << std::right << std::setw(19) << ub_C_max << std::endl
							<< " c +===============[ instance ]================+" << std::endl;
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

	int i, j, dur, mach, nJobs, nMachines, c;
	std::string tag;
	std::ifstream infile( filename, std::ios_base::in );
	
	infile >> nJobs;
	infile >> nMachines;
	
	jsp_duedate = new int[nJobs];
	//jsp_earlycost = new int[nJobs];
	//jsp_latecost = new int[nJobs];
	
	for(i=0; i<nJobs; ++i) {
		for(j=0; j<nMachines; ++j) {
			
			infile >> mach;
			infile >> dur;
			addTask(dur, i, mach);

		}
		
		infile >> jsp_duedate[i];
		infile >> c;
		jsp_earlycost.add(c);
		infile >> c;
		jsp_latecost.add(c);
	}

}

void Instance::dyn_readData( const char* filename, const int precision ) {

	DBG("Read (dyn)%s\n", "");

	int i, j, k, dur, mach, nJobs, nMachines;
	std::string tag;
	std::ifstream infile( filename, std::ios_base::in );
	
	infile >> nJobs;
	infile >> nMachines;
	
	DBG("%i jobs, %i machines\n", nJobs, nMachines);
	
	jsp_duedate = new int[nJobs];
	// jsp_earlycost = new int[nJobs];
	// jsp_latecost = new int[nJobs];
	jsp_floatcost = new double[nJobs];

	for(i=0; i<nJobs; ++i) {
		k = release_date.size();

		for(j=0; j<nMachines; ++j) {
      
			infile >> mach;
			infile >> dur;
      
			if(mach != -1 && dur != -1) addTask(dur, i, mach);
			
			DBG("  -> new task: job=%i, duration=%i, machine=%i\n", i, dur, mach);
			
		}	
		

		infile >> jsp_duedate[i];
		infile >> j;
		//infile >> jsp_latecost[i];
		infile >> jsp_floatcost[i];
		
		
		DBG("job_%i, must end at %i, penalty=%f\n", i, jsp_duedate[i], jsp_floatcost[i]);
    
		setReleaseDate(k, j);
		
		DBG("task_%i, can start at %i\n", k, j);
		
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
		jsp_earlycost.add(approx);
		jsp_latecost.add(approx);
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



SchedulingSolver::SchedulingSolver(Instance *pb, 
ParameterList *pr, 
StatisticList* st) {

	data = pb;
	params = pr;
	stats = st;
	stats->solver = this;
  
	params->initialise(this);

	//setup(); //inst, params, max_makespan);  

	//  stats->lower_bound = get_lb();//lower_bound;
	//  stats->upper_bound = get_ub();//upper_bound;
  
	//nogoods = NULL;
	//pool = new SolutionPool();


}

std::ostream& SchedulingSolver::printStats(std::ostream& os) {
	os << " c +=================[ model ]=================+" << std::endl
		<< " d " << std::left << std::setw(25)  << "LBMAKESPAN "    << std::right << std::setw(19) << lb_C_max << std::endl
			<< " d " << std::left << std::setw(25)  << "UBMAKESPAN "    << std::right << std::setw(19) << ub_C_max << std::endl
				<< " c +=================[ model ]=================+" << std::endl;
	return os;
}


void SchedulingSolver::setup() {

	int i,j,k, lb, ub, ti, tj, rki, rkj, hi, hj, aux;
  


	lb_C_max = (params->LBinit<0 ? data->getMakespanLowerBound() : params->LBinit);
	ub_C_max = (params->UBinit<0 ? data->getMakespanUpperBound(params->InitBound) : params->UBinit);
	if(params->Objective == "tardiness") {
		int max_due_date = data->getJobDueDate(0);
		for(int i=1; i<data->nJobs(); ++i) {
			if(max_due_date < data->getJobDueDate(i))
				max_due_date = data->getJobDueDate(i);
		}
		if(max_due_date < ub_C_max) max_due_date = ub_C_max;
		ub_C_max += max_due_date;
	}

	//for(i=0; i<data->nJobs(); ++i)
	//ub_C_max += data->getDuration(data->getLastTaskofJob(i));

	lb_L_sum = data->getEarlinessTardinessLowerBound(ub_C_max);
	ub_L_sum = data->getEarlinessTardinessUpperBound(ub_C_max);


	// create one variable per task
	for(i=0; i<data->nTasks(); ++i) {
		lb = data->getReleaseDate(i);
		ub = std::min(ub_C_max, data->getDueDate(i)) - data->getDuration(i);

		if(lb > ub) {
			std::cout << "INCONSISTENT" << std::endl;
			exit(1);
		}

		Variable t(lb, ub);
		tasks.add(t);
	}




	// precedence constraints
	for(i=0; i<data->nJobs(); ++i) 
	for(j=1; j<data->nTasksInJob(i); ++j) {
		ti = data->getJobTask(i,j-1);
		tj = data->getJobTask(i,j);
		add( Precedence(tasks[ti], 
		(data->getDuration(ti) + 
			(data->hasTimeLag() ? data->getMinLag(i,j-1) : 0)), 
		tasks[tj]) );
	}

	// time lags constraints
	if(data->hasTimeLag()) {
		for(i=0; i<data->nJobs(); ++i) 
		for(j=1; j<data->nTasksInJob(i); ++j) if(data->getMaxLag(i,j-1) >= 0) {
			ti = data->getJobTask(i,j-1);
			tj = data->getJobTask(i,j);
			add( Precedence(tasks[tj], 
			-(data->getDuration(ti)+data->getMaxLag(i,j-1)), 
			tasks[ti]) );
		}
	}


	// mutual exclusion constraints
	for(k=0; k<data->nMachines(); ++k) {
		for(i=0; i<data->nTasksInMachine(k); ++i) {
			for(j=i+1; j<data->nTasksInMachine(k); ++j) {
				ti = data->getMachineTask(k,i);
				tj = data->getMachineTask(k,j);
	
				if(params->OrderTasks==1) {
					rki = data->getRankInJob(ti);
					rkj = data->getRankInJob(tj);
					hi = data->getHeadInJob(ti);
					hj = data->getHeadInJob(tj);
	  
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
				//  		  << " (" << data->getDuration(ti) << "." 
				// 		  << data->getDuration(tj) << ")" 
				// 		  << std::endl; 

				first_job.add(data->getJob(ti,0));
				second_job.add(data->getJob(tj,0));

				//first_task_of_disjunct.push_back(ti);
				//second_task_of_disjunct.push_back(tj);

				disjuncts.add( ReifiedDisjunctive( tasks[ti],
				tasks[tj],
  
				data->getDuration(ti)
					+(data->hasSetupTime() ? data->getSetupTime(k,ti,tj) : 0),
				data->getDuration(tj)
					+(data->hasSetupTime() ? data->getSetupTime(k,tj,ti) : 0) ) );


				Vector< Variable > tasks_of_d;
				tasks_of_d.add(disjuncts.back().get_var());
				tasks_of_d.add(tasks[ti].get_var());
				tasks_of_d.add(tasks[tj].get_var());
				disjunct_map.add(tasks_of_d);


			}
		}
	}

  
	for(unsigned int i=0; i<disjuncts.size; ++i)
		add(Free(disjuncts[i]));
	//add( disjuncts );

	//   //exit(1);


	setup_objective();

	pool = new SolutionPool();
}

SchedulingSolver::~SchedulingSolver() {}





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

void C_max_Model::setup_objective() {
	Variable x_cmax(lb_C_max, ub_C_max);
	C_max = x_cmax;
	
	int ti;
	for(int i=0; i<data->nJobs(); ++i) {
		ti = data->getLastTaskofJob(i);
		add(Precedence(tasks[ti], data->getDuration(ti), C_max));
	}
}

Variable L_sum_Model::get_objective_var() {
	return L_sum;
}

Variable Depth_Model::get_objective_var() {
	return Depth;
}

Variable DTP_Model::get_objective_var() {
	return Weight;
}

Variable C_max_Model::get_objective_var() {
	return C_max;
}

int L_sum_Model::set_objective(const int obj) {
	//VariableInt *x_lsum = L_sum.getVariable();
	return (L_sum.set_max(obj) != FAIL_EVENT ? UNKNOWN : UNSAT);
}

void L_sum_Model::setup_objective() {
		
	int ti;
	for(int j=0; j<data->nJobs(); ++j) {
		ti = data->getLastTaskofJob(j);
					
		
		earlybool.add(tasks[ti] < data->getJobDueDate(j)-data->getDuration(ti));
		 latebool.add(tasks[ti] > data->getJobDueDate(j)-data->getDuration(ti));
	 }
	
	
	Variable EarlyCost = BoolSum(earlybool, data->getEarlyCosts());
	Variable  LateCost = BoolSum( latebool, data->getLateCosts ());
	
	
	Variable x_lsum(lb_L_sum, ub_L_sum);
	L_sum = x_lsum;

	add(L_sum == (EarlyCost + LateCost));
	
}

int Depth_Model::set_objective(const int obj) {
	//   //VariableInt *x_lsum = Depth.getVariable();
	//   Depth.print(std::cout);
	//   std::cout << std::endl;
	//   Depth.getVariable()->print(std::cout);
	//   std::cout << std::endl;
	return (Depth.set_max(obj) != FAIL_EVENT ? UNKNOWN : UNSAT);
}

void Depth_Model::setup_objective() {
	std::cout << " Depth: objective not supported!" << std::endl;
	exit(1);
}	

int DTP_Model::set_objective(const int obj) {
	return (Weight.set_max(obj) != FAIL_EVENT ? UNKNOWN : UNSAT);
}

int C_max_Model::set_objective(const int obj) {
	return (C_max.set_max(obj) != FAIL_EVENT ? UNKNOWN : UNSAT);
}

int L_sum_Model::get_objective() {
	return L_sum.get_solution_min();
}

int Depth_Model::get_objective() {
	return Depth.get_solution_min();
}

int DTP_Model::get_objective() {
	return Weight.get_solution_min();
}

void DTP_Model::setup_objective() {
	std::cout << " DTP: objective not supported!" << std::endl;
	exit(1);
}	

int C_max_Model::get_objective() {
	return C_max.get_solution_min();
}


SchedulingSolution::SchedulingSolution(SchedulingSolver *s) {
	//model = m;
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

	for(i=0; i<n; ++i) {
		all_value[i] = solver->variables[i].get_solution_int_value();
		//     std::cout << "store " ;
		//     solver->variables[i]->print(std::cout);
		//     std::cout << std::endl;
	}


	n = s->disjuncts.size;
	if(n) {
		disjunct_value = new int[n];
		for(i=0; i<n; ++i) disjunct_value[i] = s->disjuncts[i].get_solution_int_value();
	} 

	n = s->SearchVars.size;
	if(n) {
		search_value = new int[n];
		for(i=0; i<n; ++i) search_value[i] = s->SearchVars[i].get_solution_int_value();
	}

	if(s->data->hasJobDueDate()) {
		n = s->earlybool.size;
		if(n) {
			earlybool_value = new int[n];
			for(i=0; i<n; ++i) earlybool_value[i] = s->earlybool[i].get_solution_int_value();
		} 
		n = s->latebool.size;
		if(n) {
			latebool_value = new int[n];
			for(i=0; i<n; ++i) latebool_value[i] = s->latebool[i].get_solution_int_value();
		} 
		n = s->tasks.size;
		if(n) {
			task_min = new int[n];
			task_max = new int[n];
			for(i=0; i<n; ++i) {
				task_min[i] = s->tasks[i].get_solution_min();
				task_max[i] = s->tasks[i].get_solution_max();
			} 
		}

		n = s->last_tasks.size;
		if(n) {
			ltask_value = new int[n];
			for(i=0; i<n; ++i) {
				ltask_value[i] = s->last_tasks[i].get_solution_int_value();
			}
		}
	}
}


SchedulingSolution::~SchedulingSolution() {
	delete [] earlybool_value;
	delete [] latebool_value;
	delete [] task_min;
	delete [] task_max;
	delete [] ltask_value;
	delete [] disjunct_value;
	delete [] search_value;
	delete [] all_value;
}

void SchedulingSolution::guide_search() {
  
	// //   std::cout << "guide search with ";
	// //   print(std::cout);
	// //   std::cout << std::endl;

	// //   //solver->setGuidedOrdering(model->disjuncts, disjunct_value, "spl");
	// //   if( ((SchedulingSolver*)solver)->params->Type == "now2" )
	// //     solver->setGuidedBoundsOrdering(model->SearchVars, search_value);
	// //   else
    
	//     solver->setGuidedOrdering(model->SearchVars, search_value);

	//   if(model->data->hasJobDueDate()) {
	//     solver->setGuidedOrdering(model->earlybool, earlybool_value);
	//     solver->setGuidedOrdering(model->latebool, latebool_value);
	//     solver->setGuidedOrdering(model->last_tasks, ltask_value, "nbd");
	//   }
}

void SchedulingSolution::guide_search_bounds() {
  
	// //   std::cout << "guide bound search with ";
	// //   print(std::cout);
	// //   std::cout << std::endl;

	// //   //solver->setGuidedOrdering(model->disjuncts, disjunct_value, "spl");
	// //   if( ((SchedulingSolver*)solver)->params->Type == "now2" )
	// //     solver->setGuidedBoundsOrdering(model->SearchVars, search_value);
	// //   else
    
	//     solver->setGuidedOrdering(model->SearchVars, search_value);
	//   if(model->data->hasJobDueDate()) {
	//     solver->setGuidedOrdering(model->earlybool, earlybool_value);
	//     solver->setGuidedOrdering(model->latebool, latebool_value);
	//     solver->setGuidedOrdering(model->last_tasks, ltask_value, "nbd");
	//   }
}

double SchedulingSolution::distance(SchedulingSolution* s) {
  
	int i, //n = model->disjuncts.size;
	n = solver->SearchVars.size;
	double dist = 0.0;

	for(i=0; i<n; ++i) {
		//dist += (double)(disjunct_value[i] != s->disjunct_value[i]);
		dist += (double)(search_value[i] != s->search_value[i]);
	} 
	if(solver->data->hasJobDueDate()) {
		n = solver->last_tasks.size;
		for(i=0; i<n; ++i) {dist += 
			((ltask_value[i] - s->ltask_value[i])*(ltask_value[i] - s->ltask_value[i]));
		}
	}
  
	return sqrt(dist);
}

std::ostream& SchedulingSolution::print(std::ostream& os, std::string type) {
	int i, //n = solver->disjuncts.size;
	n = solver->SearchVars.size;

	//   for(i=0; i<n; ++i) {
	//     os << disjunct_value[i] ;
	//   } 
	//   os << std::endl;

	if(type == "fsp") {
		int j, m = solver->data->nJobs(), k = 0;
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

		os << " c ";

		for(i=0; i<m; ++i) {
			os << (order[i]+1) << " ";
		}

		//os << std::endl;

		delete [] rank;
		delete [] order;

	} else {
		for(i=0; i<n; ++i) {
			os << std::setw(3) << solver->SearchVars[i].id() ;
		}
		os << std::endl;
		for(i=0; i<n; ++i) {
			os << std::setw(3) << search_value[i] ;
		} 
		os << std::endl;
	}
	return os;
}



void SchedulingSolver::dichotomic_search()
{
  
	std::cout << this <<std::endl;

	stats->lower_bound = get_lb();
	stats->upper_bound = get_ub();
  
	int iteration = 1;
  
	int minfsble = stats->lower_bound;
	int maxfsble = stats->upper_bound;
  
	int objective = -1;
	int new_objective = -1;
	//   int ngd_stamp = 0;
	//   int lit_stamp = 0;
  
  
	parameters.verbosity = params->Verbose;
	parameters.time_limit = params->Cutoff;
  


	BranchingHeuristic *heu = new SchedulingWeightedDegree < TaskDomOverBoolWeight, Guided< MinValue >, 2 > (this, disjunct_map);

	//BranchingHeuristic *heu = new GenericHeuristic < NoOrder, MinValue > (this);

	RestartPolicy *pol = new Geometric();


	initialise_search(disjuncts, heu, pol);


	//propagate the bounds, with respect to the initial upper bound
	Outcome result = (propagate() ? UNKNOWN : UNSAT);

  
	////////// dichotomic search ///////////////
	while( //result == UNKNOWN && 
		minfsble<maxfsble && 
			iteration<params->Dichotomy
	) {

    
		double remaining_time = params->Optimise - stats->get_total_time();
    
		if(remaining_time < (2*params->NodeBase)) break;
    
		objective = (int)(floor(((double)minfsble + (double)maxfsble)/2));
		std::cout << "\n c +=========[ start dichotomic step ]=========+" << std::endl;
		//       setPropagsLimit(params->NodeCutoff);
    
		parameters.propagation_limit = params->NodeCutoff;
    
    
		std::cout << std::left << std::setw(30) << " c | current real range" << ":" 
			<< std::right << " " << std::setw(5) << stats->lower_bound 
				<< " to " << std::setw(5) << stats->upper_bound << " |" << std::endl;
		std::cout << std::left << std::setw(30) << " c | current dichotomic range" << ":" 
			<< std::right << " " << std::setw(5) << minfsble 
				<< " to " << std::setw(5) << maxfsble << " |" << std::endl;
		std::cout << std::left << std::setw(30) << " c | target objective" << ":"  
			<< std::right << std::setw(15) << objective << " |" << std::endl;
   
    
		statistics.start_time = get_run_time();

		save();
    
		result = set_objective(objective);
    

#ifdef _DEBUG_PRUNING
		monitor(tasks);
		monitor(C_max);
#endif
    

		result = restart_search(level);

    
		if( result == SAT ) {
			new_objective = get_objective();
      
			// 	stats->normalized_objective = solver->get_normalized_objective();
      
			// at level 0, deduce new bounds for all variables with respect to the new objective
			set_objective(new_objective);
			propagate();

			maxfsble = new_objective;
			pool->add(new SchedulingSolution(this));
      
			// 	if(nogoods) {
			// 	  for(int i=ngd_stamp; i<nogoods->base->nogood.size; ++i)
			// 	    stats->avg_nogood_size += (double)(nogoods->base->nogood[i]->size);
			// 	  if(params->Rngd>1) stats->num_nogoods = nogoods->base->nogood.size;
			// 	  else stats->num_nogoods += nogoods->base->nogood.size;
			// 	}
      
			//std::cout << std::left << std::setw(30) << " c new upper bound" << ":" << std::right << std::setw(20) << new_objective << " |" << std::endl;
			std::cout << std::left << std::setw(30) << " c | new upper bound" << ":" << std::right << std::setw(15) << new_objective << " |" << std::endl;
      
			//pool->getBestSolution()->print(std::cout);
      
		} else {

			new_objective = objective;
			minfsble = objective+1;

			if( result == UNSAT ) {
				std::cout << std::left << std::setw(30) << " c | real lower bound" << ":" << std::right << std::setw(15) << minfsble << " |" << std::endl;
			} else {
				std::cout << std::left << std::setw(30) << " c | dichotomic lower bound" << ":" << std::right << std::setw(15) << minfsble << " |" << std::endl;
			}
      

      
			// 	if(nogoods) {
			// 	  nogoods->forget(ngd_stamp);
			// 	  nogoods->reinit();
			// 	}
		}
      
		stats->add_info(new_objective, DICHO);

		std::cout << std::left << std::setw(30) << " c | cpu time" << ":" << std::right << std::setw(15) << (double)((int)((stats->time.back())*10000))/10000.0 << " |" << std::endl;
    
		//printStatistics(std::cout, ((params->Verbose ? RUNTIME : 0) + ((params->Verbose || result != UNKNOWN)  ? BTS + PPGS : 0) + OUTCOME) );
    
    
		//std::cout << "LEVEL: " << level << " " << this << std::endl;

		restore();
		statistics.initialise(this);
		pol->initialise(parameters.restart_limit);

		// std::cout << std::left << std::setw(30) << " c current dichotomic range" << ":" 
		// 	      << std::right << std::setw(6) << " " << std::setw(5) << minfsble 
		// 	      << " to " << std::setw(5) << maxfsble << " " << iteration << " " << params->Dichotomy << std::endl;
		std::cout << " c +==========[ end dichotomic step ]==========+" << std::endl;
    

		++iteration;
	} 
	//   } else if( status == SAT ) {
	//     std::cout << " c Solved during preprocessing!" << std::endl;
  
	//   } else if( status == UNSAT ) {
	//     std::cout << " c Found inconsistent during preprocessing!" << std::endl;
  
	//   }
    
	std::cout << std::endl;
}
 


void SchedulingSolver::branch_and_bound()
{
	//int ngd_stamp = 0;
	//int lit_stamp = 0;
	//resetNodeLimit();
	//resetPropagsLimit();
	parameters.propagation_limit = 0;

	statistics.start_time = get_run_time();

	//std::cout << (get_run_time() - statistics.start_time) << std::endl;

	save();
	set_objective(stats->upper_bound-1);
	addObjective();


	//std::cout << (get_run_time() - statistics.start_time) << std::endl;

	parameters.verbosity = 2;
	//setVerbosity(params->Verbose);
	//setRandomSeed( params->Seed );

	double time_limit = (params->Optimise - stats->get_total_time());

	if(time_limit > 0) {
		set_time_limit( time_limit ); 
		// //addHeuristic( params->Heuristic, params->Randomized, params->Value, params->Hlimit );
		// if(params->Value == "guided") {
		//   //function = new SolutionGuidedSearch( this, pool, stats );
		//   if(pool->size()) pool->getBestSolution()->guide_search();
		// } else {
		//   //function = new StoreStats( this, stats );
		// }
    
		std::cout << " c +=========[ start branch & bound ]==========+" << std::endl;
		std::cout << std::left << std::setw(26) << " c | current range" << ":" 
			<< std::right << std::setw(5) << " " << std::setw(5) << stats->lower_bound 
				<< " to " << std::setw(5) << objective->upper_bound << " |" << std::endl;
		std::cout << std::left << std::setw(26) << " c | run for " << ":"
			<< std::right << std::setw(18) << (time_limit) << "s |" << std::endl;
    

		//std::cout << (get_run_time() - statistics.start_time) << std::endl;
		//std::cout << C_max << " in " << C_max.get_domain() << std::endl;

		//std::cout << this << std::endl;

		//std::cout << level << std::endl;

    

    
		Outcome result = (stats->upper_bound >= stats->lower_bound ? UNKNOWN : UNSAT);

		if(result == UNKNOWN) {

			// if(nogoods) {
			// 	ngd_stamp = nogoods->base->nogood.size;
			// 	lit_stamp = sUnaryCons.size;
			// }

			//solve_and_restart(params->PolicyRestart, params->Base, params->Factor);

			//std::cout << (get_run_time() - statistics.start_time) << std::endl;

			result = restart_search(level);

			// if(nogoods) {
			// 	for(int i=ngd_stamp; i<nogoods->base->nogood.size; ++i)
			// 	  stats->avg_nogood_size += (double)(nogoods->base->nogood[i]->size);
			// 	if(params->Rngd>1) stats->num_nogoods = nogoods->base->nogood.size;
			// 	else stats->num_nogoods += nogoods->base->nogood.size;
			// }
		}

		stats->add_info(objective->upper_bound, BNB);
    
		//printStatistics(std::cout, ((params->Verbose ? RUNTIME : 0) + (params->Verbose ? BTS + PPGS : 0) + OUTCOME) );
    
		//reset(true);
		std::cout << " c +==========[ end branch & bound ]===========+" << std::endl;
    
	}
}




void SchedulingSolver::print_solution(std::ostream& os, std::string type)
{
	pool->getBestSolution()->print(os, type);
	os << std::endl;
}
