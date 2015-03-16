 
#include <mistral_solver.hpp>
#include <mistral_search.hpp>
#include <mistral_variable.hpp>
#include <iostream>
#include <fstream>
#include <string>

//#define _DEBUGSCHED true

#ifdef _DEBUGSCHED
#define DBG(fmt, args...) printf("dbg - l%d: "fmt,__LINE__,args)
#else
#define DBG(fmt, args...)
#endif


namespace Mistral {

	class SchedulingSolver;
	class StatisticList {
	private:
		int best_solution_index;
		int branch_and_bound_index;

	public:

		SchedulingSolver *solver;

		int    num_nogoods;
		int    num_solutions;
		int    lower_bound;
		int    upper_bound;

		double avg_nogood_size;
		double avg_cutoff_time;
		double avg_distance;
		double min_distance;
		double max_distance;

		double normalized_objective;
		//double normalizer;
    
		double real_start_time;

		std::vector<long unsigned int> nodes;
		std::vector<long unsigned int> backtracks;
		std::vector<long unsigned int> fails;
		std::vector<long unsigned int> propags;
		std::vector<double>            time;
		std::vector<double>            soltime;
		std::vector<int>               outcome;
		std::vector<int>               types;

		StatisticList();
		virtual ~StatisticList();

		bool solved();
		double get_total_time();
		double get_lowerbound_time();
		void start();
		void add_info(const int obj, const int tp);
		//void add_info(const int lb, const int ub);

		std::ostream& print(std::ostream& os, 
		const char* prefix=" ",
		const int start=0, 
		const int end=-1);
	};

	class Instance;
	class ParameterList {
	public:

		static const int LEX     =  0;
		static const int ANTILEX = -2;
		static const int PROMISE =  1;
		static const int ANTI    = -1;
		static const int GUIDED  =  2;
		static const int RGUIDED =  3;
		static const int RAND    =  4;

		static const int nia = 20;
		static const char* int_ident[nia];
    
		static const int nsa = 11;
		static const char* str_ident[nsa];
    
    
		const char* str_param[nsa];
		int int_param[nia];

		char* data_file;
		char* data_file_name;

		SchedulingSolver *solver;

		int UBinit; // "ub": user defined upper bound (-1 if not)
		int LBinit; // "lb": user defined lower bound (-1 if not)
		int Checked; // "check": whether the solution is checked
		int Seed; // "seed": random seed
		int Cutoff; // "cutoff": time cutoff of dichotomic steps
		unsigned long long int NodeCutoff; // "nodes": node cutoff of dichotomic steps
		int NodeBase; // "dyncutoff": node cutoff trying to mimic time cutoff
		int Dichotomy; // "dichotomy": max number of dichotomic steps
		int Base; // "base": base cutoff for restarts
		int Randomized; // "randomized": level of randomization
		int Verbose; // "verbose": level of verbosity
		int Optimise; // "optimise": total time cutoff 
		int Rngd; // "nogood": whether nogood on restart are used
		int Precision; // "-": precision when turning float weights into int weights
		int Hlimit; // "hlimit": maximum number of variables evaluated by the dvo
		int InitBound; // "init": number of iterations when initialising the upper bound
		int Neighbor; // "neighbor": number of machines to re-schedule during LNS
		int InitStep; // "initstep": whether an initial probe with high makespan should be performed
		int FixTasks; // Whether tasks' starting time should be fixed (searched)
		//int MinRank; // Whether the sum of the disjunct should be minimised
		int OrderTasks; // Whetheer tasks should be ordered within the disjuncts
		int NgdType; // nogood type for solution removal

		double Factor;
		double Decay;

		std::string Policy;
		std::string Heuristic;
		std::string Type; 
		std::string Value;
		std::string DValue;
		std::string IValue;
		std::string Objective;
		std::string Algorithm;
		std::string Presolve;

		int PolicyRestart;

		ParameterList(int length, char** commandline);
		//ParameterList();
		//ParameterList(const ParameterList& pl);
		virtual ~ParameterList() {}

		void initialise(SchedulingSolver*);
		void initialise(const ParameterList& pl);
		std::ostream& print(std::ostream& os);

	};

	class pair_ {
	public:
    
		int element;
		int rank;

		pair_(const int elt, const int rk) : element(elt), rank(rk) {}
	};

	struct Term {
		int X;
		int Y;
		int k;
	};
	class Instance {

	public:

		//dtp clauses (atoms are terms like x-y <= k)
		int dtp_nodes;
		std::vector< std::vector < Term > > dtp_clauses;
		std::vector< int > dtp_weights;

	private:

		std::vector< std::vector<int> > tasks_in_job;
		std::vector< std::vector<int> > jobs_of_task;
		std::vector< std::vector<int> > tasks_in_machine;
		std::vector< std::vector<int> > machines_of_task;
		std::vector< int >              duration;
		std::vector< int >              release_date;
		std::vector< int >              due_date;

		std::vector< std::vector< pair_ > > task_rank_in_machine;
		std::vector< std::vector< pair_ > > task_rank_in_job;

		int***  setup_time;
		int**   time_lag[2];
		int*    jsp_duedate;
		Vector<int>    jsp_latecost;
		Vector<int>    jsp_earlycost;
		double* jsp_floatcost;

		int max_makespan;


		int addTask(const int dur, const int job, const int machine) ;
		void addTaskToJob(const unsigned int index, const unsigned int j) ;
		void addTaskToMachine(const unsigned int index, const unsigned int j) ;


		void osp_readData( const char* filename );
		//void osp_readData( const char* filename );
		void sds_readData( const char* filename );
		void jtl_readData( const char* filename );
		void now_readData( const char* filename );
		void jla_readData( const char* filename );
		void tsp_readData( const char* filename );
		void fsp_readData( const char* filename );
		void jsp_readData( const char* filename );
		void jet_readData( const char* filename );
		void dtp_readData( const char* filename );
		void dyn_readData( const char* filename, const int p );


		int get_machine_rank(const int ti, const int mj) const {
			int rk=-1, mc;
			for(unsigned int i=0; i<task_rank_in_machine[ti].size(); ++i) {
				mc = task_rank_in_machine[ti][i].element;
				if(mc == mj) {
					rk = task_rank_in_machine[ti][i].rank;
					break;
				}
			}
			return rk;
		}

		int get_job_rank(const int ti, const int jj) const {
			int rk=-1, jb;
			for(unsigned int i=0; i<task_rank_in_job[ti].size(); ++i) {
				jb = task_rank_in_job[ti][i].element;
				if(jb == jj) {
					rk = task_rank_in_job[ti][i].rank;
					break;
				}
			}
			return rk;
		}


	public:

		Instance(ParameterList& params);
		virtual ~Instance();

		std::ostream& print(std::ostream& os);

		int  nJobs() const {return tasks_in_job.size();}
		int  nJobs(const int i) const {return jobs_of_task[i].size();}

		int  nMachines() const {return tasks_in_machine.size();}
		int  nMachines(const int i) const {return machines_of_task[i].size();}

		int  nTasks() const {return duration.size();}
		int  nTasksInJob(const int j) const {return tasks_in_job[j].size();}
		int  nTasksInMachine(const int j) const {return tasks_in_machine[j].size();}

		int  getJobTask(const int i, const int j) const {return tasks_in_job[i][j];}
		int  getMachineTask(const int i, const int j) const {return tasks_in_machine[i][j];}
		int  getLastTaskofJob(const int i) const {return tasks_in_job[i][nTasksInJob(i)-1];}


		int  getJob(const int i, const int j) const {return jobs_of_task[i][j];}
		int  getMachine(const int i, const int j) const {return machines_of_task[i][j];}

		int  getDuration(const int i) const {return duration[i];}
		int  getReleaseDate(const int i) const {return release_date[i];}
		void setReleaseDate(const int i, const int d) {release_date[i] = d;}
		int  getDueDate(const int i) const {return due_date[i];}

		bool hasSetupTime() const {return setup_time != NULL;}
		int  getSetupTime(const int k, const int i, const int j) const;

		bool hasTimeLag() const {return time_lag[0] != NULL;}
		int  getMinLag(const int i, const int j) const {return time_lag[0][i][j];}
		int  getMaxLag(const int i, const int j) const {return time_lag[1][i][j];}

		bool hasJobDueDate() const {return jsp_duedate != NULL;}
		int  getJobDueDate(const int i) const {return jsp_duedate[i];}

		bool hasEarlyCost() const {return !jsp_earlycost.empty(); }
		bool hasLateCost() const {return !jsp_latecost.empty(); }
		int  getJobEarlyCost(const int i) const {return jsp_earlycost[i];}
		int  getJobLateCost(const int i) const {return jsp_latecost[i];}
		const Vector<int>&  getEarlyCosts() const {return jsp_earlycost;}
		const Vector<int>&  getLateCosts() const {return jsp_latecost;}

		bool hasFloatCost() const {return jsp_floatcost != NULL; }
		double getJobFloatCost(const int i) const {return jsp_floatcost[i];}


		int getRankInJob(const int i, const int j=-1) const {return get_job_rank(i,(j==-1?getJob(i,0):j));}
		int getHeadInJob(const int i, const int j=-1) const {
			int rj = (j==-1?getJob(i,0):j);
			int rk = get_job_rank(i,rj);
			int head = 0;
			for(int k=0; k<rk; ++k) {
				head += getDuration(getJobTask(rj,k));
			}
			return head;
		}

		int getMakespanLowerBound();
		int getMakespanUpperBound(const int it=-1);

		int getEarlinessTardinessLowerBound(const int);
		int getEarlinessTardinessUpperBound(const int);

		void get_placement(const int i, const int j, Vector<int>& intervals);
		void get_placement2(const int i, const int j);

		int nDisjuncts() const;
		int nPrecedences() const;

		double getNormalizer() const;

		std::ostream& printStats(std::ostream& os);
	};



	class SchedulingSolution {
	public:
		int *earlybool_value;
		int *latebool_value;
		int *task_min;
		int *task_max;
		int *ltask_value;
		int *disjunct_value;
		int *search_value;
		int *all_value;


	public:
		//SchedulingModel  * model;
		SchedulingSolver *solver;
    
		SchedulingSolution(// SchedulingModel *m, 
		SchedulingSolver *s);
		virtual ~SchedulingSolution();

		//int& operator[](const VariableInt *x) {return all_value[x->id];}

		double distance(SchedulingSolution* s);

		void guide_search();
		void guide_search_bounds();

		std::ostream& print(std::ostream& os, std::string type="jsp");
    
	};

  
	// // used to link each disjunct with the corresponding tasks
	// class DisjunctMap : public Vector< Variable[2] > { 
    
	// public:

	// };


	template < class VarComparator, class Branching, int RAND > 
	class SchedulingWeightedDegree 
		: public GenericHeuristic< GenericDVO< VarComparator, RAND, FailureCountManager >,
	Branching > {
	public:

		Vector< Variable > *neighborhood;
		int offset;
    
		SchedulingWeightedDegree(Solver *s, Vector< Vector< Variable > >& graph) 
			: GenericHeuristic< GenericDVO< VarComparator, RAND, FailureCountManager >, 
		Branching >(s) {

			int min_idx = graph[0][0].id();
			int max_idx = min_idx;
			int idx;
      
			for(unsigned int i=1; i<graph.size; ++i) {
				idx = graph[i][0].id();
				if(min_idx > idx) min_idx = idx;
				if(max_idx < idx) max_idx = idx;
			}

			offset = min_idx;
			neighborhood = new Vector< Variable >[max_idx-min_idx+1];
			neighborhood -= offset;

			for(unsigned int i=0; i<graph.size; ++i) {
				idx = graph[i][0].id();
				for(unsigned int j=1; j<graph[i].size; ++j) {
					neighborhood[idx].add(graph[i][j]);
				}
			}
      
			for(unsigned int i=0; i<RAND; ++i) 
				GenericHeuristic< GenericDVO< VarComparator, RAND, FailureCountManager >, Branching >::var.bests[i].map = neighborhood;
			GenericHeuristic< GenericDVO< VarComparator, RAND, FailureCountManager >, Branching >::var.current.map = neighborhood;
  
			//std::cout << graph << std::endl;
			//exit(1);

		}  

		virtual ~SchedulingWeightedDegree() {
			neighborhood += offset;
			delete [] neighborhood;
		}
	};

	typedef MinNeighborDomainOverNeighborWeight TaskDomOverTaskWeight;
	typedef MinNeighborDomainOverWeight TaskDomOverBoolWeight;


	class SolutionPool;
	class SchedulingSolver : public Solver {
	public:
    
		/// DATA
		Instance *data;


		/// MODEL
		int ub_Weight;
    
		int ub_C_max;
		int lb_C_max;

		int ub_L_sum;
		int lb_L_sum;

		int ub_Depth;
		int lb_Depth;

		Vector < Vector < Variable > > disjunct_map;
    
		//std::vector< int > first_task_of_disjunct;
		//std::vector< int > second_task_of_disjunct;

		VarArray SearchVars;
    
		VarArray tasks;
		VarArray last_tasks;
		VarArray disjuncts;
		VarArray earlybool;
		VarArray latebool;
		Variable C_max;
		Variable L_sum;
		Variable Depth;
		Variable Weight;

		Vector<int> first_job;
		Vector<int> second_job;


		int lower_bound;
		int upper_bound;
		ParameterList  *params;
		StatisticList   *stats;
		SolutionPool     *pool;
		//WeighterRestartGenNogood *nogoods;

		SchedulingSolver() {}
		SchedulingSolver(Instance *pb, ParameterList *pr, StatisticList* st);
		virtual ~SchedulingSolver() ;


		virtual void setup();

		//virtual int virtual_iterative_dfs();

		// std::ostream& print_weights(std::ostream& os);
		// void decay_weights(const double decay);

		void addObjective() { 
			objective = new Goal(Goal::MINIMIZATION, get_objective_var());
		}
    

		// bool probe_ub();
		// void jtl_presolve();
		// void old_jtl_presolve();
		void dichotomic_search();
		// void all_solutions_search();
		void branch_and_bound();
		// //void large_neighborhood_search();
		// //void extract_stable(List& neighbors, Vector<VariableInt*>& stable);
		// //void repair(Solution *sol, Vector<VariableInt*>& stable);


		void print_solution(std::ostream& os, std::string tp);


		virtual int get_lb() = 0;
		virtual int get_ub() = 0;
		virtual Variable get_objective_var() = 0;
		virtual int  get_objective() = 0;
		// virtual double  get_normalized_objective() = 0;
		virtual int  set_objective(const int obj) = 0;
		
		virtual void setup_objective() = 0;
		

		std::ostream& printStats(std::ostream& os);
	};


	class C_max_Model : public SchedulingSolver {
	public:

		C_max_Model() : SchedulingSolver() {}
		C_max_Model(Instance *pb, 
		ParameterList *pr, 
		StatisticList* st) : SchedulingSolver(pb, pr, st) {}
		virtual ~C_max_Model() {}

		virtual int get_lb();
		virtual int get_ub();    
		virtual Variable get_objective_var();
		virtual int  get_objective();
		// virtual double  get_normalized_objective() {return (double)(get_objective())/data->getNormalizer();}
		virtual int  set_objective(const int obj);
		
		virtual void setup_objective();
	};

	class L_sum_Model : public SchedulingSolver {
	public:

		L_sum_Model() : SchedulingSolver() {}
		L_sum_Model(Instance* pb, ParameterList *pr, StatisticList* st// , const int L_sum
			) 
		: SchedulingSolver(pb, pr, st) { 
			//setup(prob, params, L_sum);
		}
		virtual ~L_sum_Model() {}

		virtual int get_lb();
		virtual int get_ub();    
		virtual Variable get_objective_var();
		virtual int  get_objective();
		// virtual double  get_normalized_objective() ;
		virtual int  set_objective(const int obj);
		
		virtual void setup_objective();
	};

	class Depth_Model : public SchedulingSolver {

	public:

		Depth_Model() : SchedulingSolver() {}
		Depth_Model(Instance& prob, ParameterList *params, const int Depth) { // setup(prob, params, Depth);
		}
		virtual ~Depth_Model() {}

		virtual int get_lb();
		virtual int get_ub();    
		virtual Variable get_objective_var();
		virtual int  get_objective();
		// virtual double  get_normalized_objective() {return (double)(get_objective())/(double)(disjuncts.size());}
		virtual int  set_objective(const int obj);
		
		virtual void setup_objective();
	};

	class No_wait_Model : public C_max_Model {
    
	public:
    
		int type;
		VarArray gen_disjuncts;

		Vector<int> job_size;
		Vector<int> job_index;
 
		No_wait_Model() : C_max_Model() {} 
		No_wait_Model(Instance& prob, ParameterList *params, const int C_max, const int t=1) : C_max_Model() { type=t; // setup(prob, params, C_max);
		}
		//virtual void setup(Instance& prob, ParameterList *params, const int C_max);
		virtual ~No_wait_Model() {}
	};

	class DTP_Model : public SchedulingSolver {

	public:

		DTP_Model() : SchedulingSolver() {}
		DTP_Model(Instance& prob, ParameterList *params, const int C_max) { 
			//setup(prob, params, C_max); 
		}
		virtual ~DTP_Model() {}

		virtual int get_lb();
		virtual int get_ub();   
		// virtual void setup(Instance& prob, ParameterList *params, const int C_max); 
		virtual Variable get_objective_var();
		virtual int  get_objective();
		// virtual double  get_normalized_objective() {return (double)(get_objective())/(double)(disjuncts.size());}
		virtual int  set_objective(const int obj);
		
		virtual void setup_objective();
	};



  
	class SolutionPool {
	public:
		std::vector<SchedulingSolution*> pool_;
    
		SolutionPool() {}
		virtual ~SolutionPool() {
			for(unsigned int i=0; i<pool_.size(); ++i)
				delete pool_[i];
		}

		SchedulingSolution*& operator[](const int i) { return pool_[i]; }
		void add(SchedulingSolution *s) { pool_.push_back(s); }
		SchedulingSolution* getBestSolution() { return pool_.back(); }
		unsigned int size() { return pool_.size(); }

	};








}
