 
#include <mistral_sol.h>
#include <iostream>
#include <fstream>
#include <string>

//_DEBUGSCHED = 1

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
    std::vector<bool>              types;

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
  class SchedulingModel;
  class ParameterList {
  public:

    static const int LEX     =  0;
    static const int ANTILEX = -2;
    static const int PROMISE =  1;
    static const int ANTI    = -1;
    static const int GUIDED  =  2;
    static const int RGUIDED =  3;
    static const int RAND    =  4;

    static const int nia = 21;
    static const char* int_ident[nia];
    
    static const int nsa = 12;
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
    
    int PrintSolution; // whether the solution should be printed

    double Factor;
    double Decay;
    double Skew;

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

    double getSkew();

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
    int*    jsp_latecost;
    int*    jsp_earlycost;
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

    bool hasEarlyCost() const {return jsp_earlycost != NULL; }
    bool hasLateCost() const {return jsp_latecost != NULL; }
    int  getJobEarlyCost(const int i) const {return jsp_earlycost[i];}
    int  getJobLateCost(const int i) const {return jsp_latecost[i];}

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


  class SchedulingModel : public CSP {
  public:
    
    Instance *data;

    int ub_Weight;
    
    int ub_C_max;
    int lb_C_max;

    int ub_L_sum;
    int lb_L_sum;

    int ub_Depth;
    int lb_Depth;

    std::vector< int > first_task_of_disjunct;
    std::vector< int > second_task_of_disjunct;


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

    SchedulingModel() : CSP() {}
    SchedulingModel(Instance& prob, ParameterList *params, const int C_max);
    virtual void setup(Instance& inst, ParameterList *params, const int C_max);
    virtual ~SchedulingModel();

    virtual int get_lb() = 0;
    virtual int get_ub() = 0;
    virtual VariableInt* get_objective_var() = 0;
    virtual int  get_objective() = 0;
    virtual double  get_normalized_objective() = 0;
    virtual int  set_objective(const int obj) = 0;

    std::ostream& printStats(std::ostream& os);
  };


  class C_max_Model : public SchedulingModel {
  public:

    C_max_Model() : SchedulingModel() {}
    C_max_Model(Instance& prob, ParameterList *params, const int C_max) { setup(prob, params, C_max); }
    virtual ~C_max_Model() {}

    virtual int get_lb();
    virtual int get_ub();    
    virtual VariableInt* get_objective_var();
    virtual int  get_objective();
    virtual double  get_normalized_objective() {return (double)(get_objective())/data->getNormalizer();}
    virtual int  set_objective(const int obj);
  };

  class L_sum_Model : public SchedulingModel {
  public:

    L_sum_Model() : SchedulingModel() {}
    L_sum_Model(Instance& prob, ParameterList *params, const int L_sum) { setup(prob, params, L_sum); }
    virtual ~L_sum_Model() {}

    virtual int get_lb();
    virtual int get_ub();    
    virtual VariableInt* get_objective_var();
    virtual int  get_objective();
    virtual double  get_normalized_objective() ;
    virtual int  set_objective(const int obj);
  };

  class Depth_Model : public SchedulingModel {

  public:

    Depth_Model() : SchedulingModel() {}
    Depth_Model(Instance& prob, ParameterList *params, const int Depth) { setup(prob, params, Depth); }
    virtual ~Depth_Model() {}

    virtual int get_lb();
    virtual int get_ub();    
    virtual VariableInt* get_objective_var();
    virtual int  get_objective();
    virtual double  get_normalized_objective() {return (double)(get_objective())/(double)(disjuncts.size());}
    virtual int  set_objective(const int obj);
  };

  class No_wait_Model : public C_max_Model {
    
  public:
    
    int type;
    VarArray gen_disjuncts;

    Vector<int> job_size;
    Vector<int> job_index;
 
    No_wait_Model() : C_max_Model() {} 
    No_wait_Model(Instance& prob, ParameterList *params, const int C_max, const int t=1) : C_max_Model() { type=t; setup(prob, params, C_max); }
    virtual void setup(Instance& prob, ParameterList *params, const int C_max);
    virtual ~No_wait_Model() {}
  };

  class DTP_Model : public SchedulingModel {

  public:

    DTP_Model() : SchedulingModel() {}
    DTP_Model(Instance& prob, ParameterList *params, const int C_max) { 
      setup(prob, params, C_max); 
    }
    virtual ~DTP_Model() {}

    virtual int get_lb();
    virtual int get_ub();   
    virtual void setup(Instance& prob, ParameterList *params, const int C_max); 
    virtual VariableInt* get_objective_var();
    virtual int  get_objective();
    virtual double  get_normalized_objective() {return (double)(get_objective())/(double)(disjuncts.size());}
    virtual int  set_objective(const int obj);
  };

  class Solution {
  public:
    int *earlybool_value;
    int *latebool_value;
    int *task_min;
    int *task_max;
    int *ltask_value;
    int *disjunct_value;
    int *search_value;
    int *all_value;
    int min_objective_value;
    int max_objective_value;

  public:
    SchedulingModel  * model;
    SchedulingSolver *solver;
    
    Solution(SchedulingModel *m, SchedulingSolver *s, int target=-1);
    virtual ~Solution();

    int& operator[](const VariableInt *x) {return all_value[x->id];}

    double distance(Solution* s);

    void guide_search();
    void guide_search_bounds();

    std::ostream& print(std::ostream& os, std::string type="jsp");
    
  };

  
  class SolutionPool {
  public:
    std::vector<Solution*> pool_;
    
    SolutionPool() {}
    virtual ~SolutionPool() {
      for(unsigned int i=0; i<pool_.size(); ++i)
	delete pool_[i];
    }

    Solution*& operator[](const int i) { return pool_[i]; }
    void add(Solution *s) { pool_.push_back(s); }
    Solution* getBestSolution() { return pool_.back(); }
    unsigned int size() { return pool_.size(); }

  };


class StoreStats : public SolutionMethod {

protected:
  StatisticList *stats;
  SchedulingSolver *ss;

public:

  StoreStats(Solver *s, StatisticList *t) : SolutionMethod(s) 
  {
    stats = t;
    ss = (SchedulingSolver*)s;
  }

  virtual ~StoreStats() 
  {
  }
  
  virtual void execute();

  virtual void initialise() 
  {
  }
};

class SolutionGuidedSearch : public StoreStats {

protected:
  SolutionPool *pool;

public:

  SolutionGuidedSearch(Solver *s, SolutionPool* p, StatisticList *t) 
  : StoreStats(s,t) 
  {
    pool = p;
  }

  virtual ~SolutionGuidedSearch() 
  {
  }
  
  virtual void execute() ;
//   { 
//     SchedulingSolver *ss = (SchedulingSolver *)solver;
//     pool->add(new Solution(ss->model, solver));
//     if(pool->size()) pool->getBestSolution()->guide_search();
//     StoreStats::execute();
//   }

  virtual void initialise() 
  {
    if(pool->size()) pool->getBestSolution()->guide_search();
  }

};

// class SolutionGuidedSearch : public SolutionMethod {

// protected:
//   SolutionPool *pool;

// public:

//   SolutionGuidedSearch(Solver *s, SolutionPool* p) : SolutionMethod(s) 
//   {
//     pool = p;
//   }

//   virtual ~SolutionGuidedSearch() 
//   {
//   }
  
//   virtual void execute() 
//   { 
//     if(pool->size()) pool->getBestSolution()->guide_search();
//   }

//   virtual void initialise() 
//   {
//     execute();
//   }

// };




  class SchedulingSolver : public Solver {
  public:

    //int lower_bound;
    //int upper_bound;
    ParameterList  *params;
    StatisticList   *stats;
    SchedulingModel *model;
    SolutionPool     *pool;
    WeighterRestartGenNogood *nogoods;


    SchedulingSolver(SchedulingModel* m, 
		     ParameterList* p,
		     StatisticList* s);//  : Solver(*m,m->SearchVars), params(p), stats(s) 
//     { 
//       model = m; 
//       lower_bound = m->get_lb();
//       upper_bound = m->get_ub();

//       s.lower_bound = lower_bound;
//       s.upper_bound = upper_bound;

//       pool = new SolutionPool();
//       //stats = s; //new StatisticList();

//       addHeuristic( params.Heuristic, params.Randomized, params.IValue );
//     }
    virtual ~SchedulingSolver() {
      delete pool;
    }


    virtual int virtual_iterative_dfs();

    std::ostream& print_weights(std::ostream& os);
    void decay_weights(const double decay);

    void addObjective() { goal = new MinimiseVar(this, model->get_objective_var()); }
    void removeObjective() { delete goal; goal = NULL; }
    void addHeuristic( std::string Heu, const int rdz, std::string vo, const int hlimit ) {
      

      int val_ord = ParameterList::GUIDED;
      if(  vo == "rguided")   val_ord = ParameterList::RGUIDED;
      if(  vo == "lex"    )   val_ord = ParameterList::LEX;
      if(  vo == "antilex")   val_ord = ParameterList::ANTILEX;
      if(  vo == "rand"   )   val_ord = ParameterList::RAND;
      if(  vo == "promise")   val_ord = ParameterList::PROMISE;
      if(  vo == "anti"   )   val_ord = ParameterList::ANTI;

      //std::cout << val_ord << std::endl;


      if(val_ord == ParameterList::LEX ||
	 val_ord == ParameterList::ANTILEX ||
	 val_ord == ParameterList::RAND)
      for(int i=0; i<length; ++i) {
	delete sequence[i]->branch;
	if(val_ord == ParameterList::LEX) {
	  sequence[i]->branch = new ValSelectorMin(sequence[i]);
	} else if(val_ord == ParameterList::ANTILEX) {
	  //std::cout << "HERE" << std::endl;
	  sequence[i]->branch = new ValSelectorMax(sequence[i]);
	} else {
	  sequence[i]->branch = new ValSelectorRandMinMax(sequence[i]);
	}
      }

      //exit(1);


      if( Heu == "dom" ) {
	MinDomain h(abs(rdz));
	add( h );
      }
      if( Heu == "lex" ) {
	Lexicographic h;
	add( h );
      } 
      else if( Heu == "deg") {
	MaxDegree h(abs(rdz));
	add( h );
      } 
      else if( Heu == "rand") {
	Random h;
	add( h );
      } 
      else if( Heu == "dom+deg") {
	MinDomMaxDeg h(abs(rdz));
	add( h );
      } 
      else if( Heu == "domodeg") {
	DomOverDeg h(abs(rdz));
	add( h );
      } 
      else if( Heu == "domowldeg") {
	DomOverWLDeg h(abs(rdz));
	add( h );
      }
      else if( Heu == "domowdeg") {
	DomOverWDeg h(abs(rdz));
	add( h );
      }
      else if( Heu == "neighbor") {
	Neighbor h(abs(rdz));
	add( h );
      } 
      else if( Heu == "impact") {
	Impact h(abs(rdz));
	add( h );
      }
      else if( Heu == "impactodeg") {
	ImpactOverDeg h(abs(rdz));
	add( h );
      }
      else if( Heu == "impactowdeg") {
	ImpactOverWDeg h(abs(rdz));
	add( h );
      }
      else if( Heu == "impactowldeg") {
	ImpactOverWLDeg h(abs(rdz));
	add( h );
      }
      else if( Heu == "pfsp") {
	PFSP h(abs(rdz), val_ord);
	add( h );
      }
      else if( Heu == "osp") {
	OSP h(model, abs(rdz), val_ord);
	add( h );
      }
      else if( Heu == "osp-dw") {
 	OSP h(model, abs(rdz), val_ord, OSP::DOMAIN_P_TWEIGHT);
 	add( h );
      }
      else if( Heu == "osp-d") {
	OSP h(model, abs(rdz), val_ord, OSP::DOMAIN_O_NOT);
	add( h );
      }
      else if( Heu == "osp-b") {
	OSP h(model, abs(rdz), val_ord, OSP::DOM_O_BOOLWEIGHT);
	add( h );
      }
      else if( Heu == "osp-t") {
	OSP h(model, abs(rdz), val_ord, OSP::DOM_O_TASKWEIGHT);
	add( h );
      }
      else if( Heu == "osp-bt") {
	OSP h(model, abs(rdz), val_ord, OSP::DOM_O_BOOLTASKWEIGHT);
	add( h );
      }
      else if( Heu == "osp-dr") {
	OSP h(model, abs(rdz), val_ord, OSP::DOMAIN_O_NOTTYPE);
	add( h );
      }
      else if( Heu == "osp-br") {
	OSP h(model, abs(rdz), val_ord, OSP::DOM_O_BOOLWEIGHTTYPE);
	add( h );
      }
      else if( Heu == "osp-tr") {
	OSP h(model, abs(rdz), val_ord, OSP::DOM_O_TASKWEIGHTTYPE);
	add( h );
      }
      else if( Heu == "osp-btr") {
	OSP h(model, abs(rdz), val_ord, OSP::DOM_O_BOOLTASKWEIGHTTYPE);
	add( h );
      }
      else if( Heu == "osp-jt") {
	OSP h(model, abs(rdz), val_ord, OSP::DOM_O_TASKWEIGHTPJOB);
	add( h );
      }
      else if( Heu == "osp-jtl0") {
	OSP h(model, abs(rdz), val_ord, OSP::JTL0);
	add( h );
      }
      else if( Heu == "now-d") {
	NOW h(model, abs(rdz), NOW::DOM);
	add( h );
      }
      else if( Heu == "now-dob") {
	NOW h(model, abs(rdz), NOW::DOM_O_BOOLWEIGHT);
	add( h );
      }
      else if( Heu == "now-dot") {
	NOW h(model, abs(rdz), NOW::DOM_O_TASKWEIGHT);
	add( h );
      }
      else if( Heu == "now-dobt") {
	NOW h(model, abs(rdz), NOW::DOM_O_BOOLTASKWEIGHT);
	add( h );
      }
      else if( Heu == "now-dtt") {
	NOW h(model, abs(rdz), NOW::DOM_P_TWEIGHT);
	add( h );
      }
     else if( Heu == "now-dtb") {
	NOW h(model, abs(rdz), NOW::DOM_P_BWEIGHT);
	add( h );
      }
     else if( Heu == "now-dtbt") {
	NOW h(model, abs(rdz), NOW::DOM_P_BTWEIGHT);
	add( h );
      }
      else if( Heu == "job") {
	JOB h( model );
	add( h );
      }
      else if( Heu == "osp-now") {
	OSP h(model, abs(rdz), val_ord, OSP::NOW);
	add( h );
      }
//       else if( Heu == "osp-w") {
// 	OSP h(abs(rdz), val_ord, OSP_Type);
// 	add( h );
//       }
      else {
	NoOrder h;
	add( h );
      }

      if(heuristic) heuristic->limit = hlimit;
    }
    

    bool probe_ub();
    void jtl_presolve();
    void old_jtl_presolve();
    void dichotomic_search();
    void all_solutions_search();
    void branch_and_bound();
    void large_neighborhood_search();
    void extract_stable(List& neighbors, Vector<VariableInt*>& stable);
    void repair(Solution *sol, Vector<VariableInt*>& stable);


    void print_solution(std::ostream& os, std::string tp);
  };



}
