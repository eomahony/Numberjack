/* classe des problèmes de satisfaction de contraintes CSP*/
/** Finite domain CSP class */
class CSProblem : public OpProblem
{ public :

/* nombre de contraintes */
/** constraint number */
  int nbconst;

/* tableau des domaines : chaque domaine est un vecteur d'entiers */
/** domain array : each domain is implemented by a vector of integers */
  vector<int>* tabdomains;  
/* pour chaque variable, numéro de son domaine : indice dans le tableau tabdomains */
/** for each variable, domain number : index in tabdomains array */
  int* domains;
/* tableau des connexions : pour chaque variable, vecteur des variables connectées */
/** connections table : for each variable, vector of connected variables */
  vector<int>* connections;
/* constructeur de base */
/** constructor */
  CSProblem (int nbvar, int nbconst);
/* constructeur avec borne inférieure */
/** constructor with lower bound (stopping condition when it is reached) */
  CSProblem (int nbvar, int nbconst, int lower);
  ~CSProblem();
  void move_execution(Configuration* configuration, Move* move);
/* la taille du domaine de la variable var */
/** the domain size of variable var */
  virtual int variable_domainsize (int var);
  void random_configuration(Configuration* configuration);
/* une variable choisie aléatoirement */
/** a variable randomly chosen */
  virtual int random_variable(Configuration* configuration);
/** a variable taking part to a conflict in the configuration */
  virtual int random_conflict_variable(Configuration* configuration);
/* une valeur choisie aléatoirement, si possible distincte de val : retourne l'indice de la valeur dans le domaine */
/** a value for variable var, randomly chosen in its domain, if possible distinct with val : returns the index of the value in the domain */
  virtual int random_value(int var,int val);
/* une valeur dans le domaine de la variable minimisant les conflits avec la configuration */
/** a value in the domain minimizing the conflict with the configuration (implementation of Minton min-conflict heuristics) 
returns the index of the value in the domain */
  virtual int min_conflict_value(int var,int val, Configuration * configuration);
/* initialisation des domaines par défaut : un seul domaine numéro 0 pour toutes les variables */
/** standard domain initialization : a unique domain number 0 for all variables */
  virtual void init_domains(int nbvar, int s );
/* un seul domaine par défaut : entiers de 0 à s-1 */
/** standard unique domain : integers from 0 to s-1 */
  virtual void init_tabdomains(int s);
/* calcul des variables en conflit : on reconstruit le vecteur des variables en conflit d'une configuration*/
/** compute the variables in conflict : rebuilding the vector of conflict variables of the configuration*/
  void compute_var_conflict(Configuration* configuration);
  void best_config_analysis();
  void best_config_write();
  Long move_evaluation(Configuration* configuration,Move* move);
  void init_population (Configuration** population,int populationsize);
  Configuration* create_configuration();
  Move* create_move();
  void adjust_parameters(Configuration* configuration, int & maxneighbors, int & minneighbors);
  void next_move(Configuration* configuration, Move* move, NeighborhoodSearch* nbhs);	
/* met en place les domaines et connexions d'un problème */
/** set the domains and connections of a problem */
  virtual  void set_domains_connections( int* dom, vector<int>* tabledom, vector<int> * connect );
/* initialisation des domaines : appel de init_domains et init_tabdomains */
/** initialization of the domains : call init_domains and init_tabdomains */
  virtual void init_domain_tabdomain();
  int tabuindex(Move* move, Configuration* config);
  int tabuinverseindex(Move* move, Configuration* config);
  int nbtabuindex();
};


/* CSP Binaires :  ajout du tableau des  contraintes à partir de 2 variables */
/** Binary CSPs : addition of the constraints array */ 
class BinaryCSProblem : public CSProblem
{ public :
/* pour une paire de variables (i,j) (i<j) , constraints[i][j] contient le numéro de contraintes +1 entre ces variables si
elles sont connectées, 0 sinon. On se limite à au plus une contrainte par paire de variables : dans le 
cas contraire on peut utiliser la classe  WeightExtensionBinaryCSP */
/** for a couple (i,j) of variables, (i<j) , constraints[i][j] returns the constraint number + 1 if the variables are connected, 0 si the variables are not connected. It is assumed that at most one constraint exists between two variables
(if not use WeightExtensionBinaryCSP class) */
  int** constraints;
  BinaryCSProblem( int nbvar, int nbconst);
  BinaryCSProblem (int nbvar, int nbconst, int lower);
  ~BinaryCSProblem() {;};
};




int** csp_constraintdatastructure(int nbvar);
