/* Les définitions des classes de la partie algorithme + OpProblem */
/** the definitions of classes of the algorithmic part + OpProblem */

#include <climits>
typedef long long Long;
#ifdef WINDOWS
        inline double drand48() {double(rand() / RAND_MAX); }
        inline void srand48(long seed) {return srand48(seed);}


//replace srand48(seed) by srand(seed)
//replace drand48() by (double(rand()) / RAND_MAX)


#endif



/* struct Long */
/* { */
/*   long long p; */

/*   Long() : p(0) {} */
/*   Long(const Long &l) : p(l.p) {} */
/*   Long(const long long v) : p(v) {} */
/*   double to_double() const {return (double) p;} */

/*   Long &operator=(const Long &r) { */
/* 	p = r.p; */
/* 	return *this; */
/*   } */
/*   Long &operator+=(const Long &r) { */
/* 	p += r.p; */
/* 	return *this; */
/*   } */
/*   Long &operator-=(const Long &r) { */
/* 	p -= r.p; */
/* 	return *this; */
/*   } */
/*   const Long operator-() const {return Long(-p);} */
/*   friend const Long operator+(const Long& left, const Long& right) { */
/* 	return Long(left.p + right.p); */
/*   } */
/*   friend const Long operator-(const Long& left, const Long& right) { */
/* 	return Long(left.p - right.p); */
/*   } */
/*   friend const Long operator*(const Long& left, const Long& right) { */
/* 	return Long(left.p * right.p); */
/*   } */
/*   friend const Long operator/(const Long& left, const Long& right) { */
/* 	return Long(left.p / right.p); */
/*   } */
/*   friend bool operator==(const Long& left, const Long& right) { */
/* 	return (left.p == right.p); */
/*   } */
/*   friend bool operator!=(const Long& left, const Long& right) { */
/* 	return (left.p != right.p); */
/*   } */
/*   friend bool operator<=(const Long& left, const Long& right) { */
/* 	return (left.p <= right.p); */
/*   } */
/*   friend bool operator>=(const Long& left, const Long& right) { */
/* 	return (left.p >= right.p); */
/*   } */
/*   friend bool operator<(const Long& left, const Long& right) { */
/* 	return (left.p < right.p); */
/*   } */
/*   friend bool operator>(const Long& left, const Long& right) { */
/* 	return (left.p > right.p); */
/*   } */
/*   friend ostream& operator<<(ostream& os, const Long &r) { */
/*     os << r.p; */
/* 	return os; */
/*   } */
/*   friend istream& operator>>(istream& is, Long& r) { */
/* 	is >> r.p; */
/* 	return is;  */
/*   } */
/* }; */
//#define LONG_MAX LONG_LONG_MAX

/* les classes "abstraites" utilisées dans les paramètres des méthodes */
class OpProblem;
class IncompleteAlgorithm;
class Metaheuristic;
class NeighborhoodSearch;
class Move;



/* la classe Configuration  le champ config comprend la configuration elle-même sous forme de tableau d'entiers
le champ valuation contient sa valeur pour l'évaluation */

/** the main class Configuration */

class Configuration 

{public :
  int nbvar;
 int trynumber;
/* les valeurs courantes des variables : implanté sous forme de tableau d'entiers*/
/** the current values of the variables : implemented with an array of integers*/
  int* config;
/* la valeur de la configuration */
/** the configuration value */
  Long valuation;
  int var_conflict_size;
/* les variables participant à un conflit : implanté sous forme de vecteur */
/** the variables taking part to a conflict : implemented with a vector */
   vector<int> var_conflict;	
   set<int> set_var_conflict;
/* indicateur si la configuration a été regroupée (pour GWW) */
/** indicates if the configuration has been regrouped before (for GWW) */
  int regrouped;
  virtual ~Configuration();
  Configuration();
  Configuration(int nbvar);
/* copie d'une configuration config2 dans this*/
/** copy a configuration config2 into this */
  virtual void copy_element (Configuration* config2);
/* initialisation à 0 de la structure de données des conflits */
/** initialization to 0 of the conflict datastructure */
  virtual void init_conflicts ();
/* stockage de l'augmentation des conflits de (var,val) de incr */
/** store the conflict of (var,val) incremented by incr */
  virtual void incr_conflicts (int var, int val, int index, Long incr);

/* stockage du nombre des conflits nbconf de (var,val)  */
/** store the number of conflicts nbconf of (var,val) in the conflict datastructure */
  virtual  void set_conflicts (int var, int val, int index, Long nbconf);

/* nombre de conflits de (var,val) stocké */
/** get the number of conflicts (var,val) stored in the conflict datastructure*/
  virtual  Long get_conflicts (int var, int val, int index);
/* nombre de conflits de (var,val) , au besoin recalculé */
/** get the number of conflicts of (var,val), computed if not stored */
  virtual  Long get_conflicts_problem (OpProblem* problem, int var, int val);

/* mise à jour des conflits après avoir effectué le mouvement move*/
/** update the conflict datastructure after a move is done */
  virtual void update_conflicts(OpProblem* problem, Move* move);
};

/* CSPConfiguration : pour les CSP */
/** CSPConfiguration : for the CSPs */
class CSPConfiguration: public Configuration
{
 public :
 int domainsize;

 CSPConfiguration(int nbvar, int domsize);
};

/* L'incrémentalité avec stockage de la participation à l'évaluation des valeurs courantes des 
variables de la configuration : implanté dans tabconflicts (tableau  à une dimension) */
/** Incremental evaluation with storage in the conflict datastructure
tabconflicts the participation of the current values of the configuration */
class IncrCSPConfiguration : public CSPConfiguration
{ public : 
  Long* tabconflicts;
  IncrCSPConfiguration (int nbvar);
  IncrCSPConfiguration(int nbvar, int nbcol);
  ~IncrCSPConfiguration();
  void copy_element (Configuration* config2);
  void init_conflicts ();
  void incr_conflicts (int var, int val , int index, Long incr);
  void set_conflicts (int var, int val, int index,  Long nbconf);
  Long get_conflicts (int var, int val , int index);
  Long get_conflicts_problem (OpProblem* problem, int var, int val);
  virtual  void set_variableconflicts (int var, int nbconf);
  void update_conflicts(OpProblem* problem, Move* move);		
};

/* l'incrémentalité totale : participation à l'évaluation de chaque
valeur de chaque variable stockée dans le  tableau tabconflicts à deux dimensions (variable, indice de la valeur)*/
/** Full incremental evaluation : the participation of every value of every
variable is stored in the 2 dimension array tabconflicts (variable, valueindex)
*/
class FullincrCSPConfiguration : public CSPConfiguration
{ public :
  int tabconflictsize;
  Long** tabconflicts;

  FullincrCSPConfiguration(int nbvar, int domainsize);
  ~FullincrCSPConfiguration();
  void copy_element (Configuration* config2);
  void init_conflicts ();
  void incr_conflicts (int var, int val , int index, Long incr);
  void set_conflicts (int var, int val, int index, Long nbconf);
/* nombre de conflits de (var,val) stocké : utilisation de l'indice de la valeur index*/
/** get the number of conflicts (var,val) stored in the conflict datastructure using the value index in the domain */
  Long get_conflicts (int var, int val, int index);
  Long get_conflicts_problem (OpProblem* problem, int var, int val);
  void update_conflicts(OpProblem* problem, Move* move);
};

/* classe Move */
/** root class Move */
class Move 
{public :
 Long valuation;
 Move();
 virtual ~Move() {;};
/* le test d'égalité d'un mouvement (utilisé pour la recherche d'un mouvement dans la liste taboue)*/
/** the test of equality of a move (used for searching a move in the tabu list) */
 virtual int eqmove(Move* move1);
/* copie du mouvement move1 dans this */
/** copy of move move1 into this */
 virtual void copymove(Move* move);
/* le mouvement a mettre dans la liste taboue */
/** the move to be put in the tabu list (to be implemented in the subclasses)*/
 virtual Move* computetabumove(Configuration* config){return 0;};
};

/* classe CSPMove :  un mouvement pour les CSP : variable , valeur */
/** class CSPMove : a classical move for a CSP : variable, value */
class CSPMove : public Move
{public :
 int variable;
 int value;
 CSPMove();
 ~CSPMove() {;};
 int eqmove(Move* move);
 void copymove (Move* move);
/* le mouvement stocké tabou est le mouvement inverse du mouvement effectué */
/** the move stored is the inverse of the move done */
 Move* computetabumove(Configuration* config);
};


/* classe racine des problèmes d'optimisation (minimisation) */
/** Root class of Optimization Problems (minimization) */

class OpProblem 
{public :
/* la meilleure configuration trouvée */
/** the best configuration found */
  Configuration* best_config;
/* nombre de variables */
/** the number of variables */
  int nbvar;  
/* taille maximum des domaines */
/** maximum domain size */
  int domainsize;
/* borne inférieure donnée au départ : sert de condition d'arrêt quand elle est atteinte */
/** given lower bound , is used as a stop condition when it is reached */
  Long lower_bound;
/* le mouvement courant */
/** the current move being tested */
  Move* currentmove;
/* le premier mouvement faisable essayé dans le voisinage*/
/** the first feasible move tried in the neighborhood */
  Move* firstmove;
/* le meilleur mouvement essayé */
/** the best move found in the neighborhood */
  Move* bestmove;
  OpProblem(){};
  virtual ~OpProblem(){};
/* exécution d'un mouvement (modification de la configuration courante) */
/** move execution (modification of the current configuration) */
  virtual void move_execution(Configuration* configuration, Move* move);
/* mise à jour de la structure des conflits (cas IncrCSPConfiguration) */
/** update of  the conflict data structure (case IncrCSPConfiguration) */
  virtual void incr_update_conflicts(IncrCSPConfiguration* configuration,Move* move){};
/* mise à jour de la structure des conflits (cas FullincrCSPConfiguration) */
/** update of  the conflict data structure (case FullincrCSPConfiguration) */
  virtual void fullincr_update_conflicts(FullincrCSPConfiguration* configuration,Move* move){};
/* création des 3 objets Move (currentmove,bestmove,firstmove) */
/** creation of 3 Move objects (currentmove,bestmove,firstmove) */
  virtual void allocate_moves();
/* création d'un mouvement (la classe du mouvement dépend du problène) : méthode implantée dans les sous-classes */
/** creation of 1 Move object (the class of the Move depends on the problem) : this method
is implemented in subclasses */
  virtual Move* create_move(){return 0;};
/* ajustement des paramètres du voisinage (quand la taille du voisinage est supérieure à maxneighbors) */
/** adjustment of the neighborhood parameters (when the size of the actual neighborhood is greater than maxneighbors) */
  virtual void adjust_parameters(Configuration* configuration, int & maxneighbors, int & minneighbors){};
/* prochain mouvement du voisinage à tester */
/** next move to be tested (implemented in subclasses)*/
  virtual void next_move(Configuration* configuration, Move* move, NeighborhoodSearch* nbhs){};	
/* affectation aléatoire des variables d'une configuration */
/** random assignment of the variables of a configuration */
  virtual void random_configuration(Configuration* configuration){};
/* analyse da la meilleure solution */
/** analysis of the best configuration */
  virtual void best_config_analysis(){};
/* ecriture de la meilleure solution */
/** writing the best solution */
  virtual void best_config_write(){};

/* vérification de la meilleure solution (recalcul de son coût) */
/** verification of the best solution (its cost is recomputed) */
  virtual void best_config_verification();
/* initialisation d'une population de taille populationsize */
/** initialization of the population of size populationsize */
  virtual void init_population (Configuration** population,int populationsize) {};
/* création d'une configuration (la classe exacte dépend du problème) */
/** create a configuration (the exact class depends on the problem and must defined in subclasses) */
  virtual Configuration* create_configuration(){return 0;};
/* calcul de la participation à l'évaluation de l'affectation (var,val) */
/** computation of the participation of (var,val) to the configuration evaluation */
  virtual Long compute_conflict (Configuration* configuration,int var, int val) {return 0;};
/* évaluation d'une configuration */
/** evaluation of a configuration */
  virtual Long config_evaluation(Configuration* configuration) {return 0;};
/* évaluation d'un mouvement move sur une configuration */
/** evaluation of a configuration if the move is done */
  virtual Long move_evaluation (Configuration* configuration,Move* move){return 0;};
/* passage de l'indice dans le domaine à la valeur */
/** valueindex in the domain to value */
  virtual int index2value (int index, int var) {return index;};
/* passage d'une valeur à son indice dans le domaine de la variable */
/** valueindex of value in its domain*/
  virtual int value2index(int value,int var) {return value;};
/* calcule l'ensemble des variables en conflit de la configuration*/
/** compute the variables participating to a conflict in the configuration */
  virtual void compute_var_conflict (Configuration* configuration) {};
  virtual int tabuindex(Move* move, Configuration* configuration) {return 0;};
  virtual int tabuinverseindex(Move* move, Configuration* configuration){return 0;};
  virtual int nbtabuindex(){return 0;};
};
  
/* Le voisinage paramétré d'une part par  min voisins explorés, 
max voisins explorés et épuisement voisinage et d'autre part par var_conflict et val_conflict */
/** Class NeighborhoodSearch : how the neighborhood is explored */

class NeighborhoodSearch
{public :
/* nombre minimum de voisins explorés */
/** minimum number of visited neighbors */
  int minneighbors;
/* nombre maximum de voisins explorés */
/** maximum number of explored neighbors */
  int maxneighbors;
/* indicateur de comportement quand le voisinage est épuisé sans qu'un voisin n'ait été accepté :
0 stagnation, 1 on effectue le 1er mouvement faisable, k on effectue le meilleur mouvement faisable parmi k mouvements essayés non acceptés*/
/** behavior indicator when the neighborhood is exhausted and no neighbor has been accepted :
0 stagnation, 1 the 1st feasible move is selected, k the best feasible among k tried but not accepted moves is selected */
  int finished;
/* indicateur de restriction aux variables en conflit (0 pas de restriction, 1 restriction) */
/** restriction indicator to variables participating in a conflict (0 no restriction, 1 restriction) */
  int var_conflict;
/* indicateur de restriction aux meilleures variables d'une variable (0 pas de restriction, 1 restriction) */
/** restriction indicator to best values of a variable (0 no restriction , 1 restriction) */
  int val_conflict;
  double nbhrate;
  NeighborhoodSearch( int maxneigh, int minneigh, int finish, int var_conf, int val_conf, double nbbr);
  int returnbestmove();
  void adjust_neighborhood(Configuration* configuration, OpProblem* problem, int& maxneigh, int& minneigh, int nbmoves);
  virtual void dynamicmaxneighbors(int & maxneigh, int & minneigh, int nbmoves);
  virtual void initsearch();
  virtual void spareneighboradjust(Configuration* config, Move* move){;}
};


/* Voisinage avec réglage dynamique du paramètre max-voisins*/
/** Neighborhood with dynamic parameter tuning */
class DynamicNeighborhoodSearch: public NeighborhoodSearch
{public:
 DynamicNeighborhoodSearch(int maxneigh, int minneigh, int finish, int var_conf, int val_conf, double nbbr);
/* valeur initiale du parametre maxneighbors */
/** initial value of maxneighbors parameter */
 int initmaxneighbors;
/* valeur initiale du parametre minneighbors */
/** initial value of minneighbors parameter */
 int initminneighbors;
/* période de réajustement du paramètre */
/** parameter readjustment period */
 int adjustperiod;
 void initsearch();
/* ajustement des paramètres minneighbors et maxneighbors */
/** adjust the parameters maxneighbors and minneighbors */
 void dynamicmaxneighbors(int & maxneigh, int & minneigh, int nbmoves);
};

class DynamicSpareneighbor: public NeighborhoodSearch
{public :
DynamicSpareneighbor (int maxneigh, int minneigh, int finish, int var_conf, int val_conf, double nbbr);
void spareneighboradjust(Configuration* config, Move* move);
int nbmovesdown;
};

/* Les Algorithmes
la classe mere : algo de recherche incomplet */
/** Root class of algorithms */

class IncompleteAlgorithm
{public :
  string methodname;
/* un seuil peut être utilisé pour empêcher des mouvements de coût supérieur au seuil
(utilisé dans les recherches locales des marches de GWW)*/
/** a threshold can be used to forbid moves above this threshold (used in LSAlgorithms implementing walks inside GWW)*/
  Long threshold;
  virtual ~IncompleteAlgorithm(){};
/* marche d'une particule */
/** walk for a particule */
  virtual void randomwalk (OpProblem* problem, Configuration* configuration);
  virtual void initthreshold(Configuration** population, int popsize){;};
/* exécution de l'algorithme sur une population (réduite à une particule pour une recherche locale) */
/** Run the algorithm on a population (array of configurations) */
  virtual void run (OpProblem *problem,Configuration ** population);

};

/* la classe des algos de marche aléatoire paramétrée par longueur marche
un voisinage et une metaheuristique */
/** The class of local search algorithm on one particle : the random walk is 
parameterized with the walk lengh,a neighborhood and a metaheuristics */

class LSAlgorithm: public IncompleteAlgorithm
{public :
/* longueur de la marche */
/** walk length */
  int walklength;
/* le voisinage */
/** the way the neighborhood is explored */
  NeighborhoodSearch * nbhsearch;  
/* la métaheuristique */
/** the metaheuristics used */
  Metaheuristic* mheur;
/* le nombre d'essais de mouvements (pour les stats) */
/** number of move tries  (for statistics) */
  int nhtries;
  double avgnhtries;
  double avgsqnhtries;
/* nombre de mouvements effectués */
/** number of moves done */
  int nbmoves;
  LSAlgorithm(int nbmov);
  ~LSAlgorithm();
/* faisabilité d'un mouvement  (sous ou au niveau du seuil pour marche de GWW) */
/** feasability of a move (under or at threshold level pour GWW walks) */
  virtual int isfeasible(Move* move);
  void randomwalk   (OpProblem* problem, Configuration* configuration);
/* algorithme d'exploration du voisinage pour sélectionner et effectuer un mouvement à partir de la configuration courante
Effectue le mouvement et renvoie 1 si un mvt a ete effectué et 0 si aucun mouvement ne l'a été*/
/** Neighborhood exploration algorithm for selecting and do a move from the current configuration :
returns 1 if a move has been done and 0 if no move has been done */
  virtual int configurationmove(OpProblem* problem,Configuration* configuration);
  void initthreshold(Configuration** population, int popsize);
  void run (OpProblem *problem, Configuration ** population);
/* test de meilleur trouvé (renvoie 1 si un meilleur absolu est trouvé)*/
/** test if a global best configuration has been found (returns 1 in that case) */
  int test_bestfound(Move* move);


};

class LSAlgorithmGWW: public LSAlgorithm
{public :
 LSAlgorithmGWW(int nbmov);
 int isfeasible(Move* move);
};

/* les différentes métaheuristiques */
/** Root class for Metaheuritics */
class Metaheuristic
{public :
  virtual ~Metaheuristic(){};
/* mise à jour des données de la métaheuristique juste avant qu'un mouvement soit effectué */
/** update of the metaheuristic data just before a move is performed */
  virtual void executebeforemove(Move* move, Configuration* configuration,OpProblem* problem);
/* initialisation des données de la métaheuristique */
/** initialization of the meteheuristic data at the beginning of a local search */
  virtual void reinit(OpProblem* problem);
/* condition d'acceptation d'un mouvement : renvoie 1 si le mouvement est accepté */
/** acceptance condition of a move : returns 1 if the move is accepted */
  virtual int acceptance(Move* move,Configuration*  config);
  virtual void adjustparameter (int parameter) {;};
};

/* marche avec liste taboue : parametree par longueur de la liste : cette liste de mouvements est 
implantee à l'aide d'une liste de Move* */
/** Walk with using a tabu list : this list of moves is implemented by a list<Move*> structure , the
actual class of the moves depend on the problems */
class TabuSearch: public Metaheuristic
{public :
/* longueur maximale de la liste taboue */
/** maximum length of the tabulist */
  int tabulength;
/* liste taboue : traitée comme une file */
/** tabu list : implemented FIFO */
  list<Move*> move_list;
  TabuSearch(int tabul);
/* acceptation d'un mouvement : non tabou  (le critère d'aspiration est dans l'algo de recherche du voisin) */
/** acceptance of a move : not in the tabulist (the aspiration criterion of a best is in the configurationmove algorithm) */
  int acceptance (Move* move, Configuration* config);
/* test de non présence dans la liste taboue : la présence d'un mvt est faite avec eqmove */
/** test of non presence in the tabulist (use of eqmove method) */
  int  nontabumove (Move* move);
/* mise à jour de la liste taboue qui est traitée comme une file de longueur maximale tabulength */
/** updating of the tabulist which is managed as a FIFO of maximum length tabulength */
  void executebeforemove(Move* move, Configuration* configuration, OpProblem* problem);
/* réinitialisation : la liste taboue est vidée */
/** the tabu list is cleared */
  void reinit(OpProblem* problem);
  void adjustparameter (int length);
};

class TabuGreedySearch: public TabuSearch
{public :
  TabuGreedySearch(int tabul);
 int acceptance (Move* move, Configuration* config);
};

class IncrTabuSearch: public TabuSearch
{public :
  IncrTabuSearch(int tabul);
  int nbiter;
  vector<int> tabutime;
  OpProblem* currentproblem;
  int  acceptance (Move* move, Configuration* config);
  void executebeforemove(Move* move, Configuration* configuration, OpProblem* problem);
  void reinit(OpProblem* problem);
};

class IncrTabuGreedySearch: public IncrTabuSearch
{public :
  IncrTabuGreedySearch(int tabul);
 int acceptance (Move* move, Configuration* config);
};

/* marche Metropolis : un seul paramètre = temperature */
/** Metropolis algorithm : a unique parameter - a constant temperature */
class Metropolis: public Metaheuristic
{public :
  double temperature;
  Metropolis(double temp);
/* la formule classique de Metropolis d'acceptation d'un mouvement détériorant
l'évaluation : probabilité p = exp (-evaluationdelta/temperature) */
/** the classical Metropolis formula for accepting a bad move :  probability =  exp (-evaluationdelta/temperature) */
  int acceptance (Move* move, Configuration* config);
  void adjustparameter(int parameter);
};

/* l'acceptation à seuil : un mouvement ne doit pas détériorer l'évaluation plus que le seuil courant ; 
le seuil diminue linéairement de thresholdinit à 0*/

/** Threshold accepting Metaheuristics : a move must no deteriorate the evaluation more than the
current threshod : the threshold goes down linearly from thresholdinit to 0 */
class ThresholdAccepting: public Metaheuristic
{public :
/* seuil initial */
/** initial threshold */
  double thresholdinit;
/* pas de baisse du seuil */
/** constant step to lower the threshold */
  double delta;
/* valeur courante du seuil */
/** current value of the threshold */
  double thresholdaccept; // le seuil tel que géré par TA
/* constructeur : 2 arguments seuil initial maxthreshold et nombre de pas, 
le pas constant delta de baisse du seuil est calculé*/
/** constructor : two arguments : maxthreshold the initial threshold and 
walklength , it computes a constant step for lowering the threshold */
  ThresholdAccepting(double maxthreshold, int walklength);
/* condition d'acceptation : être sous ou au niveau du  seuil */
/** acceptance condition : being under or at the threshold */
  int acceptance (Move* move,Configuration* config);
/* le seuil est diminué de delta */
/** the threshold is lowered by delta */
  void executebeforemove(Move* move, Configuration * configuration, OpProblem* problem);
/* le seuil est initialisé à thresholdinit */
/** the threshold is initialized at thresholdinit */
  void reinit(OpProblem* problem);
};

/* le recuit simulé : descente linéaire de température de inittemperature à 0 */
/** Simulated Annealing : linear temperature descent from inittemperature to 0*/
class SimulatedAnnealing: public Metaheuristic
{public :
/* temperature initiale */
/** initial temperature */
  double inittemperature;
/* pas constant de baisse de temperature */
/** constant step for lowering the temperature */ 
  double delta;
/* temperature courante */
/** current temperature */
  double temperature; 
  int walklength; 
/* Constructeur : 2 arguments : température initiale et longueur de marche */
/** Constructor : 2 parameters  : initial temperature and walk length : the fixed
temperature decrement is computed. */
  SimulatedAnnealing(double initialtemperature, int walklength);
/* acceptation en fonction de la temperature : formule classique du recuit simulé 
probablité d'acceptation d'un mouvement détériorant l'évaluation :
probabilité =  exp (-evaluationdelta/temperature) */

/** Acceptance function of the temperature : classical simulated annealing formula
for accepting a bad move :  probability =  exp (-temperature/evaluationdelta) */
  int acceptance (Move* move, Configuration* config);
/* la température est baissée de delta */
/** the temperature is lowered by delta */
  void executebeforemove(Move* move, Configuration * configuration, OpProblem* problem);
  void reinit(OpProblem* problem);
  void adjustparameter (int parameter);
};


/* marche hybride tabou + pourcentages d'acceptation selon sens des mouvements */
/** Special Tabu search with complementary acceptance condition depending on the move direction */


//                          liste taboue 

class TabuAcceptingrate: public TabuSearch
{public :
   /* probabilité d'acceptation d'un mauvais   */
  /** probability of acceptance of a worsening move */
  float Pd;       
  /* probabilité d'acceptatiion d'un mouvement de même coût que le courant */
  /** probability of acceptance of a move with same cost */
  float P0;        
  TabuAcceptingrate(int tabul, float Pd, float P0);
/* critère d'acceptation : non tabou et pourcentages d'acceptation suivant sens du mouvement (détériorant, neutre, améliorant) */
/** Acceptance condition : non tabu and probabilities depending on the move direction */
  int acceptance (Move* move, Configuration* config);
};


/* marche aléatoire : tout voisin faisable est accepté */
/** Random walk : every feasible neighbor is accepted */
class RandomSearch: public Metaheuristic
{public :
  RandomSearch();
  int acceptance (Move* move, Configuration* config);
};

/* marche gloutonne : on accepte un voisin de cout inférieur ou égal à la configuration courante*/
/** Greedy walk : a neighbor with better or same cost as the current configuration is accepted */
class GreedySearch: public Metaheuristic
{public :
  GreedySearch();
  int acceptance (Move* move, Configuration* config);
};

//-------------------------------------------------------------------------------------------------


/* les algos de type GWW 
 les différentes sous classes diffèrent par la gestion du seuil
et les regroupements de particules */

/** the GWW (Go with the winners) algorithms : the different subclasses 
differ by the way a threshold is managed and the particles are regrouped */

class GWWAlgorithm: public IncompleteAlgorithm
{public :
/* nombre de particules */
/** number of particles */
  int populationsize;
/* indicateur de marche uniquement si la particule a été regroupée 
(utile pour GWW de base, sans recherche locale, uniquement) (1 oui, 0 non) */
/** walk indicator : a walk is performed only is the particle has been regrouped : (1 yes, 0 no)
(useful for a standard GWW with random walk (and no local search)) */
  int regrouptest;
/* indicateur de baisse du seuil au dernier mouvement de la marche (pour essayer d'empecher la particule d' etre redistribuée) (1 oui, 0 non) */
/** parameter if the threshold is lowered at the last move of the walk 
(for trying to avoid the particle to be redistributed  (1 yes, 0 no)*/
  int lastmovedescent;
/* indicateur d'élitisme : remet-on le meilleur individu dans la population à chaque regroupement (1 oui, 0 non) */
/** elitism parameter : is the best particle put again in the population at each regroupment ( 1 yes, 0 no) */
  int elitism;
/* indicateur d'arrêt de la marche en cas de stagnation (1 oui, 0 non) */
/** parameter for stopping the walk in case of stagnation (1 yes, 0 no) */ 
  int nomovestop;
/* le décrément du seuil (calculé par thresholdcomputedelta) */
/** the threshold decrement (compted by thresholdcomputedelta) */
  Long thresholddelta;
/* le nombre d'iterations max : utile quand pas de seuil (NothresholdGWWAlgorithm) */
/** the maximum number of iterations : useful when no threshold is managed (NothresholdGWWAlgorithm) */
  int nbiteration;
/* le nombre de changements de seuil (pour les statistiques) */
/** number of threshold changes (for the statistics) */
  int thresholdchanges;
/* le nombre total d'essais de mouvements entre 2 regroupements (pour les statistiques)*/
/** total number of move tries between 2 regroupments (for the statistics) */
  int total_nhtries;
/* le nombre total de mouvements entre 2 regroupements (pour les statistiques)*/
/** total number of moves between 2 regroupments (for the statistics) */
  int total_nbmoves;
/* l'algorithme de recherche locale utilisé */
/** the local search algorithm used */
  LSAlgorithm* walkalgorithm;
/* destructeur */
  ~GWWAlgorithm();
/* recherche locale sur l'ensemble de la population */
/** local search on the whole population */
  virtual  void populationrandomwalk (OpProblem* problem, Configuration** population);
/* le nombre de particules au seuil (pour les statistiques),  la population étant déjà triée à l'appel */
/** the number of particles at the threshold (for statistics) , the population being yet sorted at the function call*/
  virtual int nb_threshold_population(Configuration** population);
/* une recherche locale pour une particule */
/** a local search for a particle */
  void randomwalk   (OpProblem* problem, Configuration* configuration);
/* initialisation du seuil */
/** intialization of the threshold */
  void initthreshold(Configuration** population, int popsize);
/* méthode de baisse du seuil (le delta a déjà été calculé)*/
/** method for lowering the threshold( the delta has already been computed) */
  virtual void thresholdupdate();
/* méthode de calcul du décrément du seuil */
/** method for computing the threshold decrement */
  virtual void thresholdcomputedelta(Configuration** population);
/* déroulement de l'algorithme */
/** main function for running the algorithm */
  void run (OpProblem *problem, Configuration** population);
/* regroupement des mauvais individus sur les bons */
/** regrouping of the best particles on the good ones */
  virtual  void regrouping(Configuration** population);
/* en cas d'élitisme, on remet le meilleur dans la population */
/** in case of elitism, the best particle is put into the population */
  void populationkeepbest (OpProblem* problem, Configuration** population);
/* incremente le compteur de changements de seuil (pour les statistiques) */
/** incrementing the threshold updates counter (for the statistics) */
  virtual  void thresholdchangesupdate();
};

/* Classe abstraite : GWW avec seuil */
/** Abstract class : GWW managing a threshold */
class ThresholdGWWAlgorithm : public GWWAlgorithm
{public :
  void thresholdupdate();
void thresholdchangesupdate();
void initthreshold(Configuration** population, int popsize);
int nb_threshold_population(Configuration** population);
};


/* GWW standard : descente du seuil avec taux fixe */
/** Standard GWW : threshold descent with a fixed rate */
class StandardGWWAlgorithm : public ThresholdGWWAlgorithm
{public :
/* taux de baisse du seuil */
/** threshold descent constant rate */
  double thresholddescent;
/* seuil minimum (correspond en général à une borne inférieure connue) */
/** minimum of the threshold (corresponds generally to a known lowerbound) */
  Long thresholdmin;  
  void regrouping(Configuration** population);
  StandardGWWAlgorithm(int population_size, int grtest,int lastmove, int elitisme,int stop, double thresdescent,Long threshmin );
  void thresholdcomputedelta(Configuration** population);
};

/* GWW descente du seuil au moins au niveau du pire */
/** StandardGWW with a threshold descent at least until the worst particle */
class FastStandardGWWAlgorithm: public StandardGWWAlgorithm
{public :
   void thresholdcomputedelta(Configuration** population);
   FastStandardGWWAlgorithm(int population_size, int grtest,int lastmove, int elitisme, int stop, double thresdescent,Long threshmin );	
};

/* GWW avec descente su seuil en tuant un nombre donné de particules à chaque fois */
/** GWW with a threshold descent such as a given number of particles is regrouped */
class AdaptiveGWWAlgorithm : public ThresholdGWWAlgorithm
{ public :
/* nombre de mauvaises particules à regrouper sur les bonnes */
/** number of bad particles to be regrouped on good ones */
  int nbkilled;
  AdaptiveGWWAlgorithm(int population_size, int grtest, int lastmove, int elitisme, int stop, int nbkilled );
  void regrouping(Configuration** population);
  void thresholdcomputedelta(Configuration** population);
};

/* GWW avec descente du seuil au plus bas des valeurs de seuil obtenues par AdaptiveGWWAlgorithm et FastStandardGWWAlgorithm
 (un nombre de particules et un taux) */
/** GWW with a threshold descent at the lowest value obtained by AdaptiveGWWAlgorithm et FastStandardGWWAlgorithm
using a number of particles to be redistributed and a rate */
class FastAdaptGWWAlgorithm: public AdaptiveGWWAlgorithm
{public :
/* taux de descente du seuil */
/** threshold descent rate */
  double thresholddescent;
  int nbmaxkilled;
  FastAdaptGWWAlgorithm(int population_size, int grtest, int lastmove, int elitisme, int stop, int nbkilled, int maxkilled, double thresholddescent);
  void thresholdcomputedelta(Configuration** population);
};

/* GWW avec descente du seuil en fonction du médian de la population*/
/** GWW with a descent depending on a distance between the worst and the median particle */
class MedianAdaptGWWAlgorithm: public AdaptiveGWWAlgorithm
{public :
/* taux de baisse du seuil : fraction de la distance entre la pire et la médiane (entre 0 et 1) */
/** descent rate : porcentage of the distance between the worst and the median particles (between 0 and 1) */
  double mediandescent;
  MedianAdaptGWWAlgorithm(int population_size, int grtest, int lastmove, int elitisme, int stop, double mediandescent);
  void thresholdcomputedelta(Configuration** population);
};

/* GWW avec descente du seuil en fonction du meilleur de la population*/
/** GWW with a descent depending on a distance between the worst and the best particle */
class BestAdaptGWWAlgorithm: public AdaptiveGWWAlgorithm
{public :
/* taux de baisse du seuil : fraction de la distance entre la pire et la meilleure (entre 0 et 1) */
/** descent rate : porcentage of the distance between the worst and the best particles (between 0 and 1) */

  double bestdescent;
  BestAdaptGWWAlgorithm(int population_size, int grtest, int lastmove, int elitisme, int stop, double bestdescent);
  void thresholdcomputedelta(Configuration** population);
};



/* GWW sans seuil : 2 paramètres : nombre de tués à chaque itération, nombre d'itérations */
/** GWW without threshold management : 2 parameters : number of particles redistributed at each iteration , number of iterations */
class NothresholdGWWAlgorithm : public GWWAlgorithm 
{public :
  NothresholdGWWAlgorithm(int population_size, int grtest, int lastmove, int elitisme, int stop,
	int killed,  int nbiter);
  void regrouping(Configuration** population);
/* nombre de particules à regrouper à chaque itération */
/** number of particles to be regrouped at each iteration */
  int nbkilled;
};




