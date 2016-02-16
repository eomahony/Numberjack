
 /* méthodes des classes  OpProblem, CSProblem, BinaryCSProblem  */

#include <list>
#include <vector>
#include <string>
#include <set>
#include <algorithm>


using namespace std;
#include <fstream>
#include <stdio.h>
#include "incop.h"
#include "csproblem.h"


extern ofstream* ofile;  // le fichier de sortie



void CSProblem::init_domains(int nbvar,int s)
{      for (int i=0;i<nbvar;i++)
  	domains[i]=0;
}

/* un seul domaine de taille s (valeurs de 0 à s-1) */
void CSProblem::init_tabdomains(int s)
{      tabdomains[0].clear();
      for (int i=0;i<s;i++)
	tabdomains[0].push_back(i);
}


/* calcul des variables en conflit : on reconstruit le vecteur des variables en conflit*/
void CSProblem::compute_var_conflict(Configuration* configuration)
{configuration->var_conflict.clear();
 for (int i=0;i<nbvar;i++)
  {if (configuration->get_conflicts_problem(this,i,configuration->config[i])!=0)
    configuration->var_conflict.push_back(i);
  }
 
}

/* constructeur : par defaut, la borne inférieure (condition d'arret) est 0 */

int CSProblem::variable_domainsize (int var)
{return tabdomains[domains[var]].size();}

CSProblem::CSProblem (int nvar, int nconst)
{nbvar=nvar;nbconst=nconst; lower_bound=0;} 

CSProblem::CSProblem (int nvar, int nconst,int lower)
{nbvar=nvar;nbconst=nconst; lower_bound= lower;} 

BinaryCSProblem::BinaryCSProblem(int nbvar, int nbconst) : 
  CSProblem(nbvar,nbconst) {;}
BinaryCSProblem::BinaryCSProblem (int nbvar, int nbconst, int lower) :
  CSProblem(nbvar,nbconst,lower) {;}

CSProblem::~CSProblem(){ delete currentmove; delete firstmove ; delete bestmove;}



void CSProblem::init_population(Configuration** population,int populationsize)
{
for (int i=0;i< populationsize;i++)
      population[i]= create_configuration();
 best_config=create_configuration ();
}

/* par defaut , pas d'incrementalité */
Configuration* CSProblem::create_configuration()
    { return (new Configuration(nbvar));}



void CSProblem::best_config_analysis()
{;}

void CSProblem::random_configuration(Configuration* configuration)
{
      for(int j=0;j<nbvar;j++)
	{int indice = (int) (drand48() * variable_domainsize(j));
      	configuration->config[j]= index2value(indice,j);
	//        *ofile << " j " << j << " " << configuration->config[j] << endl;
	}
}

/* choix aleatoire d'une variable */
int CSProblem::random_variable(Configuration* configuration)
{//return  (int) (drand48() * nbvar);
return (rand() % nbvar );
}


/* cas des variables en conflit : reduction de la taille du voisinage au max(voisinage, nb var conflits * taille_domaine) */
void CSProblem::adjust_parameters (Configuration* configuration, int& maxneighbors, int& minneighbors)
    {
    if 	(maxneighbors > (int) (configuration->var_conflict.size()) * (domainsize -1))
      maxneighbors = configuration->var_conflict.size() * (domainsize -1);
    if 	(minneighbors > (int) (configuration->var_conflict.size()) * (domainsize -1))
      minneighbors = configuration->var_conflict.size() * (domainsize-1);
    }



/* on choisit une variable en conflit */
int CSProblem::random_conflict_variable(Configuration* configuration)
{// *ofile << "VC " << configuration->var_conflict.size() << endl;
  // return configuration->var_conflict[(int) (drand48() * configuration->var_conflict.size())];
  return configuration->var_conflict[rand() % configuration->var_conflict.size()];
}

/* choix aléatoire d'une valeur autre que la courante : renvoie l'indice de la valeur*/
int CSProblem::random_value(int var, int val)

{ if (variable_domainsize(var) ==1) return 0; // domaine a une valeur
// int  val_changee = (int) (drand48() * (tabdomains[domains[var]].size()-1));
 int  val_changee = rand() % (tabdomains[domains[var]].size()-1);
  if (index2value(val_changee,var) >= val)
      val_changee++;
  return val_changee;}

/* une valeur du domaine de var minimisant les conflits
  autre que la courante (val)  , tirage aléatoire
entre valeurs équivalentes : renvoie l'indice de la valeur */
int CSProblem::min_conflict_value(int var, int val, Configuration* configuration)
{  if (variable_domainsize(var)==1) return 0; // domaine a une valeur
   Long minpromises=LONG_MAX; Long promises=0; int k1=1; int j=0;
   for (int i =0; i< (int) (tabdomains[domains[var]].size()); i++)
     {if (index2value(i,var)==val) i++;
     if (i ==  variable_domainsize(var)) break;
     promises= configuration->get_conflicts_problem(this,var,index2value(i,var));
     if (promises < minpromises)
       {minpromises = promises; j=i; k1=1;}
     else if (promises == minpromises)
       {k1++;
       if (drand48() < 1/(double)k1) j=i;
       }
     }
   return j;
}

/* détermination du prochain mouvement a tester et évaluation de ce mouvement */
void CSProblem::next_move(Configuration* configuration, Move* move, NeighborhoodSearch* nbhs)
{
 if (nbhs->var_conflict)
    ((CSPMove*)move)->variable = random_conflict_variable(configuration);
 else
   ((CSPMove*)move)->variable = random_variable(configuration);
// tirage d'un voisin : la valeur n'est pour le moment que l'indice dans le domaine
 if (nbhs->val_conflict)
   ((CSPMove*)move)->value = min_conflict_value(((CSPMove*)move)->variable,
						configuration->config[((CSPMove*)move)->variable], configuration);
 else
   ((CSPMove*)move)->value = random_value(((CSPMove*)move)->variable, 
					  configuration->config[((CSPMove*)move)->variable]);
 // passage de l'indice à la valeur
 ((CSPMove*)move)->value=  index2value(((CSPMove*)move)->value,((CSPMove*)move)->variable);

 move->valuation = move_evaluation (configuration,move);
 // *ofile << ((CSPMove*)move)->variable << " " << ((CSPMove*)move)->value << " " << move->valuation << endl;


}

/* exécution d'un mouvement  : fonction de base - affectation du coût du mouvement*/
void OpProblem::move_execution(Configuration* configuration,Move* move)
{configuration->valuation = move->valuation;}

/* execution d'un mouvement pour un CSP */
void CSProblem::move_execution(Configuration* configuration,Move* move)
{ 
  OpProblem::move_execution(configuration,move);
  configuration->config[((CSPMove*) move)->variable]=((CSPMove*) move)->value;
}

/* calcul incremental d'un changement de valeur  : valable pour tout type de configuration (noincr, incr, fullincr)*/
Long CSProblem::move_evaluation 
                    (Configuration* configuration,Move* move)
{ int var_changee = ((CSPMove*)move)-> variable;
  int val_changee = ((CSPMove*)move)-> value;
  //  *ofile << " move evaluation " << var_changee << " " << val_changee << " " 
  //	 << configuration->valuation << " " << 
  //    configuration->get_conflicts_problem(this,var_changee,val_changee)  
  //	 << " " <<
  //    configuration->get_conflicts_problem(this,var_changee,configuration->config[var_changee]) 
  //	 << endl;
  return (configuration->valuation+configuration->get_conflicts_problem(this,var_changee,val_changee)
	  -configuration->get_conflicts_problem(this,var_changee,configuration->config[var_changee]));
}


/* optimisation fullincr  :  
   dans le cas full incr,  tout a été préparé dans le tableau tabconflicts qu'il suffit de consulter 
cette optimisation semble inutile */
/*
int CSProblem::move_evaluation 
                    (Configuration* configuration,Move* move)
{ int var_changee = ((CSPMove*)move)-> variable;
  int val_changee = ((CSPMove*)move)-> value;
  //  *ofile << " move evaluation " << var_changee << " " << val_changee << " " 

  return (configuration->valuation+
	   ((FullincrCSPConfiguration*)configuration)->tabconflicts[var_changee][val_changee]
	  -((FullincrCSPConfiguration*)configuration)->tabconflicts[var_changee][configuration->config[var_changee]]);
}
*/

/* Ecriture de la meilleure solution trouvée (best_config) */
void CSProblem::best_config_write()
{
 *ofile<< " meilleure solution " << endl;
 for (int i = 0; i< nbvar ; i++)
  *ofile << " variable " << i << " : " << best_config->config[i] << endl;
}


/* la création des 3 mouvements utilisés par un algo de recherche locale */
void OpProblem::allocate_moves()
{firstmove = create_move();
 bestmove=  create_move();
 currentmove =  create_move();
}

/* vérification du cout de la meilleure solution en le recalculant */
void OpProblem::best_config_verification()
{
    Long value1 = best_config->valuation;
    Long value2 = config_evaluation(best_config);
    *ofile << " verification " << value1 << " " << value2 << endl;
}

/* creation d'un mouvement CSPMove */
Move* CSProblem::create_move()
{CSPMove* move = new CSPMove();
 return (Move*)move;
}

// creation du tableau des contraintes  : constraint1[i][j] = numero de la contrainte entre i et j (i<j)
int** csp_constraintdatastructure(int nbvar)
{int** constraint1= new int*[nbvar];
  for(int i=0;i<nbvar;i++) 
    {constraint1[i]=new int[nbvar];
    for(int j=0;j<nbvar;j++)
      constraint1[i][j]=0;}
  return constraint1;
}

/* mise en place des structures du CSP */
void CSProblem::set_domains_connections( int* dom, vector<int>* tabledom, vector<int> * connect )
{ 
 domains = dom;
 tabdomains = tabledom;
 connections = connect;
}
 
void CSProblem::init_domain_tabdomain()
{init_domains(nbvar,domainsize);
 init_tabdomains(domainsize);
}

/* taille du tableau de valeurs pour le tabou incremental : nb de mouvements maximum */
int CSProblem::nbtabuindex() {return (nbvar * domainsize);}

int CSProblem::tabuindex(Move* move, Configuration* config)
{return (((CSPMove*)move)->variable * domainsize + ((CSPMove*)move)->value);}


int CSProblem::tabuinverseindex(Move* move, Configuration* configuration)
{return (((CSPMove*)move)->variable * domainsize + configuration->config[((CSPMove*)move)->variable]);}
