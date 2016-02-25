


/* ------------------------------------------------------------
Les fonctions utilitaires de incop : Statistiques, 
lecture des arguments des algos , ecriture,  lancement d'un essai, 
------------------------------------------------------------ */



#include <assert.h>
#include <cerrno>
#include <stdio.h>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
using namespace std;
#include <iostream>
#include <fstream>
#include "timer.h"
#include "incop.h"
#include "incoputil.h"



#include <math.h>
#include <unistd.h>

#include <stdlib.h> 
//#include <signal.h>

ofstream* ofile = NULL;  // le fichier de sortie

Stat_GWW * Statistiques; //  l'objet pour les statistiques en variable globale 
                         // allou� dans le main() avec npb et nbessais





int TRACEMODE=0;  // variable globale : niveau de trace

//struct sigaction Action;  // trombe_ajout : pour les signaux


/* ------------------------------------ STATISTIQUES ---------------------------------------*/

Stat_GWW::Stat_GWW (int number_problems, int number_tries) {
  nb_pbs=number_problems;
  trouve=new int[number_problems]; 
  for (int i=0 ; i < number_problems ; i++) { trouve[i]=0; } 
  max_tries=number_tries;
  cost_try=new Long[number_tries];  
  nb_moves= new int[number_tries];
  nb_moves_up=new int[number_tries];
  nb_moves_down=new int[number_tries];
  nb_moves_avg=0;
  total_problem_time=new float[number_problems];
  execution_time_try=new float[number_tries]; 
  total_execution_time=0.0;
  current_pb=0;
}

void Stat_GWW::init_pb(int t)
{current_pb = t; 
 thresholdchanges=0;
}


void Stat_GWW::init_run()
{trouve[current_pb]=0;}

void Stat_GWW::init_try(int trynumber)
{current_try = trynumber;
 nb_moves[trynumber]=0;
 nb_moves_up[trynumber]=0;
 nb_moves_down[trynumber]=0;
 execution_time_try[trynumber]=0;
 costvalues.push_back(0);
 costvalues.clear();
 examinedneighbors.push_back(0);
 examinedneighbors.clear();
}


void Stat_GWW::execution_report(int nessai, Long lower_bound)
{
#ifdef WINDOWS
execution_time_try [nessai] += 1.00;
cout << " WARNING : Timer not supported under windows OS : ==> incop exec time is false " << endl;
#else
execution_time_try [nessai] += VIRTUAL_TIMELAPSE; 
#endif

    if (cost_try[nessai]== lower_bound) 
      { trouve[current_pb]++;}
    ecriture_stat_essai();
    
    total_execution_time += (float)execution_time_try[nessai];
    *ofile << " temps total execution " <<  total_execution_time << endl    ;	  
}

//void sigaction()
//  {Action.sa_handler=handler_stat;           // trombe_ajout : d�tournement du signal SIGUSR1
//  sigaction(SIGUSR1, &Action, NULL);        // trombe_ajout : d�tournement du signal SIGUSR1
//  }


// ajout_trombe : rajout de la fonction suivante pour les signaux !!!

//void handler_stat (int sig) {
//
//  *ofile  << "==========================================================================" << endl;
//  *ofile  << "Signal " << sig << " recu !!" << endl;
//
// // ecriture_stat_probleme();   // stats sur les derniers essais du probl�me courant
//
//  if (Statistiques->current_pb > 1) {
//    ecriture_statistiques_global (); // stats sur tous les probl�mes
//  }
//
//  cout << "Fin resolution (interrompue) en : " << Statistiques->total_execution_time
//       << " secondes" << endl;
//
//  sleep(1);
//
//  kill (getpid(), 9);   // suicide
//
//}



/*----------------------------------- LECTURE DES ARGUMENTS ----------------------------------*/

/* les fonctions de lecture et de v�rification de type et de valeur des donn�es */

int argument2ul(char* arg, char* message)
{char* error;
    int s = strtoul(arg,&error,10);
 if (error==arg  || errno==ERANGE)
   {cerr << "Erreur " << message << " pas un entier " <<  arg << endl;
   exit(1);
   }
 return s;
}


double argument2d(char* arg, char* message)
{char* error;
  errno=0; 
    double s = strtod (arg,&error);
 if (errno==ERANGE  || error == arg)
   {cerr << "Erreur " << message << " pas un nombre " <<  arg << endl;
   exit(1);
   }
 return s;
}

double argument2bd(char* arg, char* message, double min1, double max1)
{double s= argument2d(arg,message);
if (( s< min1) || (s > max1) )
  {cerr << "Erreur " << message << arg << " doit �tre compris entre " << min1 << " et " << max1;
  exit(1);
  }
return s;
} 

int argument2bul(char* arg, char* message, int min1, int max1)
{int s= argument2ul(arg,message);
if (( s< min1) || (s > max1) )
  {cerr << "Erreur " << message << arg << " doit �tre compris entre " << min1 << " et " << max1;
  exit(1);
  }
return s;
} 



string argument2lp(char* arg, char* message, list<string>& possibles)
{string s = arg;
 if (find(possibles.begin(),possibles.end(),arg)== possibles.end()) 
   {cerr << message << arg << endl; exit(1);}
return s;
}

/* la liste des m�thodes implant�es */

void definir_liste_methodes(list<string>& liste_methodes)
{liste_methodes.push_back("metropolis");
liste_methodes.push_back("tabu");
liste_methodes.push_back("incrtabu");

 liste_methodes.push_back("idw"); // idw 
 liste_methodes.push_back("idwa"); // idw avec min-voisins = 1 et epuisement = 1 (any)

 liste_methodes.push_back("idwb"); //  idw min-voisins =1 et epuisement=max-voisins (best)
 liste_methodes.push_back("idwbsn"); //  idw min-voisins =1 et epuisement=max-voisins dynamique
 liste_methodes.push_back("idwgra"); // idw min-voisins = max-voisins et epuisement=1 (any)

 liste_methodes.push_back("idwgrb"); // idw min-voisins=max-voisins et epuisement= max-voisins (best)

 liste_methodes.push_back("idwtabu"); // idw avec liste taboue
 liste_methodes.push_back("idwatabu"); // idwa avec liste taboue
 liste_methodes.push_back("idwbtabu"); // idwb avec liste taboue
 liste_methodes.push_back("idwgratabu"); // idwgra avec liste taboue
 liste_methodes.push_back("idwgrbtabu"); // idwgrb avec liste taboue
 liste_methodes.push_back("idwincrtabu"); // idw avec liste taboue
 liste_methodes.push_back("idwaincrtabu"); // idwa avec liste taboue
 liste_methodes.push_back("idwbincrtabu"); // idwb avec liste taboue
 liste_methodes.push_back("idwgraincrtabu"); // idwgra avec liste taboue
 liste_methodes.push_back("idwgrbincrtabu"); // idwgrb avec liste taboue

liste_methodes.push_back("simann");
liste_methodes.push_back("taburate");
liste_methodes.push_back("descente");
 liste_methodes.push_back("greedy");  // idw voisinage explor� (min=max) cf autotuning2.cc
liste_methodes.push_back("idwaminmax");
liste_methodes.push_back("idwbminmax");
liste_methodes.push_back("idwupk");
liste_methodes.push_back("idwgrupk");
 liste_methodes.push_back("grwrate");  
liste_methodes.push_back("threshold");
liste_methodes.push_back("random");
liste_methodes.push_back("gww");
liste_methodes.push_back("gww-nothreshold"); 
liste_methodes.push_back("gww-killworst");
liste_methodes.push_back("gww-adaptkillworst");
liste_methodes.push_back("gww-mediandist");
liste_methodes.push_back("gww-bestdist");
liste_methodes.push_back("gww-adapt");
}

void arguments_arret(char** argv, int& narg, int& stop)
{stop = argument2ul(argv[narg+1], " arret 1re sol");
// *ofile << " arret 1re sol " << stop << endl;
 narg++;
}

void arguments_tempscpu(char** argv, int& narg, double& maxtime)
{maxtime = argument2d( argv[narg+1], " max temps cpu" );
// *ofile << " temps cpu max " << maxtime << endl;
 narg++;}

void arguments_borneinf(char** argv, int& narg, Long& borneinf)
{borneinf= argument2ul( argv[narg+1], " borne inferieure " );
// *ofile << " borne inferieure " << borneinf << endl;
 narg++;
}

void arguments_methode(char** argv, int& narg, int& graine1, int& nbessais, string& method, list<string>& methodes_possibles)
{  graine1= argument2ul (argv[narg+1], " graine du g�n�rateur al�atoire de la m�thode ");
  nbessais = argument2ul (argv[narg+2], " nombre d'essais ");
  method = argument2lp (argv[narg+3], "methode non implantee ", methodes_possibles);

  narg=narg+3;
//  *ofile << " graine tirage al�atoire " << graine1
//	 << " nb essais " << nbessais << " m�thode " << method << endl;
}

void arguments_tracemode(char** argv, int& narg)
{ TRACEMODE = argument2ul (argv[narg+1], "indicateur trace ");
// *ofile << " mode trace " << TRACEMODE << endl;
 narg++;}


void arguments_algorithme(char** argv, int& narg, int& nbmouv)
{nbmouv= argument2ul (argv[narg+1], "longueur marche ");
// *ofile << " longueur marche " << nbmouv << endl;
 narg=narg+1;
}


void arguments_metropolis(char**argv, int& narg, double&  temp)
{ temp = argument2d (argv[narg+1], " temperature " );
// *ofile << " temperature " << temp;
 narg=narg+1;
}

void arguments_voisinage(char** argv,int& narg, int& taille_voisinage_min, int& taille_voisinage_max, int & fin_voisinage, int& var_conflit, int& val_conflit, int & dynamic )

{       
    string variables_en_conflit = argv[narg+1];

    if (variables_en_conflit == "cv") 
      var_conflit=1;

    else var_conflit=0;

    string minimum_conflit = argv[narg+2];      
    if (minimum_conflit== "mc")
      val_conflit=1;
    else val_conflit=0;
    taille_voisinage_min= argument2ul (argv[narg+3], " nombre minimum de voisins ");
    taille_voisinage_max= argument2ul (argv[narg+4], "nombre maximum de voisins ");
    //    fin_voisinage=argument2bul(argv[narg+5],"indicateur epuisement voisinage",0,taille_voisinage_max);
    fin_voisinage=argument2bul(argv[narg+5],"indicateur epuisement voisinage",0,RAND_MAX);
    dynamic = argument2bul(argv[narg+6], " voisinage dynamique", 0, 1);
//    *ofile << " variables en conflit " << variables_en_conflit << " " << var_conflit;
//    *ofile << " minimum conflit " << minimum_conflit << " " << val_conflit;
//    *ofile << " min voisins " << taille_voisinage_min ;
//    *ofile << " max voisins " << taille_voisinage_max;
//    *ofile << " epuisement voisinage " << fin_voisinage << endl;
//    *ofile << " voisinage dynamique " << dynamic << endl;
    narg= narg+6;
}

void arguments_tabu(char** argv,int & narg, int& longtabu)
{longtabu = argument2ul(argv[narg+1]," longueur liste taboue ");
// *ofile << " longueur liste taboue " << longtabu;
 narg= narg+1;
}

void arguments_recuit(char** argv,int & narg, double& inittemp)
{inittemp = argument2d(argv[narg+1]," temperature initiale ");
// *ofile << " temperature initiale " << inittemp;
 narg= narg+1;
}

void arguments_marcheseuil(char** argv,int & narg, int& seuildebut)
{seuildebut = argument2ul (argv[narg+1]," seuil debut  " );
// *ofile << " seuil debut  " << seuildebut;
 narg= narg+1;
}

void arguments_gww(char** argv, int & narg, int& taille, int& test_regroupement,
		   int & derniermouv, int & elitisme, int& stop)
{taille = argument2ul (argv[narg+1]," nombre particules " );
 test_regroupement = argument2bul (argv[narg+2]," test regroupement ",0,1 );
 derniermouv = argument2bul (argv[narg+3]," baisse dernier mouvement ",0,1 );
 elitisme= argument2bul (argv[narg+4], " indicateur elitisme ",0,1);
 stop = argument2bul (argv[narg+5], " indicateur arret stagnation ",0,1);
// *ofile << " nombre particules " << taille ;
// *ofile << " test regroupement " << test_regroupement;
// *ofile << " baisse dernier mouvement " << derniermouv;
// *ofile << " elitisme " << elitisme;
// *ofile << " arret stagnation " << stop;
narg=narg+5;}

void arguments_gww_marche(char** argv, int& narg, string& walk_method, list<string>& liste_methodes)
{
 walk_method = argument2lp (argv[narg+1], " methode marche non implant�e " , liste_methodes);
//*ofile << " methode marche " << walk_method << endl;
narg=narg+1;
}


void arguments_gww_standard(char** argv, int & narg, double& descenteseuil, int& seuilmin)
{descenteseuil= argument2bd(argv[narg+1]," facteur descente seuil ", 0, 1 );
seuilmin=argument2ul (argv[narg+2]," borne inferieure ") ;
// *ofile << " facteur descente seuil " << descenteseuil;
// *ofile << " borne-inf�rieure " << seuilmin;
 narg=narg+2;
}

void arguments_gww_descente_rapide(char** argv, int & narg, int & nb_tues, int & nb_tues_max, double& descenteseuil)
{nb_tues= argument2ul (argv[narg+1] ," nb tu�s " );
nb_tues_max=argument2ul (argv[narg+2]," nb tu�s max ");
 descenteseuil= argument2bd (argv[narg+3]," facteur descente seuil ", 0, 1 );
//*ofile << " nb tu�s " << nb_tues;
// *ofile << " nb tues max " << nb_tues_max;
// *ofile << " facteur descente seuil " << descenteseuil;
 narg=narg+3;
}


void arguments_gww_distance_median(char** argv, int & narg,  double& distance_median)
{ distance_median= argument2bd (argv[narg+1], " facteur distance median ",0, 1);
// *ofile << " facteur distance median " << distance_median;
 narg=narg+1;
}

void arguments_gww_distance_meilleur(char** argv, int & narg, double& distance_meilleur)
{
 distance_meilleur= argument2bd (argv[narg+1]," facteur distance meilleur " ,0,1);
// *ofile << " facteur distance meilleur " << distance_meilleur;
 narg=narg+1;
}



void arguments_gww_adaptatif(char** argv, int& narg, int & nb_tues)
{nb_tues = argument2ul (argv[narg+1]," nb tu�s "  );
//*ofile << " nb tu�s " << nb_tues;
narg++;
}

void arguments_gww_sans_seuil (char** argv, int & narg, int & nb_tues, int& nb_iter)
{nb_tues =  argument2ul (argv[narg+1]," nb tu�s " );
 nb_iter=   argument2ul (argv[narg+2]," nb iterations "  );
//*ofile << " nb tu�s " << nb_tues;
//*ofile << " nb iterations " << nb_iter;
 narg= narg+2;
}

void arguments_taburate(char** argv, int & narg, float & Pd, float & P0, int & longtabu) {
  Pd = argument2bd (argv[narg+1], " Pd ",0,1);
  P0 = argument2bd (argv[narg+2], " P0 ",0,1);
  longtabu = argument2ul(argv[narg+3], " longueur liste taboue ");
//  *ofile << " Pd = " << Pd << endl;
//  *ofile << " P0 =  " << P0 << endl;
//  *ofile << " longueur liste taboue = " << longtabu << endl;
  narg= narg+3;
}

void arguments_grwrate(char** argv, int & narg, double & nbr) {
  nbr = argument2bd(argv[narg+1], "taux voisinage " , 0, 1);
//  *ofile << " taux de voisinage " << nbr << endl;
  narg++;
}

/*-----------------------------------------CREATION DES OBJETS ALGORITHMES----------------------------*/
// lecture des parametres de la marche et creation de l'objet LSAlgorithm 

LSAlgorithm* algo_marche (char** argv,int& narg,string& method, int gww)
		
{LSAlgorithm* algo;
 Metaheuristic* mheuristic = new Metaheuristic();
 NeighborhoodSearch* nbhsear;
  int taille_voisinage_min, taille_voisinage_max, fin_voisinage,var_conflit,val_conflit,seuildebut;
  double temp=0;
  int longtabu=0;
  int dynamic=0;
  int nbmouv;
  double inittemp;
  double nbhr=0;  // defini pour grwrate
  float Pd, P0; 
  arguments_algorithme(argv,narg,nbmouv);
//  *ofile << " methode " << method << endl;
  if (method == "metropolis")
    arguments_metropolis(argv,narg,temp);
  else if(method == "tabu" || method == "incrtabu" || method == "idwtabu" || method =="idwatabu" || method =="idwbtabu" || method =="idwgratabu" || method =="idwgrbtabu" || method =="idwincrtabu" ||method =="idwaincrtabu" || method =="idwbincrtabu" || method =="idwgraincrtabu" || method =="idwgrbincrtabu" )
    arguments_tabu(argv,narg,longtabu);
  else if(method == "simann")
    arguments_recuit(argv,narg,inittemp);
  else if (method == "threshold")
    arguments_marcheseuil(argv,narg,seuildebut);
  else if (method == "taburate" )
    arguments_taburate(argv, narg,  Pd, P0, longtabu);
  else if (method == "grwrate")
    arguments_grwrate (argv,narg,nbhr);
//  *ofile << " arguments voisinage " << endl;
  arguments_voisinage(argv,narg,taille_voisinage_min,taille_voisinage_max,fin_voisinage,var_conflit,val_conflit,dynamic);
  if (gww)
    algo = new LSAlgorithmGWW (nbmouv);
  else 
    algo = new LSAlgorithm (nbmouv);
  if (dynamic)
    {nbhsear = new DynamicNeighborhoodSearch (taille_voisinage_min,taille_voisinage_max,fin_voisinage,var_conflit,val_conflit, nbhr);}
  else if (method == "idwbsn")
    nbhsear = new DynamicSpareneighbor (taille_voisinage_min,taille_voisinage_max,fin_voisinage,var_conflit,val_conflit, nbhr);
  else
    nbhsear = new NeighborhoodSearch (taille_voisinage_min,taille_voisinage_max,fin_voisinage,var_conflit,val_conflit, nbhr);

  algo->nbhsearch=nbhsear;
 if (method == "metropolis")
   mheuristic = new Metropolis(temp);
 else if (method == "tabu")
   mheuristic = new TabuSearch (longtabu);
 else if (method == "incrtabu")
   mheuristic = new IncrTabuSearch (longtabu);
 else if (method == "idwtabu" || method =="idwatabu" || 
	  method == "idwbtabu" ||
	  method =="idwgratabu" || method =="idwgrbtabu")
   mheuristic = new TabuGreedySearch (longtabu);
 else if (method == "idwincrtabu" || method =="idwaincrtabu" || 
	  method == "idwbincrtabu" ||
	  method =="idwgraincrtabu" || method =="idwgrbincrtabu")
   mheuristic = new IncrTabuGreedySearch (longtabu);
 else if (method == "simann")
   mheuristic = new SimulatedAnnealing (inittemp,nbmouv);
 else if (method == "threshold")
   mheuristic = new ThresholdAccepting(seuildebut,nbmouv);
 else if ((method == "descente") || (method == "greedy") || (method == "grwrate") 
	  || (method == "idwupk") || (method == "idwgrupk") || (method == "idwa") ||
	  (method == "idwb") || (method == "idwbsn") || (method == "idwgra") || (method == "idwgrb") ||
         (method == "idw" ) || (method == "idwaminmax") || (method == "idwbminmax"))
   mheuristic = new GreedySearch ();
 else if (method == "random")
   mheuristic= new RandomSearch ();
 else if (method == "taburate") 
   mheuristic = new TabuAcceptingrate (longtabu, Pd, P0);
 algo->mheur = mheuristic;
 return algo; 
}

// lecture des parametres GWW et creation de l'objet GWWAlgorithm

GWWAlgorithm* algo_gww(char** argv, int& narg, string& method , int& taille, list<string>& liste_methodes)
{GWWAlgorithm* algogww = new GWWAlgorithm();
 int testregroupement, nb_tues, nb_tues_max, nb_iter, derniermouv,elitisme, stop;
 double descenteseuil, distance_median, distance_meilleur;
 int seuilmin;
 string walk_method;
 arguments_gww(argv,narg,  taille, testregroupement, derniermouv,elitisme, stop);
 if (method == "gww-nothreshold")
   {arguments_gww_sans_seuil (argv, narg,  nb_tues, nb_iter);
   algogww = new NothresholdGWWAlgorithm 
     (taille,testregroupement, derniermouv,elitisme, stop, nb_tues,nb_iter);}
 else if (method == "gww")
   {arguments_gww_standard (argv, narg, descenteseuil, seuilmin);
   algogww = new StandardGWWAlgorithm 
     (taille,testregroupement, derniermouv,elitisme, stop, descenteseuil, seuilmin);
   }
 else if (method == "gww-killworst")
   {arguments_gww_standard (argv, narg, descenteseuil, seuilmin);
   algogww = new FastStandardGWWAlgorithm 
     (taille,testregroupement, derniermouv,elitisme, stop, descenteseuil, seuilmin);
   }

 else if (method == "gww-adapt")
   {arguments_gww_adaptatif (argv, narg, nb_tues);
   algogww = new AdaptiveGWWAlgorithm (taille, testregroupement, derniermouv,elitisme,stop, nb_tues);
   }
 else if (method == "gww-adaptkillworst")
   {arguments_gww_descente_rapide (argv, narg, nb_tues, nb_tues_max, descenteseuil);
   algogww = new FastAdaptGWWAlgorithm (taille, testregroupement, derniermouv, elitisme, stop, nb_tues, nb_tues_max, descenteseuil);}
 else if (method == "gww-mediandist")
   {arguments_gww_distance_median (argv, narg, distance_median);
   algogww = new MedianAdaptGWWAlgorithm (taille, testregroupement, derniermouv, elitisme, stop, distance_median);}
 else if (method == "gww-bestdist")
   {arguments_gww_distance_meilleur (argv, narg, distance_meilleur);
   algogww = new BestAdaptGWWAlgorithm (taille, testregroupement, derniermouv, elitisme,  stop, distance_meilleur);
   }


 arguments_gww_marche(argv,narg,walk_method,liste_methodes);
 algogww->walkalgorithm= algo_marche (argv,narg,walk_method,1);
 algogww->walkalgorithm->methodname=method; // bizarrerie a modifier (utilisee pour la trace)
 return algogww;
}


// lecture de l'argument method et appel selon l'argument de la cr�ation d'un des  2 principaux types d'algo (LS et GWW)
IncompleteAlgorithm* algo_creation(char** argv, int& narg, int& taille, int& graine1, int & nbessais)
 { IncompleteAlgorithm* algo = new IncompleteAlgorithm();
  list<string> liste_methodes;
  string method;
  definir_liste_methodes(liste_methodes);
  arguments_methode(argv,narg,graine1,nbessais,method,liste_methodes);
  if (method == "gww" || method =="gww-killworst"|| method == "gww-nothreshold" || method == "gww-adapt" || method == "gww-adaptkillworst"
      || method == "gww-mediandist" || method == "gww-bestdist")
    algo= algo_gww (argv,narg,method,taille,liste_methodes);	
  else {algo = algo_marche (argv,narg,method,0);
  taille=1;   
  }
  algo->methodname=method;
  return algo;
 }


//********************************** les utilitaires***********************************************

/* la plus mauvaise valeur de la population */
Long valeur_max(Configuration** population, int taille)
{return population[0]->valuation;
}

/* la meilleure valeur de la population */
Long valeur_min(Configuration** population, int taille)
{ return population[taille-1]->valuation; 
}

/* la valeur mediane de la population */
Long valeur_mediane(Configuration** population, int taille)
{ return population[taille/2]->valuation; 
}

// le comparateur pour le tri de la population 
static int comparepopulation(const void* e1,const void* e2)
{  return ((*(Configuration**) e1)->valuation < (*(Configuration**)e2)->valuation);}

// le tri de la population dans le sens du pire au meilleur
void populationsort(Configuration** population, int taille)
{  qsort(population,taille,sizeof(Configuration*),comparepopulation);}




/* ------------------------------ECRITURES--------------------------------------------------*/


/* le nom du fichier de sortie : pour les tests : resultatsxxx/concat�nation des arguments */
void ofile_name (char* filename, int argc, char** argv)
{ 
  sprintf(filename,"%s%s","results/",argv[2]);
  for (int i=3;i<argc;i++)
    sprintf(filename,"%s-%s",filename,argv[i]);
  cout << filename << endl;
}


/* les ecritures sur le fichier de sortie *ofile */



void ecriture_graine(int graine, int nessai)
{*ofile << " essai n. " << nessai << " graine " << graine << endl;}

void ecriture_fin_resolution(Long meilleur)
{*ofile<< "meilleure valeur --------->  " <<  meilleur << " " ;
 *ofile << " nb mouvements " << Statistiques->nb_moves[Statistiques->current_try] <<  "  " ;
 *ofile << " nb descentes " << Statistiques->nb_moves_down[Statistiques->current_try] <<  "  " ;
 *ofile << " nb remontees " << Statistiques->nb_moves_up[Statistiques->current_try] <<  "  " << endl ;
}


void ecriture_debut_resolution(Long pire, Long meilleur,string& method)
{*ofile<< "configuration initiale : " << " pire valeur  " << pire 
       << " meilleure valeur " << meilleur << endl ;
}


void ecriture_stat_essai()
{
  *ofile<<  " Temps essai: " << Statistiques->execution_time_try[Statistiques->current_try]
                             << endl;
}


void ecriture_stat_probleme() {	
  // role : �criture des stats pour le probleme courant avec les essais de 0 � current_try-1
  int current_problem;
  int nbessais;
  float temps_essais_moyen;
  float nombre_chgseuil_moyen;

  current_problem = Statistiques->current_pb;
  nbessais = Statistiques->current_try;

  *ofile << "=========================================================================="<< endl;

  if (nbessais == 0) {
    *ofile << " Interruption !! " << endl;
  }
  else {
    
    *ofile << " " << Statistiques->trouve[current_problem] << " essais avec solutions sur " 
	   <<  nbessais  << " essais. " << endl;

    double somme=0; double totalnbmoves=0; Long meilleur=LONG_MAX; Long pire=0; double ecart=0; double moyenne;

    for(int i=0; i < nbessais; i++) { 
      somme += Statistiques->cost_try[i];
	  //      somme += Statistiques->cost_try[i].to_double();
      totalnbmoves += Statistiques->nb_moves[i];
      if (Statistiques->cost_try[i] < meilleur) meilleur = Statistiques->cost_try[i];
      if (Statistiques->cost_try[i] > pire) pire = Statistiques->cost_try[i];
    }
    Statistiques-> nb_moves_avg = totalnbmoves /nbessais;
    moyenne = somme / nbessais ;
    Statistiques->cost_meanvalue = moyenne;

    for(int i=0; i < nbessais ;i++) { 
      ecart += pow((Statistiques->cost_try [i] - moyenne),2);
	  //      ecart += pow((Statistiques->cost_try [i] - moyenne).to_double(),2);
    }
    ecart=sqrt(ecart/nbessais);

    temps_essais_moyen = 0.0;
    for(int i=0; i < nbessais; i++) { 
      temps_essais_moyen += Statistiques->execution_time_try[i];
    }
    temps_essais_moyen = temps_essais_moyen / nbessais ;
    
    Statistiques->average_execution_time=temps_essais_moyen;
    nombre_chgseuil_moyen = Statistiques->thresholdchanges/nbessais;    
    *ofile   << " Meilleur essai : "    << meilleur      << endl;
    *ofile   << " Plus mauvais : "      << pire          << endl ;
    * ofile  << " Moyenne : "           << moyenne       << endl;
    * ofile  << " Ecart type : "        << ecart         << endl;
    *ofile   << " Temps essai moyen: "  << temps_essais_moyen << endl;
    if (nombre_chgseuil_moyen)  // algo de type GWWThresholdAlgorithm
      *ofile   << " Nb changements seuil moyen : " << nombre_chgseuil_moyen << endl;
    *ofile    << " Nb mouvements moyen : " << Statistiques->nb_moves_avg  << endl;
    *ofile   << "=========================================================================="
	     << endl << endl;
  }
}



void ecriture_changement_seuil (Long seuil, Long delta, Long meilleur, Long pire,Long mediane, int nbessaisvoisins, int nb_au_seuil)
{
  *ofile << " --- seuil " << seuil << " delta " << delta
	 << " meilleur " << meilleur
	 << " pire " << pire
	 << " mediane " << mediane 
	 << " au seuil " << nb_au_seuil
    	 << " voisins test�s " << nbessaisvoisins
         << endl;
}

void ecriture_nb_tues(int nb_tues)
{*ofile << " nombre particules redistribu�es " << nb_tues << endl;}


void ecriture_fin_gww( int nb_chang_seuil, int nb_mouv)
{  *ofile << " Fin GWW  : nombre de changements de seuil : " << nb_chang_seuil << endl;
   *ofile << " nombre total de mouvements : " << nb_mouv << endl;}



void ecriture_meilleure_valeur (string & method, Long valeur, Long seuil, int nbmouv, int maxvoisins)
{
*ofile << " meilleure valeur " << valeur ;
  if (method == "gww" || method == "gww-killworst" || method == "gww-adapt" || method == "gww-adaptkillworst" || method == "gww-mediandist" 
|| method == "gww-bestdist")    *ofile << " valeur seuil " << seuil ;
  else  if (method == "gww-nothreshold") 
    *ofile   << " iteration "  << nbmouv;
  else *ofile <<" nb mouvements " << nbmouv;
  *ofile       <<  " temps " <<  Statistiques->execution_time_try [Statistiques->current_try] <<
       " voisinage " << maxvoisins 
               << endl;   
}


void ecriture_statistiques_global () {

  int pbsol;      // le nombre de problemes avec au moins un essai r�ussi
  int npb;        // le nombre de problemes etudi�s
  int nb_succes;  // le nombre d'essais r�ussis sur tous les problemes

  npb = Statistiques->current_pb + 1;

  pbsol=0;
  for (int t=0 ; t < npb ; t++) {
    if (Statistiques->trouve[t] >=1) 
      pbsol++;
  }
  *ofile <<"==========================================================================" << endl;
  *ofile << pbsol << " problemes avec solution(s) sur " << npb << " problemes " << endl;

  *ofile << " Succes par probleme : | " ;
  for (int i=0 ; i < npb ; i++) {
    *ofile<< Statistiques->trouve[i] << " | " ; 
  }
  *ofile << endl;

  nb_succes=0;
  for (int i=0 ; i < npb ; i++) {
    nb_succes += Statistiques->trouve[i];
  }
  *ofile<< " Nb succes " << nb_succes << " sur " 
	<< (npb-1) * Statistiques->max_tries +  Statistiques->current_try // nb d'essais achev�s
	<< " essais " << endl;

  if (nb_succes > 0) {
    *ofile << " temps moyen pour un succes " 
	   << Statistiques->total_execution_time / nb_succes << endl;
  }

  *ofile << " Temps total : " << Statistiques->total_execution_time << endl; 
  *ofile <<"==========================================================================" <<endl; 
}


void ecriture_fin_lsrun(double avgnhtries, double avgsqnhtries)
{// *ofile << " nb moyen d'essais par mouvement " << avgnhtries << "  ecart type " << sqrt ( avgsqnhtries - avgnhtries * avgnhtries) << endl;
 if (TRACEMODE==2)
   {*ofile << " valeurs   et nb-voisins ";
   for (int i= Statistiques->costvalues.size()- 100 ; i< (int) Statistiques->costvalues.size(); i++)
     {*ofile << Statistiques->costvalues[i] << "  " ;
     *ofile << Statistiques->examinedneighbors[i] << "  " << endl;}
   }
}


/* --------------------------LANCEMENT D'UN ESSAI -------------------------------------------*/


// instanciation aleatoire
void instanciation_aleatoire(OpProblem* problem, Configuration** population, int taille)
{for (int i=0;i<taille;i++)
  problem->random_configuration(population[i]);
}


// evaluation de la population
void calcul_valeur_population(OpProblem* problem, Configuration** population,int taille)
{for (int i=0;i<taille;i++)
  population[i]->valuation = problem->config_evaluation(population[i]);
}

// execution d'un essai d'un algo sur un  probleme
void executer_essai
   (OpProblem* problem,IncompleteAlgorithm* algo, Configuration** population, int taille, int graine1, int nessai, vector<int> *initconfig)
  {
    // graine du g�n�rateur aleatoire pour l'essai
    srand48(graine1+nessai);
    srand (graine1+nessai);
//    ecriture_graine(graine1+nessai,nessai);

    Statistiques->init_try(nessai);
    // d�clenchement du chronom�tre
#ifndef WINDOWS
    start_timers(); 
#endif
    // population initiale 

    instanciation_aleatoire(problem,population,taille);
    //SdG: initial solution provided by NC/EAC supports given to INCOP
    if (initconfig && nessai==0) {
        assert(initconfig->size() == population[0]->nbvar);
        for (int i=0; i<population[0]->nbvar; i++) population[0]->config[i] = (*initconfig)[i];
    }
    //    *ofile << " population instanciee " << endl;
    // evaluation de la population
    calcul_valeur_population(problem,population,taille);
    //    *ofile << " population evaluee " << endl;
    // tri de la population
    populationsort(population,taille);
    // stockage du meilleur dans best_config
    problem->best_config->copy_element(population[taille-1]); 
    // pour les algos de type gww avec seuil , le seuil est initialis� au pire de la population, sinon � RAND_MAX
    algo->initthreshold(population,taille); 
    
    Statistiques->cost_try[nessai]=valeur_min (population,taille);  
//    ecriture_debut_resolution(valeur_max(population,taille),valeur_min(population,taille),algo->methodname);
    
    // lancement de la resolution 
    algo->run(problem,population);
    // apres resolution : arret du chronometre
#ifndef WINDOWS // time not supported under windows OS
    stop_timers(VIRTUAL);
#endif 
//    ecriture_fin_resolution(Statistiques->cost_try[nessai]);
//    problem->best_config_analysis();
//    problem->best_config_write();
    // verification de best_config en recalculant sa valeur 
//    problem->best_config_verification();
//    Statistiques->execution_report(nessai,problem->lower_bound);
}






