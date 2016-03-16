/* Réglage automatique d'un algo de recherche locale à 1 ou 2 paramètres  */

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
#include "autotuning2.h"

#include <math.h>
#include <unistd.h>

extern ofstream* ofile;  // le fichier de sortie
extern Stat_GWW * Statistiques; 

Tuning::Tuning (int pinit, int seed1, int nbessai, int maxtun, int bmin, int bmax)
 {parammin=pinit; parammax=2*pinit;
 seed= seed1; nbtries = nbessai; maxtuning = maxtun; stop=0; boundmin=bmin; boundmax=bmax; referencetime=0;}

DoubleTuning::DoubleTuning (int pinit, int pinit2, int seed1, int nbessai, int maxtun, int maxtun2, int bmin, int b2min, int bmax, int b2max) : Tuning(pinit,seed1,nbessai,maxtun,bmin,bmax) 
  { param2min=pinit2; param2max= 2*pinit2; parambest=pinit; bound2min=b2min; bound2max=b2max; maxtuning2= maxtun2;}


// arret de la mécanique : temps cpu épuisé ou tout trouvé ou arret 1er trouvé
int endsolve(double maxtime, int nbtries, int pb)
{return (Statistiques->total_execution_time > maxtime || Statistiques->trouve[pb]== nbtries ||
         (Statistiques->trouve[pb] && Statistiques->stop_trouve)); }





// algo de réglage d'un paramètre positif par dichotomie  (algorithme à un parametre)
void Tuning::run (LSAlgorithm* algo, OpProblem* problem, Configuration** population)
{ int minfound =0;
  trynumber = 0;
  double epsilon = 0.01;
  //  if (algo->methodname=="greedy" || algo->methodname =="descente") epsilon=0.1;
  onerunparam(problem,algo,population,parammin);
  if (end()) return;
  trynumber=1;
  valuemin=value;
  onerunparam(problem,algo,population,parammax);
  if (end()) return;
  valuemax=value;
  trynumber = 2;
  if ((valuemax - problem->lower_bound) > (valuemin - problem->lower_bound) * (1 - epsilon))  // epsilon sur la moyenne pour eviter d'agrandir le voisinage quand fin jamais atteinte quand le parametre a regler est une longueur de voisinage
    {parambest = parammin; valuebest= valuemin;
    if (parammin==boundmin) minfound=1;
    parammin=parammin/2;
    while (!minfound && trynumber < maxtuning && !end())
      {trynumber++;
      onerunparam(problem,algo,population,parammin);
      if (value > valuebest)
	{minfound=1;}
      else if (parammin==boundmin)
	{minfound=1; parambest=parammin; valuebest=value;}
      else
	{valuebest=value ; parammax= parambest;
	parambest=parammin; parammin = parammin/2; }
     }
    }
  else
    {parambest= parammax; valuebest = valuemax;
    parammax=parammax*2;
    while (!minfound && trynumber < maxtuning && !end() && parammax <= boundmax)
      {trynumber++;
      onerunparam(problem,algo,population,parammax);
      if ((value - problem->lower_bound)  > (valuebest - problem->lower_bound) * ( 1 - epsilon))
	{minfound=1; }
      else
	{ valuebest= value ; parammin= parambest;
	parambest=parammax; parammax = parammax*2; 
	if (parammax > boundmax) parammax = boundmax;}

      }
    }
   while (trynumber < maxtuning && parambest - parammin > 1 && parammax-parambest > 1 &&  !end() && parambest !=0)
     {trynumber++;
     *ofile << " parambest " << parambest << " parammin " << parammin << " parammax " << parammax << endl;
     paramdico1= (parambest+parammin)/2;
     onerunparam(problem,algo,population,paramdico1);
     if (end()) break;
     valuedico1=value;
     trynumber++;
     paramdico2= (parambest+parammax)/2;
     onerunparam(problem,algo,population,paramdico2);
     valuedico2=value;
     if (valuedico1 < valuebest && valuedico1 < valuedico2)
       {valuemax=valuebest; parammax=parambest;
       valuebest=valuedico1; parambest=paramdico1;}
     else if
       (valuedico2 < valuebest)
       { valuemin=valuebest; parammin=parambest;
       valuebest=valuedico2; parambest=paramdico2;}
     else {parammin=paramdico1; parammax=paramdico2;}
     }
   *ofile << " best param " << parambest << " best value " << valuebest << endl;
   *ofile << " Fin reglage " << Statistiques->total_execution_time << endl;
}

// reglage algo a 2 parametres : 2 boucles imbriquées 
// pour chaque valeur du parametre parametre param2, appel  du reglage du parametre param1 (par simplerun)
void DoubleTuning::run (LSAlgorithm* algo, OpProblem* problem, Configuration** population)
{ int minfound =0;
  trynumber1 = 0;
  int bestparam1= parambest;
  double epsilon = 0.01;
  *ofile << " reglage boucle externe " << endl;
  //  if (algo->methodname=="tabu" || algo->methodname == "tabugreedy") epsilon=0.1;
  simplerun (algo,problem,population,param2min);
  trynumber1=1;
  if (end()) return;
  int param1min=parambest;
  value2min=valuebest;
  simplerun (algo,problem,population,param2max);
  trynumber1=2;
  if (end()) return;
  int param1max=parambest;
  value2max=valuebest;
  if ((value2max - problem->lower_bound) > (value2min - problem->lower_bound) * (1 - epsilon))  // epsilon sur la moyenne pour eviter d'agrandir le voisinage quand fin jamais atteinte
    {param2best = param2min; value2best= value2min;
    if (param2min==bound2min) minfound=1;
    param2min=param2min/2; bestparam1=param1min;
    while (!minfound && trynumber1 < maxtuning2 && !end())
      {trynumber1++;
      simplerun(algo,problem,population,param2min);
      if (valuebest > value2best)
	{minfound=1;}
      else if (param2min==bound2min)
	{minfound=1; param2best=param2min; bestparam1 = parambest; value2best=valuebest;}
      else
	{value2best=valuebest ; param2max= param2best;
	param2best=param2min; bestparam1=parambest; param2min = param2min/2; }
     }
    }
  else
    {param2best= param2max; value2best = value2max;
    param2max=param2max*2; bestparam1=param1max;
    while (!minfound && trynumber1 < maxtuning2 && !end() && param2max <= bound2max)
      {trynumber1++;
      simplerun(algo,problem,population,param2max);
      if ((valuebest - problem->lower_bound) > (value2best - problem->lower_bound)  * ( 1 - epsilon))
	{minfound=1; }
      else
	{ value2best= valuebest ; param2min= param2best;
	param2best=param2max; bestparam1=parambest; param2max = param2max*2; 
	if (param2max > bound2max) param2max=bound2max;}

      }
    }
   while (trynumber1 < maxtuning2 && param2best - param2min > 1 && param2max-param2best > 1 &&  !end() && param2best !=0)
     {trynumber1++;
     *ofile << " param2best " << param2best << " param2min " << param2min << " param2max " << param2max << endl;
     param2dico1= (param2best+param2min)/2;
     simplerun(algo,problem,population,param2dico1);
     if (end()) break;
     int param1dico1=parambest;
     value2dico1=valuebest;
     //     *ofile << " valeur dico 1 " << value2dico1 << endl;
     trynumber1++;
     param2dico2= (param2best+param2max)/2;
     simplerun(algo,problem,population,param2dico2);
     int param1dico2=parambest;
     value2dico2=valuebest;
     //     *ofile << " valeur dico 2 " << value2dico2 << endl;
     if (value2dico1 < value2best && value2dico1 < value2dico2)
       {value2max=value2best; param2max=param2best;
       value2best=value2dico1; param2best=param2dico1; bestparam1=param1dico1; }
     else if
       (value2dico2 < value2best)
       { value2min=value2best; param2min=param2best;
       value2best=value2dico2; param2best=param2dico2; bestparam1=param1dico2;}
     else {param2min=param2dico1; param2max=param2dico2;}
     }
   parambest=bestparam1; // le meilleur reglage de param1 correspondant au meilleur de param2
   *ofile << " best param2 " << param2best << " best param1 " << parambest << " best value " << value2best << endl;
   *ofile << " Fin reglage " << Statistiques->total_execution_time << endl;
}


void DoubleTuning::simplerun (LSAlgorithm* algo, OpProblem* problem, Configuration** population, int param2)
{ int minfound =0;
  trynumber = 0;
  double epsilon = 0.01;

  parammin=parambest;
  if (parammin==0) parammin=1; // pour ne pas stagner en 0
  parammax=2*parammin;
  *ofile << "reglage boucle intérieure " << "  param extérieur = " << param2 << " param interieur " << parambest << endl;
  //  if (algo->methodname=="tabu" || algo->methodname == "tabugreedy") epsilon=0.1;
  //  if (algo->methodname=="idwminmax") epsilon=0.2;
  if (algo->methodname=="idwaminmax" || algo->methodname =="idwbminmax" || algo->methodname == "idwupk" ||
      algo->methodname=="idwgrupk") 
    {boundmax=param2; 
    if (parammax > boundmax) 
      {parammax=boundmax; parammin=boundmax/2;}}
  onerun2param (problem,algo,population,parammin,param2);
  trynumber=1;
  if (end()) return;
  valuemin=value;
  onerun2param(problem,algo,population,parammax,param2);
  trynumber=2;
  if (end()) return;
  valuemax=value;
  if ((valuemax - problem->lower_bound) > (valuemin - problem->lower_bound) * (1 - epsilon))  // epsilon sur la moyenne pour eviter d'agrandir le voisinage quand fin jamais atteinte
    {parambest = parammin; valuebest= valuemin;
    if (parammin <= boundmin) minfound=1;
    parammin=parammin/2;
    while (!minfound && trynumber < maxtuning && !end())
      {trynumber++;
      onerun2param(problem,algo,population,parammin,param2);
      if (value > valuebest)
	{minfound=1;}
      else if (parammin ==boundmin)
	{minfound=1; parambest=parammin; valuebest=value;}
      else
	{valuebest=value ; parammax= parambest;
	parambest=parammin; parammin = parammin/2; }
     }
    }
  else
    {parambest= parammax; valuebest = valuemax;
    if (parammax < boundmax)
      {    parammax=parammax*2;
      while (!minfound && trynumber < maxtuning && !end() && parammax <= boundmax)
	{trynumber++;
	onerun2param(problem,algo,population,parammax,param2);
	if ((value - problem->lower_bound) > (valuebest - problem->lower_bound) * ( 1 - epsilon))
	  {minfound=1; }
	else
	  { valuebest= value ; parammin= parambest;
	  parambest=parammax; parammax = parammax*2; 
	  if (parammax > boundmax) parammax=boundmax;}
	}
      }
    }
   if (parammax > boundmax) parammax=boundmax;
   while (trynumber < maxtuning && parambest - parammin > 1 && parammax-parambest > 1 &&  !end() && parambest !=0)
     {trynumber++;
     *ofile << " parambest " << parambest << " parammin " << parammin << " parammax " << parammax << endl;
     paramdico1= (parambest+parammin)/2;
     onerun2param(problem,algo,population,paramdico1,param2);
     if (end()) break;
     valuedico1=value;
     trynumber++;
     paramdico2= (parambest+parammax)/2;
     if (parammax > boundmax) parammax=boundmax;
     onerun2param(problem,algo,population,paramdico2,param2);
     valuedico2=value;
     if (valuedico1 < valuebest && valuedico1 < valuedico2)
       {valuemax=valuebest; parammax=parambest;
       valuebest=valuedico1; parambest=paramdico1;}
     else if
       (valuedico2 < valuebest)
       { valuemin=valuebest; parammin=parambest;
       valuebest=valuedico2; parambest=paramdico2;}
     else {parammin=paramdico1; parammax=paramdico2;}
     }
   *ofile << " param2 " << param2 << " best param1 " << parambest << " best value " << valuebest << endl;
   *ofile << " Fin reglage interieur " << Statistiques->total_execution_time << endl;
}

// premier essai : pour étalonner le temps
int Tuning::firsttry()
{return (trynumber==0);}

int DoubleTuning::firsttry()
{return ((trynumber ==0) && (trynumber1==0));}

// condition arret du reglage : tous essais ont trouvé ou 1 essai a trouvé (cas arret 1 solution)
int Tuning::end(){
  //  *ofile << "appel de end " << stop << endl;
   return (Statistiques->trouve[Statistiques->current_pb] == nbtries || stop);}




// un reglage : le 1er essai sert à regler la longueur de marche pour effectuer l'essai dans le 
// temps de référence (sauf pour le recuit simulé)
void Tuning::onerun (OpProblem* problem, LSAlgorithm* algo, Configuration** population)
{
  
  Statistiques->init_run();
  if (!firsttry() 
       && (algo->methodname != "simann")
      )
    {executer_essai (problem,algo,population,1,seed,0);
     float trytime=Statistiques->execution_time_try[0];
     if (Statistiques->trouve[Statistiques->current_pb]==0)
       {algo->walklength = (int) (algo->walklength * referencetime / trytime);
       }

     *ofile << " longueur marche essai "  << algo->walklength << endl;}
  
  Statistiques->init_run();     

  for (int nessai = 0;nessai< nbtries ; nessai++)
    {
     executer_essai (problem,algo,population,1,seed,nessai);
    if (Statistiques->stop_trouve && Statistiques->trouve[Statistiques->current_pb])
      {stop=1;break;}
    }
  Statistiques->current_try++; 
  //ecriture_stat_probleme();
  if (firsttry())
    //    referencetime = Statistiques->execution_time_try[0];
    referencetime = Statistiques->average_execution_time;
   value=Statistiques->cost_meanvalue;
  
}
  

// le reglage automatique d'un parametre pour une longueur de marche donnée + le lancement de l'algorithme avec
// ce réglage
int autotuning(LSAlgorithm* algo, Configuration** population, OpProblem* problem, int graine1, int nbessais, int tuninit
, int tuningwalkrate, int tuningmaxtries)
{
   int bmin=0; int bmax=RAND_MAX;
   if (algo->methodname == "idwa" || algo->methodname =="idwb" || algo->methodname =="idwbsn") bmin=1;
   else if (algo->methodname =="idwgra"   || algo->methodname == "idwgrb") bmin=5; // voisinage au moins de longueur 1
   int walklength1 = algo->walklength;
   algo->walklength= walklength1/tuningwalkrate; // marche divisee par tuningwalkrate pour le reglage
   *ofile << " reglage parametre : longueur marche " << algo->walklength << endl;
   Tuning tun (tuninit,graine1, nbessais , tuningmaxtries,bmin, bmax); // objet Tuning
   tun.run(algo,problem, population);  // reglage 
  // mise en place de la recherche avec la valeur du parametre determinee par le reglage 
   algo->walklength= walklength1;
   tun.trynumber=0;
   if (!tun.end())
     tun.onerunparam(problem,algo,population,tun.parambest);
   return tun.parambest;
}

int nb_parameters(LSAlgorithm* algo)
   {if (algo->methodname == "tabu" || algo->methodname == "incrtabu" 
	|| algo->methodname == "idwatabu" || algo->methodname == "idwbtabu"
	|| algo->methodname =="idwgratabu" || algo->methodname =="idwgrbtabu" 
	|| algo->methodname == "idwaincrtabu" || algo->methodname == "idwbincrtabu"
	|| algo->methodname =="idwgraincrtabu" || algo->methodname =="idwgrbincrtabu" 
	|| algo->methodname == "idwaminmax" || algo->methodname== "idwbminmax" 
	|| algo->methodname == "idwupk" || algo->methodname =="idwgrupk")
     return 2;
   else return 1;
   }



// algo complet 
void autosolving (LSAlgorithm* algo, Configuration** population, OpProblem* problem, int numpb, int graine1, int nbessais, double maxtime, int initwalklength)
  {
    if (nb_parameters (algo) == 2)
      autosolving2(algo,population,problem,numpb,graine1,nbessais,maxtime,initwalklength);
    else
      autosolving1(algo,population,problem,numpb,graine1,nbessais,maxtime,initwalklength);
  }

// algo complet resolution avec reglage algo à un parametre
void autosolving1 (LSAlgorithm* algo, Configuration** population, OpProblem* problem, int numpb, int graine1, int nbessais, double maxtime, int initwalklength)
{int pas = 4; // le pas de multiplication de longueur de marche
 int parameter=100 ; // valeur initiale 100
 int tuningwalkrate = 50; // la division de la marche pour le réglage
 int tuningmaxtries = 10; // le maximum d'essais pour la dichotomie
 algo->walklength=initwalklength;
 while (!endsolve(maxtime,nbessais, numpb))
   { *ofile << " ESSAI RESOLUTION :  longueur marche : " << algo->walklength << endl;
   int tuninit=parameter;
   if (tuninit==0) tuninit=1; // pour ne pas stagner avec la valeur 0
   parameter=  autotuning (algo, population, problem, graine1, nbessais, tuninit, tuningwalkrate, tuningmaxtries);
   if (algo->walklength > RAND_MAX/pas) break;
   else algo->walklength = pas * algo->walklength;
   *ofile << " Temps total utilisé " << Statistiques->total_execution_time << endl; 
   }
}

// réglage de 2 parametres  : parameter1 et parameter2 (en sortie) + lancement algo avec ce reglage
void autotuning2(LSAlgorithm* algo, Configuration** population, OpProblem* problem, int graine1, int nbessais, int&
parameter1, int& parameter2 ,int tuningwalkrate, int tuningmaxtries, int tuningmaxtries2)
{
   int bmin=1; 
   if (algo->methodname == "idwgrupk") bmin=5;
   int b2min =1;
   int bmax=RAND_MAX;
   int b2max = RAND_MAX;
   int walklength1 = algo->walklength;
   algo->walklength= walklength1/tuningwalkrate; // marche divisee par tuningwalkrate pour le reglage
   *ofile << " reglage parametre : longueur marche " << algo->walklength << endl;
   DoubleTuning tun (parameter1, parameter2,graine1, nbessais , tuningmaxtries, tuningmaxtries2, bmin,b2min, bmax, b2max); // objet DoubleTuning
   tun.run(algo,problem, population);  // reglage 
  // mise en place de la recherche avec la valeur du parametre determinee par le reglage 
   algo->walklength= walklength1;
   tun.trynumber=0;tun.trynumber1=0;
   if (!tun.end())
       {tun.onerun2param(problem,algo,population, tun.parambest, tun.param2best);}
   parameter1=tun.parambest;
   parameter2=tun.param2best;
}

// algo complet resolution avec reglage algo à deux paramètres
void autosolving2 (LSAlgorithm* algo, Configuration** population, OpProblem* problem, int numpb, int graine1, int nbessais, double maxtime, int initwalklength)
{int pas = 4; // le pas de multiplication de longueur de marche
 int parameter1=10 ; // valeur initiale 10 pour la liste taboue -> mettre une méthode ; valeur initiale min
                    // pour idwminmax et spare-neighbor pour idwupk et idwupgrk
 int parameter2=100 ; // valeur initiale 100
 int tuningwalkrate = 50; // la division de la marche pour le réglage
 int tuningmaxtries = 5; // le maximum d'essais pour les dichotomies  (parametre intérieur)
 int tuningmaxtries2 = 5; // le maximum d'essais pour les dichotomies (parametre extérieur)
 algo->walklength=initwalklength;
 while (!endsolve(maxtime,nbessais, numpb))
   { *ofile << " ESSAI RESOLUTION :  longueur marche : " << algo->walklength << endl;
   if (parameter1==0) parameter1=1; // pour ne pas stagner avec la valeur 0
   if (parameter2==0) parameter2=1; // pour ne pas stagner avec la valeur 0
   autotuning2 (algo, population, problem, graine1, nbessais, parameter1, parameter2, tuningwalkrate, tuningmaxtries, tuningmaxtries2);
   if (algo->walklength > RAND_MAX/pas) break;
   else algo->walklength = pas * algo->walklength;
   *ofile << " Temps total utilisé " << Statistiques->total_execution_time << endl; 
   }
}


void DoubleTuning::doubleparameterwrite(LSAlgorithm* algo, int parameter1, int parameter2)
{ *ofile << "methode " << algo->methodname << " " ;
 if (algo->methodname=="idwatabu" || algo->methodname=="idwbtabu" ||
     algo->methodname=="idwaincrtabu" || algo->methodname=="idwbincrtabu"  )
   *ofile << " valeur parametre voisinage " << parameter2 << " longueur liste taboue " << parameter1 << endl;
else if (algo->methodname=="idwgratabu" || algo->methodname=="idwgrbtabu" ||
	 algo->methodname=="idwgraincrtabu" || algo->methodname=="idwgrbincrtabu" ||
	 algo->methodname == "tabu" || algo->methodname =="incrtabu")
   *ofile << " valeur parametre voisinage " << parameter2/5 << " longueur liste taboue " << parameter1 << endl;
else if (algo->methodname =="idwaminmax" || algo->methodname =="idwbminmax")
  *ofile << " valeur parametre voisinage  max voisins " << parameter2 << " min voisins " << parameter1 << endl;
else if (algo->methodname =="idwupk")
    *ofile << " valeur parametre voisinage  max voisins " << parameter2 << " remontée " << parameter1 << endl;
else if (algo->methodname =="idwgrupk")
  *ofile << " valeur parametre voisinage  max voisins " << parameter2/5 << " remontée " << parameter1/5 << endl;
}
void Tuning::oneparameterwrite(LSAlgorithm* algo, int parameter)
{
  *ofile << "methode " << algo->methodname << " " ;
  if (algo->methodname =="idwa" || algo->methodname == "idwb" || algo->methodname == "idwbsn")
    *ofile << " valeur parametre voisinage " << parameter << endl;
  else if (algo->methodname =="idwgra" || algo->methodname == "idwgrb")
    *ofile << " valeur parametre voisinage " << parameter/5 << endl;
  else if (algo->methodname== "metropolis")
    *ofile << " valeur temperature " << ((double) parameter) /100  << endl;
  else if ( algo->methodname =="simann")
    *ofile << " valeur temperature initiale " << ((double) parameter /100)  << endl;}

// pour le réglage automatique d'un paramètre
// indique selon l'algo le parametre à regler et lance le reglage (onerun)
void Tuning::onerunparam (OpProblem* problem,LSAlgorithm* algo, Configuration** population, int parameter )
{ 
  if (algo->methodname =="idwa" || algo->methodname =="idwbsn")
    {algo->nbhsearch->maxneighbors = parameter;
    algo->nbhsearch->minneighbors = 1;
    algo->nbhsearch->finished = 1;}
  else if (algo->methodname =="idwb")
    {algo->nbhsearch->maxneighbors = parameter;
    algo->nbhsearch->minneighbors = 1;
    algo->nbhsearch->finished = parameter;}
  else  if (algo->methodname =="idwgra")   // minvoisins=maxvoisins - sortie any (finished=1)
    {algo->nbhsearch->maxneighbors = parameter/5; // on commence à 20
    algo->nbhsearch->minneighbors = parameter/5;
    algo->nbhsearch->finished = 1;}

  else  if (algo->methodname =="idwgrb")   // minvoisins=maxvoisins - sortie best (finished = max)
    {algo->nbhsearch->maxneighbors = parameter/5; // on commence à 20
    algo->nbhsearch->minneighbors = parameter/5;
    algo->nbhsearch->finished = parameter/5;}

  else if (algo->methodname =="metropolis") 
    {algo->mheur->adjustparameter (parameter);}
  
  else if (algo->methodname == "simann")
    {((SimulatedAnnealing*)(algo->mheur))->walklength=algo->walklength;
    algo->mheur->adjustparameter (parameter);}
  oneparameterwrite(algo,parameter);  
  onerun( problem, algo, population);
  oneparameterwrite(algo,parameter);  
}

// pour le reglage automatique de 2 paramètres 
void DoubleTuning::onerun2param (OpProblem* problem, LSAlgorithm* algo,Configuration** population, int parameter1, int parameter2 )
{ 
  if (algo->methodname=="idwatabu" || algo->methodname=="idwaincrtabu")
    {algo->nbhsearch->maxneighbors = parameter2;
    algo->mheur->adjustparameter (parameter1);
    algo->nbhsearch->finished=1;}
  else if (algo->methodname=="idwbtabu"|| algo->methodname=="idwbincrtabu")
    {algo->nbhsearch->maxneighbors = parameter2;
    algo->mheur->adjustparameter (parameter1);
    algo->nbhsearch->finished=parameter2;}
  else if (algo->methodname=="idwgratabu"|| algo->methodname=="idwgraincrtabu")
    {algo->nbhsearch->maxneighbors = parameter2/5;
    algo->nbhsearch->minneighbors = parameter2/5;
    algo->mheur->adjustparameter (parameter1);
    algo->nbhsearch->finished=1;}
  else if (algo->methodname=="idwgrbtabu"|| algo->methodname=="idwgrbincrtabu")
    {algo->nbhsearch->maxneighbors = parameter2/5;
    algo->nbhsearch->minneighbors = parameter2/5;
    algo->mheur->adjustparameter (parameter1);
    algo->nbhsearch->finished=parameter2/5;}
  else if (algo->methodname=="tabu" || algo->methodname =="incrtabu")
    {algo->nbhsearch->maxneighbors = parameter2/5; // on commence max = 20
    algo->nbhsearch->minneighbors = parameter2/5;  // min = max
    algo->nbhsearch->finished= parameter2/5;      //  sortie = max
    algo->mheur->adjustparameter (parameter1);}

  else if (algo->methodname =="idwaminmax")
    {algo->nbhsearch->maxneighbors = parameter2;
    algo->nbhsearch->minneighbors = parameter1;
    algo->nbhsearch->finished=1;}
  else if (algo->methodname =="idwbminmax")
    {algo->nbhsearch->maxneighbors = parameter2;
    algo->nbhsearch->minneighbors = parameter1;
    algo->nbhsearch->finished=parameter2;}
  else if (algo->methodname =="idwupk")
    {algo->nbhsearch->maxneighbors = parameter2;
    algo->nbhsearch->minneighbors=1;
    algo->nbhsearch->finished = parameter1;}

  else if (algo->methodname =="idwgrupk")
    {algo->nbhsearch->maxneighbors = parameter2/5;
    algo->nbhsearch->minneighbors=parameter2/5;
    algo->nbhsearch->finished = parameter1/5;}
  doubleparameterwrite(algo,parameter1,parameter2);
  onerun( problem, algo, population);
  doubleparameterwrite(algo,parameter1,parameter2);
}


