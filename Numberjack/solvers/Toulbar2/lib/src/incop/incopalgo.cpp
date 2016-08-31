/* le fichier comportant les méthodes des différents algorithmes (gww
   et recherche locale) , des métaheuristiques, des classes de
   configurations, */


#include <list>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

using namespace std;
#include <fstream>
#include <math.h>
#include <unistd.h>

#include "timer.h"
#include "incop.h"
#include "incoputil.h"


// les variables globales
extern Stat_GWW *Statistiques;  // defini dans incoputil
extern int TRACEMODE;    // indicateur de trace
extern ofstream *ofile;  // le fichier de sortie



// les constructeurs et destructeurs

LSAlgorithm::LSAlgorithm(int nbmov)
{
	walklength = nbmov;
}

LSAlgorithmGWW::LSAlgorithmGWW(int nbmov) : LSAlgorithm(nbmov) {;}

LSAlgorithm::~LSAlgorithm()
{
	delete nbhsearch;
	delete mheur;
}

NeighborhoodSearch::NeighborhoodSearch(int nbhmin, int nbhmax, int endnbh, int varconf, int valconf , double nbhr)
{
	minneighbors = nbhmin, maxneighbors = nbhmax,
	finished = endnbh, var_conflict = varconf, val_conflict = valconf;
	nbhrate = nbhr;
}

TabuSearch::TabuSearch(int tabul)
{
	list<Move *> L;
	tabulength = tabul;
	move_list = L;
}

IncrTabuSearch::IncrTabuSearch(int tabul)
	: TabuSearch(tabul)
{
	nbiter = 0;
}

IncrTabuGreedySearch::IncrTabuGreedySearch(int tabul)
	: IncrTabuSearch(tabul)
{
	;
}

TabuGreedySearch::TabuGreedySearch(int tabul)
	: TabuSearch(tabul) {;}


Metropolis::Metropolis(double temp)
{temperature = temp;}

/* constructeur : calcul du pas constant delta de baisse du seuil */
ThresholdAccepting::ThresholdAccepting(double maxthreshhold, int walklength)
{
	thresholdinit = maxthreshhold;
	delta = thresholdinit / walklength;
}

/* constructeur : calcul du pas constant delta de baisse de température */
SimulatedAnnealing::SimulatedAnnealing(double initialtemperature, int walkl)
{
	inittemperature = initialtemperature;
	walklength = walkl;
	delta = inittemperature / walklength;
}

void SimulatedAnnealing::adjustparameter(int parameter)
{
	inittemperature = ((double) parameter / 100) ;
	delta = inittemperature / walklength;
}

RandomSearch::RandomSearch() {;}

GreedySearch::GreedySearch() {;}

GWWAlgorithm::~GWWAlgorithm() {delete walkalgorithm;}

StandardGWWAlgorithm::StandardGWWAlgorithm(int popsize, int grtest, int derniermouv, int elitisme, int stop,
        double threshdescent, Long thresmin)
{
	populationsize = popsize;
	regrouptest = grtest;
	thresholddescent = threshdescent;
	thresholdmin = thresmin;
	lastmovedescent = derniermouv;
	elitism = elitisme;
	nomovestop = stop;
	nbiteration = RAND_MAX;
}

FastStandardGWWAlgorithm::FastStandardGWWAlgorithm(int popsize, int grtest, int derniermouv, int elitisme, int stop,						   double threshdescent, Long thresmin):
	StandardGWWAlgorithm(popsize, grtest, derniermouv, elitisme, stop, threshdescent, thresmin) {;}



NothresholdGWWAlgorithm::NothresholdGWWAlgorithm(int popsize, int grtest, int derniermouv, int elitisme, int stop,
        int killed, int nbiter)
{
	populationsize = popsize;
	regrouptest = grtest;
	lastmovedescent = derniermouv;
	elitism = elitisme;
	nomovestop = stop;
	nbkilled = killed;
	nbiteration = nbiter;
}

AdaptiveGWWAlgorithm::AdaptiveGWWAlgorithm(int popsize, int grtest, int derniermouv, int elitisme, int stop, int killed)
{
	populationsize = popsize;
	regrouptest = grtest;
	lastmovedescent = derniermouv;
	elitism = elitisme;
	nomovestop = stop;
	nbiteration = RAND_MAX;
	nbkilled = killed;
}


FastAdaptGWWAlgorithm::FastAdaptGWWAlgorithm(int popsize, int grtest, int derniermouv, int elitisme, int stop, int killed, int maxkilled, double threshdescent):
	AdaptiveGWWAlgorithm(popsize, grtest, derniermouv, elitisme, stop, killed)
{thresholddescent = threshdescent; nbmaxkilled = maxkilled;}

MedianAdaptGWWAlgorithm::MedianAdaptGWWAlgorithm(int popsize, int grtest, int derniermouv, int elitisme, int stop, double medianrate):
	AdaptiveGWWAlgorithm(popsize, grtest, derniermouv, elitisme, stop, 1)
{mediandescent = medianrate;}

BestAdaptGWWAlgorithm::BestAdaptGWWAlgorithm(int popsize, int grtest, int derniermouv, int elitisme, int stop, double bestrate):
	AdaptiveGWWAlgorithm(popsize, grtest, derniermouv, elitisme, stop, 1)
{bestdescent = bestrate;}



/* **********************************************************************************************************
 les différentes classes de configuration correspondant aux différentes formes d'incrémentalité
*/

Configuration::Configuration() {;}
Configuration::Configuration(int nbv)
{
	nbvar = nbv;
	config = new int[nbv];
}


CSPConfiguration::CSPConfiguration(int nbv, int domsize) : Configuration(nbv)
{
	domainsize = domsize;
	trynumber = 0;
}


IncrCSPConfiguration::IncrCSPConfiguration(int nbv, int domsize) : CSPConfiguration(nbv, domsize)
{
	tabconflicts = new Long[nbv];
	for (int i = 0; i < nbv; i++)
	{tabconflicts[i] = 0;}
}

FullincrCSPConfiguration::FullincrCSPConfiguration(int nbv, int domsize) : CSPConfiguration(nbv, domsize)
{
	tabconflictsize = domsize;
	tabconflicts = new Long*[nbv];
	for (int i = 0; i < nbv; i++) {
		tabconflicts[i] = new Long[domsize];
		for (int j = 0; j < domsize; j++)
			tabconflicts[i][j] = 0;
	}
}





Configuration::~Configuration() { delete config;}
FullincrCSPConfiguration::~FullincrCSPConfiguration()
{
	for (int i = 0; i < nbvar; i++)
		delete tabconflicts[i];
	delete tabconflicts;
}

IncrCSPConfiguration::~IncrCSPConfiguration()
{delete tabconflicts;}


void Configuration::init_conflicts() {;}
void Configuration::incr_conflicts(int var, int val, int index, Long incr) {;}
void Configuration::set_conflicts(int var, int val, int index, Long nbconf) {;}
Long Configuration::get_conflicts(int var, int val, int index) {return 0;}
Long Configuration::get_conflicts_problem(OpProblem *problem, int var, int val)
{
	return problem->compute_conflict(this, var, val);
}

void Configuration::update_conflicts(OpProblem *problem, Move *move)
{;}


/* L'implantation semi-incrementale : seuls les nb de conflits des valeurs courantes sont stockés dans
   tabconflicts (tableau a une dimension) */

void IncrCSPConfiguration::set_conflicts(int var, int val, int index, Long nbconf)
{if (config[var] == val) tabconflicts[var] = nbconf;}

// methode utile ??
Long IncrCSPConfiguration::get_conflicts(int var, int val, int index)
{
	if (config[var] == val) return tabconflicts[var];
	else return 0;
}

/* la valeur courante est stockée , pour les autres valeurs, il faut calculer */
Long IncrCSPConfiguration::get_conflicts_problem(OpProblem *problem, int var, int val)
{
	if (config[var] == val) return tabconflicts[var];
	else return  problem->compute_conflict(this, var, val);
}


void IncrCSPConfiguration::set_variableconflicts(int var, int nbconf)
{tabconflicts[var] = nbconf;}


void IncrCSPConfiguration::incr_conflicts(int var, int val, int index, Long incr)
{if (config[var] == val) tabconflicts[var] += incr;}



void IncrCSPConfiguration::init_conflicts()
{
	for (int i = 0 ; i < nbvar; i++)
		tabconflicts[i] = 0;
}


void IncrCSPConfiguration::update_conflicts(OpProblem *problem, Move *move)
{problem->incr_update_conflicts(this, move);}


/* L'implantation toute-incrementale :  les nb de conflits de toutes les valeurs sont stockés dans
   tabconflicts (tableau à deux dimensions : variable, valeur) */


void FullincrCSPConfiguration::init_conflicts()
{
	for (int i = 0 ; i < nbvar; i++)
		for (int j = 0; j < tabconflictsize; j++)
			tabconflicts[i][j] = 0;
}




Long FullincrCSPConfiguration::get_conflicts_problem(OpProblem *problem, int var, int val)
{
	return tabconflicts[var][problem->value2index(val, var)];
}

Long FullincrCSPConfiguration::get_conflicts(int var, int val, int index)
{return tabconflicts[var][index];}



void FullincrCSPConfiguration::set_conflicts(int var, int val, int index, Long nbconf)
{tabconflicts[var][index] = nbconf;}


void FullincrCSPConfiguration::incr_conflicts(int var, int val, int index, Long incr)
{
	tabconflicts[var][index] += incr;
}


void FullincrCSPConfiguration::update_conflicts(OpProblem *problem, Move *move)
{
	problem->fullincr_update_conflicts(this, move);
}



// copie la config2 dans la config1 (recopie les tableaux de conflits, l'évaluation)

void Configuration::copy_element(Configuration *config2)
{
	for (int i = 0; i < nbvar ; i++)
		config[i] = config2->config[i];
	valuation = config2->valuation;
}


void IncrCSPConfiguration::copy_element(Configuration *config2)
{
	Configuration::copy_element(config2);
	for (int i = 0; i < nbvar ; i++)
		tabconflicts[i] = ((IncrCSPConfiguration *)config2)->tabconflicts[i];
}


void FullincrCSPConfiguration::copy_element(Configuration *config2)
{
	Configuration::copy_element(config2);
	for (int i = 0; i < nbvar ; i++)
		for (int j = 0; j < tabconflictsize ; j++)
			tabconflicts[i][j] = ((FullincrCSPConfiguration *)config2)->tabconflicts[i][j];
}


// ******************************************************************************
int NeighborhoodSearch::returnbestmove()
{return ((minneighbors > 1) || (finished >= 2));}


// ********************************** LES ALGORITHMES************************************************


void IncompleteAlgorithm::randomwalk
(OpProblem *problem, Configuration *configuration)
{;}

void IncompleteAlgorithm::run(OpProblem *problem, Configuration **population)
{;}

//---------------------------LSALGORITHM----------------------------------------------------------


/*  Test pour savoir si le mouvement courant peut conduire à une configuration meilleure que toutes
celles vues jusqu'à présent.
Cette meilleure valeur est stockée dans l'objet Statistiques.
(ATTENTION : Statistiques : variable globale ; la meilleure configuration est aussi
stockée dans le champ best_config du problème, mais cette mise à jour là
ne peut s'effectuer qu'une fois le mouvement effectué )
en cas de meilleur trouvé, cette valeur est mise à jour dans l'objet Statistiques et une impression
est eventuellement faite */

int LSAlgorithm::test_bestfound(Move *move)
{
	if (move->valuation < Statistiques->cost_try[Statistiques->current_try]) {
		if (TRACEMODE) {
#ifndef WINDOWS
			stop_timers(VIRTUAL);
			Statistiques->execution_time_try[Statistiques->current_try] += VIRTUAL_TIMELAPSE;
#else // timer not supported under windows 
			Statistiques->execution_time_try[Statistiques->current_try] ;
#endif
			ecriture_meilleure_valeur(methodname, move->valuation, threshold, nbmoves, nbhsearch->maxneighbors);
#ifndef WINDOWS
			start_timers();
#endif
		}
		Statistiques->cost_try[Statistiques->current_try] = move->valuation;

		return 1;
	}

	return 0;
}

/* Analyse des evaluations des nbmove derniers mouvements : renvoie l'écart entre valeurs maximale et minimale
si la derniere valeur est inferieure a la première , renvoie -1 sinon */
Long analyselastvaluations(int nbmove)
{
	Long valmax = 0, valmin = LONG_MAX;
	for (int i = Statistiques->costvalues.size() - nbmove;
	        i < (int) Statistiques->costvalues.size(); i++) {
		if (Statistiques->costvalues[i] > valmax)
			valmax = Statistiques->costvalues[i];
		if (Statistiques->costvalues[i] < valmin)
			valmin = Statistiques->costvalues[i];
	}
	if (Statistiques->costvalues[Statistiques->costvalues.size() - 1]
	        >= Statistiques->costvalues[Statistiques->costvalues.size() - nbmove])
		return (valmax - valmin);
	else return (-1);
}



/* voisinage statique : pas de reglage dynamique */
void   NeighborhoodSearch::dynamicmaxneighbors(int &maxneigh, int &minneigh, int nbmoves)
{;}

void NeighborhoodSearch::initsearch() {;}

void DynamicNeighborhoodSearch::initsearch()
{
	maxneighbors = initmaxneighbors;
	minneighbors = initminneighbors;
}

/* Réglage dynamique du parametre maxneighbors en fonction des derniers couts*/

DynamicNeighborhoodSearch::DynamicNeighborhoodSearch(int minneigh, int maxneigh, int finish, int var_conf, int val_conf, double nbbr) : NeighborhoodSearch(minneigh, maxneigh, finish, var_conf, val_conf, nbbr)
{
	initmaxneighbors = maxneigh;
	initminneighbors = minneigh;
	adjustperiod = 50;
}

void DynamicNeighborhoodSearch::dynamicmaxneighbors(int &maxneigh, int &minneigh, int nbmoves)
{
	if ((nbmoves > 0) && (nbmoves % adjustperiod == 0)) {
		if (analyselastvaluations(10) == 0) {
			if (maxneigh <= 20)
			{if (maxneigh > 2) {maxneigh--; maxneighbors--;}}
			else {
				maxneigh = (int)(0.9 * maxneigh);
				maxneighbors = maxneigh;
			}
		} else if (analyselastvaluations(10) >= 2) {
			if (maxneigh <= 20)
			{maxneigh++; maxneighbors++;}
			else {
				if (maxneigh < 10000) {
					maxneigh = (int)(1.1 * maxneigh);
					maxneighbors = maxneigh;
				}
			}
		}

		if (minneigh > maxneigh) minneigh = maxneigh;
	}
}

DynamicSpareneighbor::DynamicSpareneighbor(int minneigh, int maxneigh, int finish, int var_conf, int val_conf, double nbbr) : NeighborhoodSearch(minneigh, maxneigh, finish, var_conf, val_conf, nbbr) { nbmovesdown = 0;}

void DynamicSpareneighbor::spareneighboradjust(Configuration *config, Move *move)
{
	if (config->valuation < move->valuation) {
		if (Statistiques->nb_moves_up[Statistiques->current_try] == 0)
		{nbmovesdown = Statistiques->nb_moves_down[Statistiques->current_try];}
		finished = (int)(finished / (1.2));
		if (finished == 0) finished = 1;
		//          *ofile << " butee spare neighbor " << finished <<  " " << move->valuation<<  endl;
	} else if (config->valuation > move->valuation) {
		double rate = 1.2 * ((double)(Statistiques->nb_moves_up[Statistiques->current_try]) /
		                     (Statistiques->nb_moves_down[Statistiques->current_try] + 1 - nbmovesdown));
		if (finished * rate - finished > 1)
			finished = (int)(finished * rate);
		else finished++;
		if (finished > maxneighbors) finished = maxneighbors;
		//    *ofile << " descente " << finished <<  " " << move->valuation << endl;
	} else {
		//  *ofile << " plateau " << finished <<  " " << move->valuation << endl;
	}

}

void NeighborhoodSearch::adjust_neighborhood(Configuration *configuration, OpProblem *problem,
        int &maxneigh, int &minneigh, int nbmoves)
{
	int nbvariables;
	minneigh = minneighbors;
	if (nbhrate) {
		if (var_conflict)
			nbvariables = configuration->var_conflict.size();
		else
			nbvariables = configuration->nbvar;
		maxneigh = (int)(nbhrate * nbvariables * (problem->domainsize - 1));
		if (minneigh > maxneigh) minneigh = maxneigh;
	} else
		maxneigh = maxneighbors;
	// reglage dynamique possible si TRACEMODE=2 (information stockée pour ce reglage)
	if (TRACEMODE == 2)
		dynamicmaxneighbors(maxneigh, minneigh, nbmoves);
	if (var_conflict)
		problem->adjust_parameters(configuration, maxneigh, minneigh);
}



/*
Algo de base qui sélectionne, teste et effectue un mouvement à partir de la configuration courante
Effectue le mouvement et renvoie 1 si un mvt a été effectué et 0 si aucun mouvement ne l'a été
*/

int LSAlgorithm::configurationmove(OpProblem *problem, Configuration *configuration)
{
	int nbtries = 0;

	Move *bestmove = problem->bestmove;
	Move *currentmove = problem->currentmove;
	Move *firstmove = problem->firstmove;

	bestmove->valuation = LONG_MAX;

	int bestfound = 0;
	int feasiblefound = 0;
	int accepted = 0;
	int nbrefused = 0;

	int maxneigh;
	int minneigh;
	// cas des variables en conflit : on calcule l'ensemble des variables en conflit
	if (nbhsearch->var_conflict)
		problem->compute_var_conflict(configuration);

	// calcul des parametres maxneigh et minneigh à partir du voisinage
	nbhsearch->adjust_neighborhood(configuration, problem, maxneigh, minneigh, nbmoves);

	// boucle sur les voisins essayés
	while ((nbtries < minneigh) || ((nbtries < maxneigh) && !accepted)) {
		problem->next_move(configuration, currentmove, nbhsearch);
		if (isfeasible(currentmove)) {
			// première valeur faisable ? (sous seuil)
			if (!feasiblefound) {
				firstmove->copymove(currentmove);
				feasiblefound = 1;
			}

			// meilleur que le meilleur trouve jusqu'à présent ?
			if (test_bestfound(currentmove))
				bestfound = 1;

			// test d'acceptation de la métaheuristique
			if (bestfound //(aspiration si meilleur_trouve)
			        || ((!accepted) && mheur->acceptance(currentmove, configuration))
			        || (accepted && minneigh > 1 && currentmove->valuation < bestmove->valuation && mheur->acceptance(currentmove, configuration))
			   ) {
				if (!accepted)  // 1er acceptable
					// remise a max de bestcost pour que le 1er acceptable soit le "meilleur du voisinage"
				{
					bestmove->valuation = LONG_MAX;
					accepted = 1;
				}
				if (minneigh > 1 && currentmove->valuation < bestmove->valuation)  // mise a jour du meilleur acceptable
				{bestmove->copymove(currentmove);}
			} else
				// mise à jour du meilleur dans le voisinage (courant non acceptable et aucun acceptable n'a été trouvé)
			{
				if (!accepted && nbhsearch->finished >= 2 && currentmove->valuation < bestmove->valuation
				        && nbrefused < nbhsearch->finished)
				{bestmove->copymove(currentmove);}
				nbrefused++;
			}
		}
		nbtries++;
	}

	nhtries += nbtries;
	avgnhtries = ((nbmoves - 1) * avgnhtries + nbtries) / nbmoves;
	avgsqnhtries = ((nbmoves - 1) * avgsqnhtries + nbtries * nbtries) / nbmoves ;
	nbmoves++;
	// choix du mouvement à exécuter
	// le meilleur accepté ou le meilleur faisable (si aucun accepté)
	if ((accepted && minneigh > 1) || (!accepted && nbhsearch->finished >= 2)) {
		currentmove = bestmove;
	}
	// le premier faisable si aucun n'est accepté
	if (!accepted && nbhsearch->finished == 1 && feasiblefound) {
		currentmove = firstmove;
	}
	if ((accepted) || ((nbhsearch->finished != 0) && feasiblefound))
		// mise à jour des structures de données de l'heuristique, du problème, puis
		// exécution du mouvement.
	{
		nbhsearch->spareneighboradjust(configuration, currentmove);
		mheur->executebeforemove(currentmove, configuration, problem);
		configuration->update_conflicts(problem, currentmove);
		problem->move_execution(configuration, currentmove);
		//	  *ofile << " valeur " << currentmove->valuation << " verif : " << problem->config_evaluation(configuration)<< endl;

		if (TRACEMODE == 2) {
			Statistiques->costvalues.push_back(currentmove->valuation);
			Statistiques->examinedneighbors.push_back(nbtries);
			//    Statistiques->examinedneighbors.push_back(maxneigh);
		}
		if (bestfound)
			// stockage du nouveau meilleur
		{
			problem->best_config->copy_element(configuration);
		}
		return 1; // mouvement effectué
	}
	return 0; // pas de mouvement effectué
}


int LSAlgorithm::isfeasible(Move *move)
{  return 1;}

/* le test de faisablité d'un mouvement dans une marche de GWW: rester sous le seuil */
int LSAlgorithmGWW::isfeasible(Move *move)
{return (move->valuation <= threshold);}


/* Une marche aléatoire : on initialise la métaheuristique et on effectue "walklength" mouvements */
void LSAlgorithm::randomwalk
(OpProblem *problem, Configuration *configuration)
{
	int res;
	mheur->reinit(problem);
	for (int i = 0; i < walklength; i++) {
		Long aval = configuration->valuation;
		res = configurationmove(problem, configuration);
		if (res) {
			Statistiques->nb_moves[Statistiques->current_try] ++;
			if (configuration->valuation < aval)
				Statistiques->nb_moves_down[Statistiques->current_try] ++;
			else if (configuration->valuation > aval)
				Statistiques->nb_moves_up[Statistiques->current_try] ++;
		}



		if (configuration->valuation == problem->lower_bound) break; // critere arret
	}
}

/* un algo simple de recherche locale n'a qu'un élément dans la population  */
void LSAlgorithm::run(OpProblem *problem, Configuration **population)

{
	nhtries = 0;
	nbmoves = 1;
	nbhsearch->initsearch();

	randomwalk(problem, population[0]);
//  ecriture_fin_lsrun(avgnhtries,avgsqnhtries);
}

//************************************** GWW ******************************************************

/* Marche aléatoire (ou recherche locale) de gww
   le dernier mouvement est avec un seuil baissé pour forcer une amélioration et éviter un regroupement  si lastmovedescent=1
   En cas de stagnation, on arrête la marche si nomovestop=1
*/

void GWWAlgorithm::randomwalk(OpProblem *problem, Configuration *configuration)
{
	int val;
	if (configuration->regrouped)
		for (int i = 0; i < walkalgorithm->walklength - 1; i++) {
			val = walkalgorithm->configurationmove(problem, configuration) ;
			if (val) Statistiques->nb_moves[Statistiques->current_try] ++;
			if (nomovestop) // stagnation : arret de la marche de la particule
				if (val == 0) return;
			if (configuration->valuation == problem->lower_bound) return;   // on a trouve une solution : arret
		}
	Long thresh1 =  walkalgorithm->threshold;
	if (lastmovedescent) thresholdupdate();  // baisse lors du dernier mouvement (on suppose la meme baisse que la
	//dernière effectuée)
	val = walkalgorithm->configurationmove(problem, configuration) ;
	if (val) Statistiques->nb_moves[Statistiques->current_try] ++;
	if (lastmovedescent)  walkalgorithm->threshold = thresh1; // remise du seuil
	if (configuration->valuation == problem->lower_bound) return;   // on a trouve une solution : arret
}


/* nombre de particules au seuil avant la baisse du seuil (pour les statistiques) : la population etant triée à l'appel */

int GWWAlgorithm::nb_threshold_population(Configuration **population) {return 0;}

int ThresholdGWWAlgorithm::nb_threshold_population(Configuration **population)
{
	int nb_threshold = 0;
	int i = 0;
	while (i < populationsize && population[i]->valuation == walkalgorithm->threshold)
	{nb_threshold++; i++;}
	return nb_threshold;
}

/* l'algo GWW : boucle sur regroupement et marche aleatoire */
void GWWAlgorithm::run(OpProblem *problem, Configuration **population)
{
	int nb_pas = nbiteration;
	Long  meilleur = valeur_min(population, populationsize);
// initialisation du champ regrouped (=0 si marche depend du regroupement, =1 sinon)
	for (int i = 0; i < populationsize; i++) {
		if (regrouptest)
			population[i]->regrouped = 0;
		else
			population[i]->regrouped = 1;
	}

	thresholdcomputedelta(population);
	thresholdupdate();
	*ofile << " premier seuil " << walkalgorithm->threshold;
	if (walkalgorithm->threshold < meilleur)
		walkalgorithm->threshold = valeur_max(population, populationsize);

	while (meilleur != problem->lower_bound && walkalgorithm->threshold >= meilleur && nb_pas > 0) {
		regrouping(population);
		total_nhtries = 0;
		total_nbmoves = 0;
		populationrandomwalk(problem, population);
		populationsort(population, populationsize);
		if (elitism) populationkeepbest(problem, population);
		meilleur = valeur_min(population, populationsize);
		thresholdcomputedelta(population);
		if (TRACEMODE)
//       ecriture_changement_seuil
			*ofile << walkalgorithm->threshold << thresholddelta <<  meilleur
			       << population[0]->valuation << population[populationsize / 2]->valuation
			       << total_nhtries << nb_threshold_population(population);
		thresholdupdate();
		thresholdchangesupdate();
		nb_pas--;
	}
	Statistiques->thresholdchanges += thresholdchanges;
// ecriture_fin_gww(thresholdchanges, total_nbmoves);

}



/* ELITISME : si le meilleur absolu est meilleur que le meilleur courant, on remplace
   le plus mauvais de la population par le meilleur absolu */

void GWWAlgorithm::populationkeepbest(OpProblem *problem, Configuration **population)
{
	if (population[populationsize - 1]->valuation > problem->best_config->valuation) {
		population[0]->copy_element(problem->best_config);
// *ofile << "best config " << problem->best_config.valuation << " evaluation " << problem->config_evaluation(problem->best_config) << endl;
		populationsort(population, populationsize);

	}
// *ofile << " meilleur " << population[populationsize-1].valuation;
}




// les méthodes de baisse du seuil

void GWWAlgorithm::thresholdcomputedelta(Configuration **population)
{;}


// baisse du seuil avec un taux thresholddescent (valeur typique 0.005) et une valeur thresholdmin


void StandardGWWAlgorithm::thresholdcomputedelta(Configuration **population)
{
	thresholddelta = (Long)(thresholddescent * (walkalgorithm->threshold - thresholdmin));
	if (thresholddelta <= 0) thresholddelta = 1;
}


/* baisse rapide du seuil avec un taux thresholddescent (valeur typique 0.005)
    à partir de la plus mauvaise valeur et une valeur thresholdmin */
void FastStandardGWWAlgorithm::thresholdcomputedelta(Configuration **population)
{
	thresholddelta = walkalgorithm-> threshold - population[0]->valuation +
	                 (Long)(thresholddescent * (population[0]->valuation - thresholdmin));
	if (thresholddelta <= 0) thresholddelta = 1;
}

// baisse pour que nbkilled soient tues
void AdaptiveGWWAlgorithm::thresholdcomputedelta(Configuration **population)
{
	Long minval = population[nbkilled - 1]->valuation;
	if (minval < walkalgorithm->threshold)  thresholddelta = walkalgorithm->threshold - minval;
	else thresholddelta = 1;
}

/* baisse jusqu'au minimum des valeurs du seuil obtenu par les 2 modes precedents */
void FastAdaptGWWAlgorithm::thresholdcomputedelta(Configuration **population)
{
	if (population[nbmaxkilled - 1]->valuation == walkalgorithm->threshold) thresholddelta = 1;
	else {
		Long descent1 =  population[nbkilled - 1]->valuation;
		Long descent2 = (Long)((1 - thresholddescent) * walkalgorithm->threshold);
		Long descent = min(descent1, descent2);

		thresholddelta = walkalgorithm->threshold - descent;
	}

}

/* baisse déterminée par un taux de la distance du pire au médian */
void MedianAdaptGWWAlgorithm::thresholdcomputedelta(Configuration **population)
{
	thresholddelta = 1 + walkalgorithm->threshold - valeur_max(population, populationsize) +
	                 (Long)(mediandescent * (valeur_max(population, populationsize) - valeur_mediane(population, populationsize)));
}

/* baisse déterminée par un taux de la distance du pire au meilleur */
void BestAdaptGWWAlgorithm::thresholdcomputedelta(Configuration **population)
{
	thresholddelta = 1 + walkalgorithm->threshold - valeur_max(population, populationsize) +
	                 (Long)(bestdescent * (valeur_max(population, populationsize) - valeur_min(population, populationsize)));
}

void GWWAlgorithm::thresholdupdate()
{;}


void ThresholdGWWAlgorithm::thresholdupdate()
{ walkalgorithm->threshold = walkalgorithm->threshold - thresholddelta;}


void LSAlgorithm::initthreshold(Configuration **population, int taille)
{threshold = LONG_MAX;}

void GWWAlgorithm::initthreshold(Configuration **population, int taille)
{
	walkalgorithm->threshold = LONG_MAX;
}

void ThresholdGWWAlgorithm::initthreshold(Configuration **population, int taille)
{
	walkalgorithm->threshold = valeur_max(population, taille);
	thresholdchanges = 0;
}


void GWWAlgorithm::thresholdchangesupdate()
{;}

void ThresholdGWWAlgorithm::thresholdchangesupdate()
{thresholdchanges++;}



/* Marche (aléatoire ou recherche locale) sur l'ensemble de la population : appel de la marche aléatoire sur chaque élément */
void GWWAlgorithm::populationrandomwalk
(OpProblem *problem,  Configuration **population)

{
	walkalgorithm->nhtries = 0;
	walkalgorithm->nbmoves = 0;
	for (int i = 0; i < populationsize; i++) {
		randomwalk(problem, population[i]);
		if (population[i]->valuation == problem->lower_bound) break; // arret : borne inferieure atteinte
	}
	total_nhtries += walkalgorithm->nhtries;
	total_nbmoves += walkalgorithm->nbmoves;
}


/* les méthodes d'interaction entre elements (travail sur la population) */

void GWWAlgorithm::regrouping(Configuration **population)
{;}

/* Elimination des configurations superieures ou egales au seuil */

void StandardGWWAlgorithm::regrouping(Configuration **population)
{
	vector<int> bonnes_config;
	vector<int> mauvaises_config;
	for (int i = 0; i < populationsize; i++) {
		if (population[i]->valuation <= walkalgorithm->threshold)
			bonnes_config.push_back(i);
		else
			mauvaises_config.push_back(i);
	}
	for (vector<int>::iterator ci = mauvaises_config.begin(); ci != mauvaises_config.end(); ci++) {
		int bon = (int)(drand48() * bonnes_config.size());
		population[*ci]->copy_element(population[bonnes_config[bon]]);
		population[*ci]->regrouped = 1;
		population[bonnes_config[bon]]->regrouped = 1;
	}
}



/* Elimination des plus mauvais au dessus du seuil  (la population est triee) */

void AdaptiveGWWAlgorithm::regrouping(Configuration **population)
{
	int nb_to_kill = 0;
	for (int i = 0; i < populationsize; i++) {
		if (population[i]->valuation <= walkalgorithm->threshold) break;
		nb_to_kill ++;
	}
	if (TRACEMODE) ecriture_nb_tues(nb_to_kill);
	for (int i = 0; i < nb_to_kill; i++) {
		int bon = ((int)(drand48() * (populationsize - nb_to_kill))) + nb_to_kill  ;
		population[i]->copy_element(population[bon]);
		population[i]->regrouped = 1;
		population[bon]->regrouped = 1;
	}
}


/* Elimination des nbkilled  plus mauvais */


void NothresholdGWWAlgorithm::regrouping(Configuration **population)

{
	for (int j = 0; j < nbkilled; j++) {
		Long val = population[0]->valuation;
		int mauvais = 0;
		for (int i = 0; i < populationsize; i++)
			if (population[i]->valuation > val)
			{mauvais = i; val = population[i]->valuation;}
		{
			int bon = (int)(drand48() * (populationsize - 1));
			if (bon >= mauvais) bon++;
			population[mauvais]->copy_element(population[bon]);
			population[mauvais]->regrouped = 1;
			population[mauvais]->regrouped = 1;
		}
	}
}



//*****************************************************************************************************


// les méthodes particulières des algorithmes de recherche locale


void Metaheuristic::reinit(OpProblem *problem)
{
	;
}

void Metaheuristic::executebeforemove(Move *move, Configuration *configuration, OpProblem *problem)
{;}

int Metaheuristic::acceptance(Move *move, Configuration *config)
{return 1;}

/* la marche aléatoire : tout mouvement est accepté */
int RandomSearch::acceptance(Move *move, Configuration  *config)
{return 1;}

/* Descente avec plateau : tout mouvement améliorant ou neutre est accepté */
int GreedySearch::acceptance(Move *move, Configuration *config)
{return (move->valuation <= config->valuation);}



/* Metropolis à temperature constante : condition d'acceptation d'un mouvement */

int Metropolis::acceptance(Move *move, Configuration *config)
{
	if ((move->valuation <= config->valuation)
	        ||
	        (drand48() < exp(-(move->valuation  - config->valuation) / temperature)))
		//     (drand48() < exp(-(move->valuation  - config->valuation).to_double()/temperature)))
		return 1;
	else return 0;
}

void Metropolis::adjustparameter(int temp)  // pour le reglage automatique qui ne marche qu'avec des paramètres entiers
{temperature = ((double)(temp)) / 100 ;}



/* la methode de recherche avec liste taboue */

/* remise a vide de la liste taboue */
void TabuSearch::reinit(OpProblem *problem)
{
	for (list<Move *>::iterator imove = move_list.begin(); imove != move_list.end(); imove++)
		delete *imove;
	move_list.clear();
}

/* remise a zero des dates de fin de temps tabou par mouvement*/
void IncrTabuSearch::reinit(OpProblem *problem)
{
	currentproblem = problem;
	tabutime.clear();
	for (int i = 0; i < problem->nbtabuindex(); i++)
		tabutime.push_back(0);
	nbiter = 0;
}

/* acceptation d'un mouvement : non tabou  (le critère d'aspiration est dans l'algo de recherche du voisin) */
int TabuSearch::acceptance(Move *move,  Configuration *config)
{
	return (nontabumove(move));
}


void TabuSearch::adjustparameter(int length)  // pour le reglage automatique de la longueur de la liste taboue
{tabulength = length;}


/* acceptation d'un mouvement : non tabou  et améliorant (le critère d'aspiration est dans l'algo de recherche du voisin) */
int TabuGreedySearch::acceptance(Move *move,  Configuration *config)
{
	return (move->valuation <= config->valuation && nontabumove(move));
}

int IncrTabuGreedySearch::acceptance(Move *move, Configuration *config)
{
	return (move->valuation <= config->valuation && IncrTabuSearch::acceptance(move, config));
}

/* la liste taboue est implantée comme une liste de pointeurs sur mouvement : test de présence avec eqmove*/
int TabuSearch::nontabumove(Move *move)
{
	for (list<Move *>::iterator imove = move_list.begin(); imove != move_list.end(); imove++)
		if ((*imove)->eqmove(move))   return 0;
	return 1;
}

/* la liste est traitée comme une file de longueur maximale tabulength : fonction de mise à jour avant mouvement */
void TabuSearch::executebeforemove(Move *move, Configuration *configuration, OpProblem *problem)
{
	Move *move1 = move->computetabumove(configuration);
	move_list.push_back(move1);
	if ((int) move_list.size() > tabulength) {
		Move *move2 = move_list.front();
		delete move2;
		move_list.pop_front();
	}
}


/* la liste tabou est implantée comme une date de fin d'état tabou pour un mouvement */
int IncrTabuSearch::acceptance(Move *move, Configuration *config)
{
	return (tabutime[currentproblem->tabuindex(move, config)] <= nbiter)
	       ;
}


void IncrTabuSearch::executebeforemove(Move *move, Configuration *config, OpProblem *problem)
{
	tabutime[currentproblem->tabuinverseindex(move, config)] = nbiter + tabulength + 1;
	nbiter++;
}



/* l'acceptation à seuil (TA) : méthode d'acceptation d'un mouvement*/

int ThresholdAccepting::acceptance(Move *move, Configuration *config)
{
	return ((move->valuation - config->valuation) <  thresholdaccept);
}



void ThresholdAccepting::executebeforemove(Move *move, Configuration *configuration, OpProblem *problem)
{thresholdaccept -= delta;}

/* remise du seuil à thresholdinit */
void ThresholdAccepting::reinit(OpProblem *problem)
{thresholdaccept = thresholdinit;}


/* le recuit simulé avec baisse de température constante à chaque mouvement*/
/* le recuit simulé : acceptation en fonction de la température */

int SimulatedAnnealing::acceptance(Move *move, Configuration *config)
{
	if ((move->valuation <= config->valuation)
	        ||
	        (drand48() < exp(-(move->valuation  - config->valuation) / temperature)))
		//     (drand48() < exp(-(move->valuation  - config->valuation).to_double()/temperature)))
		return 1;
	else return 0;
}

/* baisse de température constante de delta à chaque mouvement */
void SimulatedAnnealing::executebeforemove(Move *move, Configuration *configuration, OpProblem *problem)
{temperature -= delta;}

/* remise de la temperature initiale */
void SimulatedAnnealing::reinit(OpProblem *problem)
{temperature = inittemperature;}

/* les mouvements */

Move::Move() {valuation = 0;}
CSPMove::CSPMove() {valuation = 0; value = 0; variable = 0;}


/* copie de mouvement move1 dans this */
void Move::copymove(Move *move1)
{valuation = move1->valuation;}

void CSPMove::copymove(Move *move1)
{
	valuation = move1->valuation;
	value = ((CSPMove *)move1)->value;
	variable = ((CSPMove *)move1)->variable;
}


/* le mouvement stocké tabou est le mouvement inverse du mouvement effectué */
Move *CSPMove::computetabumove(Configuration *configuration)
{
	CSPMove *tabumove = new CSPMove();
	tabumove->variable = variable;
	tabumove->value = configuration->config[tabumove->variable];
	return tabumove;
}


/* test d'égalité de 2 mouvements : utile pour la liste taboue */
int Move::eqmove(Move *move1)
{return 0;}

/* Test que le mvt (variable,valeur) n'est pas tabou. */
int CSPMove::eqmove(Move *move1)
{
	return (
	           ((CSPMove *)move1)->value == value &&
	           ((CSPMove *)move1)->variable == variable);
}





/*-----------------------------------------------------------------------------------------*/
/*                             TabuAcceptingrate                                           */
/*-----------------------------------------------------------------------------------------*/
// implanté comme une sous-classe de TabuSearch

// le constructeur :

TabuAcceptingrate::TabuAcceptingrate(int tabul, float cPd, float cP0)
	:  TabuSearch(tabul) {  P0 = cP0; Pd = cPd;}



/* critère d'acceptation : non tabou et pourcentages d'acceptation suivant sens du mouvement (détériorant, neutre, améliorant) */
int TabuAcceptingrate::acceptance(Move *move, Configuration *config)
{
	return
	    (nontabumove(move)
	     &&
	     ((move->valuation < config->valuation) ||
	      ((move->valuation == config->valuation) &&
	       (drand48() < P0)) ||
	      ((drand48() < Pd))
	     )
	    );
}

