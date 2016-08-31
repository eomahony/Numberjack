

//---------------------------------------------------------------------------------------------
//    trombe : les stats sous forme de classe (... enfin une struct quoi !)

class Stat_GWW
{
public:
	int nb_pbs;      // Nb de problèmes différents essayés (pour les CSP aléatoires)
	int *trouve;    // trouve[i] (i ds [0,nb_pbs[) : contient le nombre de fois que
	//           le pb i est résolu (entre 0 et max_tries)
	int stop_trouve;
	int current_pb;  // Le numéro de problème courant (ds [0,nb_pbs[)
	//
	int max_tries;   // Le nombre d'essais par problème
	Long   *cost_try; // cost_try [j] (j ds [0,max_tries[): le meilleur cout trouvé à l'essai j
	float *execution_time_try;  // execution_time_try [j] (j ds [0,max_tries[): temps d'exécution
	//                                              de l'essai j
	float cost_meanvalue; // la moyenne des meilleurs couts des essais
	int current_try; // Le numéro d'essai courant (ds [0,max_tries[) pour le problème 'current_pb'
	int *nb_moves;
	int *nb_moves_up;
	int *nb_moves_down;

	double nb_moves_avg;
	float *total_problem_time;
	vector<Long> costvalues;
	vector<int> examinedneighbors;
	float total_execution_time; // temps d'execution total depuis le début
	float average_execution_time; // temps d'execution total depuis le début
	// ofstream* stat_file;           // le fichier où s'affiche les stat
	int thresholdchanges;   // nombre de changements de seuil pour l'ensemble des essais
	Stat_GWW(int number_pbs, int number_tries);
	void init_pb(int t);
	void init_run();
	void init_try(int trynumber);
	void execution_report(int ntry, Long lower_bound);
};




void sigaction();
int argument2ul(char *arg, char *message);
double argument2d(char *arg, char *message);
double argument2bd(char *arg, char *message, double min1, double max1);
int argument2bul(char *arg, char *message, int min1, int max1);
string argument2lp(char *arg, char *message, list<string> &possibles);
void handler_stat(int sig);
IncompleteAlgorithm *algo_creation(char **argv, int &narg, int &taille, int &graine1, int &nbessais);
void executer_essai
(OpProblem *problem, IncompleteAlgorithm *algo, Configuration **population, int taille, int graine1, int nessai, vector<int> *initconfig = NULL);

void ecriture_stat_probleme();
void ecriture_statistiques_global();
void arguments_tracemode(char **argv, int &narg);
void arguments_tempscpu(char **argv, int &narg, double &maxtime);
void arguments_arret(char **argv, int &narg, int &stop);
void arguments_borneinf(char **argv, int &narg, Long &borneinf);

Long valeur_max(Configuration **population, int taille);
Long valeur_min(Configuration **population, int taille);
Long valeur_mediane(Configuration **population, int taille);
void populationsort(Configuration **population, int taille);

void ecriture_changement_seuil(Long seuil, Long delta, Long meilleur, Long pire, Long mediane, int nbessaisvoisins, int nb_au_seuil);
void ecriture_nb_tues(int nb_tues);
void ecriture_fin_gww(int nb_chang_seuil, int nb_mouv);
void ecriture_meilleure_valeur(string &method, Long valeur, Long seuil, int nbmouv, int maxvoisins);
void ecriture_fin_lsrun(double avgnhtries, double avgsqnhtries);
void ecriture_stat_essai();
void ofile_name(char *filename, int argc, char **argv);
