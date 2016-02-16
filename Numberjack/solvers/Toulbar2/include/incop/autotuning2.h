/* Reglage automatique d'un algorithme à un paramètre */
/** Automatic tuning of a local search algorithm with one parameter */

class Tuning 
{public :
 int param;
 int parammin;
 int parammax;
 int parambest;
 int paramdico1;
 int paramdico2;
 int boundmin;
 int boundmax;
 double value;
 double valuemin;
 double valuemax;
 double valuebest;
 double valuedico1;
 double valuedico2;
 int seed;
 int nbtries;
 int maxtuning;
 int trynumber;
 int stop;
 float referencetime;
 virtual int firsttry();
 Tuning (int pinit, int seed1, int nbessai, int maxtun, int bmin, int bmax);
 virtual void run (LSAlgorithm* algo, OpProblem* problem, Configuration** configuration);
 virtual void onerun (OpProblem* problem, LSAlgorithm* algo, Configuration** population);
 virtual void onerunparam (OpProblem* problem, LSAlgorithm* algo, Configuration** population, int param);
 void oneparameterwrite(LSAlgorithm* algo, int parameter);
 virtual int end();
};



void autosolving (LSAlgorithm* algo, Configuration** configuration, OpProblem* problem, int npb, int graine1, int nbessais, double maxtime, int initwalklength);
int autotuning(LSAlgorithm* algo, Configuration** population, OpProblem* problem, int graine1, int nbessais);
int endoslve(double maxtime);


/* Réglage automatique d'un algorithme à deux paramètres */
/** Automatic tuning of a local search algorithm with two parameters */
class DoubleTuning : public Tuning
{public :
 int param2min;
 int param2max;
 int param2best;
 int param2dico1;
 int param2dico2;
 int bound2min;
 int bound2max;
 double value2;
 double value2min;
 double value2max;
 double value2best;
 double value2dico1;
 double value2dico2;
 int trynumber1;
 int maxtuning2;
 int firsttry();
 DoubleTuning (int pinit1, int pinit2, int seed1, int nbessai, int maxtun, int maxtun2, int bmin, int b2min,int bmax, int b2max);
 void run (LSAlgorithm* algo, OpProblem* problem, Configuration** configuration);
 void onerun2param (OpProblem* problem, LSAlgorithm* algo, Configuration** configuration, int parameter1, int parameter2);
 void simplerun (LSAlgorithm* algo, OpProblem* problem, Configuration** configuration, int parameter);
 void doubleparameterwrite(LSAlgorithm* algo, int parameter1, int parameter2);
};


void autosolving2 (LSAlgorithm* algo, Configuration** configuration, OpProblem* problem, int npb, int graine1, int nbessais, double maxtime, int initwalklength);
void autosolving1 (LSAlgorithm* algo, Configuration** configuration, OpProblem* problem, int npb, int graine1, int nbessais, double maxtime, int initwalklength);

void autotuning2(LSAlgorithm* algo, Configuration** population, OpProblem* problem, int graine1, int nbessais, int&
parameter1, int& parameter2 ,int tuningwalkrate, int tuningmaxtries);
int endoslve(double maxtime);
