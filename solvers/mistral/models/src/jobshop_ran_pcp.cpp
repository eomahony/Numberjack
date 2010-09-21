
#include <mistral_sol.h>
#include <mistral_gen.h>

using namespace std;



JSPGenerator* jsp;
char machine[20] = {'=','-','#',':','o'
		    ,'+','*','x','.','\\'
		    ,'@','%','?','z','s'
		    ,'w','v','u',';','c'};
int makespan, minfsble, maxfsble;
int upperbound(int *sol);
void print_sol(int *solution);

void store_sol(int *solution, VarArray *tasks)
{
  int k = 0;
  for(int i=0; i<jsp->Njobs; ++i) 
    for(int j=0; j<jsp->Nmach; ++j) 
      solution[k++] = tasks[i][j].min();
}


// THE MODEL ////////////////////
void setup(CSP& model, VarArray& preced, VarArray *task_in_order, int edge) {

  int x, y, z, min, max;
  int   dur_machine[jsp->Nmach][jsp->Njobs];
  VarArray task_per_machine[jsp->Nmach];
  VarArray ordering[jsp->Nmach];
 
  for( x=0; x<jsp->Nmach; ++x ) {
    VarArray mach_x;//( jsp->Njobs );
    task_per_machine[x] = mach_x;
  }

  for( y=0; y<jsp->Njobs; ++y ) {
    min=0;
    max=makespan;
    for( x=0; x<jsp->Nmach; ++x )
      max-=jsp->Durations[x][y];
    for( x=0; x<jsp->Nmach; ++x ) {
      task_in_order[y].add( Variable( min, max ) );
      min += jsp->Durations[x][y];
      max += jsp->Durations[x][y];
      //task_per_machine[jsp->Machines[x][y]][y] = task_in_order[y][x];
      task_per_machine[jsp->Machines[x][y]].add( task_in_order[y][x] );
      dur_machine[jsp->Machines[x][y]][y] = jsp->Durations[x][y];
    }
  }

  for( y=0; y<jsp->Njobs; ++y )
    for( x=1; x<jsp->Nmach; ++x )
      model.add( Precedence( 
			    task_in_order[y][x-1],
			    jsp->Durations[x-1][y], 
			    task_in_order[y][x]
			     ) );

  for( x=0; x<jsp->Nmach; ++x ) {
    for( y=0; y<jsp->Njobs; ++y )
      for( z=y+1; z<jsp->Njobs; ++z )
	
 	if( task_per_machine[x][y].min() + dur_machine[x][y] > 
 	    task_per_machine[x][z].max() ) {
	  model.add( Precedence(task_per_machine[x][z], 
				dur_machine[x][z], 
				task_per_machine[x][y]) );
	}
	else if( task_per_machine[x][z].min() + dur_machine[x][z] > 
		 task_per_machine[x][y].max() ) {
	  model.add( Precedence(task_per_machine[x][y], 
				dur_machine[x][y], 
				task_per_machine[x][z]) );
	}
	else { 
	  ordering[x].add( Disjunctive( 
				      task_per_machine[x][y],  
				      dur_machine[x][y],
				      task_per_machine[x][z],  
				      dur_machine[x][z]
				       ) );
	}
    preced.add( ordering[x] );
  }

  model.add( preced );

  //  if(edge)
  //  for( x=0; x<jsp->Nmach; ++x ) 
  //    model.add( EdgeFinder(task_per_machine[x], dur_machine[x]) );
}

int main(int argc, char** argv) {

  jsp = new JSPGenerator( (argc > 1 ? atoi(argv[1]) : 8 ),
			  (argc > 2 ? atoi(argv[2]) : 7 ) );
  int edge   = (argc > 3 ? atoi(argv[3]) : 0);
  int nbinst = (argc > 4 ? atoi(argv[4]) : 1);
  long seed  = (argc > 5 ? (long)(-1 * atoi(argv[5])) : 11041979);

  static Random_ generator;
  
  int solution[jsp->Njobs*jsp->Nmach]; 

  for(int k=0; k<nbinst; k++) {
    int act;
    int timeseed = (int)(1048576 * generator.ran2(&seed)), 
      machineseed = (int)(1048576 * generator.ran2(&seed));
    
    // GENERATE THE RANDOM INSTANCE AND PRINT IT
    jsp->jspgen( "taillard", jsp->Njobs, jsp->Nmach, 
		 timeseed, machineseed, &minfsble );
    
    cout << endl
	 << "\tJob-Shop Scheduling Problem (" 
	 << jsp->Njobs << " jobs, " 
	 << jsp->Nmach << " machine)\n\n";
     
//     for(int x=0;x<jsp->Njobs;x++) {
//       for(int y=0;y<jsp->Nmach;y++) {
// 	cout<<"[";
// 	for(int i=0;i<(jsp->Durations[y][x]-2);i++)
// 	  cout<<machine[jsp->Machines[y][x]];
// 	cout<<"] ";
//       }
//       cout<<endl;
//     }	 
//     cout<<endl;

    // COMPUTE AN UPPER BOUND WITH A GREEDY ALGORITHM 
    maxfsble = upperbound(solution);

    // maxfsble is the greatest makespan for which the instance
    // is proved to be SAT
    // minfsble is the lowest makespan for which the instance
    // is not proved to be UNSAT
    // DICHOTOMIC SEARCH
    bool solved = true;
    while( minfsble < maxfsble ) {
      makespan = (int)(floor(((double)minfsble + (double)maxfsble)/2));
      
      cout << "\t" << setw(5) << minfsble << " " << setw(5) << maxfsble 
	   << " (" << setw(5) << makespan << "):\t";
      cout.flush();
                      
      VarArray disjuncts;
      VarArray tasks[jsp->Njobs];
      CSP model;
      setup(model, disjuncts, tasks, edge);

      //Pcp heuristic;
      DomOverWDeg heuristic;
      //Impact heuristic;
      Solver s(model, disjuncts, heuristic);
      s.setTimeLimit(5);
      
      //int result = s.solve();
      //s.setVerbosity(1);
      s.setRandomized();

      //int sacres = s.sacPreprocess();

      int result = s.solve_and_restart(GEOMETRIC, 64, 1.2);


      if( result == SAT ) {
	store_sol( solution, tasks );
       	maxfsble = makespan;
      } else {
	minfsble = makespan+1;
	if( result == LIMITOUT )
	  solved = false;
      }

      s.printStatistics(std::cout);
      cout << endl;


      //exit(0);
    }

    if( !solved ) {
      VarArray disjuncts;
      VarArray tasks[jsp->Njobs];
      VarArray endlasttask( jsp->Njobs, 0, maxfsble );

      CSP model;
      setup(model, disjuncts, tasks, edge);
      for(int i=0; i<jsp->Njobs; ++i)
	model.add( tasks[i][jsp->Nmach-1] + jsp->Durations[jsp->Nmach-1][i] == endlasttask[i] );

      model.add( Minimise(Max(endlasttask)) );

      //Pcp heuristic;
      DomOverWDeg heuristic;
      Solver s(model, disjuncts, heuristic);
      s.solve();
      maxfsble = s.goal->upper_bound;
    }
    
    cout << endl << "\tbest makespan found: " << maxfsble << endl << endl;
    //print_sol(solution);
    cout << endl;
  }
  delete jsp;
}

void print_sol(int *solution) {
  int x,y,z;
  for(x = 0; x<jsp->Njobs; x++) {
    for(y = 0; y<jsp->Nmach; y++) {
      if(y == 0)
	cout<<setw(solution[x*jsp->Nmach+y]+1)<<"[";
      else 
	cout<<"]"
	    <<setw(solution[x*jsp->Nmach+y]-jsp->Durations[y-1][x]-solution[x*jsp->Nmach+y-1]+1)
	    <<"[";
      for(z=0; z<jsp->Durations[y][x]-2; z++)
	cout<<machine[jsp->Machines[y][x]];
    }
    cout <<"]"<< endl;
  }
  //delete jsp;
}

int upperbound(int *sol) {
  int currentTasks[jsp->Njobs], 
    currentDeadlinesMach[jsp->Nmach],
    currentDeadlinesJobs[jsp->Njobs],
    mkp=0, mstrt, nmkp, bi, bmkp;
  std::fill(currentTasks, currentTasks+jsp->Njobs, 0);
  std::fill(currentDeadlinesMach, currentDeadlinesMach+jsp->Nmach, 0);
  std::fill(currentDeadlinesJobs, currentDeadlinesJobs+jsp->Njobs, 0);
  
  //We assign the tasks that will increase the least the current makespan
  while(true) {
    bmkp = NOVAL; 
    for(int i=0; i<jsp->Njobs; i++) {
      if(currentTasks[i] < jsp->Nmach) {
	mstrt = currentDeadlinesMach[jsp->Machines[currentTasks[i]][i]];
	if(mstrt < currentDeadlinesJobs[i])
	  mstrt = currentDeadlinesJobs[i];
	nmkp = mstrt + jsp->Durations[currentTasks[i]][i];
	if(bmkp > nmkp) {
	  bmkp = nmkp;
	  bi = i;
	}
      }
    }
    if(bmkp == NOVAL)
      break;
    sol[bi * jsp->Nmach + currentTasks[bi]] = bmkp - jsp->Durations[currentTasks[bi]][bi];
    mkp = (mkp > bmkp ? mkp : bmkp);
    currentDeadlinesMach[jsp->Machines[currentTasks[bi]][bi]] = bmkp;
    currentDeadlinesJobs[bi] = bmkp;
    currentTasks[bi]++;
  }
  return mkp;
}




  


  




