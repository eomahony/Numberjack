
#include <mistral_sol.h>

/* 
 * a model of the bibd problem
 */

 
using namespace std;

int main(int argc, char *argv[])
{  
  int    v      = ( argc > 1 ? atoi(argv[1]) : 7     );
  int    b      = ( argc > 2 ? atoi(argv[2]) : 70   );
  int    r      = ( argc > 3 ? atoi(argv[3]) : 30   );
  int    k      = ( argc > 4 ? atoi(argv[4]) : 3     );
  int    l      = ( argc > 5 ? atoi(argv[5]) : 10    );

  
  // 7 140 60 3 20
  // 7 210 90 3 30
  // 7 280 120 3 40
  // 7 315 135 3 45
  // 7 350 150 3 50
  

  cout //<< endl
       << "//\tBIBD" 
       << setw(4) << v 
       << setw(4) << b
       << setw(4) << r 
       << setw(4) << k
       << setw(4) << l
       << " " << endl;

  // create the csp
  CSP model;
  int i, j, x, y;

  VarArray columns[v];
  for(i=0; i<v; ++i) {
    VarArray col_i( b, 2 ); 

    //cout << col_i << endl;

    columns[i] = col_i;
    model.add( Sum( columns[i] ) == r );
    if(i>0) model.add( columns[i] >= columns[i-1] );
    //if(i==1) model.add( columns[i] >= columns[i-1] );
  }

  VarArray rows[b];
  for(j=0; j<b; ++j) {
    VarArray row_j( v );
    rows[j] = row_j;
    for(i=0; i<v; ++i) {
      rows[j][i] = columns[i][j];
    }

    //cout << row_j << endl;

    model.add( Sum( rows[j] ) == k );
    if(j>0) model.add( rows[j] >= rows[j-1] );
  }

  VarArray scalar[ (v*(v-1)/2) ];
  for(j=0; j<(v*(v-1)/2); ++j) {
    VarArray scal_j( b );
    scalar[j] = scal_j;
  }
  y = 0;
  for(i=0; i<v-1; ++i)
    for(j=i+1; j<v; ++j) {
      for(x=0; x<b; ++x)
	scalar[y][x] = (columns[i][x] && columns[j][x]) ;

      //cout << scalar[y] << endl;

      model.add( Sum( scalar[y] ) == l );
      ++y;
    }

  VarArray searchVars;
  for(i=0; i<v; ++i)
    searchVars.add( columns[i] );

  //DomOverWDeg h;
  //NoOrder h;
  Lexicographic h;
  Solver s( model, searchVars, h );


//   cout << endl;
  //  s.print(cout);
//   cout << endl;

  //  s.setTimeLimit( 60 );
  s.solve();



  s.printStatistics( cout );
  cout << endl;

  for(i=0; i<v; ++i) {
    cout << "//\t";
    for(j=0; j<b; ++j) 
      cout << ( rows[j][i].value() ? "X" : " " ) ;
    //cout << rows[i][j].value() ;
    cout << endl;
  }
  cout << endl ;
}







//	BIBD   7   7   3   3   1       SAT        0 BTS       17 NDS        3 FAILS        0 BTS/s           930 CKS       0 s
//	XXX    
//	X  XX  
//	 X X X 
//	X    XX
//	 X  X X
//	  XX  X
//	  X XX 

//	BIBD   6  10   5   3   2       SAT        0 BTS       16 NDS        6 FAILS        0 BTS/s          1095 CKS       0 s
//	XXXXX     
//	XX   XXX  
//	  XX XX X 
//	X X    XXX
//	 X  XX  XX
//	   XX XX X

//	BIBD   7  14   6   3   2       SAT        0 BTS       29 NDS       12 FAILS        0 BTS/s          2289 CKS       0 s
//	XXXXXX        
//	XX    XXXX    
//	  XX  XX  XX  
//	X X     X X XX
//	   XX   XX XX 
//	 X   XX    XXX
//	    XX X XX  X

//	BIBD   9  12   4   3   1       SAT        2 BTS       35 NDS       10 FAILS      200 BTS/s          2877 CKS       0 s
//	XXXX        
//	X   XXX     
//	 X  X  XX   
//	  X  X X X  
//	 X    X  XX 
//	X      X  XX
//	  X   X X  X
//	   XX    X X
//	   X X  X X 

//	BIBD   6  20  10   3   4       SAT        2 BTS       27 NDS       12 FAILS      200 BTS/s          2672 CKS       0 s
//	XXXXXXXXXX          
//	XXXX      XXXXXX    
//	    XXXX  XXXX  XX  
//	X   X   XXXX  X X XX
//	 XX  XX       XXXXXX
//	   X   XXX  XX X XXX

//	BIBD   7  21   9   3   3       SAT        0 BTS       40 NDS       14 FAILS        0 BTS/s          3510 CKS       0 s
//	XXXXXXXXX            
//	XXX      XXXXXX      
//	   XXX   XXX   XXX   
//	   X  XX    XXXXX X  
//	X   X   X   XX X X XX
//	 X   X  XX    X X XXX
//	  X   XX  XX     XXXX

//	BIBD   6  30  15   3   6       SAT        1 BTS       35 NDS       15 FAILS      100 BTS/s          4456 CKS       0 s
//	XXXXXXXXXXXXXXX               
//	XXXXXX         XXXXXXXXX      
//	      XXXXXX   XXXXXX   XXX   
//	X     XX    XXXXXX   XX X  XXX
//	 XX     X   XXX   XXX  X XXXXX
//	   XXX   XXX         XXXXXXXXX

//	BIBD   7  28  12   3   4       SAT        4 BTS       58 NDS       27 FAILS      400 BTS/s          5468 CKS       0 s
//	XXXXXXXXXXXX                
//	XXXX        XXXXXXXX        
//	    XXXX    XXXX    XXXX    
//	    X   XXX     XXXXXXX X   
//	X    X  X  XX   XX  X  X XXX
//	 X    XX X   X    XX   XXXXX
//	  XX      XX  XX     XX XXXX

//	BIBD   9  24   8   3   2       SAT        4 BTS       70 NDS       25 FAILS      400 BTS/s          6412 CKS       0 s
//	XXXXXXXX                
//	XX      XXXXXX          
//	  XX    XX    XXXX      
//	  X X     XX  X   XXX   
//	   X X    X X  X  X  XX 
//	X    X     X    XX X X X
//	    X X X    X X    XX X
//	 X     X X      X X X XX
//	      XX    XXX  X X  X 

//	BIBD   6  40  20   3   8       SAT       15 BTS       54 NDS       35 FAILS     1500 BTS/s          8232 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXX                    
//	XXXXXXXX            XXXXXXXXXXXX        
//	        XXXXXXXX    XXXXXXXX    XXXX    
//	X       XXXX    XXX XXX     XXXXX   XXXX
//	 XXXX       XXX    X   X    XXX XXXXXXXX
//	     XXX       XXXXX    XXXX   X XXXXXXX

//	BIBD   7  35  15   3   5       SAT       27 BTS       82 NDS       54 FAILS     2700 BTS/s          9177 CKS       0 s
//	XXXXXXXXXXXXXXX                    
//	XXXXX          XXXXXXXXXX          
//	     XXXXX     XXXXX     XXXXX     
//	X    X    XXX       XXXX XXXX XX   
//	 X    X   X  XXX    XX  XXX  X  XXX
//	  X    X   XXX  XXX     X  X  XXXXX
//	   XX   XX    X    X  XX    XXXXXXX

//	BIBD   7  42  18   3   6       SAT       18 BTS       81 NDS       53 FAILS     1800 BTS/s         12088 CKS       0 s
//	XXXXXXXXXXXXXXXXXX                        
//	XXXXXX            XXXXXXXXXXXX            
//	      XXXXXX      XXXXXX      XXXXXX      
//	      X     XXXXX       XXXXXXXXXXX X     
//	X      X    XXX  XXXX   XX    X    X XXXXX
//	 X      XXXX     X   X    XXXX X    XXXXXX
//	  XXXX         XX     XX        XXXXXXXXXX

//	BIBD  10  30   9   3   2       SAT       48 BTS      128 NDS       66 FAILS     4800 BTS/s         13466 CKS       0 s
//	XXXXXXXXX                     
//	XX       XXXXXXX              
//	  XX     XX     XXXXX         
//	  X X      XX   X    XXXX     
//	X  X       X    X        XXXXX
//	 X   X       X   XX  XX  XX   
//	      XX X   X   X     XX  XX 
//	    X X       XX  XX X     X X
//	       XX   X X    XX X  X  X 
//	     X  X X    X    X  XX X  X

//	BIBD   6  50  25   3  10       SAT        5 BTS       54 NDS       27 FAILS      500 BTS/s          9109 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXXXXXX                         
//	XXXXXXXXXX               XXXXXXXXXXXXXXX          
//	          XXXXXXXXXX     XXXXXXXXXX     XXXXX     
//	X         XXXXX     XXXX XXXX      XXXXXX    XXXXX
//	 XXXXX              XXXXX    XXXXX      XXXXXXXXXX
//	      XXXX     XXXXX    X         XXXXXX XXXXXXXXX

//	BIBD   9  36  12   3   3       SAT       91 BTS      184 NDS      119 FAILS     9100 BTS/s         18699 CKS       0 s
//	XXXXXXXXXXXX                        
//	XXX         XXXXXXXXX               
//	   XXX      XXX      XXXXXX         
//	   X  XX       XXX   XX    XXXX     
//	X   X   X      X  X    XX  XX  XXX  
//	 X    X X       X  X X   XX    XX XX
//	  X  X   X  X       X    X XXX   XXX
//	         XXX      XXX XX  X  XXX    
//	       X  XX XX  X      X     X XXXX

//	BIBD  13  26   6   3   1       SAT      184 BTS      281 NDS      200 FAILS     9199 BTS/s         39565 CKS    0.02 s
//	XXXXXX                    
//	X     XXXXX               
//	 X    X    XXXX           
//	  X    X   X   XXX        
//	 X     X          XXXX    
//	  X   X           X   XXX 
//	   X    X   X  X  X      X
//	X            X  X  X  X  X
//	    X    X  X    X X   X  
//	     X   X X        X   XX
//	     X  X     X X    X X  
//	   X      X   X  X  X X   
//	    X     X  X X     X  X 

//	BIBD   7  49  21   3   7       SAT        8 BTS       77 NDS       47 FAILS      800 BTS/s         11708 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXX                            
//	XXXXXXX              XXXXXXXXXXXXXX              
//	       XXXXXXX       XXXXXXX       XXXXXXX       
//	X      X      XXXXX         XXXXXX XXXXXX XX     
//	 X      X     XXX  XXXXXX   X     XXX     X XXXXX
//	  XXXX   XX        X     X   XX      XXX XXXXXXXX
//	      X    XXX   XX X     XX   XXXX     XX XXXXXX

//	BIBD   6  60  30   3  12       SAT       25 BTS       76 NDS       50 FAILS     2500 BTS/s         13417 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                              
//	XXXXXXXXXXXX                  XXXXXXXXXXXXXXXXXX            
//	            XXXXXXXXXXXX      XXXXXXXXXXXX      XXXXXX      
//	X           XXXXX       XXXXXXXXXXXX      XXXXX X     XXXXXX
//	 XXXXX           X      XXXXXX      XXXXXX     X XXXXXXXXXXX
//	      XXXXXX      XXXXXX                  XXXXXXXXXXXXXXXXXX

//	BIBD   7  56  24   3   8       SAT       20 BTS       94 NDS       65 FAILS     2000 BTS/s         15864 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXXXXX                                
//	XXXXXXXX                XXXXXXXXXXXXXXXX                
//	        XXXXXXXX        XXXXXXXX        XXXXXXXX        
//	        X       XXXXXXX         XXXXXXXXXXXXXXX X       
//	X        X      XXXXX  XXXXXX   XX      X      X XXXXXXX
//	 XXXXXX   X          X       X    X      XXXXX XXXXXXXXX
//	       X   XXXXX      XX      XX   XXXXX      X XXXXXXXX

//	BIBD   6  70  35   3  14       SAT        0 BTS       65 NDS       31 FAILS        0 BTS/s         14618 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                   
//	XXXXXXXXXXXXXX                     XXXXXXXXXXXXXXXXXXXXX              
//	              XXXXXXXXXXXXXX       XXXXXXXXXXXXXX       XXXXXXX       
//	X             XXXXXXX       XXXXXX XXXXXX        XXXXXXXX      XXXXXXX
//	 XXXXXXX             XX     XXXX  X      XXXXX   XX     XXXXXXXXXXXXXX
//	        XXXXXX         XXXXX    XXX           XXX  XXXXX XXXXXXXXXXXXX

//	BIBD   9  48  16   3   4       SAT        8 BTS      126 NDS       46 FAILS      800 BTS/s         14737 CKS       0 s
//	XXXXXXXXXXXXXXXX                                
//	XXXX            XXXXXXXXXXXX                    
//	    XXXX        XXXX        XXXXXXXX            
//	    X   XXX         XXXX    XXX     XXXXX       
//	X    X  X  X        X   XX     XXX  XX   XXXX   
//	 X          XXX      X    XX   XX XX  XXXXX     
//	      XX X  X   X     XX  X      X    X  X XXXXX
//	  X       XX   X XX        XX     X X  X  XX XXX
//	   X         XXX   X    XX   XX    X X  X   XXXX

//	BIBD   7  63  27   3   9       SAT       28 BTS      110 NDS       72 FAILS     2799 BTS/s         19294 CKS    0.01 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXX                                    
//	XXXXXXXXX                  XXXXXXXXXXXXXXXXXX                  
//	         XXXXXXXXX         XXXXXXXXX         XXXXXXXXX         
//	X        X        XXXXXXX           XXXXXXXX XXXXXXXX XX       
//	 X        X       XXXXX  XXXXXXX    XX      XXX      X  XXXXXXX
//	  XX       XXXXXXX                    XXXXXXX  X     XXXXXXXXXX
//	    XXXXX              XXXX     XXXX            XXXXX XXXXXXXXX

//	BIBD   8  56  21   3   6       SAT       72 BTS      172 NDS      111 FAILS     7200 BTS/s         23736 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXX                                   
//	XXXXXX               XXXXXXXXXXXXXXX                    
//	      XXXXXX         XXXXXX         XXXXXXXXX           
//	      X     XXXXX          XXXXXX   XXXXX    XXXX       
//	X      X    X    XXX       XX    XXXX    XXXXXX  XXX    
//	 XXX    X    XX      X       X   X   X   XXX   XXXX XXXX
//	         XX      XXXX XXX     XX  X   X      X XX  XXXXX
//	    XX     X   XX   X    XX     X  X   XX   X X  XXXXXXX

//	BIBD   6  80  40   3  16       SAT       53 BTS      109 NDS       83 FAILS     5300 BTS/s         24683 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                        
//	XXXXXXXXXXXXXXXX                        XXXXXXXXXXXXXXXXXXXXXXXX                
//	                XXXXXXXXXXXXXXXX        XXXXXXXXXXXXXXXX        XXXXXXXX        
//	X               XXXXXXX         XXXXXXXXXXXXXXXX        XXXXXXX X       XXXXXXXX
//	 XXXXXXXX              XXXX     XXXX            XXXX    XXX    XXXXXXXXXXXXXXXXX
//	         XXXXXXX           XXXXX    XXXX            XXXX   XXXXX XXXXXXXXXXXXXXX

//	BIBD   7  70  30   3  10       SAT       48 BTS      128 NDS       98 FAILS     4800 BTS/s         25210 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                        
//	XXXXXXXXXX                    XXXXXXXXXXXXXXXXXXXX                    
//	          XXXXXXXXXX          XXXXXXXXXX          XXXXXXXXXX          
//	          X         XXXXXXXXX           XXXXXXXXXXXXXXXXXXX X         
//	X          X        XXXXXXX  XXXXXXXX   XX        X        X XXXXXXXXX
//	 XX         XXXXXXX        X         X    XXXXXXX  X       XXXXXXXXXXX
//	   XXXXXXX         X        XX        XX         X  XXXXXXX XXXXXXXXXX

//	BIBD  15  35   7   3   1       SAT      429 BTS      566 NDS      450 FAILS    21450 BTS/s         93852 CKS    0.02 s
//	XXXXXXX                            
//	X      XXXXXX                      
//	 X     X     XXXXX                 
//	X            X    XXXXX            
//	 X      X         X    XXXX        
//	  X    X          X        XXXX    
//	  X      X   X         X       XXX 
//	   X    X     X    X       X   X  X
//	    X     X    X    X   X  X    X  
//	     X     X  X      X  X   X    X 
//	   X        X   X    X   X   X  X  
//	      X     X  X      XX    X     X
//	      X   X      X X      X  X   X 
//	     X   X      X   X     X   X   X
//	    X      X     X    X  X    XX   

//	BIBD  12  44  11   3   2       SAT      455 BTS      614 NDS      495 FAILS    22750 BTS/s         84457 CKS    0.02 s
//	XXXXXXXXXXX                                 
//	XX         XXXXXXXXX                        
//	  XX       XX       XXXXXXX                 
//	  X X        XX     X      XXXXXX           
//	X  X         X      X            XXXXXXX    
//	    XX         XX    XX    X     XX     XX  
//	 X   X           X   X X    XX     XX     XX
//	      XX      X   X     XX X       X X  X X 
//	        XX X   X        X   X X       XXX  X
//	        XX        XX   X  X    XXXX       X 
//	      X   X     X  X     XX  XX     X X  X  
//	       X  X X    X    X        XX    X X X X

//	BIBD   7  77  33   3  11       SAT       57 BTS      149 NDS      110 FAILS     5700 BTS/s         28278 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                            
//	XXXXXXXXXXX                      XXXXXXXXXXXXXXXXXXXXXX                      
//	           XXXXXXXXXXX           XXXXXXXXXXX           XXXXXXXXXXX           
//	X          X          XXXXXXXXX             XXXXXXXXXX XXXXXXXXXX XX         
//	 X          X         XXXXXXX  XXXXXXXXX    XX        XXX        X  XXXXXXXXX
//	  XXXXXX     XXX             XX         XX    XX      X  XXXXX   XXXXXXXXXXXX
//	        XXX     XXXXXX         XX         XX    XXXXXX        XXX XXXXXXXXXXX

//	BIBD   9  60  20   3   5       SAT      144 BTS      266 NDS      182 FAILS     7199 BTS/s         38442 CKS    0.02 s
//	XXXXXXXXXXXXXXXXXXXX                                        
//	XXXXX               XXXXXXXXXXXXXXX                         
//	     XXXXX          XXXXX          XXXXXXXXXX               
//	     X    XXXX           XXXXX     XXXX      XXXXXX         
//	      X   X   XXX        X    XXXX     XXXX  XXX   XXX      
//	XX     X         XX       X   X   XX   X   XX   XXXXXXXX    
//	  X     XX X       X       X   XX X     XX X X  XX    X XXXX
//	              XXXX XXXX     XX      X       X X   XX   XXXXX
//	   XX       XX    X    XX        X   XX   X    X    XXXXXXXX

//	BIBD   7  84  36   3  12       SAT       28 BTS      130 NDS       90 FAILS     2800 BTS/s         27496 CKS       0 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                
//	XXXXXXXXXXXX                        XXXXXXXXXXXXXXXXXXXXXXXX                        
//	            XXXXXXXXXXXX            XXXXXXXXXXXX            XXXXXXXXXXXX            
//	            X           XXXXXXXXXXX             XXXXXXXXXXXXXXXXXXXXXXX X           
//	X            X          XXXXXXXXX  XXXXXXXXXX   XX          X          X XXXXXXXXXXX
//	 XXXXXXXXXXX                     X           X               XXXXXXXXXXXXXXXXXXXXXXX
//	              XXXXXXXXXX          XX          XX  XXXXXXXXXX            XXXXXXXXXXXX

//	BIBD  10  60  18   3   4       SAT      543 BTS      707 NDS      590 FAILS    18099 BTS/s        103090 CKS    0.03 s
//	XXXXXXXXXXXXXXXXXX                                          
//	XXXX              XXXXXXXXXXXXXX                            
//	    XXXX          XXXX          XXXXXXXXXX                  
//	    X   XXX           XXXX      XXX       XXXXXXX           
//	X       X  XX         X   XX       XXXX   XX     XXXXX      
//	 X   XX    X              X XX     X   X    XXXX X    XXXX  
//	            XXXX       XX   XX  X   X   XX      XXX   X   XX
//	  X    X     XX            X  XX XX  X    X X      X   XXXXX
//	   X     X      XX       X    XX      XXXX X X      XXXX  X 
//	          X    XXXXXXX                        XXX XXXX  XX X

//	BIBD  11  55  15   3   3       SAT     3496 BTS     3682 NDS     3536 FAILS    20564 BTS/s        626705 CKS    0.17 s
//	XXXXXXXXXXXXXXX                                        
//	XXX            XXXXXXXXXXXX                            
//	   XXX         XXX         XXXXXXXXX                   
//	   X  XX          XXX      XX       XXXXXXX            
//	        XXX          XXX     XXX    XXX    XXX         
//	    X   X  X      XXX        X  X          X  XXXXXX   
//	            XXXX     X  X  X    X   X  X   X  X     XXX
//	      X     XX  X        XX      XX  X  X   XXXXX      
//	     X     X  X       X  XX X      X  X  X  X    XX XX 
//	X        XX      X      X        X X   XX X  X   XXX  X
//	 XX    X               X      XX  X      XX    XX  XXXX

//	BIBD   7  91  39   3  13       SAT      107 BTS      221 NDS      174 FAILS     5349 BTS/s         45012 CKS    0.02 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                    
//	XXXXXXXXXXXXX                          XXXXXXXXXXXXXXXXXXXXXXXXXX                          
//	             XXXXXXXXXXXXX             XXXXXXXXXXXXX             XXXXXXXXXXXXX             
//	X            X            XXXXXXXXXXX               XXXXXXXXXXXX XXXXXXXXXXXX XX           
//	 X            X           XXXXXXXXX  XXXXXXXXXXX    XX          XXX          X  XXXXXXXXXXX
//	  XXXXXXXX     XXXX                  X          X     XXX       X  XXXXXXXX   XXXXXXXXXXXXX
//	          XXX      XXXXXXX         XX X          XXX     XXXXXXX           XXXXXXXXXXXXXXXX

//	BIBD   9  72  24   3   6       SAT       15 BTS      137 NDS       64 FAILS     1499 BTS/s         24217 CKS    0.01 s
//	XXXXXXXXXXXXXXXXXXXXXXXX                                                
//	XXXXXX                  XXXXXXXXXXXXXXXXXX                              
//	      XXXXXX            XXXXXX            XXXXXXXXXXXX                  
//	      X     XXXXX             XXXXXX      XXXXX       XXXXXXX           
//	X      X    X    XXX          X     XXXX       XXXXX  XXXX   XXXX       
//	 X           X      XXXX       X    XX  XXX    XXXX X     XXX    XXXX   
//	  X     X     XXX   X   X       X     X XX X       XXXX      XXX XX  XXX
//	         XX      XXX X   XX      XXX   X    X        X X  X     XXXXXXXX
//	   XXX     X          XX   XXX               XX         XX XXXXXX  XXXXX

//	BIBD  13  52  12   3   2       SAT      354 BTS      554 NDS      391 FAILS    11799 BTS/s         85296 CKS    0.03 s
//	XXXXXXXXXXXX                                        
//	XX          XXXXXXXXXX                              
//	  XX        XX        XXXXXXXX                      
//	  X X         XX      X       XXXXXXX               
//	X  X            X     X       X      XXXXXXX        
//	 X   X      X          X       XX    XX     XXXX    
//	      XX      XX       XX              XX   X   XXX 
//	     X  X        XX     XX       XX      XX X      X
//	    X X          X X      XX  X            X XX X  X
//	         XX     X  X     X  X  X X         X   X XX 
//	       X   X      X X     X X      XXX X       X   X
//	         X X X       X       X    XX  X  X   X  XX  
//	        X X         XX     X X  X   X   X X   X   X 

//	BIBD   9  84  28   3   7       SAT      154 BTS      302 NDS      215 FAILS     7700 BTS/s         65701 CKS    0.02 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                        
//	XXXXXXX                     XXXXXXXXXXXXXXXXXXXXX                                   
//	       XXXXXXX              XXXXXXX              XXXXXXXXXXXXXX                     
//	       X      XXXXXX               XXXXXXX       XXXXXX        XXXXXXXX             
//	        X     X     XXXXX          X      XXXXXX       XXXXXX  XXXXX   XXX          
//	X        X     X    X    XXXXXX     X     X     X      X     XXXX   XXXXX XXXXX     
//	 X        XXX   XXX            X     X     XXX  X       X    XX  XXX      XXXX XXXXX
//	  XXX              X X   XX     X     XX      X  XX      XXXX       XX   XXX  XXXXXX
//	     XX      X        XXX  X     XX     XX     X   XXXX               XXXX  XXXXXXXX

//	BIBD   9  96  32   3   8       SAT      455 BTS      605 NDS      518 FAILS    15166 BTS/s        105046 CKS    0.03 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                                
//	XXXXXXXX                        XXXXXXXXXXXXXXXXXXXXXXXX                                        
//	        XXXXXXXX                XXXXXXXX                XXXXXXXXXXXXXXXX                        
//	        X       XXXXXXX                 XXXXXXXX        XXXXXXX         XXXXXXXXX               
//	X        X      X      XXXXX            X       XXXXXX         XXXXXXX  XXXXXX   XXXX           
//	 X        X      X     X    XXXXXXXX     X      X     X        X      XXXXX   XXXXX  XXXXXX     
//	  XXX      X      XXXX              X     X      X    XX        XXXXXX        XXX  X XXXX  XXXXX
//	            X           XXXXXXX      XX    XXXX   X    XXXXX          X          X XXX   XXXXXXX
//	     XXX     XXX      X        X       X       X   XXX      XXX        X   XXX    X X XXXXXXXXXX

//	BIBD  10  90  27   3   6       SAT       17 BTS      215 NDS       88 FAILS     1699 BTS/s         35355 CKS    0.01 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXX                                                               
//	XXXXXX                     XXXXXXXXXXXXXXXXXXXXX                                          
//	      XXXXXX               XXXXXX               XXXXXXXXXXXXXXX                           
//	      X     XXXXX                XXXXXX         XXXXX          XXXXXXXXXX                 
//	X      X    X    XXX             X     XXXX          XXXXX     XXXX      XXXXXXX          
//	 X               X  XXXX          X    X   XXX       X    XXXXX    XXXXX XXX    XX        
//	  XXX   X         X     X          X    X  X    X     X   XXX  X   XX   X   XX    XXXXXXX 
//	         X   X      X   XXX              XX XXXX XXXX        X          X   XXXXXXXX     X
//	              X      XXX XXXXX      XXX                XX     X X    X   X    XX  X XXXXXX
//	     X    XX   XX  X          XXX             XX         X       XX   XX  XX    XX XXXXXXX

//	BIBD   9 108  36   3   9       SAT      144 BTS      313 NDS      216 FAILS     7200 BTS/s         72890 CKS    0.02 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                                        
//	XXXXXXXXX                           XXXXXXXXXXXXXXXXXXXXXXXXXXX                                             
//	         XXXXXXXXX                  XXXXXXXXX                  XXXXXXXXXXXXXXXXXX                           
//	         X        XXXXXXXX                   XXXXXXXXX         XXXXXXXX          XXXXXXXXXX                 
//	          X       X       XXXXXXX            X        XXXXXXXX         XXXXXXXX  XXXXXXX   XXX              
//	X          X       X      XXX    XXXXXXXX     X       X       X        X       XXXXXX   XXX   XXXXXXXXX     
//	 XXXXX      X       X        X   X                     XXXX    XXXXXX   X      X     X  X  XX XXXXXX   XXXXX
//	      XXX    XX      XXX      X          XXX   XX          X             XXXX         XX XX  XXXXX  XXXXXXXX
//	               XXX      XX     XX XX        X    XXXXX      XXX      XX      XX X          XXX    XXXXXXXXXX

//	BIBD  13  78  18   3   3       SAT     1392 BTS     1642 NDS     1449 FAILS    15466 BTS/s        408543 CKS    0.09 s
//	XXXXXXXXXXXXXXXXXX                                                            
//	XXX               XXXXXXXXXXXXXXX                                             
//	   XXX            XXX            XXXXXXXXXXXX                                 
//	   X  XX             XXX         XX          XXXXXXXXXX                       
//	X   X   X            X  X        XX                    XXXXXXXXXXX            
//	      XX X        X      XX        XX        X         XXX        XXXXXX      
//	     X  XX         X    X  X       X          XXX         X       X     XXXXXX
//	 X        XX        X    X           XX          XXX       XXX     XX   XXX   
//	            XXX       XX   X           XXX       X         X  XX   X XX    XX 
//	               XXX          XXX      X    XX X    XX      X   XX  X  X       X
//	          X XX              XX X    X X     X XX    X       X   XX    XX     X
//	  X           XX          X     X      X  XX    X    XX      X  XX  X  XX  X  
//	           X    XX            XXX       XX  X       XXXXXX               XX X 

//	BIBD  15  70  14   3   2       SAT     2321 BTS     2566 NDS     2375 FAILS    16578 BTS/s        493780 CKS    0.14 s
//	XXXXXXXXXXXXXX                                                        
//	XX            XXXXXXXXXXXX                                            
//	  XX          XX          XXXXXXXXXX                                  
//	  X X           XX        X         XXXXXXXXX                         
//	X  X            X          X        X        XXXXXXXXX                
//	 X  X         X            X        X                 XXXXXXXXX       
//	     XX        X X          X        X       XX       XX       XXXX   
//	     XX           XX         XX       XX       XX       XX         XX 
//	       XX         XX           XX    X  X        XX       XX   X     X
//	         XX         XX      X    X       XX      X X    X X     X  X  
//	       X X          X X   X       X     X           XX   X  X    XX X 
//	           XX        X X     X     X  X  X          XXX      X X     X
//	        X X            XX      X X         XX  XX            XX  XX   
//	           X X          XX      X X    X  X  XX             X X    X X
//	            XX        X  X    X    X       XX     XX   X   X    X   X 

//	BIBD  12  88  22   3   4       SAT    51462 BTS    51723 NDS    51530 FAILS    16984 BTS/s      10396572 CKS       3 s
//	XXXXXXXXXXXXXXXXXXXXXX                                                                  
//	XXXX                  XXXXXXXXXXXXXXXXXX                                                
//	    XXXX              XXXX              XXXXXXXXXXXXXX                                  
//	    X   XXX               XXXX          XXX           XXXXXXXXXXX                       
//	X    X  X  X              XX  X            XXX        X          XXXXXXXXXXX            
//	 X    X  X  X         X        XX       X  X           XX        XXX        XXXXXXXX    
//	  X          XXX            XX X         X    XXX      X            XXXX    XX      XXXX
//	       X  X  X  X                XXXX         XX X       XXX        X   XXX   XXXX      
//	                 XXXX         X X    XX   X XX    X     XXX                X  XX    XXXX
//	              XXXX               XXX X            XXXXX     XXX  XX  X            XXX   
//	           X      XX X XX           X  X        XX          XX XX  X  X X       X XX XX 
//	   X        X       XX   X            XX           XXX     X  XXX      X XXXXX   X     X

//	BIBD   9 120  40   3  10       SAT      115 BTS      301 NDS      198 FAILS     5750 BTS/s         70373 CKS    0.02 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                                                
//	XXXXXXXXXX                              XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                  
//	          XXXXXXXXXX                    XXXXXXXXXX                    XXXXXXXXXXXXXXXXXXXX                              
//	          X         XXXXXXXXX                     XXXXXXXXXX          XXXXXXXXX           XXXXXXXXXXX                   
//	X          X        X        XXXXXXX              X         XXXXXXXX           XXXXXXXXX  XXXXXXXX   XXXX               
//	 X          X        X       XXXX   XXX XXXXXX     X                XX         X        XXXXXXX   XXX    XXXXXXXXXX     
//	  X          XXXXX                  XXXX            XXXXXXXX        X           XXXXX          X  X  XXXXXXXXX     XXXXX
//	   XX             X   XXXXXXX                 XXXX          XXX      XX              XX XX      X  X XXXXX    XXXXXXXXXX
//	     XXXXX         X             XXX   X                       XXXXX   XXXXXXXX        X         X  X     XXXXXXXXXXXXXX

//	BIBD  19  57   9   3   1       SAT     2562 BTS     2841 NDS     2601 FAILS    11139 BTS/s        847752 CKS    0.23 s
//	XXXXXXXXX                                                
//	X        XXXXXXXX                                        
//	 X       X       XXXXXXX                                 
//	X                X      XXXXXXX                          
//	 X        X             X      XXXXXX                    
//	  X      X               X     X     XXXXX               
//	  X       X      X                        XXXXXX         
//	   X       X      X     X            X    X     XXX      
//	    X       X      X     X      X         X        XXX   
//	     X       X      X     X     X     X    X    X     X  
//	      X     X        X    X      X     X    X    X     X 
//	       X     X     X       X      X    X     X    X     X
//	        X     X      X      X      X    X    X  X  X     
//	      X    X          X      X      X   X  X        X   X
//	       X      X   X           X  X       X    X     X X  
//	     X         X      X    X       X X        X      X X 
//	        X       X      X      X     X X     X     X  X   
//	    X          X       X    X  X               X X    X X
//	   X            X   X        X    X      X     X   X   X 

//	BIBD  10 120  36   3   8       SAT      661 BTS      885 NDS      745 FAILS    11016 BTS/s        182643 CKS    0.06 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                                                    
//	XXXXXXXX                            XXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                        
//	        XXXXXXXX                    XXXXXXXX                    XXXXXXXXXXXXXXXXXXXX                                    
//	        X       XXXXXXX                     XXXXXXXX            XXXXXXX             XXXXXXXXXXXXX                       
//	X        X      X      XXXXX                X       XXXXXX             XXXXXXX      XXXXXX       XXXXXXXX               
//	 X        X            X    XXXXX                   X     XXXXXX       X      XXXXXXX     XXXXXXXXXXX    XX             
//	  XXXX                  XX  XX      X        X       X    X     XXXXXX        X           X      X   XXXXXXXXXXXXXXX    
//	           XXX   XXXX     X          XX               XXXX XX           XX     X           XXXX   X        XXXXXXXX XXXX
//	              X      X     X  XX XXX   X      XXXXXX         X            XXX   XXX            X   X XXX   XXXX    XXXXX
//	      XX       X      X         XXXX    XXXX                  XX      X      X     X XXXXX      X   X   XXX    XXXXXXXXX

//	BIBD  11 110  30   3   6       SAT      354 BTS      628 NDS      428 FAILS    11800 BTS/s        120199 CKS    0.03 s
//	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                                                
//	XXXXXX                        XXXXXXXXXXXXXXXXXXXXXXXX                                                        
//	      XXXXXX                  XXXXXX                  XXXXXXXXXXXXXXXXXX                                      
//	      X     XXXXX                   XXXXXX            XXXXX             XXXXXXXXXXXXX                         
//	            X    XXXXX              X     XXXXX            XXXXXX       XXXX         XXXXXXXX                 
//	             X   X    XXXX                X    XXXXX       X     XXXXX      XXXXX    XXX     XXXX             
//	XX            XX  X       X          XX        XX           XXXX XX         XX          X        XXXXXXXXXXX  
//	       X        X  XX XX      XX       X   XX       X           X     XX         XXXX   X    XXXXXXXXX      X 
//	        X            X  X  XXX  XXXX         X      X X                 XX    XXX        XX  XX  XX   XXXX  XX
//	         XXX             XXX            XX       XXX X             X  XX  XX     XX      XXXX  X   X  XX  XX X
//	  XXXX                      XX                X      X XXXX         XX             XXXXX   XX   X   XX  XXXXXX

//	BIBD  16  80  15   3   2       SAT     6701 BTS     7042 NDS     6779 FAILS    15229 BTS/s       1486766 CKS    0.44 s
//	XXXXXXXXXXXXXXX                                                                 
//	XX             XXXXXXXXXXXXX                                                    
//	  XX           XX           XXXXXXXXXXX                                         
//	X   X            X          XX         XXXXXXXXXX                               
//	 X  X            X            XX                 XXXXXXXXXX                     
//	  X  X            XX        X          X         XX        XXXXXXX              
//	      XX          X X        XX         X        X                XXXXXXX       
//	   X X              XX          X        XX        XX      X      X      XXXX   
//	        XX            XX       X X      XX         X        XX     X         XXX
//	          XX          X X         XX      XX        XX        XX    XX       X  
//	          X X  X       X          X    X    X         XX        X     XX XX   X 
//	             XX X        X          X       XX    X  X      X       X X    XX  X
//	            XX          XX      XX            XX        XX X     X     XX    X  
//	        X  X         X    X         XX        X X     X   X   X X  X    X  X    
//	      X  X         X       X         XX    X    X      XX        X   X   X  X  X
//	       X      X           XX       X  X      X X         XX  X X  X       X   X 

//	BIBD  13 104  24   3   4       SAT     2676 BTS     3015 NDS     2744 FAILS    14084 BTS/s        647641 CKS    0.19 s
//	XXXXXXXXXXXXXXXXXXXXXXXX                                                                                
//	XXXX                    XXXXXXXXXXXXXXXXXXXX                                                            
//	    XXXX                XXXX                XXXXXXXXXXXXXXXX                                            
//	    X   XXX                 XXXX            XXX             XXXXXXXXXXXXX                               
//	X          XXX              XX  X           X  XXX          X            XXXXXXXXXXXX                   
//	     XX XX              X       XXX               X          XX          XXX         XXXXXXXXXX         
//	 X     X      XX              XX   X         XX   X                         XXXX     XXX       XXXXXXX  
//	                XXXX                XXXX       X   XXX      XXXX            X   X       XX     XXX    X 
//	  XX      X   X          X              X          X  XX        XXX             XXXX    XXXX      XXX  X
//	                    XXXX         XX      XX     XX  XX          X  XXX       X   X   X    X    X     XXX
//	               XX   XX              XX   X X          X XXX      X    XXXXX       X X XX        X      X
//	           X     XX   X            X  X X  X           XXX X   X   XXX        XX    X      XXXX   X     
//	            XX     X   X  XX           X  X               XX      X   XXX  X       X        XXX  X XXXX 

//	BIBD   8  14   7   4   3       SAT       17 BTS       51 NDS       25 FAILS     1700 BTS/s          4250 CKS       0 s
//	XXXXXXX       
//	XXX    XXXX   
//	   XXX XXX X  
//	X  X  XX  XXX 
//	 X  X X X XX X
//	X    XX XX  XX
//	 XX  X X   XXX
//	  XXX    XX XX

//	BIBD  11  11   5   5   2       SAT       19 BTS       56 NDS       30 FAILS     1900 BTS/s          5167 CKS       0 s
//	XXXXX      
//	XX   XXX   
//	  XX XX X  
//	X X    XXX 
//	 X X   XX X
//	  X XX X  X
//	X   X X X X
//	 X  XX  XX 
//	 XX   X  XX
//	X  X X   XX
//	   XX XX X 

//	BIBD  10  15   6   4   2       SAT      133 BTS      180 NDS      151 FAILS    13300 BTS/s         22456 CKS       0 s
//	XXXXXX         
//	XX    XXXX     
//	  XX  XX  XX   
//	  X X   XXX X  
//	X  X    X  XXX 
//	X    X   XXX  X
//	  X  XX X    XX
//	 X   X X  X XX 
//	   XX  X X   XX
//	 X  X X    XX X

//	BIBD   9  18   8   4   3       SAT       76 BTS      128 NDS       90 FAILS     7600 BTS/s         10996 CKS       0 s
//	XXXXXXXX          
//	XXX     XXXXX     
//	   XXX  XXX  XX   
//	X  X  X    XXXXX  
//	 X  X X X  X X  XX
//	X    X X X  XX  XX
//	 X   X X  XX  XXX 
//	  X X  XX   X XX X
//	  XX  X  XX    XXX

//	BIBD  13  13   4   4   1       SAT        4 BTS       68 NDS       13 FAILS      400 BTS/s          6661 CKS       0 s
//	XXXX         
//	X   XXX      
//	 X  X  XX    
//	  X  X X X   
//	 X    X  XX  
//	X      X  XX 
//	   X X  X X  
//	  X X     X X
//	X       XX  X
//	   XX    X X 
//	 X   X     XX
//	  X   X X  X 
//	   X  XX    X

//	BIBD  10  18   9   5   4       SAT     2043 BTS     2101 NDS     2067 FAILS    25537 BTS/s        254501 CKS    0.08 s
//	XXXXXXXXX         
//	XXXX     XXXXX    
//	    XXXX XXXX X   
//	X   XX  XXX  X XX 
//	XX    X X  XX XXX 
//	 XX    XXXX   XX X
//	X  XX  X   X XXX X
//	 XX XX      XXX XX
//	   X  XXXX  XX  XX
//	  XX XX   XX   XXX

//	BIBD   8  28  14   4   6       SAT       44 BTS       95 NDS       66 FAILS     4400 BTS/s         11713 CKS       0 s
//	XXXXXXXXXXXXXX              
//	XXXXXX        XXXXXXXX      
//	      XXXXXX  XXXXXX  XX    
//	X     XXX   XXXXX   XX  XXX 
//	 X       XXXXX   XXXXX  XX X
//	  XXX    X  XXXX X    XXX XX
//	  XX XXX  X       X XXXX XXX
//	XX  XX  X  X    X  X  XXXXXX

//	BIBD  15  15   7   7   3       SAT      444 BTS      528 NDS      469 FAILS    22199 BTS/s         68466 CKS    0.02 s
//	XXXXXXX        
//	XXX    XXXX    
//	   XXX XXX X   
//	X  X  XX  XXX  
//	X   X X XX  XX 
//	 X  XX X  X XX 
//	  X  XX X XX X 
//	X X X  X   X XX
//	 X X  XXX    XX
//	  XXX   X X X X
//	XX   X  X  XX X
//	  X  XXX X  X X
//	 X  X X  XXX  X
//	 XXX     X XXX 
//	X  X X   XX  XX

//	BIBD  11  22  10   5   4       SAT     9367 BTS     9457 NDS     9399 FAILS    26019 BTS/s       1206778 CKS    0.36 s
//	XXXXXXXXXX            
//	XXXX      XXXXXX      
//	    XXXX  XXXX  XX    
//	X   X   XXX   XXXXX   
//	 X   X  XX XX X X  XX 
//	  XX XX   X    XX XXX 
//	 X  X X X  X X X  XX X
//	X X    X X X X  X X XX
//	 X X   X XX X    XXX X
//	   X X XX    XXX X  XX
//	X X X X     X X  X XXX

//	BIBD  16  16   6   6   2       SAT      368 BTS      463 NDS      392 FAILS    18399 BTS/s         79714 CKS    0.02 s
//	XXXXXX          
//	XX    XXXX      
//	  XX  XX  XX    
//	  X X   XXX X   
//	X  X    X X  XX 
//	X X      X X X X
//	 X X    X  XX  X
//	    XXX X  X X  
//	  X  X XX     XX
//	   X X X X  XX  
//	X    XX   X X  X
//	X   X  X   XX X 
//	   XX X  X    XX
//	 X  X  X  X  X X
//	 X   X   XXX  X 
//	 XX   X     XXX 

//	BIBD  12  22  11   6   5       SAT     1640 BTS     1745 NDS     1673 FAILS    27333 BTS/s        216342 CKS    0.06 s
//	XXXXXXXXXXX           
//	XXXXX      XXXXXX     
//	     XXXXX XXXXX X    
//	XX   XX   XXX   XXXX  
//	X X    XX X  XX XXX X 
//	 X XX  X X X X   XXXX 
//	X XX X   X  X  X XX XX
//	  X XX  X XX X X  XX X
//	X  X  XX  XX  XX   XXX
//	 XX  XX  X   XX X  XXX
//	 X  X XXX   X  XX X XX
//	   XX   XXX X X XX X X

//	BIBD  10  30  12   4   4       SAT      403 BTS      502 NDS      437 FAILS    20149 BTS/s         65943 CKS    0.02 s
//	XXXXXXXXXXXX                  
//	XXXX        XXXXXXXX          
//	    XXXX    XXXX    XXXX      
//	    X   XXX     XXXXXXX X     
//	X    X  X  XX   XX  X  X XXX  
//	 X   X   X X X    XX X X X  XX
//	X   X X  X  X X   X     X XXXX
//	  XX  X   X  X     X  XXXXXX  
//	  XX   X   X   XX   XX  X X XX
//	 X     XX X   XX X    X  X XXX

