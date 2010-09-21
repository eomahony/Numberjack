

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <sol.h>


using namespace std;


static const int north = 0;
static const int east = 1;
static const int south = 2;
static const int west = 3;


typedef int itile;

int getNorth(itile t) {
  return (t & 0xff);
}

int getEast(itile t) {
  return ((t & 0xff00) >> 8);
}

int getSouth(itile t) {
  return ((t & 0xff0000) >> 16);
}

int getWest(itile t) {
  return ((t & 0xff000000) >> 24);
}

int setNorth(itile t, int val) {
  return (t |= val);
}

int setEast(itile t, int val) {
  return (t |= (val << 8));
}

int setSouth(itile t, int val) {
  return (t |= (val << 16));
}

int setWest(itile t, int val) {
  return (t |= (val << 24));
}



struct tile {
  int cols[4];
};

// get the tile at position i,j in the solution
tile getTile( int*** colors, int i, int j, int K, int N )
{
  tile sq;
  sq.cols[north] = ( i     ? colors[0][j][i-1] : K );
  sq.cols[east]  = ( j<N-1 ? colors[1][i][j]   : K );
  sq.cols[south] = ( i<N-1 ? colors[0][j][i]   : K );
  sq.cols[west]  = ( j     ? colors[1][i][j-1] : K );
  return sq;
}

// get a unique id for a tile and an angle (0,1,2,3)
int getIndex( tile t, int angle, int K )
{
  int index=0, k=0;
  for(int i=angle; k<4; i=(++i%4)) 
    index += (t.cols[i] * pow((K+1),k++));
  return index;
}

int main(int argc, char *argv[])
{  
  long int Seed = 11041979;
  srand( Seed );

  //itile t;



  // struct
  int i, j, k, N = ( argc > 1 ? atoi(argv[1]) : 4 ), 
    K = ( argc > 2 ? atoi(argv[2]) : 4 );
  int allcells[2*N*(N-1)];
  int nTiles = (K+1)*(K+1)*(K+1)*(K+1);
  int tiles[nTiles];
  std::fill( tiles, tiles+nTiles, 0 );
  for( i=0; i<(2*N*(N-1)); ++i )  
    allcells[i] = i;
  int **colors[2];
  for( j=0; j<2; ++j) {
    colors[j] = new int*[N];
    for( i=0; i<N; ++i )
      colors[j][i] = new int[N-1]; 
  }


  // 
  int int2ind[nTiles];
  int ind2int[nTiles];
  int nbt=0;
  int int2ind_r[nTiles];
  int ind2int_r[nTiles];
  int nbt_r=0;

  
  ///////////////////////////////////////
  // shuffle the cells
  for( i=0; i<(2*N*(N-1)); ++i) {
    j = (rand() % ((2*N*(N-1)) - i));
    k = allcells[j];
    allcells[j] = allcells[i];
    allcells[i] = k;
  }

  // assign colors
  k=0;
  int x, y, z;  
  for( i=0; i<(2*N*(N-1)); ++i) {
    x = (allcells[i] % (N*(N-1)));
    y = ( x % N );
    z = ( x / N );    
    colors[(x == allcells[i])][y][z] = k;
    k = ((k+1)%K);
  }

  // get the set of tiles
  vector<tile> pieces;
  for( i=0; i<N; ++i ) 
    for( j=0; j<N; ++j ) {
      pieces.push_back( getTile( colors, i, j, K, N ) );
      for(k=0; k<4; ++k) {
	x = getIndex( pieces.back(), k, K );
	if( !k ) {
	  ++tiles[x];
	  int2ind[x] = nbt;
	  ind2int[nbt++] = x;
	}
	int2ind_r[x] = nbt_r;
	ind2int_r[nbt_r++] = x;
      }
    }
  ///////////////////////////////////////


  // display tiles and puzzle 
  int r0, r1, r2, r3;
  for( i=0; i<pieces.size(); ++i ) {
    r0 = getIndex( pieces[i], 0, K );
    r1 = getIndex( pieces[i], 1, K );
    r2 = getIndex( pieces[i], 2, K );
    r3 = getIndex( pieces[i], 3, K );
    cout << "   " << setw(6) << pieces[i].cols[north] << endl
	 << setw(6) << pieces[i].cols[west] << setw(6) << pieces[i].cols[east] 
	 << " => " << r0 << " (" << int2ind_r[r0] << ") "
	 << r1 << " (" << int2ind_r[r1] << ") "
	 << r2 << " (" << int2ind_r[r2] << ") "
	 << r3 << " (" << int2ind_r[r3] << ") " << endl
	 << "   " << setw(6) << pieces[i].cols[south] << endl << endl;
  }

  cout << "   " ;
  for(i=0; i<N; ++i) 
    cout << setw(6) << K; 
  cout << endl;
  for(i=0; i<N; ++i) {
    for(k=1; k>=0; --k) {
      if( k ) cout << setw(6) << K; 
      else cout << "   ";
      if(i<N-1 || k) for(j=0; j<N-k; ++j) { 	  
	  if( k )
	    cout << setw(6) << colors[k][i][j];
	  else
	    cout << setw(6) << colors[k][j][i];
	}
      if( k ) cout << setw(6) << K; 
      if(i<N-1 || k) cout << endl;
    }
  }
  for(i=0; i<N; ++i) 
    cout << setw(6) << K; 
  cout << endl;
  cout << endl;


  int nEdges = 2*N*(N-1);
  int nPieces = pieces.size();

  CSP model;

  VarArray edges( nEdges, 0, K-1 );
  //VarArray tiles( N*N, 0, nPieces-1 );
  
  








  for( j=0; j<2; ++j) 
    delete [] colors[j];
}




