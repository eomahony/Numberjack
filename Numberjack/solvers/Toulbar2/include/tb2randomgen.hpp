/* \file tb2randomgen.hpp
 * \brief Random WCSP generator.
 */

#ifndef TB2RANDOMGEN_H_
#define TB2RANDOMGEN_H_

#include "tb2wcsp.hpp"

class naryRandom {
public:

    WCSP& wcsp;

    naryRandom(WCSP* wcspin, int seed = 0) : wcsp(*wcspin) { mysrand(seed);  }
    ~naryRandom() {}

    int n,m;

    bool connected();
    void generateNaryCtr( vector<int>& indexs, long nogoods, Cost costMin = SMALL_COST, Cost costMax = MEDIUM_COST);
    void generateTernCtr( int i, int j, int k, long p, Cost costMin = SMALL_COST, Cost costMax = MEDIUM_COST );
    void generateBinCtr( int i, int j, long p, Cost costMin = SMALL_COST, Cost costMax = MEDIUM_COST );
    void generateSubModularBinCtr( int i, int j, Cost costMin = SMALL_COST, Cost costMax = MEDIUM_COST );
    void Input( int in_n, int in_m, vector<int>& p, bool forceSubModular = false );

    void ini( vector<int>& index, int arity );
    long toIndex( vector<int>& index );
    int inc( vector<int>& index, int i );
    bool inc( vector<int>& index );

};






#endif /*TB2RANDOMGEN_H_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

