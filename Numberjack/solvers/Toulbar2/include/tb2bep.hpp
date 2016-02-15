/** \file tb2bep.hpp
 *  \brief BEP benchmark: selecting and scheduling earth observations for agile satellite
 * 
 */

#ifndef TB2BEP_HPP_
#define TB2BEP_HPP_

#include "tb2wcsp.hpp"

class BEP {
public:
    int size;
    vector<int> duration;
    vector<int> earliest;
    vector<int> latest;
    vector<int> revenue;
    vector<int> delay;

    BEP() : size(0) {}

    void read(const char *fileName, WCSP *wcsp);
    void printSolution(WCSP *wcsp);
};

#endif /*TB2BEP_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

