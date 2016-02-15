/*
 * **************** Read BEP files **************************
 * 
 */

#include "toulbar2lib.hpp"
#include "tb2bep.hpp"

const bool BAC = true;

void BEP::read(const char *fileName, WCSP *wcsp)
{
    Cost top = MIN_COST;

    // open the file
    ifstream file(fileName);
    if (!file) {
        cerr << "Could not open file " << fileName << endl;
        exit(EXIT_FAILURE);
    }
    string str;
    getline(file,str);
    getline(file,str);
    file >> size;
    duration.resize(size);
    earliest.resize(size);
    latest.resize(size);
    revenue.resize(size);
    delay.resize(size*size);
    getline(file,str);
    for (int i=0; i<size; i++) {
        int pos;
        file >> pos;
        pos--;
        int x;
        file >> x;
        int y;
        file >> y;
        file >> duration[pos];
        file >> earliest[pos];
        earliest[pos] = max(0,earliest[pos]);
        file >> latest[pos];
        latest[pos] -= duration[pos];
        file >> revenue[pos];
        top += revenue[pos] * (size-1);
    }
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            file >> delay[i*size+j];
        }
    }

    wcsp->updateUb(top);

    /* create variables */
    for (int i=0; i<size; i++) {
        string varname = "t";
        varname += to_string(i+1);
        if (BAC) wcsp->makeIntervalVariable(varname, earliest[i],latest[i]+1);
        else {
            wcsp->makeEnumeratedVariable(varname, 0, latest[i]+1);
            wcsp->increase(i, earliest[i]);
        }
    }

    /* create binary special disjunction */
    for (int i=0; i<size; i++) {
        for (int j=i+1; j<size; j++) {
            if (BAC) {
                wcsp->postSpecialDisjunction(i,j,duration[i]+delay[i*size+j],duration[j]+delay[j*size+i],latest[i]+1,latest[j]+1,revenue[i],revenue[j]);
            } else {
                vector<Cost> costs((latest[i]+2) * (latest[j]+2), 0);
                for (int a=earliest[i];a<=latest[i]+1;a++) {
                    for (int b=earliest[j];b<=latest[j]+1;b++) {
                        if (a>latest[i] && b>latest[j]) costs[a * (latest[j]+2) + b] = revenue[i]+revenue[j];
                        else if (a>latest[i]) costs[a * (latest[j]+2) + b] = revenue[i];
                        else if (b>latest[j]) costs[a * (latest[j]+2) + b] = revenue[j];
                        else costs[a * (latest[j]+2) + b] = (((a>=b+duration[j]+delay[j*size+i])||(b>=a+duration[i]+delay[i*size+j]))?0:wcsp->getUb()*MEDIUM_COST);
                    }
                }
                wcsp->postBinaryConstraint(i,j,costs);
            }
        }
    }
    wcsp->sortConstraints();
    wcsp->histogram();

    if (ToulBar2::verbose >= 0) {
        cout << "Read BEP with " << size << " photographs and total gain " << top/(size-1) << endl;
    }
}

void BEP::printSolution(WCSP *wcsp)
{
    cout << "Id \t\tRev\t\tTime\t\tMin\t\tMax\t\tDuration\tDelay\t\tSlack" << endl;
    int cost = 0;
    int nbphotos = 0;
    Value lasttime = -1;
    int lastcurr = -1;
    for (int i=0; i<size; i++) {
        Value time = MAX_VAL;
        int curr = -1;
        for (int j=0; j<size; j++) {
            if (wcsp->getValue(j) <= ToulBar2::bep->latest[j] && wcsp->getValue(j) < time && wcsp->getValue(j) > lasttime) {
                curr = j;
                time = wcsp->getValue(j);
            }
        }
        if (curr >= 0) {
            cout << curr+1 << "\t\t" << ToulBar2::bep->revenue[curr] << "\t\t" << wcsp->getValue(curr) << "\t\t" << ToulBar2::bep->earliest[curr] << "\t\t" << ToulBar2::bep->latest[curr] << "\t\t" << ToulBar2::bep->duration[curr];
            if (lastcurr>=0) {
                cout << "\t\t" << ToulBar2::bep->delay[lastcurr*ToulBar2::bep->size+curr];
                int slack = time - max(ToulBar2::bep->earliest[curr],(lasttime + ToulBar2::bep->duration[lastcurr] + ToulBar2::bep->delay[lastcurr*ToulBar2::bep->size+curr]));
                cout << "\t\t" << slack;
                if (slack < 0) cout << " " << "**********";
            }
            cout << endl;
            cost += ToulBar2::bep->revenue[curr];
            nbphotos++;
            lasttime = time;
            lastcurr = curr;
        } else break;
    }
    cout << "Gain = " << cost << "\t\tNbPhotos = " << nbphotos << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

