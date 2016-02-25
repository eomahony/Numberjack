/*
 * **************** Read wcsp format files **************************
 *
 */

#include "tb2wcsp.hpp"
#include "tb2enumvar.hpp"
#include "tb2pedigree.hpp"
#include "tb2haplotype.hpp"
#include "tb2bep.hpp"
#include "tb2naryconstr.hpp"
#include "tb2randomgen.hpp"
#include "tb2globaldecomposable.hpp"
#include <list>

typedef struct {
    EnumeratedVariable *var;
    vector<Cost> costs;
} TemporaryUnaryConstraint;

/**
 * \defgroup wcspformat Weighted Constraint Satisfaction Problem file format (wcsp)
 *
 * It is a text format composed of a list of numerical and string terms separated by spaces.
 * Instead of using names for making reference to variables, variable
 * indexes are employed. The same for domain values. All indexes start at
 * zero.
 *
 * Cost functions can be defined in intention (see below) or in extension, by their list of
 * tuples. A default cost value is defined per function in order to
 * reduce the size of the list. Only tuples with a different cost value
 * should be given (not mandatory). All the cost values must be positive. The arity of a cost function in extension may be
 * equal to zero. In this case, there is no tuples and the default cost
 * value is added to the cost of any solution. This can be used to represent
 * a global lower bound constant of the problem.
 *
 * The wcsp file format is composed of three parts: a problem header, the list of
 * variable domain sizes, and the list of cost functions.
 *
 * - Header definition for a given problem:
 * \verbatim
 <Problem name>
 <Number of variables (N)>
 <Maximum domain size>
 <Number of cost functions>
 <Initial global upper bound of the problem (UB)>
 \endverbatim
 * The goal is to find an assignment of all the variables with minimum total cost,
 * strictly lower than UB.
 * Tuples with a cost greater than or equal to UB are forbidden (hard constraint).
 *
 * - Definition of domain sizes
 * \verbatim
 <Domain size of variable with index 0>
 ...
 <Domain size of variable with index N - 1>
 \endverbatim
 * \note domain values range from zero to \e size-1
 * \note a negative domain size is interpreted as a variable with an interval domain in \f$[0,-size-1]\f$
 * \warning variables with interval domains are restricted to arithmetic and disjunctive cost functions in intention (see below)
 * - General definition of cost functions
 *   - Definition of a cost function in extension
 * \verbatim
 <Arity of the cost function>
 <Index of the first variable in the scope of the cost function>
 ...
 <Index of the last variable in the scope of the cost function>
 <Default cost value>
 <Number of tuples with a cost different than the default cost>
 \endverbatim
 * followed by for every tuple with a cost different than the default cost:
 * \verbatim
 <Index of the value assigned to the first variable in the scope>
 ...
 <Index of the value assigned to the last variable in the scope>
 <Cost of the tuple>
 \endverbatim
 * \note Shared cost function: A cost function in extension can be shared by several cost functions with the same arity (and same domain sizes) but different scopes. In order to do that, the cost function to be shared must start by a negative scope size. Each shared cost function implicitly receives an occurrence number starting from 1 and incremented at each new shared definition. New cost functions in extension can reuse some previously defined shared cost functions in extension by using a negative number of tuples representing the occurrence number of the desired shared cost function. Note that default costs should be the same in the shared and new cost functions. Here is an example of 4 variables with domain size 4 and one AllDifferent hard constraint decomposed into 6 binary constraints.
 *   - Shared CF used inside a small example in wcsp format:
 * \code
 AllDifferentDecomposedIntoBinaryConstraints 4 4 6 1
 4 4 4 4
 -2 0 1 0 4
 0 0 1
 1 1 1
 2 2 1
 3 3 1
 2 0 2 0 -1
 2 0 3 0 -1
 2 1 2 0 -1
 2 1 3 0 -1
 2 2 3 0 -1
 \endcode
 *   - Definition of a cost function in intension by replacing the default cost value by -1 and by giving its keyword name and its K parameters
 * \verbatim
 <Arity of the cost function>
 <Index of the first variable in the scope of the cost function>
 ...
 <Index of the last variable in the scope of the cost function>
 -1
 <keyword>
 <parameter1>
 ...
 <parameterK>
 \endverbatim
 *   .
 * .
 * Possible keywords of cost functions defined in intension followed by their specific parameters:
 * - >= \e cst \e delta to express soft binary constraint \f$x \geq y + cst\f$ with associated cost function \f$max( (y + cst - x \leq delta)?(y + cst - x):UB , 0 )\f$
 * - > \e cst \e delta to express soft binary constraint \f$x > y + cst\f$ with associated cost function  \f$max( (y + cst + 1 - x \leq delta)?(y + cst + 1 - x):UB , 0 )\f$
 * - <= \e cst \e delta to express soft binary constraint \f$x \leq y + cst\f$ with associated cost function  \f$max( (x - cst - y \leq delta)?(x - cst - y):UB , 0 )\f$
 * - < \e cst \e delta to express soft binary constraint \f$x < y + cst\f$ with associated cost function  \f$max( (x - cst + 1 - y \leq delta)?(x - cst + 1 - y):UB , 0 )\f$
 * - = \e cst \e delta to express soft binary constraint \f$x = y + cst\f$ with associated cost function  \f$(|y + cst - x| \leq delta)?|y + cst - x|:UB\f$
 * - disj \e cstx \e csty \e penalty to express soft binary disjunctive constraint \f$x \geq y + csty \vee y \geq x + cstx\f$ with associated cost function \f$(x \geq y + csty \vee y \geq x + cstx)?0:penalty\f$
 * - sdisj \e cstx \e csty \e xinfty \e yinfty \e costx \e costy to express a special disjunctive constraint with three implicit hard constraints \f$x \leq xinfty\f$ and \f$y \leq yinfty\f$ and \f$x < xinfty \wedge y < yinfty \Rightarrow (x \geq y + csty \vee y \geq x + cstx)\f$ and an additional cost function \f$((x = xinfty)?costx:0) + ((y= yinfty)?costy:0)\f$
 * - Global cost functions using a flow-based propagator:
 *     - salldiff var|dec|decbi \e cost to express a soft alldifferent constraint with either variable-based (\e var keyword) or decomposition-based (\e dec and \e decbi keywords) cost semantic with a given \e cost per violation (\e decbi decomposes into a binary cost function complete network)
 *     - sgcc var|dec|wdec \e cost \e nb_values (\e value \e lower_bound \e upper_bound (\e shortage_weight \e excess_weight)?)* to express a soft global cardinality constraint with either variable-based (\e var keyword) or decomposition-based (\e dec keyword) cost semantic with a given \e cost per violation and for each value its lower and upper bound (if \e wdec then violation cost depends on each value shortage or excess weights)
 *     - ssame \e cost \e list_size1 \e list_size2 (\e variable_index)* (\e variable_index)* to express a permutation constraint on two lists of variables of equal size (implicit variable-based cost semantic)
 *     - sregular var|edit \e cost \e nb_states \e nb_initial_states (\e state)* \e nb_final_states (\e state)* \e nb_transitions (\e start_state \e symbol_value \e end_state)* to express a soft regular constraint with either variable-based (\e var keyword) or edit distance-based (\e edit keyword) cost semantic with a given \e cost per violation followed by the definition of a deterministic finite automaton with number of states, list of initial and final states, and list of state transitions where symbols are domain values
 *     .
 * - Global cost functions using a dynamic programming DAG-based propagator:
 *     - sregulardp var \e cost \e nb_states \e nb_initial_states (\e state)* \e nb_final_states (\e state)* \e nb_transitions (\e start_state \e symbol_value \e end_state)* to express a soft regular constraint with a variable-based (\e var keyword) cost semantic with a given \e cost per violation followed by the definition of a deterministic finite automaton with number of states, list of initial and final states, and list of state transitions where symbols are domain values
 *     - sgrammar|sgrammardp var|weight \e cost \e nb_symbols \e nb_values \e start_symbol \e nb_rules ((0 \e terminal_symbol \e value)|(1 \e nonterminal_in \e nonterminal_out_left \e nonterminal_out_right)|(2 \e terminal_symbol \e value \e weight)|(3 \e nonterminal_in \e nonterminal_out_left \e nonterminal_out_right \e weight))* to express a soft/weighted grammar in Chomsky normal form
 *     - samong|samongdp var \e cost \e lower_bound \e upper_bound \e nb_values (\e value)* to express a soft among constraint to restrict the number of variables taking their value into a given set of values
 *     - max|smaxdp \e defCost \e nbtuples (\e variable \e value \e cost)* to express a weighted max cost function to find the maximum cost over a set of unary cost functions associated to a set of variables (by default, \e defCost if unspecified)
 *     - smst|smstdp hard to express a spanning tree hard constraint where each variable is assigned to its parent variable index in order to build a spanning tree (the root being assigned to itself)
 *     .
 * - Global cost functions using a cost function network-based propagator:
 *     - wregular \e nb_states \e nb_initial_states (\e state and cost enum)* \e nb_final_states (\e state and cost enum)* \e nb_transitions (\e start_state \e symbol_value \e end_state \e cost)* to express a wregular constraint with variable-based cost semantic with a given \e cost per violation followed by the definition of a deterministic finite automaton with number of states, list of initial and final states, and list of state transitions where symbols are domain values
 *     - walldiff hard|lin|quad \e cost to express a soft alldifferent constraint as a set of wamong hard constraint (\e hard keyword) or decomposition-based (\e lin and \e quad keywords) cost semantic with a given \e cost per violation
 *     - wgcc hard|lin|quad \e cost \e nb_values (\e value \e lower_bound \e upper_bound)* to express a soft global cardinality constraint as either a hard constraint (\e hard keyword) or with decomposition-based (\e lin and \e quad keyword) cost semantic with a given \e cost per violation and for each value its lower and upper bound
 *     - wsame hard|lin|quad \e cost to express a permutation constraint on two lists of variables of equal size (implicitly concatenated in the scope) using implicit decomposition-based cost semantic
 *     - wsamegcc hard|lin|quad \e cost \e nb_values (\e value \e lower_bound \e upper_bound)* to express the combination of a soft global cardinality constraint and a permutation constraint
 *     - wamong hard|lin|quad \e cost \e nb_values (\e value)* \e lower_bound \e upper_bound to express a soft among constraint to restrict the number of variables taking their value into a given set of values
 *     - wvaramong hard \e cost \e nb_values (\e value)* to express a hard among constraint to restrict the number of variables taking their value into a given set of values to be equal to the last variable in the scope
 *     - woverlap hard|lin|quad \e cost \e comparator \e righthandside  overlaps between two sequences of variables X, Y (i.e. set the fact that Xi and Yi take the same value (not equal to zero))
 *     - wsum hard|lin|quad \e cost \e comparator \e righthandside to express a soft sum constraint with unit coefficients to test if the sum of a set of variables matches with a given comparator and right-hand-side value
 *     - wvarsum hard \e cost \e comparator to express a hard sum constraint to restrict the sum to be \e comparator to the value of the last variable in the scope
 *
 *       Let us note <> the comparator, K the right-hand-side value associated to the comparator, and Sum the result of the sum over the variables. For each comparator, the gap is defined according to the distance as follows:
 *       -	if <> is == : gap = abs(K - Sum)
 *       -  if <> is <= : gap = max(0,Sum - K)
 *       -  if <> is < : gap = max(0,Sum - K - 1)
 *       -	if <> is != : gap = 1 if Sum != K and gap = 0 otherwise
 *       -  if <> is > : gap = max(0,K - Sum + 1);
 *       -	if <> is >= : gap = max(0,K - Sum);
 *       .
 *     .
 * .
 *
 * \warning The decomposition of wsum and wvarsum may use an exponential size (sum of domain sizes).
 * \warning  \e list_size1 and \e list_size2 must be equal in \e ssame.
 * \warning  Cost functions defined in intention cannot be shared.
 *
 * \note More about network-based global cost functions can be found here https://metivier.users.greyc.fr/decomposable/
 *
 * Examples:
 * - quadratic cost function \f$x0 * x1\f$ in extension with variable domains \f$\{0,1\}\f$ (equivalent to a soft clause \f$\neg x0 \vee \neg x1\f$): \code 2 0 1 0 1 1 1 1 \endcode
 * - simple arithmetic hard constraint \f$x1 < x2\f$: \code 2 1 2 -1 < 0 0 \endcode
 * - hard temporal disjunction\f$x1 \geq x2 + 2 \vee x2 \geq x1 + 1\f$: \code 2 1 2 -1 disj 1 2 UB \endcode
 * - soft_alldifferent({x0,x1,x2,x3}): \code 4 0 1 2 3 -1 salldiff var 1 \endcode
 * - soft_gcc({x1,x2,x3,x4}) with each value \e v from 1 to 4 only appearing at least v-1 and at most v+1 times: \code 4 1 2 3 4 -1 sgcc var 1 4 1 0 2 2 1 3 3 2 4 4 3 5 \endcode
 * - soft_same({x0,x1,x2,x3},{x4,x5,x6,x7}): \code 8 0 1 2 3 4 5 6 7 -1 ssame 1 4 4 0 1 2 3 4 5 6 7 \endcode
 * - soft_regular({x1,x2,x3,x4}) with DFA (3*)+(4*): \code 4 1 2 3 4 -1 sregular var 1 2 1 0 2 0 1 3 0 3 0 0 4 1 1 4 1 \endcode
 * - soft_grammar({x0,x1,x2,x3}) with hard cost (1000) producing well-formed parenthesis expressions: \code 4 0 1 2 3 -1 sgrammardp var 1000 4 2 0 6 1 0 0 0 1 0 1 2 1 0 1 3 1 2 0 3 0 1 0 0 3 1 \endcode
 * - soft_among({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^4(x_i \in \{1,2\}) < 1\f$ or \f$\sum_{i=1}^4(x_i \in \{1,2\}) > 3\f$: \code 4 1 2 3 4 -1 samongdp var 1000 1 3 2 1 2 \endcode
 * - soft max({x0,x1,x2,x3}) with cost equal to \f$\max_{i=0}^3((x_i!=i)?1000:(4-i))\f$: \code 4 0 1 2 3 -1 smaxdp 1000 4 0 0 4 1 1 3 2 2 2 3 3 1 \endcode
 * - wregular({x0,x1,x2,x3}) with DFA (a(ba)*c*): \code 4 0 1 2 3 -1 wregular 3 1 0 0 1 2 0 9 0 0 1 0 0 1 1 1 0 2 1 1 1 1 0 0 1 0 0 1 1 2 0 1 1 2 2 0 1 0 2 1 1 1 2 1 \endcode
 * - wamong ({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^4(x_i \in \{1,2\}) < 1\f$ or \f$\sum_{i=1}^4(x_i \in \{1,2\}) > 3\f$: \code 4 1 2 3 4 -1 wamong hard 1000 2 1 2 1 3 \endcode
 * - wvaramong ({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^3(x_i \in \{1,2\}) \neq x_4\f$: \code 4 1 2 3 4 -1 wvaramong hard 1000 2 1 2 \endcode
 * - woverlap({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^2(x_i = x_{i+2}) \geq 1\f$: \code 4 1 2 3 4 -1 woverlap hard 1000 < 1\endcode
 * - wsum ({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^4(x_i) \neq 4\f$: \code 4 1 2 3 4 -1 wsum hard 1000 == 4 \endcode
 * - wvarsum ({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^3(x_i) \neq x_4\f$: \code 4 1 2 3 4 -1 wvarsum hard 1000 == \endcode
 * .
 *
 * Latin Square 4 x 4 crisp CSP example in wcsp format:
 * \code
 latin4 16 4 8 1
 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
 4 0 1 2 3 -1 salldiff var 1
 4 4 5 6 7 -1 salldiff var 1
 4 8 9 10 11 -1 salldiff var 1
 4 12 13 14 15 -1 salldiff var 1
 4 0 4 8 12 -1 salldiff var 1
 4 1 5 9 13 -1 salldiff var 1
 4 2 6 10 14 -1 salldiff var 1
 4 3 7 11 15 -1 salldiff var 1
 \endcode
 *
 * 4-queens binary weighted CSP example with random unary costs in wcsp format:
 * \code
 4-WQUEENS 4 4 10 5
 4 4 4 4
 2 0 1 0 10
 0 0 5
 0 1 5
 1 0 5
 1 1 5
 1 2 5
 2 1 5
 2 2 5
 2 3 5
 3 2 5
 3 3 5
 2 0 2 0 8
 0 0 5
 0 2 5
 1 1 5
 1 3 5
 2 0 5
 2 2 5
 3 1 5
 3 3 5
 2 0 3 0 6
 0 0 5
 0 3 5
 1 1 5
 2 2 5
 3 0 5
 3 3 5
 2 1 2 0 10
 0 0 5
 0 1 5
 1 0 5
 1 1 5
 1 2 5
 2 1 5
 2 2 5
 2 3 5
 3 2 5
 3 3 5
 2 1 3 0 8
 0 0 5
 0 2 5
 1 1 5
 1 3 5
 2 0 5
 2 2 5
 3 1 5
 3 3 5
 2 2 3 0 10
 0 0 5
 0 1 5
 1 0 5
 1 1 5
 1 2 5
 2 1 5
 2 2 5
 2 3 5
 3 2 5
 3 3 5
 1 0 0 2
 1 1
 3 1
 1 1 0 2
 1 1
 2 1
 1 2 0 2
 1 1
 2 1
 1 3 0 2
 0 1
 2 1
 \endcode
 **/

void WCSP::read_wcsp(const char *fileName)
{
    char *Nfile2;
    Nfile2 = strdup(fileName);
    name = to_string(basename(Nfile2));

    if (ToulBar2::haplotype) {
        ToulBar2::haplotype->read(fileName, this);
        return;
    } else if (ToulBar2::pedigree) {
        if (!ToulBar2::bayesian) ToulBar2::pedigree->read(fileName, this);
        else ToulBar2::pedigree->read_bayesian(fileName, this);
        return;
    } else if (ToulBar2::uai) {
        read_uai2008(fileName);
        return;
    } else if (ToulBar2::xmlflag) {
        read_XML(fileName);
        return;
    } else if (ToulBar2::bep) {
        ToulBar2::bep->read(fileName, this);
        return;
    } else if (ToulBar2::wcnf) {
        read_wcnf(fileName);
        return;
    } else if (ToulBar2::qpbo) {
        read_qpbo(fileName);
        return;
    }
    string pbname;
    int nbvar,nbval,nbconstr;
    int nbvaltrue = 0;
    Cost top;
    int i,j,k,t, ic;
    string varname;
    int domsize;
    unsigned int a;
    unsigned int b;
    unsigned int c;
    Cost defval;
    Cost cost;
    int ntuples;
    int arity;
    string funcname;
    Value funcparam1;
    Value funcparam2;
    vector<TemporaryUnaryConstraint> unaryconstrs;
    Cost inclowerbound = MIN_COST;
    int maxarity = 0;
    vector< int > sharedSize;
    vector< vector<Cost> > sharedCosts;
    vector< vector<String> > sharedTuples;
    vector<String> emptyTuples;

    // open the file
    ifstream file(fileName);
    if (!file) {
        cerr << "Could not open file " << fileName << endl;
        exit(EXIT_FAILURE);
    }

    // read problem name and sizes
    file >> pbname;
    file >> nbvar;
    file >> nbval;
    file >> nbconstr;
    file >> top;
    if (ToulBar2::verbose >= 1) cout << "Read problem: " << pbname << endl;
    ToulBar2::nbvar= nbvar;


    assert(vars.empty());
    assert(constrs.empty());

    double K = ToulBar2::costMultiplier;
    if(top < MAX_COST / K)	top = top * K;
    else top = MAX_COST;
    updateUb(top);

    // read variable domain sizes
    for (i = 0; i < nbvar; i++) {
        string varname;
        varname = to_string(i);
        file >> domsize;
        if(domsize > nbvaltrue) nbvaltrue = domsize;
        if (ToulBar2::verbose >= 3) cout << "read variable " << i << " of size " << domsize << endl;
        int theindex = -1;

        if (domsize >= 0) theindex = makeEnumeratedVariable(varname,0,domsize-1);
        else theindex = makeIntervalVariable(varname,0,-domsize-1);
        assert(theindex == i);
    }

    // read each constraint
    for (ic = 0; ic < nbconstr; ic++) {
        file >> arity;
        if (!file) {
            cerr << "Warning: EOF reached before reading all the cost functions (initial number of cost functions too large?)" << endl;
            break;
        }
        bool shared = (arity < 0);
        if (shared) arity = -arity;
        if (arity > 3) {
            maxarity = max(maxarity,arity);
            if (ToulBar2::verbose >= 3) cout << "read " << arity << "-ary cost function " << ic << " on";
            int scopeIndex[arity]; // replace arity by MAX_ARITY in case of compilation problem
            for(i=0;i<arity;i++) {
                file >> j;
                if (ToulBar2::verbose >= 3) cout << " " << j;
                scopeIndex[i] = j;
            }
            if (ToulBar2::verbose >= 3) cout << endl;
            file >> defval;
            if (defval == -1) {
                string gcname;
                file >> gcname;
                if (gcname.substr(0,1) == "w") { // global cost functions decomposed into a cost function network
                    DecomposableGlobalCostFunction* decomposableGCF = DecomposableGlobalCostFunction::FactoryDGCF(gcname, arity, scopeIndex, file);
                    decomposableGCF->addToCostFunctionNetwork(this);
                }
                else { // monolithic global cost functions
                    postGlobalConstraint(scopeIndex, arity, gcname, file, &nbconstr);
                }

            } else {
                if(arity > MAX_ARITY)  { cerr << "Nary cost functions of arity > " << MAX_ARITY << " not supported" << endl; exit(EXIT_FAILURE); }
                file >> ntuples;
                int reusedconstr = -1;
                bool reused = (ntuples < 0);
                if (reused) {
                    reusedconstr = -ntuples-1;
                    if (reusedconstr >= (int) sharedSize.size()) {
                        cerr << "Shared cost function number " << reusedconstr << " not already defined! Cannot reuse it!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    ntuples = sharedSize[reusedconstr];
                }
                if((defval != MIN_COST) || (ntuples > 0))
                {
                    Cost tmpcost = defval*K;
                    if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb()) && getUb()<(MAX_COST / MEDIUM_COST)) tmpcost *= MEDIUM_COST;
                    int naryIndex = postNaryConstraintBegin(scopeIndex,arity,tmpcost);
                    NaryConstraint *nary = (NaryConstraint *) constrs[naryIndex];

                    Char buf[MAX_ARITY];
                    vector<String> tuples;
                    vector<Cost> costs;
                    for (t = 0; t < ntuples; t++) {
                        if (!reused) {
                            for(i=0;i<arity;i++) {
                                file >> j;
                                buf[i] = j + CHAR_FIRST;
                            }
                            buf[i] = '\0';
                            file >> cost;
                            Cost tmpcost = cost * K;
                            if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb()) && getUb()<(MAX_COST / MEDIUM_COST)) tmpcost *= MEDIUM_COST;
                            String tup = buf;
                            if (shared) {
                                tuples.push_back(tup);
                                costs.push_back(tmpcost);
                            }
                            postNaryConstraintTuple(naryIndex, tup, tmpcost);
                        } else {
                            postNaryConstraintTuple(naryIndex, sharedTuples[reusedconstr][t], sharedCosts[reusedconstr][t]);
                        }
                    }
                    if (shared) {
                        assert(ntuples == (int) costs.size());
                        sharedSize.push_back(costs.size());
                        sharedCosts.push_back(costs);
                        sharedTuples.push_back(tuples);
                    }

                    if (ToulBar2::preprocessNary>0) {
                        Cost minc = nary->getMinCost();
                        if (minc > MIN_COST) {
                            Cost defcost = nary->getDefCost();
                            if (CUT(defcost, minc)) nary->setDefCost(defcost - minc);
                            String tuple;
                            Cost cost;
                            nary->first();
                            while (nary->next(tuple,cost)) {
                                nary->setTuple(tuple, cost-minc, NULL);
                            }
                            if (ToulBar2::verbose >= 2) cout << "IC0 performed for cost function " << nary << " with initial minimum cost " << minc << endl;
                            inclowerbound += minc;
                        }
                    }
                    postNaryConstraintEnd(naryIndex);
                }
            }
        } else if (arity == 3) {
            maxarity = max(maxarity,arity);
            file >> i;
            file >> j;
            file >> k;
            if ((i == j) || (i == k) || (k == j)) {
                cerr << "Error: ternary cost function!" << endl;
                exit(EXIT_FAILURE);
            }
            file >> defval;
            if (defval >= MIN_COST) {
                assert(vars[i]->enumerated());
                assert(vars[j]->enumerated());
                assert(vars[k]->enumerated());
                EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
                EnumeratedVariable *y = (EnumeratedVariable *) vars[j];
                EnumeratedVariable *z = (EnumeratedVariable *) vars[k];
                if (ToulBar2::verbose >= 3) cout << "read ternary cost function " << ic << " on " << i << "," << j << "," << k << endl;
                file >> ntuples;
                if (ntuples < 0) {
                    int reusedconstr = -ntuples-1;
                    if (reusedconstr >= (int) sharedSize.size()) {
                        cerr << "Shared cost function number " << reusedconstr << " not already defined! Cannot reuse it!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    ntuples = sharedSize[reusedconstr];
                    assert(ntuples == (int) (x->getDomainInitSize() * y->getDomainInitSize() * z->getDomainInitSize()));
                    if((defval != MIN_COST) || (ntuples > 0)) postTernaryConstraint(i,j,k,sharedCosts[reusedconstr]);
                    continue;
                }
                vector<Cost> costs;
                for (a = 0; a < x->getDomainInitSize(); a++) {
                    for (b = 0; b < y->getDomainInitSize(); b++) {
                        for (c = 0; c < z->getDomainInitSize(); c++) {
                            Cost tmpcost = defval*K;
                            if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
                            costs.push_back(tmpcost);
                        }
                    }
                }
                for (t = 0; t < ntuples; t++) {
                    file >> a;
                    file >> b;
                    file >> c;
                    file >> cost;
                    Cost tmpcost = cost*K;
                    if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb()) && getUb()<(MAX_COST / MEDIUM_COST)) tmpcost *= MEDIUM_COST;
                    assert(a>= 0 && a < x->getDomainInitSize());
                    assert(b>= 0 && b < y->getDomainInitSize());
                    assert(c>= 0 && c < z->getDomainInitSize());
                    costs[a * y->getDomainInitSize() * z->getDomainInitSize() + b * z->getDomainInitSize() + c] = tmpcost;
                }
                if (shared) {
                    sharedSize.push_back(costs.size());
                    sharedCosts.push_back(costs);
                    sharedTuples.push_back(emptyTuples);
                }
                if((defval != MIN_COST) || (ntuples > 0)) postTernaryConstraint(i,j,k,costs);
            } else if (defval == -1) {
                int scopeIndex[3];
                scopeIndex[0] = i;
                scopeIndex[1] = j;
                scopeIndex[2] = k;
                string gcname;
                file >> gcname;
                postGlobalConstraint(scopeIndex, arity, gcname, file, &nbconstr);
            }
        } else if (arity == 2) {
            maxarity = max(maxarity,arity);
            file >> i;
            file >> j;
            if (ToulBar2::verbose >= 3) cout << "read binary cost function " << ic << " on " << i << "," << j << endl;
            if (i == j) {
                cerr << "Error: binary cost function with only one variable in its scope!" << endl;
                exit(EXIT_FAILURE);
            }
            file >> defval;
            if (defval >= MIN_COST) {
                assert(vars[i]->enumerated());
                assert(vars[j]->enumerated());
                EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
                EnumeratedVariable *y = (EnumeratedVariable *) vars[j];
                file >> ntuples;
                if (ntuples < 0) {
                    int reusedconstr = -ntuples-1;
                    if (reusedconstr >= (int) sharedSize.size()) {
                        cerr << "Shared cost function number " << reusedconstr << " not already defined! Cannot reuse it!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    ntuples = sharedSize[reusedconstr];
                    assert(ntuples == (int) (x->getDomainInitSize() * y->getDomainInitSize()));
                    if((defval != MIN_COST) || (ntuples > 0)) postBinaryConstraint(i,j,sharedCosts[reusedconstr]);
                    continue;
                }
                vector<Cost> costs;
                for (a = 0; a < x->getDomainInitSize(); a++) {
                    for (b = 0; b < y->getDomainInitSize(); b++) {
                        Cost tmpcost = defval*K;
                        if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb()) && getUb()<(MAX_COST / MEDIUM_COST)) tmpcost *= MEDIUM_COST;
                        costs.push_back(tmpcost);
                    }
                }
                for (k = 0; k < ntuples; k++) {
                    file >> a;
                    file >> b;
                    file >> cost;
                    Cost tmpcost = cost*K;
                    if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb()) && getUb()<(MAX_COST / MEDIUM_COST)) tmpcost *= MEDIUM_COST;
                    assert(a>= 0 && a < x->getDomainInitSize());
                    assert(b>= 0 && b < y->getDomainInitSize());
                    costs[a * y->getDomainInitSize() + b] = tmpcost;
                }
                if (shared) {
                    sharedSize.push_back(costs.size());
                    sharedCosts.push_back(costs);
                    sharedTuples.push_back(emptyTuples);
                }
                if((defval != MIN_COST) || (ntuples > 0)) postBinaryConstraint(i,j,costs);
            } else {
                file >> funcname;
                if (funcname == ">=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i,j,funcparam1,funcparam2);
                } else if (funcname == ">") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i,j,funcparam1 + 1,funcparam2);
                } else if (funcname == "<=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(j,i, -funcparam1,funcparam2);
                } else if (funcname == "<") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(j,i, -funcparam1 + 1,funcparam2);
                } else if (funcname == "=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i,j,funcparam1,funcparam2);
                    postSupxyc(j,i,-funcparam1,funcparam2);
                } else if (funcname == "disj") {
                    Cost funcparam3;
                    file >> funcparam1;
                    file >> funcparam2;
                    file >> funcparam3;
                    postDisjunction(i,j,funcparam1,funcparam2,funcparam3);
                } else if (funcname == "sdisj") {
                    Value funcparam3;
                    Value funcparam4;
                    Cost funcparam5;
                    Cost funcparam6;
                    file >> funcparam1;
                    file >> funcparam2;
                    file >> funcparam3;
                    file >> funcparam4;
                    file >> funcparam5;
                    file >> funcparam6;
                    postSpecialDisjunction(i,j,funcparam1,funcparam2,funcparam3,funcparam4,funcparam5,funcparam6);
                } else {
                    cerr << "Error: function " << funcname << " not implemented!" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        } else if (arity == 1) {
            maxarity = max(maxarity,arity);
            file >> i;
            if (ToulBar2::verbose >= 3) cout << "read unary cost function " << ic << " on " << i << endl;
            if (vars[i]->enumerated()) {
                EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
                file >> defval;
                file >> ntuples;
                TemporaryUnaryConstraint unaryconstr;
                unaryconstr.var = x;
                if (ntuples < 0) {
                    int reusedconstr = -ntuples-1;
                    if (reusedconstr >= (int) sharedSize.size()) {
                        cerr << "Shared cost function number " << reusedconstr << " not already defined! Cannot reuse it!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    ntuples = sharedSize[reusedconstr];
                    assert(ntuples == (int) x->getDomainInitSize());
                    unaryconstr.costs = sharedCosts[reusedconstr];
                    unaryconstrs.push_back(unaryconstr);
                    continue;
                }
                for (a = 0; a < x->getDomainInitSize(); a++) {
                    Cost tmpcost = defval*K;
                    if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb()) && getUb()<(MAX_COST / MEDIUM_COST)) tmpcost *= MEDIUM_COST;
                    unaryconstr.costs.push_back(tmpcost);
                }
                for (k = 0; k < ntuples; k++) {
                    file >> a;
                    file >> cost;
                    Cost tmpcost = cost*K;
                    if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb()) && getUb()<(MAX_COST / MEDIUM_COST)) tmpcost *= MEDIUM_COST;
                    assert(a>= 0 && a < x->getDomainInitSize());
                    unaryconstr.costs[a] = tmpcost;
                }
                if (shared) {
                    sharedSize.push_back(x->getDomainInitSize());
                    sharedCosts.push_back(unaryconstr.costs);
                    sharedTuples.push_back(emptyTuples);
                }
                unaryconstrs.push_back(unaryconstr);
            } else {
                file >> defval;
                if (defval == MIN_COST) {
                    cerr << "Error: unary cost function with zero penalty cost!" << endl;
                    exit(EXIT_FAILURE);
                }
                file >> ntuples;
                Value *dom = new Value[ntuples];
                for (k = 0; k < ntuples; k++) {
                    file >> dom[k];
                    file >> cost;
                    if (cost != MIN_COST) {
                        cerr << "Error: unary cost function with non-zero cost tuple!" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                postUnary(i,dom,ntuples,defval);
                delete [] dom;
            }
        } else if (arity == 0) {
            file >> defval;
            file >> ntuples;
            if (ToulBar2::verbose >= 3) cout << "read global lower bound contribution " << ic << " of " << defval << endl;
            if (ntuples > 1) {
                cerr << "Error: global lower bound contribution with several tuples!" << endl;
                exit(EXIT_FAILURE);
            }
            if (ntuples == 1) file >> cost;
            else cost = defval;
            inclowerbound += cost*K;
        }
    }

    file >> funcname;
    if (file) {
        cerr << "Warning: EOF not reached after reading all the cost functions (initial number of cost functions too small?)" << endl;
    }

    // merge unarycosts if they are on the same variable
    vector<int> seen(nbvar, -1);
    vector<TemporaryUnaryConstraint> newunaryconstrs;
    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        if (seen[unaryconstrs[u].var->wcspIndex] == -1) {
            seen[unaryconstrs[u].var->wcspIndex] = newunaryconstrs.size();
            newunaryconstrs.push_back(unaryconstrs[u]);
        }
        else {
            for (unsigned int i=0; i<unaryconstrs[u].var->getDomainInitSize(); i++) {
                if (newunaryconstrs[seen[unaryconstrs[u].var->wcspIndex]].costs[i] < getUb()) {
                    if (unaryconstrs[u].costs[i] < getUb()) newunaryconstrs[seen[unaryconstrs[u].var->wcspIndex]].costs[i] += unaryconstrs[u].costs[i];
                    else newunaryconstrs[seen[unaryconstrs[u].var->wcspIndex]].costs[i] = getUb();
                }
            }
        }
    }
    unaryconstrs = newunaryconstrs;
    if (ToulBar2::sortDomains) {
        if (maxarity > 2) {
            cout << "Warning! Cannot sort domains in preprocesing with non-binary cost functions." << endl;
        } else {
            ToulBar2::sortedDomains.clear();
            for (unsigned int u=0; u<unaryconstrs.size(); u++) {
                ToulBar2::sortedDomains[unaryconstrs[u].var->wcspIndex] = unaryconstrs[u].var->sortDomain(unaryconstrs[u].costs);
            }
        }
    }

    sortConstraints();
    // apply basic initial propagation AFTER complete network loading
    increaseLb(inclowerbound);

    // unary cost functions are delayed for compatibility issue (same lowerbound found) with old toolbar solver
    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        postUnary(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
    }
    histogram();

    cout << "Read " << nbvar << " variables, with " << nbvaltrue << " values at most, and " << nbconstr << " cost functions, with maximum arity " << maxarity  << "." << endl;
}


void WCSP::read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular )
{
    naryRandom randwcsp(this,seed);
    randwcsp.Input(n,m,p,forceSubModular);
    ToulBar2::nbvar= n;

    unsigned int nbconstr = numberOfConstraints();
    sortConstraints();
    histogram();

    if (ToulBar2::verbose >= 0) {
        cout << "Generated random problem " << n << " variables, with " << m << " values, and " << nbconstr << " cost functions." << endl;
    }
}




void WCSP::read_uai2008(const char *fileName)
{
    // Compute the factor that enables to capture the difference in log for probability (1-10^resolution):
    ToulBar2::NormFactor = (-1.0/Log1p(-Exp10(-(TLogProb)ToulBar2::resolution)));
    if (ToulBar2::NormFactor > (Pow((TProb)2., (TProb)INTEGERBITS)-1)/(TLogProb)ToulBar2::resolution) {
        cerr << "This resolution cannot be ensured on the data type used to represent costs." << endl;
        exit(EXIT_FAILURE);
    } else if (ToulBar2::verbose >= 1) {
        cout << "NormFactor= " << ToulBar2::NormFactor << endl;
    }

    // Cost inclowerbound = MIN_COST;
    string uaitype;
    ifstream file(fileName);
    if (!file) { cerr << "Could not open file " << fileName << endl; exit(EXIT_FAILURE); }

    Cost inclowerbound = MIN_COST;
    updateUb( (MAX_COST-UNIT_COST)/MEDIUM_COST/MEDIUM_COST/MEDIUM_COST/MEDIUM_COST );
    Cost upperbound = UNIT_COST;

    int nbval = 0;
    int nbvar,nbconstr;
    int i,j,k,ic;
    string varname;
    int domsize;
    EnumeratedVariable *x;
    EnumeratedVariable *y;
    EnumeratedVariable *z;
    unsigned int a;
    unsigned int b;
    unsigned int c;
    Cost cost;
    int ntuples;
    int arity;
    int maxarity = 0;
    vector<TemporaryUnaryConstraint> unaryconstrs;

    list<int> lctrs;

    file >> uaitype;

    if (ToulBar2::verbose >= 3) cout << "Reading " << uaitype << "  file." << endl;


    bool markov = uaitype == string("MARKOV");
    //bool bayes = uaitype == string("BAYES");


    file >> nbvar;
    ToulBar2::nbvar= nbvar;
    // read variable domain sizes
    for (i = 0; i < nbvar; i++) {
        string varname;
        varname = to_string(i);
        file >> domsize;
        if (ToulBar2::verbose >= 3) cout << "read variable " << i << " of size " << domsize << endl;
        if(domsize > nbval) nbval = domsize;
        int theindex = -1;
        if (domsize >= 0) theindex = makeEnumeratedVariable(varname,0,domsize-1);
        else theindex = makeIntervalVariable(varname,0,-domsize-1);
        assert(theindex == i);
    }


    file >> nbconstr;
    // read each constraint
    for (ic = 0; ic < nbconstr; ic++) {
        file >> arity;
        maxarity = max(maxarity,arity);

        if(arity > MAX_ARITY)  { cerr << "Nary cost functions of arity > " << MAX_ARITY << " not supported" << endl; exit(EXIT_FAILURE); }
        if (!file) {
            cerr << "Warning: EOF reached before reading all the cost functions (initial number of cost functions too large?)" << endl;
            break;
        }

        if (arity > 3) {
            int scopeIndex[MAX_ARITY];
            if (ToulBar2::verbose >= 3) cout << "read nary cost function on ";

            for(i=0;i<arity;i++) {
                file >> j;
                scopeIndex[i] = j;
                if (ToulBar2::verbose >= 3) cout << j << " ";
            }
            if (ToulBar2::verbose >= 3) cout << endl;
            lctrs.push_back( postNaryConstraintBegin(scopeIndex,arity,MIN_COST) );
            assert(lctrs.back() >= 0);
        }
        else if (arity == 3) {
            file >> i;
            file >> j;
            file >> k;
            if ((i == j) || (i == k) || (k == j)) {
                cerr << "Error: ternary cost function!" << endl;
                exit(EXIT_FAILURE);
            }
            x = (EnumeratedVariable *) vars[i];
            y = (EnumeratedVariable *) vars[j];
            z = (EnumeratedVariable *) vars[k];
            if (ToulBar2::verbose >= 3) cout << "read ternary cost function " << ic << " on " << i << "," << j << "," << k << endl;
            vector<Cost> costs;
            for (a = 0; a < x->getDomainInitSize(); a++) {
                for (b = 0; b < y->getDomainInitSize(); b++) {
                    for (c = 0; c < z->getDomainInitSize(); c++) {
                        costs.push_back(MIN_COST);
                    }
                }
            }
            lctrs.push_back( postTernaryConstraint(i,j,k,costs) );
            assert(lctrs.back() >= 0);
        }
        else if (arity == 2) {
            file >> i;
            file >> j;
            if (ToulBar2::verbose >= 3) cout << "read binary cost function " << ic << " on " << i << "," << j << endl;
            if (i == j) {
                cerr << "Error: binary cost function with only one variable in its scope!" << endl;
                exit(EXIT_FAILURE);
            }
            x = (EnumeratedVariable *) vars[i];
            y = (EnumeratedVariable *) vars[j];
            vector<Cost> costs;
            for (a = 0; a < x->getDomainInitSize(); a++) {
                for (b = 0; b < y->getDomainInitSize(); b++) {
                    costs.push_back(MIN_COST);
                }
            }
            lctrs.push_back( postBinaryConstraint(i,j,costs) );
            assert(lctrs.back() >= 0);
        }
        else if (arity == 1) {
            file >> i;
            if (ToulBar2::verbose >= 3) cout << "read unary cost function " << ic << " on " << i << endl;
            x = (EnumeratedVariable *) vars[i];
            TemporaryUnaryConstraint unaryconstr;
            unaryconstr.var = x;
            unaryconstrs.push_back(unaryconstr);
            lctrs.push_back(-1);
        } else if (arity == 0) {
            lctrs.push_back(-2);
        }
    }

    int iunaryctr = 0;
    int ictr = 0;
    Constraint* ctr = NULL;
    TernaryConstraint* tctr = NULL;
    BinaryConstraint* bctr = NULL;
    NaryConstraint* nctr = NULL;
    String s;

    ToulBar2::markov_log = 0;   // for the MARKOV Case

    int ntuplesarray[lctrs.size()];
    vector< vector<Cost> > costs;
    costs.resize(lctrs.size());
    list<int>::iterator it = lctrs.begin();
    while(it !=  lctrs.end()) {
        file >> ntuples;
        ntuplesarray[ictr] = ntuples;

        TProb p;
        vector<TProb> costsProb;

        TProb maxp = 0.;
        for (k = 0; k < ntuples; k++) {
            file >> p;
            assert(ToulBar2::uai>1 || (p >= 0. && (markov || p <= 1.)));
            costsProb.push_back( p );
            if(p > maxp) maxp = p;
        }
        if ( ToulBar2::uai==1 && maxp == 0. ) THROWCONTRADICTION;

        Cost minc = MAX_COST;
        Cost maxc = MIN_COST;
        for (k = 0; k < ntuples; k++) {
            p = costsProb[k];
            Cost cost;
            // ToulBar2::uai is 1 for .uai and 2 for .LG (log domain)
            if (markov) cost = ((ToulBar2::uai>1)?LogProb2Cost((TLogProb)(p - maxp)):Prob2Cost(p / maxp));
            else        cost = ((ToulBar2::uai>1)?LogProb2Cost((TLogProb)p):Prob2Cost(p));
            costs[ictr].push_back(cost);
            if(cost < minc) minc = cost;
            if(cost > maxc && cost < getUb()) maxc = cost;
        }
        upperbound += maxc;

        if(ToulBar2::preprocessNary>0 && minc > MIN_COST) {
            for (k = 0; k < ntuples; k++) {
                costs[ictr][k] -= minc;
            }
            if (ToulBar2::verbose >= 2) cout << "IC0 performed for cost function " << ictr << " with initial minimum cost " << minc << endl;
            inclowerbound += minc;
        }

        if(markov) ToulBar2::markov_log += ((ToulBar2::uai>1)?maxp:Log( maxp ));

        ictr++;
        ++it;
    }
    updateUb( upperbound );

    ictr = 0;
    it = lctrs.begin();
    while(it !=  lctrs.end()) {
        ntuples = ntuplesarray[ictr];
        for (k = 0; k < ntuples; k++) {
            if(CUT(costs[ictr][k], getUb())) costs[ictr][k] = getUb() * MEDIUM_COST;
        }

        int arity;
        if(*it == -1) { ctr = NULL; arity = 1; }
        else if(*it == -2) { ctr = NULL; arity = 0; }
        else { assert(*it >= 0); ctr = getCtr(*it); arity = ctr->arity(); }
        switch(arity) {
        case 0: {
            inclowerbound += costs[ictr][0];
            break; }

        case 1: {
            unaryconstrs[iunaryctr].costs.clear();
            for (a = 0; a < unaryconstrs[iunaryctr].var->getDomainInitSize(); a++) {
                unaryconstrs[iunaryctr].costs.push_back(costs[ictr][a]);
            }
            iunaryctr++;
            if (ToulBar2::verbose >= 3) cout << "read unary costs."  << endl;
            break; }

        case 2: {
            bctr = (BinaryConstraint*) ctr;
            x = (EnumeratedVariable*) bctr->getVar(0);
            y = (EnumeratedVariable*) bctr->getVar(1);
            postBinaryConstraint( x->wcspIndex, y->wcspIndex, costs[ictr] );
            if (ToulBar2::verbose >= 3) cout << "read binary costs."  << endl;
            break; }

        case 3: {
            tctr = (TernaryConstraint*) ctr;
            x = (EnumeratedVariable*) tctr->getVar(0);
            y = (EnumeratedVariable*) tctr->getVar(1);
            z = (EnumeratedVariable*) tctr->getVar(2);
            postTernaryConstraint( x->wcspIndex, y->wcspIndex, z->wcspIndex, costs[ictr] );
            if (ToulBar2::verbose >= 3) cout << "read ternary costs." << endl;
            break; }

        default: {
            nctr = (NaryConstraint*) ctr;
            j = 0;
            nctr->firstlex();
            while(nctr->nextlex(s,cost)) {
                //					  if (costs[j]>MIN_COST) nctr->setTuple(s, costs[j]);
                postNaryConstraintTuple(nctr->wcspIndex, s, costs[ictr][j]);
                j++;
            }
            if (ToulBar2::verbose >= 3) cout << "read arity " << arity << " table costs."  << endl;
            postNaryConstraintEnd(nctr->wcspIndex);
            break; }
        }
        ictr++;
        ++it;
    }
    if (ToulBar2::verbose >= 1) {
        cout << "MarkovShiftingValue= " << ToulBar2::markov_log << endl;
    }

    sortConstraints();
    // apply basic initial propagation AFTER complete network loading
    increaseLb(inclowerbound);

    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        postUnary(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
    }
    histogram();
    cout << "Read " << nbvar << " variables, with " << nbval << " values at most, and " << nbconstr << " cost functions, with maximum arity " << maxarity  << "." << endl;

    int nevi = 0;
    ifstream fevid(ToulBar2::evidence_file.c_str());
    if (!fevid)	{
        string strevid(string(fileName) + string(".evid"));
        fevid.open(strevid.c_str());
        cerr << "No evidence file specified. Trying " << strevid << endl;
        if(!fevid) cerr << "No evidence file. " << endl;
    }
    if(fevid) {
        vector<int> variables;
        vector<Value> values;
        fevid >> nevi;
        bool firstevid =true;
        if (nevi == 0) return;
        if (nevi == 1) fevid >> nevi;  // UAI 2010 evidence file format assumes possible multiple evidence samples, but toulbar2 will search for the first evidence sample only!
        while(nevi) {
            if(!fevid) {
                cerr << "Error: incorrect number of evidences." << endl;
                exit(EXIT_FAILURE);
            }
            fevid >> i;
            fevid >> j;
            if (firstevid && !fevid) { // old UAI 2008 evidence format
                variables.push_back(nevi);
                values.push_back(i);
                break;
            } else firstevid = false;
            variables.push_back(i);
            values.push_back(j);
            nevi--;
        }
        assignLS(variables, values);
    }
}


void WCSP::solution_UAI(Cost res, bool opt)
{
    if (!ToulBar2::uai && !ToulBar2::uaieval) return;
    if (ToulBar2::isZ) return;
    // UAI 2012 Challenge output format
    if (ToulBar2::uai_firstoutput && !ToulBar2::uaieval) ToulBar2::uai_firstoutput = false;
    else {
        //	    ToulBar2::solution_file << "-BEGIN-" << endl;
        ToulBar2::solution_file.close();
        ToulBar2::solution_file.open(ToulBar2::solution_uai_filename.c_str());
        ToulBar2::solution_file << "MPE" << endl;
    }
    //	ToulBar2::solution_file << "1" << endl; // we assume a single evidence sample
    if (ToulBar2::showSolutions && !ToulBar2::uaieval) {
        cout << "t " << cpuTime() - ToulBar2::startCpuTime << endl;
        cout << "s " << (Cost2LogProb(res) + ToulBar2::markov_log)/Log(10.) << endl;
        cout << numberOfVariables();
        printSolution(cout);
    }
    ToulBar2::solution_file << numberOfVariables();
    printSolution(ToulBar2::solution_file);
    ToulBar2::solution_file << endl;
    //	if (opt) {
    //	  if (ToulBar2::showSolutions) cout << " LU" << endl;
    //	  ToulBar2::solution_file << " LU" << endl;
    //	} else {
    //	  if (ToulBar2::showSolutions) cout << " L" << endl;
    //	  ToulBar2::solution_file << " L" << endl;
    //	}
    ToulBar2::solution_file.flush();
}



#ifdef XMLFLAG
#include "./xmlcsp/xmlcsp.h"
#endif


void WCSP::read_XML(const char *fileName)
{
#ifdef XMLFLAG
    MyCallback xmlCallBack;
    xmlCallBack.wcsp = this;
    xmlCallBack.fname = string(fileName);
    xmlCallBack.convertWCSP = true;
    try {
        XMLParser_libxml2<> parser( xmlCallBack );
        parser.setPreferredExpressionRepresentation(INFIX_C);
        parser.parse(fileName);
    } catch (exception &e) {
        cout.flush();
        cerr << "\n\tUnexpected exception in XML parsing\n";
        cerr << "\t" << e.what() << endl;
        exit(1);
    }
#else
    cerr << "\nXML format without including in Makefile flag XMLFLAG and files ./xmlcsp\n" << endl;
    exit(1);
#endif
}


void WCSP::solution_XML(bool opt)
{
#ifdef XMLFLAG
    if (!ToulBar2::xmlflag) return;

    if(opt)  cout << "s OPTIMUM FOUND" << endl;

    //ofstream fsol;
    ifstream sol;
    sol.open("sol");
    //if(!sol) { cout << "cannot open solution file to translate" << endl; exit(1); }
    //fsol.open("solution");
    //fsol << "SOL ";


    cout << "v ";
    for (unsigned int i=0; i<vars.size(); i++) {
        int value;
        sol >> value;
        int index = ((EnumeratedVariable*) getVar(i))->toIndex(value);
        cout << Doms[varsDom[i]][ index ] << " ";
    }
    cout << endl;

    //fsol << endl;
    //fsol.close();
    sol.close();
#endif
}

void WCSP::read_wcnf(const char *fileName)
{
    ifstream file(fileName);
    if (!file) { cerr << "Could not open file " << fileName << endl; exit(EXIT_FAILURE); }

    double K = ToulBar2::costMultiplier;
    Cost inclowerbound = MIN_COST;
    updateUb( (MAX_COST-UNIT_COST)/MEDIUM_COST/MEDIUM_COST );

    int maxarity = 0;
    vector<TemporaryUnaryConstraint> unaryconstrs;

    int nbvar,nbclauses;
    string dummy,sflag;

    file >> sflag;
    while (sflag == "c") {
        getline( file, dummy );
        file >> sflag;
    }
    if (sflag != "p") {
        cerr << "Wrong wcnf format in " << fileName << endl;
        exit(EXIT_FAILURE);
    }

    string format,strtop;
    Cost top;
    file >> format;
    file >> nbvar;
    ToulBar2::nbvar= nbvar;
    file >> nbclauses;
    if (format == "wcnf") {
        getline( file, strtop );
        if (string2Cost((char*) strtop.c_str())>0) {
            cout << "c (Weighted) Partial Max-SAT input format" << endl;
            top = string2Cost((char*) strtop.c_str());
            if(top < MAX_COST / K)	top = top * K;
            else top = MAX_COST;
            updateUb(top);
        } else {
            cout << "c Weighted Max-SAT input format" << endl;
        }
    } else {
        cout << "c Max-SAT input format" << endl;
        updateUb((nbclauses+1)*K);
    }

    // create Boolean variables
    for (int i = 0; i<nbvar; i++) {
        string varname;
        varname = to_string(i);
        int theindex = -1;
        theindex = makeEnumeratedVariable(varname,0,1);
        assert(theindex == i);
    }

    // Read each clause
    for (int ic = 0; ic < nbclauses; ic++) {

        int scopeIndex[MAX_ARITY];
        Char buf[MAX_ARITY];
        int arity = 0;
        if (ToulBar2::verbose >= 3) cout << "read clause on ";

        int j = 0;
        Cost cost = UNIT_COST;
        if (format == "wcnf") file >> cost;
        bool tautology = false;
        do {
            file >> j;
            if (j != 0 && !tautology) {
                scopeIndex[arity] = abs(j) - 1;
                buf[arity] = ((j>0)?0:1) + CHAR_FIRST;
                int k = 0;
                while (k<arity) {
                    if (scopeIndex[k]==scopeIndex[arity]) {
                        break;
                    }
                    k++;
                }
                if (k<arity) {
                    if (buf[k]!=buf[arity]) {
                        tautology = true;
                        if (ToulBar2::verbose >= 3) cout << j << " is a tautology! skipped.";
                    }
                    continue;
                }
                arity++;
                if (ToulBar2::verbose >= 3) cout << j << " ";
            }
        } while (j != 0);
        if (ToulBar2::verbose >= 3) cout << endl;
        if (tautology) continue;
        buf[arity] = '\0';
        maxarity = max(maxarity,arity);

        if (arity > 3) {
            int index = postNaryConstraintBegin(scopeIndex,arity,MIN_COST);
            String tup = buf;
            postNaryConstraintTuple(index, tup, cost*K);
            postNaryConstraintEnd(index);
        } else if (arity == 3) {
            vector<Cost> costs;
            for (int a = 0; a < 2; a++) {
                for (int b = 0; b < 2; b++) {
                    for (int c = 0; c < 2; c++) {
                        costs.push_back(MIN_COST);
                    }
                }
            }
            costs[(buf[0] - CHAR_FIRST)*4 + (buf[1] - CHAR_FIRST)*2 + (buf[2] - CHAR_FIRST)] = cost*K;
            postTernaryConstraint(scopeIndex[0], scopeIndex[1], scopeIndex[2], costs);
        } else if (arity == 2) {
            vector<Cost> costs;
            for (int a = 0; a < 2; a++) {
                for (int b = 0; b < 2; b++) {
                    costs.push_back(MIN_COST);
                }
            }
            costs[(buf[0] - CHAR_FIRST)*2 + (buf[1] - CHAR_FIRST)] = cost*K;
            postBinaryConstraint(scopeIndex[0], scopeIndex[1], costs);
        } else if (arity == 1) {
            EnumeratedVariable *x = (EnumeratedVariable *) vars[scopeIndex[0]];
            TemporaryUnaryConstraint unaryconstr;
            unaryconstr.var = x;
            if ((buf[0] - CHAR_FIRST)==0) {
                unaryconstr.costs.push_back(cost*K);
                unaryconstr.costs.push_back(MIN_COST);
            } else {
                unaryconstr.costs.push_back(MIN_COST);
                unaryconstr.costs.push_back(cost*K);
            }
            unaryconstrs.push_back(unaryconstr);
        } else if (arity == 0) {
            inclowerbound += cost*K;
        } else {
            cerr << "Wrong clause arity " << arity << " in " << fileName << endl;
            exit(EXIT_FAILURE);
        }
    }

    file >> dummy;
    if (file) {
        cerr << "Warning: EOF not reached after reading all the clauses (initial number of clauses too small?)" << endl;
    }

    sortConstraints();
    // apply basic initial propagation AFTER complete network loading
    increaseLb(inclowerbound);

    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        postUnary(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
    }
    histogram();
    cout << "c Read " << nbvar << " variables, with 2 values at most, and " << nbclauses << " clauses, with maximum arity " << maxarity  << "." << endl;
}

/// \brief minimizes/maximizes \f$ X^t \times W \times X = \sum_{i=1}^N \sum_{j=1}^N W_{ij} \times X_i \times X_j \f$
/// where W is expressed by its M non-zero half squared matrix costs (can be positive or negative float numbers)
/// \note Costs for \f$ i \neq j \f$ are multiplied by 2 by this method (symmetric N*N squared matrix)
/// \note If N is positive, then variable domain values are {0,1}
/// \note If N is negative, then variable domain values are {1,-1} with value 1 having index 0 and value -1 having index 1 in the output solutions
/// \note If M is positive then minimizes the quadratic function, else maximizes it
/// \warning It does not allow infinite costs (no forbidden assignments)
void WCSP::read_qpbo(const char *fileName)
{
    ifstream file(fileName);
    if (!file) { cerr << "Could not open file " << fileName << endl; exit(EXIT_FAILURE); }

    int n = 0;
    file >> n;
    bool booldom = (n >= 0); // n positive means variable domains {0,1} else {1,-1}
    if (!booldom) n = -n;
    int m = 0;
    file >> m;
    if (n == 0 || m == 0) return;
    bool minimize = (m >= 0); // m positive means minimize the quadratic function, else maximize it
    if (!minimize) m = -m;
    int e = 0;
    int dummy;

    vector<int> posx(m, 0);
    vector<int> posy(m, 0);
    vector<double> cost(m, 0.);
    for (e=0; e<m; e++) {
        file >> posx[e];
        if (!file) {
            cerr << "Warning: EOF reached before reading all the cost sparse matrix (number of nonzero costs too large?)" << endl;
            break;
        }
        if (posx[e]>n) {
            cerr << "Warning: variable index too large!" << endl;
            break;
        }
        file >> posy[e];
        if (posy[e]>n) {
            cerr << "Warning: variable index too large!" << endl;
            break;
        }
        file >> cost[e];
    }
    file >> dummy;
    if (file) {
        cerr << "Warning: EOF not reached after reading all the cost sparse matrix (wrong number of nonzero costs too small?)" << endl;
    }
    m = e;

    // create Boolean variables
    ToulBar2::nbvar= n;
    for (int i=0; i<n; i++) {
        makeEnumeratedVariable(to_string(i), 0, 1);
    }

    vector<Cost> unaryCosts0(n, 0);
    vector<Cost> unaryCosts1(n, 0);

    // find total cost
    Double sumcost = 0.;
    for (int e=0; e<m; e++) {
        sumcost += 2. * abs(cost[e]);
    }
    Double multiplier = Exp10((Double) ToulBar2::resolution);
    if (multiplier * sumcost >= (Double) MAX_COST) {
        cerr << "This resolution cannot be ensured on the data type used to represent costs! (see option -precision)" << endl;
        exit(EXIT_FAILURE);
    }
    updateUb((Cost) multiplier * sumcost  +1);

    // create weighted binary clauses
    for (int e=0; e<m; e++) {
        if (posx[e] != posy[e]) {
            vector<Cost> costs(4, 0);
            if (booldom) {
                if (cost[e] > 0) {
                    if (minimize) {
                        costs[3] = (Cost) (multiplier * 2. * cost[e]);
                    } else {
                        costs[0] = (Cost) (multiplier * 2. * cost[e]);
                        costs[1] = costs[0];
                        costs[2] = costs[0];
                    }
                } else {
                    if (minimize) {
                        costs[0] = (Cost) (multiplier * -2. * cost[e]);
                        costs[1] = costs[0];
                        costs[2] = costs[0];
                    } else {
                        costs[3] = (Cost) (multiplier * -2. * cost[e]);
                    }
                }
            } else {
                if (cost[e] > 0) {
                    if (minimize) {
                        costs[0] = (Cost) (multiplier * 2. * cost[e]);
                        costs[3] = costs[0];
                    } else {
                        costs[1] = (Cost) (multiplier * 2. * cost[e]);
                        costs[2] = costs[1];
                    }
                } else {
                    if (minimize) {
                        costs[1] = (Cost) (multiplier * -2. * cost[e]);
                        costs[2] = costs[1];
                    } else {
                        costs[0] = (Cost) (multiplier * -2. * cost[e]);
                        costs[3] = costs[0];
                    }
                }
            }
            postBinaryConstraint(posx[e] - 1, posy[e] - 1, costs);
        } else {
            if (booldom) {
                if (cost[e] > 0) {
                    if (minimize) {
                        unaryCosts1[posx[e] - 1] += (Cost) (multiplier * cost[e]);
                    } else {
                        unaryCosts0[posx[e] - 1] += (Cost) (multiplier * cost[e]);
                    }
                } else {
                    if (minimize) {
                        unaryCosts0[posx[e] - 1] += (Cost) (multiplier * -cost[e]);
                    } else {
                        unaryCosts1[posx[e] - 1] += (Cost) (multiplier * -cost[e]);
                    }
                }
            } else {
                if (cost[e] > 0) {
                    if (minimize) {
                        unaryCosts0[posx[e] - 1] += (Cost) (multiplier * cost[e]);
                    } else {
                        unaryCosts1[posx[e] - 1] += (Cost) (multiplier * cost[e]);
                    }
                } else {
                    if (minimize) {
                        unaryCosts1[posx[e] - 1] += (Cost) (multiplier * -cost[e]);
                    } else {
                        unaryCosts0[posx[e] - 1] += (Cost) (multiplier * -cost[e]);
                    }
                }
            }
        }
    }

    sortConstraints();

    // create weighted unary clauses
    for (int i=0; i<n; i++) {
        if (unaryCosts0[i] > 0 || unaryCosts1[i] > 0) {
            vector<Cost> costs(2, 0);
            costs[0] = unaryCosts0[i];
            costs[1] = unaryCosts1[i];
            postUnary(i, costs);
        }
    }
    histogram();
    cout << "Read " << n << " variables, with " << 2 << " values at most, and " << m << " nonzero matrix costs." << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

