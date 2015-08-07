#include "Gecode.hpp"
#include <limits>

/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

// typedef enum {NOTVAR, BOOLVAR, INTVARNVAL, INTVARBOUNDS, INTVARVALUES} VariableType;


Gecode_Expression::Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a Boolean expression" << std::endl;
#endif

    _solver = NULL;
    _gcvarid = -1;
    _vartype = BOOLVAR;
}

Gecode_Expression::Gecode_Expression(const int nval) {

#ifdef _DEBUGWRAP
    std::cout << "creating a variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif

    _solver = NULL;
    _gcvarid = -1;
    _vartype = INTVARNVAL;
    _nval = nval;
}

Gecode_Expression::Gecode_Expression(const int lb, const int ub) {

#ifdef _DEBUGWRAP
    std::cout << "creating a variable [" << lb << ".." << ub << "]" << std::endl;
#endif

    _solver = NULL;
    _gcvarid = -1;
    _vartype = INTVARBOUNDS;
    _lb = lb; _ub = ub;
}

Gecode_Expression::Gecode_Expression(GecodeIntArray& vals) {

#ifdef _DEBUGWRAP
    std::cout << "creating a domain variable" << std::endl;
#endif

    _solver = NULL;
    _gcvarid = -1;
    _vartype = INTVARVALUES;
    for(int i=0; i<vals.size(); i++) {
        _vals.add(vals.get_item(i));
    }
}


int Gecode_Expression::getVariableId() const {

#ifdef _DEBUGWRAP
    std::cout << "return identity of expression" << std::endl;
#endif

    return _gcvarid;
}

int Gecode_Expression::get_value() const {

// #ifdef _DEBUGWRAP
//     std::cout << "return value of expression" << std::endl;
// #endif

    if(_solver != NULL && _gcvarid >= 0){
        // std::cout << "get value for variable id " << _gcvarid << std::endl;
        try {
            // std::cout << "var:" << _solver->gecodespace->getVar(_gcvarid) << std::endl;
            return _solver->gecodespace->getVar(_gcvarid).val();
        } catch(Gecode::Int::ValOfUnassignedVar e){
            std::cout << "ERROR asked for value of unassigned var. FIXME" << std::endl;
        }
    }
    return -1;
}

int Gecode_Expression::get_size() const {

#ifdef _DEBUGWRAP
    std::cout << "return size of expression" << std::endl;
#endif

    // FIXME add support for if the variable has been loaded, to check its value
    // in Gecode

    if(_vartype == BOOLVAR){
        return 2;
    } else if(_vartype == INTVARNVAL){
        return _nval;
    } else if(_vartype == INTVARBOUNDS){
        return _ub - _lb + 1;
    } else if(_vartype == INTVARVALUES){
        return _vals.size();
    }
    return 0;
}

int Gecode_Expression::get_max() const {

#ifdef _DEBUGWRAP
    std::cout << "return max of expression" << std::endl;
#endif
    if(_vartype == BOOLVAR){
        return 1;
    } else if(_vartype == INTVARNVAL){
        return _nval-1;
    } else if(_vartype == INTVARBOUNDS){
        return _ub;
    } else if(_vartype == INTVARVALUES){
        if(_vals.size()) return _vals.get_item(_vals.size()-1);
    }
    return 0;
}

int Gecode_Expression::get_min() const {

#ifdef _DEBUGWRAP
    std::cout << "return min of expression" << std::endl;
#endif
    if(_vartype == BOOLVAR){
        return 0;
    } else if(_vartype == INTVARNVAL){
        return 0;
    } else if(_vartype == INTVARBOUNDS){
        return _lb;
    } else if(_vartype == INTVARVALUES){
        if(_vals.size()) return _vals.get_item(0);
    }
    return 0;
}

bool Gecode_Expression::contain(const int v) const {

#ifdef _DEBUGWRAP
    std::cout << "check expression contains " << v << std::endl;
#endif
    // return _self.contain(v);
    return 0;
}

Gecode_Expression::~Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "delete expression" << std::endl;
#endif
}

bool Gecode_Expression::has_been_added() const {
    return _solver != NULL;
}

Gecode::IntVar Gecode_Expression::getGecodeVar() const {
    assert(_solver != NULL);
    if(_gcintrepr.size() > 0) {
        // std::cout << "Returning _gcintrepr" << _gcintrepr[0] << std::endl;
        return _gcintrepr[0];
    }
    assert(_gcvarid >= 0);
    return _solver->gecodespace->getVar(_gcvarid);
}

Gecode_Expression* Gecode_Expression::add(GecodeSolver *solver, bool top_level) {

#ifdef _DEBUGWRAP
    std::cout << "add expression to model" << std::endl;
    if(top_level) std::cout << "\tAdding at top level" << std::endl;
    else std::cout << "\tAdding within tree" <<std::endl;
#endif

    if(!has_been_added()) {
        _solver = solver;
        if(_vartype == BOOLVAR){
            _gcvarid = solver->gecodespace->createIntVar(0, 1);
        } else if(_vartype == INTVARNVAL){
            _gcvarid = solver->gecodespace->createIntVar(0, _nval-1);
        } else if(_vartype == INTVARBOUNDS){
            _gcvarid = solver->gecodespace->createIntVar(_lb, _ub);
        } else if(_vartype == INTVARVALUES){
            _gcvarid = solver->gecodespace->createIntVar(_vals);
        } else {
            std::cerr << "Error trying to create variable, unknown variable type " << _vartype << std::endl;
            exit(1);
        }

        if(top_level) {
            _solver->add(this);
        }
    }

    return this;
}

// /* Binary operators */

Gecode_binop::Gecode_binop(Gecode_Expression *var1,
                           Gecode_Expression *var2)
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a binary operator" << std::endl;
#endif
    _vars[0] = var1;
    _vars[1] = var2;
}

Gecode_binop::Gecode_binop(Gecode_Expression *var1, int constant)
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a binary (constant) operator" << std::endl;
#endif
    _vars[0] = var1;
    _vars[1] = NULL;
    _constant = constant;
}

Gecode_binop::Gecode_binop(int(*opfunc)(int, int), Gecode_Expression *var1, Gecode_Expression *var2) :
        binaryopfunc(opfunc){
    _vars[0] = var1;
    _vars[1] = var2;
    initialise();
}

Gecode_binop::Gecode_binop(int(*opfunc)(int, int), Gecode_Expression *var1, int constant) :
        binaryopfunc(opfunc) {
    _vars[0] = var1;
    _vars[1] = NULL;
    _constant = constant;
    initialise();
}

Gecode_binop::~Gecode_binop() {

#ifdef _DEBUGWRAP
    std::cout << "delete binary operator" << std::endl;
#endif
}

void Gecode_binop::initialise() {
    int lb1 = _vars[0]->get_min(), ub1 = _vars[0]->get_max(), lb2, ub2;
    if(_vars[1] != NULL){
        lb2 = _vars[1]->get_min();
        ub2 = _vars[1]->get_max();
    } else {
        lb2 = ub2 = _constant;
    }
    std::vector<int> bounds;
    bounds.push_back(binaryopfunc(lb1, lb2));
    bounds.push_back(binaryopfunc(lb1, ub2));
    bounds.push_back(binaryopfunc(ub1, lb2));
    bounds.push_back(binaryopfunc(ub1, ub2));
    _lb = *(std::min_element(bounds.begin(), bounds.end()));
    _ub = *(std::max_element(bounds.begin(), bounds.end()));
#ifdef _DEBUGWRAP
    std::cout << "Initialised bounds of binary operator to " << _lb << ".." << _ub << std::endl;
#endif
}

/**
 * Constraints
 */

Gecode_Min::Gecode_Min( GecodeExpArray& vars )
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating an Min constraint" << std::endl;
#endif

    _vars = vars;
    initialise();
}

Gecode_Min::Gecode_Min( Gecode_Expression *var1, Gecode_Expression *var2 )
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a binary Min constraint" << std::endl;
#endif

    _vars.add(var1);
    _vars.add(var2);
    initialise();
}

void Gecode_Min::initialise(){
    int vmin, vmax;
    _lb = std::numeric_limits<int>::max();
    _ub = std::numeric_limits<int>::max();
    for(unsigned int i=0; i<_vars.size(); ++i) {
        vmin = _vars.get_item(i)->get_min();
        vmax = _vars.get_item(i)->get_max();

        if(vmin < _lb) _lb = vmin;
        if(vmax < _ub) _ub = vmax;
    }
#ifdef _DEBUGWRAP
    std::cout << "init bounds of minimum to " << _lb << ".." << _ub << std::endl;
#endif
}

Gecode_Min::~Gecode_Min() {

#ifdef _DEBUGWRAP
    std::cout << "delete Min" << std::endl;
#endif
}

Gecode_Expression* Gecode_Min::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {

#ifdef _DEBUGWRAP
        std::cout << "add a Min constraint to solver" << std::endl;
#endif
        if(top_level){
            std::cerr << "Error: minimum expression not valid at the top-level, it should be used as a sub-expression." << std::endl;
            exit(1);
        }

        _solver = solver;

        int i, n=_vars.size();
        Gecode::IntVarArgs scope(n);
        for(i=0; i<n; ++i){
            _vars.set_item(i, _vars.get_item(i)->add(_solver, false));
            scope[i] = _vars.get_item(i)->getGecodeVar();
        }

        _gcintrepr << expr(*(solver->gecodespace), Gecode::min(scope));
    }
    return this;
}


Gecode_Max::Gecode_Max( GecodeExpArray& vars )
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a Max constraint" << std::endl;
#endif

    _vars = vars;
    initialise();
}

Gecode_Max::Gecode_Max( Gecode_Expression *var1, Gecode_Expression *var2 )
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a binary Max constraint" << std::endl;
#endif

    _vars.add(var1);
    _vars.add(var2);
    initialise();
}

void Gecode_Max::initialise(){
    int vmin, vmax;
    _lb = std::numeric_limits<int>::min();
    _ub = std::numeric_limits<int>::min();
    for(unsigned int i=0; i<_vars.size(); ++i) {
        vmin = _vars.get_item(i)->get_min();
        vmax = _vars.get_item(i)->get_max();

        if(vmin > _lb) _lb = vmin;
        if(vmax > _ub) _ub = vmax;
    }
#ifdef _DEBUGWRAP
    std::cout << "init bounds of maximum to " << _lb << ".." << _ub << std::endl;
#endif
}

Gecode_Max::~Gecode_Max() {

#ifdef _DEBUGWRAP
    std::cout << "delete Max" << std::endl;
#endif
}

Gecode_Expression* Gecode_Max::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {

#ifdef _DEBUGWRAP
        std::cout << "add an Max constraint" << std::endl;
#endif

        if(top_level){
            std::cerr << "Error: maximum expression not valid at the top-level, it should be used as a sub-expression." << std::endl;
            exit(1);
        }
        _solver = solver;

        int i, n=_vars.size();
        Gecode::IntVarArgs scope(n);
        for(i=0; i<n; ++i){
            _vars.set_item(i, _vars.get_item(i)->add(_solver, false));
            scope[i] = _vars.get_item(i)->getGecodeVar();
        }

        _gcintrepr << expr(*(solver->gecodespace), Gecode::max(scope));
    }
    return this;
}

Gecode_AllDiff::Gecode_AllDiff( GecodeExpArray& vars )
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating an alldiff constraint" << std::endl;
#endif
    _vars = vars;
}

Gecode_AllDiff::Gecode_AllDiff( Gecode_Expression *var1, Gecode_Expression *var2 )
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating a binary alldiff constraint" << std::endl;
#endif
    _vars.add(var1);
    _vars.add(var2);
}

Gecode_AllDiff::~Gecode_AllDiff() {
#ifdef _DEBUGWRAP
    std::cout << "delete alldiff" << std::endl;
#endif
}

Gecode_Expression* Gecode_AllDiff::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add an alldiff constraint" << std::endl;
#endif

        if(!top_level){
            std::cerr << "AllDiff sub-expression not supported by this solver." << std::endl;
            exit(1);
        }

        _solver = solver;
        int n = _vars.size();
        Gecode::IntVarArgs gcvars(n);
        for(int i=0; i<n; ++i){
            _vars.set_item(i, _vars.get_item(i)->add(solver, false));
            gcvars[i] = _vars.get_item(i)->getGecodeVar();
        }

        if(n==2) {
            Gecode::rel(*(solver->gecodespace), _vars.get_item(0)->getGecodeVar() != _vars.get_item(1)->getGecodeVar());
        } else {
            Gecode::distinct(*(solver->gecodespace), gcvars);
        }
    }
    return this;
}

Gecode_Gcc::Gecode_Gcc(GecodeExpArray& vars,
                       GecodeIntArray& vals,
                       GecodeIntArray& lb_card,
                       GecodeIntArray& ub_card)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating a gcc constraint" << std::endl;
#endif

    _vars = vars;
    _vals = vals;
    _lb_card = lb_card;
    _ub_card = ub_card;

}

Gecode_Gcc::~Gecode_Gcc() {
#ifdef _DEBUGWRAP
    std::cout << "delete gcc" << std::endl;
#endif
}

Gecode_Expression* Gecode_Gcc::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add gcc constraint" << std::endl;
#endif

        _solver = solver;
        int i, n=_vars.size(), m=_vals.size();
        Gecode::IntVarArgs scope(n);
        Gecode::IntSetArgs c(m);
        Gecode::IntArgs v(m);

        for(i=0; i<n; i++){
            _vars.set_item(i, _vars.get_item(i)->add(solver, false));
            scope[i] = _vars.get_item(i)->getGecodeVar();
        }

        for(i=0; i<m; i++){
            c[i] = Gecode::IntSet(_lb_card.get_item(i), _ub_card.get_item(i));
            v[i] = _vals.get_item(i);
        }

        count(*(solver->gecodespace), scope, c, v);
    }

    return this;
}

Gecode_Element::Gecode_Element( GecodeExpArray& vars )
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating element" << std::endl;
#endif

    _vars = vars;
    initialise();
}

Gecode_Element::~Gecode_Element() {
#ifdef _DEBUGWRAP
    std::cout << "delete element" << std::endl;
#endif
}

void Gecode_Element::initialise() {
    int indexlb, indexub, vmin, vmax, n = _vars.size();
    indexlb = _vars.get_item(n-1)->get_min();
    indexub = _vars.get_item(n-1)->get_max();
    indexlb = indexlb >= 0 ? indexlb : 0;
    indexub = indexub < n-1 ? indexub : n-1;

    _lb = std::numeric_limits<int>::max();
    _ub = std::numeric_limits<int>::min();
    for(unsigned int i=indexlb; i<indexub; ++i) {
        vmin = _vars.get_item(i)->get_min();
        vmax = _vars.get_item(i)->get_max();

        if(vmin < _lb) _lb = vmin;
        if(vmax > _ub) _ub = vmax;
    }
    std::cout << "init bounds of element to " << _lb << ".." << _ub << std::endl;
}

Gecode_Expression* Gecode_Element::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add element constraint" << std::endl;
#endif

        if(top_level) {
            std::cerr << "Error: element expression not valid at the top-level, it should be used as a sub-expression." << std::endl;
            exit(1);
        }
        _solver = solver;
        int n = _vars.size();
        Gecode::IntVarArgs shared(n - 1);
        for(int i=0; i<n; ++i){
            _vars.set_item(i, _vars.get_item(i)->add(solver, false));
            if(i<n-1) shared[i] = _vars.get_item(i)->getGecodeVar();
        }

        Gecode_Expression *aux = new Gecode_Expression(_lb, _ub);
        aux = aux->add(solver, false);
        element(*(solver->gecodespace), shared, _vars.get_item(n-1)->getGecodeVar(), aux->getGecodeVar());
        return aux;
    }

    return this;
}

Gecode_LeqLex::Gecode_LeqLex( GecodeExpArray& vars )
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating lexleq" << std::endl;
#endif
    _vars = vars;
}

Gecode_LeqLex::~Gecode_LeqLex() {
#ifdef _DEBUGWRAP
    std::cout << "delete leqlex" << std::endl;
#endif
}

Gecode_Expression* Gecode_LeqLex::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add leqlex constraint" << std::endl;
#endif
        _solver = solver;
        int n = _vars.size();
        Gecode::IntVarArgs left(_vars.size() / 2), right(_vars.size() / 2);
        for(int i=0; i<n; ++i){
            _vars.set_item(i, _vars.get_item(i)->add(solver, false));

            if(i<n / 2) left[i] = _vars.get_item(i)->getGecodeVar();
            else right[i - n / 2] = _vars.get_item(i)->getGecodeVar();
        }

        if(top_level) {
            lex(*(solver->gecodespace), left, Gecode::IRT_LQ, right);
        } else {
            std::cerr << "Error: reified lexicographic constraint not supported with this solver." << std::endl;
            exit(1);
        }
    }
    return this;
}


Gecode_LessLex::Gecode_LessLex( GecodeExpArray& vars )
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating lexless" << std::endl;
#endif
    _vars = vars;
}

Gecode_LessLex::~Gecode_LessLex() {
#ifdef _DEBUGWRAP
    std::cout << "delete leslex" << std::endl;
#endif
}

Gecode_Expression* Gecode_LessLex::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add lesslex constraint" << std::endl;
#endif
        _solver = solver;
        int n = _vars.size();
        Gecode::IntVarArgs left(_vars.size() / 2), right(_vars.size() / 2);
        for(int i=0; i<n; ++i){
            _vars.set_item(i, _vars.get_item(i)->add(solver, false));

            if(i<n / 2) left[i] = _vars.get_item(i)->getGecodeVar();
            else right[i - n / 2] = _vars.get_item(i)->getGecodeVar();
        }

        if(top_level) {
            lex(*(solver->gecodespace), left, Gecode::IRT_LE, right);
        } else {
            std::cerr << "Error: reified lexicographic constraint not supported with this solver." << std::endl;
            exit(1);
        }
    }
    return this;
}


Gecode_Sum::Gecode_Sum(GecodeExpArray& vars,
                       GecodeIntArray& weights,
                       const int offset)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating sum A" << std::endl;
#endif
    _vars = vars;
    _weights = weights;
    _offset = offset;
    initialise();
}

Gecode_Sum::Gecode_Sum(Gecode_Expression *arg1,
                       Gecode_Expression *arg2,
                       GecodeIntArray& w,
                       const int offset)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating sum B" << std::endl;
#endif
    _vars.add(arg1);
    _vars.add(arg2);
    _weights = w;
    _offset = offset;
    initialise();
}

Gecode_Sum::Gecode_Sum(Gecode_Expression *arg,
                       GecodeIntArray& w,
                       const int offset)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating sum C" << std::endl;
#endif
    _vars.add(arg);
    _weights = w;
    _offset = offset;
    initialise();
}

Gecode_Sum::Gecode_Sum()
    : Gecode_Expression() {
    _offset = 0;
#ifdef _DEBUGWRAP
    std::cout << "creating sum D" << std::endl;
#endif
}

Gecode_Sum::~Gecode_Sum() {
#ifdef _DEBUGWRAP
    std::cout << "delete sum" << std::endl;
#endif
}

void Gecode_Sum::addVar(Gecode_Expression* v) {
#ifdef _DEBUGWRAP
    std::cout << "adding variable" << std::endl;
#endif
    _vars.add(v);
}

void Gecode_Sum::addWeight(const int w) {
#ifdef _DEBUGWRAP
    std::cout << "adding weight" << std::endl;
#endif
    _weights.add(w);
}

void Gecode_Sum::initialise(){
    int w;
    _lb = _offset;
    _ub = _offset;
    for(unsigned int i = 0; i < _vars.size(); ++i) {
        w = _weights.get_item(i);
        if( w > 0 ) {
            _lb += (w * _vars.get_item(i)->get_min());
            _ub += (w * _vars.get_item(i)->get_max());
        } else {
            _ub += (w * _vars.get_item(i)->get_min());
            _lb += (w * _vars.get_item(i)->get_max());
        }
    }
#ifdef _DEBUGWRAP
    std::cout << "init bounds of sum to " << _lb << ".." << _ub << std::endl;
#endif
}

Gecode_Expression* Gecode_Sum::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add sum constraint" << std::endl;
#endif
        _solver = solver;

        int i, n=_vars.size();
        Gecode::IntVarArgs x(_offset == 0 ? n : n + 1);
        Gecode::IntArgs w(_offset == 0 ? n : n + 1);
        bool weighted = false;

        for(i=0; i<n; ++i) {
            _vars.set_item(i, _vars.get_item(i)->add(_solver,false));
            if(_weights.size() > i) w[i] = _weights.get_item(i);
            else w[i] = 1;
            x[i] = _vars.get_item(i)->getGecodeVar();

            if(w[i] != 1) weighted = true;
        }

        if(n == 1){
            if(_offset == 0 && !weighted){
                // std::cout << "simplify single item sum 1" << std::endl;
                return _vars.get_item(0);
            } else if(_offset != 0 && !weighted){
                // std::cout << "simplify single item sum 2, offset:" << _offset << std::endl;
                // std::cout << "gecodevar:" << _vars.get_item(0)->getGecodeVar() << std::endl;
                _gcintrepr << expr(*(solver->gecodespace), _vars.get_item(0)->getGecodeVar() + _offset);
                return this;
            } else if(_offset == 0 && weighted){
                // std::cout << "simplify single item sum 3, w:" << w[0] << std::endl;
                _gcintrepr << expr(*(solver->gecodespace), _vars.get_item(0)->getGecodeVar() * w[0]);
                return this;
            } else if(_offset != 0 && weighted){
                // std::cout << "simplify single item sum 4, offset:" << _offset << " w:" << w[0] << std::endl;
                _gcintrepr << expr(*(solver->gecodespace), _vars.get_item(0)->getGecodeVar() * w[0] + _offset);
                return this;
            }
        } else if(n == 2){
            if(_offset == 0 && !weighted){
                // std::cout << "simplify single item sum 5" << std::endl;
                _gcintrepr << expr(*(solver->gecodespace), _vars.get_item(0)->getGecodeVar() + _vars.get_item(1)->getGecodeVar());
                return this;
            }
        }

        if(_offset!=0){  // Create a dummy constant variable to model the offset
            Gecode::IntVar dummy1(*(solver->gecodespace), 1, 1);
            x[n] = dummy1;
            w[n] = _offset;
            weighted = true;
        }

        Gecode_Expression *aux = new Gecode_Expression(_lb, _ub);
        aux = aux->add(solver, false);
        // std::cout << "lb:" << _lb << " _ub:" << _ub << std::endl;
        // std::cout << "_offset:" << _offset << std::endl;
        // std::cout << "w:" << w << std::endl;
        // std::cout << "x:" << x << std::endl;
        if(weighted){
            Gecode::linear(*(solver->gecodespace), w, x, Gecode::IRT_EQ, aux->getGecodeVar());
        } else {
            Gecode::linear(*(solver->gecodespace), x, Gecode::IRT_EQ, aux->getGecodeVar());
        }
        return aux;
    }

    return this;
}

int Gecode_Sum::get_value() const {
    int val=_offset, i, n=_vars.size(), w;
    for(i=0; i<n; ++i) {
        if(_weights.size() > i) w = _weights.get_item(i);
        else w = 1;
        val += w * _vars.get_item(i)->getGecodeVar().val();
    }
    return val;
}

/**
 * Binary constraints
 */

int operatorintmul(int a, int b){ return a * b;}

Gecode_mul::Gecode_mul(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(operatorintmul, var1, var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating mul constraint between two variables" << std::endl;
#endif
}

Gecode_mul::Gecode_mul(Gecode_Expression *var1, int constant)
    : Gecode_binop(operatorintmul, var1, constant) {
#ifdef _DEBUGWRAP
    std::cout << "creating mul constraint between variable and constant" << std::endl;
#endif
}

Gecode_mul::~Gecode_mul() {
#ifdef _DEBUGWRAP
    std::cout << "delete mul constraint" << std::endl;
#endif
}

Gecode_Expression* Gecode_mul::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add mul constraint" << std::endl;
#endif

        if(top_level) {
            std::cerr << "Error: mulitply expression not valid at the top-level, it should be used as a sub-expression." << std::endl;
            exit(1);
        }

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(_vars[1] != NULL) {
            _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() * _vars[1]->getGecodeVar());
        } else {
            _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() * _constant);
        }
    }
    return this;
}

int operatorintdiv(int a, int b){ return a * b;}

Gecode_div::Gecode_div(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(operatorintdiv, var1, var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating div constraint" << std::endl;
#endif
}

Gecode_div::Gecode_div(Gecode_Expression *var1, int constant)
    : Gecode_binop(operatorintdiv, var1, constant) {
#ifdef _DEBUGWRAP
    std::cout << "creating div constraint" << std::endl;
#endif
}


Gecode_div::~Gecode_div() {
#ifdef _DEBUGWRAP
    std::cout << "delete div" << std::endl;
#endif
}

Gecode_Expression* Gecode_div::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add div constraint" << std::endl;
#endif

        if(top_level) {
            std::cerr << "Error: division expression not valid at the top-level, it should be used as a sub-expression." << std::endl;
            exit(1);
        }

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(_vars[1] != NULL) {
            _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() / _vars[1]->getGecodeVar());
        } else {
            _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() / _constant);
        }
    }
    return this;
}

int operatorintmod(int a, int b){ return a * b;}

Gecode_mod::Gecode_mod(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(operatorintmod, var1, var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating mod predicate" << std::endl;
#endif
}

Gecode_mod::Gecode_mod(Gecode_Expression *var1, int constant)
    : Gecode_binop(operatorintmod, var1, constant) {
#ifdef _DEBUGWRAP
    std::cout << "creating mod predicate" << std::endl;
#endif
}

Gecode_mod::~Gecode_mod() {
#ifdef _DEBUGWRAP
    std::cout << "delete mod" << std::endl;
#endif
}

Gecode_Expression* Gecode_mod::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add mod predicate" << std::endl;
#endif

        if(top_level) {
            std::cerr << "Error: modulo expression not valid at the top-level, it should be used as a sub-expression." << std::endl;
            exit(1);
        }

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(_vars[1] != NULL) {
            _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() % _vars[1]->getGecodeVar());
        } else {
            _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() % _constant);
        }
    }
    return this;
}

Gecode_and::Gecode_and(Gecode_Expression *var1,
                       Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating and predicate" << std::endl;
#endif
}

Gecode_and::Gecode_and(Gecode_Expression *var1,
                       int constant)
    : Gecode_binop(var1,constant) {
#ifdef _DEBUGWRAP
    std::cout << "creating and predicate" << std::endl;
    std::cout << "You should consider a different model if you have a conjunction with a constant." << std::endl;
#endif
}

Gecode_and::~Gecode_and() {
#ifdef _DEBUGWRAP
    std::cout << "delete and" << std::endl;
#endif
}

Gecode_Expression* Gecode_and::add(GecodeSolver *solver,
                                   bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add and constraint" << std::endl;
#endif
        _solver = solver;

        if(top_level) {
            if(_vars[1] != NULL) {
                _vars[0] = _vars[0]->add(_solver, true);
                _vars[1] = _vars[1]->add(_solver, true);
            } else if(_constant) {
                _vars[0]->add(_solver, true);
            } else {
                // Fail, model unsatisfiable
                solver->gecodespace->fail();
            }
        } else {
            _vars[0] = _vars[0]->add(_solver, false);

            if(_vars[1] != NULL) {
                _vars[1] = _vars[1]->add(_solver, false);
                _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() * _vars[1]->getGecodeVar());
            
            } else {
                _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() * _constant);
            }
        }

    }
    return this;
}

Gecode_or::Gecode_or(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating or predicate" << std::endl;
#endif
}

Gecode_or::Gecode_or(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
#ifdef _DEBUGWRAP
    std::cout << "creating or predicate" << std::endl;
#endif
}

Gecode_or::~Gecode_or() {
#ifdef _DEBUGWRAP
    std::cout << "delete or" << std::endl;
#endif
}

Gecode_Expression* Gecode_or::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add or predicate" << std::endl;
#endif
        _solver = solver;

        if(top_level) {
            if(_vars[1] != NULL) {
                _vars[0] = _vars[0]->add(_solver, false);
                _vars[1] = _vars[1]->add(_solver, false);
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() + _vars[1]->getGecodeVar()) >= 1);

            } else if(_constant) { // OR trivalially satisfied
            } else {
                _vars[0]->add(_solver, true);
            }
        } else {
            _vars[0] = _vars[0]->add(_solver, false);

            if(_vars[1] != NULL) {
                _vars[1] = _vars[1]->add(_solver, false);
                _gcintrepr << expr(*(solver->gecodespace), _vars[0]->getGecodeVar() + _vars[1]->getGecodeVar() - _vars[0]->getGecodeVar() * _vars[1]->getGecodeVar());
            } else if(_constant){
                // trivially satisfied
                _gcintrepr << Gecode::IntVar(*(solver->gecodespace), 1, 1);
            } else {
                // satisfied if _vars[0] is
                _gcintrepr << _vars[0]->getGecodeVar();
            }
        }
    }
    return this;
}

Gecode_eq::Gecode_eq(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating equality" << std::endl;
#endif
}

Gecode_eq::Gecode_eq(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
#ifdef _DEBUGWRAP
    std::cout << "creating equality" << std::endl;
#endif
}

Gecode_eq::~Gecode_eq() {
#ifdef _DEBUGWRAP
    std::cout << "delete eq" << std::endl;
#endif
}

Gecode_Expression* Gecode_eq::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add equality constraint (" << _constant << ")" << std::endl;
#endif

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(top_level){
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() == _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() == _constant);
            }
        } else{  // Reified
            Gecode_Expression *r = new Gecode_IntVar();
            r = r->add(_solver, false);
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() == _vars[1]->getGecodeVar()) == r->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() == _constant) == r->getGecodeVar());
            }
            return r;
        }
    }

    return this;
}

/* Disequality operator */

Gecode_ne::Gecode_ne(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating notequal" << std::endl;
#endif
}

Gecode_ne::Gecode_ne(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
#ifdef _DEBUGWRAP
    std::cout << "creating notequal" << std::endl;
#endif
}

Gecode_ne::~Gecode_ne() {
#ifdef _DEBUGWRAP
    std::cout << "delete notequal" << std::endl;
#endif
}

Gecode_Expression* Gecode_ne::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add notequal constraint" << std::endl;
#endif

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(top_level){
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() != _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() != _constant);
            }
        } else{  // Reified
            Gecode_Expression *r = new Gecode_IntVar();
            r = r->add(_solver, false);
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() != _vars[1]->getGecodeVar()) == r->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() != _constant) == r->getGecodeVar());
            }
            return r;
        }
    }
    return this;
}

/* Disjunctive constraint operator */

Gecode_NoOverlap::Gecode_NoOverlap(Gecode_Expression *var1, Gecode_Expression *var2, int constant, int bonstant)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating nooverlap" << std::endl;
#endif
    _constant = constant;
    _bonstant = bonstant;
}

Gecode_NoOverlap::~Gecode_NoOverlap() {
#ifdef _DEBUGWRAP
    std::cout << "delete nooverlap" << std::endl;
#endif
}

Gecode_Expression* Gecode_NoOverlap::add(GecodeSolver *solver, bool top_level) {

    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add nooverlap constraint" << std::endl;
#endif

        if(!top_level){
            std::cerr << "NoOverlap sub-expression not supported by this solver." << std::endl;
            exit(1);
        }

        _solver = solver;

        _vars[0] = _vars[0]->add(_solver, false);
        _vars[1] = _vars[1]->add(_solver, false);

        Gecode::IntVarArgs scope(2);
        Gecode::IntArgs durations(2);
        scope[0] = _vars[0]->getGecodeVar();
        scope[1] = _vars[1]->getGecodeVar();
        durations[0] = _constant;
        durations[1] = _bonstant;

        unary(*(solver->gecodespace), scope, durations);
    }

    return this;
}


/* Leq operator */

Gecode_le::Gecode_le(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "Creating less than or equal" << std::endl;
#endif
}

Gecode_le::Gecode_le(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
#ifdef _DEBUGWRAP
    std::cout << "Creating less than or equal" << std::endl;
#endif
}


Gecode_le::~Gecode_le() {
#ifdef _DEBUGWRAP
    std::cout << "delete lessequal" << std::endl;
#endif
}

Gecode_Expression* Gecode_le::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add lessequal constraint" << std::endl;
#endif

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(top_level){
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() <= _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() <= _constant);
            }
        } else{  // Reified
            Gecode_Expression *r = new Gecode_IntVar();
            r = r->add(_solver, false);
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() <= _vars[1]->getGecodeVar()) == r->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() <= _constant) == r->getGecodeVar());
            }
            return r;
        }
    }
    return this;
}

/* Geq operator */

Gecode_ge::Gecode_ge(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "Creating a greater or equal constraint" << std::endl;
#endif
}

Gecode_ge::Gecode_ge(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
#ifdef _DEBUGWRAP
    std::cout << "Creating a greater or equal constraint" << std::endl;
#endif
}

Gecode_ge::~Gecode_ge() {
#ifdef _DEBUGWRAP
    std::cout << "delete greaterequal" << std::endl;
#endif
}

Gecode_Expression* Gecode_ge::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add greaterequal constraint" << std::endl;
#endif

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(top_level){
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() >= _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() >= _constant);
            }
        } else{  // Reified
            Gecode_Expression *r = new Gecode_IntVar();
            r = r->add(_solver, false);
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() >= _vars[1]->getGecodeVar()) == r->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() >= _constant) == r->getGecodeVar());
            }
            return r;
        }

    }

    return this;
}

/* Lt object */

Gecode_lt::Gecode_lt(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "Creating a less than constraint" << std::endl;
#endif
}

Gecode_lt::Gecode_lt(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
#ifdef _DEBUGWRAP
    std::cout << "Creating a less than constraint" << std::endl;
#endif
}

Gecode_lt::~Gecode_lt() {
#ifdef _DEBUGWRAP
    std::cout << "delete lessthan" << std::endl;
#endif
}

Gecode_Expression* Gecode_lt::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add lessthan constraint" << std::endl;
#endif

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(top_level){
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() < _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() < _constant);
            }
        } else{  // Reified
            Gecode_Expression *r = new Gecode_IntVar();
            r = r->add(_solver, false);
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() < _vars[1]->getGecodeVar()) == r->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() < _constant) == r->getGecodeVar());
            }
            return r;
        }

    }
    return this;
}

/* Gt object */

Gecode_gt::Gecode_gt(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "Creating a greater than constraint" << std::endl;
#endif
}

Gecode_gt::Gecode_gt(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
#ifdef _DEBUGWRAP
    std::cout << "Creating a greater than constraint" << std::endl;
#endif
}

Gecode_gt::~Gecode_gt() {
#ifdef _DEBUGWRAP
    std::cout << "delete greaterthan" << std::endl;
#endif
}

Gecode_Expression* Gecode_gt::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add greaterthan constraint" << std::endl;
#endif

        _solver = solver;
        _vars[0] = _vars[0]->add(_solver, false);
        if(_vars[1] != NULL) _vars[1] = _vars[1]->add(_solver, false);

        if(top_level){
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() > _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() > _constant);
            }
        } else{  // Reified
            Gecode_Expression *r = new Gecode_IntVar();
            r = r->add(_solver, false);
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() > _vars[1]->getGecodeVar()) == r->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), expr(*(solver->gecodespace), _vars[0]->getGecodeVar() > _constant) == r->getGecodeVar());
            }
            return r;
        }

    }
    return this;
}


/* Minimise object */

Gecode_Minimise::Gecode_Minimise(Gecode_Expression *var)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "Creating a minimise objective" << std::endl;
#endif
    _exp = var;
}

Gecode_Minimise::~Gecode_Minimise() {
#ifdef _DEBUGWRAP
    std::cout << "delete minimise" << std::endl;
#endif
}

Gecode_Expression* Gecode_Minimise::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add minimise objective" << std::endl;
#endif

        _solver = solver;
        _solver->isoptimisation = true;
        _exp = _exp->add(solver, false);
        solver->gecodespace->setGecodeCostVar(_exp->getGecodeVar());

        // std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        // exit(1);

        // //std::cout << 11 << std::endl;

        // // This will be the objective function
        // _exp = _exp->add(solver, false);

        // //std::cout << 22 << std::endl;

        // //std::cout << _solver->solver << std::endl;

        // if(top_level) {
        //     _solver->solver->minimize(_exp->_self);
        //     //_search_goal = new Mistral::Goal(Mistral::Goal::MINIMIZATION, _exp);
        // }

    }
    return this;
}

/* Maximise object */

Gecode_Maximise::Gecode_Maximise(Gecode_Expression *var)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "Creating a maximise objective" << std::endl;
#endif

    _exp = var;
}

Gecode_Maximise::~Gecode_Maximise() {
#ifdef _DEBUGWRAP
    std::cout << "delete maximise" << std::endl;
#endif
}

Gecode_Expression* Gecode_Maximise::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add maximise objective" << std::endl;
#endif
        _solver = solver;
        _solver->isoptimisation = true;
        _exp = _exp->add(solver, false);  
        solver->gecodespace->setGecodeCostVar(expr(*(solver->gecodespace), _exp->getGecodeVar() * -1));

        // FIXME should be using a minus view for this instead.
        // Gecode::Int::MinusView objvarview(_exp->getGecodeVar());
        // solver->gecodespace->setGecodeCostVar(objvarview);
    }
    return this;
}

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

GecodeSolver::GecodeSolver() {
#ifdef _DEBUGWRAP
    std::cout << "c (wrapper) creating solver" << std::endl;
#endif

    gecodespace = originalspace = new NJGecodeSpace();
    lastsolvestatus = UNKNOWN;
    isoptimisation = false;
    // solver = new Mistral::Solver();
    // solver->parameters.verbosity = 2;

    // //_search_goal = NULL;

    // _restart_policy_str = "geom";
    // _heuristic_randomization = 2;
    // _var_heuristic_str = "dom/wdeg";
    // _val_heuristic_str = "minval+guided";
}

GecodeSolver::~GecodeSolver() {
#ifdef _DEBUGWRAP
    std::cout << "c (wrapper) delete solver" << std::endl;
#endif
    delete gecodespace;
}

void GecodeSolver::add(Gecode_Expression* arg) {
#ifdef _DEBUGWRAP
    std::cout << "adding expression" << std::endl;
#endif

    arg->add(this, true);
}

void GecodeSolver::initialise(GecodeExpArray& arg) {
#ifdef _DEBUGWRAP
    std::cout << "Initialising solver with array of expressions" << std::endl;
#endif
    // FIXME add support for specifying decision variables.
}

void GecodeSolver::initialise() {
#ifdef _DEBUGWRAP
    std::cout << "Initialising solver with no expressions" << std::endl;
#endif
}

int GecodeSolver::num_vars() {
#ifdef _DEBUGWRAP
    std::cout << "return number of variables" << std::endl;
#endif
    return 0;
}

int GecodeSolver::solveAndRestart(const int policy,
                                  const unsigned int base,
                                  const double factor,
                                  const double decay,
                                  const int reinit) {
#ifdef _DEBUGWRAP
    std::cout << " SOLVE!!" << std::endl;
#endif

    return solve();
}

int GecodeSolver::solve() {

#ifdef _DEBUGWRAP
    std::cout << " SOLVE!! " << std::endl;
#endif

    gecodespace->closeNJ();
    // Gecode::branch(*gecodespace, gecodespace->getVar(0), Gecode::INT_VAL_MIN());
    // Gecode::branch(*gecodespace, gecodespace->getGecodeSearchVariables(), Gecode::INT_VAR_MIN_MIN(), Gecode::INT_VAL_SPLIT_MIN());
    std::cout << "first status:" << gecodespace->status() << " solved:" << Gecode::SS_SOLVED << " branch:" << Gecode::SS_BRANCH << std::endl;
    gecodespace->print();
    NJGecodeSpace *sol = NULL;
    if(isoptimisation){
        sol = Gecode::bab(gecodespace);
    } else {
        Gecode::DFS<NJGecodeSpace> gecodedfs(gecodespace);
        sol = gecodedfs.next();
    }

    std::cout << "gecodespace:" << gecodespace <<" sol:" << sol << std::endl;
    std::cout << "solve finished. gecodespace status:" << gecodespace->status() << std::endl;
    if(sol != NULL){
        lastsolvestatus = SAT;
        std::cout << "solve finished. sol status:" << sol->status() << std::endl;
        sol->print();
        gecodespace = sol;
    } else {
        lastsolvestatus = UNSAT;
    }

    return (is_sat());
}

int GecodeSolver::startNewSearch() {
#ifdef _DEBUGWRAP
    std::cout << "starting new search" << std::endl;
#endif
    return 0;
}

int GecodeSolver::getNextSolution() {
#ifdef _DEBUGWRAP
    std::cout << "getting next solution" << std::endl;
#endif
    return 0;
}

// int GecodeSolver::sacPreprocess(const int type) {
// #ifdef _DEBUGWRAP
//     std::cout << "running sac preprocessing" << std::endl;
// #endif
//     return 0;
// }

bool GecodeSolver::propagate() {
#ifdef _DEBUGWRAP
    std::cout << "propagating" << std::endl;
#endif
    return 0;
}

int GecodeSolver::get_level() {
#ifdef _DEBUGWRAP
    std::cout << "return current level of solver" << std::endl;
#endif
    return 0;
}

bool GecodeSolver::undo(const int nlevel) {
#ifdef _DEBUGWRAP
    std::cout << "Going back nlevel levels" << std::endl;
#endif
    return 0;
}

void GecodeSolver::save() {
#ifdef _DEBUGWRAP
    std::cout << "saving state" << std::endl;
#endif
}

int GecodeSolver::next(Gecode_Expression* x, int v) {
#ifdef _DEBUGWRAP
    std::cout << "get next" << std::endl;
#endif
    return 0;
}

void GecodeSolver::post(const char* op, Gecode_Expression* x, int v) {
#ifdef _DEBUGWRAP
    std::cout << "Posting a constraint" << std::endl;
#endif


}

void GecodeSolver::deduce(const char* op, Gecode_Expression* x, int v) {
#ifdef _DEBUGWRAP
    std::cout << "deducing" << std::endl;
#endif
}

void GecodeSolver::deduce() {
#ifdef _DEBUGWRAP
    std::cout << "deducing" << std::endl;
#endif
}

bool GecodeSolver::branch_right() {
#ifdef _DEBUGWRAP
    std::cout << "branching right" << std::endl;
#endif
    return 0;
}

void GecodeSolver::store_solution() {
#ifdef _DEBUGWRAP
    std::cout << "storing solution" << std::endl;
#endif
}

void GecodeSolver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand) {
#ifdef _DEBUGWRAP
    std::cout << "Setting heuristics" << std::endl;
#endif

    // _heuristic_randomization = rand;
    // _var_heuristic_str = var_heuristic;
    // _val_heuristic_str = val_heuristic;

}

void GecodeSolver::addNogood(GecodeExpArray& vars,
                             GecodeIntArray& vals) {
#ifdef _DEBUGWRAP
    std::cout << "Adding nogood constraint" <<std::endl;
#endif
}

void GecodeSolver::guide(GecodeExpArray& vars,
                         GecodeIntArray& vals,
                         GecodeDoubleArray& probs) {
#ifdef _DEBUGWRAP
    std::cout << "Adding nogood constraint" <<std::endl;
#endif
}

void GecodeSolver::forceFiniteDomain(GecodeExpArray& vars) {
#ifdef _DEBUGWRAP
    std::cout << "Forcing finite domain" <<std::endl;
#endif
}

void GecodeSolver::backtrackTo(const int level) {
#ifdef _DEBUGWRAP
    std::cout << "backtracking to level" <<std::endl;
#endif
}

void GecodeSolver::presolve() {
#ifdef _DEBUGWRAP
    std::cout << "presolving" <<std::endl;
#endif
}

void GecodeSolver::increase_init_level(const int i) {
#ifdef _DEBUGWRAP
    std::cout << "increasing initial level" <<std::endl;
#endif
}

void GecodeSolver::decrease_init_level(const int i) {
#ifdef _DEBUGWRAP
    std::cout << "decreasing initial level" <<std::endl;
#endif
}

void GecodeSolver::assign(Gecode_Expression *X, const int v) {
#ifdef _DEBUGWRAP
    std::cout << "Setting domain of expression X to v" <<std::endl;
#endif
}

void GecodeSolver::upOneLevel() {
#ifdef _DEBUGWRAP
    std::cout << "stepping up one level" <<std::endl;
#endif
}

void GecodeSolver::reset(bool full) {
#ifdef _DEBUGWRAP
    std::cout << "resetting solver" <<std::endl;
#endif

    // solver->restore(solver->search_root);
}

void GecodeSolver::setLowerBounds(GecodeExpArray& vars,
                                  GecodeIntArray& vals) {
#ifdef _DEBUGWRAP
    std::cout << "stepping up one level" <<std::endl;
#endif
}

void GecodeSolver::setUpperBounds(GecodeExpArray& vars,
                                  GecodeIntArray& vals) {
#ifdef _DEBUGWRAP
    std::cout << "stetting upper bounds" <<std::endl;
#endif
}

void GecodeSolver::setRestartNogood() {
#ifdef _DEBUGWRAP
    std::cout << "setting restart no good" <<std::endl;
#endif
}

void GecodeSolver::setFailureLimit(const int cutoff) {
#ifdef _DEBUGWRAP
    std::cout << "setting failure limit" <<std::endl;
#endif

    // solver->parameters.fail_limit = cutoff;

}

void GecodeSolver::setNodeLimit(const int cutoff) {
#ifdef _DEBUGWRAP
    std::cout << "setting node limit" <<std::endl;
#endif

    // solver->parameters.node_limit = cutoff;

}

void GecodeSolver::setTimeLimit(const int cutoff) {
#ifdef _DEBUGWRAP
    std::cout << "setting time limit" <<std::endl;
#endif

    // solver->parameters.time_limit = cutoff;

}

void GecodeSolver::setVerbosity(const int degree) {
#ifdef _DEBUGWRAP
    std::cout << "setting verbosity to " << degree <<std::endl;
#endif

    // solver->parameters.verbosity = degree;

}

void GecodeSolver::setRandomized(const int degree) {
#ifdef _DEBUGWRAP
    std::cout << "setting randomised" <<std::endl;
#endif

    // solver->parameters.randomization = degree;

}

void GecodeSolver::setRandomSeed(const int seed) {
#ifdef _DEBUGWRAP
    std::cout << "setting random seed" <<std::endl;
#endif

    // solver->parameters.seed = seed;

}

bool GecodeSolver::is_opt() {
// #ifdef _DEBUGWRAP
//     std::cout << "returning is_opt()" <<std::endl;
// #endif
    return false;  // FIXME
}

bool GecodeSolver::is_sat() {
// #ifdef _DEBUGWRAP
//     std::cout << "returning is satisfied?" <<std::endl;
// #endif
    return gecodespace->status() == Gecode::SS_SOLVED;
}

bool GecodeSolver::is_unsat() {
// #ifdef _DEBUGWRAP
//     std::cout << "returning is NOT satisfied?" <<std::endl;
// #endif
    return lastsolvestatus == UNSAT;
}

void GecodeSolver::printStatistics() {
#ifdef _DEBUGWRAP
    std::cout << "printing statistics" <<std::endl;
#endif
    // FIXME
}

int GecodeSolver::getBacktracks() {
#ifdef _DEBUGWRAP
    std::cout << "return number of backtracks" <<std::endl;
#endif
    // return solver->statistics.num_backtracks;
    return 0;
}

int GecodeSolver::getNodes() {
#ifdef _DEBUGWRAP
    std::cout << "return number of nodes" <<std::endl;
#endif
    // return solver->statistics.num_nodes;
    return 0;
}

int GecodeSolver::getFailures() {
#ifdef _DEBUGWRAP
    std::cout << "return number of failures" <<std::endl;
#endif
    // return solver->statistics.num_failures;
    return 0;
}

int GecodeSolver::getChecks() {
#ifdef _DEBUGWRAP
    std::cout << "return number of checks" <<std::endl;
#endif
    // return solver->statistics.num_propagations;
    return 0;
}

int GecodeSolver::getPropags() {
#ifdef _DEBUGWRAP
    std::cout << "return amount of propagation" <<std::endl;
#endif
    // return solver->statistics.num_propagations;
    return 0;
}

double GecodeSolver::getTime() {
#ifdef _DEBUGWRAP
    std::cout << "return duration" <<std::endl;
#endif
    // return solver->statistics.end_time - solver->statistics.start_time;
    return 0;
}

int GecodeSolver::getNumVariables() {
#ifdef _DEBUGWRAP
    std::cout << "return number of variables" <<std::endl;
#endif
    return gecodespace->getNumVariables();
}

int GecodeSolver::getNumConstraints() {
#ifdef _DEBUGWRAP
    std::cout << "return number of constraints" <<std::endl;
#endif
    // return solver->posted_constraints.size;
    return 0;
}

int GecodeSolver::getRandomNumber() {
    // why is this here? deprecate
    return 42;
}
