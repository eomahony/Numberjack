#include "Gecode.hpp"

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
    // Mistral::Variable x(0,1);
    // _self = x;
}

Gecode_Expression::Gecode_Expression(const int nval) {

#ifdef _DEBUGWRAP
    std::cout << "creating a variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif

    _solver = NULL;
    _gcvarid = -1;
    _vartype = INTVARNVAL;
    _nval = nval;
    // Mistral::Variable x(0,nval-1);
    // _self = x;
}

Gecode_Expression::Gecode_Expression(const int lb, const int ub) {

#ifdef _DEBUGWRAP
    std::cout << "creating a variable [" << lb << ".." << ub << "]" << std::endl;
#endif

    _solver = NULL;
    _gcvarid = -1;
    _vartype = INTVARBOUNDS;
    _lb = lb; _ub = ub;
    // Mistral::Variable x(lb,ub);
    // _self = x;
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

    // return _self.id() ;
    return 0;
}

int Gecode_Expression::get_value() const {

#ifdef _DEBUGWRAP
    std::cout << "return value of expression" << std::endl;
#endif

    if(_solver != NULL && _gcvarid >= 0){
        std::cout << "get value for variable id " << _gcvarid << std::endl;
        try {
            std::cout << "var:" << _solver->gecodespace->getVar(_gcvarid) << std::endl;
            return _solver->gecodespace->getVar(_gcvarid).val();
        } catch(Gecode::Int::ValOfUnassignedVar e){
            std::cout << "ignore unassigned var" << std::endl;
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
    return _solver->gecodespace->getVar(_gcvarid);
}

Gecode_Expression* Gecode_Expression::add(GecodeSolver *solver, bool top_level) {

#ifdef _DEBUGWRAP
    std::cout << "add expression to model" << std::endl;
#endif

    if(top_level) {
#ifdef _DEBUGWRAP
        std::cout << "\tAdding at top level" << std::endl;
#endif
    } else {
#ifdef _DEBUGWRAP
        std::cout << "\tAdding within tree" <<std::endl;
#endif
    }

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

// Gecode_binop::Gecode_binop(int constant, Gecode_Expression *var1)
//   : Gecode_Expression()
// {

// #ifdef _DEBUGWRAP
//   std::cout << "creating a binary (constant) operator" << std::endl;
// #endif
//   _vars[0] = NULL;
//   _vars[1] = var1;
//   _constant = constant;
// }


Gecode_binop::~Gecode_binop() {

#ifdef _DEBUGWRAP
    std::cout << "delete binary operator" << std::endl;
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

}

Gecode_Min::Gecode_Min( Gecode_Expression *var1, Gecode_Expression *var2 )
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a binary Min constraint" << std::endl;
#endif

    _vars.add(var1);
    _vars.add(var2);
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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // int i, n=_vars.size();
        // for(i=0; i<n; ++i)
        //     _vars.get_item(i)->add(_solver,false);

        // if(n == 2) {
        //     _self = Min(_vars.get_item(0)->_self,
        //                 _vars.get_item(1)->_self);
        // } else if(n > 2) {
        //     Mistral::VarArray scope(n);
        //     for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
        //     _self = Min(scope);
        // }

        // if( top_level )
        //     _solver->solver->add( _self );

    }
    return this;
}


Gecode_Max::Gecode_Max( GecodeExpArray& vars )
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a Max constraint" << std::endl;
#endif

    _vars = vars;
}

Gecode_Max::Gecode_Max( Gecode_Expression *var1, Gecode_Expression *var2 )
    : Gecode_Expression() {

#ifdef _DEBUGWRAP
    std::cout << "creating a binary Max constraint" << std::endl;
#endif

    _vars.add(var1);
    _vars.add(var2);
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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // int i, n=_vars.size();
        // for(i=0; i<n; ++i)
        //     _vars.get_item(i)->add(_solver,false);

        // if(n == 2) {
        //     _self = Max(_vars.get_item(0)->_self,
        //                 _vars.get_item(1)->_self);
        // } else if(n > 2) {
        //     Mistral::VarArray scope(n);
        //     for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
        //     _self = Max(scope);
        // }

        // if( top_level )
        //     _solver->solver->add( _self );

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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // int i, n=_vars.size(), m=_vals.size();
        // int min_val=_vals.get_item(0);
        // int max_val=_vals.get_item(m-1);
        // int M = (max_val - min_val + 1);
        // Mistral::VarArray scope(n);

        // int *tmp_lb = new int[M];
        // int *tmp_ub = new int[M];

        // for(i=0; i<n; ++i) {
        //     _vars.get_item(i)->add(_solver,false);
        //     scope[i] = _vars.get_item(i)->_self;
        // }

        // //std::cout << min_val << " to " << max_val << std::endl;

        // for(i=0; i<M; ++i) {
        //     tmp_lb[i] = _lb_card.get_item(i);
        //     tmp_ub[i] = _ub_card.get_item(i);

        //     //std::cout << " " << i+min_val << ": in [" << tmp_lb[i] << ".." << tmp_ub[i] << "]\n";
        // }

        // _self = Occurrences(scope, min_val, max_val, tmp_lb, tmp_ub);

        // if( top_level )
        //     _solver->solver->add( _self );

        // delete [] tmp_lb;
        // delete [] tmp_ub;

    }

    return this;
}

Gecode_Element::Gecode_Element( GecodeExpArray& vars )
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating element" << std::endl;
#endif

    _vars = vars;

}

Gecode_Element::~Gecode_Element() {
#ifdef _DEBUGWRAP
    std::cout << "delete element" << std::endl;
#endif
}

Gecode_Expression* Gecode_Element::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add element constraint" << std::endl;
#endif
        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // int i, n=_vars.size();
        // Mistral::VarArray scope(n-1);



        // for(i=0; i<n-1; ++i) {
        //     _vars.get_item(i)->add(_solver,false);
        //     scope[i] = _vars.get_item(i)->_self;
        // }

        // _vars.get_item(n-1)->add(_solver,false);
        // Mistral::Variable index = _vars.get_item(n-1)->_self;
        // _self = scope[index];

    }

    return this;
}

Gecode_LeqLex::Gecode_LeqLex( GecodeExpArray& vars )
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating lexleq" << std::endl;
#endif
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


        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        for(int i = 0; i < _vars.size(); ++i)
            _vars.set_item(i, _vars.get_item(i)->add(solver, false));

        if(top_level) {
#ifdef _DEBUGWRAP
            std::cout << "\tAdding at top level" << std::endl;
#endif
        } else {
#ifdef _DEBUGWRAP
            std::cout << "\tAdding within tree" <<std::endl;
#endif
        }

    }
    return this;
}


Gecode_LessLex::Gecode_LessLex( GecodeExpArray& vars )
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating lexless" << std::endl;
#endif
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

        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        for(int i = 0; i < _vars.size(); ++i)
            _vars.set_item(i, _vars.get_item(i)->add(solver, false));

        if(top_level) {
#ifdef _DEBUGWRAP
            std::cout << "\tAdding at top level" << std::endl;
#endif
        } else {
#ifdef _DEBUGWRAP
            std::cout << "\tAdding within tree" <<std::endl;
#endif
        }

    }
    return this;
}


Gecode_Sum::Gecode_Sum(GecodeExpArray& vars,
                       GecodeIntArray& weights,
                       const int offset)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating sum" << std::endl;
#endif
    _vars = vars;
    _weights = weights;
    _offset = offset;
}

Gecode_Sum::Gecode_Sum(Gecode_Expression *arg1,
                       Gecode_Expression *arg2,
                       GecodeIntArray& w,
                       const int offset)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating sum" << std::endl;
#endif
    _vars.add(arg1);
    _vars.add(arg2);
    _weights = w;
    _offset = offset;
}

Gecode_Sum::Gecode_Sum(Gecode_Expression *arg,
                       GecodeIntArray& w,
                       const int offset)
    : Gecode_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating sum" << std::endl;
#endif
    _vars.add(arg);
    _weights = w;
    _offset = offset;
}

Gecode_Sum::Gecode_Sum()
    : Gecode_Expression() {
    _offset = 0;
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
}

void Gecode_Sum::addWeight(const int w) {
#ifdef _DEBUGWRAP
    std::cout << "adding weight" << std::endl;
#endif
}

Gecode_Expression* Gecode_Sum::add(GecodeSolver *solver, bool top_level) {
    if(!has_been_added()) {
#ifdef _DEBUGWRAP
        std::cout << "add sum constraint" << std::endl;
#endif
        _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // _v * _w + o = _x
        // _v * _w - _x = -o


        // int i, n=_vars.size();
        // Mistral::VarArray scope(n);
        // Mistral::Vector<int> w(n);
        // bool weighted = false;
        // bool boolean = true;


        // for(i=0; i<n; ++i) {

        //     _vars.get_item(i)->add(_solver,false);
        //     if(_weights.size() > i)
        //         w[i] = _weights.get_item(i);
        //     else
        //         w[i] = 1;
        //     scope[i] = _vars.get_item(i)->_self;

        //     if(w[i] != 1) weighted = true;
        //     if(!scope[i].is_boolean()) boolean = false;
        // }
        // //w[n] = _weights.get_item(n);


        // if(boolean) {
        //     if(weighted)
        //         _self = BoolSum(scope, w);
        //     else
        //         _self = BoolSum(scope, w);
        // } else {
        //     if(weighted)
        //         _self = Sum(scope, w, -Mistral::INFTY, Mistral::INFTY, _offset);
        //     else
        //         _self = Sum(scope, w, -Mistral::INFTY, Mistral::INFTY, _offset);
        // }

        // if( top_level ) {
        //     _solver->solver->add( _self );
        // }
    }

    return this;
}

/**
 * Binary constraints
 */
Gecode_mul::Gecode_mul(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating mul constraint between two variables" << std::endl;
#endif
}

Gecode_mul::Gecode_mul(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // _vars[0]->add(_solver,false);
        // if(_vars[1]) {
        //     _vars[1]->add(_solver,false);
        //     _self = (_vars[0]->_self *
        //              _vars[1]->_self);
        // } else {
        //     _self = (_vars[0]->_self *
        //              _constant);
        // }
        // if( top_level )
        //     _solver->solver->add( _self );

    }
    return this;
}

Gecode_div::Gecode_div(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating div constraint" << std::endl;
#endif
}

Gecode_div::Gecode_div(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // if(_vars[0]) {
        //     _vars[0]->add(_solver,false);

        //     if(_vars[1]) {
        //         _vars[1]->add(_solver,false);
        //         _self = (_vars[0]->_self /
        //                  _vars[1]->_self);
        //     } else {
        //         _self = (_vars[0]->_self /
        //                  _constant);
        //     }
        // } else {

        //     std::cout << "Cannot divide constants" << std::endl;
        //     exit(1);

        // }
        // if( top_level )
        //     _solver->solver->add( _self );

    }
    return this;
}

Gecode_mod::Gecode_mod(Gecode_Expression *var1, Gecode_Expression *var2)
    : Gecode_binop(var1,var2) {
#ifdef _DEBUGWRAP
    std::cout << "creating mod predicate" << std::endl;
#endif
}

Gecode_mod::Gecode_mod(Gecode_Expression *var1, int constant)
    : Gecode_binop(var1,constant) {
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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // _vars[0]->add(_solver,false);
        // if(_vars[1]) {
        //     _vars[1]->add(_solver,false);
        //     _self = (_vars[0]->_self % _vars[1]->_self);
        // } else {
        //     _self = (_vars[0]->_self % _constant);
        // }
        // if( top_level )
        //     _solver->solver->add( _self );

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
    std::cout << "I DON'T THINK I SHOULD BE HERE" << std::endl;
#endif

    /**
     * Should never be in this constructor???
     */

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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // _vars[0]->add(_solver,false);
        // if(_vars[1]) {
        //     _vars[1]->add(_solver,false);
        //     _self = (_vars[0]->_self && _vars[1]->_self);
        // } else if(_constant) {
        //     _self = (_vars[0]->_self == 1);
        // } else {
        //     _solver->solver->fail();
        // }
        // if( top_level )
        //     _solver->solver->add( _self );

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
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // bool used = false;
        // _vars[0]->add(_solver,false);
        // if(_vars[1]) {
        //     _vars[1]->add(_solver,false);
        //     _self = (_vars[0]->_self || _vars[1]->_self);
        //     used = true;
        // } else if(_constant == 0) {
        //     _self = (_vars[0]->_self == 1);
        //     used = true;
        // }
        // if( top_level && used )
        //     _solver->solver->add( _self );

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
            std::cout << "creating rel eq" << std::endl;
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() == _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() == _constant);
            }
        } else{  // Reified
            std::cout << "creating reified rel eq" << std::endl;
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
            std::cout << "creating rel ne" << std::endl;
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() != _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() != _constant);
            }
        } else{  // Reified
            std::cout << "creating reified rel ne" << std::endl;
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
        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // _vars[0]->add(_solver,false);
        // _vars[1]->add(_solver,false);

        // _self = Disjunctive(_vars[0]->_self, _vars[1]->_self, _constant, _bonstant);

        // if( top_level )
        //     _solver->solver->add( _self );
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
            std::cout << "creating rel le" << std::endl;
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() <= _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() <= _constant);
            }
        } else{  // Reified
            std::cout << "creating reified rel le" << std::endl;
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
            std::cout << "creating rel ge" << std::endl;
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() >= _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() >= _constant);
            }
        } else{  // Reified
            std::cout << "creating reified rel ge" << std::endl;
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
            std::cout << "creating rel lt" << std::endl;
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() < _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() < _constant);
            }
        } else{  // Reified
            std::cout << "creating reified rel lt" << std::endl;
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
            std::cout << "creating rel gt" << std::endl;
            if(_vars[1] != NULL) {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() > _vars[1]->getGecodeVar());
            } else {
                Gecode::rel(*(solver->gecodespace), _vars[0]->getGecodeVar() > _constant);
            }
        } else{  // Reified
            std::cout << "creating reified rel gt" << std::endl;
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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

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

        // _solver = solver;
        std::cerr << "Error constraint not supported with this solver, yet." << std::endl;
        exit(1);

        // // This will be the objective function
        // _exp = _exp->add(solver, false);

        // if(top_level) {
        //     _solver->solver->maximize(_exp->_self);
        //     //_search_goal = new Mistral::Goal(Mistral::Goal::MAXIMIZATION, _exp);
        // }

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
    // delete solver;
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
    Gecode::DFS<NJGecodeSpace> gecodedfs(gecodespace);
    NJGecodeSpace *sol = gecodedfs.next();
    // gecodespace->choice();
    // NJGecodeSpace *sol = gecodespace;
    std::cout << "gecodespace:" << gecodespace <<" sol:" << sol << std::endl;
    std::cout << "solve finished. status:" << sol->status() << std::endl;
    sol->print();

    gecodespace = sol;

//     solver->consolidate();


// #ifdef _DEBUGWRAP
//     std::cout << solver << std::endl;
//     solver->parameters.verbosity = 2;
// #endif

//     _branching_heuristic = solver->heuristic_factory(_var_heuristic_str, _val_heuristic_str, _heuristic_randomization);

//     //if(!_search_goal)
//     //_search_goal = new Mistral::Goal(Mistral::Goal::SATISFACTION);
//     _restart_policy = solver->restart_factory(_restart_policy_str);

//     //Mistral::Outcome result =
//     solver->depth_first_search(solver->variables, _branching_heuristic, _restart_policy, NULL, false); //, _search_goal);

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
#ifdef _DEBUGWRAP
    std::cout << "returning is_opt()" <<std::endl;
#endif
    // return solver->statistics.outcome == 3;
    return false;
}

bool GecodeSolver::is_sat() {
#ifdef _DEBUGWRAP
    std::cout << "returning is satisfied?" <<std::endl;
#endif
    // return solver->statistics.num_solutions > 0; //solver->statistics.outcome == Mistral::SAT || solver->statistics.outcome == Mistral::OPT;
    // return false;
    return gecodespace->status() == Gecode::SS_SOLVED;
}

bool GecodeSolver::is_unsat() {
#ifdef _DEBUGWRAP
    std::cout << "returning is NOT satisfied?" <<std::endl;
#endif
    // return solver->statistics.outcome == 0;
    return false;
}

void GecodeSolver::printStatistics() {
#ifdef _DEBUGWRAP
    std::cout << "printing statistics" <<std::endl;
#endif

    // solver->statistics.print_full(std::cout);

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
    // return solver->variables.size;
    return 0;
}

int GecodeSolver::getNumConstraints() {
#ifdef _DEBUGWRAP
    std::cout << "return number of constraints" <<std::endl;
#endif
    // return solver->posted_constraints.size;
    return 0;
}

int GecodeSolver::getRandomNumber() {
    // return Mistral::randint(0xffffffff);
    return 42;
}
