#include "Osi.hpp"
#include <coin/OsiCbcSolverInterface.hpp>

class OsiCbcSolver: public OsiSolver {
private:
public:
	OsiCbcSolver();
	virtual ~OsiCbcSolver();
};
