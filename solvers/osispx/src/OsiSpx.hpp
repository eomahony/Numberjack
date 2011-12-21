#include "Osi.hpp"
#include <coin/OsiSpxSolverInterface.hpp>

class OsiSpxSolver: public OsiSolver {
private:
public:
	OsiSpxSolver();
	virtual ~OsiSpxSolver();
};
