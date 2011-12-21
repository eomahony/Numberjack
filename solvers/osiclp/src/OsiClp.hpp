#include "Osi.hpp"
#include <OsiClpSolverInterface.hpp>

class OsiClpSolver: public OsiSolver {
private:
public:
	OsiClpSolver();
	virtual ~OsiClpSolver();
};
