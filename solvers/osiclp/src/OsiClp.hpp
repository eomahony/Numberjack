#include "Osi.hpp"
#include <OsiClpSolverInterface.hpp>

class OsiClpSolver: public OsiSolver {
private:
	OsiClpSolverInterface* clpSi;
public:
	OsiClpSolver();
	virtual ~OsiClpSolver();
};
