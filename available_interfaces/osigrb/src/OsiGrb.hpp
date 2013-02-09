#include "Osi.hpp"
#include <OsiGrbSolverInterface.hpp>

class OsiGrbSolver: public OsiSolver {
private:
	OsiGrbSolverInterface* clpSi;
public:
	OsiGrbSolver();
	virtual ~OsiGrbSolver();
};
