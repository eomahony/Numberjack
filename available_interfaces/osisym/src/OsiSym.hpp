#include "Osi.hpp"
#include <OsiSymSolverInterface.hpp>

class OsiSymSolver: public OsiSolver {
private:
public:
	OsiSymSolver();
	virtual ~OsiSymSolver();
};
