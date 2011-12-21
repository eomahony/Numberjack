#include "Osi.hpp"
#include <OsiCbcSolverInterface.hpp>

class OsiCbcSolver: public OsiSolver {
private:
public:
	int getNodes();
	OsiCbcSolver();
	virtual ~OsiCbcSolver();
};
