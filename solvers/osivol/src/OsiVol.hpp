#include "Osi.hpp"
#include <OsiVolSolverInterface.hpp>

class OsiVolSolver: public OsiSolver {
private:
public:
	int solve();
	OsiVolSolver();
	virtual ~OsiVolSolver();
};
