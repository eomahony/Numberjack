#include "Osi.hpp"
#include <OsiDylpSolverInterface.hpp>

class OsiDylpSolver: public OsiSolver {
private:
public:
	int solve();
	OsiDylpSolver();
	virtual ~OsiDylpSolver();
};
