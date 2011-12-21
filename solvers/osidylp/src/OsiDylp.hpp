#include "Osi.hpp"
#include <OsiDylpSolverInterface.hpp>

class OsiDylpSolver: public OsiSolver {
private:
public:
	OsiDylpSolver();
	virtual ~OsiDylpSolver();
};
