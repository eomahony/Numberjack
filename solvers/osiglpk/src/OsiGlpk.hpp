#include "Osi.hpp"
#include <OsiGlpkSolverInterface.hpp>

class OsiGlpkSolver: public OsiSolver {
private:
public:
	OsiGlpkSolver();
	virtual ~OsiGlpkSolver();
};
