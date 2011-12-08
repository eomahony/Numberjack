#include "Osi.hpp"
#include <coin/OsiGlpkSolverInterface.hpp>

class OsiGlpkSolver: public OsiSolver {
private:
public:
	OsiGlpkSolver();
	virtual ~OsiGlpkSolver();
};
