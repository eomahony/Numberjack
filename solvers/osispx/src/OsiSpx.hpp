#include "Osi.hpp"
#include <OsiSpxSolverInterface.hpp>

class OsiSpxSolver: public OsiSolver {
private:
	OsiSpxSolverInterface* spxS;
public:
	OsiSpxSolver();
	virtual ~OsiSpxSolver();
	void setTimeLimit(const int cutoff);
};
