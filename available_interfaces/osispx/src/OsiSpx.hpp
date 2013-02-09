#include "Osi.hpp"
#include <OsiSpxSolverInterface.hpp>

class OsiSpxSolver: public OsiSolver {
private:
	OsiSpxSolverInterface* spxS;
public:
	OsiSpxSolver();
	virtual ~OsiSpxSolver();
	int solve();
	void setTimeLimit(const int cutoff);
};
