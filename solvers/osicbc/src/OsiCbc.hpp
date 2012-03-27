#include "Osi.hpp"
#include <OsiCbcSolverInterface.hpp>

class OsiCbcSolver: public OsiSolver {
private:
public:
	int getNodes();
    void setLPRelaxationSolver(OsiSolver njSolver); 
	OsiCbcSolver();
	virtual ~OsiCbcSolver();
};
