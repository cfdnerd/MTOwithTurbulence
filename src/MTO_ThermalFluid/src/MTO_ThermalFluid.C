//Author: Yu Minghao    Updated: May 2020 

static char help[] = "topology optimization of fluid problem\n";
#define MPI_NO_CPPBIND 1
#include "mpi.h"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"//
#include "MMA/MMA.h"
#include <cmath>
#include <diff.c>
#include "solver_status.H"
#include "solver_diagnostics.H"

int main(int argc, char *argv[])
{
    // Do not call MPI_Init here - OpenFOAM initializes MPI during setRootCase/startup.
    // Calling MPI_Init explicitly causes "MPI_Init called twice" when run with mpirun.

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFvOptions.H"//
    #include "createFields.H"
    #include "readTransportProperties.H" 
    #include "initContinuityErrs.H"
    #include "readThermalProperties.H" 
    #include "opt_initialization.H"
    MTOdiag::setConfig(diagLevel, diagEvery, diagRank, diagOptStart);
    while (simple.loop(runTime))
    {
        #include "update.H"
        turbulence->correct();
        #include "Primal_U.H"
        #include "Primal_T.H"
        #include "AdjointHeat_Tb.H"
        const bool doTurbAdjSolveThisOpt =
        (
            adjTurbLagEvery <= 1
         || opt <= adjTurbRampEnd
         || ((opt % adjTurbLagEvery) == 0)
        );
        MTOdiag::logMetric(runTime, opt, "AdjointLoop", "doTurbAdjSolve", doTurbAdjSolveThisOpt ? 1.0 : 0.0);
        for (label pass = 0; pass < nAdjTurbPasses; pass++)
        {
            if (doTurbAdjSolveThisOpt || pass == 0)
            {
                #include "Adjoint_kOmegaSST.H"
            }
            #include "AdjointHeat_Ub.H"
            #include "AdjointFlow_Ua.H"
        }
        #include "costfunction.H"              
        #include "sensitivity.H"
        #include "writeSolverStatus.H"
    }
    #include "finalize.H"
    return 0;
}
