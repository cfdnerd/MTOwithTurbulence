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
    #include "setRootCase.H"
    // Ensure MPI is initialized before opt_initialization (MMA constructor calls MPI_Allreduce).
    // Some run configurations (job schedulers, non-mpirun) may not init MPI; guard against MPI_Allreduce-before-MPI_Init.
    {
        int mpiAlreadyInit = 0;
        MPI_Initialized(&mpiAlreadyInit);
        if (!mpiAlreadyInit)
        {
            MPI_Init(&argc, &argv);
        }
    }
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
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "loop_enter", 1.0);
        #include "update.H"
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "after_update", 1.0);
        turbulence->correct();
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "after_turbulence_correct", 1.0);
        #include "Primal_U.H"
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "after_primal_U", 1.0);
        #include "Primal_T.H"
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "after_primal_T", 1.0);
        #include "AdjointHeat_Tb.H"
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "after_adjoint_Tb", 1.0);
        const bool doTurbAdjSolveThisOpt =
        (
            adjTurbLagEvery <= 1
         || opt <= adjTurbRampEnd
         || ((opt % adjTurbLagEvery) == 0)
        );
        MTOdiag::logMetric(runTime, opt, "AdjointLoop", "doTurbAdjSolve", doTurbAdjSolveThisOpt ? 1.0 : 0.0);
        for (label pass = 0; pass < nAdjTurbPasses; pass++)
        {
            MTOdiag::logMetric(runTime, opt, "Heartbeat", word("adj_pass_") + Foam::name(pass) + "_enter", 1.0);
            if (doTurbAdjSolveThisOpt || pass == 0)
            {
                #include "Adjoint_kOmegaSST.H"
                MTOdiag::logMetric(runTime, opt, "Heartbeat", word("adj_pass_") + Foam::name(pass) + "_after_kOmega", 1.0);
            }
            #include "AdjointHeat_Ub.H"
            MTOdiag::logMetric(runTime, opt, "Heartbeat", word("adj_pass_") + Foam::name(pass) + "_after_Ub", 1.0);
            #include "AdjointFlow_Ua.H"
            MTOdiag::logMetric(runTime, opt, "Heartbeat", word("adj_pass_") + Foam::name(pass) + "_after_Ua", 1.0);
        }
        #include "costfunction.H"              
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "after_costfunction", 1.0);
        #include "sensitivity.H"
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "after_sensitivity", 1.0);
        #include "writeSolverStatus.H"
        MTOdiag::logMetric(runTime, opt, "Heartbeat", "after_writeSolverStatus", 1.0);
    }
    #include "finalize.H"
    return 0;
}
