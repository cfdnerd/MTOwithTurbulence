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
#include <diff.c>

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
    while (simple.loop(runTime))
    {
        #include "update.H"
        turbulence->correct();
        #include "Primal_U.H"
        #include "Primal_T.H"
        #include "AdjointHeat_Tb.H"
        #include "Adjoint_kOmegaSST.H"
        #include "AdjointHeat_Ub.H"
        #include "AdjointFlow_Ua.H"
        #include "costfunction.H"              
        #include "sensitivity.H"
    }
    #include "finalize.H"
    return 0;
}
