/*---------------------------------------------------------------------------*\
  solver_status.H - Solver performance logging for MTO_ThermalFluid
  Logs residuals, iteration counts, and convergence status for debugging
\*---------------------------------------------------------------------------*/

#ifndef solver_status_H
#define solver_status_H

#include "OFstream.H"

//- Log solver performance to solver_status.log (master rank only)
//  Format: optIter,timeStep,block,outerIter,equation,initRes,finalRes,nIter,converged,singular
template<class Type>
inline void logSolverStats
(
    const word& blockName,
    const label outerIter,
    const word& eqnName,
    const SolverPerformance<Type>& sp,
    const label optIter,
    const word& timeName
)
{
    if (!Pstream::master()) return;

    const fileName logPath(runTime.path()/"solver_status.log");
    OFstream logFile(logPath, IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED, true);
    if (!logFile.good()) return;

    const scalar initRes = sp.initialResidual();
    const scalar finRes = sp.finalResidual();

    logFile << optIter << ","
            << timeName << ","
            << blockName << ","
            << outerIter << ","
            << eqnName << ","
            << initRes << ","
            << finRes << ","
            << sp.nIterations() << ","
            << (sp.converged() ? 1 : 0) << ","
            << (sp.singular() ? 1 : 0)
            << "\n";
}

#endif
