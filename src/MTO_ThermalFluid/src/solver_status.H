/*---------------------------------------------------------------------------*\
  solver_status.H - Solver performance logging for MTO_ThermalFluid
  Accumulates residuals/iterations per opt, writes summary at end of each opt
\*---------------------------------------------------------------------------*/

#ifndef solver_status_H
#define solver_status_H

#include "OFstream.H"
#include "HashTable.H"

namespace MTOsolverStatus
{
    static HashTable<scalar> maxInitRes;
    static HashTable<scalar> maxFinRes;
    static HashTable<label> maxNIter;
    static HashTable<label> nSolves;
    static HashTable<label> nConverged;
    static HashTable<label> nSingular;
}

//- Accumulate solver stats (master rank only)
inline void logSolverStats
(
    const Time& runTime,
    const word& blockName,
    const label outerIter,
    const word& eqnName,
    const SolverPerformance<scalar>& sp,
    const label optIter,
    const word& timeName
)
{
    if (!Pstream::master()) return;

    const word key(eqnName);
    const scalar initRes = sp.initialResidual();
    const scalar finRes = sp.finalResidual();
    const label nIter = sp.nIterations();

    using namespace MTOsolverStatus;
    if (!maxInitRes.found(key))
    {
        maxInitRes.insert(key, initRes);
        maxFinRes.insert(key, finRes);
        maxNIter.insert(key, nIter);
        nSolves.insert(key, 1);
        nConverged.insert(key, sp.converged() ? 1 : 0);
        nSingular.insert(key, sp.singular() ? 1 : 0);
    }
    else
    {
        maxInitRes[key] = Foam::max(maxInitRes[key], initRes);
        maxFinRes[key] = Foam::max(maxFinRes[key], finRes);
        maxNIter[key] = Foam::max(maxNIter[key], nIter);
        nSolves[key]++;
        if (sp.converged()) nConverged[key]++;
        if (sp.singular()) nSingular[key]++;
    }
}

inline void logSolverStats
(
    const Time& runTime,
    const word& blockName,
    const label outerIter,
    const word& eqnName,
    const SolverPerformance<vector>& sp,
    const label optIter,
    const word& timeName
)
{
    if (!Pstream::master()) return;

    const word key(eqnName);
    const scalar initRes = mag(sp.initialResidual());
    const scalar finRes = mag(sp.finalResidual());
    const label nIter = sp.nIterations();

    using namespace MTOsolverStatus;
    if (!maxInitRes.found(key))
    {
        maxInitRes.insert(key, initRes);
        maxFinRes.insert(key, finRes);
        maxNIter.insert(key, nIter);
        nSolves.insert(key, 1);
        nConverged.insert(key, sp.converged() ? 1 : 0);
        nSingular.insert(key, sp.singular() ? 1 : 0);
    }
    else
    {
        maxInitRes[key] = Foam::max(maxInitRes[key], initRes);
        maxFinRes[key] = Foam::max(maxFinRes[key], finRes);
        maxNIter[key] = Foam::max(maxNIter[key], nIter);
        nSolves[key]++;
        if (sp.converged()) nConverged[key]++;
        if (sp.singular()) nSingular[key]++;
    }
}

//- Write accumulated solver summary to solver_status.log and clear
inline void writeSolverStatusSummary(const Time& runTime, const label optIter)
{
    if (!Pstream::master()) return;

    using namespace MTOsolverStatus;
    if (maxInitRes.empty()) return;

    const fileName logPath(runTime.rootPath()/runTime.globalCaseName()/"solver_status.log");
    OFstream logFile(logPath, IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED, true);
    if (!logFile.good()) return;

    logFile.setf(std::ios::scientific);
    logFile.precision(6);

    if (optIter == 1)
    {
        logFile << "# ================================================================================\n";
        logFile << "# MTO_ThermalFluid Solver Status (critical residuals per optimization iteration)\n";
        logFile << "# ================================================================================\n";
        logFile << "# Primal: U,p,T | Adjoint: Ua,pa,Ub,pb,Tb,ka,omegaa | Filter: xp,fsensMeanT,gsensPowerDiss,gsensVol\n";
        logFile << "# (k,omega from turbulence->correct() not logged)\n";
        logFile << "# --------------------------------------------------------------------------------\n";
        logFile << "#  optIter | equation      | maxInitRes [-] | maxFinRes [-] | maxNIter | nSolves | conv | sing\n";
        logFile << "# ---------+---------------+----------------+---------------+----------+---------+------+-----\n";
    }

    // Ordered list of equations for readable output
    const wordList eqnOrder
    (
        {"U", "p", "T", "Ua", "pa", "Ub", "pb", "Tb", "ka", "omegaa",
         "xp", "fsensMeanT", "gsensPowerDiss", "gsensVol"}
    );

    forAll(eqnOrder, i)
    {
        const word& eqn = eqnOrder[i];
        if (!maxInitRes.found(eqn)) continue;

        logFile << "  ";
        logFile.width(7);  logFile << optIter << " | ";
        logFile.width(13); logFile << eqn << " | ";
        logFile.width(14); logFile << maxInitRes[eqn] << " | ";
        logFile.width(13); logFile << maxFinRes[eqn] << " | ";
        logFile.width(8);  logFile << maxNIter[eqn] << " | ";
        logFile.width(7);  logFile << nSolves[eqn] << " | ";
        logFile.width(4);  logFile << nConverged[eqn] << " | ";
        logFile.width(4);  logFile << nSingular[eqn] << "\n";
    }
    forAllConstIter(HashTable<scalar>, maxInitRes, iter)
    {
        const word& eqn = iter.key();
        bool inOrder = false;
        forAll(eqnOrder, j) { if (eqnOrder[j] == eqn) { inOrder = true; break; } }
        if (inOrder) continue;  // already printed

        logFile << "  ";
        logFile.width(7);  logFile << optIter << " | ";
        logFile.width(13); logFile << eqn << " | ";
        logFile.width(14); logFile << maxInitRes[eqn] << " | ";
        logFile.width(13); logFile << maxFinRes[eqn] << " | ";
        logFile.width(8);  logFile << maxNIter[eqn] << " | ";
        logFile.width(7);  logFile << nSolves[eqn] << " | ";
        logFile.width(4);  logFile << nConverged[eqn] << " | ";
        logFile.width(4);  logFile << nSingular[eqn] << "\n";
    }

    logFile << "# ---------\n";

    // Clear for next opt
    maxInitRes.clear();
    maxFinRes.clear();
    maxNIter.clear();
    nSolves.clear();
    nConverged.clear();
    nSingular.clear();
}

#endif
