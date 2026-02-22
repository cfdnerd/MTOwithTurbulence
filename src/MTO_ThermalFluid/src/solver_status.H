/*---------------------------------------------------------------------------*\
  solver_status.H - Solver performance logging for MTO_ThermalFluid
  Accumulates residuals/iterations per opt, writes summary at end of each opt
\*---------------------------------------------------------------------------*/

#ifndef solver_status_H
#define solver_status_H

#include "OFstream.H"
#include "HashTable.H"
#include <sstream>

namespace MTOsolverStatus
{
    static HashTable<scalar> maxInitRes;
    static HashTable<scalar> maxFinRes;
    static HashTable<label> maxNIter;
    static HashTable<label> nSolves;
    static HashTable<label> nConverged;
    static HashTable<label> nSingular;
    static std::string meshFilterHeaderContent;
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
    const label nIter = cmptMax(sp.nIterations());

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

//- Store mesh/filter header for deferred write (no file I/O at init; avoids stall on clusters)
inline void storeMeshFilterHeaderContent
(
    const label nCell,
    const scalar volTot,
    const scalar Uref,
    const scalar area,
    const scalar N,
    const scalar len,
    const scalar filterR,
    const scalar filterR2,
    const scalar GeoDim,
    const scalar nuVal
)
{
    if (!Pstream::master()) return;

    const scalar hChar = Foam::pow(Foam::max(volTot/nCell, 1e-30), 1.0/3.0);
    const scalar lenCell = (GeoDim >= 2) ? std::pow(Foam::max(area/N, 1e-30), 1.0/GeoDim) : hChar;
    const scalar radius1 = filterR * len;
    const scalar radius2 = filterR2 * len;
    const scalar cellsInRadius1 = Foam::max(radius1/lenCell, 0.5);
    const scalar cellsInRadius2 = Foam::max(radius2/lenCell, 0.5);
    const scalar filterR_suggested_min = Foam::max(1.5, 1.5 * lenCell / Foam::max(len, 1e-30));
    const scalar filterR_suggested_max = Foam::max(2.5, 2.5 * lenCell / Foam::max(len, 1e-30));
    const scalar Re_est = Uref * len / nuVal;

    std::ostringstream os;
    os.precision(6);
    os.setf(std::ios::scientific);
    os << "# ================================================================================\n";
    os << "# MTO_ThermalFluid - Mesh & Filter Estimates (deferred write at opt 1)\n";
    os << "# ================================================================================\n";
    os << "Mesh:  nCells " << nCell << "  area " << area << "  N_design " << N
       << "  h_char " << hChar << "  len_cell " << lenCell << "\n";
    os << "Filter:  len " << len << "  filterR " << filterR << "  filterR2 " << filterR2
       << "  GeoDim " << GeoDim << "\n";
    os << "  (len = characteristic filter length scale [m]; physical radius = filterR*len)\n";
    os << "Filter radius (physical):  stage1 " << radius1 << "  stage2 " << radius2
       << "  (~" << cellsInRadius1 << ", " << cellsInRadius2 << " cells)\n";
    os << "Suggested filterR for mesh (less gray):  " << filterR_suggested_min
       << " ... " << filterR_suggested_max << "  (rule-of-thumb: 2-3 cell radii)\n";
    os << "# --------------------------------------------------------------------------------\n";
    os << "# Turbulence mesh criteria (RAS k-omega SST):\n";
    os << "  Re_est " << Re_est
       << "  | U_ref*len/nu at design-cell scale; bulk Re uses channel Dh instead\n";
    os << "  y+ < 5 wall-resolved:  first cell delta_y < 50*nu/U_bulk (approx)\n";
    os << "    | delta_y in wall-normal dir; for fully resolved BL use h < l_t/2 below\n";
    os << "  Turbulent length scale l_t ~ 0.07*L (L=BL thickness or channel half-height)\n";
    os << "    | Resolve boundary layers with cell size h < l_t/2 to capture shear\n";
    os << "# ================================================================================\n\n";

    MTOsolverStatus::meshFilterHeaderContent = os.str();
}

//- Write accumulated solver summary to solver_status.log and clear
inline void writeSolverStatusSummary(const Time& runTime, const label optIter)
{
    if (!Pstream::master()) return;

    using namespace MTOsolverStatus;
    if (maxInitRes.empty()) return;

    const fileName logPath(runTime.rootPath()/runTime.globalCaseName()/"solver_status.log");
    // opt 1 with stored mesh header: overwrite file (header at top, then solver table); otherwise append
    const bool appendMode = (optIter > 1) || meshFilterHeaderContent.empty();
    OFstream logFile(logPath, IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED, appendMode);
    if (!logFile.good()) return;

    logFile.setf(std::ios::scientific);
    logFile.precision(6);

    // Write mesh/filter header first if stored (deferred from init to avoid cluster stall)
    if (optIter == 1 && !meshFilterHeaderContent.empty())
    {
        logFile << meshFilterHeaderContent;
        meshFilterHeaderContent.clear();
    }

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
