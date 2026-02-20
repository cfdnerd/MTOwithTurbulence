/*---------------------------------------------------------------------------*\
  solver_diagnostics.H - Per-rank runtime diagnostics for crash forensics
\*---------------------------------------------------------------------------*/

#ifndef solver_diagnostics_H
#define solver_diagnostics_H

#include "OFstream.H"
#include "Pstream.H"
#include <cmath>

namespace MTOdiag
{
inline fileName diagLogPath(const Time& runTime)
{
    return runTime.rootPath()/runTime.globalCaseName()
        /("solver_diagnostics_rank" + Foam::name(Pstream::myProcNo()) + ".log");
}

inline OFstream openDiag(const Time& runTime)
{
    return OFstream
    (
        diagLogPath(runTime),
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED,
        true
    );
}

inline void writeHeaderIfNeeded(const Time& runTime)
{
    static bool wroteHeader = false;
    if (wroteHeader) return;
    OFstream os = openDiag(runTime);
    if (!os.good()) return;
    os << "# ================================================================================\n";
    os << "# MTO_ThermalFluid diagnostics (per-rank)\n";
    os << "# rank=" << Pstream::myProcNo() << "\n";
    os << "# columns: optIter stage metric value\n";
    os << "# ================================================================================\n";
    wroteHeader = true;
}

inline void logMetric(const Time& runTime, const label optIter, const word& stage, const word& metric, const scalar value)
{
    writeHeaderIfNeeded(runTime);
    OFstream os = openDiag(runTime);
    if (!os.good()) return;
    os.setf(std::ios::scientific);
    os.precision(6);
    os << optIter << " " << stage << " " << metric << " " << value << "\n";
}

inline void scalarStats(const UList<scalar>& f, label& nBad, scalar& minVal, scalar& maxVal)
{
    nBad = 0;
    minVal = GREAT;
    maxVal = -GREAT;
    bool anyFinite = false;
    forAll(f, i)
    {
        const scalar v = f[i];
        if (std::isfinite(v))
        {
            anyFinite = true;
            minVal = Foam::min(minVal, v);
            maxVal = Foam::max(maxVal, v);
        }
        else
        {
            nBad++;
        }
    }
    if (!anyFinite)
    {
        minVal = 0.0;
        maxVal = 0.0;
    }
}

inline void vectorMagStats(const UList<vector>& f, label& nBad, scalar& minMag, scalar& maxMag)
{
    nBad = 0;
    minMag = GREAT;
    maxMag = -GREAT;
    bool anyFinite = false;
    forAll(f, i)
    {
        const vector& v = f[i];
        if (std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z()))
        {
            const scalar m = mag(v);
            anyFinite = true;
            minMag = Foam::min(minMag, m);
            maxMag = Foam::max(maxMag, m);
        }
        else
        {
            nBad++;
        }
    }
    if (!anyFinite)
    {
        minMag = 0.0;
        maxMag = 0.0;
    }
}

inline void logScalarFieldHealth(const Time& runTime, const label optIter, const word& stage, const word& name, const UList<scalar>& fld)
{
    label nBad = 0;
    scalar minVal = 0.0, maxVal = 0.0;
    scalarStats(fld, nBad, minVal, maxVal);
    logMetric(runTime, optIter, stage, name + "_min", minVal);
    logMetric(runTime, optIter, stage, name + "_max", maxVal);
    logMetric(runTime, optIter, stage, name + "_nBad", scalar(nBad));
}

inline void logVectorFieldHealth(const Time& runTime, const label optIter, const word& stage, const word& name, const UList<vector>& fld)
{
    label nBad = 0;
    scalar minMag = 0.0, maxMag = 0.0;
    vectorMagStats(fld, nBad, minMag, maxMag);
    logMetric(runTime, optIter, stage, name + "_magMin", minMag);
    logMetric(runTime, optIter, stage, name + "_magMax", maxMag);
    logMetric(runTime, optIter, stage, name + "_nBad", scalar(nBad));
}
}

#endif
