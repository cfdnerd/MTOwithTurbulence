#include <iomanip>
Info<<"cost function"<<endl;
MeanT=fvc::domainIntegrate(T).value()/area;
V=fvc::domainIntegrate(xh).value()/area-(voluse+0*Foam::max(0,0.5-opt*0.015));

// Actual (raw) power dissipation before normalization
scalar PowerDissActual = 0;
for (i=0; i<nObjPatch; i++)
{
    PowerDissActual = PowerDissActual - sum(
        phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
        + 0.5*magSqr(U.boundaryField()[conPatchList[i]])));
}
reduce(PowerDissActual, sumOp<scalar>());

PowerDiss=PowerDissActual/PowerDiss0;
if(GeoDim==2)
{
    PowerDiss=PowerDiss/std::pow(mesh.V()[0],1.0/3);
}
reduce(PowerDiss, sumOp<scalar>());

if (Pstream::master())
{
    const scalar elapsedTime = runTime.elapsedClockTime();
    ofstream optLog("optimization.hst", std::ios::app);
    optLog.setf(std::ios::scientific);
    optLog.precision(6);
    if (opt == 1)
    {
        optLog << "# ================================================================================\n";
        optLog << "# MTO_ThermalFluid Optimization History\n";
        optLog << "# ================================================================================\n";
        optLog << "#  optIter  | PowerDissActual [W] | PowerDiss [-] |  MeanT [K]  |  maxT [K]  |     V [-]     | elapsedTime [s]\n";
        optLog << "# ---------+---------------------+---------------+-------------+------------+---------------+------------------\n";
    }
    optLog << std::setw(9) << opt << " | "
           << std::setw(19) << PowerDissActual << " | "
           << std::setw(13) << PowerDiss << " | "
           << std::setw(11) << MeanT << " | "
           << std::setw(10) << maxT << " | "
           << std::setw(15) << V << " | "
           << std::setw(16) << elapsedTime << "\n";
    optLog.close();
}
if(runTime.writeTime())
{
   T.write();
   //U.write();
}