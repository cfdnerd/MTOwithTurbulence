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
    // Use global domain area (not mesh.V()[0]) for parallel consistency and stable scaling
    PowerDiss=PowerDiss/std::pow(area,1.0/3);
}

if (Pstream::master())
{
    const scalar elapsedTime = runTime.elapsedClockTime();
    const fileName histFile(runTime.rootPath()/runTime.globalCaseName()/"optimization.hst");
    ofstream optLog(histFile.c_str(), std::ios::app);
    if (!optLog.good())
    {
        FatalErrorInFunction
            << "Cannot open " << histFile << " for writing"
            << exit(FatalError);
    }
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
    optLog.width(9);  optLog << opt << " | ";
    optLog.width(19); optLog << PowerDissActual << " | ";
    optLog.width(13); optLog << PowerDiss << " | ";
    optLog.width(11); optLog << MeanT << " | ";
    optLog.width(10); optLog << maxT << " | ";
    optLog.width(15); optLog << V << " | ";
    optLog.width(16); optLog << elapsedTime << "\n";
    optLog.close();
}
if(runTime.writeTime())
{
   T.write();
   //U.write();
}