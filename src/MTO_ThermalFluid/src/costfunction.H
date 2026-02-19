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
    if (opt == 1)
    {
        optLog << "# optIter PowerDissActual PowerDiss MeanT maxT V elapsedTime[s]\n";
    }
    optLog << opt << " " << PowerDissActual << " " << PowerDiss << " " << MeanT
           << " " << maxT << " " << V << " " << elapsedTime << "\n";
    optLog.close();
}
if(runTime.writeTime())
{
   T.write();
   //U.write();
}