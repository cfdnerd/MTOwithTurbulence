fsenshMeanT=fsenshMeanT/gMax(mag(fsenshMeanT.primitiveField()));
gsenshPowerDiss=gsenshPowerDiss/gMax(mag(gsenshPowerDiss.primitiveField()));
gsenshVol.primitiveFieldRef() = 1.0;
//***********************************//Heaviside_filter
if(opt>1)
{
  for(i=0;i<n;i++)
  {
     if(xp[i]<=eta5)
     {
       drho[i]=del*Foam::exp(-del*(1-xp[i]/eta5))+Foam::exp(-del); 
     }
     else
     {
       drho[i]=del*Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+Foam::exp(-del); 
     } 
  }
  fsenshMeanT=fsenshMeanT*drho;
  gsenshPowerDiss=gsenshPowerDiss*drho;
  gsenshVol=gsenshVol*drho;
}
//***********************************//PDE filter
{
    fvScalarMatrix fsensEqn(fvm::laplacian(fsensMeanT)-fvm::Sp(b,fsensMeanT)+fsenshMeanT*b);
    SolverPerformance<scalar> spFsens = fsensEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 0, "fsensMeanT", spFsens, opt, runTime.timeName());
}
{
    fvScalarMatrix gsensEqn(fvm::laplacian(gsensPowerDiss)-fvm::Sp(b,gsensPowerDiss)+gsenshPowerDiss*b);
    SolverPerformance<scalar> spGsensP = gsensEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 0, "gsensPowerDiss", spGsensP, opt, runTime.timeName());
}
{
    fvScalarMatrix gvolEqn(fvm::laplacian(gsensVol)-fvm::Sp(b,gsensVol)+gsenshVol*b);
    SolverPerformance<scalar> spGsensV = gvolEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 0, "gsensVol", spGsensV, opt, runTime.timeName());
}

fsensMeanT.primitiveFieldRef()=fsensMeanT.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensPowerDiss.primitiveFieldRef()=gsensPowerDiss.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensVol.primitiveFieldRef()=gsensVol.primitiveFieldRef()*mesh.V()/gMax(mesh.V());

