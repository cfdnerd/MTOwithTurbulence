// Normalize by max over DESIGN DOMAIN only; fixed-zone sensitivities can be orders
// of magnitude larger (e.g. inlet/outlet) and would otherwise squash design updates
{
    volScalarField tmpF(mag(fsenshMeanT));
    volScalarField tmpG(mag(gsenshPowerDiss));
    if(solid_area) { forAll(cells_solid,i) { tmpF[cells_solid[i]]=0; tmpG[cells_solid[i]]=0; } }
    if(fluid_area) { forAll(cells_fluid,i) { tmpF[cells_fluid[i]]=0; tmpG[cells_fluid[i]]=0; } }
    if(test_area)  { forAll(cells_test,i)  { tmpF[cells_test[i]]=0;  tmpG[cells_test[i]]=0;  } }
    scalar maxF = Foam::max(gMax(tmpF), VSMALL);
    scalar maxG = Foam::max(gMax(tmpG), VSMALL);
    fsenshMeanT /= maxF;
    gsenshPowerDiss /= maxG;
}
gsenshVol.primitiveFieldRef() = 1.0;
// Zero sensitivities in fixed non-design zones so MMA does not update them
if(solid_area)
{
   forAll(cells_solid, i)
   {
      fsenshMeanT[cells_solid[i]] = 0;
      gsenshPowerDiss[cells_solid[i]] = 0;
      gsenshVol[cells_solid[i]] = 0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      fsenshMeanT[cells_fluid[i]] = 0;
      gsenshPowerDiss[cells_fluid[i]] = 0;
      gsenshVol[cells_fluid[i]] = 0;
   }
}
if(test_area)
{
   forAll(cells_test, i)
   {
      fsenshMeanT[cells_test[i]] = 0;
      gsenshPowerDiss[cells_test[i]] = 0;
      gsenshVol[cells_test[i]] = 0;
   }
}
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

