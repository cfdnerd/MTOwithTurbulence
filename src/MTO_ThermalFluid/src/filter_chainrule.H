// Cascaded sensitivity chain rule: dJ/dx via two PDE adjoint solves
// Normalize by max over DESIGN DOMAIN only; fixed-zone sensitivities can be orders
// of magnitude larger (e.g. inlet/outlet) and would otherwise squash design updates
{
    volScalarField tmpF(mag(fsenshMeanT));
    volScalarField tmpG(mag(gsenshPowerDiss));
    if(solid_area) { forAll(cells_solid,i) { tmpF[cells_solid[i]]=0; tmpG[cells_solid[i]]=0; } }
    if(fluid_area) { forAll(cells_fluid,i) { tmpF[cells_fluid[i]]=0; tmpG[cells_fluid[i]]=0; } }
    if(test_area)  { forAll(cells_test,i)  { tmpF[cells_test[i]]=0;  tmpG[cells_test[i]]=0;  } }
    const scalar sensFloor = 1e-12;
    scalar maxF = Foam::max(gMax(tmpF), sensFloor);
    scalar maxG = Foam::max(gMax(tmpG), sensFloor);
    fsenshMeanT /= maxF;
    gsenshPowerDiss /= maxG;
}
gsenshVol.primitiveFieldRef() = 1.0;
// Zero sensitivities in fixed non-design zones
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
//***********************************//Heaviside chain rule: dJ/dxh -> dJ/dxp
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
//***********************************//Cascaded PDE adjoints: dJ/dxp -> dJ/dxp1 -> dJ/dx
// Reuse xp1 as temporary (must be in fvSolution solvers)
// Stage 2 adjoint: (laplacian - bFilter2)*xp1 = -bFilter2*dJ/dxp
{
    xp1 = fsenshMeanT;
    fvScalarMatrix sens1Eqn(fvm::laplacian(xp1)-fvm::Sp(bFilter2,xp1)+fsenshMeanT*bFilter2);
    SolverPerformance<scalar> spS1 = sens1Eqn.solve();
    logSolverStats(runTime, "filter_chainrule", 0, "fsensStage2", spS1, opt, runTime.timeName());
    fsenshMeanT = xp1;
}
// Stage 1 adjoint: (laplacian - b)*fsensMeanT = -b*dJ/dxp1
{
    fvScalarMatrix fsensEqn(fvm::laplacian(fsensMeanT)-fvm::Sp(b,fsensMeanT)+fsenshMeanT*b);
    SolverPerformance<scalar> spFsens = fsensEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 1, "fsensMeanT", spFsens, opt, runTime.timeName());
}
// gsensPowerDiss: same cascaded chain
{
    xp1 = gsenshPowerDiss;
    fvScalarMatrix sens2Eqn(fvm::laplacian(xp1)-fvm::Sp(bFilter2,xp1)+gsenshPowerDiss*bFilter2);
    SolverPerformance<scalar> spS2 = sens2Eqn.solve();
    logSolverStats(runTime, "filter_chainrule", 0, "gsensPStage2", spS2, opt, runTime.timeName());
    gsenshPowerDiss = xp1;
}
{
    fvScalarMatrix gsensEqn(fvm::laplacian(gsensPowerDiss)-fvm::Sp(b,gsensPowerDiss)+gsenshPowerDiss*b);
    SolverPerformance<scalar> spGsensP = gsensEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 1, "gsensPowerDiss", spGsensP, opt, runTime.timeName());
}
// gsensVol: same cascaded chain
{
    xp1 = gsenshVol;
    fvScalarMatrix sens3Eqn(fvm::laplacian(xp1)-fvm::Sp(bFilter2,xp1)+gsenshVol*bFilter2);
    SolverPerformance<scalar> spS3 = sens3Eqn.solve();
    logSolverStats(runTime, "filter_chainrule", 0, "gsensVStage2", spS3, opt, runTime.timeName());
    gsenshVol = xp1;
}
{
    fvScalarMatrix gvolEqn(fvm::laplacian(gsensVol)-fvm::Sp(b,gsensVol)+gsenshVol*b);
    SolverPerformance<scalar> spGsensV = gvolEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 1, "gsensVol", spGsensV, opt, runTime.timeName());
}

fsensMeanT.primitiveFieldRef()=fsensMeanT.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensPowerDiss.primitiveFieldRef()=gsensPowerDiss.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensVol.primitiveFieldRef()=gsensVol.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
