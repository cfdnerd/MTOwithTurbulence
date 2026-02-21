// Cascaded sensitivity chain rule: dJ/dx via two PDE adjoint solves
// Normalize by max over DESIGN DOMAIN only; fixed-zone sensitivities can be orders
// of magnitude larger (e.g. inlet/outlet) and would otherwise squash design updates
{
    volScalarField tmpF(mag(fsenshMeanT));
    volScalarField tmpG(mag(gsenshPowerDiss));
    if(solid_area) { forAll(cells_solid,i) { tmpF[cells_solid[i]]=0; tmpG[cells_solid[i]]=0; } }
    if(fluid_area) { forAll(cells_fluid,i) { tmpF[cells_fluid[i]]=0; tmpG[cells_fluid[i]]=0; } }
    // Reference: do not exclude test_area from normalization; test zone participates in design
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
// Reference: do not zero test_area sensitivities; test zone participates in MMA update
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
// Stage 2: linear PDE uses fsenshMeanT as RHS; nonlinear fW-mean uses modified RHS and chain rule
// sensTmp1/2/3 must be in fvSolution solvers
{
    volScalarField sensTmp1("sensTmp1", fsenshMeanT);
    if(filterStage2Nonlinear)
    {
        const scalar pMean = Foam::max(filterPMean, 1.01);
        const scalar xpFloor = 1e-12;
        volScalarField dJdxpRaw("dJdxpRaw", fsenshMeanT * (1.0/pMean) * Foam::pow(Foam::max(xp, scalar(xpFloor)), 1.0 - pMean));
        fvScalarMatrix sens1Eqn(fvm::laplacian(sensTmp1)-fvm::Sp(bFilter2,sensTmp1)+dJdxpRaw*bFilter2);
        SolverPerformance<scalar> spS1 = sens1Eqn.solve();
        logSolverStats(runTime, "filter_chainrule", 0, "fsensStage2_fW", spS1, opt, runTime.timeName());
        fsenshMeanT = sensTmp1 * pMean * Foam::pow(Foam::max(xp1, scalar(xpFloor)), pMean - 1.0);
    }
    else
    {
        fvScalarMatrix sens1Eqn(fvm::laplacian(sensTmp1)-fvm::Sp(bFilter2,sensTmp1)+fsenshMeanT*bFilter2);
        SolverPerformance<scalar> spS1 = sens1Eqn.solve();
        logSolverStats(runTime, "filter_chainrule", 0, "fsensStage2", spS1, opt, runTime.timeName());
        fsenshMeanT = sensTmp1;
    }
}
{
    fvScalarMatrix fsensEqn(fvm::laplacian(fsensMeanT)-fvm::Sp(b,fsensMeanT)+fsenshMeanT*b);
    SolverPerformance<scalar> spFsens = fsensEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 1, "fsensMeanT", spFsens, opt, runTime.timeName());
}
{
    volScalarField sensTmp2("sensTmp2", gsenshPowerDiss);
    if(filterStage2Nonlinear)
    {
        const scalar pMean = Foam::max(filterPMean, 1.01);
        const scalar xpFloor = 1e-12;
        volScalarField dJdxpRaw("dJdxpRaw2", gsenshPowerDiss * (1.0/pMean) * Foam::pow(Foam::max(xp, scalar(xpFloor)), 1.0 - pMean));
        fvScalarMatrix sens2Eqn(fvm::laplacian(sensTmp2)-fvm::Sp(bFilter2,sensTmp2)+dJdxpRaw*bFilter2);
        SolverPerformance<scalar> spS2 = sens2Eqn.solve();
        logSolverStats(runTime, "filter_chainrule", 0, "gsensPStage2_fW", spS2, opt, runTime.timeName());
        gsenshPowerDiss = sensTmp2 * pMean * Foam::pow(Foam::max(xp1, scalar(xpFloor)), pMean - 1.0);
    }
    else
    {
        fvScalarMatrix sens2Eqn(fvm::laplacian(sensTmp2)-fvm::Sp(bFilter2,sensTmp2)+gsenshPowerDiss*bFilter2);
        SolverPerformance<scalar> spS2 = sens2Eqn.solve();
        logSolverStats(runTime, "filter_chainrule", 0, "gsensPStage2", spS2, opt, runTime.timeName());
        gsenshPowerDiss = sensTmp2;
    }
}
{
    fvScalarMatrix gsensEqn(fvm::laplacian(gsensPowerDiss)-fvm::Sp(b,gsensPowerDiss)+gsenshPowerDiss*b);
    SolverPerformance<scalar> spGsensP = gsensEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 1, "gsensPowerDiss", spGsensP, opt, runTime.timeName());
}
{
    volScalarField sensTmp3("sensTmp3", gsenshVol);
    if(filterStage2Nonlinear)
    {
        const scalar pMean = Foam::max(filterPMean, 1.01);
        const scalar xpFloor = 1e-12;
        volScalarField dJdxpRaw("dJdxpRaw3", gsenshVol * (1.0/pMean) * Foam::pow(Foam::max(xp, scalar(xpFloor)), 1.0 - pMean));
        fvScalarMatrix sens3Eqn(fvm::laplacian(sensTmp3)-fvm::Sp(bFilter2,sensTmp3)+dJdxpRaw*bFilter2);
        SolverPerformance<scalar> spS3 = sens3Eqn.solve();
        logSolverStats(runTime, "filter_chainrule", 0, "gsensVStage2_fW", spS3, opt, runTime.timeName());
        gsenshVol = sensTmp3 * pMean * Foam::pow(Foam::max(xp1, scalar(xpFloor)), pMean - 1.0);
    }
    else
    {
        fvScalarMatrix sens3Eqn(fvm::laplacian(sensTmp3)-fvm::Sp(bFilter2,sensTmp3)+gsenshVol*bFilter2);
        SolverPerformance<scalar> spS3 = sens3Eqn.solve();
        logSolverStats(runTime, "filter_chainrule", 0, "gsensVStage2", spS3, opt, runTime.timeName());
        gsenshVol = sensTmp3;
    }
}
{
    fvScalarMatrix gvolEqn(fvm::laplacian(gsensVol)-fvm::Sp(b,gsensVol)+gsenshVol*b);
    SolverPerformance<scalar> spGsensV = gvolEqn.solve();
    logSolverStats(runTime, "filter_chainrule", 1, "gsensVol", spGsensV, opt, runTime.timeName());
}

fsensMeanT.primitiveFieldRef()=fsensMeanT.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensPowerDiss.primitiveFieldRef()=gsensPowerDiss.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensVol.primitiveFieldRef()=gsensVol.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
