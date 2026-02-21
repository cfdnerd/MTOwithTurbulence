//***********************************// Cascaded PDE filter
// Stage 1: x -> xp1 (linear Helmholtz). Stage 2: xp1 -> xp (linear or nonlinear fW-mean). Heaviside: xp -> xh
//***********************************//Stage 1: PDE filter x -> xp1
{
    fvScalarMatrix xp1Eqn(fvm::laplacian(xp1)-fvm::Sp(b,xp1)+x*b);
    SolverPerformance<scalar> spXp1 = xp1Eqn.solve();
    logSolverStats(runTime, "filter_x", 0, "xp1", spXp1, opt, runTime.timeName());
}
//***********************************//Stage 2: xp1 -> xp (linear PDE or nonlinear fW-mean for gray suppression)
if(filterStage2Nonlinear)
{
    // fW-mean: xp = (W*xp1^p)^{1/p}. PDE: solve (L-bFilter2)*tmp = bFilter2*xp1^p, then xp = tmp^{1/p}
    const scalar pMean = Foam::max(filterPMean, 1.01);
    const scalar xp1Floor = 1e-12;
    volScalarField xp1Pow("xp1Pow", Foam::pow(Foam::max(xp1, scalar(xp1Floor)), pMean));
    {
        fvScalarMatrix xpRawEqn(fvm::laplacian(xp)-fvm::Sp(bFilter2,xp)+xp1Pow*bFilter2);
        SolverPerformance<scalar> spXpR = xpRawEqn.solve();
        logSolverStats(runTime, "filter_x", 1, "xp_fWmean", spXpR, opt, runTime.timeName());
    }
    xp.primitiveFieldRef() = Foam::pow(Foam::max(xp.primitiveField(), xp1Floor), 1.0/pMean);
    xp.correctBoundaryConditions();
}
else
{
    fvScalarMatrix xpEqn(fvm::laplacian(xp)-fvm::Sp(bFilter2,xp)+xp1*bFilter2);
    SolverPerformance<scalar> spXp = xpEqn.solve();
    logSolverStats(runTime, "filter_x", 1, "xp", spXp, opt, runTime.timeName());
}
//***********************************//Heaviside projection (turbulence-tuned: delayed sharpening, gentler ramp, lower cap)
del = (opt > heavisideDelStart)
    ? Foam::min(heavisideDelFactor*(opt - heavisideDelStart), heavisideDelCap)
    : heavisideDelInit;
eta0=0.0001;
eta1=1.0;
y0=diff(xp,mesh.V(),del,eta0,n);
reduce(y0, sumOp<scalar>());
do
{
  eta5=(eta0+eta1)/2.0;
  y5=diff(xp,mesh.V(),del,eta5,n);
  reduce(y5, sumOp<scalar>());
  if(y0*y5<0)
  { 
     eta1=eta5;
  }
  else
  { 
     eta0=eta5;
     y0=y5;
  }
} while ((eta1-eta0)>0.0001);

for(i=0;i<n;i++)
{
   if(xp[i]<=eta5)
   {
     xh[i]=eta5*(Foam::exp(-del*(1-xp[i]/eta5))-(1-xp[i]/eta5)*Foam::exp(-del)); 
   }
   else
   {
     xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+(xp[i]-eta5)*Foam::exp(-del)/(1-eta5));
   }   
}
// Overwrite fixed zones (reference MTO: solid/fluid only; xh[test] set at init in createFields, not per iteration)
if(solid_area) { forAll(cells_solid, i) { xh[cells_solid[i]] = 0; } }
if(fluid_area) { forAll(cells_fluid, i) { xh[cells_fluid[i]] = 1.0; } }
